
#include "mutate.hpp"
#include "sequence.hpp"
#include "hibf_util.hpp"
//#include "../megax_util/minimiser.hpp"

#include <cmath>
#include <regex>
#include <zlib.h>
#include <omp.h>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <filesystem>
#include <mutex>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp> 
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

#include "../megax_util/minimiser_hashes.hpp"

#include <seqan3/search/fm_index/all.hpp>


using namespace seqan3::literals;
using namespace std::string_literals;


//==================================================================
// Member functions.
//==================================================================

// Definition of the member constructor Mutate::Mutate. 
Mutate::Mutate(const std::string& matrixFilePath){
    
    //loadMatrix(matrixFilePath);
    double start {omp_get_wtime()};
    processMatrix(matrixFilePath);
    double end {omp_get_wtime()};
    std::cout << "[INFO] Elapsed total time for loading and processing the input matrix:" << end - start << "s" << '\n';
}


// Definition of the member destructor Mutate::~Mutate. 
Mutate::~Mutate(){
std::cout<<"[INFO-DEV] Delete constructer of Mutate functionality....."<<std::endl;
}

// Definition of the member destructor Mutate::compareAminoAcids. 
bool Mutate::compareAminoAcids(const std::pair<std::string, int>& aa_1, const std::pair<std::string, int>& aa_2){

    return aa_1.second > aa_2.second;
}


// Definition of the member Mutate::loadMatrix. 
void Mutate::loadMatrix(const std::string& matrixFilePath){

    std::ifstream matrixFile(matrixFilePath);
    if (!matrixFile.good())  throw std::runtime_error("[ERROR] Failed to open matrix file.");

        std::string column;
        std::vector<char> aminoAcid;
        bool knownAminoAcid = false; // avoid reading the whole matrix 

        while (std::getline(matrixFile, column)) {
            if (!column.empty()) {
                std::istringstream iss(column);
                char aminoAcidTemp;
                if (!knownAminoAcid) {
                    while (iss >> aminoAcidTemp) {
                        aminoAcid.push_back(aminoAcidTemp);
                    }
                    knownAminoAcid = true;
                } else {
                    char rowAminoAcid;
                    iss >> rowAminoAcid; 
                    for (char columnAminoAcid : aminoAcid) {
                        int aaScore;
                        iss >> aaScore;
                        int aa_1 = AminoAcidsProcessing::mapCharToIndex(rowAminoAcid);
                        int aa_2 = AminoAcidsProcessing::mapCharToIndex(columnAminoAcid);
                        //std::cout << "row: " << rowAminoAcid << " column: "<< columnAminoAcid << " Score: " << aaScore << std::endl;
                        //std::cout << "row: " << aa_1 << " column: "<< aa_2 << " Score: " << aaScore << std::endl;
                        this->scoreSetter(aa_1, aa_2, aaScore);
                            
                    }
                }
            }
        }

        matrixFile.close();

    }    


// Definition of the member Mutate::processMatrix. 
void Mutate::processMatrix(const std::string& matrixFilePath){

    std::cout<<"\n[INFO] Start loading the input matrix....."<<std::endl;
    std::cout<<"\n[WARNING] We are assuming symmetric matrix as input!"<<std::endl;

    double start {omp_get_wtime()};

    loadMatrix(matrixFilePath);

    double end {omp_get_wtime()};
    std::cout << "[INFO] Elapsed time for loading the matrix:" << end - start << "s" << std::endl;

    std::cout<<"\n[INFO] Start processing the input matrix....."<<std::endl;

    // resize vector size to ensure no memory issues
    this->sortedSubstitutionMatrix.resize(AminoAcidsProcessing::TOTAL_RESIDUES);

    for (unsigned short int i = 0; i < AminoAcidsProcessing::TOTAL_RESIDUES; i++){ // row

        for (unsigned short int j = 0; j < AminoAcidsProcessing::TOTAL_RESIDUES; j++){ // column

            std::pair<std::string, int> rowHolder
            {};
            rowHolder.first.push_back(AminoAcidsProcessing::convertIndexToAminoAcid(j));// get Amino Acid
            rowHolder.second = this->scoreGetter(i, j);
            this->sortedSubstitutionMatrix[i].push_back(rowHolder);
            //std::cout << rowHolder.first << " " << rowHolder.second << std::endl;
        }
    }

    // Sorting the matrix to avoid any low score neighbors
    for (unsigned short int i = 0; i < AminoAcidsProcessing::TOTAL_RESIDUES; i++){

        std::sort(this->sortedSubstitutionMatrix[i].begin(), this->sortedSubstitutionMatrix[i].end(),
          [this](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
              return this->compareAminoAcids(a, b);
          });
    }

}


// Definition of the member destructor Mutate::scoreSetter. 
void Mutate::scoreSetter(int& aa_1, int& aa_2, int score){
    this->substitutionMatrix[std::make_pair(aa_1, aa_2)] = score;
    this->substitutionMatrix[std::make_pair(aa_2, aa_1)] = score;  // Only one scan, matrix is symmetric
}


// Definition of the member Mutate::scoreGetter. 
int Mutate::scoreGetter(int aa_1, int aa_2){
    return (this->substitutionMatrix.find(std::make_pair(aa_1, aa_2)) != this->substitutionMatrix.end()) ?
           this->substitutionMatrix.find(std::make_pair(aa_1, aa_2))->second : 0;
}


// Definition of the member Mutate::calculateQmerScore. 
int Mutate::calculateQmerScore(const std::string &qMer, const std::string &neighbor)
{
    int score {0};
    for(size_t i = 0; i < qMer.size(); i++)
    {
        score += this->scoreGetter(AminoAcidsProcessing::mapCharToIndex(qMer[i]), AminoAcidsProcessing::mapCharToIndex(neighbor[i]));
        //std::cout << "check "<< score << " " << mapCharToIndex(qMer[i]) << " " <<  mapCharToIndex(neighbor[i]) << std::endl;
        // 5 score += this->scoreGetter(qMer[i], neighbor[i]);

    }

    return score;
}


// Definition of the member Mutate::writeDBToFasta. 
void Mutate::writeDBToFasta(uint8_t qMerSize, std::string& sequenceFilePath, int minScore, uint8_t threads, std::string& outputDir){

    std::filesystem::path fileName{sequenceFilePath};

    const std::string outputDB = outputDir + "mutated_" + fileName.filename().string();
    seqan3::sequence_file_output fout{outputDB};

    std::mutex uniqueSetMutex;
    std::unordered_set<std::string> uniqueSet;

    using types = seqan3::type_list<std::vector<seqan3::aa27>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;
    using traits_type = seqan3::sequence_file_input_default_traits_aa;

    seqan3::sequence_file_input<traits_type> fin{sequenceFilePath};
    
    seqan3::debug_stream << "[LOG] Reference File Name: "<< sequenceFilePath << '\n';
    seqan3::debug_stream << "[LOG] Output File Name   : "<< outputDB << '\n';
    seqan3::debug_stream << "[LOG] Output Dir         : "<< outputDir << '\n';
    seqan3::debug_stream << "[LOG] K-mer Size         : "<< unsigned(qMerSize) << '\n';
    seqan3::debug_stream << "[LOG] Mutation Score     : "<< unsigned(minScore) << '\n';
    seqan3::debug_stream << "[LOG] Number of Threads  : "<< unsigned(threads) << '\n';

    seqan3::debug_stream << "[INFO] Start mutating input sequences, you can check the sequences in real time in the output file! "<< '\n';

    unsigned seqCounter {0u};
    unsigned seqCounterOr {0u};

    for (auto & record : fin)
    {
        seqCounterOr++;

        unsigned counter {0u};

        auto const& seq = record.sequence();
        auto const& header = record.id();

        std::string sequence = "";
        for (auto const & c : seq) sequence += seqan3::to_char(c); 

        if (sequence.length() < qMerSize){

            seqan3::debug_stream << "[STEP-ERROR] Sequence is shorter than K-mer, we will skip the sequence: " << sequence << '\n';
            seqan3::debug_stream << "[STEP-ERROR] ID: " << record.id() << '\n';
            continue;
        }

        const unsigned numberOfQMers = (sequence.length() - qMerSize + 1);

        #pragma omp parallel for 
        for(unsigned i = 0; i<numberOfQMers; i++){
            MatrixSearchResults neighboringResults;
            const std::string_view qMer = std::string_view{sequence}.substr(i, qMerSize);

            if (qMer.find('*') != std::string_view::npos) {
                continue;
            }

            neighboringResults = this->scoring(qMer, minScore);
            neighboringResults.neighbors.emplace_back(std::make_pair(qMer,100));

            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                seqan3::aa27_vector qMerHolder{};
                for (const char& c : neighbor.first){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    qMerHolder.push_back(aminoAcid);
                }       

                std::string newSequenceSTD = sequence;
                newSequenceSTD.replace(i, qMerSize, neighbor.first);

                seqan3::aa27_vector newSequence{};
                for (const char& c : newSequenceSTD){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    newSequence.push_back(aminoAcid);
                }  

                if (uniqueSet.find(newSequenceSTD) == uniqueSet.end())
                {
                
                    std::lock_guard<std::mutex> lock(uniqueSetMutex);
                    uniqueSet.insert(newSequenceSTD);
                    counter++;
                    const std::string newHeader = header + " (mut_" + std::to_string(counter) + ")";

                    sequence_record_type record{std::move(newSequence), std::move(newHeader)};
                    seqCounter++;
                    fout.push_back(record);

                }
            }
        }
    }

    seqan3::debug_stream << "[INFO] Finished mutating input sequences! "<< '\n';
    seqan3::debug_stream << "[INFO] Total number of loaded sequences: " << seqCounterOr << '\n';
    seqan3::debug_stream << "[INFO] Total number of written sequences: " << seqCounter << '\n';
}


struct ScoreChecker
{
    bool operator()(const std::pair<std::string, int>& seq1, const std::pair<std::string, int>& seq2) const
    {
        return seq1.second > seq2.second;
    }
};

// Definition of the member Mutate::mutateToSequenceLength
void Mutate::mutateToSequenceLength(uint8_t numberOfTopSequences, std::string& sequenceFilePath, std::string& outputDir, int minScore){


    std::filesystem::path fileName{sequenceFilePath};

    const std::string outputDB = outputDir + "mutated_" + fileName.filename().string();
    seqan3::sequence_file_output fout{outputDB};

    std::mutex uniqueSetMutex;
    std::unordered_set<std::string> uniqueSet;

    using types = seqan3::type_list<std::vector<seqan3::aa27>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;
    using traits_type = seqan3::sequence_file_input_default_traits_aa;

    seqan3::sequence_file_input<traits_type> fin{sequenceFilePath};
    
    seqan3::debug_stream << "[LOG] Reference File Name   : "<< sequenceFilePath << '\n';
    seqan3::debug_stream << "[LOG] Output File Name      : "<< outputDB << '\n';
    seqan3::debug_stream << "[LOG] Output Dir            : "<< outputDir << '\n';
    seqan3::debug_stream << "[LOG] Sequences per Protein : "<< unsigned(numberOfTopSequences) << '\n';
    seqan3::debug_stream << "[INFO] Start mutating input sequences, you can check the sequences in real time in the output file! "<< '\n';


    for (auto & record : fin)
    {
        std::vector<std::pair<std::string, int>> topScoredSequences;
        auto const& seq = record.sequence();
        auto const& header = record.id();

        std::string sequence = "";
        for (auto const & c : seq) sequence += seqan3::to_char(c); 

        if (sequence.length() < 5){ // general case, skip all sequences shorter than 5

            seqan3::debug_stream << "[STEP-ERROR] Sequence is shorter than 5, we will skip the sequence: " << sequence << '\n';
            seqan3::debug_stream << "[STEP-ERROR] ID: " << record.id() << '\n';
            continue;
        }

            MatrixSearchResults neighboringResults;
            const std::string_view qMer = std::string_view{sequence};
            neighboringResults = this->scoring(qMer, minScore);            
            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                topScoredSequences.emplace_back(std::make_pair(neighbor.first, neighbor.second));
    
            }

            std::sort(topScoredSequences.begin(), topScoredSequences.end(), ScoreChecker());

            for (unsigned i = 0; i < numberOfTopSequences; i++)
            {
                const std::string newHeader = header + " (mut_" + std::to_string(topScoredSequences[i].second) + ")";
                seqan3::aa27_vector newSequence{};
                for (const char& c :topScoredSequences[i].first){
                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    newSequence.push_back(aminoAcid);
                }  
                sequence_record_type record{std::move(newSequence), std::move(newHeader)};
                fout.push_back(record);
            }
     
    }

    seqan3::debug_stream << "[INFO] Finished mutating input sequences! "<< '\n';
}

// Definition of the member Mutate::printSequenceStatstics. 
void Mutate::printSequenceStatstics(const std::string& inputFastaFile, uint8_t kMerSize, std::string& outputDir){

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{inputFastaFile};
    auto hashAdaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});

    uint64_t numberOfSequences {0u};
    uint64_t totalNumberOfKMersRaw {0u};

    // output file for mapping sequence name to the target index of counting
    // index starts by 0! 
    const std::string outputMapper = outputDir + "reference_id_map.log";
    const std::string sequenceHistogram = outputDir + "reference_hist.log";
    const std::string outputLengths = outputDir + "lengths.log";

    std::ofstream logFile(outputMapper);
    std::ofstream histFile(sequenceHistogram);
    std::ofstream lengthFile(outputLengths);

    if (!logFile.is_open()) 
      throw std::runtime_error("[ERROR] Couldn't open the file: " + outputMapper);

    if (!histFile.is_open()) 
      throw std::runtime_error("[ERROR] Couldn't open the file: " + sequenceHistogram);

    if (!lengthFile.is_open()) 
      throw std::runtime_error("[ERROR] Couldn't open the file: " + outputLengths);

    std::map<std::string, uint64_t> histogram;

    uint64_t referenceID {0u};
    
    for (auto & record : fin)
    {
        std::vector<uint64_t> hashesBin;

        numberOfSequences += 1; 

        auto const seq = record.sequence();
        auto const id_ = record.id();
        //if (id_ == "[Chlorocebus sabaeus]"){

        //    seqan3::debug_stream << seq << std::endl;
        //}

        if (seq.size() < kMerSize){

            seqan3::debug_stream << "[STEP-ERROR] Sequence is shorter than K-mer, we will skip the sequence: " << seq << '\n';
            seqan3::debug_stream << "[STEP-ERROR] ID: " << record.id() << '\n';
            continue;
        }
        histogram[id_] = seq.size();
        
        auto hashes = (seq | hashAdaptor);
        totalNumberOfKMersRaw += hashes.size();

        lengthFile << seq.size() << '\n';
        logFile << referenceID << ' ' << id_ << '\n';

        referenceID += 1;
    }

    lengthFile.close();
    logFile.close();

    std::vector<std::pair<std::string, uint64_t>> histogramLambda(histogram.begin(), histogram.end());
    std::sort(histogramLambda.begin(), histogramLambda.end(), [](const auto& left, const auto& right) {
        return left.second > right.second;
    });

    seqan3::debug_stream << "[INFO] Writing histogram data to: " << sequenceHistogram << '\n';
    for (const auto& value : histogramLambda) {
        histFile << value.first << ' ' << value.second << '\n';
    }

    histFile.close();
    

    std::size_t tMax = std::ceil(std::sqrt(numberOfSequences)/ 64) * 64;

    std::cout << "###########################################################################" << '\n';
    std::cout << "# Input sequence                   : " << inputFastaFile << '\n';
    std::cout << "# Mapper file name                 : " << outputMapper << '\n';
    std::cout << "# K-mer size                       : "<< unsigned(kMerSize) << '\n';
    std::cout << "# Total number of k-mers (original): " << unsigned(totalNumberOfKMersRaw) << '\n';
    std::cout << "# Total number of sequences        : " << unsigned(numberOfSequences) << '\n';
    std::cout << "# Number of technical bins         : " << unsigned(tMax) << '\n';
    std::cout << "###########################################################################" << '\n';

}

// Definition of the member Mutate::printMutationStatstics. 
void Mutate::printMutationStatstics(const std::string& inputFastaFile, uint8_t kMerSize, int minScore, bool minimiser){

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{inputFastaFile};

    auto hashAdaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});


    uint64_t numberOfSequences {0u};
    uint64_t totalNumberOfKMersRaw {0u};
    uint64_t totalNumberOfKMersMutated {0u};

    
    for (auto & record : fin)
    {
        auto const seq = record.sequence();
        auto const id_ = record.id();
        numberOfSequences += 1; 

        std::vector<uint64_t> hashesOr; // set of original hash values
        std::vector<uint64_t> hashesMu; // set of mutated hash values

        std::string sequence;
        for (auto const & c : seq) sequence += seqan3::to_char(c);

        const unsigned numberOfQMers = (sequence.length() - kMerSize + 1);

        for(unsigned i = 0; i<numberOfQMers; i++){

            MatrixSearchResults neighboringResults;
            const std::string_view qMer = std::string_view{sequence}.substr(i, kMerSize);
            const std::string kMer = std::string{sequence}.substr(i, kMerSize);
            
            seqan3::aa27_vector kmerHolder{};
            for (const char& c : kMer){
                seqan3::aa27 aminoAcid = {};
                aminoAcid.assign_char(c);
                kmerHolder.push_back(aminoAcid);
            }

            auto hashValue = this->hashQMer(kmerHolder, kMerSize);
            hashesOr.push_back(hashValue);

            neighboringResults = this->scoring(qMer, minScore);

            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                /// Convert neighbor string to seqan3::aa27 
                seqan3::aa27_vector qMerHolder{};
                for (const char& c : neighbor.first){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    qMerHolder.push_back(aminoAcid);
                }
                auto hashValue = this->hashQMer(qMerHolder, kMerSize);
                hashesMu.push_back(hashValue);
            }
         }

         if(minimiser){// ToDo: check if minimiser might give also the same results in our scenario! 

            seqan3::debug_stream << "[NOTE] The minimiser is computed for each user-bin! " << std::endl;

            size_t windowSize = 6 - kMerSize + 1; 
            seqan3::debug_stream << "[INFO] Using the window size 6 as default value! " << std::endl;
            seqan3::debug_stream << hashesOr << std::endl;
            auto minimisers = hashesOr | seqan3::views::minimiser(windowSize);
            hashesOr.clear();
            seqan3::debug_stream << hashesOr << std::endl;
            for (auto hash : minimisers){
                seqan3::debug_stream << uint64_t(hash) << std::endl;
            }
            //seqan3::debug_stream << typeid(minimisers).name()  << std::endl;
            seqan3::debug_stream << hashesOr << std::endl;
         }

         std::sort(hashesOr.begin(), hashesOr.end());
         auto uniqueValuesOr = std::unique(hashesOr.begin(), hashesOr.end());
         hashesOr.erase(uniqueValuesOr, hashesOr.end());

         std::sort(hashesMu.begin(), hashesMu.end());
         auto uniqueValuesMu = std::unique(hashesMu.begin(), hashesMu.end());
         hashesMu.erase(uniqueValuesMu, hashesMu.end());

        // Remove elements from hashesMu that are present in hashesOr
        auto it = std::remove_if(hashesMu.begin(), hashesMu.end(),
            [&hashesOr](uint64_t hashMu) {
                return std::binary_search(hashesOr.begin(), hashesOr.end(), hashMu);
            });

        hashesMu.erase(it, hashesMu.end());

        totalNumberOfKMersRaw += hashesOr.size();
        totalNumberOfKMersMutated += hashesMu.size();

        //seqan3::debug_stream << hashesOr << std::endl;
        //seqan3::debug_stream << hashesMu << std::endl;
    }

    uint64_t sumOfOrMu = totalNumberOfKMersRaw + totalNumberOfKMersMutated;

    std::cout << "###########################################################################" << '\n';
    std::cout << "# Input sequence                          : " << inputFastaFile << '\n';
    std::cout << "# K-mer size                              : "<< unsigned(kMerSize) << '\n';
    std::cout << "# Mutation Score                          : "<< unsigned(minScore) << '\n';
    std::cout << "# Total number of original k-mers (unique): " << totalNumberOfKMersRaw << '\n';
    std::cout << "# Total number of mutated k-mers (unique) : " << totalNumberOfKMersMutated << '\n';
    std::cout << "# Total number of k-mers                  : " << sumOfOrMu << '\n';
    std::cout << "# Total number of sequences               : " << unsigned(numberOfSequences) << '\n';
    std::cout << "###########################################################################" << '\n';

}

// Definition of the member Mutate::scoring. 
MatrixSearchResults Mutate::scoring(const std::string_view& qMer, int minScore){

    MatrixSearchResults neighboringResults; // output struct to store the results

    // get index of first amino acid in the input qMer
    unsigned indexFirstAA = AminoAcidsProcessing::mapCharToIndex(qMer[0]);
    if (qMer.size() == 1){

        for (unsigned short int i = 0; i < AminoAcidsProcessing::TOTAL_RESIDUES; i++){

            int score = this->sortedSubstitutionMatrix[indexFirstAA][i].second; // get score of first amino acid from the qMer
            if(score >= minScore){
                neighboringResults.neighbors.emplace_back(this->sortedSubstitutionMatrix[indexFirstAA][i]);
            }
            else{

                break;
            }
        }
    }

    else{// other amino acids

        for (unsigned short int i = 0; i < AminoAcidsProcessing::TOTAL_RESIDUES; i++){

            int score = this->sortedSubstitutionMatrix[indexFirstAA][i].second;
            
            std::string qMerHolder 
            {qMer};
            qMerHolder.erase(0,1); // get other amino acids MNA -> NA

            int newScore = minScore - score;
            MatrixSearchResults recursionCall = this->scoring(qMerHolder, newScore);

            unsigned recursionOutputSize = recursionCall.neighbors.size();
            if (recursionOutputSize > 0){

                for (unsigned j = 0; j<(recursionOutputSize); ++j){
                    const auto neighbors = this->sortedSubstitutionMatrix[indexFirstAA][i].first + recursionCall.neighbors[j].first;
                    auto scores = score + recursionCall.neighbors[j].second;
                    //std::cout << neighbors << " " << scores << std::endl;
                    neighboringResults.neighbors.emplace_back(std::make_pair(neighbors,scores));
                }
            }
            else{

                break;
            }
        }

    }

    return neighboringResults;
}

// Definition of the member Mutate::mutateSequencesRec. 
void Mutate::mutateSequencesRec(uint8_t qMerSize, std::string& sequenceFilePath, int minScore, uint8_t threads, std::string& outputDir, bool minimiser, uint8_t windowSize){

    QMersHashes hashes{};
    
    if(minimiser){
        uint8_t minimiserWindow = windowSize - qMerSize + 1;
        seqan3::debug_stream << "[INFO] Minimiser compuation is applied to the user bins!" << '\n';
        seqan3::debug_stream << "[INFO] Minimiser window size: " << unsigned(minimiserWindow) << '\n';
    }

    const std::string hashesVect = outputDir + "hashes_vect_" + std::to_string(qMerSize) + "_" + std::to_string(minScore) + ".vect";
    std::ofstream  os(hashesVect, std::ios::binary); 
	cereal::BinaryOutputArchive archive(os); 

    seqan3::debug_stream << "[INFO] Start writing binary data base to file (without duplications) ...."<< '\n';
    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{sequenceFilePath};

    omp_set_num_threads(threads);
    std::mutex hashesMutex;

    for (auto & record : fin)
    {
        std::vector <uint64_t> userBin{};

        this->numberOfInputSequences += 1;
        auto const & seq = record.sequence();
        std::string sequence = "";
        for (auto const & c : seq) sequence += seqan3::to_char(c); 

        if (sequence.length() < qMerSize){

            seqan3::debug_stream << "[STEP-ERROR] Sequence is shorter than K-mer, we will skip the sequence: " << sequence << '\n';
            seqan3::debug_stream << "[STEP-ERROR] ID: " << record.id() << '\n';
            continue;
        }

        const unsigned numberOfQMers = (sequence.length() - qMerSize + 1);

        this->totalNumberOfOriginalQMers += numberOfQMers;

        #pragma omp parallel for // run each qMer on safe thread
        for(unsigned i = 0; i<numberOfQMers; i++){
            //const std::string qMer = sequence.substr(i, qMerSize);
            MatrixSearchResults neighboringResults; // each qmer has n neighbors stored in this structure
            const std::string_view qMer = std::string_view{sequence}.substr(i, qMerSize);

            if (qMer.find('*') != std::string_view::npos) {
                continue;
            }
            //std::cout << qMer << std::endl;

            neighboringResults = this->scoring(qMer, minScore);
            // ToDo if QMer contains * then skip! 
            
            //std::cout << "QMer: " << std::string(qMer) << std::endl;
            // Ensure the insertion of original q-mers with score of 100 for the original ones
            neighboringResults.neighbors.emplace_back(std::make_pair(qMer,100));
            
            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                /// Convert neighbor string to seqan3::aa27 
                seqan3::aa27_vector qMerHolder{};
                for (const char& c : neighbor.first){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    qMerHolder.push_back(aminoAcid);
                }
                auto hashValue = this->hashQMer(qMerHolder, qMerSize);
                std::lock_guard<std::mutex> lock(hashesMutex);
                userBin.push_back(hashValue);
                //hashes[userBins].push_back(hashValue);
                // std::cout << "Neighbor: "<< neighbor.first << " | Hash_value: "<< this->hashQMer(qMerHolder, qMerSize) << std::endl;
                // std::cout << "QMer: " << std::string(qMer) << " | Neighbor: "<< neighbor.first << " | Score: "<< neighbor.second << std::endl;
                        
            }
        }

        // Remove duplicates from user bins, which will reduce the size of the assigend user bin 
        // We are sorting the hash values inside the bin and not changing the order of the user bins
        std::sort(userBin.begin(), userBin.end());
        auto uniqueValues = std::unique(userBin.begin(), userBin.end());
        userBin.erase(uniqueValues, userBin.end());
        this->totalNumberOfQMers += userBin.size();

        if(minimiser){

            std::vector<uint64_t> minimiser = computeMinimiserSingleBin(qMerSize, windowSize, userBin);
            archive(minimiser);
        }
        else{

            archive(userBin);
        }
        
        //hashes.push_back(userBin);

    }

    seqan3::debug_stream << "[INFO] Finished mutating input sequences! "<< '\n';
    //seqan3::debug_stream << "[INFO] Start writing histogram log file...." << '\n';

    //uint64_t counter = 0; 
    //const std::string hashesHist = outputDir + "hashes_hist_" + std::to_string(qMerSize) + "_" + std::to_string(minScore) + ".log";

	/*std::ofstream outputFileHis(hashesHist);
	for (auto const& hashVector : hashes){
			
		counter += 1; 
		outputFileHis << counter << " " << hashVector.size() << std::endl;
	}*/

    //seqan3::debug_stream << "[INFO] Finished writing histogram log file!" << '\n';
    
	//archive(hashes);

    seqan3::debug_stream << "[INFO] Finished writing binary dataset to file! "<< '\n';

}


// Definition of the member Mutate::mutateSingleSequence. 
std::string Mutate::mutateSingleSequence(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads){

    std::string mutatedQMers;
    //seqan3::debug_stream << "[ERROR] " << qMerSize << std::endl;
    const unsigned numberOfQMers = (sequence.length() - qMerSize + 1);
    omp_set_num_threads(threads);
    std::mutex sequenceMutex;

    #pragma omp parallel for // run each qMer on safe thread
    for(unsigned i = 0; i<numberOfQMers; i++){
            MatrixSearchResults neighboringResults;
            const std::string_view qMer = std::string_view{sequence}.substr(i, qMerSize);
            if (qMer.find('*') != std::string_view::npos) {
                continue;
            }
            //std::cout << "qMer: " << qMer << std::endl;
            neighboringResults = this->scoring(qMer, minScore);

            // Ensure the insertion of original q-mers with score of 100 for the original ones
            neighboringResults.neighbors.emplace_back(std::make_pair(qMer,100));
            
            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                #pragma omp critical
                {
                    mutatedQMers += '$';
                    mutatedQMers += neighbor.first; 
                } 
            }
        }

    std::unordered_set<std::string> storedString;
    std::string mutatedQMersUnique;

    for (size_t i = 0; i < mutatedQMers.size(); ++i) {
        if (!isalpha(static_cast<unsigned char>(mutatedQMers[i]))) { // for non valuable char: $
            continue;
        }

        std::string qmer;
        for (size_t j = i; j < mutatedQMers.size() && j - i < qMerSize; ++j) {
            if (isalpha(static_cast<unsigned char>(mutatedQMers[j]))) {
                qmer += mutatedQMers[j];
            }
        }

        if (qmer.size() == qMerSize && storedString.insert(qmer).second) {
            if (!mutatedQMersUnique.empty()) {
                mutatedQMersUnique += '$'; // Insert a '$' separator if not the first substring.
            }
            mutatedQMersUnique += qmer;
        }
        i += qmer.size() - 1;
    }

    return mutatedQMersUnique;
}

// Definition of the member Mutate::mutateSingleSequenceHashes. 
std::vector<uint64_t> Mutate::mutateSingleSequenceHashes(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads, bool minimiser, uint8_t windowSize){

    if(minimiser){
        uint8_t minimiserWindow = windowSize - qMerSize + 1;
        seqan3::debug_stream << "[INFO] Minimiser compuation is applied to the user bins!" << '\n';
        seqan3::debug_stream << "[INFO] Minimiser window size: " << unsigned(minimiserWindow) << '\n';
    }
    std::vector<uint64_t> userBin;
    const unsigned numberOfQMers = (sequence.length() - qMerSize + 1);
    omp_set_num_threads(threads);
    std::mutex hashesMutex;

    #pragma omp parallel for // run each qMer on safe thread
    for(unsigned i = 0; i<numberOfQMers; i++){
            MatrixSearchResults neighboringResults;
            const std::string_view qMer = std::string_view{sequence}.substr(i, qMerSize);
            if (qMer.find('*') != std::string_view::npos) {
                continue;
            }
            neighboringResults = this->scoring(qMer, minScore);

            neighboringResults.neighbors.emplace_back(std::make_pair(qMer,100));
            
            for (const auto& neighbor : neighboringResults.neighbors)
            { 
                /// Convert neighbor string to seqan3::aa27 
                seqan3::aa27_vector qMerHolder{};
                for (const char& c : neighbor.first){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    qMerHolder.push_back(aminoAcid);
                }
                auto hashValue = this->hashQMer(qMerHolder, qMerSize);
                std::lock_guard<std::mutex> lock(hashesMutex);
                userBin.push_back(std::move(hashValue));
            }
        }

    std::sort(userBin.begin(), userBin.end());
        auto uniqueValues = std::unique(userBin.begin(), userBin.end());
        userBin.erase(uniqueValues, userBin.end());
        this->totalNumberOfQMers += userBin.size();

        if(minimiser){

            std::vector<uint64_t> minimiser = computeMinimiserSingleBin(qMerSize, windowSize, userBin);
            return std::move(minimiser);
        }
        else{

            return std::move(userBin);
        }
}
