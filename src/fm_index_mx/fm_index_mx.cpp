
#include "fm_index_mx.hpp"
#include "mutate.hpp"


//==================================================================
// Member functions.
//==================================================================


// Definition of the member constructor FM_Index_MX::FM_Index_MX. 
FM_Index_MX::FM_Index_MX(){}


// Definition of the member destructor FM_Index_MX::~FM_Index_MX. 
FM_Index_MX::~FM_Index_MX(){
std::cout<<"[INFO-DEV] Delete constructer of FM_Index_MX functionality....."<<std::endl;
}

struct FMParams
{
    /// FM index params struct
    uint8_t   kMerSize{2};
    uint8_t threads {1};
    uint64_t numberOfReferenceSequences {1u};
    std::string referenceFileName {""};
    std::string OutputFileName {""};
    
};

// Definition of the member function FM_Index_MX::initialiseParamsIndex. 
void FM_Index_MX::initialiseParamsIndex(FMParams& userParams){

    seqan3::debug_stream << "[INFO] Start loading FM index parameters..." << '\n';

    if (userParams.kMerSize < 2){

        this->kMerSize = 2; 
    }
    else {this->kMerSize = userParams.kMerSize;}

    if (userParams.threads < 1){

        this->threads = 1; 
    }
    else {this->threads = userParams.threads;}

    this->OutputFileName = userParams.OutputFileName;
    this->referenceFileName = userParams.referenceFileName;
    this->numberOfReferenceSequences = userParams.numberOfReferenceSequences;
}

// Definition of the member function FM_Index_MX::printInputParams.
void FM_Index_MX::printInputParams(){

    seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
	seqan3::debug_stream << "            Input Parameters for Building FM Index              " << " " << '\n';
	seqan3::debug_stream << "Input reference file               : " << this->referenceFileName << '\n';
    seqan3::debug_stream << "Output FM index file               : " << this->OutputFileName << '\n';
    seqan3::debug_stream << "Number of building threads         : " << unsigned(this->threads) << '\n';
    seqan3::debug_stream << "K-mer size                         : " << this->kMerSize << '\n';
    seqan3::debug_stream << "Number of input reference sequences: " << this->numberOfReferenceSequences << '\n';
    seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
}

// Definition of the member function FM_Index_MX::buildFM. 
void FM_Index_MX::buildFM(const std::string& matrixFilePath, int minScore){

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{this->referenceFileName};

    std::vector<std::string> referenceSequences{};

    Mutate mutate(matrixFilePath);

    for (auto & record : fin)
    {
        auto const & seq = record.sequence();
        std::string sequence;
        for (auto const & c : seq) sequence += seqan3::to_char(c);

       const std::string mutatedSequence = mutate.mutateSingleSequence(this->kMerSize, sequence,  minScore, this->threads);
       referenceSequences.emplace_back(mutatedSequence);
       //seqan3::debug_stream << mutatedSequence << std::endl;
    }

    seqan3::debug_stream << "[INFO] Finished mutating input sequences, total number of sequences: " << referenceSequences.size() << '\n';
    seqan3::debug_stream << "[INFO] Start building FM index....." << '\n';

    seqan3::fm_index index{referenceSequences};

    seqan3::debug_stream << "[INFO] Finished building FM index." << '\n';
    seqan3::debug_stream << "[INFO] Writing FM index to the output file: "<< this->OutputFileName << '\n';

    {
        std::ofstream os{this->OutputFileName, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }

}

// Definition of the member function FM_Index_MX::buildFM. 
void FM_Index_MX::countFM(const std::string& IndexFileName, const std::string& queryFileName, uint8_t kMerSize, uint64_t numberOfReferenceSequences, uint8_t threads){

    std::vector<std::vector<std::string>> queries; 

    seqan3::debug_stream << "[INFO] Start querying the de-novo sequences." << '\n';
    seqan3::debug_stream << "[INFO] Loading input FM index...." << '\n';

    seqan3::fm_index<char, seqan3::text_layout::collection> index;
    {
        std::ifstream is{IndexFileName, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }

    seqan3::debug_stream << "[INFO] Finished loading input FM index!" << '\n';
    seqan3::debug_stream << "[INFO] Loading input sequences...." << '\n';

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{queryFileName};

    for (auto & record : fin)
    {
        auto const & seq = record.sequence();
        //seqan3::debug_stream << seq << std::endl;
        std::string sequence;
        for (auto const & c : seq) sequence += seqan3::to_char(c);
        
        if(sequence.length() < kMerSize)
            continue;
        
        //seqan3::debug_stream << sequence << std::endl;
        const unsigned numberOfKMers = (sequence.length() - kMerSize + 1);
        std::vector<std::string> tempKMers {};
        
        for(unsigned i = 0; i<numberOfKMers; i++){

            const std::string kmer = sequence.substr(i, kMerSize);
            tempKMers.emplace_back(kmer);
        }
        queries.emplace_back(tempKMers);
        
    }

    std::vector<std::vector<uint8_t>> resultsVector(queries.size(), std::vector<uint8_t>(numberOfReferenceSequences, 0));

    seqan3::debug_stream << "[INFO] Finished loading de-novo peptide sequences! " << "(" << queries.size() << ")" << '\n';
    
    assert(threads >= 1);
    omp_set_num_threads(threads);

    int queryCounter = 0;
    std::mutex fm; 
    #pragma omp parallel for
    //#pragma omp parallel for reduction(+:resultsVector[:queries.size()][:numberOfReferenceSequences])
    for (unsigned queryIdx = 0; queryIdx < queries.size(); ++queryIdx) {

        std::vector<uint8_t> counts(numberOfReferenceSequences, 0);
        seqan3::debug_stream << "[INFO] Processing query: " << queryCounter << std::endl;
        for (const auto& kmerSet : queries[queryIdx]) {
            // seqan3::debug_stream << kmerSet << std::endl;
            seqan3::algorithm_result_generator_range results = search(kmerSet, index);

            for (auto && result : results){
                //seqan3::debug_stream << result << std::endl;
                uint64_t refernceID = result.reference_id();
                //seqan3::debug_stream << refernceID << std::endl;
                if(refernceID >= 0){
                    #pragma omp critical
                    {
                    // std::lock_guard<std::mutex> lock(fm);
                    resultsVector[queryIdx][refernceID] += 1;
                    // std::cout << resultsVector[queryIdx][refernceID] << std::endl;
                    }       
                }
                
            }
        }

        queryCounter += 1;
    }

    seqan3::debug_stream << "[INFO] Finished processing all queries!" << '\n';

    std::filesystem::path filePath(IndexFileName);
    std::filesystem::path queryFile(queryFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / (queryFile.filename().stem().string() + "_results_FM_index.log");
    //std::ofstream logFile(resultsFile);

    //if (!logFile.is_open()) 
    //  throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());

    seqan3::debug_stream << "[INFO] Start writing results to the output file: " << resultsFile.string() << '\n';

    std::ofstream  outputLog(resultsFile.string() , std::ios::binary); 
	cereal::BinaryOutputArchive archiveOutput(outputLog); 
	archiveOutput(resultsVector);

    /*
    for (size_t i = 0; i < resultsVector.size(); ++i) {
        for (size_t j = 0; j < resultsVector[i].size(); ++j) {
            std::cout << unsigned(resultsVector[i][j]) << " ";
        }
        std::cout << std::endl; // Move to the next row
    }*/

    /*for (const auto& row : resultsVector) {
        for (const auto& element : row) {
            logFile << element << ',';
        }
        logFile << '\n';

    }
    
    logFile.close();*/

    resultsVector.clear();

}