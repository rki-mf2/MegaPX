
#include <fstream>
#include <vector>
#include <cereal/archives/binary.hpp>
#include <unordered_set>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
 
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
//#include "../megax_util/minimiser_hashes.hpp"
 


void test_cereal() {
    // Serialization
    {
        std::ofstream outputFile("vector_of_vectors.vect", std::ios::binary);
        cereal::BinaryOutputArchive archive(outputFile);

        for (int i = 0; i < 3; ++i) {
            std::vector<uint64_t> innerVector;
            for (int j = 1; j <= 4; ++j) {
                innerVector.push_back(i * 10 + j);
            }

            std::cout << "Values from innerVector" << std::endl;
            for (const auto& value : innerVector){

                std::cout << value << std::endl;
            }
            archive(innerVector);
        }
    }

    // Deserialization
    {
        std::ifstream inputFile("vector_of_vectors.vect", std::ios::binary);
        cereal::BinaryInputArchive archive(inputFile);

        std::vector<std::vector<uint64_t>> loadedData;

        while (true) {
            std::vector<uint64_t> innerVector;
            try {
                archive(innerVector);
                loadedData.push_back(innerVector);
            } catch (const cereal::Exception& e) {
                break; // End of file reached
            }
        }

        std::cout << "Loaded Data:" << std::endl;
        for (const auto& innerVec : loadedData) {
            for (const auto& num : innerVec) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }
    }

}

void test_unordered_set(){

   std::vector<std::vector<uint64_t>> values{{10, 20, 30, 40}, {15, 25, 35, 45}, {5, 15, 25, 35}};
   std::vector<std::set<uint64_t>> hashSet{values.size()};

    // Insert elements from vectors into the unordered_set
    for (size_t i = 0; i < values.size(); ++i) {
        const auto& vec = values[i];
        for (const auto& elem : vec) {
            hashSet[i].insert(elem);
        }
    }

    std::cout << "Elements in hashSet: ";
    for (const auto& vect : hashSet) {

        for (const auto& element : vect) std::cout << element << std::endl;
    }
    std::cout << std::endl;

    std::vector<uint64_t> singleVector = {5, 10, 15, 20, 25};
    /*
    std::cout << "Elements found in hashSet:" << std::endl;
    for (uint64_t element : singleVector) {
        auto it = hashSet.find(element);
        if (it != hashSet.end()) {
            size_t index = std::distance(hashSet.begin(), it);
            std::cout << element << " found at index " << index << std::endl;
        }
    }*/

    int counter = 1; 
    for (auto const& searchValue : singleVector){
        for (auto const& targetVector : hashSet){

            if (auto iter = targetVector.find(searchValue); iter != targetVector.end())
            std::cout << "Found: " << searchValue << " In: " << counter << '\n';
        }
        counter++; 
    }


}

 

void write_db()
{
    using namespace seqan3::literals;
 
    //seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};
     seqan3::sequence_file_output fout{"file.fa"};

    // Define a set to store unique sequences
    std::unordered_set<std::string> unique_sequences;

    using types = seqan3::type_list<std::vector<seqan3::aa27>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    for (int i = 0; i < 5; ++i)
    {
        std::string id = "test_id_" + std::to_string(i);
        seqan3::aa27_vector sequence {"ANML"_aa27};

        // Convert the sequence to a string
        std::string sequence_str;
        for (auto const & base : sequence)
        {
            sequence_str += seqan3::to_char(base);
        }

        // Check if the sequence is already in the set
        if (unique_sequences.find(sequence_str) == unique_sequences.end())
        {
            // If not, add it to the set and write the record to the file
            unique_sequences.insert(sequence_str);
            sequence_record_type record{std::move(sequence), std::move(id)};
            fout.push_back(record);
        }
        // If the sequence is already in the set, skip writing the record to the file
    }
    

}


 

void test_minimiser()
{
    using namespace seqan3::literals;
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};
 
    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
    // auto minimisers = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
     
    std::vector<uint64_t> hashesSet = {254, 841, 3264, 1, 0, 12, 1846, 2, 15, 14};
    //auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);
    auto minimisers = hashesSet | seqan3::views::minimiser(4); // [1, 0, 2]
    auto minimisers_3 = hashesSet | seqan3::views::minimiser(3);
    // results in: [10322096095657499240, 10322096095657499142, 10322096095657499224]
    // representing the k-mers [GTAC, TCGA, GACG]
    seqan3::debug_stream << "Four minimizers: " <<  minimisers << '\n';
    seqan3::debug_stream << "Three minimizers: " <<  minimisers_3 << '\n'; // [254, 1, 0, 2]
    std::string s = "AILM";
    seqan3::aa27_vector aa{};
    for (const char& c : s){

                    seqan3::aa27 aminoAcid = {};
                    aminoAcid.assign_char(c);
                    aa.push_back(aminoAcid);
                }
     
    auto hash  = aa | seqan3::views::kmer_hash(seqan3::ungapped{4});
    seqan3::debug_stream << seqan3::to_rank('A'_aa27) << '\n';
    seqan3::debug_stream << seqan3::to_rank('I'_aa27) << '\n';
    seqan3::debug_stream << seqan3::to_rank('L'_aa27) << '\n';
    seqan3::debug_stream << seqan3::to_rank('M'_aa27) << '\n';
    seqan3::debug_stream << "aa: "<< aa << " " << aa.size() << "(8 * 27^2 + 11 * 27^1 + 12 * 27^0) " <<  hash << '\n';

    //AAHILNMY 
    uint8_t kMerSize = 3;
    uint8_t windowSize = 5;
    uint8_t minimiserWindow = windowSize - kMerSize + 1;
    seqan3::aa27_vector testPeptideSeq{"AAHILNMY"_aa27};
    // [AAH, AHI, HIL, ILN, LNM, NMY]
    // 3: [7,197,5330,6142,8382,9825]
    // 4: [197,5330,143923,165846,226338]
    // 5: [5330,143923,3885933,4477866]

    // (3,5): [7,197,5330,6142] --> 3
    // (4,5): [197,5330,143923,165846] --> 2
    // (5,6): [5330,143923,3885933] --> 2 better to use (5,8) --> 4 for long sequences
    auto hashPeptide  = testPeptideSeq | seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
    seqan3::debug_stream << "Hashes of AAHILNMY: " << unsigned(minimiserWindow) << '\n';
    seqan3::debug_stream << hashPeptide << '\n';

    auto minimisersPeptide = hashPeptide | seqan3::views::minimiser(minimiserWindow); 

    seqan3::debug_stream << minimisersPeptide << '\n';
    std::vector<uint64_t> hashVector = {7,197,5330,6142,8382,9825};
    auto minimisersVect = hashVector | seqan3::views::minimiser(minimiserWindow); 
    //std::vector<uint64_t> resultsVect = computeMinimiser(kMerSize, windowSize, hashVector);
    std::vector<uint64_t> resultsVect;
    for(const auto& hash:minimisersVect){
        resultsVect.emplace_back(hash);
    }
    seqan3::debug_stream << resultsVect << '\n';
    }
