#include <iostream>
#include <cstddef> // size_t
#include <cinttypes> // uint64_t
#include <cassert> // assert
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp> 
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/io/sequence_file/all.hpp>


/*
* @fn binarySearch
* @brief function for assigning the specific IBF parameters with building IBF.
* @signature void binarySearch(const std::string& vectorFileName, const std::string& queryFileName, uint8_t kMerSize);
* @param vectorFileName: name and path of the input vector. 
* @param queryFileName: path to query file name.
* @param kMerSize k-mer size for input query hashing.
* @param threads number of input threads.
* @throws None.
* @return None.
*/
void binarySearch(const std::string& vectorFileName, const std::string& queryFileName, uint8_t kMerSize, uint8_t threads) {

    omp_set_num_threads(threads);

    seqan3::debug_stream << "[INFO] Loading input user bins from file..." << '\n';
    std::vector<std::vector<uint64_t>> userBins;

    std::ifstream inputVect(vectorFileName, std::ios::binary); 
	cereal::BinaryInputArchive archiveInput(inputVect); 

    while (true) {
        std::vector<uint64_t> userBin;
        try {
            archiveInput(userBin);
            userBins.push_back(userBin);
        } catch (const cereal::Exception& e) {
            break; // End of file reached
        }
    }

    seqan3::debug_stream << "[INFO] Finished loading user bins." << '\n';
    seqan3::debug_stream << "[INFO] Number of input user bins: " << userBins.size() << '\n';

    seqan3::debug_stream << "[INFO] Loading input sequences...." << '\n';

    std::vector<std::vector<uint64_t>> queries; 

    uint64_t numberOfQuerySequences {0};

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{queryFileName};

    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
    seqan3::debug_stream << "[WARNING] K-mer size should be the same while buildling vector index!" << '\n';

    for (auto & record : fin)
    {
        uint16_t numberOfHits {0};

        std::vector<uint64_t> queryHashes; 
        numberOfQuerySequences += 1 ;

        auto const & seq = record.sequence();

        if(seq.size() < kMerSize)
            continue;

        auto hashValue = (seq | hash_adaptor);

        for (auto const & hash : hashValue)
            queryHashes.emplace_back(hash);

        queries.emplace_back(queryHashes);
    }

    std::vector<std::vector<uint8_t>> countsQueries(queries.size(), std::vector<uint8_t>(userBins.size(), 0));

    assert(threads >= 1);
    omp_set_num_threads(threads);

    seqan3::debug_stream << "[INFO] Number of peptide sequences: " << queries.size() << '\n';
    seqan3::debug_stream << "[INFO] Start querying the de-novo sequences." << '\n';
    seqan3::debug_stream << "[INFO] Number of query threads: "<< threads << '\n';

    int queryCounter = 0;

    #pragma omp parallel for
    for (unsigned queryIdx = 0; queryIdx < queries.size(); ++queryIdx) {

        seqan3::debug_stream << "[INFO] Processing query: " << queryCounter << std::endl;
        std::vector<uint8_t> counts(userBins.size(), 0);

        for (const auto& kmer : queries[queryIdx]) {
            for (unsigned userBinIdx = 0; userBinIdx < userBins.size(); ++userBinIdx) {
                if (std::binary_search(userBins[userBinIdx].begin(), userBins[userBinIdx].end(), kmer)) {
                    counts[userBinIdx]++;
                }
            }
        }
        
        //#pragma omp critical
        //{
        // Store the counts in the correct order
            countsQueries[queryIdx] = counts;
        //}

        queryCounter += 1;
    }

    // Print the counts after processing all queries
    /*
    for (size_t i = 0; i < countsQueries.size(); ++i) {
        for (size_t j = 0; j < countsQueries[i].size(); ++j) {
            std::cout << countsQueries[i][j] << " ";
        }
        std::cout << std::endl; // Move to the next row
    }
    */

    std::filesystem::path filePath(vectorFileName);
    std::filesystem::path queryFile(queryFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / (queryFile.filename().stem().string() + "_results_BS.log");
    
    //std::ofstream logFile(resultsFile);

    //if (!logFile.is_open()) 
     // throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());
    
    seqan3::debug_stream << "[INFO] Start writing results to the output file: " << resultsFile.string() << '\n';

    std::ofstream  outputLog(resultsFile.string() , std::ios::binary); 
	cereal::BinaryOutputArchive archiveOutput(outputLog); 
	archiveOutput(countsQueries);

    /*
    for (const auto& row : resultsVector) {
        for (const auto& element : row) {
            logFile << element << ' ';
        }
        logFile << '\n';
    }*/
    countsQueries.clear();
    
}

