#ifndef MINIMISER_HASHES_HPP
#define MINIMISER_HASHES_HPP

#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

typedef std::vector<uint64_t> hashesType;


/*
* @fn computeMinimiserSingleBin
* @brief function for computing minimisers from input hashes bin.
* @signature hashesType computeMinimiserSingleBin(uint8_t kMerSize, uint8_t windowSize, hashesType hashes);
* @param hashes hashes vector.
* @param kMerSize k-mer size for input query hashing.
* @param windowSize range of window size.
* @throws None.
* @return vector of minimisers in type hashesType.
*/
hashesType computeMinimiserSingleBin(uint8_t kMerSize, uint8_t windowSize, hashesType hashes)
{
    using namespace seqan3::literals;
    
    uint8_t minimiserWindow = windowSize - kMerSize + 1;
    auto minimisersVect = hashes | seqan3::views::minimiser(minimiserWindow); 
    hashesType resultsVect;

    for(const auto& hash:minimisersVect) resultsVect.emplace_back(hash);
    
    return resultsVect;

    }
#endif