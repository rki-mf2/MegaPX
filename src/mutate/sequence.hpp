
#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <cmath>
#include <map>
#include <algorithm>
#include "mutate.hpp"



//==================================================================
// Classes, Tags, Enums, Structs
//==================================================================


/*
 * @namespace AminoAcidsProcessing
 * @inherits None
 * @brief Main functionalities of Amino Acids processing as namespace
 */
namespace AminoAcidsProcessing
{
    const unsigned short int TOTAL_RESIDUES = 20;

    const char INDEX_TO_AMINO_ACID[TOTAL_RESIDUES] = {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    }; 

    char convertIndexToAminoAcid(unsigned i)
    {
    if (i >= AminoAcidsProcessing::TOTAL_RESIDUES)   return 'X';
    return AminoAcidsProcessing::INDEX_TO_AMINO_ACID[i];}

    // Source:  a4_util.h fu berlin algorithms course 
    const unsigned CHAR_TO_INT[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 0,
    // A, B, C, D, E, F, G, H, I, J, K, L,  M,  N, O
    12, 13, 14, 15, 16, 0, 17, 18, 0, 19, 0, 0, 0, 0, 0, 0,
    // P,  Q,  R,  S,  T, U,  V,  W, X,  Y, Z
    0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 0,
    //   a, b, c, d, e, f, g, h, i, j, k, l,  m, n, o
    12, 13, 14, 15, 16, 0, 17, 18, 0, 19, 0, 0, 0, 0, 0, 0,
    // p,  q,  r,  s,  t, u,  v,  w, x,  y, z
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    unsigned mapCharToIndex(char aminoAcid)
{
    return CHAR_TO_INT[static_cast<unsigned>(aminoAcid)];
}

}


//==================================================================
// Global functions.
//==================================================================

/*
 * @fn qMerizeSequence
 * @brief convert the input sequences to q-mers.
 * @signature std::vector<std::string> qMerizeSequence(const std::string& sequence, int qMerSize);
 * @param sequence the input sequence to q-merize.
 * @param qMerSize integer value corresponds to the length of each q-mer.
 * @return vector of strings filled with q-mers.
 * @throws None.
 */
std::vector<std::string> qMerizeSequence(const std::string& sequence, int qMerSize) {
    
    //std::cout << "[INFO] Generating q-mers...." << std::endl;
    std::vector<std::string> substrings;
    for (size_t i = 0; i + qMerSize <= sequence.size(); ++i) {
        substrings.push_back(sequence.substr(i, qMerSize));
    }

    return substrings;
}


/*
 * @fn mapSequenceToIndices
 * @brief map each char to index integer.
 * @signature std::vector<int> mapSequenceToIndices(const std::string& proteinSequence);
 * @param proteinSequence target protein sequence.
 * @return vector of type int containing the converted char.
 * @throws None.
 */
std::vector<int> mapSequenceToIndices(const std::string& proteinSequence) {
    
    // C sytle for locating memory 
    const int numAminoAcids = 20;
    int aaTOInteger['Z' + 1]; // ASCII characters

    for (int i = 0; i < 'Z' + 1; ++i) {
        aaTOInteger[i] = -1; // Default to unknown amino acid
    }

    // Update the known amino acids
    aaTOInteger['A'] = 0;
    aaTOInteger['C'] = 1;
    aaTOInteger['D'] = 2;
    aaTOInteger['E'] = 3;
    aaTOInteger['F'] = 4;
    aaTOInteger['G'] = 5;
    aaTOInteger['H'] = 6;
    aaTOInteger['I'] = 7;
    aaTOInteger['K'] = 8;
    aaTOInteger['L'] = 9;
    aaTOInteger['M'] = 10;
    aaTOInteger['N'] = 11;
    aaTOInteger['P'] = 12;
    aaTOInteger['Q'] = 13;
    aaTOInteger['R'] = 14;
    aaTOInteger['S'] = 15;
    aaTOInteger['T'] = 16;
    aaTOInteger['V'] = 17;
    aaTOInteger['W'] = 18;
    aaTOInteger['Y'] = 19;

    std::vector<int> convertedSequence;
    convertedSequence.reserve(proteinSequence.size());

    for (char aminoAcid : proteinSequence) {
        int index = aaTOInteger[aminoAcid];
        convertedSequence.push_back(index);
    }

    return convertedSequence;
}

#endif /* SEQUENCE_HPP_ */
