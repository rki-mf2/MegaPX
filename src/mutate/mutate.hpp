
#ifndef MUTATE_HPP_
#define MUTATE_HPP_

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <functional>
#include <map>
#include <unordered_map>
#include <memory>
#include <unordered_set>
#include <cctype>
#include <filesystem>

#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

//==================================================================
// Forwards
//==================================================================

typedef std::vector<std::pair<std::string, std::string>> mutatedDB;
typedef std::vector<std::vector<uint64_t>> QMersHashes;
struct IBFParams;
//==================================================================
// Classes, Tags, Enums, Structs
//==================================================================


/*
 * @class PairHash
 * @inherits None
 * @brief Struct to provide a custom hash function for std::pair (C++11)
 */
struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& stdPairCPlusPlus_11) const {

        auto h1 = std::hash<T1>{}(stdPairCPlusPlus_11.first); 
        auto h2 = std::hash<T2>{}(stdPairCPlusPlus_11.second); 
        return h1 ^ (h2 << 1); // XOR
    }
};

/*
 * @struct MatrixSearchResults
 * @inherits None
 * @brief struct contains the target qMer with the neighbors 
 */
struct MatrixSearchResults
{
  std::string qMer;
  std::vector< std::pair <std::string, int> > neighbors;
};

/*
 * @class Mutate
 * @inherits None
 * @brief Process and Mutate sequences according to given substitution matrix 
 */
class Mutate
{
	public:

        /*
         * @fn Mutate::Mutate
         * @brief Constructor.
         */
		Mutate(const std::string& matrixFilePath);

        /*
         * @fn Mutate::~Mutate
         * @brief Destructor.
         */
		~Mutate();

        /*
         * @fn Mutate::loadMatrix
         * @brief loads the scoring matrix from path.
		 * @signature void loadMatrix(const std::string& matrixFilePath);
		 * @param matrixFilePath the path, where the input matrix is stored.
		 * @throws std::runtime_error if the matrix not found. 
		 * @return None
         */
		void loadMatrix(const std::string& matrixFilePath);

		/*
         * @fn Mutate::processMatrix
         * @brief process the loaded matrix with sorting.
		 * @signature void processMatrix(const std::string& matrixFilePath);
		 * @param matrixFilePath the path, where the input matrix is stored.
		 * @throws None. 
		 * @return None
         */
		void processMatrix(const std::string& matrixFilePath);

        /*
         * @fn Mutate::scoreSetter
         * @brief set the input score to the pair of amino acids.
		 * @signature void scoreSetter(char& aa_1, char& aa_2, int score);
		 * @param aa_1 first amino acid.
		 * @param aa_2 second amino acid.
		 * @param score score to be set to the pair of first and second amino acid.
		 * @throws None. 
		 * @return None
         */
		//2void scoreSetter(char& aa_1, char& aa_2, int score);
		void scoreSetter(int& aa_1, int& aa_2, int score);
        /*
         * @fn Mutate::scoreGetter
         * @brief get the input score of pair of amino acids.
		 * @signature int scoreGetter(const char& aa_1, const char& aa_2) const;
		 * @param aa_1 first amino acid.
		 * @param aa_2 second amino acid.
		 * @throws None. 
		 * @return integer value corresponding the score of both amino acids. 
         */
    	//4int scoreGetter(const char& aa_1, const char& aa_2) const;
		int scoreGetter(int aa_1, int aa_2);

		/*
         * @fn Mutate::calculateQmerScore
         * @brief compute the score between q-mer and neighbor.
		 * @signature calculateQmerScore(const std::string &qMer, const std::string &neighbor);
		 * @param qMer query q-mer string.
		 * @param neighbor the q-mer from the matrix .
		 * @throws None. 
		 * @return integer value corresponding the score of comparing both. 
         */
		int calculateQmerScore(const std::string &qMer, const std::string &neighbor);

		/*
         * @fn Mutate::mutateSequencesRec
         * @brief mutate input sequences given into FASTA file.
		 * @signature QMersHashes mutateSequencesRecmake(int qMerSize, std::string& sequenceFile, int minScore, int threads);
		 * @param qMerSize integer value of the user-defined q-mer size [default 1].
		 * @param sequenceFilePath path to the input FASTA file.
		 * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @param threads number of input threads [default: 1].
		 * @param outputDir path to the output directory.
		 * @param minimiser bool value to use minimiser. 
		 * @param windowSize used window size for the minimiser computation.
		 * @throws std::runtime_error if the sequences file not found.
		 * @return None. 
         */
		void mutateSequencesRec(uint8_t qMerSize, std::string& sequenceFilePath, int minScore, uint8_t threads, std::string& outputDir, bool minimiser, uint8_t windowSize);

		/*
         * @fn Mutate::scoring
         * @brief score different qMers against the input matrix.
		 * @signature MatrixSearchResults scoring(const std::string& qMer, int minScore);
		 * @param qMer database sequence query.
		 * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @return MatrixSearchResults avector contains pair of the results and related scores. 
         */
		MatrixSearchResults scoring(const std::string_view& qMer, int minScore);

		/*
         * @fn Mutate::writeDBToFasta
         * @brief write mutated sequences to output file.
		 * @signature writeDBToFasta(uint8_t qMerSize, std::string& sequenceFilePath, int minScore, uint8_t threads, std::string& outputDir);
		 * @param qMerSize integer value of the user-defined q-mer size [default 1].
		 * @param sequenceFilePath path to the input FASTA file.
		 * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @param threads number of input threads [default: 1].
		 * @param outputDir path to the output directory.
		 * @throws std::runtime_error if the sequences file not found.
         */
		void writeDBToFasta(uint8_t qMerSize, std::string& sequenceFilePath, int minScore, uint8_t threads, std::string& outputDir);

		/*
         * @fn Mutate::mutateToSequenceLength
         * @brief write top N mutated sequences to output file.
		 * @signature mutateToSequenceLength(uint8_t numberOfTopSequences, std::string& sequenceFilePath, uint8_t threads, std::string& outputDir);
		 * @param numberOfTopSequences integer value of the number of top generated sequences pro protein[default 1].
		 * @param sequenceFilePath path to the input FASTA file.
		 * @param threads number of input threads [default: 1].
		 * @param outputDir path to the output directory.
		 * @param minScore mutation score. 
		 * @throws std::runtime_error if the sequences file not found.
         */
		void mutateToSequenceLength(uint8_t numberOfTopSequences, std::string& sequenceFilePath, std::string& outputDir, int minScore);

		/*
         * @fn Mutate::mutateSingleSequence
         * @brief mutate one single sequence.
		 * @signature void mutateSingleSequence(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads);
		 * @param sequence target sequence to mutate.
		 * @param qMerSize integer value of the user-defined q-mer size [default 1].
		 * @param sequenceFilePath path to the input FASTA file.
		 * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @param threads number of input threads [default: 1].
		 * @throws None.
		 * @return std::string mutated q-mers as one sequence with $ between the q-mers. 
         */
		std::string mutateSingleSequence(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads);

		/*
         * @fn Mutate::mutateSingleSequenceHashes
         * @brief mutate one single sequence.
		 * @signature void mutateSingleSequenceHashes(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads);
		 * @param sequence target sequence to mutate.
		 * @param qMerSize integer value of the user-defined q-mer size [default 1].
		 * @param sequenceFilePath path to the input FASTA file.
		 * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @param threads number of input threads [default: 1].
		 * @param minimiser bool value to enable minimiser. 
		 * @param windowSize used window size for the minimiser computation.
		 * @throws None.
		 * @return std::vector<uint64_t> mutated q-mers as one vector. 
         */
		std::vector<uint64_t> mutateSingleSequenceHashes(uint8_t qMerSize, const std::string& sequence,  int minScore, uint8_t threads, bool minimiser, uint8_t windowSize);

		/*
		* @fn printSequenceStatstics
		* @brief print sequence statistics.
		* @signature void printSequenceStatstics(const std::string& inputFastaFile, uint8_t kMerSize, std::string& outputDir);
		* @param inputFastaFile target protein sequence.
		* @param kMerSize k-mer size. 
		* @param outputDir path to output directory for bin mapping files.
		* @return None.
		* @throws std::runtime_error if the sequences file not found.
		*/
		void printSequenceStatstics(const std::string& inputFastaFile, uint8_t kMerSize, std::string& outputDir);

		/*
		* @fn printMutationStatstics
		* @brief print mutation statistics.
		* @signature void printMutationStatstics(const std::string& inputFastaFile, uint8_t kMerSize);
		* @param inputFastaFile target protein sequence.
		* @param kMerSize k-mer size. 
		* @param minScore minimum mutation score.
		* @param minimiser True/False value
		* @return None.
		* @throws std::runtime_error if the sequences file not found.
		*/
		void printMutationStatstics(const std::string& inputFastaFile, uint8_t kMerSize, int minScore, bool minimiser);

		/*
         * @var numberOfInputSequences
         * @brief number of input sequences.
		 * @signature uint64_t numberOfInputSequences {0};
         */
		uint64_t numberOfInputSequences 
		{0};

		/*
         * @var totalNumberOfOriginalQMers
         * @brief total number of original q_mers.
		 * @signature uint64_t totalNumberOfOriginalQMers {0};
         */
		uint64_t totalNumberOfOriginalQMers
		{0};

		/*
         * @var totalNumberOfQMers
         * @brief total number of q_mers.
		 * @signature uint64_t totalNumberOfQMers {0};
         */
		uint64_t totalNumberOfQMers
		{0};


		/*
         * @fn hashQMer
         * @brief hash provided input q_mer string.
		 * @signature uint64_t inline hashQMer(seqan3::aa27_vector kmer, u_int8_t qmerSize);
		 * @param QMer q_mer string of type seqan3::aa27_vector.
		 * @param qmerSize the size of used qmer for hashing adaptor.
		 * @throws None
		 * @return uint64_t hash value. 
         */
		uint64_t inline hashQMer(seqan3::aa27_vector QMer, u_int8_t qmerSize){

			auto hashAdaptor = seqan3::views::kmer_hash(seqan3::ungapped{qmerSize});
			auto hash = (QMer | hashAdaptor);

			return hash[0];
		}

		/*
         * @fn writeVector
         * @brief function to write vector of hashes to file.
		 * @signature void inline writeVector(QMersHashes hashesVector, const std::string& outputPath);
		 * @param hashesVector target vector.
		 * @param outputPath output path.
		 * @throws None.
		 * @return None. 
         */
		void inline writeVector(QMersHashes hashesVector, const std::string& outputPath){

			std::ofstream  os(outputPath, std::ios::binary); 
			cereal::BinaryOutputArchive archive(os); 
			archive(hashesVector);

		}
		
	// Protected member variables and methods:
	protected:

		/*
         * @var substitutionMatrix
         * @brief substitution matrix loaded as pair of char and mapped with the corresponding score.
		 * @signature std::map<std::pair<char, char>, int> substitutionMatrix;
         */
		std::unordered_map<std::pair<int, int>, int, PairHash> substitutionMatrix;

		/*
         * @var sortedSubstitutionMatrix
         * @brief sorted substitution matrix.
		 * @signature std::vector <std::vector <std::pair<std::string, int>>> sortedSubstitutionMatrix;
         */
		std::vector <std::vector <std::pair<std::string, int>>> sortedSubstitutionMatrix;

	// Private member variables and methods:
	private:

		/*
         * @fn compareAminoAcids
         * @brief function for comparing the scores odf the internal amino acids in the matrix to sort them.
		 * @signature compareAminoAcids(const std::pair<std::string, int>& aa_1, const std::pair<std::string, int>& aa_2);
		 * @param aa_1 first amino acid pair of type std::pair<std::string, int>.
		 * @param aa_2 second amino acid pair of type std::pair<std::string, int>.
		 * @throws None.
		 * @return bool value if aa_1.score < aa_2.score. 
         */
		bool compareAminoAcids(const std::pair<std::string, int>& aa_1, const std::pair<std::string, int>& aa_2);

		/*
         * @var mutatedSequences
         * @brief mutateDB vector contains original and mutated sequences (might be removed due to memory issues!).
		 * @signature mutatedDB mutatedSequences;
         */
		mutatedDB mutatedSequences;

};


#endif /* MUTATE_HPP_ */
