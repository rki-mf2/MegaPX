#ifndef HIBF_HPP_
#define HIBF_HPP_


#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp> 
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "mutate.hpp"
#include <omp.h>
#include <filesystem>
#include <bitset>

//==================================================================
// Forwards
//==================================================================

struct HIBFParams;
struct IBFParams;

//==================================================================
// Classes, Tags, Enums, Structs
//==================================================================

/*
 * @class HIBF
 * @inherits None
 * @brief Process and build sequences according into interleaved Bloom filters
 */
class HIBF
{
	public:

        /*
         * @fn HIBF::HIBF
         * @brief Constructor.
         */
		HIBF();

        /*
         * @fn HIBF::~HIBF
         * @brief Destructor.
         */
		~HIBF();
		
        void run_test();
        
        /*
         * @var bins
         * @brief Number of input bins [default: 1].
		 * @signature uint32_t bins{1};
         */
        uint64_t  bins{1};

        /*
         * @var maxIBFSize
         * @brief maximal IBF sizess [default: 0].
		 * @signature double maxIBFSize {0u};
         */
        double maxIBFSize {0u};

        /*
         * @var kMerSize 
         * @brief K-mer size [default: 2].
		 * @signature uint32_t kMerSize{2};
         */
        uint8_t   kMerSize {2}; /// only for IBF 

        /*
         * @var binSize
         * @brief Bin size [default: 0].
		 * @signature uint32_t kMerSize{2};
         */
        size_t    binSize { 0u }; /// only for IBF 

        /*
         * @var falsePositiveRate
         * @brief false positive rate [unchangeable value].
		 * @signature double   falsePositiveRate {0.01};
         */
        double   falsePositiveRate {0.01}; /// only for IBF 

        /*
         * @var sequenceLength
         * @brief length of longest sequence.
		 * @signature uint64_t    sequenceLength; 
         */
        uint64_t    sequenceLength;  /// only for IBF 

        /*
         * @var numberOfHashFunctions
         * @brief Number of hash functions [default: 3].
		 * @signature uint32_t numberOfHashFunctions{3};
         */
        uint8_t   numberOfHashFunctions {2};

        /*
         * @var IBFFileName
         * @brief HIBF output file [default: ""].
		 * @signature std::string IBFFileName {""};
         */
        std::string IBFFileName {""};

        /*
         * @var threads
         * @brief Number of input threads [default: 1].
		 * @signature uint8_t threads{1};
         */
        uint8_t threads {1}; 

        /*
         * @var sketchBits
         * @brief Number of sketch bits [default: 12].
		 * @signature uint8_t sketchBits{12};
         */
        uint8_t sketchBits {12};

        /*
         * @var alpha
         * @brief Value of alpha [default: 1.2].
		 * @signature double alpha{1.2};
         */
        double alpha {1.2};

        /*
         * @var maxUserBins
         * @brief Maximum number of user bins to insert in each HIBF [default: 1].
		 * @signature uint64_t maxUserBins {1u};
         */
        uint64_t maxUserBins {1u};

        /*
         * @var maximumFalsePositiveRate
         * @brief Maximum false positive rate [default: 0.01].
		 * @signature double maximumFalsePositiveRate {0.01};
         */
        double maximumFalsePositiveRate {0.01};

        /*
         * @var maxRearrangementRatio
         * @brief Maximum rearrangement ratio [default: 0.5].
		 * @signature double maxRearrangementRatio {0.5};
         */
        double maxRearrangementRatio {0.5};

        /*
         * @var disableEstimationRatio
         * @brief Disable estimation ratio [default: false].
		 * @signature bool disableEstimationRatio = false;
         */
        bool disableEstimationRatio = false;

        /*
         * @var disableRearrangement
         * @brief Disable rearrangement [default: false].
		 * @signature bool disableRearrangement = false;
         */
        bool disableRearrangement = false;

        /*
         * @var vectFileName
         * @brief Path to input vector file.
		 * @signature std::string vectFileName {""};
         */
        std::string vectFileName {""};

        /*
         * @var filterFileName
         * @brief Path to output hibf file.
		 * @signature std::string filterFileName {""};
         */
        std::string filterFileName {""};

        /*
         * @fn buildIBF
         * @brief function for assigning the specific IBF parameters with building IBF.
		 * @signature void buildIBF(uint8_t numberOfHashFunctions, std::string vectFileName);
		 * @param numberOfHashFunctions (uint8_t).
         * @param vectorFileName: name and path of the input vector. 
         * @param outputDir: path to output dir.
		 * @throws None.
		 * @return None.
         */
        void buildIBF(uint8_t numberOfHashFunctions, std::string& vectorFileName, std::string& outputDir);

         /*
         * @fn initialiseParamsHIBF
         * @brief function for assigning the specific HIBF parameters.
		 * @signature void initialiseParamsHIBF(HIBFParams& userParams);
		 * @param userParams Struct contains the user input arguments.
         * @param inputFastaFilePath path of the input fasta file. 
		 * @throws std::runtime_error.
		 * @return None.
         */
        void initialiseParamsHIBF(HIBFParams& userParams);

        /*
         * @fn printInputParams
         * @brief function for printing HIBF params.
		 * @signature void printInputParams();
		 * @param None.
		 * @throws None.
		 * @return None.
         */
        void printInputParams();

        /*
         * @fn buildHIBF
         * @brief Build hibf function.
		 * @signature void buildHIBF();
		 * @param None.
		 * @throws None.
		 * @return None.
         */
        void buildHIBF();

        /*
         * @fn computeBinSize
         * @brief function for computing the bin size.
		 * @signature size_t computeBinSize(uint64_t numberOfKmers);
		 * @param uint64_t numberOfKmers number of input k-mers in the longest sequence.
		 * @throws None.
		 * @return bin size.
         */
        uint64_t computeBinSizeIBF(uint64_t numberOfKmers); /// only for IBF 

        /*
         * @fn countIBF
         * @brief function for counting k-mers into IBF.
		 * @signature void countIBF(const std::string& IBFFileName, const std::string& queryFileName);
		 * @param IBFFileName path to the input IBF file name.
         * @param queryFileName path to the query file, which contains the target de-novo sequences.
         * @param kMerSize size of input k-mer for counting. 
		 * @throws None.
		 * @return atm None.
         */
        void countIBF(const std::string& IBFFileName, const std::string& queryFileName, uint8_t kMerSize);

         /*
         * @fn buildHIBFRef
         * @brief Build HIBF directly from input reference sequences.
		 * @signature void buildHIBFRef(const std::string& inputFastaFile, uint8_t kMerSize, uint8_t threads);
		 * @param inputFastaFile input fasta file.
         * @param kMerSize k-mer size for construction. 
		 * @throws None.
		 * @return None.
         */
        void buildHIBFRef(const std::string& inputFastaFile, uint8_t kMerSize);

        /*
         * @fn countHIBF
         * @brief function for counting k-mers into HIBF.
		 * @signature void countHIBF(const std::string& IBFFileName, const std::string& queryFileName);
		 * @param HIBFFileName file contain the paths of the splitted HIBFs.
         * @param queryFileName path to the query file, which contains the target de-novo sequences.
         * @param kMerSize size of input k-mer for counting. 
         * @param hibfThreshold The membership_for function takes the query and a threshold. Here, a threshold of two means that
         * at least (>=) 2 values of the query must be found within a user bin to be a hit.
		 * @throws None.
		 * @return atm None.
         */
        void countHIBF(const std::string& HIBFFileName, const std::string& queryFileName, uint8_t kMerSize, uint8_t hibfThreshold);

        /*
         * @fn buildSingleIBF
         * @brief function for building one IBF.
		 * @signature void buildSingleIBF(std::vector<std::vector<uint64_t>> userBins, uint64_t IBFSize);
		 * @param userBins vector of user bins hashes.
         * @param numberOfHashFunctionsIn number of input hash functions. 
		 * @throws None.
		 * @return IBF.
         */
        seqan::hibf::interleaved_bloom_filter buildSingleIBF(std::vector<std::vector<uint64_t>> userBins, uint64_t numberOfHashFunctionsIn);

        /*
         * @fn multiIndexing
         * @brief function for building many IBF in real time.
		 * @signature void multiIndexing(const std::string& matrixFilePath, const std::string& inputFastaFile, const std::string& queryFileName, uint8_t kMerSize, int minScore, uint8_t threads, bool minimiser, uint8_t windowSize, uint64_t numberOfHashFunctionsIn, int splitSize);
		 * @param matrixFilePath path to blosum/pam matrix. 
         * @param inputFastaFile input reference file.
         * @param queryFileName input query file.
         * @param kMerSize k_mer size.
         * @param minScore minimum mutation score.
         * @param threads number of building and mutating threads.
         * @param minimiser bool value for enabling minimiser.
         * @param windowSize minimiser assigned window size.
         * @param numberOfHashFunctionsIn number of input hash functions. 
         * @param splitSize number of input references in each IBF. 
         * @param outputFile output file for sroting the results. 
         * @param blackList input file contains the blacklist references.
         * @param threshold assignment threshold.
		 * @throws std::runtime_error.
		 * @return None.
         */
        void multiIndexing(const std::string& matrixFilePath, const std::string& blackList, const std::string& inputFastaFile, const std::string& queryFileName, uint8_t kMerSize, int minScore, uint8_t threads, bool minimiser, uint8_t windowSize, uint64_t numberOfHashFunctionsIn, uint64_t splitSize, const std::string& outputFile, double threshold);
        
	// Protected member variables and methods:
	protected:

	// Private member variables and methods:
	private:

		/*
        * @fn countSequencesRef
        * @brief counting the number of sequences in the input fasta file.
        * @signature countSequencesRef(const std::string& inputFastaFilePath);
        * @param inputFastaFilePath path to the input fasta file.
        * @return uint32_t type as number of sequences.
        * @throws std::runtime_error if the file is not openable.
        */
        uint32_t countSequencesRef(const std::string& inputFastaFilePath); /// only for IBF 

        /*
        * @fn optimalNumberOfBins
        * @brief computes optimal number of bins for HIBF (multiple 64 bits).
        * @signature inline uint64_t optimalNumberOfBins().
        * @param class intern variables.
        * @return uint64_t computed optimal number of bins.
        * @throws None.
        */
        inline uint64_t optimalNumberOfBins() /// only for IBF 
        {
            return std::ceil( this->bins / 64.0 ) * 64;
        }

        /*
        * @fn loadBlackList
        * @brief get blacklist outliers from the txt file.
        * @signature inline void loadBlackList().
        * @param blackList path to the blacklist file.
        * @param outliers vector of string to store outliers in.
        * @return None.
        * @throws None.
        */
        inline void loadBlackList(const std::string& blackList, std::vector<std::string>& outliers) 
        {
            std::ifstream inputFile(blackList);
            if (!inputFile.is_open()) {
                std::cerr << "[Error] Unable to open blacklist file: " << blackList << std::endl;
                return;
            }

            std::string line;
            while (std::getline(inputFile, line)) {
                outliers.push_back(line);
            }

            inputFile.close();

            std::cout << "[INFO] Total number of outliers: "<< outliers.size() << std::endl;
        }
};


#endif /* HIBF_HPP_ */