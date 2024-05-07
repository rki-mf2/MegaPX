#ifndef FM_INDEX_MX_HPP_
#define FM_INDEX_MX_HPP_


#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp> 
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

#include <seqan3/core/debug_stream.hpp> // pretty printing
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search_result.hpp>
#include <seqan3/search/search.hpp>

#include <seqan3/io/sequence_file/all.hpp>

#include <filesystem>

//==================================================================
// Forwards
//==================================================================

struct FMParams;

//==================================================================
// Classes, Tags, Enums, Structs
//==================================================================

/*
 * @class FM_Index_MX
 * @inherits None
 * @brief Build and count into FM Index
 */
class FM_Index_MX
{
	public:

        /*
         * @fn FM_Index_MX::FM_Index_MX
         * @brief Constructor.
         */
		FM_Index_MX();

        /*
         * @fn FM_Index_MX::~FM_Index_MX
         * @brief Destructor.
         */
		~FM_Index_MX();
		

        /*
         * @var kMerSize 
         * @brief K-mer size [default: 2].
		 * @signature uint32_t kMerSize{2};
         */
        uint8_t   kMerSize {2}; 

        /*
         * @var numberOfReferenceSequences 
         * @brief Number of input reference sequences [default: 1].
		 * @signature uint64_t numberOfReferenceSequences {1u}
         */
        uint64_t numberOfReferenceSequences {1u};


        /*
         * @var OutputFileName
         * @brief Index output file [default: ""].
		 * @signature std::string OutputFileName {""};
         */
        std::string OutputFileName {""};

        /*
         * @var referenceFileName
         * @brief Path to reference input file.
		 * @signature std::string referenceFileName {""};
         */
        std::string referenceFileName {""};

        /*
         * @var queryFileName
         * @brief Path to query input file.
		 * @signature std::string queryFileName {""};
         */
        std::string queryFileName {""};

        /*
         * @var threads
         * @brief Number of input threads [default: 1].
		 * @signature uint8_t threads{1};
         */
        uint8_t threads {1};

        /*
         * @fn buildFM
         * @brief function for building the fm index.
		 * @signature void buildFM(const std::string& matrixFilePath);
		 * @param matrixFilePath path to mutation matrix.
         * @param minScore minimum score of the matched q-mer, to assign the neighbor as mutation [default: 0].
		 * @throws None.
		 * @return None.
         */
        void buildFM(const std::string& matrixFilePath, int minScore);

         /*
         * @fn initialiseParamsIndex
         * @brief function for assigning the specific FM index parameters.
		 * @signature void initialiseParamsIndex(FMParams& userParams);
		 * @param userParams Struct contains the user input arguments.
		 * @throws std::runtime_error.
		 * @return None.
         */
        void initialiseParamsIndex(FMParams& userParams);

         /*
         * @fn printInputParams
         * @brief function for printing FM index params.
		 * @signature void printInputParams();
		 * @param None.
		 * @throws None.
		 * @return None.
         */
        void printInputParams();

        /*
         * @fn countFM
         * @brief function for counting k-mers into FM index.
		 * @signature void countFM(const std::string& IndexFileName, const std::string& queryFileName, uint8_t kMerSize);
		 * @param IndexFileName path to the input index file name.
         * @param queryFileName path to the query file, which contains the target de-novo sequences.
         * @param kMerSize size of input k-mer for counting. 
         * @param numberOfReferenceSequences number of input sequences for vectorising the output. 
         * @param threads number of input threads.
		 * @throws std::runtime_error.
		 * @return atm None.
         */
        void countFM(const std::string& IndexFileName, const std::string& queryFileName, uint8_t kMerSize, uint64_t numberOfReferenceSequences, uint8_t threads);

	// Protected member variables and methods:
	protected:

	// Private member variables and methods:
	private:

};


#endif /* FM_INDEX_MX */