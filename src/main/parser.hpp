
#ifndef PARSER_HPP_
#define PARSER_HPP_

// MegaX headers
#include "mutate.hpp"
#include "hibf_util.hpp"
#include "../megax_util/binary_search.hpp"
#include "../megax_util/evaluation.hpp"
//#include "../megax_util/taxor_tax_profile.hxx"
#include "../megax_util/tax_classification.hxx"
#include "../megax_util/peptide_simulation.hxx"
#include "fm_index_mx.hpp"

// Testing
//#include "test_functions.hxx"

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>


using namespace seqan3;


//==================================================================
// Classes, Tags, Enums, Structs
//==================================================================

/*
 * @struct CmdArguments
 * @inherits None
 * @brief Main user arguments
 */
struct CmdArguments
{
	/// Mutation params struct
	uint8_t  numberOfThreads{1u};
	int      minScore{0};
	uint8_t  qMerSize{2};
	uint8_t windowSize{0u};
	uint8_t numberOfTopSequences{1u};

	bool minimiser = false; 
	uint64_t numberOfReferenceSequences {1u};

	std::string outputDir 
	{""};

    /// HIBF params struct
    uint8_t   numberOfHashFunctions {2};
	uint8_t hibfThreshold {1u};
	uint8_t threads {1u}; 
	uint8_t sketchBits {12};
	uint64_t maxUserBins {1u};
	double alpha {1.2};
	double maximumFalsePositiveRate {0.01};
	double maxRearrangementRatio {0.5};
	bool disableEstimationRatio = false;
	bool disableRearrangement = false;
	std::string vectFileName
	{""};
    std::string filterFileName
    {""};

	std::string indexFileName
    {""};

	std::string queryFileName
	{""};

    /// General params struct
	std::string sequenceFilePath
    {""};
	std::string matrixFilePath
    {""}; 
	std::string mode
    {""};
	std::string countingResults
    {""};
	std::string referenceMapping
    {""}; 
	std::string evaluationResults
    {""};
	std::string taxClassificationFile
	{""};
	std::string blackListFile
	{""};

	double threshold {0.75};

	
};

/*
 * @struct HIBFParams
 * @inherits None
 * @brief Main hibf arguments
 */
struct HIBFParams
{
    /// HIBF params struct
    uint8_t   numberOfHashFunctions {2};
	//uint8_t hibfThreshold {1u};
	uint8_t threads {1u}; 
	uint8_t sketchBits {12};
	uint64_t maxUserBins {1u};
	double alpha {1.2};
	double maximumFalsePositiveRate {0.01};
	double maxRearrangementRatio {0.5};
	bool disableEstimationRatio = false;
	bool disableRearrangement = false;
	std::string vectFileName {""};
    std::string filterFileName {""};
	
};

/*
 * @struct FMParams
 * @inherits None
 * @brief Main fm index arguments
 */
struct FMParams
{
    /// FM index params struct
    uint8_t   kMerSize{2};
    uint8_t threads {1};
	uint64_t numberOfReferenceSequences {1u};
    std::string referenceFileName {""};
	std::string OutputFileName {""};
};

//==================================================================
// Global functions.
//==================================================================


/*
 * @fn initializeMainArgumentParser
 * @brief get user defined arguments.
 * @signature void initializeMainArgumentParser(argument_parser &parser, CmdArguments &args);
 * @param MegaXParser seqan3 parser argument.
 * @param args user parameters struct.
 * @return None.
 * @throws None.
 */
void initializeMainArgumentParser(argument_parser &MegaXParser, CmdArguments &args)
{
	MegaXParser.info.author = "Ahmad Lutfi";
	MegaXParser.info.short_description = "MegaX: protein and sequence level taxonomic classification using mutated databases.";
	MegaXParser.info.version = "0.0.0";
	MegaXParser.info.date = "2024";
	MegaXParser.info.email = "ahmad.lutfi@fu-berlin.de";

	std::vector<std::string> description
	{ "MegaX builds and counts mutations from and in datasets with the classification of unknown samples." };
	MegaXParser.info.description = description;

	std::vector<std::string> synopsis
	{"[build_vect, build_ibf, count_ibf, build_hibf, count_hibf, binary_search, build_fm, count_fm, counting, stat, mutate_stat, evaluate, hibf_ref, classification, ibf_stat, write_db, mutate_seq_len, simulate_peptides, profile, test, multi_indexing] [OPTIONS]"};
	MegaXParser.info.synopsis = synopsis;


	MegaXParser.add_positional_option(args.mode, "Modus to run MegaX : ", value_list_validator
	{ "build_vect", "build_ibf", "count_ibf", "build_hibf", "count_hibf", "binary_search", "build_fm", "count_fm", "counting", "stat", "mutate_stat", "evaluate", "hibf_ref", "classification", "ibf_stat", "write_db", "mutate_seq_len", "simulate_peptides", "profile", "test", "multi_indexing"});

}


/*
 * @fn initializeArgumentParser
 * @brief get user defined arguments.
 * @signature void initializeArgumentParser(argument_parser &MegaXParser, CmdArguments &args);
 * @param MegaXParser seqan3 parser argument.
 * @param args user parameters struct.
 * @return None.
 * @throws None.
 */
void initializeArgumentParser(argument_parser &MegaXParser, CmdArguments &args)
{
	MegaXParser.info.author = "Ahmad Lutfi";
	MegaXParser.info.short_description = "MegaX builds and counts mutations from and in datasets with the classification of unknown samples.";
	MegaXParser.info.version = "0.0.0";
	MegaXParser.info.date = "2024";
	MegaXParser.info.email = "ahmad.lutfi@fu-berlin.de";

    /// Mutations params
	MegaXParser.add_option(args.numberOfThreads, 't', "Threads", "Number of input threads.");
	MegaXParser.add_option(args.minScore, 's', "Score", "MinimumScore for each mutation.");
	MegaXParser.add_option(args.qMerSize, 'q', "QMerSize", "QMerSize for each k-mer.");
	MegaXParser.add_option(args.sequenceFilePath, 'i', "InputFasta", "Input fasta file");
	MegaXParser.add_option(args.matrixFilePath, 'm', "Matrix", "Path to input matrix");
	MegaXParser.add_option(args.outputDir, 'o', "outputDir", "Output directory for writing logs");
	MegaXParser.add_option(args.minimiser, 'Z', "minimiser", "Use minimiser in one level.");
	MegaXParser.add_option(args.windowSize, 'w', "windowSize", "Window size for minimiser computation, valid only with -Z true.");
	MegaXParser.add_option(args.numberOfTopSequences, 'N', "numberOfTopSequences", "Number of top scored written sequences of each protein.");

    /// HIBF params
    MegaXParser.add_option(args.alpha, 'p', "alpha", "Alpha value.");
    MegaXParser.add_option(args.sketchBits, 'S', "sketchBits", "HyperLogLog sketch bits.");
	MegaXParser.add_option(args.maximumFalsePositiveRate, 'r', "maximumFalsePositiveRate", "Maximum false positive rate.");
	MegaXParser.add_option(args.maxRearrangementRatio, 'g', "maxRearrangementRatio", "Maximum rearrangement ratio.");
	MegaXParser.add_option(args.disableEstimationRatio, 'E', "disableEstimationRatio", "Disable estimation ratio.");
	MegaXParser.add_option(args.disableRearrangement, 'R', "disableRearrangement", "Disable rearrangement ratio.");
    MegaXParser.add_option(args.numberOfHashFunctions, 'a', "numberOfHashFunctions", "Number of used hashfunctions.");
	MegaXParser.add_option(args.maxUserBins, 'M', "maxUserBins", "Maximum number of user bins in each HIBF.");
	MegaXParser.add_option(args.hibfThreshold, 'd', "hibfThreshold", "Threshold to assign query to user bin.");
	MegaXParser.add_option(args.filterFileName, 'F', "filterFileName", "HIBF output file name.");
	MegaXParser.add_option(args.vectFileName, 'v', "vectFileName", "Mutated DB input file name.");

	/// FM Index params
	MegaXParser.add_option(args.indexFileName, 'X', "indexFileName", "Output index file name.");
	MegaXParser.add_option(args.numberOfReferenceSequences, 'n', "numberOfReferenceSequences", "Number of reference sequences.");

	/// query params
	MegaXParser.add_option(args.queryFileName, 'f', "queryFileName", "Query file name.");

	/// evaluation params
	MegaXParser.add_option(args.countingResults, 'C', "countingResults", "Counting results file name.");
	MegaXParser.add_option(args.referenceMapping, 'I', "referenceMapping", "Reference ID mapping file.");
	MegaXParser.add_option(args.evaluationResults, 'V', "evaluationResults", "Output file name for evaluation results.");
	MegaXParser.add_option(args.threshold, 'D', "threshold", "Mapping threshold to assign query as part of sequence.");
	MegaXParser.add_option(args.taxClassificationFile, 'T', "taxClassificationFile", "Output file name for strain level classification.");
	MegaXParser.add_option(args.blackListFile, 'b', "blackListFile", "Path to black list file.");
	
}

/*
 * @fn runProgram
 * @brief run the whole tool.
 * @signature void runProgram(mdArguments &args);
 * @param args user parameters struct.
 * @return None.
 * @throws None.
 */
void runProgram(CmdArguments &args)
{
    

	if (std::string("build_vect").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath);
        seqan3::debug_stream << '\n';
        seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
		seqan3::debug_stream << "      Input Parameters for Mutating DB             " << " " << '\n';
		seqan3::debug_stream << "Viruses input file: " << args.sequenceFilePath << '\n';
		seqan3::debug_stream << "Matrix File: " << args.matrixFilePath << '\n';
		seqan3::debug_stream << "(Threads: " << unsigned(args.numberOfThreads) << ", QMerSize: " << unsigned(args.qMerSize) << ", minScore: " << unsigned(args.minScore) << ")" << '\n';
        seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
		seqan3::debug_stream << '\n';

		double start {omp_get_wtime()};
		mutate.mutateSequencesRec(args.qMerSize, args.sequenceFilePath, args.minScore, args.numberOfThreads, args.outputDir, args.minimiser, args.windowSize);
		double end {omp_get_wtime()};
		
		long int mutatedQMers = (mutate.totalNumberOfQMers - mutate.totalNumberOfOriginalQMers);
		std::cout << "[INFO] Elapsed time to mutate the input sequences: " << end - start << "s" << std::endl;
        std::cout<< "Number of original sequences: " << mutate.numberOfInputSequences << std::endl;
		std::cout<< "Number of original q_mers: " << mutate.totalNumberOfOriginalQMers << std::endl;
		std::cout<< "Number of mutated q_mers: " << mutatedQMers << std::endl;
        std::cout<< "Total number of q_mers: " << mutate.totalNumberOfQMers << std::endl;
        seqan3::debug_stream << "###################################################################################################################################" << '\n';
		std::cout << std::endl;

    }


	else if (std::string("build_hibf").compare(args.mode) == 0)
	{
		HIBF HIBFClassObject;
		HIBFParams Params = {args.numberOfHashFunctions, args.numberOfThreads, args.sketchBits, args.maxUserBins, args.alpha,
							 args.maximumFalsePositiveRate, args.maxRearrangementRatio, args.disableEstimationRatio,
							 args.disableRearrangement, args.vectFileName, args.filterFileName};

		HIBFClassObject.initialiseParamsHIBF(Params);
		HIBFClassObject.printInputParams();

		auto start = std::chrono::high_resolution_clock::now();
		HIBFClassObject.buildHIBF();
		auto finish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time (for building HIBF)           : " << elapsed.count() << " s\n";
		std::cout << std::endl;
    }


	else if (std::string("count_hibf").compare(args.mode) == 0)
	{
		std::cout << "[WARNING] HIBF counting results give only the bins, where the query is reported and without yielding the nmber of hits!" << std::endl;
		HIBF HIBFCountintObject;
		HIBFCountintObject.countHIBF(args.filterFileName, args.queryFileName, args.qMerSize, args.hibfThreshold);
    }
	
	else if (std::string("build_ibf").compare(args.mode) == 0)
	{
		HIBF IBFClassObject;
		IBFClassObject.buildIBF(args.numberOfHashFunctions, args.vectFileName, args.outputDir);
    }

	else if (std::string("count_ibf").compare(args.mode) == 0)
	{
		HIBF IBFCountintObject;
		IBFCountintObject.countIBF(args.filterFileName, args.queryFileName, args.qMerSize);
    }


	else if (std::string("hibf_ref").compare(args.mode) == 0)
	{
		HIBF HIBFClassObject; 
		HIBFParams Params = {args.numberOfHashFunctions, args.numberOfThreads, args.sketchBits, args.maxUserBins, args.alpha,
							 args.maximumFalsePositiveRate, args.maxRearrangementRatio, args.disableEstimationRatio,
							 args.disableRearrangement, args.vectFileName, args.filterFileName};

		HIBFClassObject.initialiseParamsHIBF(Params);
		HIBFClassObject.printInputParams();

		auto start = std::chrono::high_resolution_clock::now();
		HIBFClassObject.buildHIBFRef(args.sequenceFilePath, args.qMerSize);
		auto finish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time (for building HIBF)           : " << elapsed.count() << " s\n";
		std::cout << std::endl;
    }


	else if (std::string("binary_search").compare(args.mode) == 0)
	{
		std::cout << "[INFO] Start binary search algorithm!" << std::endl;
		binarySearch(args.vectFileName, args.queryFileName, args.qMerSize, args.numberOfThreads);
    }

	else if (std::string("build_fm").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath);

		FM_Index_MX FMIdxClassObject;
		FMParams Params = {args.qMerSize, args.numberOfThreads, args.numberOfReferenceSequences, args.sequenceFilePath, args.indexFileName};

		FMIdxClassObject.initialiseParamsIndex(Params);
		FMIdxClassObject.printInputParams();

		auto start = std::chrono::high_resolution_clock::now();
		FMIdxClassObject.buildFM(args.matrixFilePath, args.minScore);
		auto finish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time (for building/writing FM index)           : " << elapsed.count() << " s\n";
		std::cout << std::endl;
    }

	else if (std::string("count_fm").compare(args.mode) == 0)
	{
		FM_Index_MX FMIdxClassObject;

		auto start = std::chrono::high_resolution_clock::now();
		FMIdxClassObject.countFM(args.indexFileName, args.queryFileName, args.qMerSize, args.numberOfReferenceSequences, args.numberOfThreads);
		auto finish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time (for counting in FM index)        : " << elapsed.count() << " s\n";
		std::cout << std::endl;
    }


	else if (std::string("stat").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath);
        mutate.printSequenceStatstics(args.sequenceFilePath, args.qMerSize, args.outputDir);
    }

	else if (std::string("mutate_stat").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath);
        mutate.printMutationStatstics(args.sequenceFilePath, args.qMerSize, args.minScore, args.minimiser);
    }

	else if (std::string("evaluate").compare(args.mode) == 0)
	{
        resultsMapping(args.countingResults, args.referenceMapping, args.evaluationResults, args.queryFileName, args.threshold, args.qMerSize);
    }

	else if (std::string("classification").compare(args.mode) == 0)
	{
        taxaClassification(args.evaluationResults, args.taxClassificationFile);
    }
	else if (std::string("ibf_stat").compare(args.mode) == 0)
	{
        Mutate mutate(args.matrixFilePath);

		mutate.mutateSequencesRec(args.qMerSize, args.sequenceFilePath, args.minScore, args.numberOfThreads, args.outputDir, args.minimiser, args.windowSize);
		std::vector<std::vector<uint64_t>> userBins;

		std::ifstream inputVect(args.vectFileName, std::ios::binary); 
		cereal::BinaryInputArchive archiveInput(inputVect);              
		archiveInput(userBins); 

		uint64_t bins = userBins.size();
		seqan3::debug_stream << "[INFO] Finished loading user bins." << '\n';
		seqan3::debug_stream << "[INFO] Number of input user bins: " << bins << '\n';

		std::size_t maxSize = 0;
		for (const auto& userBin : userBins) {
			if (userBin.size() > maxSize) {
				maxSize = userBin.size();
			}
		}
		uint64_t numberOfKmers = uint64_t(maxSize);

		uint64_t binSize = ceil(-1 / (pow(1 - pow((double) 0.01, 1.0 / (double) args.numberOfHashFunctions), 1.0 / ((double) (args.numberOfHashFunctions * numberOfKmers))) - 1));

		auto IBFMB = (double)(binSize * bins)/ static_cast< double >( 8388608u );
		auto IBFGB = IBFMB / 1024.0;

		std::cout << "Main Params: (" << unsigned(args.qMerSize) << ", " << unsigned(args.minScore) << ", " << args.sequenceFilePath <<
		 ")" << std::endl;
		std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;
		std::cout << "[INFO-DEV] Bloom Filter Stat   : "<< std::endl;
		std::cout << "Number of user bins            : " << unsigned(bins) << std::endl;
		std::cout << "Max k-mer size                 : " << unsigned(numberOfKmers) << std::endl;
		std::cout << "Bin size                       : " << unsigned(binSize) << std::endl;
		std::cout << "IBF size in MB                 : " <<  IBFMB << " MBytes"<< std::endl;
		std::cout << "IBF size in GB                 : " << std::fixed << std::setprecision(2) << IBFGB << " GBytes" << std::endl;
		std::cout << "Number of hash functions       : " << unsigned(args.numberOfHashFunctions) << std::endl;
		std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;
    }

	else if (std::string("write_db").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath); 
        mutate.writeDBToFasta(args.qMerSize, args.sequenceFilePath, args.minScore, args.numberOfThreads, args.outputDir);
    }

	else if (std::string("simulate_peptides").compare(args.mode) == 0)
	{
		uint8_t errors {args.hibfThreshold};

		simulate(args.numberOfReferenceSequences, errors, 0, 5, 50, args.sequenceFilePath, args.outputDir);
    }

	else if (std::string("test").compare(args.mode) == 0)
	{

	 //test_unordered_set();
		//write_db();
		//test_minimiser();
    }

	else if (std::string("mutate_seq_len").compare(args.mode) == 0)
	{
		Mutate mutate(args.matrixFilePath);
        mutate.mutateToSequenceLength(args.numberOfTopSequences, args.sequenceFilePath, args.outputDir, args.minScore);
    }

	// ./megaX multi_indexing -m blosum62 -b black_list.txt -i refSeqViral_monkeypox.fasta -f monkeypox_sample/E02292_MonkeyPox_SP3_DDA_1.fasta -q 5 -s 100 -t 40 -Z 0 -a 2 -M 1000 -F monckey.log -D 0.8
	else if (std::string("multi_indexing").compare(args.mode) == 0)
	{
		HIBF MultiIndexingObj;
		MultiIndexingObj.multiIndexing(args.matrixFilePath, args.blackListFile, args.sequenceFilePath, args.queryFileName, args.qMerSize, args.minScore,
		 args.numberOfThreads, args.minimiser, args.windowSize, args.numberOfHashFunctions, args.maxUserBins, args.filterFileName, args.threshold);
    }
    
}

#endif /* PARSER_HPP_ */
