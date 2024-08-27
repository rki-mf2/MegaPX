
#include "hibf_util.hpp"


//==================================================================
// Member functions.
//==================================================================

// Definition of the member constructor HIBF::HIBF. 
HIBF::HIBF(){}


// Definition of the member destructor HIBF::~HIBF. 
HIBF::~HIBF(){
std::cout<<"[INFO-DEV] Delete constructer of HIBF functionality....."<<std::endl;
}

struct FMParams
{
    /// IBF params struct
    uint32_t  bins{1};
    uint8_t   kMerSize{2};
    uint64_t    sequenceLength; 

    uint8_t   numberOfHashFunctions {3};

    std::string IBFFileName {""};
};


struct HIBFParams
{
    /// HIBF params struct
    uint8_t   numberOfHashFunctions {2};
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

/// only for IBF 
// Definition of the member function HIBF::countSequencesRef. 
uint32_t HIBF::countSequencesRef(const std::string& inputFastaFilePath) {
    std::ifstream sequenceFile(inputFastaFilePath);

    if (!sequenceFile.good())
        throw std::runtime_error("[ERROR] Failed to open sequence file.");

    uint32_t numberOfSequences = 0;
    std::string line;

    while (std::getline(sequenceFile, line)) {
        if (!line.empty() && line[0] == '>') {
            numberOfSequences++;
        }
    }

    sequenceFile.close();
    return numberOfSequences;
}

// Definition of the member function HIBF::computeBinSize.
// Source: https://github.com/JensUweUlrich/ReadBouncer/blob/1.1.0/src/IBF/IBFBuild.cpp#L400
uint64_t HIBF::computeBinSizeIBF(uint64_t numberOfKmers){

    uint64_t BinSizeBits = ceil(-1 / (pow(1 - pow((double) this->falsePositiveRate, 1.0 / (double) this->numberOfHashFunctions), 1.0 / ((double) (this->numberOfHashFunctions * numberOfKmers))) - 1));
    return BinSizeBits;
}


// Definition of the member function HIBF::buildIBF. 
void HIBF::buildIBF(uint8_t numberOfHashFunctionsIn, std::string& vectorFileName, std::string& outputDir){

    seqan3::debug_stream << "[INFO] Loading input user bins from file..." << '\n';

    std::vector<std::vector<uint64_t>> userBins;

    std::ifstream inputVect(vectorFileName, std::ios::binary); 
	cereal::BinaryInputArchive archiveInput(inputVect); 

    /**/
    while (true) {
        std::vector<uint64_t> userBin;
        try {
            archiveInput(userBin);
            userBins.push_back(userBin);
        } catch (const cereal::Exception& e) {
            break; // End of file reached
        }
    }


    this->bins = userBins.size();
    seqan3::debug_stream << "[INFO] Finished loading user bins." << '\n';
    seqan3::debug_stream << "[INFO] Number of input user bins: " << this->bins << '\n';
    
    if (numberOfHashFunctionsIn < 1 || numberOfHashFunctionsIn > 5){

        this-> numberOfHashFunctions = 2;
    }
    else{

        this->numberOfHashFunctions = numberOfHashFunctionsIn;
    }

    std::size_t maxSize = 0;
    for (const auto& userBin : userBins) {
        //seqan3::debug_stream << " Size of the user bin: "<< userBin.size() << std::endl;
        if (userBin.size() > maxSize) {
            maxSize = userBin.size();
        }
    }

    uint64_t numberOfKmers = uint64_t(maxSize);

    // Get file name without extension     
    auto fileName = std::filesystem::path(vectorFileName).stem();
    this->IBFFileName = outputDir + fileName.string() + "_" + std::to_string(this-> numberOfHashFunctions) + ".ibf";
    this->binSize = this->computeBinSizeIBF(numberOfKmers);

    //uint64_t optimalNumberOfBinsMegaX = optimalNumberOfBins();

    auto IBFMB = (double)(this->binSize * this->bins)/ static_cast< double >( 8388608u );
    auto IBFGB = IBFMB / 1024.0;

    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;
    std::cout << "[INFO-DEV] Bloom Filter Stat   : "<< std::endl;
    std::cout << "Number of user bins            : " << unsigned(this->bins) << std::endl;
    //std::cout << "Computed optimal number of bins: "<< unsigned(optimalNumberOfBinsMegaX) << std::endl;
    std::cout << "Max k-mer size                 : " << unsigned(numberOfKmers) << std::endl;
    std::cout << "Bin size                       : " << unsigned(this->binSize) << std::endl;
    std::cout << "IBF size in MB                 : " <<  IBFMB << " MBytes"<< std::endl;
    std::cout << "IBF size in GB                 : " << std::fixed << std::setprecision(2) << IBFGB << " GBytes" << std::endl;
    std::cout << "Number of hash functions       : " << unsigned(this->numberOfHashFunctions) << std::endl;
    std::cout << "IBF file name                  : " << this->IBFFileName << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;

    seqan3::debug_stream << "[INFO] Start building IBF configuration." << '\n';

    //@ToDo seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{this->bins}, // number of bins
                                              seqan::hibf::bin_size{binSize}, // each bin using this size
                                              seqan::hibf::hash_function_count{this->numberOfHashFunctions}};
    
    seqan3::debug_stream << "[INFO] Start building IBF..." << '\n';
    
    uint64_t userBinID = 0; 

    for (const auto& userBin : userBins){
        
        //seqan3::debug_stream << "Start inserting into bin: "<< userBinID << '\n';

        for (const auto& hashValue : userBin){

            ibf.emplace(hashValue, seqan::hibf::bin_index{userBinID});

        }
        userBinID += 1; 
    }

    seqan3::debug_stream << "[INFO] Finsihed building IBF." << '\n';
    seqan3::debug_stream << "[INFO] Start writing IBF to the file: " << this->IBFFileName << '\n';

    std::ofstream  outputIBF(this->IBFFileName , std::ios::binary); 
	cereal::BinaryOutputArchive archiveOutput(outputIBF); 
	archiveOutput(ibf);

    seqan3::debug_stream << "[INFO] Finsihed writing IBF to the target output file." << '\n';

}

// Definition of the member function HIBF::countIBF.
void HIBF::countIBF(const std::string& IBFFileName, const std::string& queryFileName, uint8_t kMerSize){

    seqan3::debug_stream << "[INFO] Start querying the de-novo sequences." << '\n';
    seqan3::debug_stream << "[INFO] Loading input IBF...." << '\n';

    seqan::hibf::interleaved_bloom_filter ibf;

    std::ifstream is(IBFFileName, std::ios::binary); 
    cereal::BinaryInputArchive archive(is);
    archive(ibf); 

    seqan3::debug_stream << "[INFO] Finished loading input IBF. " << "(" << ibf.bin_count() << ")" << '\n';
    seqan3::debug_stream << "[INFO] Loading input sequences...." << '\n';

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{queryFileName};

    uint64_t numberOfQueries {0u};
    uint64_t skippedQueries {0u};

    for (auto & record : fin)
    {
       numberOfQueries += 1;
    }

    seqan3::debug_stream << "[INFO] Number of query input sequences: " << numberOfQueries << '\n';

    std::vector<std::vector<uint8_t>> resultsVector(numberOfQueries, std::vector<uint8_t>(ibf.bin_count(), 0));

    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});

    auto agent = ibf.counting_agent();

    seqan3::sequence_file_input<traits_type> fin1{queryFileName}; // query loop! @ToDo make it faster without two jumps into cache and read file twice. 

    uint64_t queryIDX {0u};

    for (auto & record : fin1)
    {
        auto const & seq = record.sequence();

        if(seq.size() < kMerSize){
            skippedQueries++;
            continue;
        }
            
        //auto numberOfK = (seq | hash_adaptor);
        //seqan3::debug_stream << "Size of the query: "<< numberOfK.size() << '\n';
        auto counts = agent.bulk_count(seq | hash_adaptor);

        uint64_t referenceID {0u};
        for (const auto & count : counts){
            resultsVector[queryIDX][referenceID] = count;
            referenceID += 1;
        }
        queryIDX += 1;        
    }


    std::filesystem::path filePath(IBFFileName);
    std::filesystem::path queryFile(queryFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / (queryFile.filename().stem().string() + "_results_IBF.log");
    
    //std::ofstream logFile(resultsFile);

    //if (!logFile.is_open()) 
     // throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());
    
    seqan3::debug_stream << "[INFO] Number of skipped sequences (size < k): " << skippedQueries << '\n';
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
    }*
    /*
    for (const auto& row : resultsVector) {
        for (const auto& element : row) {
            logFile << element << ' ';
        }
        logFile << '\n';
    }*/
    resultsVector.clear();
    //logFile.close();

}


// Definition of the member function HIBF::buildSingleIBF. 
seqan::hibf::interleaved_bloom_filter HIBF::buildSingleIBF(std::vector<std::vector<uint64_t>> userBins, uint64_t numberOfHashFunctionsIn){

    this->bins = userBins.size();
    seqan3::debug_stream << "[INFO] Number of input user bins (single IBF): " << this->bins << '\n';
    
    if (numberOfHashFunctionsIn < 1 || numberOfHashFunctionsIn > 5){

        this-> numberOfHashFunctions = 2;
    }
    else{

        this->numberOfHashFunctions = numberOfHashFunctionsIn;
    }

    std::size_t maxSize = 0;
    for (const auto& userBin : userBins) {
        //seqan3::debug_stream << " Size of the user bin: "<< userBin.size() << std::endl;
        if (userBin.size() > maxSize) {
            maxSize = userBin.size();
        }
    }

    uint64_t numberOfKmers = uint64_t(maxSize);
    this->binSize = this->computeBinSizeIBF(numberOfKmers);
    //this->binSize = ceil(-1 / (pow(1 - pow((double) 0.09, 1.0 / (double) this->numberOfHashFunctions), 1.0 / ((double) (this->numberOfHashFunctions * numberOfKmers))) - 1));
    //this->binSize = 1048576;

    auto IBFMB = (double)(this->binSize * this->bins)/ static_cast< double >( 8388608u );
    auto IBFGB = IBFMB / 1024.0;
    if (IBFGB > this->maxIBFSize) this->maxIBFSize = IBFGB;

    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;
    std::cout << "[INFO-DEV] Bloom Filter Stat   : "<< std::endl;
    std::cout << "Number of user bins            : " << unsigned(this->bins) << std::endl;
    //std::cout << "Computed optimal number of bins: "<< unsigned(optimalNumberOfBinsMegaX) << std::endl;
    std::cout << "Max k-mer size                 : " << unsigned(numberOfKmers) << std::endl;
    std::cout << "Bin size                       : " << unsigned(this->binSize) << std::endl;
    std::cout << "IBF size in MB                 : " <<  IBFMB << " MBytes"<< std::endl;
    std::cout << "IBF size in GB                 : " << std::fixed << std::setprecision(2) << IBFGB << " GBytes" << std::endl;
    std::cout << "Number of hash functions       : " << unsigned(this->numberOfHashFunctions) << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<< std::endl;

    seqan3::debug_stream << "[INFO] Start building IBF configuration." << '\n';

    //@ToDo seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{this->bins}, // number of bins
                                              seqan::hibf::bin_size{binSize}, // each bin using this size
                                              seqan::hibf::hash_function_count{this->numberOfHashFunctions}};
    
    seqan3::debug_stream << "[INFO] Start building IBF..." << '\n';
    
    uint64_t userBinID = 0; 

    for (const auto& userBin : userBins){
        
        //seqan3::debug_stream << "Start inserting into bin: "<< userBinID << '\n';

        for (const auto& hashValue : userBin){

            ibf.emplace(hashValue, seqan::hibf::bin_index{userBinID});

        }
        userBinID += 1; 
    }

    seqan3::debug_stream << "[INFO] Finsihed building IBF." << '\n';

    return ibf;
}

// Definition of the member function HIBF::initialiseParamsHIBF. 
void HIBF::initialiseParamsHIBF(HIBFParams& userParams){

    seqan3::debug_stream << "[INFO] Start loading HIBF parameters..." << '\n';

    if (userParams.numberOfHashFunctions < 1){

        this->numberOfHashFunctions = 2; 
    }
    else {this->numberOfHashFunctions = userParams.numberOfHashFunctions;}

    if (userParams.threads < 1){

        this->threads = 1; 
    }
    else {this->threads = userParams.threads;}

    if (userParams.sketchBits < 1){

        this->sketchBits = 12; 
    }
    else {this->sketchBits = userParams.sketchBits;}


    if (userParams.alpha < 0){

        this->alpha = 1.2; 
    }
    else {this->alpha = userParams.alpha;}

    if (userParams.maxUserBins <= 0){

        this->maxUserBins = 1; 
    }
    else {this->maxUserBins = userParams.maxUserBins;}

    if (userParams.maximumFalsePositiveRate < 0){

        this->maximumFalsePositiveRate = 0.05; 
    }
    else {this->maximumFalsePositiveRate = userParams.maximumFalsePositiveRate;}

    if (userParams.maxRearrangementRatio < 0){

        this->maxRearrangementRatio = 0.5; 
    }
    else {this->maxRearrangementRatio = userParams.maxRearrangementRatio;}

    this->disableEstimationRatio = userParams.disableEstimationRatio;
    this->disableRearrangement = userParams.disableRearrangement;
    this->vectFileName = userParams.vectFileName;
    this->filterFileName = userParams.filterFileName;

}

// Definition of the member function HIBF::printInputParams.
void HIBF::printInputParams(){

    seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
	seqan3::debug_stream << "            Input Parameters for Building HIBF              " << " " << '\n';
	seqan3::debug_stream << "Input vector file               : " << this->vectFileName << '\n';
    seqan3::debug_stream << "Output HIBF file                : " << this->filterFileName << '\n';
    seqan3::debug_stream << "Disable estimation ratio        : " << this->disableEstimationRatio << '\n';
    seqan3::debug_stream << "Disable max rearrangement ratio : " << this->disableRearrangement << '\n';
    seqan3::debug_stream << "Max rearrangement ratio         : " << this->maxRearrangementRatio << '\n';
    seqan3::debug_stream << "Max false positive rate         : " << this->maximumFalsePositiveRate << '\n';
    seqan3::debug_stream << "Value of alpha                  : " << this->alpha << '\n';
    seqan3::debug_stream << "Maximum number of user bins HIBF: " << this->maxUserBins << '\n';
    seqan3::debug_stream << "Sketch bits                     : " << this->sketchBits << '\n';
    seqan3::debug_stream << "Number of building threads      : " << unsigned(this->threads) << '\n';
    seqan3::debug_stream << "Number of hash functions        : " << this->numberOfHashFunctions << '\n';
    seqan3::debug_stream << "----------------------------------------------------------------" << '\n';
}

// Definition of the member function HIBF::buildHIBF.
void HIBF::buildHIBF(){

    seqan3::debug_stream << "[INFO] Start building HIBF using the input parameters..." << '\n';
    seqan3::debug_stream << "[INFO] Loading input user bins from file: "<< this->vectFileName << '\n';

    std::vector<std::vector<uint64_t>> userBins;

    std::ifstream inputVect(this->vectFileName, std::ios::binary); 
	cereal::BinaryInputArchive archiveInput(inputVect); 

    /**/
    while (true) {
        std::vector<uint64_t> userBin;
        try {
            archiveInput(userBin);
            userBins.push_back(userBin);
        } catch (const cereal::Exception& e) {
            break; // End of file reached
        }
    }
    /*
    std::vector<std::vector<uint64_t>> userBins;
    std::ifstream inputVect(this->vectFileName, std::ios::binary); 
	cereal::BinaryInputArchive archiveInput(inputVect);              
	archiveInput(userBins); */

    seqan3::debug_stream << "[INFO] Finished loading user bins." << '\n';
    seqan3::debug_stream << "[INFO] Number of input user bins: " << userBins.size() << '\n';
    seqan3::debug_stream << "[INFO] Start splitting input user bins into different HIBF's. "<< '\n';

    size_t numberOfHIBFs = std::ceil(double(userBins.size())/double(this->maxUserBins));
    std::vector<std::vector<std::vector<uint64_t>>> splits;
    splits.resize(numberOfHIBFs);

    seqan3::debug_stream << "[INFO] Computed number of HIBF's: " << unsigned(numberOfHIBFs) << '\n';

    // Insertion into collection of HIBF's 
    uint8_t targetBin {0u};

    for (auto& userBin : userBins) {
        // seqan3::debug_stream << "User Bin Size DEV: " << userBin.size() << std::endl;

        // Check if the split is not full and there are remaining user bins
        if (splits[targetBin].size() < this->maxUserBins) {
            splits[targetBin].emplace_back(std::move(userBin));
        }
        
        if (splits[targetBin].size() == this->maxUserBins) {
            seqan3::debug_stream << "[INFO] Finished insertion into: " << targetBin << '\n';
            targetBin++;
        }
    }

    userBins.clear();

    seqan3::debug_stream << "[INFO] Finished all insertions."<< '\n';

    int splitCounts {0u};
    for (const auto& hibfSplit : splits){

        splitCounts++;
        seqan3::debug_stream << "[INFO] Split number: "<< splitCounts << " | size: " << hibfSplit.size()<< '\n';

    }

    std::vector<std::string> outputHIBFsPaths{};

    // Start split insertions!
    int userBinIDX = 0; 
    for (const auto& split : splits){

        seqan3::debug_stream << "[INFO] Start constructing HIBF number: " << userBinIDX << '\n';
        seqan3::debug_stream << "##########################################" << '\n';

        auto inputUserBins = [&](size_t const userBinID, seqan::hibf::insert_iterator it)
        {
            for (auto const hash : split[userBinID])
                it = hash;
        };

        std::size_t tMax = std::ceil(std::sqrt(split.size())/ 64) * 64;

        seqan3::debug_stream << "[INFO] MegaX computes the number of technical bins according to input user bins."<< '\n';
        seqan3::debug_stream << "[INFO] Number of technical bins in each HIBF layer: " << tMax << '\n';

        std::filesystem::path fPath(this->filterFileName);
        const std::string hibfTempName = fPath.parent_path() / (fPath.filename().stem().string() + "_split_" + std::to_string(userBinIDX) + ".hibf");

        outputHIBFsPaths.emplace_back(hibfTempName);
        seqan3::debug_stream << "[INFO] HIBF will be written to the file: " << hibfTempName << '\n';

        seqan3::debug_stream << "[INFO] Start building HIBF configuration." << '\n';
        seqan::hibf::config config{.input_fn = inputUserBins,     // required
                                .number_of_user_bins = split.size(), // required
                                .number_of_hash_functions = this->numberOfHashFunctions,
                                .maximum_fpr = this->maximumFalsePositiveRate, // recommended to adapt
                                .threads = this->threads,                        // recommended to adapt
                                .sketch_bits = this->sketchBits,
                                .tmax = tMax, // number of technical bins! if tmax is 0 this will be autoamtically computed. 
                                .alpha = this->alpha,
                                .max_rearrangement_ratio = this->maxRearrangementRatio,
                                .disable_estimate_union = this->disableEstimationRatio,
                                .disable_rearrangement = this->disableRearrangement};

        seqan3::debug_stream << "[INFO] Start building HIBF." << '\n';

        seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

        seqan3::debug_stream << "[INFO] Finsihed building HIBF." << '\n';
        seqan3::debug_stream << "[INFO] Start writing HIBF to the file: " << hibfTempName << '\n';
        
        std::ofstream  outputHIBF(hibfTempName , std::ios::binary); 
        cereal::BinaryOutputArchive archiveOutput(outputHIBF); 
        archiveOutput(hibf);

        seqan3::debug_stream << "[INFO] Finsihed writing HIBF to the target output file." << '\n';
        seqan3::debug_stream << "[INFO] Finished constructing HIBF number: " << userBinIDX << '\n';
        seqan3::debug_stream << "##########################################" << '\n';

        userBinIDX += 1;
    }

    std::filesystem::path filePath(this->filterFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / ((filePath.stem()).string() + "_HIBF_files.log");
    std::ofstream logFile(resultsFile);

    if (!logFile.is_open()) 
      throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());
    
    seqan3::debug_stream << "[INFO] Start writing hibfs paths to the output file: " << resultsFile.string() << '\n';

    for (const auto& path : outputHIBFsPaths){

        logFile << path;
        logFile << '\n';
    }

    logFile.close();

}

/// Definition of the function HIBF::buildHIBFRef
void HIBF::buildHIBFRef(const std::string& inputFastaFile, uint8_t kMerSize){

    seqan3::debug_stream << "[INFO] Start building HIBF using the input parameters..." << '\n';
    seqan3::debug_stream << "[INFO] Loading input user bins from file..." << '\n';

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{inputFastaFile};
    std::vector<seqan3::aa27_vector> userBins; 

    for (auto & record : fin)
    {
        userBins.push_back(record.sequence());
    }

    seqan3::debug_stream << "[INFO] Finished loading user bins." << '\n';
    seqan3::debug_stream << "[INFO] Number of input user bins: " << userBins.size() << '\n';
    seqan3::debug_stream << "[INFO] Start splitting input user bins into different HIBF's. "<< '\n';

    size_t numberOfHIBFs = std::ceil(double(userBins.size())/double(this->maxUserBins));
    std::vector<std::vector<seqan3::aa27_vector>> splits;
    splits.resize(numberOfHIBFs);

    seqan3::debug_stream << "[INFO] Computed number of HIBF's: " << unsigned(numberOfHIBFs) << '\n';

     // Insertion into collection of HIBF's 
    uint8_t targetBin {0u};

    for (auto& userBin : userBins) {
        // seqan3::debug_stream << "User Bin Size DEV: " << userBin.size() << std::endl;

        // Check if the split is not full and there are remaining user bins
        if (splits[targetBin].size() < this->maxUserBins) {
            splits[targetBin].emplace_back(std::move(userBin));
        }
        
        if (splits[targetBin].size() == this->maxUserBins) {
            seqan3::debug_stream << "[INFO] Finished insertion into: " << targetBin << '\n';
            targetBin++;
        }
    }

    userBins.clear();
    seqan3::debug_stream << "[INFO] Finished all insertions."<< '\n';

    int splitCounts {0u};
    for (const auto& hibfSplit : splits){

        splitCounts++;
        seqan3::debug_stream << "[INFO] Split number: "<< splitCounts << " | size: " << hibfSplit.size()<< '\n';

    }

    std::vector<std::string> outputHIBFsPaths{};

    // Start split insertions!
    int userBinIDX = 0; 
    for (const auto& split : splits){

        seqan3::debug_stream << "[INFO] Start constructing HIBF number: " << userBinIDX << '\n';
        seqan3::debug_stream << "##########################################" << '\n';

       auto inputUserBins = [&](size_t const userBinID, seqan::hibf::insert_iterator it)
        {
            auto const & seq = split[userBinID];

            auto hashAdaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
            auto hashes = (seq | hashAdaptor);

            //seqan3::debug_stream << "Inserting into the user bin number: " <<  userBinID << std::endl;
            for (auto const & hash : hashes)
                it = hash; 

        };

        std::size_t tMax = std::ceil(std::sqrt(split.size())/ 64) * 64;

        seqan3::debug_stream << "[INFO] MegaX computes the number of technical bins according to input user bins."<< '\n';
        seqan3::debug_stream << "[INFO] Number of technical bins in each HIBF layer: " << tMax << '\n';

        std::filesystem::path fPath(this->filterFileName);
        const std::string hibfTempName = fPath.parent_path() / (fPath.filename().stem().string() + "_split_" + std::to_string(userBinIDX) + ".hibf");

        outputHIBFsPaths.emplace_back(hibfTempName);
        seqan3::debug_stream << "[INFO] HIBF will be written to the file: " << hibfTempName << '\n';

        seqan3::debug_stream << "[INFO] Start building HIBF configuration." << '\n';
        seqan::hibf::config config{.input_fn = inputUserBins,     // required
                                .number_of_user_bins = split.size(), // required
                                .number_of_hash_functions = this->numberOfHashFunctions,
                                .maximum_fpr = this->maximumFalsePositiveRate, // recommended to adapt
                                .threads = this->threads,                        // recommended to adapt
                                .sketch_bits = this->sketchBits,
                                .tmax = tMax, // number of technical bins! if tmax is 0 this will be autoamtically computed. 
                                .alpha = this->alpha,
                                .max_rearrangement_ratio = this->maxRearrangementRatio,
                                .disable_estimate_union = this->disableEstimationRatio,
                                .disable_rearrangement = this->disableRearrangement};

        seqan3::debug_stream << "[INFO] Start building HIBF." << '\n';

        seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

        seqan3::debug_stream << "[INFO] Finsihed building HIBF." << '\n';
        seqan3::debug_stream << "[INFO] Start writing HIBF to the file: " << hibfTempName << '\n';
        
        std::ofstream  outputHIBF(hibfTempName , std::ios::binary); 
        cereal::BinaryOutputArchive archiveOutput(outputHIBF); 
        archiveOutput(hibf);

        seqan3::debug_stream << "[INFO] Finsihed writing HIBF to the target output file." << '\n';
        seqan3::debug_stream << "[INFO] Finished constructing HIBF number: " << userBinIDX << '\n';
        seqan3::debug_stream << "##########################################" << '\n';

        userBinIDX += 1;
    }

    std::filesystem::path filePath(this->filterFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / ((filePath.stem()).string() + "_HIBF_files.log");
    std::ofstream logFile(resultsFile);

    if (!logFile.is_open()) 
      throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());
    
    seqan3::debug_stream << "[INFO] Start writing hibfs paths to the output file: " << resultsFile.string() << '\n';

    for (const auto& path : outputHIBFsPaths){

        logFile << path;
        logFile << '\n';
    }

    logFile.close();

}

// Definition of the function HIBF::countHIBF
void HIBF::countHIBF(const std::string& HIBFFileName, const std::string& queryFileName, uint8_t kMerSize, uint8_t hibfThreshold){

    seqan3::debug_stream << "[INFO] Start querying the de-novo sequences." << '\n';

    std::vector<std::string> outputHIBFsCollection{};

    std::string path;

    std::ifstream pathsFile (HIBFFileName);
    std::filesystem::path filePath_(HIBFFileName);
    //const std::string pathToHIBF = (filePath_.parent_path()).string();

    // std::ifstream pathsFile (HIBFFileName);
    
    if (pathsFile.is_open())
    {
        while (std::getline (pathsFile,path) )
        {
            auto path_ = filePath_.parent_path() / path;
            outputHIBFsCollection.emplace_back(path_);
        }
        pathsFile.close();
    }

    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{queryFileName};

    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
    
    int queryCounter {1u};

    std::vector<std::vector<uint8_t>> resultsVector{};
    seqan3::debug_stream << "[WARNING] We use the counting agent from IBF counting, this process will be slower than only using threshold!" << '\n';

    seqan3::debug_stream << "[INFO] Processing query sequences..." << '\n';

    for (auto & record : fin)
    {
        queryCounter += 1;

        std::vector<uint8_t> results;
        std::vector<uint64_t> queryHashes; 

        auto const & seq = record.sequence();

        if(seq.size() < kMerSize)
            continue;

        auto hashValue = (seq | hash_adaptor);

        for (auto const & hash : hashValue)
            queryHashes.emplace_back(hash);

        for (const auto& HIBFFileName : outputHIBFsCollection){

            //seqan3::debug_stream << "[INFO-DEV] Querying HIBF: " << HIBFFileName << '\n'; 
            seqan::hibf::hierarchical_interleaved_bloom_filter hibf;
            std::ifstream is(HIBFFileName, std::ios::binary); 
            cereal::BinaryInputArchive archive(is);
            archive(hibf); 

            // In this case we will assign results using threshold! we can use counting but will use more memory! 
            // ToDo Add to thesis describtion 
            //auto agent = hibf.membership_agent(); 

            // The membership_for function takes the query and a threshold. Here, a threshold of two means that
            // at least (>=) 2 values of the query must be found within a user bin to be a hit.
            //auto & result = agent.membership_for(queryHashes, hibfThreshold);

            auto agent = hibf.template counting_agent<uint64_t>();
            auto & result = agent.bulk_count(queryHashes, hibfThreshold);
 
            if (!result.empty())
            {
                for (size_t i = 0u; i < result.size() - 1u; ++i)
                {
                    results.emplace_back(result[i]);
                }
                results.emplace_back(result.back());
            }

        }

        resultsVector.emplace_back(std::move(results));
        
    }

    /*
    for (size_t i = 0; i < resultsVector.size(); ++i) {
        for (size_t j = 0; j < resultsVector[i].size(); ++j) {
            std::cout << resultsVector[i][j] << " ";
        }
        std::cout << std::endl; // Move to the next row
    }*/

    std::filesystem::path filePath(HIBFFileName);
    std::filesystem::path queryFile(queryFileName);
    std::filesystem::path resultsFile = filePath.parent_path() / (queryFile.filename().stem().string() + "_results_HIBF.log");
    //std::ofstream logFile(resultsFile);

    //if (!logFile.is_open()) 
    //  throw std::runtime_error("[ERROR] Couldn't open the file: " + resultsFile.string());
    
    seqan3::debug_stream << "[INFO] Start writing results to the output file: " << resultsFile.string() << '\n';

    std::ofstream  outputLog(resultsFile.string() , std::ios::binary); 
	cereal::BinaryOutputArchive archiveOutput(outputLog); 
	archiveOutput(resultsVector);

    resultsVector.clear();

    /*for (const auto& row : resultsVector) {
        for (const auto& element : row) {
            logFile << element << ' ';
        }
        logFile << '\n';
    }

    logFile.close();*/

}

//Definition of the function HIBF::multiIndexing
void HIBF::multiIndexing(const std::string& matrixFilePath, const std::string& blackList, const std::string& inputFastaFile, const std::string& queryFileName, uint8_t kMerSize,
                         int minScore, uint8_t threads, bool minimiser, uint8_t windowSize, uint64_t numberOfHashFunctionsIn, uint64_t splitSize, const std::string& outputFile, double threshold)
{

    seqan3::debug_stream << "[INFO] Loading blacklist file..." << '\n';
    std::vector<std::string> outliers;
    this->loadBlackList(blackList, outliers);
    Mutate mutate(matrixFilePath);
    std::vector<std::vector<uint64_t>> userBins; 
    seqan3::debug_stream << "[INFO] Loading input user bins from file..." << '\n';

    using traits_type = seqan3::sequence_file_input_default_traits_aa;

    seqan3::sequence_file_input<traits_type> fin_query_count{queryFileName};
    seqan3::sequence_file_input<traits_type> fin_ref{inputFastaFile};
    seqan3::sequence_file_input<traits_type> fin_ref_count{inputFastaFile};

    uint64_t numberOfReferences {0u};
    uint64_t numberOfQueries {0u};
    uint64_t skippedQueries {0u};

    seqan3::debug_stream << "[INFO] Split size: " << splitSize << '\n';

    seqan3::debug_stream << "[INFO] Counting number of input peptides... " << '\n';
    for (auto & record : fin_query_count) numberOfQueries += 1; 
    seqan3::debug_stream << "[INFO] Counting number of input proteins... " << '\n';
    for (auto & record : fin_ref_count) numberOfReferences += 1; 

    seqan3::debug_stream << "[INFO] Number of queries: " << numberOfQueries << '\n';
    seqan3::debug_stream << "[INFO] Number of references: " << numberOfReferences << '\n';

    int numberOfOutliers {0u};
    int currentIBF {0u};
    uint64_t sequencesCounter {0u};
    uint64_t restSequences {numberOfReferences};
    uint64_t currentInsertedSequences{0u};
    uint64_t beforeInsertedSequences{0u};
    uint64_t currentReferencesCounter{0u};

    std::vector<std::vector<uint64_t>> queries; 
    seqan3::sequence_file_input<traits_type> fin{queryFileName};

    std::ofstream assignedPeptides;
    assignedPeptides.open(outputFile);

    for (auto & record : fin)
    {

        auto const & seq = record.sequence();
        if(seq.size() < kMerSize){
                skippedQueries++;
                continue;
        }
        std::vector<uint64_t> queryHashes;
        auto hashAdaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
        auto hashes = (seq | hashAdaptor);
        for (auto const & hash : hashes) queryHashes.push_back(hash);
        queries.push_back(queryHashes);
    }

    seqan3::debug_stream << "[INFO] Finished loading queries! " << numberOfQueries << '\n';
    std::string type = "RVDB"; 
    std::vector<std::string> headers;
    for (auto & record : fin_ref)
    {

        if ((record.sequence()).size() < kMerSize){
            numberOfReferences -= 1;
            seqan3::debug_stream << "[STEP-ERROR] Sequence is shorter than K-mer, we will skip the sequence: " << record.sequence() << '\n';
            seqan3::debug_stream << "[STEP-ERROR] ID: " << record.id() << '\n';
            continue;
        }
        if (headers.size() == splitSize) headers.clear();
        sequencesCounter += 1; 
        currentReferencesCounter += 1;
        auto const & seq = record.sequence();
        auto id = record.id(); 

        // Outliers detection
        bool outlierDetection = false;
        for (const auto& outlier : outliers) {
            if (id.find(outlier) != std::string::npos) {
                outlierDetection = true;
                break;
            }
        }
        if (outlierDetection){
            numberOfOutliers++;
            //seqan3::debug_stream << "[INFO] Detected outlier in the header: "<< id << '\n';
            continue;
        }

        headers.push_back(id);
        
        std::string sequence;
        for (auto const & c : seq) sequence += seqan3::to_char(c);
        std::vector<uint64_t> userBin = mutate.mutateSingleSequenceHashes(kMerSize, sequence, minScore, threads, minimiser, windowSize);
        userBins.push_back(userBin);

        if (sequencesCounter == splitSize || (currentReferencesCounter == numberOfReferences)){
            std::vector<std::vector<uint8_t>> resultsVector(numberOfQueries, std::vector<uint8_t>(userBins.size(), 0));
            restSequences -= splitSize;
            currentIBF += 1;
            sequencesCounter = 0;
            currentInsertedSequences = userBins.size();
            seqan3::debug_stream << "[INFO] Building IBF number: "<< currentIBF << " || (" << currentInsertedSequences << " bins)" << '\n';
            seqan::hibf::interleaved_bloom_filter IBF = this->buildSingleIBF(userBins, numberOfHashFunctionsIn);

            userBins.clear();
            seqan3::debug_stream << "[INFO] Start querying the de-novo sequences." << '\n';
            int referenceCheckerLog {0u};
            if (currentIBF > 1) referenceCheckerLog += (beforeInsertedSequences);
            seqan3::debug_stream << "[INFO] Query assignment will start with: " << referenceCheckerLog << '\n';
            auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{kMerSize});
            auto agent = IBF.counting_agent();
            uint64_t queryIDX {0u};

            for (auto& query : queries)
            {

                uint64_t numberOfKMers = (query.size());
	            uint64_t content = static_cast<int>(std::floor(threshold * numberOfKMers));
                auto counts = agent.bulk_count(query);
                uint64_t referenceID {0u};

                for (const auto & count : counts){
                    if (count >= content)
                        {
                            resultsVector[queryIDX][referenceID] = 1;
                        }
                        else
                        {
                            resultsVector[queryIDX][referenceID] = 0;
                        }
                    referenceID += 1;
                }
                queryIDX += 1;        
            }
            beforeInsertedSequences += currentInsertedSequences;

            seqan3::debug_stream << "[INFO] Processing results vector .... " << '\n';

            std::vector<int> bits(currentInsertedSequences, 0);
            auto numberOfPeptides = resultsVector.size();

            // Checkpoint 1
            if (currentInsertedSequences > resultsVector[0].size() || currentInsertedSequences > headers.size()) {
                // Handle the error, e.g., log and exit or throw an exception
                std::cerr << "Error: currentInsertedSequences exceeds the size of resultsVector or headers\n";
                
            }

            for (int col = 0; col < currentInsertedSequences; ++col)
            {
                for (int row = 0; row < numberOfPeptides; ++row)
                {
                    bits[col] += resultsVector[row][col];
                }
                //assignedPeptides << (headers[col]) << " | " << static_cast<double>(bits[col]) << "/" << numberOfPeptides << " (" << (static_cast<double>(bits[col])/numberOfPeptides)*100 << " %) \n";
                assignedPeptides << (headers[col]) << "\t" << static_cast<double>(bits[col]) << "\t" << numberOfPeptides << "\t" << (static_cast<double>(bits[col])/numberOfPeptides)*100 << " %) \n";
            }

            seqan3::debug_stream << "[INFO] Writing results vector to file: " << outputFile << '\n';
            bits.clear();
            resultsVector.clear();
            seqan3::debug_stream << "[INFO] Finished processing IBF number: " << currentIBF << std::endl;
        }
        
    }
    seqan3::debug_stream << "[INFO] Maximum IBF Size: " << this->maxIBFSize << " GBytes" << std::endl;
    seqan3::debug_stream << "[INFO] Total number of detected outliers: " <<numberOfOutliers << std::endl;

}

//Definition of the function HIBF::screening
void HIBF::screening(const std::string& inputFastaFile, const std::string& outputScreen, uint8_t kMerSize, double threshold){


}