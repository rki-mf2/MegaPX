#include <regex>


// Global name holder
std::string assignedPeptidesFile;


/*
 * @fn computeContent
 * @brief compute k-mers content of user-defined threshold.
 * @signature uint64_t computeContent(const std::string& query, double threshold, uint8_t qMerSize);
 * @param query query string.
 * @param threshold classification threshold [default:0.75].
 * @param qMerSize k-mer size for classification [default:2].
 * @return number of k-mers content for tax.classification.
 * @throws None.
 */
uint64_t computeContent(const std::string& query, double threshold, uint8_t qMerSize){

	uint64_t numberOfKMers = (query.length() - qMerSize + 1);
	uint64_t thresholdNumber = static_cast<int>(std::floor(threshold * numberOfKMers));

	return thresholdNumber;
}

/*
 * @fn estimation
 * @brief estimate k-mer contents using ML.
 * @signature std::vector<double> estimation(const std::vector<std::vector<uint8_t>>& resultsVector);
 * @param resultsVector lookup table.
 * @return std::vector<double> vector contains the estimations of each bin.
 * @throws None.
 */
std::vector<double> estimation(const std::vector<std::vector<uint8_t>>& resultsVector){
    
    std::vector<double> estimations;
    std::vector<size_t> bits(resultsVector[0].size(), 0);
    auto numberOfPeptides = resultsVector.size();

    for (size_t col = 0; col < resultsVector[0].size(); ++col)
    {
        for (size_t row = 0; row < numberOfPeptides; ++row)
        {
            bits[col] += resultsVector[row][col];
        }
    }

    estimations.resize(resultsVector[0].size());

    std::ofstream assignedPeptides;
    assignedPeptides.open (assignedPeptidesFile);

    for (size_t col = 0; col < estimations.size(); ++col)
    {
        estimations[col] = static_cast<double>(bits[col]) / numberOfPeptides;
        assignedPeptides << "Index " << col << " | " << static_cast<double>(bits[col]) << "/" << numberOfPeptides << " (" << (static_cast<double>(bits[col])/numberOfPeptides)*100 << " %) \n";
        //std::cout << estimations[col] << std::endl;
    }

    //resultsVector.clear();

    return estimations;
}

/*
 * @fn resultsMapping
 * @brief mapping counting results to reference ID.
 * @signature void resultsMapping(const std::string& countingResults, const std::string& referenceMapping, const std::string& evaluationResults, const std::string& queryFilePath, double threshold, uint8_t qMerSiz);
 * @param countingResults path to counting results file.
 * @param referenceMapping path to reference ID's file.
 * @param evaluationResultsFile path to output classification results file.
 * @param queryFilePath path to the query sequence.
 * @param threshold classification threshold [default:0.75].
 * @param qMerSize k-mer size for classification [default:2].
 * @return None.
 * @throws std::run_time.
 */
void resultsMapping(const std::string& countingResults, const std::string& referenceMapping, const std::string& evaluationResultsFile, const std::string& queryFilePath, double threshold, uint8_t qMerSize){


    std::string outputLengths = referenceMapping;
    assignedPeptidesFile = referenceMapping;;

    outputLengths = std::regex_replace(outputLengths, std::regex("reference_id_map"), "lengths");
    assignedPeptidesFile = std::regex_replace(assignedPeptidesFile, std::regex("reference_id_map"), "assigned_peptides_idx");

	std::cout << "[INFO] Start evaluating counting results..." << std::endl;
	std::cout << "Query file              : " << queryFilePath << std::endl;
	std::cout << "Counts file             : " << countingResults << std::endl;
	std::cout << "Mapping file            : " << referenceMapping << std::endl;
    std::cout << "Lengths file            : " << outputLengths << std::endl;
    std::cout << "Assignemnts file        : " << assignedPeptidesFile << std::endl;
	std::cout << "Classification threshold: " << threshold << std::endl;
	std::cout << "K-mer size              : " << static_cast<unsigned>(qMerSize) << std::endl;
	std::cout << '\n';

	// Read reference names from reference file
    std::vector<uint64_t> referenceLengths; 
    std::map<int, std::string> referenceMap;
    std::ifstream referenceStream(referenceMapping);
    std::ifstream lengthStream(outputLengths);

	if (referenceStream.is_open()) {
        std::string line;
        while (std::getline(referenceStream, line)) {
            std::istringstream iss(line);
            int id;
            std::string name;
            if (iss >> id >> std::ws && std::getline(iss, name)) {
                referenceMap[id] = name;
            } else {
                std::cerr << "Error parsing line: " << line << std::endl;
            }
        }
        referenceStream.close();
    } else {
        std::cerr << "Error opening reference file" << std::endl;
    }

    std::cout << "[INFO] Finshed loading reference headers!" << std::endl;
    std::cout << "[INFO] Start loading binary results ..." << std::endl;

	std::ifstream is(countingResults, std::ios::binary);
    cereal::BinaryInputArchive archive(is);
    std::vector<std::vector<uint8_t>> resultsVector;
    archive(resultsVector);

    std::cout << "[INFO] Finished loading binary archive!" << std::endl; 
	using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{queryFilePath};

    /*
    for (auto const &row : resultsVector)
    {
        seqan3::debug_stream << row << '\n';
    }*/
    
	std::cout << "[INFO] Start constructing lookup table...." << std::endl;

	size_t rowIndex = 0;
 
	for (auto & record : fin)
    {
        auto const & seq = record.sequence();
        std::string sequence;
        for (auto const & c : seq) sequence += seqan3::to_char(c);
        
        if(sequence.length() < qMerSize){

			++rowIndex;
			continue;

		}

        uint64_t content = computeContent(sequence, threshold, qMerSize);
        for (size_t i = 0; i < resultsVector[rowIndex].size(); ++i)
        {
            if (resultsVector[rowIndex][i] >= content)
            {
                resultsVector[rowIndex][i] = 1;
            }
            else
            {
                resultsVector[rowIndex][i] = 0;
            }
        }
        ++rowIndex;
    }
    
    /*
	for (auto const &row : resultsVector)
    {
        seqan3::debug_stream << row << '\n';
    }*/

    std::cout << "[INFO] Finished constructing lookup table." << std::endl;
    std::cout << "[INFO] Start estimating peptide Counts." << std::endl;

    std::vector<double> estimations = estimation(resultsVector);

    resultsVector.clear();

    uint64_t length;
    while (lengthStream >> length) {
        referenceLengths.push_back(length);
    }

    //seqan3::debug_stream << referenceLengths << std::endl;

    /*
    seqan3::debug_stream << "Estimations: ";
    for (auto const &estimation : estimations)
    {
        seqan3::debug_stream << estimation << " ";
    }
    seqan3::debug_stream << '\n';*/

    std::vector<std::pair<double,std::string>> results;

    if (estimations.size() == referenceMap.size()) {
        for (size_t i = 0; i < estimations.size(); ++i) {
            //results.push_back(std::make_pair(estimations[i]/referenceLengths[i], referenceMap[i]));
            //results.push_back(std::make_pair((estimations[i]/log(referenceLengths[i])), referenceMap[i]));
            results.push_back(std::make_pair(estimations[i], referenceMap[i]));
        }
    } else {
        std::cerr << "[ERROR] Counting and mapping results should have same size!" << std::endl;
    } 

    std::ofstream outputResults(evaluationResultsFile);

     if (outputResults.is_open()) {
            for (const auto& resultPair : results) {
                outputResults << resultPair.first << ' ' << resultPair.second << '\n';
            }

            outputResults.close();
            std::cout << "[INFO] Finished writing the results to: "<< evaluationResultsFile << std::endl;
    }
    else {
            std::cerr << "[ERROR] Error opening the file: "<< evaluationResultsFile << std::endl;
    }

}