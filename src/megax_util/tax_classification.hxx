#include <vector>
#include <map>
#include <regex>
#include <fstream>



/*
 * @fn taxaClassification
 * @brief tax.classification of target strains.
 * @signature void taxaClassification(const std::string& countingResultsFile, const std::string& taxClassificationFile);
 * @param countingResultsFile path to counting results file.
 * @param taxClassificationFile path to output classification file.
 * @return None.
 * @throws std::run_time.
 */
void taxaClassification(const std::string& countingResultsFile, const std::string& taxClassificationFile) {
    
    std::ifstream countingResultsFileStream(countingResultsFile);
    if (!countingResultsFileStream.is_open()) {
        throw std::runtime_error("[ERROR] Error opening counting results file!");
    }

    std::vector<std::pair<double, std::string>> countingResults;
    double maximumLikelihood;
    std::string proteinName;

    std::cout << "[INFO] Start parsing counting results...." << std::endl;

    while (countingResultsFileStream >> maximumLikelihood >> std::ws && std::getline(countingResultsFileStream, proteinName)) {
        countingResults.emplace_back(maximumLikelihood, proteinName);
    }

    countingResultsFileStream.close();

    std::cout << "[INFO] Finished parsing counting results." << std::endl;

    std::map<std::string, double> tmp;

    std::cout << "[INFO] Start data normalization..." << std::endl;

    // Multiply likelihood values for each strain
    for (auto& classificationLevel : countingResults) {
        std::smatch match;
        if (std::regex_search(classificationLevel.second, match, std::regex("\\[([^\\]]+)\\]"))) {
            // Extract the strain name from regex
            const std::string strain = match[1].str();

            // If the classification on the protein level is 0 then replace with 1 to avoid removing the whole strain
            if (classificationLevel.first == 0.0) {
                //classificationLevel.first += 1.0;
                continue; // skip the classification level for correct normalization 
            }
            
            // Initialize with 1.0 if not already present
            if (tmp.find(strain) == tmp.end()) {
                tmp[strain] = 1.0;
            }
            tmp[strain] *= classificationLevel.first;
        }
    }

    countingResults.clear();

    // Normalization sum
    double sum = 0.0;
    for (const auto& class_ : tmp) {
        sum += class_.second;
    }
    // Normalize the values
    std::map<std::string, double> result;
    for (const auto& class_ : tmp) {
        result[class_.first] = class_.second / sum;
    }

    std::vector<std::pair<std::string, double>> taxClassificationResults(result.begin(), result.end());
    std::sort(taxClassificationResults.begin(), taxClassificationResults.end(), [](const std::pair<std::string, double>& lCmp, const std::pair<std::string, double>& rCmp) {
        return lCmp.second > rCmp.second; 
    });

    result.clear();

    std::cout << "[INFO] Finished data normalization." << std::endl;
    std::cout << "[INFO] Start writing the classification results to the output file: " << taxClassificationFile << std::endl;


    std::ofstream output(taxClassificationFile); 

    if (!output.is_open()) {
        throw std::runtime_error("[ERROR] Error opening results file!");
    }

    for (const auto& tax : taxClassificationResults) {
        output << "[" << tax.first << "] " << tax.second << '\n';
    }

    std::cout << "[INFO] Finished writing the classification results to the output file: " << taxClassificationFile << std::endl;
}

/*
int main() {
    std::vector<std::pair<double, std::string>> test = {
        {0.2, "YP_010771932.1 MAG: hypothetical protein QIT33_gp32 [unclassified virus PBV305]"},
        {0.2, "YP_010771933.1 MAG: hypothetical protein QIT33_gp33 [Salmonella phage SP069]"},
        {0.15, "YP_010771941.1 MAG: transcription initiation factor [Salmonella phage L13]"},
        {0.3, "YP_010771949.1 MAG: hypothetical protein QIT33_gp49 [Salmonella phage SP069]"},
        {0.25, "YP_010772626.1 MAG: hypothetical protein QIT47_gp08 [Salmonella phage L13]"},
        {0.1, "YP_010772622.1 MAG: hypothetical protein QIT47_gp04 [unclassified virus PBV305]"},
        {0.0, "YP_009617910.1 hypothetical protein QII00_sBgp17 [unclassified virus PBV305]"},
        {0, "YP_008058256.1 tail assembly chaperone [Salmonella phage L13]"},
        {0, "YP_010773461.1 hypothetical protein QIT85_gp71 [Salmonella phage L13]"},
        {0.1, "YP_010773455.1 recombination protein [Salmonella phage L13]]"},
        {0.05, "YP_010768515.1 Phage tail protein (Tail_P2_I) [unclassified virus PBV305]"},
    };

    std::map<std::string, double> result = taxaClassification(test);
    
    std::vector<std::pair<std::string, double>> taxClassificationResults(result.begin(), result.end());
    std::sort(taxClassificationResults.begin(), taxClassificationResults.end(), [](const std::pair<std::string, double>& lCmp, const std::pair<std::string, double>& rCmp) {
        return lCmp.second > rCmp.second; 
    });

    for (const auto& entry : taxClassificationResults) {
        std::cout << "[" << entry.first << "]: " << entry.second << '\n';
    }

    return 0;
}
*/