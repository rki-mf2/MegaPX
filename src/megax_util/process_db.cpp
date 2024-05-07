#include <iostream>
#include <fstream>
#include <regex>
#include <map>



void sequenceConc(const std::string& referenceFilePath, const std::string& outputDir) {
    std::map<std::string, std::string> sequencesDict;
    std::string currentSpecies = "";
    std::string currentSequence = "";
    int count = 0;

    std::ifstream fastaFile(referenceFilePath);
    if (!fastaFile.is_open()) {
        std::cerr << "Error opening file: " << referenceFilePath << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fastaFile, line)) {
        if (line[0] == '>') {
            count++;
            //std::cout << "Processing protein number: " << count << std::endl;
            std::smatch match;
            std::regex_search(line, match, std::regex("\\[([^\\[\\]]+(?:\\[[^\\]]+\\])?)\\]"));
            
            if (!currentSequence.empty()) {
                if (sequencesDict.find(currentSpecies) != sequencesDict.end()) {
                    currentSequence = '*' + currentSequence;
                    sequencesDict[currentSpecies] += currentSequence;
                } else {
                    sequencesDict[currentSpecies] = currentSequence;
                }
            }
            
            currentSequence = "";
            
            if (match.size() > 1) {
                currentSpecies = '[' + match[1].str() + ']';
            }
        } else {
            currentSequence += line;
        }
    }

    fastaFile.close();

    size_t lastSlash = referenceFilePath.find_last_of("/\\");
    std::string path = referenceFilePath.substr(0, lastSlash);
    std::string DBFileName = referenceFilePath.substr(lastSlash + 1);
    size_t lastDot = DBFileName.find_last_of(".");
    std::string file_ = DBFileName.substr(0, lastDot);

    std::string outputDB = outputDir + file_ + "_tmp_cpp.fasta";
    std::string outputHist = outputDir + file_ + "_sequence_distribution_cpp.hist";

    std::ofstream DB(outputDB);
    if (!DB.is_open()) {
        std::cerr << "Error creating output file: " << outputDB << std::endl;
        return;
    }

    for (const auto& entry : sequencesDict) {
        DB << '>' << entry.first << '\n';
        DB << entry.second << '\n';
    }

    DB.close();

    std::cout << "Start constructing protein distribution data: " << std::endl;
    std::ofstream histFile(outputHist);
    if (!histFile.is_open()) {
        std::cerr << "Error creating output file: " << outputHist << std::endl;
        return;
    }

    for (const auto& entry : sequencesDict) {
        size_t seqLength = std::count(entry.second.begin(), entry.second.end(), '*');
        histFile << entry.first << ": " << (seqLength + 1) << " (Proteins) | "
                  << (entry.second.length() - seqLength) << " (AA)\n";
    }

    histFile.close();
}

int main() {
    std::string db_file = "/home/lutfia/scratch/megaX/build/main/Monkeypox_Chlorocebus_2022.fasta";
    //std::string db_file = "refSeqViral.fasta";
    std::string output_path = "./";

    size_t lastSlash = db_file.find_last_of("/\\");
    std::string path = db_file.substr(0, lastSlash);
    std::string db_file_name = db_file.substr(lastSlash + 1);
    size_t lastDot = db_file_name.find_last_of(".");
    std::string file = db_file_name.substr(0, lastDot);
    std::string extension = db_file_name.substr(lastDot);
    std::cout << "Start conatenating sequences" << std::endl;
    sequenceConc(db_file, output_path);
    std::cout << "Finish conatenating sequences" << std::endl;
    std::cout << file << std::endl;
    std::string tmp_db = output_path + file + "_tmp_cpp" + extension;
    std::string output_db = output_path + file + "_concatenated_cpp" + extension;
    std::cout << "Start reordering sequences" << std::endl;
    std::ifstream input_file(tmp_db);
    std::ofstream output_file(output_db);
    std::cout << tmp_db << '\n';
    std::cout << output_db << '\n';
    if (input_file.is_open() && output_file.is_open()) {
        std::string line;
        while (std::getline(input_file, line)) {
            if (line[0] == '>') {
                output_file << '\n' << line << '\n';
            } else {
                output_file << line << '\n';
            }
        }
        std::cout << "Output DB: " << output_db << std::endl;
        input_file.close();
        output_file.close();
        std::remove(tmp_db.c_str());
    } else {
        std::cerr << "Error opening files for final processing." << std::endl;
    }

    return 0; // sed -i '/^$/d' : remove empty lines 
    // count number of > in file: tr -cd '>' < U-RVDBv26.0-prot.fasta| wc -c
}
