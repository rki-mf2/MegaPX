#include <string>
#include <vector>
#include <random>

/*
 * @fn simulatePeptide
 * @brief Simulate peptides of the protein sequence.
 * @signature std::string simulatePeptide(const std::string& sequence, uint8_t errors, uint8_t seed, uint8_t rangeBegin, uint8_t rangeEnd);
 * @param sequence input protein sequence.
 * @param errors number of errors [default: 0].
 * @param seed seeds value [default: 0], the default value produces new sequences for each run.
 * @param rangeBegin starting length [default: 5].
 * @param rangeEnd ending length [default: 50].
 * @return Simulated peptide sequence.
 * @throws None.
 */
std::string simulatePeptide(const std::string& sequence, uint8_t errors, uint8_t seed, uint8_t rangeBegin, uint8_t rangeEnd) {
    std::random_device randomObject;

    std::mt19937 generateRandomNumber;
    if (seed == 0) {
        generateRandomNumber = std::mt19937(randomObject());
    } else {
        generateRandomNumber = std::mt19937(seed);
    }

    std::uniform_int_distribution<int> lengthValue(rangeBegin, rangeEnd);
    std::uniform_int_distribution<int> startPosition(0, sequence.length() - rangeBegin);

    size_t length = lengthValue(generateRandomNumber);
    
    size_t start = startPosition(generateRandomNumber);

    std::string subPeptideSeq = sequence.substr(start, length);

    const std::string aminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    std::uniform_int_distribution<int> errorPosition(0, subPeptideSeq.length() - 1);

    std::string mutatedPeptide(subPeptideSeq);

    if (mutatedPeptide.size() < 5){

            std::cout<< mutatedPeptide.size() << std::endl;
        }

    for (size_t i = 0; i < errors; ++i) {
        size_t index = errorPosition(generateRandomNumber);
        char newAminoAcid = aminoAcids[generateRandomNumber() % aminoAcids.length()];

        while (newAminoAcid == subPeptideSeq[index]) {
            newAminoAcid = aminoAcids[generateRandomNumber() % aminoAcids.length()];
        }

        mutatedPeptide[index] = newAminoAcid;
    }

    return mutatedPeptide;
}

/*
 * @fn simulate
 * @brief Simulate peptides of the protein sequence.
 * @signature void simulate(unsigned numberOfSequences, uint8_t errors, uint8_t seed, uint8_t rangeBegin, uint8_t rangeEnd, std::string& sequenceFilePath, std::string& outputDir);
 * @param numberOfSequences number of target sequences for the simulation.
 * @param errors number of errors [default: 0].
 * @param seed seeds value [default: 0], the default value produces new sequences for each run.
 * @param rangeBegin starting length [default: 5].
 * @param rangeEnd ending length [default: 50].
 * @param sequenceFilePath path to input sequence file.
 * @param outputDir output dir for the target peptide sequences
 * @return None.
 * @throws None.
 */
void simulate(uint64_t numberOfSequences, uint8_t errors, uint8_t seed, uint8_t rangeBegin, uint8_t rangeEnd, std::string& sequenceFilePath, std::string& outputDir){

    std::filesystem::path fileName{sequenceFilePath};

    const std::string outputDB = outputDir + "simulated_" + fileName.filename().string();
    seqan3::sequence_file_output fout{outputDB};

    std::mutex uniqueSetMutex;
    std::unordered_set<std::string> uniqueSet;

    using types = seqan3::type_list<std::vector<seqan3::aa27>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;
    using traits_type = seqan3::sequence_file_input_default_traits_aa;

    seqan3::sequence_file_input<traits_type> fin{sequenceFilePath};
    
    seqan3::debug_stream << "[LOG] Reference File Name        : "<< sequenceFilePath << '\n';
    seqan3::debug_stream << "[LOG] Output File Name           : "<< outputDB << '\n';
    seqan3::debug_stream << "[LOG] Output Dir                 : "<< outputDir << '\n';
    seqan3::debug_stream << "[LOG] Number of Sequences        : "<< numberOfSequences << '\n';
    seqan3::debug_stream << "[LOG] Number of errors           : "<< unsigned(errors) << '\n';
    seqan3::debug_stream << "[LOG] Seed                       : "<< unsigned(seed) << '\n';
    seqan3::debug_stream << "[LOG] Range Begin                : "<< unsigned(rangeBegin) << '\n';
    seqan3::debug_stream << "[LOG] Range End                  : "<< unsigned(rangeEnd) << '\n';

    unsigned seqCounter {0u};
    unsigned seqCounterOr {0u};

    for (auto & record : fin)
    {
        seqCounterOr++;

        unsigned counter {0u};

        auto const& seq = record.sequence();
        auto const& header = record.id();

        std::string sequence = "";
        for (auto const & c : seq) sequence += seqan3::to_char(c); 

        for (size_t i = 0; i < numberOfSequences; ++i) {

            std::string newPeptide = simulatePeptide(sequence, errors, seed, rangeBegin, rangeEnd);
            const std::string newHeader = header + " (sim_" + std::to_string(i) + ")";

            seqan3::aa27_vector newSequence{};
            for (const char& c : newPeptide){

                seqan3::aa27 aminoAcid = {};
                aminoAcid.assign_char(c);
                newSequence.push_back(aminoAcid);
            } 

            sequence_record_type record{std::move(newSequence), std::move(newHeader)};
            seqCounter++;
            fout.push_back(record);
        }

    }

    seqan3::debug_stream << "[INFO] Number of target sequences   : " << seqCounterOr << '\n';
    seqan3::debug_stream << "[INFO] Number of simulated sequences: " << seqCounter << '\n';

}

