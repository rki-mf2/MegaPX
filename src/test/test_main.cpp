#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include "../mutate/mutate.hpp"
#include "../megax_util/evaluation.hpp"
#include "../megax_util/peptide_simulation.hxx"
#include "../megax_util/tax_classification.hxx"

namespace AminoAcidsProcessing
{
char convertIndexToAminoAcid(unsigned i);
unsigned mapCharToIndex(char aminoAcid);
}

std::vector<std::string> qMerizeSequence(std::string const & sequence, int qMerSize);
std::vector<int> mapSequenceToIndices(std::string const & proteinSequence);

using hashesType = std::vector<uint64_t>;
hashesType computeMinimiserSingleBin(uint8_t kMerSize, uint8_t windowSize, hashesType hashes);

namespace
{

std::filesystem::path temp_root()
{
    auto root = std::filesystem::temp_directory_path() / "megapx_tests";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    return root;
}

void require(bool condition, std::string const & message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template <typename left_t, typename right_t>
void require_equal(left_t const & actual, right_t const & expected, std::string const & message)
{
    if (!(actual == expected))
    {
        std::ostringstream stream;
        stream << message << " (actual=" << actual << ", expected=" << expected << ")";
        throw std::runtime_error(stream.str());
    }
}

std::vector<std::string> split_dollar(std::string const & value)
{
    std::vector<std::string> parts;
    std::stringstream stream{value};
    std::string item;

    while (std::getline(stream, item, '$'))
    {
        if (!item.empty())
            parts.push_back(item);
    }

    return parts;
}

seqan3::aa27_vector to_aa27(std::string const & value)
{
    seqan3::aa27_vector result{};
    for (char c : value)
    {
        seqan3::aa27 aa{};
        aa.assign_char(c);
        result.push_back(aa);
    }
    return result;
}

std::size_t hamming_distance(std::string const & lhs, std::string const & rhs)
{
    require_equal(lhs.size(), rhs.size(), "Hamming distance requires equal-length strings");
    std::size_t distance{0};
    for (std::size_t i = 0; i < lhs.size(); ++i)
        distance += lhs[i] != rhs[i];
    return distance;
}

std::size_t count_fasta_records(std::filesystem::path const & file)
{
    std::ifstream input{file};
    require(input.is_open(), "Failed to open fasta file for counting");

    std::size_t count{0};
    std::string line;
    while (std::getline(input, line))
        count += !line.empty() && line[0] == '>';
    return count;
}

std::vector<std::string> read_lines(std::filesystem::path const & file)
{
    std::ifstream input{file};
    require(input.is_open(), "Failed to open file");

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
        lines.push_back(line);
    return lines;
}

void test_sequence_helpers()
{
    require_equal(AminoAcidsProcessing::convertIndexToAminoAcid(0u), 'A', "Index 0 should map to A");
    require_equal(AminoAcidsProcessing::convertIndexToAminoAcid(19u), 'Y', "Index 19 should map to Y");
    require_equal(AminoAcidsProcessing::convertIndexToAminoAcid(20u), 'X', "Out-of-range index should map to X");

    require_equal(AminoAcidsProcessing::mapCharToIndex('A'), 0u, "A should map to 0");
    require_equal(AminoAcidsProcessing::mapCharToIndex('Y'), 19u, "Y should map to 19");
    require_equal(AminoAcidsProcessing::mapCharToIndex('a'), 0u, "Lowercase amino acids should map");

    auto qmers = qMerizeSequence("ABCDE", 3);
    require_equal(qmers.size(), std::size_t{3}, "Expected three 3-mers");
    require_equal(qmers[0], std::string{"ABC"}, "Unexpected first q-mer");
    require_equal(qmers[2], std::string{"CDE"}, "Unexpected last q-mer");
    require(qMerizeSequence("ABC", 5).empty(), "q-mers should be empty when q > sequence length");

    auto indices = mapSequenceToIndices("ACDYZ");
    require_equal(indices.size(), std::size_t{5}, "Sequence should map element-wise");
    require_equal(indices[0], 0, "A should map to 0");
    require_equal(indices[1], 1, "C should map to 1");
    require_equal(indices[2], 2, "D should map to 2");
    require_equal(indices[3], 19, "Y should map to 19");
    require_equal(indices[4], -1, "Unknown amino acids should map to -1");
}

void test_minimiser_and_evaluation_helpers(std::filesystem::path const & root)
{
    hashesType hashes{7, 197, 5330, 6142, 8382, 9825};
    auto minimisers = computeMinimiserSingleBin(3, 5, hashes);
    require_equal(minimisers.size(), std::size_t{4}, "Unexpected minimiser count");
    require_equal(minimisers[0], uint64_t{7}, "Unexpected first minimiser");
    require_equal(minimisers[3], uint64_t{6142}, "Unexpected last minimiser");

    require_equal(computeContent("ACDEF", 0.75, 2), uint64_t{3}, "Unexpected threshold content");
    require_equal(computeContent("ACDEFG", 0.50, 3), uint64_t{2}, "Unexpected threshold content");

    assignedPeptidesFile = (root / "assigned_peptides_idx.log").string();
    std::vector<std::vector<uint8_t>> lookup_table{{1, 0, 1}, {1, 1, 0}};
    auto estimations = estimation(lookup_table);
    require_equal(estimations.size(), std::size_t{3}, "Unexpected estimation size");
    require(std::abs(estimations[0] - 1.0) < 1e-9, "Unexpected estimation for bin 0");
    require(std::abs(estimations[1] - 0.5) < 1e-9, "Unexpected estimation for bin 1");
    require(std::filesystem::exists(root / "assigned_peptides_idx.log"), "Assigned peptides file was not written");
}

void test_results_mapping(std::filesystem::path const & root)
{
    auto mapping_file = root / "reference_id_map.log";
    auto lengths_file = root / "lengths.log";
    auto results_file = root / "counts.bin";
    auto query_file = root / "queries.fasta";
    auto output_file = root / "evaluation.log";

    {
        std::ofstream mapping{mapping_file};
        mapping << "0 Reference Alpha\n";
        mapping << "1 Reference Beta\n";
    }

    {
        std::ofstream lengths{lengths_file};
        lengths << "100\n";
        lengths << "120\n";
    }

    {
        std::ofstream output{results_file, std::ios::binary};
        cereal::BinaryOutputArchive archive{output};
        std::vector<std::vector<uint8_t>> counts{{3, 0}, {1, 3}};
        archive(counts);
    }

    {
        std::ofstream query{query_file};
        query << ">q1\nABCDE\n";
        query << ">q2\nFGHIK\n";
    }

    resultsMapping(results_file.string(),
                   mapping_file.string(),
                   output_file.string(),
                   query_file.string(),
                   0.75,
                   2);

    auto lines = read_lines(output_file);
    require_equal(lines.size(), std::size_t{2}, "Expected one mapped result per reference");
    require(lines[0].find("0.5 Reference Alpha") == 0, "Unexpected first evaluation line");
    require(lines[1].find("0.5 Reference Beta") == 0, "Unexpected second evaluation line");
    require(std::filesystem::exists(root / "assigned_peptides_idx.log"), "resultsMapping should write assigned peptide stats");
}

void test_peptide_simulation(std::filesystem::path const & root)
{
    std::string protein{"ACDEFGHIKLMNPQRSTVWY"};
    auto original = simulatePeptide(protein, 0, 42, 5, 5);
    auto mutated = simulatePeptide(protein, 1, 42, 5, 5);

    require_equal(original.size(), std::size_t{5}, "Peptide length should match fixed range");
    require_equal(mutated.size(), std::size_t{5}, "Mutated peptide length should match fixed range");
    require(protein.find(original) != std::string::npos, "Zero-error peptide should be a substring of the input protein");
    require_equal(hamming_distance(original, mutated), std::size_t{1}, "Single-error peptide should differ by one amino acid");

    auto fasta_file = root / "protein.fasta";
    auto output_dir = (root / "simulated").string() + "/";
    std::filesystem::create_directories(output_dir);

    {
        std::ofstream fasta{fasta_file};
        fasta << ">protein_1\nACDEFGHIKLMNPQRSTVWY\n";
    }

    std::string fasta_path = fasta_file.string();
    simulate(2, 0, 42, 5, 5, fasta_path, output_dir);

    auto simulated_fasta = std::filesystem::path{output_dir} / ("simulated_" + fasta_file.filename().string());
    require(std::filesystem::exists(simulated_fasta), "Simulation output fasta was not created");
    require_equal(count_fasta_records(simulated_fasta), std::size_t{2}, "Expected two simulated records");
}

void test_tax_classification(std::filesystem::path const & root)
{
    auto counting_file = root / "counting_results.log";
    auto output_file = root / "tax_classification.log";

    {
        std::ofstream output{counting_file};
        output << "0.5 protein A [Strain Alpha]\n";
        output << "0.5 protein B [Strain Alpha]\n";
        output << "0.2 protein C [Strain Beta]\n";
        output << "0.0 protein D [Strain Beta]\n";
    }

    taxaClassification(counting_file.string(), output_file.string());

    auto lines = read_lines(output_file);
    require_equal(lines.size(), std::size_t{2}, "Expected one normalized line per strain");
    require(lines[0].find("[Strain Alpha]") == 0, "Highest-probability strain should be listed first");
    require(lines[1].find("[Strain Beta]") == 0, "Second strain should be listed second");
}

void test_mutate_behaviour()
{
    std::filesystem::path matrix_file = std::filesystem::path{MEGAPX_REPO_DIR} / "data/matrix/blosum62";
    Mutate mutate{matrix_file.string()};

    int aa_a = static_cast<int>(AminoAcidsProcessing::mapCharToIndex('A'));
    int aa_r = static_cast<int>(AminoAcidsProcessing::mapCharToIndex('R'));
    int aa_s = static_cast<int>(AminoAcidsProcessing::mapCharToIndex('S'));

    require_equal(mutate.scoreGetter(aa_a, aa_a), 4, "A/A score should match BLOSUM62");
    require_equal(mutate.scoreGetter(aa_a, aa_r), -1, "A/R score should match BLOSUM62");
    require_equal(mutate.scoreGetter(aa_a, aa_s), 1, "A/S score should match BLOSUM62");
    require_equal(mutate.calculateQmerScore("AA", "AR"), 3, "Unexpected q-mer score");

    auto single = mutate.scoring("A", 1);
    std::set<std::string> single_neighbors;
    for (auto const & neighbor : single.neighbors)
        single_neighbors.insert(neighbor.first);
    require_equal(single_neighbors.size(), std::size_t{2}, "Expected A and S neighbors for threshold 1");
    require(single_neighbors.count("A") == 1, "A should be retained as a neighbor");
    require(single_neighbors.count("S") == 1, "S should be retained as a neighbor");

    auto recursive = mutate.scoring("AA", 8);
    require_equal(recursive.neighbors.size(), std::size_t{1}, "Only AA should satisfy score 8");
    require_equal(recursive.neighbors[0].first, std::string{"AA"}, "Unexpected recursive scoring result");
    require_equal(recursive.neighbors[0].second, 8, "Unexpected recursive scoring score");

    auto mutated_qmers = split_dollar(mutate.mutateSingleSequence(1, "A", 1, 1));
    std::set<std::string> mutated_set{mutated_qmers.begin(), mutated_qmers.end()};
    require_equal(mutated_set.size(), std::size_t{2}, "Expected unique mutated q-mers");
    require(mutated_set.count("A") == 1, "Original q-mer should be present");
    require(mutated_set.count("S") == 1, "Mutated q-mer should be present");

    auto hashed = mutate.mutateSingleSequenceHashes(1, "A", 1, 1, false, 0);
    std::vector<uint64_t> expected{
        mutate.hashQMer(to_aa27("A"), 1),
        mutate.hashQMer(to_aa27("S"), 1)
    };
    std::sort(expected.begin(), expected.end());
    require_equal(hashed.size(), expected.size(), "Unexpected hash count");
    require(std::equal(hashed.begin(), hashed.end(), expected.begin()), "Unexpected mutated q-mer hashes");
}

void test_mutate_statistics(std::filesystem::path const & root)
{
    std::filesystem::path matrix_file = std::filesystem::path{MEGAPX_REPO_DIR} / "data/matrix/blosum62";
    Mutate mutate{matrix_file.string()};

    auto fasta_file = root / "references.fasta";
    auto output_dir = (root / "stats").string() + "/";
    std::filesystem::create_directories(output_dir);

    {
        std::ofstream fasta{fasta_file};
        fasta << ">ref_a\nACDEF\n";
        fasta << ">ref_b\nLMNPQ\n";
    }

    mutate.printSequenceStatstics(fasta_file.string(), 3, output_dir);

    auto mapping_lines = read_lines(std::filesystem::path{output_dir} / "reference_id_map.log");
    auto lengths_lines = read_lines(std::filesystem::path{output_dir} / "lengths.log");
    auto histogram_lines = read_lines(std::filesystem::path{output_dir} / "reference_hist.log");

    require_equal(mapping_lines.size(), std::size_t{2}, "Expected mapping file to contain both references");
    require_equal(lengths_lines.size(), std::size_t{2}, "Expected length file to contain both references");
    require_equal(histogram_lines.size(), std::size_t{2}, "Expected histogram file to contain both references");
}

} // namespace

int main()
{
    auto root = temp_root();

    struct named_test
    {
        std::string name;
        std::function<void()> run;
    };

    std::vector<named_test> tests{
        {"sequence_helpers", [] { test_sequence_helpers(); }},
        {"minimiser_and_evaluation_helpers", [&] { test_minimiser_and_evaluation_helpers(root); }},
        {"results_mapping", [&] { test_results_mapping(root); }},
        {"peptide_simulation", [&] { test_peptide_simulation(root); }},
        {"tax_classification", [&] { test_tax_classification(root); }},
        {"mutate_behaviour", [] { test_mutate_behaviour(); }},
        {"mutate_statistics", [&] { test_mutate_statistics(root); }}
    };

    std::size_t failed{0};

    for (auto const & test : tests)
    {
        try
        {
            test.run();
            std::cout << "[PASS] " << test.name << '\n';
        }
        catch (std::exception const & ex)
        {
            ++failed;
            std::cerr << "[FAIL] " << test.name << ": " << ex.what() << '\n';
        }
    }

    std::filesystem::remove_all(root);
    return failed == 0 ? 0 : 1;
}
