#include "CascadiaMzmlReader.h"
#include "CascadiaOutput.h"
#include "CascadiaPeptideTokenizer.h"
#include <torch/torch.h>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  try {
    auto tokenizer = CascadiaPeptideTokenizer::from_bundle("#cascadia-bundle-v1\n$\t0\nA\t71.037113805\n");
    if (tokenizer.stop_token_id() != 1 || tokenizer.vocab_size_with_padding() != 3) throw std::runtime_error("vocabulary order changed");
    auto ids = torch::tensor({{2, 1, 2}}, torch::kInt64);
    if (tokenizer.detokenize(ids)[0] != "A") throw std::runtime_error("stop token ignored");
    if (std::abs(tokenizer.peptide_mz(ids[0], 2) - 45.53111571938) > 1e-6) throw std::runtime_error("incorrect peptide mass");
    bool rejected = false;
    try { CascadiaPeptideTokenizer::from_bundle(""); } catch (const std::exception&) { rejected = true; }
    if (!rejected) throw std::runtime_error("missing metadata accepted");
    auto details = cascadia_peptide_details("[Acetyl]-AM[Oxidation]C[Carbamidomethyl]D-[Amidated]");
    if (details.sequence != "AMCD" || !details.valid || details.modifications != "Acetyl@N-term, Oxidation@2, Carbamidomethyl@3, Amidated@C-term") throw std::runtime_error("Modification locations changed");
    auto header = cascadia_fasta_header({"demo.mzML", "0", "2", "AM", "UNKNOWN", "0.9", "1", "0", "2"}, details.modifications);
    if (header.find("file=demo.mzML_scan=0_charge=2_score-type=UNKNOWN_score=0.9_retention-time=1_start-time=0_end-time=2 [modification:") != 0) throw std::runtime_error("FASTA header fields lost");
    auto modified = CascadiaPeptideTokenizer::from_bundle("#cascadia-bundle-v1\n$\t0\nM\t131.040484645\nM[Oxidation]\t147.035399275\n[+25.980265]-\t25.980265\n");
    if (modified.ssl_sequence(torch::tensor({4,3,1}, torch::kInt64)) != "[+25.980265]-M[+15.99]") throw std::runtime_error("Python modification formatting changed");
    CascadiaMzmlReader fixture_reader;
    auto fixture = fixture_reader.read(std::filesystem::path(CASCADIA_TEST_DATA_DIR) / "tiny.mzML");
    if (fixture.size() != 3 || std::abs(fixture[0].retention_time_seconds - 60.0) > 1e-6 || std::abs(fixture[1].retention_time_seconds - 60.2) > 1e-6) throw std::runtime_error("mzML scan time units changed");
    if (fixture[0].mz[1] != 509.5 || fixture[0].intensity[1] != 25 || fixture[1].mz[0] != 100.5 || fixture[1].intensity[1] != 16) throw std::runtime_error("mzML float/zlib decoding changed");
    auto augmented_fixture = fixture_reader.build_augmented_spectra(fixture, {});
    auto fixture_batch = fixture_reader.to_tensors(augmented_fixture);
    if (augmented_fixture.size() != 3 || fixture_batch.spectra.size(1) != 4 || augmented_fixture.front().charge != 2 || augmented_fixture.back().charge != 4) throw std::runtime_error("DIA candidate augmentation changed");
    if (std::abs(fixture_batch.spectra.index({0, torch::indexing::Slice(), 1}).sum().item<float>() - 1.0F) > 1e-6) throw std::runtime_error("Combined intensities are not normalized");
    if (std::abs(fixture_batch.precursors[0][0].item<float>() - 1993.95634213248) > 1e-3) throw std::runtime_error("Legacy precursor transform changed");
    if (argc == 3) {
      CascadiaMzmlReader reader;
      auto raw = reader.read(argv[1]);
      auto augmented = reader.build_augmented_spectra(raw, {});
      auto batch = reader.to_tensors(augmented);
      torch::serialize::OutputArchive archive;
      archive.write("spectra", batch.spectra);
      archive.write("precursors", batch.precursors);
      archive.save_to(argv[2]);
      std::cout << "Prepared " << augmented.size() << " candidates, " << batch.spectra.sizes() << '\n';
    }
    std::cout << "Native vocabulary/mass validation passed\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
