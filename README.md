# MegaPX: Fast Peptide Assignment Method Using Index Databases with Various Counting Algorithms

<p align="center">
<img src="img/logo_1.png" alt="MegaPX logo" width="40%">
</p>
<div align="center">

[![bioRxiv](https://img.shields.io/badge/biorxiv-2025.04.14.648734v1-b31b1b)](https://www.biorxiv.org/content/10.1101/2025.04.14.648734v1)&nbsp;
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-10.1093%2Fbioinformatics%2Fbtag134-blue)](https://doi.org/10.1093/bioinformatics/btag134)
[![code](https://img.shields.io/badge/code-MegaPX-blue)](https://github.com/lutfia95/MegaPX)&nbsp;

[![CMake tests](https://github.com/lutfia95/MegaPX/actions/workflows/cmake-single-platform.yml/badge.svg?branch=master)](https://github.com/lutfia95/MegaPX/actions/workflows/cmake-single-platform.yml)
[![Sequencer and sanitizers CI](https://github.com/lutfia95/MegaPX/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/lutfia95/MegaPX/actions/workflows/ci.yml)
[![Doxygen docs](https://github.com/lutfia95/MegaPX/actions/workflows/docs.yml/badge.svg?branch=master)](https://github.com/lutfia95/MegaPX/actions/workflows/docs.yml)
![C++](https://img.shields.io/badge/C%2B%2B-23%20%7C%2017-00599C?logo=cplusplus&logoColor=white)
![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.20-064F8C?logo=cmake&logoColor=white)
![OpenMP](https://img.shields.io/badge/OpenMP-search%20required-5C2D91)
![LibTorch](https://img.shields.io/badge/LibTorch-optional%20sequencing-EE4C2C?logo=pytorch&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?logo=python&logoColor=white)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-FCC624?logo=linux&logoColor=black)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)


</div>



<!-- ## Citation  -->

## Table of contents

* [Description](#description)
* [Installation](#installation)
* [Tests and CI](#tests)
* [Cascadia C++ loader](#cascadia-mzml-sequencing)
* [All loader settings](src_libtorch/README.md#complete-command-line-reference)
* [Python/C++ agreement benchmarks](src_libtorch/README.md#reproduce-all-three-pythonc-benchmarks)
* [Commands](#commands)
* [Use-cases](#use-cases)
* [Parameters](#parameters )

## <a name="description"></a>Description
MegaPX is designed to mutate protein databases, enabling users to search for de novo peptides in both the original and mutated versions. The tool employs a _k_-mer-based search, providing rapid and highly accurate peptide assignment capability. Beyond this primary function, the tool has numerous other use cases, including metaproteomics classification, virus assignment, and peptide classification.   
MegaPX utilizes four distinct counting algorithms: [IBF](https://github.com/seqan/hibf/), [HBF](https://github.com/seqan/hibf/), [FM-Index](https://github.com/seqan/seqan3), and binary search. As we assessed these four algorithms in terms of runtime and performance, it became imperative to incorporate them into a single tool, allowing users the freedom to choose based on their specific use case. Mutation generation relies on a substitution matrix containing values proportional to the likelihood of amino acid `i` mutating into amino acid `j` across all possible pairs of amino acids. Each reference amino acid sequence undergoes _k_-merization, and mutated _k_-mers are generated based on the user-defined score, which is the sum of _k_-mer alignment scores. The tool is designed in step commands, which will be included as rules in the Snakemake pipeline. 
## <a name="installation"></a>Installation   
The main MegaPX search tool targets Linux x86_64 and has also been tested under WSL 2. The standalone Cascadia sequencer additionally supports macOS; see its build instructions below. We recommend Conda for the search dependencies:
```
conda create --name MegaPX_env OR  conda create --prefix path/to/MegaPX_env  
conda activate MegaPX_env  OR conda activate path/to/MegaPX_env     
conda install -c conda-forge cxx-compiler   
conda install conda-forge::gcc  
conda install conda-forge::clang
conda install -c conda-forge zlib
```
For optional Cascadia sequencing on Linux x86_64 or ARM64, CMake can download the CPU LibTorch SDK automatically. See the [sequencer build instructions and local SDK options](src_libtorch/README.md#build-the-standalone-sequencer).

Building MegaPX from source:   
```
git clone https://github.com/lutfia95/MegaPX.git   
cd MegaPX  
mkdir build && cd build   
cmake ../src   
make    
```
The executable will be built in `path/to/MegaPX/build/main/`   
Users can also download the binary pre-built version: [Linux x86_64](https://github.com/lutfia95/MegaPX/releases/download/v.0.5.4/megapx-Linux_x64.tar.gz).    
Run: `tar -xzvf megapx-Linux_x64.tar.gz`    
## <a name="tests"></a>Tests
Build the test target and run the full test suite with CTest from the build directory:
```
mkdir -p build
cd build
cmake ../src
cmake --build . --target megapx_tests
ctest --output-on-failure
```
If you want to run only the MegaPX unit test binary directly:
```
./test/megapx_tests
```

### Sequencer tests and CI

The [native unit tests](tests/cascadia_native_tests.cpp) cover vocabulary order,
stop tokens, mass calculations, SSL modification formatting, FASTA headers,
and mixed-precision/zlib mzML decoding with a [tiny fixture](tests/fixtures/tiny.mzML).
The [benchmark unit tests](tests/test_cascadia_benchmark.py) verify agreement
failure detection and independent peak-memory measurements using only Python's
standard library. The [pipeline tests](tests/test_cascadia_pipeline.py) run real
Cascadia inference and check outputs, overwrite protection and the optional
search subprocess contract; their search executable is a test stub.

```bash
# No model or Torch installation needed for these Python tests.
python tests/test_cascadia_benchmark.py

# After building the standalone native loader:
ctest --test-dir build_libtorch --output-on-failure

# Requires a local model bundle and demo mzML (kept out of Git):
python tests/test_cascadia_pipeline.py \
  --sequencer build_libtorch/megapx_sequence \
  --bundle models/cascadia.pt --mzml ../demo.mzML
```

[Sequencer CI](.github/workflows/ci.yml) runs benchmark unit tests and native
CTest with and without Address/UndefinedBehavior sanitizers. CMake downloads the CPU LibTorch SDK. The prebuilt Torch dependency is not
itself sanitizer-instrumented; leak checking of its process-lifetime caches is
disabled. Full three-checkpoint agreement/performance tests are run locally
using the [reproduction script](scripts/reproduce_benchmark.sh).

Read the [published test record and anonymized logs](docs/benchmarks/2026-11-21/README.md)
for exactly what was tested, all three checkpoint comparisons, runtime/peak-RAM
measurements, commands, environment details and local test results.

Linux sanitizer builds run AddressSanitizer and UndefinedBehaviorSanitizer.
With the prebuilt macOS Torch package, local builds use UndefinedBehaviorSanitizer
because a local ASan build stalled during startup; Linux CI is configured to cover ASan.
Configure a separate sanitizer build. On macOS, include the local Torch package path described in the [sequencer build instructions](src_libtorch/README.md#build-the-standalone-sequencer).

```bash
cmake -S src_libtorch -B build_sanitizers \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCASCADIA_ENABLE_SANITIZERS=ON
cmake --build build_sanitizers -j 2
ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build_sanitizers --output-on-failure
```

Build API docs with `doxygen Doxyfile` and open `build-docs/html/index.html`.
The [documentation workflow](.github/workflows/docs.yml) provides the generated
HTML as a GitHub Actions artifact; a Pages site is not required. Workflow
badges reflect GitHub's status once these files are pushed and Actions runs.

### Cascadia mzML sequencing

The optional LibTorch integration provides `--sequence` to sequence an mzML
file, write matched `.ssl` and peptide FASTA outputs, and optionally forward
those peptides to `multi_indexing`. The standalone sequencer can be built
without the MegaPX search dependencies. The main search tool uses C++23 on
Linux; the standalone sequencer uses C++17 and has also been tested on macOS.
OpenMP is required for search. LibTorch and zlib are required for sequencing.
See [Cascadia build, export and usage](src_libtorch/README.md), including the
[complete loader option table](src_libtorch/README.md#complete-command-line-reference)
and [benchmark/measurement commands](src_libtorch/README.md#reproduce-all-three-pythonc-benchmarks).

```bash
./build_libtorch/megapx_sequence \
  --sequence ../demo.mzML --model models/cascadia.pt --out results/demo
```

The exporter now keeps compact scripted methods by default (`--no-freeze`).
Full-file tests with all three public checkpoints matched Python's predictions;
native peak RAM fell from 3,172–3,259 MiB to 848–963 MiB on the tested Mac.
See the [measured runtime/RAM table and limitations](src_libtorch/README.md#measured-memory-and-runtime).
The [benchmark record](docs/benchmarks/2026-11-21/README.md) links the actual logs
and explains the compact and frozen exports of the same model weights.
Re-export existing frozen bundles to obtain the memory reduction.

#### Loader settings

| Option | Default | Meaning |
|---|---|---|
| `--sequence PATH` | Required | DIA mzML input; charges and scan context are expanded into candidates. |
| `--model PATH` | Required | Exported Cascadia TorchScript `.pt` bundle with its vocabulary/masses. |
| `--out DIR`, `-o DIR` | Required | Results directory; existing generated result files are protected. |
| `--help`, `-h` | — | Show command help and exit. |
| `--device auto\|cpu\|cuda` | `auto` | CUDA when available, otherwise CPU; explicit CUDA fails when unavailable. |
| `--threads N` | `4` | Positive CPU intra-op thread count. |
| `--batch-size N` | `32` | Positive number of augmented candidates decoded per batch. |
| `--score-threshold X` | `0.8` | Retain confidence strictly greater than X; allowed range is 0–1. |
| `--max-charge N` | `4` | Generate candidates with charges 2 through N; N must be at least 2. |
| `--scan-width N` | `1` | Nearest MS1 and matching MS2 scans included in augmentation. |
| `--top-peaks N` | `150` | Most intense peaks selected from each contributing scan. |
| `--max-length N` | `32` | Maximum decoding iterations; 32 reproduces Python's default. |
| `--output-mode python\|token` | `python` | Original Python scoring/filtering, or optional token-based scoring and per-window filtering. |
| `--jit-optimize on\|off` | `on` | Runtime TorchScript profiling/optimization; `off` is available for memory/performance comparison. |
| `--limit N` | All | Process only the first N augmented candidates for a smoke test. |
| `--dump-tensors PATH` | Disabled | Save the first prepared batch in a Torch archive for diagnostics. |
| `--search` | Disabled | Run MegaPX `multi_indexing` after sequencing and closing output files. |
| `--megapx PATH` | Required with search | Existing executable used for the search subprocess. |
| `--search-mode MODE` | `multi_indexing` | Currently only `multi_indexing` is supported. |
| `--` | — | Everything afterward is passed as search options; sequencing manages `-f`, `-F`, and `-o`. |

Builds, model weights, raw spectra, predictions, benchmarks and logs are excluded
by [.gitignore](.gitignore). Source, workflows, configuration examples, tests and
small reference FASTA/matrix fixtures remain trackable.

## <a name="commands"></a>Commands  
|Subcommand                                                                |Description                                                     |
|:-------------------------------------------------------------------------|:---------------------------------------------------------------|
|[**stat**](#stat)                                                         |Print reference statistics generating mapping and length file   |
|[**mutate_stat**](#mutate_stat)                                           |Print mutation statistics with number of mutated sequences      |
|[**ibf_stat**](#ibf_stat)                                                 |Print approximate IBF size with the number of technical bins    |
|[**write_db**](#write_db)                                                 |Write mutated database to output fasta file (big files)         |
|[**simulate_peptides**](#simulate_peptides)                               |Simulate peptides from given reference database and error rate  |
|[**build_vect**](#build_vect)                                             |Build first serialized index containing k_mer hashes `uint64_t` |
|[**build_hibf**](#build_hibf)                                             |Build HIBF from input index with number of included bins        |
|[**count_hibf**](#count_hibf)                                             |Counting query peptide sequences into pre-built HIBF(s)         |
|[**build_ibf**](#build_ibf)                                               |Build IBF from input index with number of included bins         |
|[**build_ibf**](#count_ibf)                                               |Counting query peptide sequences into pre-built IBF             |
|[**hibf_ref**](#hibf_ref)                                                 |Build HIBF from input database without mutation assignment      |
|[**binary_search**](#binary_search)                                       |Run binary search based classification on the serialized index  |
|[**build_fm**](#build_fm)                                                 |Build FM index from target reference database                   |
|[**count_fm**](#count_fm)                                                 |Count query peptide sequences into pre-built FM index           |
|[**evaluate**](#evaluate)                                                 |Run results evaluation and generate assignment report           |
|[**classification**](#classification)                                     |Normalize sequences on species level (used for refSeqViral)     |
|[**multi_indexing**](#multi_indexing)                                     |Run multi-indexing search, builds set of IBFs for the input     |

### Mulit-Indexing
If running the easy-assignment of peptides without any other algorithm, we recommend using `multi_indexing`
User parameters: 
```
-m Path to input matrix.
-b Path to blacklist file.
-i Input fasta file (reference).
-f Query file name.
-q K-mer size.
-s Minimum mutation score.
-t Number of building threads.
-Z Use minimizer in one level (bool value).
-w Window size is for minimizer computation and is valid only with -Z true.
-a Number of hash functions.
-M Maximum number of user bins in each filter.
-F Results file name.
-D Mapping threshold is used to assign a query as part of the sequence.
-x Length of peptide before dropping (L = K + x). Default: 1.

```
Example use case: 
```
./megapx multi_indexing -m blosum62 -b blacklist.txt -i input_references.fasta -f input_peptides.fasta -q 5 -s 100 -t 15 -Z 0 -a 2 -M 1000 -F output.log -D 0.85
```
### <a name="stat"></a>MegaPX stat
Print reference statistics by generating mapping and length files. The last evaluation step uses the mapping file to map each assignment score to the target reference name. 
```
./megapx stat -q kmer_size -i path/to/reference.fasta -m path/to/blosum62 -o path/to/output_dir
```

### <a name="mutate_stat"></a>MegaPX mutate_stat
Print mutation statistics with a number of mutated sequences. 
```
./megapx mutate_stat -q kmer_size -s min_mutation_score -i path/to/reference.fasta -m path/to/blosum62
```

### <a name="ibf_stat"></a>MegaPX ibf_stat
Print approximate IBF size with the number of technical bins (HIBF). 

```
./megapx ibf_stat -i path/to/file.fasta -m path/to/blosum62 -t threads -s min_mutation_score -q kmer_size -o path/to/output_dir -v path/to/output_dir/hashes.vect -a number_of_hash_functions

```
### <a name="write_db"></a>MegaPX write_db
Write the mutated database to an output FASTA file. However, utilizing this option for mutation is not advisable, as MegaPX will generate a FASTA file containing mutated sequences (one mutation per sequence). This approach consumes significant hard disk storage and results in slow runtime.
```
./megapx write_db -q kmer_size -s min_mutation_score -i test1.fasta -m path/to/blosum62 -o path/to/output_dir -t threads
```

### <a name="simulate_peptides"></a>MegaPX simulate_peptides
Simulate peptides from the given reference database and error rate. The number of errors sequence-wise should be assigned. 

```
./megapx simulate_peptides -i hendra_protein.fasta -n number_of_peptides -o path/to/output_dir -d number_of_errors
```

### <a name="build_vect"></a>MegaPX build_vect
For counting peptides into index, we recommend directly using multi-indexing command. 
Build first serialized index containing k_mer hashes `uint64_t`. The output is written to `path/to/output_dir/hashes_vect_q_s.vect` where `q`is the _k_-mer size and `s`is the user-definied minimum mutation score. 

```
./megapx build_vect -i path/to/file.fasta -m path/to/blosum62  -t threads -s min_mutation_score -q kmer_size -o path/to/output_dir -w window_size -Z bool_minimiser
```

### <a name="build_hibf"></a>MegaPX build_hibf
Build HIBF from the input index with several included bins. The output file is the log file, which contains the path to all generated HIBF(s): `path/to/output_dir/filter_HIBF_files.log`.
```
./megapx build_hibf -v path/to/output_dir/hashes_vect_q_s.vect -F path/to/output_dir/filter.hibf -E disable_estimation_ratio -R disable_rearrangement -g max_rearrangement_ratio -r maximum_false_positive_rate -p alpha -S sketch_bits -t threads -a number_of_hash_functions -M split_size

```

### <a name="count_hibf"></a>MegaPX count_hibf
Counting query peptide sequences into pre-built HIBF. _we highly recommend using the matching_threshold of 0_for the evaluation step. 
```
./megapx count_hibf -F filter_HIBF_files.log -f path/to/query.fasta -q kmer_size -d matching_threshold

```

### <a name="build_ibf"></a>MegaPX build_ibf
Build IBF from the input index with a number of included bins. 
```
./megapx build_ibf -v path/to/output_dir/hashes_vect_q_s.vect -o path/to/output_dir -a number_of_hash_functions

```

### <a name="count_ibf"></a>MegaPX count_ibf
Counting query peptide sequences into pre-built IBF.
```
./megapx count_ibf -F path/to/output_dir/hashes_vect_q_s.ibf -f path/to/query.fasta -q kmer_size

```

### <a name="hibf_ref"></a>MegaPX hibf_ref
Build HIBF from the input database without mutation assignment.
```
./megapx hibf_ref -q kmer_size -i path/to/file.fasta -F path/to/output_dir/filter.hibf -E disable_estimation_ratio -R disable_rearrangement -g max_rearrangement_ratio -r maximum_false_positive_rate -p alpha -S sketch_bits -t threads -a number_of_hash_functions -M split_size

```

### <a name="binary_search"></a>MegaPX binary_search
Run binary search-based classification on the serialized index.
```
./megapx binary_search -v path/to/output_dir/hashes_vect_q_s.vect -f path/to/query.fasta -q kmer_size -t threads

```

### <a name="build_fm"></a>MegaPX build_fm
Build FM index from input reference sequences, running the `stat`sub-commands is important to know the total number of sequences. 
```
./megapx build_fm -q kmer_size -t threads -s min_mutation_score -n total_number_of_sequences -i path/to/file.fasta -X path/to/output_dir/FM.idx -m path/to/blosum62 

```

### <a name="count_fm"></a>MegaPX count_fm
Count query sequences in the pre-built FM index.
```
./megapx count_fm -q kmer_size -n total_number_of_sequences -t threads -X path/to/output_dir/FM.idx  -f path/to/query.fasta

```

### <a name="evaluate"></a>MegaPX evaluate
Evaluate counting results, this step will generate a file containing two columns (assignment_score, sequence_name). The counting results file is named based on the query file name and the building/counting algorithm type. E.g. if the user build the index using IBF, the counting file will be: `path/to/output_dir/query_results_IBF.log`. The name combination makes the follow-rule simple for Snakemake pipeline (`query_file_name` + `_results_` + `algorithm_type` + `.log`). 
```
./megapx evaluate -C path/to/output_dir/query_results_algorithm.log -I path/to/output_dir/reference_id_map.log -V path/to/output_dir/evaluation_results.log -f path/to/query.fasta -D peptide_assignment_score -q kmer_size

```

### <a name="classification"></a>MegaPX classification
Species or Strain level normalization of counting results, this sub-command will be used for `refSeqViral` or `RVDB` viruses databases, where the sequences headers are named based on protein and species name (e.g. `> id protein_name [species/strain name]`).
```
./megapx classification -V path/to/output_dir/evaluation_results.log -T path/to/output_dir/classification_report.txt 

```
## <a name="use-cases "></a>Use-cases
### IBF Classification
IBF classification of a set of peptides against reference database using different commands. 
```
./megapx build_vect -i path/to/file.fasta -m path/to/blosum62  -t threads -s min_mutation_score -q kmer_size -o path/to/output_dir -w window_size -Z bool_minimiser
./megapx build_ibf -v path/to/output_dir/hashes_vect_q_s.vect -o path/to/output_dir -a number_of_hash_functions
./megapx count_ibf -F path/to/output_dir/hashes_vect_q_s.ibf -f path/to/query.fasta -q kmer_size
./megapx evaluate -C path/to/output_dir/query_results_algorithm.log -I path/to/output_dir/reference_id_map.log -V path/to/output_dir/evaluation_results.log -f path/to/query.fasta -D peptide_assignment_score -q kmer_size
./megapx classification -V path/to/output_dir/evaluation_results.log -T path/to/output_dir/classification_report.txt 
```

## <a name="parameters "></a>Parameters  
```
MegaPX - MegaPX builds and counts mutations from and in datasets with the classification of unknown samples.
============================================================================================================

SYNOPSIS
    [build_vect, build_ibf, count_ibf, build_hibf, count_hibf, binary_search, build_fm, count_fm, counting, stat,
    mutate_stat, evaluate, hibf_ref, classification, ibf_stat, write_db, mutate_seq_len, simulate_peptides, profile,
    test, multi_indexing] [OPTIONS]

DESCRIPTION
    MegaPX builds and counts mutations from and in datasets with the classification of unknown samples.

POSITIONAL ARGUMENTS
    ARGUMENT-1 (std::string)
          Modus to run MegaPX : Value must be one of
          [build_vect,build_ibf,count_ibf,build_hibf,count_hibf,binary_search,build_fm,count_fm,counting,stat,mutate_stat,evaluate,hibf_ref,classification,ibf_stat,write_db,mutate_seq_len,simulate_peptides,profile,test,multi_indexing].

OPTIONS

  Basic options:
    -h, --help
          Prints the help page.
    -hh, --advanced-help
          Prints the help page including advanced options.
    --version
          Prints the version information.
    --copyright
          Prints the copyright/license information.
    --export-help (std::string)
          Export the help page information. Value must be one of [html, man].
    --version-check (bool)
          Whether to check for the newest app version. Default: true.
    -t, --Threads (unsigned 8 bit integer)
          Number of input threads. Default: 1.
    -s, --Score (signed 32 bit integer)
          MinimumScore for each mutation. Default: 0.
    -q, --QMerSize (unsigned 8 bit integer)
          QMerSize for each k-mer. Default: 2.
    -i, --InputFasta (std::string)
          Input fasta file Default: .
    -m, --Matrix (std::string)
          Path to input matrix Default: .
    -o, --outputDir (std::string)
          Output directory for writing logs Default: .
    -Z, --minimiser (bool)
          Use minimiser in one level. Default: 0.
    -w, --windowSize (unsigned 8 bit integer)
          Window size for minimiser computation, valid only with -Z true. Default: 0.
    -N, --numberOfTopSequences (unsigned 8 bit integer)
          Number of top scored written sequences of each protein. Default: 1.
    -p, --alpha (double)
          Alpha value. Default: 1.2.
    -S, --sketchBits (unsigned 8 bit integer)
          HyperLogLog sketch bits. Default: 12.
    -r, --maximumFalsePositiveRate (double)
          Maximum false positive rate. Default: 0.01.
    -g, --maxRearrangementRatio (double)
          Maximum rearrangement ratio. Default: 0.5.
    -E, --disableEstimationRatio (bool)
          Disable estimation ratio. Default: 0.
    -R, --disableRearrangement (bool)
          Disable rearrangement ratio. Default: 0.
    -a, --numberOfHashFunctions (unsigned 8 bit integer)
          Number of used hashfunctions. Default: 2.
    -M, --maxUserBins (unsigned 64 bit integer)
          Maximum number of user bins in each HIBF. Default: 1.
    -d, --hibfThreshold (unsigned 8 bit integer)
          Threshold to assign query to user bin. Default: 1.
    -F, --filterFileName (std::string)
          HIBF output file name. Default: .
    -v, --vectFileName (std::string)
          Mutated DB input file name. Default: .
    -X, --indexFileName (std::string)
          Output index file name. Default: .
    -n, --numberOfReferenceSequences (unsigned 64 bit integer)
          Number of reference sequences. Default: 1.
    -f, --queryFileName (std::string)
          Query file name. Default: .
    -C, --countingResults (std::string)
          Counting results file name. Default: .
    -I, --referenceMapping (std::string)
          Reference ID mapping file. Default: .
    -V, --evaluationResults (std::string)
          Output file name for evaluation results. Default: .
    -D, --threshold (double)
          Mapping threshold to assign query as part of sequence. Default: 0.75.
    -T, --taxClassificationFile (std::string)
          Output file name for strain level classification. Default: .
    -b, --blackListFile (std::string)
          Path to black list file. Default: .

VERSION
    Last update: 2024
    MegaPX version: 0.5.4
    SeqAn version: 3.4.0-rc.4

LEGAL
    Author: Ahmad Lutfi
    Contact: mutht@rki.de & ahmad.lutfi.op@gmail.com
    SeqAn Copyright: 2006-2025 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.
```
















