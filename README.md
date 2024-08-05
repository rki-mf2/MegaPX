# MegaPX: Fast Peptide Assignment Method Using Index Databases with Various Counting Algorithms

<!-- ## Citation  -->

## Table of contents

* [Description](#description)
* [Installation](#installation)
* [Commands](#commands)
* [Use-cases](#use-cases)
* [Parameters](#parameters )

## <a name="description"></a>Description
MegaPX is designed to mutate protein databases, enabling users to search for de novo peptides in both the original and mutated versions. The tool employs a k-mer-based search, providing rapid and highly accurate peptide assignment capability. Beyond this primary function, the tool has numerous other use cases, including metaproteomics classification, virus assignment, and peptide classification.   
MegaPX utilizes four distinct counting algorithms: [IBF](https://github.com/seqan/hibf/), [HBF](https://github.com/seqan/hibf/), [FM-Index](https://github.com/seqan/seqan3), and binary search. As we assessed these four algorithms in terms of runtime and performance, it became imperative to incorporate them into a single tool, allowing users the freedom to choose based on their specific use case. Mutation generation relies on a substitution matrix containing values proportional to the likelihood of amino acid `i` mutating into amino acid `j` across all possible pairs of amino acids. Each reference amino acid sequence undergoes k-merization, and mutated k-mers are generated based on the user-defined score, which is the sum of k-mer alignment scores. The tool is designed in step commands, which will be included as rules in the Snakemake pipeline. 
## <a name="installation"></a>Installation   
MegaPX runs only under Linux x86_64; we recommend the Conda installations for all dependencies by:   
```
conda create --name MegaPX_env OR  conda create --prefix path/to/MegaPX_env  
conda activate MegaPX_env  OR conda activate path/to/MegaPX_env     
conda install -c conda-forge cxx-compiler   
conda install conda-forge::gcc  
conda install conda-forge::clang
```
Building MegaPX from source:   
```
git clone https://github.com/lutfia95/MegaPX.git   
cd MegaPX  
mkdir build && cd build   
cmake ../src   
make
```
The executable will be built in `path/to/MegaPX/build/main/`   
Users can also download the binary pre-built version: [Linux x86_64](https://github.com/lutfia95/MegaPX/releases/download/v.0.0.0/MegaPX-0.0.0-Linux_x64.tar.gz)
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
Build first serialized index containing k_mer hashes `uint64_t`. The output is written to `path/to/output_dir/hashes_vect_q_s.vect` where `q`is the k-mer size and `s`is the user-definied minimum mutation score. 

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

### Mulit-Indexing 
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

```
Example use case: 
```
./megapx multi_indexing -m blosum62 -b black_list.txt -i refSeqViral_monkeypox.fasta -f monkeypox_sample/E02292_MonkeyPox_SP3_DDA_1.fasta -q 5 -s 100 -t 40 -Z 0 -a 2 -M 1000 -F monckey.log -D 0.8
```


## <a name="parameters "></a>Parameters  
```
DESCRIPTION
    MegaPX builds and counts mutations from and in datasets with the classification of unknown samples.

POSITIONAL ARGUMENTS
    ARGUMENT-1 (std::string)
          Modus to run MegaPX: Value must be one of
          [build_vect, build_ibf, count_ibf, build_hibf, count_hibf, binary_search, build_fm, count_fm, counting, stat,
            mutate_stat, evaluate, hibf_ref, classification, ibf_stat, write_db, mutate_seq_len, simulate_peptides, profile,
            test, multi_indexing].

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
    MegaPX version: 0.0.0
    SeqAn version: 3.4.0-rc.1

LEGAL
    Author: Ahmad Lutfi
    Contact: ahmad.lutfi@fu-berlin.de
    SeqAn Copyright: 2006-2023 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.
```

















