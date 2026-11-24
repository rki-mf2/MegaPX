# Cascadia sequencing in C++

`megapx_sequence` reads DIA mzML, runs Cascadia with LibTorch, and writes SSL
and peptide FASTA files. An optional second step runs MegaPX `multi_indexing`.
Python is needed once per checkpoint for export; sequencing itself runs in C++.
The main MegaPX executable also accepts `--sequence` when built with
`MEGAPX_WITH_CASCADIA=ON`.

## Build the standalone sequencer

**Linux x86_64 and ARM64 automatic setup.** Run this block from the MegaPX repository root. The standalone build
requires CMake 3.18+, a C++17 compiler and zlib development files.

```bash
cmake -S src_libtorch -B build_libtorch -DCMAKE_BUILD_TYPE=Release \
  -U Torch_DIR -DCASCADIA_DOWNLOAD_LIBTORCH=ON
cmake --build build_libtorch -j 4
ctest --test-dir build_libtorch --output-on-failure
```

CMake first looks for an installed Torch package. On Linux x86_64 or ARM64, if it
cannot find one, it downloads Torch 2.0.1 for CPU, checks its SHA-256 and
reuses it in the build directory. On x86_64 it uses the C++11 ABI SDK.
On ARM64 (`uname -m` reports `aarch64` or `arm64`), CMake extracts the
C++ SDK from the official ARM64 wheel without running Python or pip. The first
download needs internet access. Building and sequencing an exported model
require no Python Torch installation. Python is used for checkpoint export.

The configure command also recovers an existing build directory after an
offline attempt by clearing the cached Torch path and enabling downloads.

**Alternative for an installed SDK.** Use this configuration only when you
already have LibTorch installed, for example for macOS, CUDA or offline use.
Replace the example path with the directory containing your `TorchConfig.cmake`:

```bash
cmake -S src_libtorch -B build_libtorch -DCMAKE_BUILD_TYPE=Release \
  -DTorch_DIR=/absolute/path/to/libtorch/share/cmake/Torch \
  -DCASCADIA_DOWNLOAD_LIBTORCH=OFF
```

You can also use an existing Python environment's PyTorch installation.
Activate the environment used for export, then run:

```bash
torch_prefix="$(python -c 'import torch; print(torch.utils.cmake_prefix_path)')" &&
cmake -S src_libtorch -B build_libtorch -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$torch_prefix"
```

If this prints `ModuleNotFoundError: No module named 'torch'`, use the
automatic Linux build above, supply a standalone SDK, or activate the
environment containing PyTorch.

To include sequencing in the main MegaPX executable, enable Cascadia.
The same package discovery and download options apply:

```bash
cmake -S src -B build -DMEGAPX_WITH_CASCADIA=ON
cmake --build build -j 4
```

The main build still needs MegaPX's existing SeqAn3, HIBF, OpenMP and other
search dependencies. Use a LibTorch build with the C++ standard-library ABI
matching the search dependencies (normally the cxx11-ABI build on Linux).
The ARM64 wheel uses the old libstdc++ ABI. For an integrated ARM64 build
with C++11 ABI search dependencies, supply a matching SDK via `Torch_DIR`.
The standalone sequencer can also call an already-built MegaPX executable.

## Export each checkpoint once

```bash
python src_libtorch/export_cascadia_torchscript.py \
  --checkpoint ../cascadia.ckpt --output models/cascadia.pt
```

Repeat for `cascadia_pretrained.ckpt` and `cascadia_astral_tuned.ckpt` to get
separate `.pt` bundles. A Lightning `.ckpt` contains weights rather than an
executable C++ model; renaming it to `.pt` does not convert it.

To keep a frozen export alongside the default compact bundle, give it a
separate filename:

```bash
python src_libtorch/export_cascadia_torchscript.py \
  --checkpoint ../cascadia.ckpt --output models/cascadia_frozen.pt --freeze
```

The same naming convention applies to `cascadia_pretrained_frozen.pt` and
`cascadia_astral_tuned_frozen.pt`. Select either export using the native
`--model` option; no separate freezing flag is needed during sequencing.

The exporter infers embedding width, layer count, feedforward width, charge
embedding size and vocabulary size from the checkpoint and restores weights
strictly. It compares eager and exported forwards at several batch/peak/token
lengths, including empty token prefixes and padding, before saving. The bundle
contains the exact ordered token vocabulary, residue/modification masses,
architecture metadata and the source checkpoint SHA-256. The C++ runtime reads
that vocabulary instead of assuming a fixed list.

The public Cascadia models use the canonical MassiveKB vocabulary and eight
attention heads. Checkpoints without tokenizer metadata cannot reveal token
meanings from their tensor shapes. For a custom model, supply its training
metadata rather than choosing a vocabulary just because its size matches:

```bash
python src_libtorch/export_cascadia_torchscript.py \
  --checkpoint custom.ckpt --output models/custom.pt \
  --n-head 8 --modifications training_ptms.json
```

`training_ptms.json` is a list of ProForma peptide examples containing the PTMs
used during training. Alternatively, `--vocabulary training_vocabulary.json`
accepts `{"tokens": [null, "$", "A", ...], "masses": {"A": 71.037113805, ...}}`
with the **exact training token order**, including padding at index zero.
Reverse tokenizers and unrelated architectures require a matching exporter;
this runtime supports Cascadia's forward tokenizer and AugmentedSpec2Pep model.
It cannot guarantee compatibility with arbitrary future architectures.

The exporter normalizes the empty modification lists returned by Pyteomics 5,
which otherwise create two false zero-mass tokens in Cascadia 0.0.7. It also
avoids PyTorch 2.0 attention fast-path differences and nested-tensor trimming.
Export errors stop the process; weights are never truncated or padded to fit.

## Sequence an mzML file

From the MegaPX repository:

```bash
./build_libtorch/megapx_sequence \
  --sequence ../demo.mzML --model models/cascadia.pt --out results/demo
```

The output directory contains:

- `predictions.ssl`: accepted, scored predictions; modification masses are
  written as signed shifts. The default mode preserves Python's time-unit
  heuristic: divide seconds by 60 when any candidate retention time exceeds
  500 seconds, otherwise retain seconds. `--output-mode token` always writes minutes.
- `peptides.fasta`: the same accepted predictions with modifications removed
  from the sequence. Every other SSL field is copied into its header as
  `field=value`, separated by underscores. Named modifications and their
  1-based residue positions follow the identifier. Terminal positions are
  `N-term` or `C-term`; unmodified peptides have `modification: none`.
- `candidates.tsv`: every charge candidate, including rejected predictions,
  unmodified sequence, named modification positions, confidence, calculated
  precursor m/z and rejection reason.
- `model.json`: the model bundle's embedded provenance.
- `memory.txt`: the same usage report printed to the terminal, with total
  elapsed time, sequencing time, CPU time (user + system), and peak resident
  memory in MiB. With `--search`, it includes separate search time, CPU and
  child-process peak memory; total CPU includes both processes. Peak memory
  is reported per process, rather than adding peaks that may occur at different
  times. These measurements describe host RAM. Failed runs also report usage.

Example FASTA record:

```text
>file=demo.mzML_scan=12_charge=2_score-type=UNKNOWN_score=0.95_retention-time=24_start-time=23.9_end-time=24.1 [modification: Oxidation@3, Carbamidomethyl@5]
AMMDCGR
```

Whitespace in header field values is percent-encoded. Each candidate keeps its
own record even if another charge candidate predicts the same peptide.
Existing result files are protected from overwriting; choose a fresh `--out`.

Useful options:

```text
--device auto|cpu|cuda    auto selects CUDA when available, otherwise CPU
--threads 4              CPU intra-op threads
--batch-size 32
--score-threshold 0.8
--max-length 32           original Python decoding iterations
--output-mode python      default: reproduce Cascadia 0.0.7 scoring/filtering
--output-mode token       optional token-based score/per-window mass filter
--max-charge 4           candidates have charges 2 through this value
--scan-width 1
--top-peaks 150
--limit 12               first N candidates, for smoke tests
```

C++ caches the spectrum encoder once per batch, then reuses its output during
peptide decoding. This removes the repeated encoder work in Python's original
`predict_step`; speed still depends on hardware, batch size and the model.

Input preparation follows Cascadia 0.0.7, including its legacy precursor
transformation, m/z sorting, scan selection and float32 normalization. The
**default** output mode also reproduces Python's behavior: 32 autoregressive
steps, continued predictions after the first stop, confidence summed using the
length of the peptide string (including modification labels), and the final
MS2 isolation-window width used for every mass filter. This preserves the
original tool's conventions even where they are surprising. The default does
not reject a candidate simply because no stop appeared within 32 steps.

`--output-mode token` selects the earlier alternative: score through the first
stop, reject unfinished sequences, and use each candidate's isolation width.
It intentionally gives different SSL rows/scores and is not used for Python
agreement benchmarks. The SSL file path/number formatting can differ between
languages; agreement requires identical reported candidates, order and peptide
strings, with explicitly stated tolerances for floating-point columns.

The mzML reader supports ordinary 32/64-bit little-endian float arrays with
zlib or no compression, and scan times in seconds or minutes. MS-Numpress and
other encodings need conversion first. DIA MS2 isolation windows and nearby
MS1 spectra are required. This is not a general replacement for OpenMS.

## Sequence, then run multi_indexing

Only `multi_indexing` is supported by this handoff. Pass its ordinary options
after `--`; the pipeline supplies the peptide query and results paths:

```bash
./build_libtorch/megapx_sequence \
  --sequence ../demo.mzML --model models/cascadia.pt --out results/demo_search \
  --search --megapx /absolute/path/to/megapx -- \
  -i reference.fasta -m data/matrix/blosum62 \
  -q 5 -s 100 -t 4 -Z 0 -a 2 -M 1000 -D 0.85
```

An integrated build uses the same syntax with `./build/main/megapx --sequence`.
The optional subprocess runs without a shell. Paths passed with `-i`, `-m`
and `-b` are resolved before changing to the result directory. Omit `-f`, `-F`
and `-o`; those paths are managed by the pipeline. Sequencing files are closed
before search starts. Search writes `search.tsv`, MegaPX's additional reports,
and `search.log`. A search failure returns a nonzero status and preserves all
sequencing outputs. A result with zero accepted peptides stops before search.

For a plumbing check, `peptides.fasta` may itself be used as the reference.
This only tests handoff and self-matching; use an independent protein database
for meaningful classification. Model inference uses I-to-L replacement as in
Cascadia, which should be considered when preparing a search reference.

## Validate against Python

The native unit test can also dump all prepared demo tensors:

```bash
./build_libtorch/cascadia_native_tests ../demo.mzML results/native_inputs.pt
python tests/validate_cascadia.py \
  --index ../cascadia_11-21-14:35:52/index.hdf5 \
  --inputs results/native_inputs.pt --checkpoint ../cascadia.ckpt \
  --candidates results/demo/candidates.tsv --batch-size 32
```

The HDF5 index must come from Python augmentation of the same input and options.
The validator checks every prepared spectrum/precursor input, then compares a
batch of decoded peptides with original Cascadia `predict_step`. Use `--batches`
to check more batches. `--dump-tensors PATH` on the sequencer saves only its
first batch; the native test utility saves all batches for this comparison.

`torch_loader` remains available for checkpoint inspection and older TOML-based
experiments. Use `megapx_sequence` for the SSL/FASTA/search workflow.

The subprocess/output integration tests use real Cascadia inference with a
**stub search executable**, so they can run without SeqAn/HIBF dependencies:

```bash
python tests/test_cascadia_pipeline.py \
  --sequencer build_libtorch/megapx_sequence \
  --bundle models/cascadia.pt --mzml ../demo.mzML
```


## Reproduce all three Python/C++ benchmarks

From the parent `cascadia` folder, with the existing environment active:

```bash
python MegaPX/scripts/benchmark_cascadia.py run \
  --mzml demo.mzML \
  --model cascadia.ckpt MegaPX/models/cascadia.pt \
  --model cascadia_pretrained.ckpt MegaPX/models/cascadia_pretrained.pt \
  --model cascadia_astral_tuned.ckpt MegaPX/models/cascadia_astral_tuned.pt \
  --cpp MegaPX/build_libtorch/megapx_sequence \
  --python "$(command -v python)" \
  --threads 4 --batch-size 32 --out MegaPX/results/benchmark_reproduced
```

To rebuild the native tool, export all checkpoints, and run the same benchmark
in one command:

```bash
CASCADIA_PYTHON="$(command -v python)" \
  bash MegaPX/scripts/reproduce_benchmark.sh "$PWD" "$PWD/MegaPX/results/benchmark_fresh"
```

The benchmark runs **six full sequencing processes sequentially**, C++ then
Python for each checkpoint, on the CPU with four threads, batch size 32, score
threshold 0.8 and maximum charge 4. Export time is excluded. Each measured run
includes process startup/imports, model loading, mzML preprocessing, sequencing
and writing results. No predictions are substituted or read from earlier runs.
The reference wrapper calls installed Cascadia's original `sequence()` and
`predict_step` with unchanged scoring/filtering. It normalizes the Pyteomics 5
empty-modification lists, uses a writable Numba cache, and sets zero data-loader
workers to avoid the macOS pickle problem. It also supplies the Unimod XML
snapshot bundled with psims to the original Pyteomics relational database
loader, avoiding its network download. That same database is reused by the
unmodified output writer in the process. The database is initialized before
inference (included in each measured process's runtime), and a separate writer
preflight checks ordinary, named-modified and numeric-terminal-modified peptides
before any long benchmark starts. The snapshot's SHA-256 is logged. CPU/thread settings are
controlled equally for both implementations. The benchmark records those
settings rather than comparing a CPU run to a GPU run.

Results are saved under the chosen output directory:

- `benchmark.tsv`: one row per backend/model, with elapsed seconds, CPU seconds
  (user + system), peak RSS in MiB, row count and agreement result.
- `results.json`: complete measurements and executed commands.
- `manifest.json`: input/checkpoint/bundle/executable hashes and run settings;
  the bundle's embedded checkpoint hash is checked before running.
- `<model>/python.log` and `cpp.log`: complete stdout/stderr logs plus the end-to-end
  elapsed/CPU/peak-RSS summary.
- `<model>/python.metrics.json` and `cpp.metrics.json`: independently measured
  process runtime and memory, captured by the same `wait4` mechanism.
- `<model>/agreement.json`: exact candidate/sequence/order checks, missing and
  additional rows, maximum numeric errors, and individual mismatches.
- `<model>/python/predictions.ssl` and `cpp/predictions.ssl`: actual outputs.

The parent benchmark process does not import Torch or use cumulative
`RUSAGE_CHILDREN` peak values. Each child's peak is collected separately;
macOS bytes and Linux KiB are converted to MiB. Peak RSS describes host memory,
not GPU VRAM. These are individual end-to-end runs, not repeated-trial averages;
OS caching and other machine activity can affect timings. Do not run the
backends concurrently for a speed comparison. Automated `run` commands hold a
temporary macOS idle-sleep assertion that is released at exit. Elapsed time uses
Python's monotonic performance counter; metrics also record UTC start/end and
wall-clock duration, so sleep/clock gaps can be distinguished. On this macOS host, the performance
counter excludes system sleep while the progress bar uses wall time. Avoid closing
the laptop or suspending the machine during a performance run.

Default agreement tolerances are absolute score error <= `1e-6` and time error
<= `1e-5` minutes. Candidate keys, row order, modification strings and score
labels must match exactly. File names are compared by basename, while exact
path equality and byte equality are reported separately. Score and retention-time
values are also compared after round-tripping to Python's original float32
precision, which distinguishes harmless decimal formatting from numeric changes. A failed agreement
makes the benchmark exit nonzero; it is never silently called equal because
peptide sets merely overlap. Tighten tolerances with `--score-atol` and
`--time-atol`. `--resume` reuses completed measurements only when the stored
input/model/executable/settings manifest and commands are unchanged.

### Measure either tool independently

From the parent `cascadia` folder:

```bash
CUDA_VISIBLE_DEVICES="" python MegaPX/scripts/benchmark_cascadia.py measure \
  --metrics python.metrics.json --stdout python.log -- \
  python MegaPX/scripts/run_cascadia_python.py sequence \
  demo.mzML cascadia.ckpt -o python_results -b 32 -t 0.8 -c 4

python MegaPX/scripts/benchmark_cascadia.py measure \
  --metrics cpp.metrics.json --stdout cpp.log -- \
  ./MegaPX/build_libtorch/megapx_sequence \
  --sequence demo.mzML --model MegaPX/models/cascadia.pt --out cpp_results \
  --device cpu --threads 4 --batch-size 32 --score-threshold 0.8 --max-charge 4

python MegaPX/scripts/benchmark_cascadia.py compare \
  --python-ssl python_results.ssl --cpp-ssl cpp_results/predictions.ssl \
  --report agreement.json
```

Use fresh output names. `measure` prints the runtime/memory summary, writes its
JSON metrics file and captures tool output in the requested log. The C++ tool
also prints its own stage report and writes `<out>/memory.txt`. `compare` can
be run on any existing pair of SSL files, including after independent manual
runs, and returns exit status 1 when agreement fails.


## Complete command-line reference

### Measured memory and runtime

Full `demo.mzML` runs on this macOS ARM64 host, CPU with four threads, batch
size 32, Cascadia 0.0.7 and Torch 2.0.1. Each model processed all 2,880
augmented candidates with the original Python scoring/filtering defaults.
Times include process startup, model loading, preprocessing, inference and
output; checkpoint export is excluded. Memory is peak process RSS in MiB.

| Checkpoint | SSL rows, both tools | Python seconds | Python MiB | Compact C++ seconds | Compact C++ MiB | Frozen C++ MiB |
|---|---:|---:|---:|---:|---:|---:|
| `cascadia` | 220 | 1189.7 | 2293.7 | 362.6 | 847.9 | 3258.6 |
| `cascadia_pretrained` | 40 | 1690.7 | 2204.7 | 362.0 | 870.3 | 3248.0 |
| `cascadia_astral_tuned` | 35 | 1665.9 | 2304.4 | 363.4 | 963.4 | 3172.2 |

All three comparisons passed: identical candidate order and peptide/modification
strings, with scores and retention times identical after float32 round-trip.
SSL files are not byte-identical because numeric formatting differs.
Compact exports reduced peak native RAM by **69.6–74.0%** versus frozen exports.
Frozen C++ runs took 254.8, 257.3 and 366.5 seconds respectively: the lower-RAM
export has a speed tradeoff on these individual runs, though it remained faster
than Python here. These are measurements, not repeated-trial averages or
performance guarantees.

Python figures come from the completed six-process baseline; compact C++ was
measured again and checked against those same saved Python predictions. The
pretrained Python run experienced host sleep: its 1690.7 seconds is the measured
performance-counter duration excluding sleep, not wall-clock turnaround. New
benchmark runs prevent macOS idle sleep and record both clocks. Local logs,
commands, CPU times, hashes and per-model agreement reports are published in the
[anonymized benchmark record](../docs/benchmarks/2026-11-21/README.md).
Redundant generated results and archived binaries were removed after retaining
that evidence. The record contains reports rather than raw predictions; use the
commands above to regenerate SSL/FASTA output and fresh reference runs.

### Memory tuning

Export uses `--no-freeze` by default. Keeping shared scripted module methods
avoids the large inlined TorchScript graphs produced by `--freeze`. Re-export
existing frozen bundles to get this change; rebuilding the C++ executable
alone does not change an old bundle. Both forms retain the original weights
and use the same inference/scoring settings.

`--jit-optimize off` disables runtime JIT profiling/optimization for comparison.
Reducing `--batch-size` can also reduce tensor memory, but may change floating-point
results; keep batch size 32 for the tested Python agreement configuration.
Host peak RSS includes model loading and inference, and is not GPU VRAM.
The reader currently prepares candidates for the full mzML file in memory;
larger inputs can therefore increase peak RAM even at a fixed batch size.

To measure another native build or export against completed Python reference
outputs without repeating Python inference, from the repository root:

```bash
python scripts/benchmark_cpp_variant.py \
  --reference-run results/benchmark_reproduced \
  --cpp build_libtorch/megapx_sequence --models models \
  --out results/benchmark_variant --jit-optimize on
```

This checks input/checkpoint hashes, runs C++ again for all three bundles,
and compares every output against Python. Its manifest and Python metrics
explicitly identify the reused Python measurements. Use the full reproduction
script above when fresh measurements of both implementations are needed.

### Native sequencing CLI

All options below apply to `megapx_sequence`, and to an integrated `megapx
--sequence` build. Positive integer options reject zero and negative values.
Architecture and vocabulary settings belong to export, not runtime overrides.

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

Charge candidates must fit the model's trained charge embedding. Increasing
`--max-length`, `--batch-size`, `--top-peaks` or `--scan-width` can increase
memory/runtime and can change results. Matching the original Python defaults
requires batch size 32, length 32, score threshold 0.8, charges 2–4, scan width
1, top 150 peaks and `--output-mode python`. CPU and CUDA runs can have small
floating-point differences; benchmark both tools on the same device/settings.

### Checkpoint exporter

| Option | Default | Meaning |
|---|---|---|
| `--checkpoint PATH` | Required | Trusted original Lightning checkpoint. |
| `--output PATH` | Required | TorchScript `.pt` bundle; also writes a provenance `.pt.json` sidecar. |
| `--device DEVICE` | `cpu` | Torch device used during export and parity validation. |
| `--threads N` | `4` | CPU thread count during export. |
| `--freeze`, `--no-freeze` | `--no-freeze` | Keep compact scripted module methods by default; freezing inlines the graph and increased host RAM in local measurements. |
| `--d-model N` | Inferred | Embedding width; an explicit mismatch is rejected. |
| `--n-layers N` | Inferred | Decoder/encoder depth supported by this exporter; mismatch rejected. |
| `--dim-feedforward N` | Inferred | Feedforward width; mismatch rejected. |
| `--max-charge N` | Inferred | Charge-embedding configuration; mismatch rejected. |
| `--n-tokens N` | Inferred | Vocabulary size excluding padding; mismatch rejected. |
| `--n-head N` | Checkpoint metadata or `8` | Attention heads; cannot be recovered from matrix shapes alone. |
| `--dropout X` | `0.0` | Constructor dropout setting; inference runs in eval mode. |
| `--rt-width X` | `2.0` | Constructor retention-time encoding width. |
| `--modifications PATH` | Canonical MassiveKB | Training PTMs as JSON ProForma examples; retains sorted tokenizer semantics. |
| `--vocabulary PATH` | Checkpoint/canonical tokenizer | Explicit ordered training tokens and masses JSON; takes precedence over PTM examples. |
| `--example-batch N` | `2` | Validation batch size for the configurable example. |
| `--example-peaks N` | `16` | Validation peak count for the configurable example. |
| `--example-sequence-length N` | `4` | Validation prefix length; additional empty-prefix/padded cases are always checked. |

An unrelated architecture or unknown/custom vocabulary must be exported with
its matching model implementation and training metadata. Tensor size alone
does not establish token meaning.

### Inspection and TOML experiments

`torch_loader MODEL` prints checkpoint/TorchScript metadata.
`torch_loader CONFIG.toml` uses the settings in
[`cascadia_sequence.example.toml`](cascadia_sequence.example.toml) for the older
SSL-only experiment. Its architecture labels do not rewrite an exported model.
Use the native sequencing CLI above for the full SSL/FASTA/search workflow.

### Benchmark CLI

Run `python scripts/benchmark_cascadia.py --help` and each subcommand's `--help`
for the complete benchmark arguments. `run` accepts repeated `--model CKPT PT`
pairs, `--mzml`, `--cpp`, `--python`, `--out`, `--threads`, `--batch-size`,
`--max-charge`, `--score-threshold`, `--score-atol`, `--time-atol` and `--resume`.
`measure` accepts `--metrics`, `--stdout`, optional `--cwd`, then `-- COMMAND`.
`compare` accepts `--python-ssl`, `--cpp-ssl`, `--report`, `--score-atol`, and
`--time-atol`. Full examples and artifact descriptions appear in the benchmark
section above.
