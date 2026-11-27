# Cascadia verification record

This is the retained test evidence for the native Cascadia loader. Logs and
metrics below are checked into the repository so readers can inspect what was
actually run. Local paths have been replaced with `${WORKSPACE}`, `${CASCADIA_ENV}`,
`${HOME}`, `${TMP}` and `${SYSTEM_PYTHON}`. Measurements, predictions printed in logs, and model/input
SHA-256 values are unchanged. Any hashes of original manifest text refer to the
original files before path redaction; they are not hashes of this edited archive.
Historical command/output paths describe the original runs, not files still
present in this cleaned repository.

## What “compact” and “frozen” mean

These are two TorchScript exports of the **same model weights**, using the
**same C++ loader**. They are not different model sizes or quantized weights.

- **Frozen export (earlier baseline):** `torch.jit.freeze` replaces parameter
  lookups with constants and inlines model methods into a large graph, enabling
  further optimization. It used more loader RAM in these tests.
- **Compact export (current default):** `--no-freeze` preserves shared scripted
  methods. The original weights and inference settings are retained.

Layer count, learned weight values and numerical precision are unchanged.
“Compact” describes the execution representation, not weight compression or a
smaller neural network.

The compact export reduced host peak RAM by 69.6–74.0% on this machine.
Frozen export was faster in two of these three individual measurements; the
current default prioritizes the lower measured RAM use.

## Full-file agreement and performance

Input: `demo.mzML`, SHA-256
`fc0a6391a7f0e947a106864e1a8bb6dc7a538f8ecc4f74289c6cfc3cd0cdb225`.
All 2,880 augmented candidates were processed per run, using CPU, four threads,
batch size 32, 32 decoding steps, charges 2–4 and score threshold 0.8.
See the [environment](baseline/environment.json), [baseline manifest](baseline/manifest.json)
and [compact manifest](compact/manifest.json) for versions, commands, settings
and checkpoint hashes. Python used the compatibility wrapper described in the
[loader documentation](../../../src_libtorch/README.md#reproduce-all-three-pythonc-benchmarks).

| Checkpoint | Accepted rows, both tools | Python seconds | Python MiB | Current C++ seconds | Current C++ MiB | Earlier frozen C++ MiB |
|---|---:|---:|---:|---:|---:|---:|
| `cascadia` | 220 | 1189.7 | 2293.7 | 362.6 | 847.9 | 3258.6 |
| `cascadia_pretrained` | 40 | 1690.7 | 2204.7 | 362.0 | 870.3 | 3248.0 |
| `cascadia_astral_tuned` | 35 | 1665.9 | 2304.4 | 363.4 | 963.4 | 3172.2 |

All three reports passed: candidate keys/order, peptide and modification strings,
and score types match. Scores and retention times are also identical after
float32 round-trip. Decimal formatting means SSL files are not byte-identical.
Agreement reports include explicit tolerances and maximum numeric differences.

| Checkpoint | Current C++ log | Python log | Agreement report |
|---|---|---|---|
| Cascadia | [C++](compact/cascadia/cpp.log) | [Python](baseline/cascadia/python.log) | [JSON](compact/cascadia/agreement.json) |
| Pretrained | [C++](compact/cascadia_pretrained/cpp.log) | [Python](baseline/cascadia_pretrained/python.log) | [JSON](compact/cascadia_pretrained/agreement.json) |
| Astral tuned | [C++](compact/cascadia_astral_tuned/cpp.log) | [Python](baseline/cascadia_astral_tuned/python.log) | [JSON](compact/cascadia_astral_tuned/agreement.json) |

Each log has an adjacent `.metrics.json` file with its command, exit status,
elapsed time, CPU time and per-process peak RSS. Browse the
[current summary TSV](compact/benchmark.tsv), [baseline summary TSV](baseline/benchmark.tsv)
and [memory comparison JSON](compact/memory_comparison.json). Earlier frozen C++
logs are retained in each model's `baseline/` folder.

Python was measured during the original full six-process benchmark. The compact
variant reran all three C++ processes and compared against the saved Python
predictions; **Python timings were reused, not measured again**. Export time is
excluded. These are individual end-to-end runs, not repeated-trial averages.
Peak RSS is host RAM, not GPU VRAM. The pretrained Python run included host sleep;
its reported elapsed value excludes sleep on this Mac. See the
[timing observation](baseline/timing_observations.json). New benchmark runs
prevent idle sleep and record UTC/wall-clock duration as well.

## Additional checks

These short checks were rerun while preparing this record; they are separate
from the full-file timings above. See [exact commands](tests/commands.json).

| Check | Scope | Retained log |
|---|---|---|
| Native CTest | Vocabulary, masses, modifications, FASTA header, tiny mzML fixture; benchmark tests | [Passed, 2/2](tests/native-ctest.log) |
| macOS UBSan CTest | Same checks with UndefinedBehaviorSanitizer | [Passed, 2/2](tests/macos-ubsan-ctest.log) |
| Pipeline tests | Real model inference, SSL/FASTA output, overwrite protection, search success/failure handoff | [Passed, 3/3](tests/pipeline.log) |
| Benchmark tests | Agreement failures, duplicate keys, independent per-child RSS | [Passed, 3/3](tests/benchmark-tests.log) |
| Doxygen | Generate API documentation without warnings | [Passed](tests/doxygen.log) |

The pipeline test uses a stub search executable: it checks the handoff contract,
not classification accuracy. Linux ASan/UBSan is configured in CI; this local
record demonstrates macOS UBSan, not an executed Linux CI run. The normal local
executable and all three working bundles were hash-verified as byte-identical
to the full-file tested artifacts; see [production hashes](compact/production_paths.json).

## Reproduce

From the directory containing `MegaPX/` and the three original checkpoints:

```bash
conda activate cascadia_env
bash MegaPX/scripts/reproduce_benchmark.sh
```

See [all benchmark commands and comparison rules](../../../src_libtorch/README.md#reproduce-all-three-pythonc-benchmarks).
Generated predictions, temporary smoke runs, failed trials and archived model
binaries were removed during cleanup. This archive keeps logs and reports,
not raw SSL/FASTA predictions; rerun the benchmark to regenerate those files or
to create a reference directory for `benchmark_cpp_variant.py`.
