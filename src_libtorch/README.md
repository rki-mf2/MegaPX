# src_libtorch

Minimal C++/CMake smoke test for the local LibTorch install. The program can
print a random tensor or summarize a pretrained checkpoint/model file.

## Prerequisites

- CMake 3.18 or newer
- A C++20 compiler
- LibTorch extracted at the repository root as `libtorch/`

## Build

From the repository root:

```bash
cmake -S src_libtorch -B build_libtorch \
  -DCMAKE_PREFIX_PATH="$PWD/libtorch" \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build_libtorch --config Release
```

Or, from an existing build directory:

```bash
cmake -DCMAKE_PREFIX_PATH=../libtorch ../src_libtorch
cmake --build . --config Release
```

If you previously configured the build while a Conda toolchain was being picked
up, remove the old build directory or reconfigure into a fresh one. Stale
`CMakeCache.txt` entries can keep pointing at Conda's `ld`, `ar`, or `nm`.

## Run

```bash
./build_libtorch/torch_loader
```

Expected output is a random tensor similar to:

```text
 0.3030  0.1593  0.9083
 0.9048  0.6349  0.7786
[ CPUFloatType{2,3} ]
```

To inspect a checkpoint:

```bash
./build_libtorch/torch_loader model/cascadia.ckpt
```

You can also pass a Cascadia sequence TOML config:

```bash
./build_libtorch/torch_loader src_libtorch/cascadia_sequence.example.toml
```

## Running Inference

LibTorch can only run a TorchScript export, not the original PyTorch-Lightning
`.ckpt` file. First export the checkpoint from a Python environment with
Cascadia and PyTorch installed:

```bash
python src_libtorch/export_cascadia_torchscript.py \
  --checkpoint model/cascadia.ckpt \
  --output model/cascadia_sequence.pt
```

Then update the TOML:

```toml
[model]
path = "../model/cascadia_sequence.pt"
```

Now the same command will read `demo.mzML`, build the Cascadia tensors, run
greedy sequence decoding, print the first decoded candidate peptides, and write
filtered SSL results to the configured `output_path`:

```bash
./build_libtorch/torch_loader src_libtorch/cascadia_sequence.example.toml
```

`CheckpointModelLoader` supports two cases:

- TorchScript files saved with `torch.jit.save(...)`: loads the module with
  LibTorch and reports methods, parameters, buffers, tensor shapes, dtypes, and
  devices.
- Python `torch.save(...)` or PyTorch-Lightning `.ckpt` zip checkpoints:
  reports archive metadata, checkpoint keys, PyTorch-Lightning version when
  present, sample `state_dict` names, tensor-storage payload counts, and the
  largest storage entries. These files need the original Python model class to
  restore weights, so they cannot be executed directly by LibTorch until
  exported as TorchScript.

`InferenceParameterInspector` then reports the parameters needed for inference:

- For TorchScript models, it reads the `forward` schema and prints required
  inputs, optional inputs, return types, and runtime tensor requirements.
- For Python `.ckpt` files, it prints the required conversion/runtime items
  that are still missing from C++ inference, plus best-effort hints from
  `state_dict` names.

`CascadiaModelConfig` reads the sequence-inference TOML file. The config keeps
the model path and the Python sequence defaults from
`cascadia/cascadia/cascadia.py`:

```toml
[model]
path = "../model/cascadia.ckpt"
d_model = 512
n_layers = 9
n_head = 8
dim_feedforward = 1024
dropout = 0
rt_width = 2
max_charge = 10
tokenizer = "massivekb"

[runtime]
device = "cpu"

[sequence]
batch_size = 32
augmentation_width = 2
max_sequence_length = 64
score_threshold = 0.8
output_path = "demo_results"
```

The C++ runner follows Cascadia's `--out` convention: `output_path =
"demo_results"` writes `demo_results.ssl`; setting a path that already ends in
`.ssl` uses that exact file name.

`CascadiaSequenceForward` is the tensor-level C++ inference wrapper. It expects
an exported TorchScript model whose `forward` accepts:

- `spectra`: float tensor `[batch, peaks, 4]` containing `(m/z, intensity,
  retention_time, ms_level)`
- `precursors`: float tensor `[batch, 2]` containing `(neutral_precursor_mass,
  charge)`
- `partial_sequence_tokens`: integer tensor `[batch, sequence_length]`

It returns token logits, and optionally the precursor and fragment predictions
if the TorchScript export returns the same tuple as Python `_forward_step`.

## mzML Input

`CascadiaMzmlReader` reads centroided mzML files like `cascadia/demo.mzML` and
mirrors the tensor preparation used by `cascadia/cascadia/augment.py`, without
writing an intermediate ASF file.

For each spectrum it extracts:

- MS level
- scan start time, converted from minutes to seconds
- m/z and intensity arrays from mzML binary data
- for MS2 scans: isolation window target m/z, lower/upper offsets, selected ion
  m/z, and precursor charge

The reader then builds augmented spectra by combining nearby MS2 peaks with
nearby MS1 peaks. It applies the same top-N peak selection and intensity
normalization pattern used in the Python augmentation code:

- MS2 intensity: square root, then normalize by max intensity
- MS1 intensity: square root twice, then normalize by max intensity

The resulting tensors are:

- `spectra`: `[candidate_count, padded_peak_count, 4]` containing `(m/z,
  intensity, retention_time_delta, ms_level)`
- `precursors`: `[candidate_count, 2]` containing `(neutral_precursor_mass,
  charge)`

The local reader uses standard C++ plus zlib and supports the ProteoWizard-style
zlib-compressed 32/64-bit float arrays in the included demo file. For a broader
production mzML implementation, prefer a dedicated C++ MS library such as
OpenMS `MzMLFile` or ProteoWizard `msdata`.

## Notes

- `library kineto not found` is a warning emitted by this LibTorch package. It
  does not stop this sample from building or running.
- `undefined reference to ... GLIBC_PRIVATE` during linking usually means the
  system compiler is being mixed with Conda's linker/sysroot. The CMake target
  prefers `/usr/bin/ld` for system GNU compiler builds to avoid that mismatch.
