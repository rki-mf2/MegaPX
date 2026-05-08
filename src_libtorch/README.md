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
./build_libtorch/torch_loader src_libtorch/model/cascadia.ckpt
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

## Notes

- `library kineto not found` is a warning emitted by this LibTorch package. It
  does not stop this sample from building or running.
- `undefined reference to ... GLIBC_PRIVATE` during linking usually means the
  system compiler is being mixed with Conda's linker/sysroot. The CMake target
  prefers `/usr/bin/ld` for system GNU compiler builds to avoid that mismatch.
