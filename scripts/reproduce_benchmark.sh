#!/usr/bin/env bash
# Usage: CASCADIA_PYTHON=/path/to/cascadia_env/bin/python bash scripts/reproduce_benchmark.sh [data_dir] [output_dir]
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "$script_dir/.." && pwd)"
data_dir="${1:-$(cd "$repo_dir/.." && pwd)}"
benchmark_dir="${2:-$repo_dir/results/benchmark_reproduced}"
cascadia_python="${CASCADIA_PYTHON:-python}"

torch_prefix="$("$cascadia_python" -c 'import torch; print(torch.utils.cmake_prefix_path)')"
cmake -S "$repo_dir/src_libtorch" -B "$repo_dir/build_libtorch" \
  -DCMAKE_PREFIX_PATH="$torch_prefix" -DCMAKE_BUILD_TYPE=Release
cmake --build "$repo_dir/build_libtorch" -j 4
ctest --test-dir "$repo_dir/build_libtorch" --output-on-failure

for model in cascadia cascadia_pretrained cascadia_astral_tuned; do
  "$cascadia_python" "$script_dir/../src_libtorch/export_cascadia_torchscript.py" \
    --checkpoint "$data_dir/$model.ckpt" --output "$repo_dir/models/$model.pt"
done

exec "$cascadia_python" "$script_dir/benchmark_cascadia.py" run \
  --mzml "$data_dir/demo.mzML" \
  --model "$data_dir/cascadia.ckpt" "$repo_dir/models/cascadia.pt" \
  --model "$data_dir/cascadia_pretrained.ckpt" "$repo_dir/models/cascadia_pretrained.pt" \
  --model "$data_dir/cascadia_astral_tuned.ckpt" "$repo_dir/models/cascadia_astral_tuned.pt" \
  --cpp "$repo_dir/build_libtorch/megapx_sequence" --python "$cascadia_python" \
  --threads 4 --batch-size 32 --out "$benchmark_dir"
