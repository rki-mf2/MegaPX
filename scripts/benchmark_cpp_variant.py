#!/usr/bin/env python3
"""Measure a C++/bundle variant against completed Python reference runs.

No Python inference is rerun. Original Python timing/command provenance is kept
and marked as reused; every C++ run and every agreement check is performed anew.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import zipfile
from benchmark_cascadia import compare, keep_machine_awake, measure, save_json, sha256, write_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-run", type=Path, required=True)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--models", type=Path, required=True, help="Directory of <checkpoint-stem>.pt bundles")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jit-optimize", choices=["on", "off"], default="on")
    args = parser.parse_args()
    reference, out, cpp = args.reference_run.resolve(), args.out.resolve(), args.cpp.resolve()
    manifest = json.loads((reference/"manifest.json").read_text())
    assert sha256(manifest["mzml"]) == manifest["mzml_sha256"], "Input mzML changed"
    if out.exists():
        raise ValueError("Choose a fresh variant --out directory")
    bundles = []
    for model in manifest["models"]:
        name = Path(model["checkpoint"]).stem
        bundle = (args.models/name).with_suffix(".pt").resolve()
        with zipfile.ZipFile(bundle) as archive:
            key = next(n for n in archive.namelist() if n.endswith("/extra/cascadia.metadata.json"))
            meta = json.loads(archive.read(key))
        assert meta["checkpoint_sha256"] == model["checkpoint_sha256"], f"Wrong weights in {bundle}"
        bundles.append({"name": name, "path": str(bundle), "sha256": sha256(bundle), "metadata": meta})
    out.mkdir(parents=True)
    save_json(out/"manifest.json", {"reference_run": str(reference), "reference_manifest_sha256": sha256(reference/"manifest.json"),
                                 "cpp": str(cpp), "cpp_sha256": sha256(cpp), "jit_optimize": args.jit_optimize,
                                 "models": bundles, "settings": manifest,
                                 "python_measurements": "Reused unchanged from the completed reference run, not measured again"})
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(manifest["threads"]), MKL_NUM_THREADS=str(manifest["threads"]))
    results = []
    all_agree = True
    for bundle in bundles:
        name = bundle["name"]
        case = out/name
        case.mkdir()
        command = [str(cpp), "--sequence", manifest["mzml"], "--model", bundle["path"], "--out", str(case/"cpp"),
                   "--device", "cpu", "--threads", str(manifest["threads"]), "--batch-size", str(manifest["batch_size"]),
                   "--max-charge", str(manifest["max_charge"]), "--score-threshold", str(manifest["score_threshold"]),
                   "--jit-optimize", args.jit_optimize]
        measured = measure(command, case/"cpp.metrics.json", case/"cpp.log", env=env)
        if measured["exit_code"]:
            raise RuntimeError(f"C++ failed for {name}")
        ref_ssl = reference/name/"python/predictions.ssl"
        agreement = compare(ref_ssl, case/"cpp/predictions.ssl", manifest["score_atol"], manifest["time_atol"])
        save_json(case/"agreement.json", agreement)
        original = json.loads((reference/name/"python.metrics.json").read_text())
        if original["exit_code"] != 0:
            raise RuntimeError("Reference Python run was unsuccessful")
        copied = case/"python"
        copied.mkdir()
        shutil.copy2(ref_ssl, copied/"predictions.ssl")
        original.update(model=name, backend="python", agreement=agreement["agreement"], rows=agreement["python_rows"],
                        measurement_reused_from=str(reference/name/"python.metrics.json"))
        save_json(case/"python.reference.metrics.json", original)
        measured.update(model=name, backend="cpp", agreement=agreement["agreement"], rows=agreement["cpp_rows"])
        results.extend([measured, original])
        write_summary(out, results)
        all_agree &= agreement["agreement"]
        print(f"{name}: {'PASS' if agreement['agreement'] else 'FAIL'}, rows={agreement['cpp_rows']}, float32 scores identical={agreement['score_float32_roundtrip_equal']}", flush=True)
    save_json(out/"results.json", {"all_agree": all_agree, "runs": results})
    return 0 if all_agree else 1


if __name__ == "__main__":
    try:
        with keep_machine_awake():
            raise SystemExit(main())
    except Exception as error:
        print(f"Variant benchmark error: {error}", file=sys.stderr)
        raise SystemExit(1)
