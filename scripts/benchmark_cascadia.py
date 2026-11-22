#!/usr/bin/env python3
"""Measure Python/C++ Cascadia in separate processes and check SSL agreement.

Only the Python standard library is needed by the benchmark driver. The chosen
--python interpreter must have Cascadia installed. Measurements include process
startup, imports, model loading, mzML preprocessing, inference and output I/O.
Checkpoint export is excluded. Runs are sequential, with identical CPU thread,
batch, charge and score settings. Peak RSS uses wait4, with OS-specific units.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess
import struct
import sys
import time
import zipfile

REPO = Path(__file__).resolve().parents[1]
NUMERIC = {"score": 1e-6, "retention-time": 1e-5, "start-time": 1e-5, "end-time": 1e-5}
FIELDS = ["file", "scan", "charge", "sequence", "score-type", *NUMERIC]


def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


@contextmanager
def keep_machine_awake():
    keeper = None
    if sys.platform == "darwin" and Path("/usr/bin/caffeinate").exists():
        keeper = subprocess.Popen(["/usr/bin/caffeinate", "-i", "-w", str(os.getpid())],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("Preventing macOS idle sleep for this benchmark only", flush=True)
    try:
        yield
    finally:
        if keeper is not None:
            if keeper.poll() is None:
                keeper.terminate()
            keeper.wait()


def append_usage_summary(path, metrics):
    marker = "=== End-to-end benchmark measurement ==="
    path = Path(path)
    if marker in path.read_text(errors="replace"):
        return
    with path.open("a") as handle:
        handle.write(f"\n{marker}\n")
        handle.write(f"Elapsed time: {metrics['elapsed_seconds']:.6f} sec\n")
        handle.write(f"CPU time (user + system): {metrics['cpu_seconds']:.6f} sec\n")
        handle.write(f"Peak process RSS: {metrics['peak_rss_mib']:.6f} MiB\n")
        handle.write(f"Exit code: {metrics['exit_code']}\n")


def measure(command, metrics_path, stdout_path, cwd=None, env=None):
    """Measure precisely this child, without cumulative RUSAGE_CHILDREN peaks."""
    metrics_path, stdout_path = Path(metrics_path), Path(stdout_path)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    print("Running: " + shlex.join(command), flush=True)
    wall_started = time.time()
    started = time.perf_counter()
    with stdout_path.open("w") as log:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
        while True:
            try:
                _, status, usage = os.wait4(process.pid, 0)
                break
            except InterruptedError:
                continue
        code = os.waitstatus_to_exitcode(status)
        process.returncode = code  # wait4 already reaped this exact child.
    elapsed = time.perf_counter() - started
    wall_completed = time.time()
    data = {"command": command, "shell_command": shlex.join(command),
            "cwd": str(Path(cwd or os.getcwd()).resolve()),
            "elapsed_seconds": elapsed,
            "wall_clock_elapsed_seconds": wall_completed - wall_started,
            "started_at_utc": datetime.fromtimestamp(wall_started, timezone.utc).isoformat(),
            "completed_at_utc": datetime.fromtimestamp(wall_completed, timezone.utc).isoformat(),
            "elapsed_clock": time.get_clock_info("perf_counter").implementation,
            "cpu_seconds": usage.ru_utime + usage.ru_stime,
            "peak_rss_mib": usage.ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024),
            "exit_code": code, "stdout_log": str(stdout_path.resolve()),
            "memory_measurement": "peak process RSS from wait4; host RAM, not GPU VRAM"}
    save_json(metrics_path, data)
    append_usage_summary(stdout_path, data)
    print(f"Finished: exit={code}, elapsed={data['elapsed_seconds']:.3f}s, CPU={data['cpu_seconds']:.3f}s, peak RSS={data['peak_rss_mib']:.3f} MiB", flush=True)
    return data


def read_ssl(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Unexpected SSL columns in {path}: {reader.fieldnames}")
        rows = list(reader)
    indexed = {}
    order = []
    for row in rows:
        key = (Path(row["file"]).name, int(row["scan"]), int(row["charge"]))
        if key in indexed:
            raise ValueError(f"Duplicate spectrum/charge key {key} in {path}")
        for field in NUMERIC:
            if not math.isfinite(float(row[field])):
                raise ValueError(f"Nonfinite {field} for {key} in {path}")
        indexed[key] = row
        order.append(key)
    return indexed, order


def compare(python_ssl, cpp_ssl, score_atol=1e-6, time_atol=1e-5):
    a, ao = read_ssl(python_ssl)
    b, bo = read_ssl(cpp_ssl)
    common = a.keys() & b.keys()
    missing, extra = sorted(a.keys() - b.keys()), sorted(b.keys() - a.keys())
    categorical = []
    numeric = []
    maximum = {field: 0.0 for field in NUMERIC}
    for key in sorted(common):
        for field in ["sequence", "score-type"]:
            if a[key][field] != b[key][field]:
                categorical.append({"key": key, "field": field, "python": a[key][field], "cpp": b[key][field]})
        for field in NUMERIC:
            delta = abs(float(a[key][field]) - float(b[key][field]))
            maximum[field] = max(maximum[field], delta)
            tolerance = score_atol if field == "score" else time_atol
            if delta > tolerance:
                numeric.append({"key": key, "field": field, "python": a[key][field], "cpp": b[key][field], "absolute_error": delta})
    agreement = not missing and not extra and not categorical and not numeric and ao == bo
    return {"agreement": agreement, "python_rows": len(a), "cpp_rows": len(b),
            "shared_rows": len(common), "python_only": missing, "cpp_only": extra,
            "categorical_mismatches": categorical, "numeric_mismatches": numeric,
            "row_order_equal": ao == bo, "maximum_absolute_errors": maximum,
            "score_absolute_tolerance": score_atol, "time_absolute_tolerance_minutes": time_atol,
            "file_paths_text_equal": all(a[k]["file"] == b[k]["file"] for k in common),
            "score_float32_roundtrip_equal": not missing and not extra and all(struct.pack("f", float(a[k]["score"])) == struct.pack("f", float(b[k]["score"])) for k in common),
            "retention_time_float32_roundtrip_equal": not missing and not extra and all(struct.pack("f", float(a[k]["retention-time"])) == struct.pack("f", float(b[k]["retention-time"])) for k in common),
            "byte_identical": Path(python_ssl).read_bytes() == Path(cpp_ssl).read_bytes(),
            "comparison_rule": "Exact candidate keys, order, peptide/modification strings and score types; numeric absolute tolerances; file column matched by basename"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def write_summary(out, results):
    fields = ["model", "backend", "elapsed_seconds", "cpu_seconds", "peak_rss_mib", "exit_code", "agreement", "rows", "stdout_log"]
    with (out / "benchmark.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in results)


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / "manifest.json").exists() and not args.resume:
        raise ValueError("Benchmark already exists; use a fresh --out or --resume")
    mzml, cpp, interpreter = args.mzml.resolve(), args.cpp.resolve(), args.python.resolve()
    pairs = [(Path(c).resolve(), Path(b).resolve()) for c, b in args.model]
    names = [c.stem for c, _ in pairs]
    if len(names) != len(set(names)):
        raise ValueError("Checkpoint stems must be unique")
    env = dict(os.environ, CUDA_VISIBLE_DEVICES="", CASCADIA_TORCH_THREADS=str(args.threads),
               CASCADIA_NUM_WORKERS="0", OMP_NUM_THREADS=str(args.threads), MKL_NUM_THREADS=str(args.threads))
    worker = REPO / "scripts/run_cascadia_python.py"
    manifest = {"mzml": str(mzml), "mzml_sha256": sha256(mzml), "cpp": str(cpp),
                "cpp_sha256": sha256(cpp), "python": str(interpreter), "python_runner_sha256": sha256(worker),
                "platform": platform.platform(), "threads": args.threads, "batch_size": args.batch_size,
                "score_threshold": args.score_threshold, "max_charge": args.max_charge,
                "device": "cpu", "python_data_loader_workers": 0,
                "execution": "sequential; end-to-end process timings; export excluded; one measured run per backend/model",
                "score_atol": args.score_atol, "time_atol": args.time_atol, "models": []}
    for checkpoint, bundle in pairs:
        checkpoint_hash = sha256(checkpoint)
        with zipfile.ZipFile(bundle) as archive:
            key = next(n for n in archive.namelist() if n.endswith("/extra/cascadia.metadata.json"))
            metadata = json.loads(archive.read(key))
        if metadata["checkpoint_sha256"] != checkpoint_hash:
            raise ValueError(f"Bundle {bundle} does not come from {checkpoint}")
        manifest["models"].append({"checkpoint": str(checkpoint), "checkpoint_sha256": checkpoint_hash,
                                   "bundle": str(bundle), "bundle_sha256": sha256(bundle)})
    preflight = subprocess.run([str(interpreter), str(worker), "--preflight"], env=env, text=True, capture_output=True)
    (out / "python_preflight.log").write_text(preflight.stdout + preflight.stderr)
    if preflight.returncode != 0:
        raise RuntimeError(f"Python output preflight failed; inspect {out / 'python_preflight.log'}")
    manifest["offline_unimod"] = json.loads(preflight.stdout.splitlines()[-1])
    if args.resume and (out / "manifest.json").exists():
        previous = json.loads((out / "manifest.json").read_text())
        if previous != manifest:
            raise ValueError("Resume manifest differs: inputs, executable, model or benchmark settings changed")
    # Capture the reference environment outside the measured processes.
    version_code = "import importlib.metadata as m,json,platform; packages=['cascadia','torch','pytorch-lightning','lightning','numpy','numba','pyteomics','spectrum-utils','psims']; print(json.dumps({'python':platform.python_version(),'platform':platform.platform(),'packages':{p:m.version(p) for p in packages}},indent=2))"
    versions = subprocess.run([str(interpreter), "-c", version_code], check=True, text=True, capture_output=True)
    save_json(out / "environment.json", json.loads(versions.stdout))
    save_json(out / "manifest.json", manifest)
    results = []
    all_agree = True
    for checkpoint, bundle in pairs:
        case = out / checkpoint.stem
        case.mkdir(exist_ok=True)
        python_dir, cpp_dir = case / "python", case / "cpp"
        python_dir.mkdir(exist_ok=True)
        commands = {
            "python": [str(interpreter), str(worker), "sequence", str(mzml), str(checkpoint),
                       "-o", str(python_dir / "predictions"), "-b", str(args.batch_size),
                       "-t", str(args.score_threshold), "-c", str(args.max_charge)],
            "cpp": [str(cpp), "--sequence", str(mzml), "--model", str(bundle), "--out", str(cpp_dir),
                    "--threads", str(args.threads), "--batch-size", str(args.batch_size),
                    "--score-threshold", str(args.score_threshold), "--max-charge", str(args.max_charge),
                    "--device", "cpu"]}
        case_results = []
        for backend in ["cpp", "python"]:
            metrics = case / f"{backend}.metrics.json"
            ssl = (python_dir if backend == "python" else cpp_dir) / "predictions.ssl"
            if args.resume and metrics.exists():
                result = json.loads(metrics.read_text())
                if result["command"] != commands[backend] or result["exit_code"] != 0 or not ssl.exists():
                    raise ValueError(f"Cannot resume incomplete/changed run: {metrics}; use a fresh --out")
                print(f"Reusing completed {checkpoint.stem}/{backend}", flush=True)
            else:
                result = measure(commands[backend], metrics, case / f"{backend}.log", cwd=python_dir if backend == "python" else case, env=env)
            result = dict(result, model=checkpoint.stem, backend=backend)
            results.append(result)
            case_results.append(result)
            write_summary(out, results)
            if result["exit_code"] != 0:
                raise RuntimeError(f"{backend} failed for {checkpoint.name}; inspect {result['stdout_log']}")
        agreement = compare(python_dir / "predictions.ssl", cpp_dir / "predictions.ssl", args.score_atol, args.time_atol)
        save_json(case / "agreement.json", agreement)
        for result in case_results:
            result["agreement"] = agreement["agreement"]
            result["rows"] = agreement[f"{result['backend']}_rows"]
        write_summary(out, results)
        all_agree &= agreement["agreement"]
        print(f"Agreement {checkpoint.name}: {'PASS' if agreement['agreement'] else 'FAIL'}; Python={agreement['python_rows']}, C++={agreement['cpp_rows']}; max deltas={agreement['maximum_absolute_errors']}", flush=True)
    save_json(out / "results.json", {"all_agree": all_agree, "runs": results})
    return 0 if all_agree else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    measure_parser = sub.add_parser("measure", help="Measure any one Python or C++ command")
    measure_parser.add_argument("--metrics", type=Path, required=True)
    measure_parser.add_argument("--stdout", type=Path, required=True)
    measure_parser.add_argument("--cwd", type=Path)
    measure_parser.add_argument("command", nargs=argparse.REMAINDER)
    compare_parser = sub.add_parser("compare", help="Compare existing SSL outputs; exit nonzero on disagreement")
    compare_parser.add_argument("--python-ssl", type=Path, required=True)
    compare_parser.add_argument("--cpp-ssl", type=Path, required=True)
    compare_parser.add_argument("--report", type=Path, required=True)
    run_parser = sub.add_parser("run", help="Benchmark Python and C++ sequentially for every model pair")
    run_parser.add_argument("--mzml", type=Path, required=True)
    run_parser.add_argument("--model", nargs=2, action="append", required=True, metavar=("CHECKPOINT", "BUNDLE"))
    run_parser.add_argument("--cpp", type=Path, required=True)
    run_parser.add_argument("--python", type=Path, default=Path(sys.executable))
    run_parser.add_argument("--out", type=Path, required=True)
    run_parser.add_argument("--threads", type=int, default=4)
    run_parser.add_argument("--batch-size", type=int, default=32)
    run_parser.add_argument("--max-charge", type=int, default=4)
    run_parser.add_argument("--score-threshold", type=float, default=0.8)
    run_parser.add_argument("--resume", action="store_true")
    for target in [compare_parser, run_parser]:
        target.add_argument("--score-atol", type=float, default=1e-6)
        target.add_argument("--time-atol", type=float, default=1e-5)
    args = parser.parse_args()
    if args.action == "measure":
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        if not command:
            parser.error("measure needs a command after --")
        return 0 if measure(command, args.metrics, args.stdout, args.cwd)["exit_code"] == 0 else 1
    if args.action == "compare":
        result = compare(args.python_ssl, args.cpp_ssl, args.score_atol, args.time_atol)
        save_json(args.report, result)
        print(json.dumps(result, indent=2))
        return 0 if result["agreement"] else 1
    if args.threads < 1 or args.batch_size < 1 or args.max_charge < 2 or not 0 <= args.score_threshold <= 1:
        parser.error("Invalid thread, batch, charge or threshold setting")
    with keep_machine_awake():
        return run(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Benchmark error: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1)
