#!/usr/bin/env python3
"""Real-model sequencing/output tests with a stubbed search subprocess.

This validates the handoff contract; it does not test MegaPX classification.
"""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="cascadia pipeline ")
        cls.root = Path(cls.temp.name)
        cls.reference = cls.root / "reference.fasta"
        cls.reference.write_text(">self_reference\nHGESAWNLENR\n")
        cls.matrix = cls.root / "test_matrix"
        cls.matrix.write_text("test fixture only\n")
        cls.fake = cls.root / "fake_megapx.py"
        cls.fake.write_text("#!" + sys.executable + "\n" + '''
import csv, json, pathlib, re, sys
args = sys.argv[1:]
assert args[0] == "multi_indexing"
query = pathlib.Path(args[args.index("-f") + 1])
output = pathlib.Path(args[args.index("-F") + 1])
assert query.is_absolute() and output.is_absolute()
assert query.read_text().startswith(">file=")
assert " [modification: " in query.read_text()
assert all(re.fullmatch("[A-Z]+", line) for line in query.read_text().splitlines() if not line.startswith(">"))
with query.with_name("predictions.ssl").open() as handle:
    assert len(list(csv.DictReader(handle, delimiter="\\t"))) == 1
query.with_name("received_arguments.json").write_text(json.dumps(args))
if "--fail-test" in args:
    sys.exit(7)
output.write_text("stub_handoff_only\\n")
''')
        cls.fake.chmod(0o755)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def command(self, name):
        return [OPTIONS.sequencer, "--sequence", OPTIONS.mzml, "--model", OPTIONS.bundle,
                "--out", str(self.root / name), "--limit", "1", "--batch-size", "1"]

    def search(self):
        return ["--search", "--megapx", str(self.fake), "--", "-i", str(self.reference), "-m", str(self.matrix)]

    def test_sequence_outputs_and_search_handoff(self):
        result = subprocess.run(self.command("success") + self.search(), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = self.root / "success"
        with (out / "predictions.ssl").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        header, sequence = (out / "peptides.fasta").read_text().splitlines()
        self.assertEqual(len(rows), 1)
        self.assertEqual(sequence, re.sub(r"\[.*?\]|-", "", rows[0]["sequence"]))
        for key, value in rows[0].items():
            if key != "sequence":
                self.assertIn(key + "=" + value, header)
        self.assertTrue((out / "search.tsv").exists())
        report = (out / "memory.txt").read_text()
        self.assertIn(report, result.stdout)
        self.assertIn("Status                 : complete", report)
        for label in ["Real time (total)", "Real time (sequencing)", "CPU time (total)", "Peak RSS (sequencer)", "Real time (search)", "Peak RSS (search)"]:
            self.assertIn(label, report)
        self.assertRegex(report, r"Peak RSS \(sequencer\)\s+: [1-9][0-9]*\.[0-9]+ MiB")
        retry = subprocess.run(self.command("success"), text=True, capture_output=True)
        self.assertNotEqual(retry.returncode, 0)
        self.assertIn("Output already exists", retry.stderr)

    def test_search_failure_keeps_sequencing(self):
        result = subprocess.run(self.command("failure") + self.search() + ["--fail-test"], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sequencing outputs retained", result.stderr)
        for name in ["peptides.fasta", "predictions.ssl", "candidates.tsv", "search.log", "memory.txt"]:
            self.assertTrue((self.root / "failure" / name).exists())

        report = (self.root / "failure" / "memory.txt").read_text()
        self.assertIn("Status                 : failed", report)
        self.assertIn("Peak RSS (search)", report)

    def test_ckpt_requires_export(self):
        cmd = self.command("unexported")
        cmd[cmd.index("--model") + 1] = "unexported.ckpt"
        result = subprocess.run(cmd, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Export it once", result.stderr)
        self.assertFalse((self.root / "unexported").exists())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequencer", required=True)
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--mzml", required=True)
    OPTIONS = parser.parse_args()
    for key in vars(OPTIONS):
        setattr(OPTIONS, key, str(Path(getattr(OPTIONS, key)).resolve()))
    unittest.main(argv=[sys.argv[0]], verbosity=2)
