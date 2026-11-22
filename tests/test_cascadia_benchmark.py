import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("benchmark", Path(__file__).parents[1] / "scripts/benchmark_cascadia.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def ssl(self, name, score="0.9", sequence="AM[+15.99]K", scan="0"):
        path = self.root / name
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(benchmark.FIELDS)
            writer.writerow(["demo.mzML", scan, "2", sequence, "UNKNOWN", score, "24.0", "23.9", "24.1"])
        return path

    def test_agreement_checks_numeric_and_peptide_results(self):
        a = self.ssl("python.ssl")
        self.assertTrue(benchmark.compare(a, self.ssl("rounding.ssl", score="0.90000005"))["agreement"])
        self.assertFalse(benchmark.compare(a, self.ssl("score.ssl", score="0.95"))["agreement"])
        self.assertFalse(benchmark.compare(a, self.ssl("peptide.ssl", sequence="AMK"))["agreement"])
        missing = benchmark.compare(a, self.ssl("candidate.ssl", scan="1"))
        self.assertFalse(missing["agreement"])
        self.assertEqual(len(missing["python_only"]), 1)
        self.assertEqual(len(missing["cpp_only"]), 1)

    def test_duplicate_keys_are_rejected(self):
        a = self.ssl("duplicate.ssl")
        with a.open("a") as handle:
            handle.write(a.read_text().splitlines()[1] + "\n")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            benchmark.read_ssl(a)

    def test_peak_rss_belongs_to_each_child(self):
        first = benchmark.measure([sys.executable, "-c", "x=bytearray(60*1024*1024); print(len(x))"], self.root/"large.json", self.root/"large.log")
        second = benchmark.measure([sys.executable, "-c", "print('small')"], self.root/"small.json", self.root/"small.log")
        self.assertEqual(first["exit_code"], 0)
        self.assertEqual(second["exit_code"], 0)
        self.assertGreater(first["peak_rss_mib"] - second["peak_rss_mib"], 20)
        self.assertGreater(first["elapsed_seconds"], 0)
        self.assertEqual((self.root/"small.log").read_text().splitlines()[0], "small")
        self.assertIn("Peak process RSS:", (self.root/"small.log").read_text())


if __name__ == "__main__":
    unittest.main(verbosity=2)
