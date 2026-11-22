#!/usr/bin/env python3
"""Run the original Cascadia CLI with dependency/macOS compatibility fixes.

Scoring, filtering and predict_step are unchanged. CPU thread and data-loader
worker counts are controlled by CASCADIA_TORCH_THREADS (4) and
CASCADIA_NUM_WORKERS (0) for reproducible comparisons with the native runner.
"""
import os
from pathlib import Path
import sys
import tempfile
import gzip
import importlib.resources
import json
import hashlib

cache = Path(tempfile.gettempdir()) / "cascadia-benchmark-cache"
cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(cache / "numba"))
os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))

import torch
from cascadia.depthcharge.primitives import Peptide
from cascadia.depthcharge.data.spectrum_datasets import AnnotatedSpectrumDataset

_original_split = Peptide.split
_original_loader = AnnotatedSpectrumDataset.loader


def compatible_split(self):
    if self.modifications is not None:
        self.modifications = [mods or None for mods in self.modifications]
    return _original_split(self)


def compatible_loader(self, *args, **kwargs):
    kwargs["num_workers"] = int(os.environ.get("CASCADIA_NUM_WORKERS", "0"))
    return _original_loader(self, *args, **kwargs)


Peptide.split = compatible_split
AnnotatedSpectrumDataset.loader = compatible_loader
# Use the same Unimod relational schema from the psims packaged XML snapshot.
# Build it before inference so a missing output dependency fails immediately.
from pyteomics import mass
xml_resource = importlib.resources.files("psims.controlled_vocabulary.vendor").joinpath("unimod_tables.xml.gz")
xml_path = cache / "unimod_tables.xml"
with xml_resource.open("rb") as handle:
    xml = gzip.decompress(handle.read())
if not xml_path.exists() or xml_path.read_bytes() != xml:
    xml_path.write_bytes(xml)
mass.unimod._unimod_xml_download_url = xml_path.as_uri()
_original_unimod = mass.unimod.Unimod
_unimod_database = _original_unimod()
# Resolve all canonical named PTMs now; no network call is needed at output.
for name in ["Acetyl", "Carbamyl", "Ammonia-loss", "Oxidation", "Deamidated", "Carbamidomethyl"]:
    _unimod_database.by_name(name).monoisotopic_mass


def offline_unimod(path=None):
    return _unimod_database if path is None else _original_unimod(path)


mass.unimod.Unimod = offline_unimod

# Keep the original Python computation, including its output conventions.
torch.set_num_threads(int(os.environ.get("CASCADIA_TORCH_THREADS", "4")))

if __name__ == "__main__":
    if sys.argv[1:] == ["--preflight"]:
        import csv
        import math
        from cascadia.utils import write_results
        with tempfile.TemporaryDirectory(prefix="cascadia-output-preflight-") as directory:
            prefix = str(Path(directory) / "check")
            precursors = [((mz - 1.00727646688) * 2 - 1.007276) * 2 for mz in [659.0, 720.0, 659.0, 659.0]]
            predictions = [(["HGESAWNLENR", "SVEEYANC[Carbamidomethyl]HLAR", "[Acetyl]-HGESAWNLENR", "[+25.980265]-HGESAWNLENR"], [],
                            torch.tensor([math.log(0.95), math.log(0.9), math.log(0.9), math.log(0.9)]), torch.ones(4, 32),
                            torch.tensor([1434.0, 1435.0, 1436.0, 1437.0]), torch.tensor([[value, 2.0] for value in precursors]))]
            write_results(predictions, prefix, "preflight.mzML", 31.0, 0.8, 8.5)
            with open(prefix + ".ssl") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            assert len(rows) == 4 and rows[1]["sequence"] == "SVEEYANC[+57.02]HLAR"
            assert rows[2]["sequence"] == "[+42.01]-HGESAWNLENR"
            assert rows[3]["sequence"] == "[+25.980265]-HGESAWNLENR"
        print(json.dumps({"preflight": "passed", "unimod_source": str(xml_resource), "unimod_xml_sha256": hashlib.sha256(xml).hexdigest(), "writer_rows": len(rows)}))
    else:
        from cascadia.cascadia import main
        main()
