#!/usr/bin/env python3
"""Compare native preparation/decoding with installed Cascadia on real spectra.

Requires the Python environment used to train/export the checkpoint. --index
is the Python-generated augmented HDF5 index for the same mzML/options.
Native input tensors can be generated with cascadia_native_tests mzML dump.pt.
"""
import argparse
import csv
import os
from pathlib import Path
import torch

os.environ.setdefault("NUMBA_CACHE_DIR", "/private/tmp/cascadia-validation-numba")
from cascadia.depthcharge.primitives import Peptide
from cascadia.depthcharge.tokenizers import PeptideTokenizer
from cascadia.depthcharge.data.spectrum_datasets import AnnotatedSpectrumDataset
from cascadia.depthcharge.data.preprocessing import scale_intensity, scale_to_unit_norm
from cascadia.model import AugmentedSpec2Pep

_original_split = Peptide.split

def compatible_split(self):
    if self.modifications is not None:
        self.modifications = [mods or None for mods in self.modifications]
    return _original_split(self)

Peptide.split = compatible_split


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--index", required=True)
    p.add_argument("--inputs", required=True)
    p.add_argument("--checkpoint", required=True)
    p.add_argument("--candidates")
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--batches", type=int, default=1)
    args = p.parse_args()
    torch.set_num_threads(4)
    t = PeptideTokenizer.from_massivekb(reverse=False, replace_isoleucine_with_leucine=True)
    dataset = AnnotatedSpectrumDataset(t, index_path=args.index, preprocessing_fn=[scale_intensity(scaling="root"), scale_to_unit_norm])
    native = torch.jit.load(args.inputs)
    assert len(dataset) == native.spectra.shape[0], "Candidate count differs"
    max_error = 0
    for begin, batch in enumerate(dataset.loader(batch_size=args.batch_size, num_workers=0)):
        start = begin * args.batch_size
        actual = native.spectra[start:start+len(batch[0]), :batch[0].shape[1]]
        try:
            torch.testing.assert_close(actual, batch[0], rtol=1e-6, atol=2e-7)
        except AssertionError:
            print(f"Input mismatch beginning at candidate {start}", flush=True)
            raise
        torch.testing.assert_close(native.precursors[start:start+len(batch[0])], batch[1], rtol=1e-6, atol=1e-4)
        max_error = max(max_error, (actual-batch[0]).abs().max().item())
    print(f"PASS: all {len(dataset)} native spectrum/precursor inputs match Python (max spectrum error {max_error:g})", flush=True)
    if not args.candidates:
        return
    model = AugmentedSpec2Pep.load_from_checkpoint(args.checkpoint, d_model=512, n_layers=9, n_head=8, dim_feedforward=1024, dropout=0, rt_width=2, tokenizer=t, max_charge=10, map_location="cpu", strict=True).eval()
    with open(args.candidates) as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    checked = 0
    for k, batch in enumerate(dataset.loader(batch_size=args.batch_size, num_workers=0)):
        if k >= args.batches:
            break
        prediction = model.predict_step(batch)
        actual = [row["sequence"] for row in candidates[checked:checked+len(batch[0])]]
        assert prediction[0] == actual, f"Prediction mismatch at batch {k}: Python={prediction[0]}, C++={actual}"
        checked += len(batch[0])
    print(f"PASS: {checked} native decoded peptides match original Cascadia predict_step", flush=True)

if __name__ == "__main__":
    main()
