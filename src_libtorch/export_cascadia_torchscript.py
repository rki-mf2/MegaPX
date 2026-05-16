#!/usr/bin/env python3
"""Export Cascadia's Lightning checkpoint to a LibTorch-loadable file.

Run from the repository root with a Python environment that has Cascadia,
PyTorch, and PyTorch-Lightning installed:

    python export_cascadia_torchscript.py --checkpoint cascadia.ckpt --output cascadia_sequence.pt
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import torch

_torch_load = torch.load


def _torch_load_trusted_checkpoint(*args, **kwargs):
    kwargs.setdefault("weights_only", False)
    return _torch_load(*args, **kwargs)


torch.load = _torch_load_trusted_checkpoint


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--d-model", default=512, type=int)
    parser.add_argument("--n-layers", default=9, type=int)
    parser.add_argument("--n-head", default=8, type=int)
    parser.add_argument("--dim-feedforward", default=1024, type=int)
    parser.add_argument("--dropout", default=0.0, type=float)
    parser.add_argument("--rt-width", default=2.0, type=float)
    parser.add_argument("--max-charge", default=10, type=int)
    parser.add_argument("--n-tokens", default=28, type=int)
    parser.add_argument("--example-batch", default=2, type=int)
    parser.add_argument("--example-peaks", default=16, type=int)
    parser.add_argument("--example-sequence-length", default=4, type=int)
    return parser.parse_args()


class CascadiaSequenceModule(torch.nn.Module):
    """Small TorchScript-facing wrapper around Cascadia's _forward_step."""

    def __init__(self, model: torch.nn.Module) -> None:
        super().__init__()
        self.peak_encoder = model.spectrum_encoder.peak_encoder
        self.latent_spectrum = model.spectrum_encoder.latent_spectrum
        self.spectrum_transformer_encoder = model.spectrum_encoder.transformer_encoder
        self.aa_encoder = model.decoder.aa_encoder
        self.mass_encoder = model.decoder.mass_encoder
        self.charge_encoder = model.decoder.charge_encoder
        self.peptide_positional_encoder = model.decoder.positional_encoder
        self.transformer_decoder = model.decoder.transformer_decoder
        self.final = model.decoder.final
        self.prec_layer = model.prec_layer
        self.frag_layer = model.frag_layer

    def forward(
        self,
        spectra: torch.Tensor,
        precursors: torch.Tensor,
        partial_sequence_tokens: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        peak_padding_mask = ~spectra.sum(dim=2).to(torch.bool)
        spectrum_token_mask = torch.zeros(
            (spectra.shape[0], 1),
            dtype=torch.bool,
            device=spectra.device,
        )
        memory_key_padding_mask = torch.cat(
            [spectrum_token_mask, peak_padding_mask],
            dim=1,
        )
        peaks = self.peak_encoder(spectra)
        latent_spectra = self.latent_spectrum.expand(peaks.shape[0], -1, -1)
        peaks = torch.cat([latent_spectra, peaks], dim=1)
        memory = self.spectrum_transformer_encoder(
            peaks,
            src_key_padding_mask=memory_key_padding_mask,
        )
        encoded_tokens = self.aa_encoder(partial_sequence_tokens)
        masses = self.mass_encoder(precursors[:, None, 0])
        charges = self.charge_encoder(precursors[:, 1].to(torch.int64) - 1)
        encoded_precursors = masses + charges[:, None, :]
        target = torch.cat([encoded_precursors, encoded_tokens], dim=1)
        target_key_padding_mask = target.sum(dim=2) == 0
        target = self.peptide_positional_encoder(target)
        target_length = target.shape[1]
        target_mask = ~torch.triu(
            torch.ones(
                (target_length, target_length),
                dtype=torch.bool,
                device=target.device,
            )
        ).transpose(0, 1)
        decoded = self.transformer_decoder(
            tgt=target,
            memory=memory,
            tgt_mask=target_mask,
            tgt_key_padding_mask=target_key_padding_mask,
            memory_key_padding_mask=memory_key_padding_mask.to(target.device),
        )
        token_logits = self.final(decoded)
        spectrum_representation = memory[:, 0]
        precursor_prediction = self.prec_layer(spectrum_representation)[:, 0]
        fragment_logits = self.frag_layer(memory)
        return token_logits, precursor_prediction, fragment_logits


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "cascadia"))

    from cascadia.depthcharge.encoders.sinusoidal import FloatEncoder, PositionalEncoder
    from cascadia.depthcharge.transformers.peptides import PeptideTransformerDecoder
    from cascadia.model import AugmentedSpec2Pep

    def float_encoder_forward(self: FloatEncoder, x: torch.Tensor) -> torch.Tensor:
        sin_mz = torch.sin(x[:, :, None] / self.sin_term)
        cos_mz = torch.cos(x[:, :, None] / self.cos_term)
        return torch.cat([sin_mz, cos_mz], dim=-1)

    def positional_encoder_forward(
        self: PositionalEncoder,
        x: torch.Tensor,
    ) -> torch.Tensor:
        positions = torch.arange(x.shape[1], dtype=self.sin_term.dtype, device=x.device)
        positions = positions.unsqueeze(0).expand(x.shape[0], -1)
        sin_in = positions.unsqueeze(2).expand(-1, -1, self.sin_term.shape[0])
        cos_in = positions.unsqueeze(2).expand(-1, -1, self.cos_term.shape[0])
        sin_pos = torch.sin(sin_in / self.sin_term)
        cos_pos = torch.cos(cos_in / self.cos_term)
        return torch.cat([sin_pos, cos_pos], dim=2) + x

    FloatEncoder.forward = float_encoder_forward
    PositionalEncoder.forward = positional_encoder_forward

    def decoder_forward(
        self: PeptideTransformerDecoder,
        tokens: torch.Tensor,
        precursors: torch.Tensor,
        memory: torch.Tensor,
        memory_key_padding_mask: torch.Tensor,
    ) -> torch.Tensor:
        encoded_tokens = self.aa_encoder(tokens)
        masses = self.mass_encoder(precursors[:, None, 0])
        charges = self.charge_encoder(precursors[:, 1].to(torch.int64) - 1)
        encoded_precursors = masses + charges[:, None, :]
        target = torch.cat([encoded_precursors, encoded_tokens], dim=1)
        target_key_padding_mask = target.sum(dim=2) == 0
        target = self.positional_encoder(target)
        target_length = target.shape[1]
        target_mask = ~torch.triu(
            torch.ones(
                (target_length, target_length),
                dtype=torch.bool,
                device=target.device,
            )
        ).transpose(0, 1)
        decoded = self.transformer_decoder(
            tgt=target,
            memory=memory,
            tgt_mask=target_mask,
            tgt_key_padding_mask=target_key_padding_mask,
            memory_key_padding_mask=memory_key_padding_mask.to(target.device),
        )
        return self.final(decoded)

    PeptideTransformerDecoder.forward = decoder_forward

    model = AugmentedSpec2Pep.load_from_checkpoint(
        str(args.checkpoint),
        d_model=args.d_model,
        n_layers=args.n_layers,
        n_head=args.n_head,
        dim_feedforward=args.dim_feedforward,
        dropout=args.dropout,
        rt_width=args.rt_width,
        tokenizer=args.n_tokens,
        max_charge=args.max_charge,
        map_location=args.device,
    )
    model.eval().to(args.device)

    wrapper = CascadiaSequenceModule(model).eval().to(args.device)
    spectra = torch.zeros(
        (args.example_batch, args.example_peaks, 4),
        dtype=torch.float32,
        device=args.device,
    )
    precursors = torch.ones(
        (args.example_batch, 2),
        dtype=torch.float32,
        device=args.device,
    )
    precursors[:, 0] = 1000.0
    precursors[:, 1] = 2.0
    partial_tokens = torch.zeros(
        (args.example_batch, args.example_sequence_length),
        dtype=torch.long,
        device=args.device,
    )

    with torch.no_grad():
        traced = torch.jit.script(wrapper)
        traced = torch.jit.freeze(traced)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    traced.save(str(args.output))
    print(f"Saved TorchScript model to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
