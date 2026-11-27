#!/usr/bin/env python3
"""Export Cascadia's Lightning checkpoint to a LibTorch-loadable file.

Run from the repository root with a Python environment that has Cascadia,
PyTorch, and PyTorch-Lightning installed:

    python export_cascadia_torchscript.py --checkpoint cascadia.ckpt --output cascadia_sequence.pt
"""

from __future__ import annotations

import argparse
from pathlib import Path
import json
import hashlib
import os
from typing import Optional, Tuple
import copy
import math

os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(__file__).resolve().parents[1] / ".numba_cache"))

import torch

_torch_load = torch.load


def _torch_load_trusted_checkpoint(*args, **kwargs):
    if "weights_only" in __import__("inspect").signature(_torch_load).parameters:
        kwargs.setdefault("weights_only", False)
    return _torch_load(*args, **kwargs)


torch.load = _torch_load_trusted_checkpoint


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--d-model", type=int)
    parser.add_argument("--n-layers", type=int)
    parser.add_argument("--n-head", type=int)
    parser.add_argument("--dim-feedforward", type=int)
    parser.add_argument("--dropout", default=0.0, type=float)
    parser.add_argument("--rt-width", default=2.0, type=float)
    parser.add_argument("--max-charge", type=int)
    parser.add_argument("--n-tokens", type=int)
    parser.add_argument("--example-batch", default=2, type=int)
    parser.add_argument("--example-peaks", default=16, type=int)
    parser.add_argument("--example-sequence-length", default=4, type=int)
    parser.add_argument("--modifications", type=Path, help="Training PTMs as a JSON list of ProForma peptides")
    parser.add_argument("--vocabulary", type=Path, help="Exact training vocabulary JSON: ordered tokens (padding=null) and masses mapping")
    parser.add_argument("--freeze", action=argparse.BooleanOptionalAction, default=False, help="Freeze/in-line the scripted graph (higher loader memory); default --no-freeze keeps compact module methods")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


class StableAttention(torch.nn.MultiheadAttention):
    """Use the public functional attention path, matching Cascadia predict_step.

    PyTorch 2.0's eval fast path can differ from the grad-enabled path used by
    Cascadia. Calling the functional implementation also avoids its repeated
    mask-conversion warnings in a frozen C++ module.
    """
    def forward(self, query: torch.Tensor, key: torch.Tensor, value: torch.Tensor,
                key_padding_mask: Optional[torch.Tensor] = None,
                need_weights: bool = True, attn_mask: Optional[torch.Tensor] = None,
                average_attn_weights: bool = True, is_causal: bool = False
                ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        result, weights = torch.nn.functional.multi_head_attention_forward(
            query.transpose(0, 1), key.transpose(0, 1), value.transpose(0, 1),
            self.embed_dim, self.num_heads, self.in_proj_weight, self.in_proj_bias,
            self.bias_k, self.bias_v, self.add_zero_attn, 0.0,
            self.out_proj.weight, self.out_proj.bias, training=False,
            key_padding_mask=key_padding_mask, need_weights=need_weights,
            attn_mask=attn_mask, average_attn_weights=average_attn_weights,
            is_causal=is_causal)
        return result.transpose(0, 1), weights


def stable_attention(model):
    for name, child in list(model.named_children()):
        if isinstance(child, torch.nn.MultiheadAttention):
            replacement = StableAttention(child.embed_dim, child.num_heads, dropout=0.0, batch_first=True)
            replacement.load_state_dict(child.state_dict(), strict=True)
            setattr(model, name, replacement.eval())
        else:
            stable_attention(child)
        if isinstance(child, torch.nn.TransformerEncoderLayer):
            child.activation_relu_or_gelu = 0

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

    @torch.jit.export
    def encode_spectra(self, spectra: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        peak_padding_mask = ~spectra.sum(dim=2).to(torch.bool)
        mask = torch.cat([torch.zeros((spectra.shape[0], 1), dtype=torch.bool, device=spectra.device), peak_padding_mask], dim=1)
        peaks = self.peak_encoder(spectra)
        peaks = torch.cat([self.latent_spectrum.expand(peaks.shape[0], -1, -1), peaks], dim=1)
        return self.spectrum_transformer_encoder(peaks, src_key_padding_mask=mask), mask

    @torch.jit.export
    def decode_tokens(self, memory: torch.Tensor, memory_key_padding_mask: torch.Tensor, precursors: torch.Tensor, partial_sequence_tokens: torch.Tensor) -> torch.Tensor:
        encoded_tokens = self.aa_encoder(partial_sequence_tokens)
        masses = self.mass_encoder(precursors[:, None, 0])
        charges = self.charge_encoder(precursors[:, 1].to(torch.int64) - 1)
        target = torch.cat([masses + charges[:, None, :], encoded_tokens], dim=1)
        padding_mask = target.sum(dim=2) == 0
        target = self.peptide_positional_encoder(target)
        length = target.shape[1]
        mask = ~torch.triu(torch.ones((length, length), dtype=torch.bool, device=target.device)).transpose(0, 1)
        return self.final(self.transformer_decoder(tgt=target, memory=memory, tgt_mask=mask, tgt_key_padding_mask=padding_mask, memory_key_padding_mask=memory_key_padding_mask))

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
    torch.set_num_threads(args.threads)
    torch.manual_seed(0)
    from cascadia.depthcharge.primitives import Peptide
    from cascadia.depthcharge.tokenizers import PeptideTokenizer
    original_split = Peptide.split
    def compatible_split(self):
        if self.modifications is not None:
            self.modifications = [mods or None for mods in self.modifications]
        return original_split(self)
    Peptide.split = compatible_split

    checkpoint = torch.load(str(args.checkpoint), map_location="cpu")
    state = checkpoint["state_dict"]
    hp = checkpoint.get("hyper_parameters", {})
    inferred = {
        "d_model": state["decoder.aa_encoder.weight"].shape[1],
        "n_layers": len({k.split(".")[3] for k in state if k.startswith("decoder.transformer_decoder.layers.")}),
        "dim_feedforward": state["decoder.transformer_decoder.layers.0.linear1.weight"].shape[0],
        "max_charge": state["decoder.charge_encoder.weight"].shape[0] - 1,
        "n_tokens": state["decoder.aa_encoder.weight"].shape[0] - 1,
    }
    for name, value in inferred.items():
        supplied = getattr(args, name)
        if supplied is not None and supplied != value:
            raise ValueError(f"--{name.replace('_', '-')}={supplied} disagrees with checkpoint ({value})")
        setattr(args, name, value)
    if args.n_head is None:
        args.n_head = hp.get("n_head", 8)
    if args.n_head < 1 or args.d_model % args.n_head:
        raise ValueError("Attention heads must divide the embedding width")
    if args.vocabulary:
        vocab = json.loads(args.vocabulary.read_text())
        tokens, masses = vocab["tokens"], vocab["masses"]
        if tokens[0] is not None or len(set(tokens)) != len(tokens) or tokens.count("$") != 1:
            raise ValueError("Vocabulary must start with null padding and contain unique tokens and one '$'")
    else:
        tokenizer = hp.get("tokenizer")
        if args.modifications:
            tokenizer = PeptideTokenizer.from_proforma(json.loads(args.modifications.read_text()), reverse=False, replace_isoleucine_with_leucine=True)
        elif tokenizer is None or isinstance(tokenizer, int):
            if args.n_tokens != 28:
                raise ValueError("Nonstandard vocabulary: supply --vocabulary or --modifications from training; token meanings cannot be inferred from weight shapes")
            tokenizer = PeptideTokenizer.from_massivekb(reverse=False, replace_isoleucine_with_leucine=True)
        if getattr(tokenizer, "reverse", False):
            raise ValueError("Reverse tokenizers require an explicit forward-order vocabulary export; this runtime supports Cascadia forward tokenizers")
        tokens, masses = tokenizer.reverse_index, dict(tokenizer.residues)
    if any(not isinstance(t, str) or not t or any(c in t for c in "\t\n\r") for t in tokens[1:]):
        raise ValueError("Vocabulary tokens must be nonempty strings without tabs/newlines")
    if len(tokens) != args.n_tokens + 1:
        raise ValueError(f"Tokenizer has {len(tokens)} entries, checkpoint requires {args.n_tokens + 1}")
    if any(not math.isfinite(float(m)) for m in masses.values()):
        raise ValueError("Vocabulary masses must be finite")
    if any(t not in masses for t in tokens[1:] if t != "$"):
        raise ValueError("Every residue/modification token needs its training mass")
    print("Architecture:", inferred, "attention heads:", args.n_head, flush=True)
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
    model.spectrum_encoder.transformer_encoder.enable_nested_tensor = False
    model.spectrum_encoder.transformer_encoder.use_nested_tensor = False
    model.eval().to(args.device)
    del checkpoint, state

    reference = copy.deepcopy(model)
    stable_attention(model)
    wrapper = CascadiaSequenceModule(model).eval().to(args.device)
    with torch.no_grad():
        traced = torch.jit.script(wrapper)
        if args.freeze:
            traced = torch.jit.freeze(traced, preserved_attrs=["encode_spectra", "decode_tokens"])

    # Compare to the original eager forward over empty and nonempty token prefixes.
    # This catches exporter changes that silently change the model computation.
    with torch.no_grad():
        for batch, peaks, length in [(1, 7, 0), (args.example_batch, args.example_peaks, args.example_sequence_length), (1, 19, 1)]:
            spectra = torch.rand(batch, peaks, 4, device=args.device)
            spectra[:, :, 0] *= 1000
            spectra[:, :, 3] = 2
            if peaks > 7:
                spectra[:, -2:] = 0
            precursors = torch.tensor([[2000., 2.]], device=args.device).expand(batch, -1)
            prefix = torch.full((batch, length), 2, dtype=torch.long, device=args.device)
            with torch.enable_grad():
                expected = tuple(x.detach() for x in reference._forward_step(spectra, precursors, prefix))
            actual = traced(spectra, precursors, prefix)
            for a, b in zip(actual, expected):
                torch.testing.assert_close(a, b, rtol=2e-4, atol=2e-4)
            memory, mask = traced.encode_spectra(spectra)
            torch.testing.assert_close(traced.decode_tokens(memory, mask, precursors, prefix), expected[0], rtol=2e-4, atol=2e-4)
    metadata = dict(format="cascadia-bundle-v1", frozen=args.freeze, architecture=inferred, n_head=args.n_head,
                    rt_width=args.rt_width, torch_version=torch.__version__,
                    checkpoint=str(args.checkpoint.name), tokenizer_direction="forward",
                    preprocessing="cascadia-0.0.7", vocabulary_source="explicit vocabulary" if args.vocabulary else "custom PTMs" if args.modifications else "checkpoint tokenizer" if hp.get("tokenizer") is not None else "canonical MassiveKB (public Cascadia models)", parity="passed: 3 variable input shapes")
    digest = hashlib.sha256()
    with args.checkpoint.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    metadata["checkpoint_sha256"] = digest.hexdigest()
    vocabulary = "#cascadia-bundle-v1\n" + "".join(f"{t}\t{float(masses.get(t, 0)):.15g}\n" for t in tokens[1:])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    extra = {"cascadia.vocab.tsv": vocabulary, "cascadia.metadata.json": json.dumps(metadata, indent=2)}
    torch.jit.save(traced, str(args.output), _extra_files=extra)
    restored = torch.jit.load(str(args.output), map_location=args.device)
    with torch.no_grad():
        torch.testing.assert_close(restored(spectra, precursors, prefix)[0], actual[0])
    args.output.with_suffix(args.output.suffix + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Saved validated TorchScript bundle to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
