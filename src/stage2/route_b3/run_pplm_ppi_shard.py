#!/usr/bin/env python3
"""Resumable PPLM-PPI scoring for PLM-interact-selected B3 pairs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import tempfile
import time
from pathlib import Path

import torch


OUTPUT_FIELDS = [
    "pair_index", "candidate_index", "anchor_index", "isoform_accession",
    "canonical_accession", "candidate_gene", "candidate_length",
    "anchor_accession", "anchor_gene", "anchor_length", "aa_length_sum",
    "plm_interact_probability", "pplm_ppi_probability", "pplm_gt_0p9",
    "shard_id", "num_shards", "host", "gpu_name",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def clean_sequence(value: str) -> str:
    return "".join(value.split()).upper()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temp_name = handle.name
    os.replace(temp_name, path)


def already_processed(path: Path) -> set[int]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != OUTPUT_FIELDS:
            raise ValueError(f"Unexpected output header in {path}")
        return {int(row["pair_index"]) for row in reader}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-pairs", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--pplm-repo", type=Path, required=True)
    parser.add_argument("--pplm-model", type=Path, required=True)
    parser.add_argument("--ppi-models", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--progress", type=Path, required=True)
    parser.add_argument("--shard-id", type=int, required=True)
    parser.add_argument("--num-shards", type=int, required=True)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    if not (0 <= args.shard_id < args.num_shards):
        raise ValueError("shard-id must be in [0, num-shards)")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required")

    sys.path.insert(0, str(args.pplm_repo))
    from pplm import Alphabet, PPLM  # noqa: PLC0415
    from pplm_ppi import PPLM_PPI  # noqa: PLC0415

    device = torch.device("cuda:0")
    torch.manual_seed(2)
    torch.cuda.manual_seed_all(2)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False

    candidates = read_tsv(args.candidates)
    anchors = read_tsv(args.anchors)
    selected = read_tsv(args.selected_pairs)
    if len(candidates) != 497 or len(anchors) != 2368:
        raise ValueError(f"Unexpected frozen inputs: {len(candidates)} x {len(anchors)}")

    candidate_sequences = [clean_sequence(row["sequence"]) for row in candidates]
    anchor_sequences = [clean_sequence(row["sequence"]) for row in anchors]
    for index, (row, sequence) in enumerate(zip(candidates, candidate_sequences, strict=True)):
        if len(sequence) != int(row["length"]):
            raise ValueError(f"Candidate length mismatch at row {index}")
    for index, (row, sequence) in enumerate(zip(anchors, anchor_sequences, strict=True)):
        if len(sequence) != int(row["anchor_length"]):
            raise ValueError(f"Anchor length mismatch at row {index}")

    pairs = []
    seen = set()
    for row in selected:
        pair_index = int(row["pair_index"])
        if pair_index in seen:
            raise ValueError(f"Duplicate selected pair_index: {pair_index}")
        seen.add(pair_index)
        if float(row["plm_interact_probability"]) <= 0.99:
            raise ValueError(f"Selected pair does not satisfy PLM >0.99: {pair_index}")
        candidate_index, anchor_index = divmod(pair_index, len(anchors))
        if int(row["candidate_index"]) != candidate_index or int(row["anchor_index"]) != anchor_index:
            raise ValueError(f"Pair index mapping mismatch: {pair_index}")
        if int(row["aa_length_sum"]) > 1600:
            raise ValueError(f"Length-excluded pair present: {pair_index}")
        if pair_index % args.num_shards == args.shard_id:
            pairs.append((int(row["aa_length_sum"]), pair_index, candidate_index, anchor_index, row))
    pairs.sort()
    if args.limit is not None:
        pairs = pairs[: args.limit]

    processed = already_processed(args.output)
    pending = [item for item in pairs if item[1] not in processed]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_exists = args.output.exists() and args.output.stat().st_size > 0

    alphabet = Alphabet.from_architecture()
    batch_converter = alphabet.get_batch_converter()
    model_data = torch.load(args.pplm_model, map_location="cpu", weights_only=False)
    model_param = model_data["param"]
    model = PPLM(
        num_layers=model_param["encoder_layers"],
        embed_dim=model_param["encoder_embed_dim"],
        attention_heads=model_param["encoder_attention_heads"],
        token_dropout=False,
        alphabet=alphabet,
    )
    incompatible = model.load_state_dict(model_data["model"], strict=False)
    model.eval().to(device)

    ensemble_weights = torch.load(args.ppi_models, map_location="cpu", weights_only=False)
    mean_models = []
    max_models = []
    for state in ensemble_weights["mean"]:
        classifier = PPLM_PPI()
        classifier.load_state_dict(state, strict=True)
        classifier.eval().to(device)
        mean_models.append(classifier)
    for state in ensemble_weights["max"]:
        classifier = PPLM_PPI()
        classifier.load_state_dict(state, strict=True)
        classifier.eval().to(device)
        max_models.append(classifier)

    host = os.uname().nodename
    gpu_name = torch.cuda.get_device_name(0)
    start = time.time()
    completed_this_run = 0
    high_this_run = 0

    with args.output.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        if not output_exists:
            writer.writeheader()
            handle.flush()

        for aa_sum, pair_index, candidate_index, anchor_index, source_row in pending:
            seq_a = candidate_sequences[candidate_index]
            seq_b = anchor_sequences[anchor_index]
            _, _, tokens_a = batch_converter([("seqA", seq_a)])
            _, _, tokens_b = batch_converter([("seqB", seq_b)])
            tokens = torch.cat([tokens_a, tokens_b], dim=-1).to(device)
            token_length = int(tokens.shape[-1])
            if token_length != aa_sum + 4:
                raise RuntimeError(f"Unexpected PPLM token length for pair {pair_index}")

            inter_chain_mask = torch.ones((token_length, token_length), device=device)
            split = len(seq_a) + 2
            inter_chain_mask[:split, :split] = 0
            inter_chain_mask[split:, split:] = 0

            with torch.inference_mode():
                out = model(
                    tokens,
                    inter_chain_mask,
                    repr_layers=[33],
                    need_head_weights=True,
                    return_contacts=False,
                )
                embed_a = out["representations"][33][0, 1 : len(seq_a) + 1, :]
                embed_b = out["representations"][33][0, -(len(seq_b) + 1) : -1, :]
                attentions = out["attentions"].squeeze(0).reshape(33 * 20, token_length, token_length)
                attn_aa = attentions[:, 1 : len(seq_a) + 1, 1 : len(seq_a) + 1]
                attn_ab = attentions[:, 1 : len(seq_a) + 1, -(len(seq_b) + 1) : -1]
                attn_ba = attentions[:, -(len(seq_b) + 1) : -1, 1 : len(seq_a) + 1]
                attn_bb = attentions[:, -(len(seq_b) + 1) : -1, -(len(seq_b) + 1) : -1]
                inter_attn = (attn_ab + attn_ba.transpose(1, 2)) / 2

                mean_features = (
                    inter_attn.mean(dim=(1, 2)),
                    attn_aa.mean(dim=(1, 2)),
                    attn_bb.mean(dim=(1, 2)),
                    embed_a.mean(dim=0),
                    embed_b.mean(dim=0),
                )
                max_features = (
                    torch.amax(inter_attn, dim=(1, 2)),
                    torch.amax(attn_aa, dim=(1, 2)),
                    torch.amax(attn_bb, dim=(1, 2)),
                    torch.amax(embed_a, dim=0),
                    torch.amax(embed_b, dim=0),
                )
                predictions = []
                for classifier in mean_models:
                    forward = classifier(*mean_features)
                    reverse = classifier(
                        mean_features[0], mean_features[2], mean_features[1],
                        mean_features[4], mean_features[3],
                    )
                    predictions.append((forward + reverse) / 2)
                for classifier in max_models:
                    forward = classifier(*max_features)
                    reverse = classifier(
                        max_features[0], max_features[2], max_features[1],
                        max_features[4], max_features[3],
                    )
                    predictions.append((forward + reverse) / 2)
                score = float(torch.stack(predictions).mean().item())

            is_high = int(score > 0.9)
            high_this_run += is_high
            candidate = candidates[candidate_index]
            anchor = anchors[anchor_index]
            writer.writerow({
                "pair_index": pair_index,
                "candidate_index": candidate_index,
                "anchor_index": anchor_index,
                "isoform_accession": candidate["isoform_accession"],
                "canonical_accession": candidate["canonical_accession"],
                "candidate_gene": candidate["gene_primary"],
                "candidate_length": candidate["length"],
                "anchor_accession": anchor["anchor_accession"],
                "anchor_gene": anchor["anchor_gene"],
                "anchor_length": anchor["anchor_length"],
                "aa_length_sum": aa_sum,
                "plm_interact_probability": source_row["plm_interact_probability"],
                "pplm_ppi_probability": f"{score:.9f}",
                "pplm_gt_0p9": is_high,
                "shard_id": args.shard_id,
                "num_shards": args.num_shards,
                "host": host,
                "gpu_name": gpu_name,
            })
            handle.flush()
            completed_this_run += 1
            elapsed = time.time() - start
            rate = completed_this_run / elapsed if elapsed > 0 else 0.0
            remaining = len(pending) - completed_this_run
            write_json_atomic(args.progress, {
                "status": "running" if remaining else "complete",
                "shard_id": args.shard_id,
                "num_shards": args.num_shards,
                "host": host,
                "gpu_name": gpu_name,
                "pairs_in_shard": len(pairs),
                "already_processed_before_run": len(processed),
                "completed_this_run": completed_this_run,
                "completed_total": len(processed) + completed_this_run,
                "remaining": remaining,
                "pplm_gt_0p9_this_run": high_this_run,
                "last_token_length": token_length,
                "elapsed_seconds": elapsed,
                "pairs_per_second": rate,
                "naive_eta_seconds": remaining / rate if rate > 0 else None,
                "pplm_model_missing_keys": list(incompatible.missing_keys),
                "pplm_model_unexpected_keys": list(incompatible.unexpected_keys),
                "ensemble_mean_models": len(mean_models),
                "ensemble_max_models": len(max_models),
                "output": str(args.output),
            })


if __name__ == "__main__":
    main()
