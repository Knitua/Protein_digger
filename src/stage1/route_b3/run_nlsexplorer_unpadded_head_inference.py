#!/usr/bin/env python
import argparse
import csv
import sys
import time
from itertools import groupby
from pathlib import Path

import esm
import torch
from torch.nn.utils.rnn import pad_sequence

SHARED_MODEL_SCRIPT_ROOT = Path("/root/autodl-tmp/Agent_analysis/04_Nuclear_Models/scripts")
sys.path.insert(0, str(SHARED_MODEL_SCRIPT_ROOT))

from benchmark_nlsexplorer_runtime import LSTMTagger, make_windows, read_fasta


def predict_prob_only(model, embding):
    batch_size = embding.size(0)
    vec_store = []
    t_emb = embding

    for i, fig in enumerate(model.Att_config):
        vec_s = []
        att_s = []
        for k in range(fig):
            vec, att = model.Att_li[i][k](t_emb)
            vec_s.append(torch.relu(model.pro_li[i][k](vec)))
            att_s.append(att)

        att_s = torch.stack(att_s).squeeze(3)
        sum_att = model.project_li[i](att_s.sum(0).unsqueeze(2))
        t_emb = model.norm_li[i]((t_emb * sum_att) + t_emb)

        vec_s = torch.stack(vec_s).transpose(0, 1)
        z = vec_s.permute(0, 2, 1)
        output = z.reshape(z.size(0), -1)
        output = torch.relu(model.FI_li[i](output)) + output
        output = torch.relu(model.Set_li[i](output))
        vec_store.append(output)

    ott = torch.stack(vec_store).transpose(0, 1).reshape(batch_size, -1)
    summed, _ = model.AAt(t_emb)
    logits = model.hidden2p(model.dropout(torch.cat((summed, ott), 1)))
    return torch.sigmoid(logits).view(-1)


def format_seconds(seconds):
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m{s:02d}s"
    return f"{m}m{s:02d}s"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--nls-weights", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=1022)
    parser.add_argument("--overlap", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    windows_tsv = outdir / "nlsexplorer_window_scores.tsv"
    proteins_tsv = outdir / "nlsexplorer_protein_scores.tsv"
    manifest = outdir / "MANIFEST.md"

    records = read_fasta(args.fasta)
    protein_lengths = {acc: len(seq) for acc, _, seq in records}
    windows = make_windows(records, args.window, args.overlap)
    # The classifier head has no padding-mask argument. Build batches only from
    # windows of exactly the same length, so vectorized head inference is both
    # efficient and independent of unrelated batch neighbours.
    windows.sort(key=lambda row: (row["length"], row["accession"], row["start0"]))
    window_batches = []
    for _, same_length_iter in groupby(windows, key=lambda row: row["length"]):
        same_length = list(same_length_iter)
        window_batches.extend(
            same_length[start : start + args.batch_size]
            for start in range(0, len(same_length), args.batch_size)
        )

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    if device.type == "cuda":
        torch.cuda.set_device(device)
    dtype = torch.float16 if device.type == "cuda" else torch.bfloat16

    t0 = time.perf_counter()
    esm_model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    esm_model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()

    nls_model = LSTMTagger().eval()
    nls_model.load_state_dict(torch.load(args.nls_weights, map_location="cpu"), strict=True)
    nls_model.to(device)
    load_seconds = time.perf_counter() - t0

    protein_stats = {}
    processed = 0
    infer_start = time.perf_counter()

    with windows_tsv.open("w", newline="") as wf:
        writer = csv.writer(wf, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "window_id",
                "accession",
                "protein_length",
                "window_start_1based",
                "window_end_1based",
                "window_length",
                "nls_prob",
            ]
        )

        with torch.inference_mode():
            window_id = 0
            next_progress = args.progress_every
            for batch in window_batches:
                data = [(w["accession"], w["sequence"]) for w in batch]
                _, _, toks = batch_converter(data)
                toks = toks.to(device)

                with torch.autocast(
                    device_type=device.type, dtype=dtype, enabled=device.type == "cuda"
                ):
                    reps = esm_model(toks, repr_layers=[33])["representations"][33]
                    emb_list = [
                        reps[j, 1 : 1 + len(w["sequence"])] for j, w in enumerate(batch)
                    ]
                    # All lengths are equal by construction; pad_sequence is an
                    # exact stack here and adds no zero-padding residues.
                    unpadded_batch = pad_sequence(emb_list, batch_first=True)
                    probs = predict_prob_only(nls_model, unpadded_batch)

                probs_cpu = probs.detach().float().cpu().tolist()
                for offset, (w, prob) in enumerate(zip(batch, probs_cpu)):
                    acc = w["accession"]
                    row = [
                        window_id,
                        acc,
                        protein_lengths[acc],
                        w["start0"] + 1,
                        w["end0"],
                        w["length"],
                        f"{prob:.8f}",
                    ]
                    writer.writerow(row)
                    window_id += 1

                    stat = protein_stats.setdefault(
                        acc,
                        {
                            "n_windows": 0,
                            "sum_prob": 0.0,
                            "max_prob": -1.0,
                            "max_start": None,
                            "max_end": None,
                        },
                    )
                    stat["n_windows"] += 1
                    stat["sum_prob"] += prob
                    if prob > stat["max_prob"]:
                        stat["max_prob"] = prob
                        stat["max_start"] = w["start0"] + 1
                        stat["max_end"] = w["end0"]

                processed += len(batch)
                if processed >= next_progress or processed == len(windows):
                    elapsed = time.perf_counter() - infer_start
                    rate = processed / elapsed if elapsed else 0.0
                    remaining = (len(windows) - processed) / rate if rate else 0.0
                    print(
                        f"processed_windows={processed}/{len(windows)} "
                        f"rate={rate:.3f}/s eta={format_seconds(remaining)}",
                        flush=True,
                    )
                    while next_progress <= processed:
                        next_progress += args.progress_every

                del toks, reps, emb_list, unpadded_batch, probs
                if device.type == "cuda":
                    torch.cuda.empty_cache()

    if device.type == "cuda":
        torch.cuda.synchronize(device)
    inference_seconds = time.perf_counter() - infer_start
    total_seconds = time.perf_counter() - t0

    with proteins_tsv.open("w", newline="") as pf:
        writer = csv.writer(pf, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "accession",
                "protein_length",
                "n_windows",
                "nls_prob_max",
                "max_window_start_1based",
                "max_window_end_1based",
                "nls_prob_mean",
            ]
        )
        for acc, _, _ in records:
            stat = protein_stats[acc]
            writer.writerow(
                [
                    acc,
                    protein_lengths[acc],
                    stat["n_windows"],
                    f"{stat['max_prob']:.8f}",
                    stat["max_start"],
                    stat["max_end"],
                    f"{stat['sum_prob'] / stat['n_windows']:.8f}",
                ]
            )

    manifest.write_text(
        f"""# NLSExplorer Window Inference

Input FASTA: `{args.fasta}`

Settings:

```text
window = {args.window}
overlap = {args.overlap}
batch_size = {args.batch_size}
device = {device}
nls_head_evaluation = exact-length batches; no zero-padding is introduced
```

Counts:

```text
proteins = {len(records)}
windows = {len(windows)}
```

Runtime:

```text
model_load_seconds = {load_seconds:.3f}
inference_seconds = {inference_seconds:.3f}
total_seconds = {total_seconds:.3f}
```

Outputs:

- `nlsexplorer_window_scores.tsv`: one row per sequence window.
- `nlsexplorer_protein_scores.tsv`: one row per protein, using max window probability as the protein-level NLS score.
"""
    )

    print(f"proteins\t{len(records)}")
    print(f"windows\t{len(windows)}")
    print(f"load_seconds\t{load_seconds:.3f}")
    print(f"inference_seconds\t{inference_seconds:.3f}")
    print(f"total_seconds\t{total_seconds:.3f}")
    print(f"windows_tsv\t{windows_tsv}")
    print(f"proteins_tsv\t{proteins_tsv}")
    print(f"manifest\t{manifest}")


if __name__ == "__main__":
    main()
