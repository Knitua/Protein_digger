#!/usr/bin/env python3
"""Prepare and validate the current 14,625-protein B1 inference universe."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        token = header.split()[0]
        parts = token.split("|")
        accession = parts[1] if len(parts) >= 2 and parts[0] in {"sp", "tr"} else token
        records.append((accession, header, "".join(chunks)))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    flush()
    return records


def fasta_bytes(records: list[tuple[str, str, str]]) -> bytes:
    stream = io.StringIO()
    for _, header, sequence in records:
        stream.write(f">{header}\n")
        for start in range(0, len(sequence), 80):
            stream.write(sequence[start : start + 80] + "\n")
    return stream.getvalue().encode("utf-8")


def safe_write(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() == content:
            return
        raise FileExistsError(path)
    partial = path.with_name(path.name + ".partial")
    with partial.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-tsv", type=Path, required=True)
    parser.add_argument("--current-fasta", type=Path, required=True)
    parser.add_argument("--previous-tsv", type=Path, required=True)
    parser.add_argument("--new69-tsv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(f"output directory exists: {args.out}")

    current_rows = read_tsv(args.current_tsv)
    previous_rows = read_tsv(args.previous_tsv)
    new_rows = read_tsv(args.new69_tsv)
    current_set = {row["Entry"].strip() for row in current_rows}
    previous_set = {row["Entry"].strip() for row in previous_rows}
    new_set = {row["Entry"].strip() for row in new_rows}
    if len(current_rows) != 14625 or len(current_set) != 14625:
        raise ValueError("current B1 universe must contain 14,625 unique accessions")
    if len(previous_rows) != 14556 or len(previous_set) != 14556:
        raise ValueError("previous B1 universe must contain 14,556 unique accessions")
    if len(new_rows) != 69 or len(new_set) != 69:
        raise ValueError("new set must contain 69 unique accessions")
    if previous_set & new_set:
        raise ValueError("previous B1 and new69 overlap")
    if previous_set | new_set != current_set:
        raise ValueError("14,625 universe is not exactly 14,556 + 69")

    fasta_records = parse_fasta(args.current_fasta)
    fasta_set = {record[0] for record in fasta_records}
    if len(fasta_records) != 14625 or len(fasta_set) != 14625 or fasta_set != current_set:
        raise ValueError("current FASTA does not exactly match current TSV")
    sequence_by_accession = {record[0]: record for record in fasta_records}
    new_records = [sequence_by_accession[row["Entry"].strip()] for row in new_rows]

    old_records = [record for record in fasta_records if record[0] in previous_set]
    old_records.sort(key=lambda record: (len(record[2]), record[0]))
    reference_indices = sorted({round(i * (len(old_records) - 1) / 31) for i in range(32)})
    reference_records = [old_records[index] for index in reference_indices]
    if len(reference_records) != 32:
        raise ValueError("reference selection did not yield 32 records")

    safe_write(args.out / "input/new69_canonical.fasta", fasta_bytes(new_records))
    safe_write(args.out / "input/reference32_existing.fasta", fasta_bytes(reference_records))
    safe_write(
        args.out / "input/new69_accessions.txt",
        ("\n".join(row["Entry"].strip() for row in new_rows) + "\n").encode("utf-8"),
    )
    manifest_rows = [
        {
            "role": "INPUT",
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in [args.current_tsv, args.current_fasta, args.previous_tsv, args.new69_tsv]
    ]
    manifest_rows.extend(
        {
            "role": "OUTPUT",
            "path": str(path.relative_to(args.out)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in [
            args.out / "input/new69_canonical.fasta",
            args.out / "input/reference32_existing.fasta",
            args.out / "input/new69_accessions.txt",
        ]
    )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=["role", "path", "bytes", "sha256"], delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(manifest_rows)
    safe_write(args.out / "audit/input_manifest.tsv", stream.getvalue().encode("utf-8"))
    summary = (
        "metric\tvalue\n"
        "current_common_input\t14625\n"
        "previous_validated_input\t14556\n"
        "new_isoform_aware_canonical_input\t69\n"
        "set_relation\tcurrent=previous_disjoint_union_new69\n"
        "reference_existing_proteins\t32\n"
    )
    safe_write(args.out / "audit/input_counts.tsv", summary.encode("utf-8"))
    safe_write(args.out / "scripts/prepare_b1_14625.py", Path(__file__).read_bytes())
    print("PASS\t14625 = 14556 disjoint-union 69")
    print(f"output\t{args.out}")


if __name__ == "__main__":
    main()
