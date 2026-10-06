#!/usr/bin/env python3
"""Read-only validation of the five primary 02_anchors datasets."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def accession_set(rows: list[dict[str, str]]) -> set[str]:
    if not rows:
        return set()
    key = "accession" if "accession" in rows[0] else "anchor_accession"
    return {row[key].strip() for row in rows}


def require_dataset(path: Path, expected: int) -> list[dict[str, str]]:
    rows = read_tsv(path)
    accessions = accession_set(rows)
    if len(rows) != expected or len(accessions) != expected or "" in accessions:
        raise ValueError(
            f"{path}: rows={len(rows)}, unique_nonblank_accessions={len(accessions - {''})}, expected={expected}"
        )
    return rows


def fasta_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    accession: str | None = None
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                accession = line[1:].split("|")[0].split()[0]
                if accession in records:
                    raise ValueError(f"duplicate FASTA accession in {path}: {accession}")
                records[accession] = []
            else:
                if accession is None:
                    raise ValueError(f"sequence before FASTA header in {path}")
                records[accession].append(line)
    return {key: "".join(parts) for key, parts in records.items()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-root", type=Path, required=True)
    parser.add_argument("--reference-metadata", type=Path, required=True)
    args = parser.parse_args()

    root = args.anchor_root
    source = root / "01_source_filtered"
    merged = root / "02_merged"

    lambert = require_dataset(source / "lambert_tf_filtered_1632.tsv", 1632)
    epi = require_dataset(source / "epifactors_filtered_794.tsv", 794)
    animal = require_dataset(source / "animaltfdb_accepted_72.tsv", 72)
    lambert_animal = require_dataset(
        merged / "lambert_animaltfdb_deduplicated_1696.tsv", 1696
    )
    all_three = require_dataset(
        merged / "lambert_epifactors_animaltfdb_deduplicated_2368.tsv", 2368
    )
    stable = require_dataset(root / "working_regulatory_anchor_2368.tsv", 2368)

    l_set = accession_set(lambert)
    e_set = accession_set(epi)
    a_set = accession_set(animal)
    la_set = accession_set(lambert_animal)
    all_set = accession_set(all_three)

    if Counter(row["selection_status"] for row in animal) != {
        "ACCEPTED_STRICT_TF": 37,
        "ACCEPTED_REGULATORY_EXTENSION": 35,
    }:
        raise ValueError("AnimalTFDB 72 selection_status must be 37 strict + 35 regulatory")
    if Counter(row["target_dataset_role"] for row in animal) != {
        "STRICT_TF": 37,
        "REGULATORY_EXTENSION_NOT_TF_OR_EPIFACTOR": 35,
    }:
        raise ValueError("AnimalTFDB 72 target_dataset_role must preserve both evidence levels")

    if len(l_set & a_set) != 8:
        raise ValueError("Lambert/AnimalTFDB accepted overlap must be 8 proteins")
    if len(l_set & e_set) != 122 or a_set & e_set:
        raise ValueError("expected Lambert/EpiFactors overlap=122 and Animal/EpiFactors overlap=0")
    if la_set != l_set | a_set:
        raise ValueError("1,696 table is not the exact Lambert/AnimalTFDB accession union")
    if all_set != l_set | a_set | e_set:
        raise ValueError("2,368 table is not the exact three-source accession union")
    if accession_set(stable) != all_set:
        raise ValueError("stable 2,368 interface differs from the three-source union")

    expected_roles = {
        "is_tf_anchor": {"1": 1661, "0": 707},
        "is_epifactor_anchor": {"1": 794, "0": 1574},
        "is_regulatory_extension": {"1": 35, "0": 2333},
    }
    for field, expected in expected_roles.items():
        observed = dict(Counter(row[field] for row in all_three))
        if observed != expected:
            raise ValueError(f"unexpected {field} composition: {observed}")

    reference = read_tsv(args.reference_metadata)
    ref_sequences = {row["Entry"]: row["Sequence"] for row in reference}
    if len(reference) != 20416 or len(ref_sequences) != 20416:
        raise ValueError("reference is not the frozen 20,416-protein universe")
    for rows in (lambert, epi, animal, lambert_animal, all_three, stable):
        acc_key = "accession" if "accession" in rows[0] else "anchor_accession"
        for row in rows:
            accession = row[acc_key]
            if ref_sequences.get(accession) != row["sequence"]:
                raise ValueError(f"canonical sequence mismatch: {accession}")

    fasta_expectations = {
        source / "lambert_tf_filtered_1632.fasta": 1632,
        source / "epifactors_filtered_794.fasta": 794,
        source / "animaltfdb_accepted_72.fasta": 72,
        merged / "lambert_animaltfdb_deduplicated_1696.fasta": 1696,
        merged / "lambert_epifactors_animaltfdb_deduplicated_2368.fasta": 2368,
        root / "working_regulatory_anchor_2368.fasta": 2368,
    }
    for path, expected in fasta_expectations.items():
        observed = fasta_count(path)
        if observed != expected:
            raise ValueError(f"{path}: FASTA records={observed}, expected={expected}")

    if sha256(merged / "lambert_epifactors_animaltfdb_deduplicated_2368.tsv") != sha256(
        root / "working_regulatory_anchor_2368.tsv"
    ):
        raise ValueError("stable TSV is not byte-identical to the merged 2,368 table")
    if read_fasta(merged / "lambert_epifactors_animaltfdb_deduplicated_2368.fasta") != read_fasta(
        root / "working_regulatory_anchor_2368.fasta"
    ):
        raise ValueError("stable FASTA differs from the merged 2,368 accession-sequence mapping")

    expected_source_files = {
        "README.md",
        "lambert_tf_filtered_1632.tsv",
        "lambert_tf_filtered_1632.fasta",
        "epifactors_filtered_794.tsv",
        "epifactors_filtered_794.fasta",
        "animaltfdb_accepted_72.tsv",
        "animaltfdb_accepted_72.fasta",
    }
    expected_merged_files = {
        "README.md",
        "lambert_animaltfdb_deduplicated_1696.tsv",
        "lambert_animaltfdb_deduplicated_1696.fasta",
        "lambert_epifactors_animaltfdb_deduplicated_2368.tsv",
        "lambert_epifactors_animaltfdb_deduplicated_2368.fasta",
    }
    if {p.name for p in source.iterdir() if p.is_file()} != expected_source_files:
        raise ValueError("01_source_filtered contains files outside the three primary source datasets")
    if {p.name for p in merged.iterdir() if p.is_file()} != expected_merged_files:
        raise ValueError("02_merged contains files outside the two primary merged datasets")

    print("five primary anchor datasets: PASS")
    print("Lambert=1632; EpiFactors=794; AnimalTFDB accepted=72 (37+35)")
    print("Lambert|AnimalTFDB=1696; all-three union=2368")


if __name__ == "__main__":
    main()
