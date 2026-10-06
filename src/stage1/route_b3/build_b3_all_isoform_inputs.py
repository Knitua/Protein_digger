#!/usr/bin/env python3
"""Build the full alternative-isoform input and direct UniProt nucleus evidence for B3."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
from collections import defaultdict
from pathlib import Path


TAB = "\t"
NUCLEUS_RE = re.compile(r"\bNucleus\b")
ISOFORM_LOCATION_RE = re.compile(r"^\[Isoform ([^\]]+)\]:")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=TAB))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=TAB, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_label(value: str) -> str:
    value = re.sub(r"\s*\{[^{}]*\}\s*", " ", value)
    return re.sub(r"\s+", " ", value.strip())


def normalized_label(value: str) -> str:
    return clean_label(value).casefold()


def parse_alternative_products(value: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for raw_token in value.split(";"):
        token = raw_token.strip()
        if token.startswith("Name="):
            if current is not None:
                records.append(current)
            current = {
                "name": clean_label(token[5:]),
                "synonyms": [],
                "isoform_ids": [],
                "displayed": False,
            }
        elif current is not None and token.startswith("Synonyms="):
            current["synonyms"] = [clean_label(item) for item in token[9:].split(",")]
        elif current is not None and token.startswith("IsoId="):
            current["isoform_ids"] = [item.strip() for item in token[6:].split(",")]
        elif current is not None and token.startswith("Sequence="):
            current["displayed"] = token[9:].strip() == "Displayed"
    if current is not None:
        records.append(current)
    return records


def location_sections(value: str) -> list[str]:
    value = re.sub(r"^SUBCELLULAR LOCATION:\s*", "", value)
    return re.split(r";\s*SUBCELLULAR LOCATION:\s*", value) if value else []


def formal_location_part(section: str) -> str:
    return section.split("Note=", 1)[0].strip()


def fasta_from_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f">{row['isoform_accession']}|canonical_parent={row['canonical_accession']}"
                f"|gene={row['gene_primary']}\n{row['sequence']}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alternative-metadata", type=Path, required=True)
    parser.add_argument("--alternative-fasta", type=Path, required=True)
    parser.add_argument("--alternative-mapping", type=Path, required=True)
    parser.add_argument("--frozen-annotations", type=Path, required=True)
    parser.add_argument("--alternative-products", type=Path, required=True)
    parser.add_argument("--legacy-strict-nucleus-accessions", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()

    out_root = args.out_root
    input_dir = out_root / "00_input"
    direct_dir = out_root / "01_nuclear_evidence/01_uniprot_direct"
    audit_dir = out_root / "03_audit"
    input_dir.mkdir(parents=True, exist_ok=True)
    direct_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    alt_rows = read_tsv(args.alternative_metadata)
    if len(alt_rows) != 22131:
        raise RuntimeError(f"expected 22,131 alternative isoforms, found {len(alt_rows)}")
    alt_by_id = {row["isoform_accession"]: row for row in alt_rows}
    if len(alt_by_id) != 22131:
        raise RuntimeError("alternative isoform accessions are not unique")

    annotation_rows = read_tsv(args.frozen_annotations)
    annotation_by_accession = {row["Entry"]: row for row in annotation_rows}
    products_rows = read_tsv(args.alternative_products)
    products_by_accession = {row["Entry"]: row for row in products_rows}
    if len(annotation_by_accession) != 20416 or len(products_by_accession) != 20416:
        raise RuntimeError("canonical annotation/product resources do not each cover 20,416 accessions")

    legacy_strict = {
        line.strip()
        for line in args.legacy_strict_nucleus_accessions.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if len(legacy_strict) != 5742:
        raise RuntimeError(f"expected 5,742 legacy strict-nucleus accessions, found {len(legacy_strict)}")

    direct_evidence: dict[str, dict[str, object]] = {}
    canonical_nucleus_status: dict[str, bool] = {}
    parent_direct_isoforms: dict[str, set[str]] = defaultdict(set)

    for accession, annotation in annotation_by_accession.items():
        product_row = products_by_accession.get(accession)
        product_records = parse_alternative_products(
            product_row["Alternative products (isoforms)"] if product_row is not None else ""
        )
        label_lookup: dict[str, list[dict[str, object]]] = defaultdict(list)
        displayed_labels: set[str] = set()
        for product in product_records:
            labels = [str(product["name"]), *[str(x) for x in product["synonyms"]]]
            for label in labels:
                label_lookup[normalized_label(label)].append(product)
                if bool(product["displayed"]):
                    displayed_labels.add(normalized_label(label))

        canonical_has_nucleus = False
        for section in location_sections(annotation["Subcellular location [CC]"]):
            formal = formal_location_part(section)
            if not NUCLEUS_RE.search(formal):
                continue
            tagged = ISOFORM_LOCATION_RE.match(section)
            if tagged is None:
                canonical_has_nucleus = True
                continue
            label = tagged.group(1)
            norm_label = normalized_label(label)
            if norm_label in displayed_labels:
                canonical_has_nucleus = True
            products = label_lookup.get(norm_label, [])
            if not products:
                raise RuntimeError(f"cannot map nuclear isoform label: {accession} / {label}")
            for product in products:
                if bool(product["displayed"]):
                    continue
                for isoform_id in product["isoform_ids"]:
                    isoform_id = str(isoform_id)
                    if isoform_id not in alt_by_id:
                        raise RuntimeError(f"nuclear isoform missing from alternative metadata: {isoform_id}")
                    parent_direct_isoforms[accession].add(isoform_id)
                    evidence = direct_evidence.setdefault(
                        isoform_id,
                        {
                            "isoform_accession": isoform_id,
                            "canonical_accession": accession,
                            "gene_primary": alt_by_id[isoform_id]["gene_primary"],
                            "isoform_number": alt_by_id[isoform_id]["isoform_number"],
                            "isoform_label": label,
                            "isoform_protein_name": alt_by_id[isoform_id]["isoform_protein_name"],
                            "length": alt_by_id[isoform_id]["length"],
                            "uniprot_location_sections": [],
                        },
                    )
                    evidence["uniprot_location_sections"].append(formal)
        canonical_nucleus_status[accession] = canonical_has_nucleus

    corrected_parents = {
        accession
        for accession in legacy_strict
        if not canonical_nucleus_status[accession] and parent_direct_isoforms.get(accession)
    }
    if len(direct_evidence) != 266 or len(parent_direct_isoforms) != 203:
        raise RuntimeError(
            "unexpected direct isoform-nucleus counts: "
            f"isoforms={len(direct_evidence)}, parents={len(parent_direct_isoforms)}"
        )
    if len(corrected_parents) != 73:
        raise RuntimeError(f"expected 73 corrected parents, found {len(corrected_parents)}")
    corrected_direct_count = sum(
        len(parent_direct_isoforms[accession]) for accession in corrected_parents
    )
    if corrected_direct_count != 82:
        raise RuntimeError(
            f"expected 82 direct nuclear isoforms among corrected parents, found {corrected_direct_count}"
        )

    frozen_only_accessions = sorted(set(annotation_by_accession) - set(products_by_accession))
    api_only_accessions = sorted(set(products_by_accession) - set(annotation_by_accession))

    all_tsv = input_dir / "b3_all_alternative_isoforms_22131.tsv"
    all_fasta = input_dir / "b3_all_alternative_isoforms_22131.fasta"
    all_mapping = input_dir / "b3_all_alternative_isoforms_to_canonical_22131.tsv"
    shutil.copy2(args.alternative_metadata, all_tsv)
    shutil.copy2(args.alternative_fasta, all_fasta)
    shutil.copy2(args.alternative_mapping, all_mapping)
    (input_dir / "b3_all_alternative_isoforms_22131.accessions.txt").write_text(
        "".join(f"{row['isoform_accession']}\n" for row in alt_rows), encoding="utf-8"
    )

    direct_rows: list[dict[str, object]] = []
    for isoform_id in sorted(direct_evidence):
        evidence = direct_evidence[isoform_id]
        parent = str(evidence["canonical_accession"])
        direct_rows.append(
            {
                **{key: value for key, value in evidence.items() if key != "uniprot_location_sections"},
                "parent_canonical_isoform_aware_nucleus": int(canonical_nucleus_status[parent]),
                "legacy_parent_in_strict_nucleus_5742": int(parent in legacy_strict),
                "corrected_parent_group_73": int(parent in corrected_parents),
                "uniprot_direct_nucleus": 1,
                "uniprot_location_sections": " || ".join(evidence["uniprot_location_sections"]),
            }
        )
    direct_fields = [
        "isoform_accession",
        "canonical_accession",
        "gene_primary",
        "isoform_number",
        "isoform_label",
        "isoform_protein_name",
        "length",
        "parent_canonical_isoform_aware_nucleus",
        "legacy_parent_in_strict_nucleus_5742",
        "corrected_parent_group_73",
        "uniprot_direct_nucleus",
        "uniprot_location_sections",
    ]
    write_tsv(direct_dir / "b3_uniprot_isoform_specific_nucleus_266.tsv", direct_fields, direct_rows)
    (direct_dir / "b3_uniprot_isoform_specific_nucleus_266.accessions.txt").write_text(
        "".join(f"{row['isoform_accession']}\n" for row in direct_rows), encoding="utf-8"
    )
    fasta_from_rows(
        direct_dir / "b3_uniprot_isoform_specific_nucleus_266.fasta",
        [alt_by_id[str(row["isoform_accession"])] for row in direct_rows],
    )

    alt_counts: dict[str, int] = defaultdict(int)
    for row in alt_rows:
        alt_counts[row["canonical_accession"]] += 1
    correction_rows = []
    for parent in sorted(corrected_parents):
        annotation = annotation_by_accession[parent]
        correction_rows.append(
            {
                "canonical_accession": parent,
                "gene_primary": annotation["Gene Names (primary)"],
                "alternative_isoform_count": alt_counts[parent],
                "direct_nuclear_alternative_isoform_count": len(parent_direct_isoforms[parent]),
                "direct_nuclear_alternative_isoforms": ";".join(sorted(parent_direct_isoforms[parent])),
                "canonical_isoform_aware_nucleus": 0,
                "legacy_strict_nucleus_hit": 1,
                "reason": "NUCLEUS_ONLY_ON_NON_DISPLAYED_ALTERNATIVE_ISOFORM",
            }
        )
    write_tsv(
        audit_dir / "b3_legacy_parent_isoform_misattribution_73.tsv",
        [
            "canonical_accession",
            "gene_primary",
            "alternative_isoform_count",
            "direct_nuclear_alternative_isoform_count",
            "direct_nuclear_alternative_isoforms",
            "canonical_isoform_aware_nucleus",
            "legacy_strict_nucleus_hit",
            "reason",
        ],
        correction_rows,
    )

    write_tsv(
        audit_dir / "b3_input_counts.tsv",
        ["metric", "value"],
        [
            {"metric": "alternative_isoform_input", "value": 22131},
            {"metric": "unique_parent_canonical", "value": len({r["canonical_accession"] for r in alt_rows})},
            {"metric": "unique_gene_primary_nonblank", "value": len({r["gene_primary"] for r in alt_rows if r["gene_primary"]})},
            {"metric": "uniprot_isoform_specific_nucleus", "value": 266},
            {"metric": "uniprot_isoform_specific_nucleus_parents", "value": 203},
            {"metric": "corrected_legacy_parent_group", "value": 73},
            {"metric": "all_alternative_isoforms_of_corrected_parents", "value": sum(alt_counts[p] for p in corrected_parents)},
            {"metric": "direct_nuclear_isoforms_of_corrected_parents", "value": 82},
            {"metric": "frozen_annotation_accessions_absent_from_20260902_api", "value": len(frozen_only_accessions)},
            {"metric": "api_accessions_absent_from_frozen_annotations", "value": len(api_only_accessions)},
        ],
    )

    write_tsv(
        audit_dir / "b3_alternative_products_accession_drift.tsv",
        ["accession", "status"],
        [
            *[
                {"accession": accession, "status": "FROZEN_ANNOTATION_ONLY"}
                for accession in frozen_only_accessions
            ],
            *[
                {"accession": accession, "status": "20260902_API_ONLY"}
                for accession in api_only_accessions
            ],
        ],
    )

    manifest_rows = []
    for role, path in [
        ("source_alternative_metadata", args.alternative_metadata),
        ("source_alternative_fasta", args.alternative_fasta),
        ("source_alternative_mapping", args.alternative_mapping),
        ("source_frozen_annotations", args.frozen_annotations),
        ("source_alternative_products", args.alternative_products),
        ("source_legacy_strict_nucleus_accessions", args.legacy_strict_nucleus_accessions),
        ("output_all_metadata", all_tsv),
        ("output_all_fasta", all_fasta),
        ("output_all_mapping", all_mapping),
        ("output_uniprot_direct", direct_dir / "b3_uniprot_isoform_specific_nucleus_266.tsv"),
        ("output_corrected_parent_audit", audit_dir / "b3_legacy_parent_isoform_misattribution_73.tsv"),
    ]:
        manifest_rows.append(
            {"role": role, "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
        )
    write_tsv(audit_dir / "b3_input_manifest.tsv", ["role", "path", "bytes", "sha256"], manifest_rows)

    print("alternative_isoform_input\t22131")
    print("uniprot_isoform_specific_nucleus\t266")
    print("uniprot_isoform_specific_nucleus_parents\t203")
    print("corrected_legacy_parent_group\t73")
    print("direct_nuclear_isoforms_of_corrected_parents\t82")


if __name__ == "__main__":
    main()
