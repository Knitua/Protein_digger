#!/usr/bin/env python3
"""Build the isoform-aware canonical A/B partition and Stage1-A filtered views."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


TAB = "\t"
NUCLEUS_RE = re.compile(r"\bNucleus\b")
ISOFORM_LOCATION_RE = re.compile(r"^\[Isoform ([^\]]+)\]:")
RUN_ID = "stage0_isoform_aware_canonical_partition_20260903"
ANCHOR_VERSION = "working_regulatory_anchor_2368"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-annotations", type=Path, required=True)
    parser.add_argument("--alternative-products", type=Path, required=True)
    parser.add_argument("--anchor-table", type=Path, required=True)
    parser.add_argument("--classification-master", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=TAB)
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {path}")
        return list(reader.fieldnames), list(reader)


def atomic_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite: {path}")
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"stale partial exists: {partial}")
    with partial.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def tsv_payload(fields: list[str], rows: list[dict[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter=TAB,
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    atomic_bytes(path, tsv_payload(fields, rows))


def write_fasta(path: Path, rows: list[dict[str, object]]) -> None:
    payload = "".join(
        f">sp|{row['Entry']}|{row.get('Entry Name', '')} {row.get('Protein names', '')}\n"
        f"{row['Sequence']}\n"
        for row in rows
    ).encode("utf-8")
    atomic_bytes(path, payload)


def write_accessions(path: Path, rows: list[dict[str, object]]) -> None:
    atomic_bytes(path, "".join(f"{row['Entry']}\n" for row in rows).encode("utf-8"))


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
            current = {"name": clean_label(token[5:]), "synonyms": [], "isoform_ids": [], "displayed": False}
        elif current is not None and token.startswith("Synonyms="):
            current["synonyms"] = [clean_label(item) for item in token[9:].split(",")]
        elif current is not None and token.startswith("IsoId="):
            current["isoform_ids"] = [item.strip() for item in token[6:].split(",") if item.strip()]
        elif current is not None and token.startswith("Sequence="):
            current["displayed"] = token[9:].strip() == "Displayed"
    if current is not None:
        records.append(current)
    return records


def location_sections(value: str) -> list[str]:
    value = re.sub(r"^SUBCELLULAR LOCATION:\s*", "", value or "")
    return re.split(r";\s*SUBCELLULAR LOCATION:\s*", value) if value else []


def formal_location_part(section: str) -> str:
    return section.split("Note=", 1)[0].strip()


def atomic_genes(value: str) -> set[str]:
    return {item.strip().upper() for item in re.split(r"[;|]", value or "") if item.strip()}


def classify_canonical_location(
    annotation: dict[str, str], product_value: str
) -> tuple[bool, list[str], list[str], list[str]]:
    products = parse_alternative_products(product_value)
    label_lookup: dict[str, list[dict[str, object]]] = defaultdict(list)
    displayed_labels: set[str] = set()
    for product in products:
        labels = [str(product["name"]), *[str(x) for x in product["synonyms"]]]
        for label in labels:
            label_lookup[normalized_label(label)].append(product)
            if bool(product["displayed"]):
                displayed_labels.add(normalized_label(label))

    canonical_sections: list[str] = []
    alternative_sections: list[str] = []
    alternative_ids: set[str] = set()
    for section in location_sections(annotation.get("Subcellular location [CC]", "")):
        formal = formal_location_part(section)
        if not NUCLEUS_RE.search(formal):
            continue
        tagged = ISOFORM_LOCATION_RE.match(section)
        if tagged is None:
            canonical_sections.append(formal)
            continue
        label = tagged.group(1)
        norm_label = normalized_label(label)
        products_for_label = label_lookup.get(norm_label, [])
        if not products_for_label:
            raise RuntimeError(f"cannot map nuclear isoform label: {annotation['Entry']} / {label}")
        if norm_label in displayed_labels:
            canonical_sections.append(formal)
        for product in products_for_label:
            if bool(product["displayed"]):
                continue
            alternative_sections.append(formal)
            alternative_ids.update(str(item) for item in product["isoform_ids"])
    return bool(canonical_sections), canonical_sections, alternative_sections, sorted(alternative_ids)


def append_partition_fields(
    row: dict[str, str],
    partition: str,
    canonical_sections: list[str],
    alternative_sections: list[str],
    alternative_ids: list[str],
) -> dict[str, object]:
    return {
        **row,
        "stage0_localization_partition": partition,
        "canonical_nucleus_annotation": 1 if partition == "A_CANONICAL_NUCLEUS" else 0,
        "canonical_nucleus_evidence_sections": " || ".join(canonical_sections),
        "non_displayed_isoform_nucleus_sections": " || ".join(alternative_sections),
        "non_displayed_nuclear_isoform_accessions": ";".join(alternative_ids),
        "partition_rule_version": RUN_ID,
    }


def apply_anchor_filter(
    rows: list[dict[str, object]], anchor_accessions: set[str], anchor_genes: set[str]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    retained: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    for row in rows:
        accession = str(row["Entry"])
        genes = atomic_genes(str(row.get("Gene Names (primary)", "")))
        accession_hit = accession in anchor_accessions
        matched_genes = sorted(genes & anchor_genes)
        if accession_hit or matched_genes:
            reasons = []
            if accession_hit:
                reasons.append("ANCHOR_ACCESSION")
            if matched_genes:
                reasons.append("ANCHOR_ATOMIC_GENE")
            excluded.append(
                {
                    **row,
                    "anchor_accession_hit": int(accession_hit),
                    "matched_anchor_atomic_genes": ";".join(matched_genes),
                    "anchor_exclusion_reason": ";".join(reasons),
                    "anchor_version": ANCHOR_VERSION,
                }
            )
        else:
            retained.append(row)
    return retained, excluded


def build_a_classification(
    output_root: Path,
    current_master: list[dict[str, str]],
    current_fields: list[str],
    route_a_non_anchor: list[dict[str, object]],
) -> list[Path]:
    allowed = {str(row["Entry"]) for row in route_a_non_anchor}
    master = [dict(row) for row in current_master if row["candidate_accession"] in allowed]
    if len(master) != 3418 or {row["candidate_accession"] for row in master} != allowed:
        raise RuntimeError("Stage1-A evidence master does not exactly cover new Route A non-anchor set")
    for row in master:
        row["source_universe"] = "isoform_aware_canonical_nucleus_5669"
        row["pass_rule"] = (
            "UniProt Nucleus annotation applies to canonical/displayed sequence; "
            "non-anchor by accession and atomic gene"
        )
        row["run_id"] = RUN_ID
        if "working_stage1_run_id" in row:
            row["working_stage1_run_id"] = RUN_ID
        if "classification_scope" in row:
            row["classification_scope"] = "ROUTE_A_ISOFORM_AWARE_NON_ANCHOR_3418"
        if "classification_method_version" in row:
            row["classification_method_version"] = "A1_FUNCTIONAL_TIERING_ISOFORM_AWARE_20260903"
    master.sort(key=lambda row: row["candidate_accession"])

    stage1 = output_root / "stage1_A"
    outputs: list[Path] = []
    base_fields = current_fields[: current_fields.index("regulatory_evidence_status")]
    nonanchor_path = stage1 / "02_non_anchor/route_A_uniprot_canonical_nucleus_non_anchor_3418.tsv"
    write_tsv(nonanchor_path, base_fields, master)
    outputs.append(nonanchor_path)
    master_path = stage1 / "route_A_candidate_classification_3418.tsv"
    write_tsv(master_path, current_fields, master)
    outputs.append(master_path)

    l1_names = {
        "KNOWN_CURATED": "01_known_curated_541.tsv",
        "KNOWN_FUNCTIONAL_DIRECT": "02_known_functional_direct_236.tsv",
        "KNOWN_COMPLEX_ASSOCIATED": "03_known_complex_associated_440.tsv",
        "NO_KNOWN_REGULATORY_EVIDENCE": "04_no_known_regulatory_evidence_2201.tsv",
    }
    l1_expected = {"KNOWN_CURATED": 541, "KNOWN_FUNCTIONAL_DIRECT": 236, "KNOWN_COMPLEX_ASSOCIATED": 440, "NO_KNOWN_REGULATORY_EVIDENCE": 2201}
    for status, filename in l1_names.items():
        subset = [row for row in master if row["regulatory_evidence_status"] == status]
        if len(subset) != l1_expected[status]:
            raise RuntimeError(f"unexpected {status} count: {len(subset)}")
        path = stage1 / "03_level1_regulatory_evidence" / filename
        write_tsv(path, current_fields, subset)
        outputs.append(path)

    no_known = [row for row in master if row["regulatory_evidence_status"] == "NO_KNOWN_REGULATORY_EVIDENCE"]
    l2_names = {
        "HK_STRONG": "01_hk_strong_0.tsv",
        "HK_PROBABLE": "02_hk_probable_170.tsv",
        "UNDER_CHARACTERIZED": "03_under_characterized_924.tsv",
        "MIXED_MOONLIGHTING": "04_mixed_moonlighting_1093.tsv",
        "EVIDENCE_INCOMPLETE": "05_evidence_incomplete_14.tsv",
    }
    l2_expected = {"HK_STRONG": 0, "HK_PROBABLE": 170, "UNDER_CHARACTERIZED": 924, "MIXED_MOONLIGHTING": 1093, "EVIDENCE_INCOMPLETE": 14}
    for status, filename in l2_names.items():
        subset = [row for row in no_known if row["level2_no_known_functional_tier"] == status]
        if len(subset) != l2_expected[status]:
            raise RuntimeError(f"unexpected {status} count: {len(subset)}")
        path = stage1 / "04_level2_no_known_regulatory" / filename
        write_tsv(path, current_fields, subset)
        outputs.append(path)

    summary = [
        *({"classification_level": "LEVEL1", "class": key, "proteins": value} for key, value in l1_expected.items()),
        *({"classification_level": "LEVEL2_NO_KNOWN", "class": key, "proteins": value} for key, value in l2_expected.items()),
    ]
    path = stage1 / "05_audit/route_A_classification_count_summary.tsv"
    write_tsv(path, ["classification_level", "class", "proteins"], summary)
    outputs.append(path)
    return outputs


def main() -> None:
    args = parse_args()
    if args.output_root.exists():
        raise FileExistsError(f"output root exists: {args.output_root}")
    args.output_root.mkdir(parents=True)

    annotation_fields, annotations = read_tsv(args.canonical_annotations)
    product_fields, products = read_tsv(args.alternative_products)
    anchor_fields, anchors = read_tsv(args.anchor_table)
    class_fields, class_master = read_tsv(args.classification_master)
    if len(annotations) != 20416 or len({row["Entry"] for row in annotations}) != 20416:
        raise RuntimeError("canonical annotation universe must contain 20,416 unique accessions")
    if len(products) != 20416 or len({row["Entry"] for row in products}) != 20416:
        raise RuntimeError("alternative-products resource must cover 20,416 unique accessions")
    if len(anchors) != 2368 or len({row["anchor_accession"] for row in anchors}) != 2368:
        raise RuntimeError("anchor table must contain 2,368 unique accessions")
    if len(class_master) != 3487:
        raise RuntimeError("frozen evidence master must contain 3,487 rows before set projection")

    product_by = {row["Entry"]: row for row in products}
    route_a: list[dict[str, object]] = []
    route_b: list[dict[str, object]] = []
    isoform_only: list[dict[str, object]] = []
    for annotation in annotations:
        accession = annotation["Entry"]
        product_row = product_by.get(accession, {})
        product_value = product_row.get("Alternative products (isoforms)", "")
        canonical, canonical_sections, alternative_sections, alternative_ids = classify_canonical_location(
            annotation, product_value
        )
        partition = "A_CANONICAL_NUCLEUS" if canonical else "B_NO_CANONICAL_NUCLEUS"
        extended = append_partition_fields(
            annotation, partition, canonical_sections, alternative_sections, alternative_ids
        )
        (route_a if canonical else route_b).append(extended)
        if not canonical and alternative_sections:
            isoform_only.append(
                {
                    "canonical_accession": accession,
                    "gene_primary": annotation.get("Gene Names (primary)", ""),
                    "protein_name": annotation.get("Protein names", ""),
                    "alternative_isoform_count": sum(
                        len(product["isoform_ids"])
                        for product in parse_alternative_products(product_value)
                        if not bool(product["displayed"])
                    ),
                    "nuclear_alternative_isoform_count": len(alternative_ids),
                    "nuclear_alternative_isoforms": ";".join(alternative_ids),
                    "nuclear_location_sections": " || ".join(alternative_sections),
                    "stage0_partition": partition,
                    "partition_reason": "NUCLEUS_RESTRICTED_TO_NON_DISPLAYED_ALTERNATIVE_ISOFORM",
                }
            )
    route_a.sort(key=lambda row: str(row["Entry"]))
    route_b.sort(key=lambda row: str(row["Entry"]))
    isoform_only.sort(key=lambda row: str(row["canonical_accession"]))
    if len(route_a) != 5669 or len(route_b) != 14747 or len(isoform_only) != 73:
        raise RuntimeError(
            f"unexpected partition counts A={len(route_a)} B={len(route_b)} isoform_only={len(isoform_only)}"
        )
    if {str(row["Entry"]) for row in route_a} & {str(row["Entry"]) for row in route_b}:
        raise RuntimeError("A and B partition overlap")

    anchor_accessions = {row["anchor_accession"] for row in anchors}
    anchor_genes = {gene for row in anchors for gene in atomic_genes(row["anchor_gene"])}
    if len(anchor_genes) != 2372:
        raise RuntimeError(f"unexpected atomic anchor gene count: {len(anchor_genes)}")
    a_nonanchor, a_excluded = apply_anchor_filter(route_a, anchor_accessions, anchor_genes)
    b_nonanchor, b_excluded = apply_anchor_filter(route_b, anchor_accessions, anchor_genes)
    if (len(a_nonanchor), len(a_excluded), len(b_nonanchor), len(b_excluded)) != (3418, 2251, 14625, 122):
        raise RuntimeError(
            "unexpected anchor partition counts: "
            f"A={len(a_nonanchor)}+{len(a_excluded)} B={len(b_nonanchor)}+{len(b_excluded)}"
        )

    partition_fields = annotation_fields + [
        "stage0_localization_partition",
        "canonical_nucleus_annotation",
        "canonical_nucleus_evidence_sections",
        "non_displayed_isoform_nucleus_sections",
        "non_displayed_nuclear_isoform_accessions",
        "partition_rule_version",
    ]
    exclusion_fields = partition_fields + [
        "anchor_accession_hit",
        "matched_anchor_atomic_genes",
        "anchor_exclusion_reason",
        "anchor_version",
    ]
    pdir = args.output_root / "stage0/01_canonical_partition"
    edir = args.output_root / "stage0/02_anchor_exclusion"
    adir = args.output_root / "stage0/03_audit"
    outputs: list[Path] = []

    datasets = [
        (pdir / "route_A_canonical_nucleus_5669.tsv", route_a, partition_fields),
        (pdir / "route_B_without_canonical_nucleus_14747.tsv", route_b, partition_fields),
        (edir / "route_A_anchor_excluded_2251.tsv", a_excluded, exclusion_fields),
        (edir / "route_B_anchor_excluded_122.tsv", b_excluded, exclusion_fields),
        (edir / "route_A_non_anchor_3418.tsv", a_nonanchor, partition_fields),
        (edir / "route_B_non_anchor_14625.tsv", b_nonanchor, partition_fields),
    ]
    for path, rows, fields in datasets:
        write_tsv(path, fields, rows)
        outputs.append(path)
        if "anchor_excluded" not in path.name:
            fasta = path.with_suffix(".fasta")
            accessions = path.with_suffix(".accessions.txt")
            write_fasta(fasta, rows)
            write_accessions(accessions, rows)
            outputs.extend([fasta, accessions])

    isoform_audit_path = adir / "canonical_nucleus_restricted_to_non_displayed_isoforms_73.tsv"
    isoform_audit_fields = [
        "canonical_accession",
        "gene_primary",
        "protein_name",
        "alternative_isoform_count",
        "nuclear_alternative_isoform_count",
        "nuclear_alternative_isoforms",
        "nuclear_location_sections",
        "stage0_partition",
        "partition_reason",
    ]
    write_tsv(isoform_audit_path, isoform_audit_fields, isoform_only)
    outputs.append(isoform_audit_path)

    isoform_only_ids = {str(row["canonical_accession"]) for row in isoform_only}
    b_nonanchor_69 = [row for row in b_nonanchor if str(row["Entry"]) in isoform_only_ids]
    if len(b_nonanchor_69) != 69:
        raise RuntimeError(f"expected 69 non-anchor isoform-only canonical records, found {len(b_nonanchor_69)}")
    b69_path = adir / "route_B_isoform_only_nucleus_non_anchor_69.tsv"
    write_tsv(b69_path, partition_fields, b_nonanchor_69)
    outputs.append(b69_path)

    counts = [
        {"stage": "REFERENCE", "dataset": "reviewed_human_canonical", "proteins": 20416, "definition": "UniProt reviewed human canonical universe"},
        {"stage": "PARTITION", "dataset": "route_A_canonical_nucleus", "proteins": 5669, "definition": "Nucleus annotation applies to canonical/displayed sequence"},
        {"stage": "PARTITION", "dataset": "route_B_without_canonical_nucleus", "proteins": 14747, "definition": "exact canonical complement of Route A"},
        {"stage": "ANCHOR_EXCLUSION", "dataset": "route_A_anchor_excluded", "proteins": 2251, "definition": "anchor accession or atomic gene"},
        {"stage": "ANCHOR_EXCLUSION", "dataset": "route_B_anchor_excluded", "proteins": 122, "definition": "anchor accession or atomic gene"},
        {"stage": "STAGE1_INPUT", "dataset": "route_A_non_anchor", "proteins": 3418, "definition": "Stage1-A canonical input"},
        {"stage": "STAGE1_INPUT", "dataset": "route_B_non_anchor", "proteins": 14625, "definition": "Stage1-B1/B2 common canonical input"},
        {"stage": "AUDIT", "dataset": "non_displayed_isoform_only_nucleus_parents", "proteins": 73, "definition": "canonical belongs to B; alternative isoforms handled independently by B3"},
        {"stage": "AUDIT", "dataset": "non_displayed_isoform_only_nucleus_non_anchor", "proteins": 69, "definition": "members of the 73 remaining after anchor exclusion"},
    ]
    counts_path = adir / "stage0_partition_counts.tsv"
    write_tsv(counts_path, ["stage", "dataset", "proteins", "definition"], counts)
    outputs.append(counts_path)

    # Stage1-A evidence rows are a deterministic projection onto the new A set.
    a_outputs = build_a_classification(args.output_root, class_master, class_fields, a_nonanchor)
    outputs.extend(a_outputs)

    # Stage1-A and Stage1-B interfaces copied from the validated Stage0 datasets.
    interface_files = [
        (pdir / "route_A_canonical_nucleus_5669.tsv", args.output_root / "stage1_A/01_nuclear_annotation/route_A_uniprot_canonical_nucleus_5669.tsv"),
        (pdir / "route_A_canonical_nucleus_5669.fasta", args.output_root / "stage1_A/01_nuclear_annotation/route_A_uniprot_canonical_nucleus_5669.fasta"),
        (edir / "route_A_non_anchor_3418.tsv", args.output_root / "stage1_A/02_non_anchor/route_A_uniprot_canonical_nucleus_non_anchor_3418_metadata.tsv"),
        (pdir / "route_B_without_canonical_nucleus_14747.tsv", args.output_root / "stage1_B_common/01_without_nuclear_annotation/route_B_uniprot_without_canonical_nucleus_14747.tsv"),
        (pdir / "route_B_without_canonical_nucleus_14747.fasta", args.output_root / "stage1_B_common/01_without_nuclear_annotation/route_B_uniprot_without_canonical_nucleus_14747.fasta"),
        (edir / "route_B_non_anchor_14625.tsv", args.output_root / "stage1_B_common/02_non_anchor/route_B_uniprot_without_canonical_nucleus_non_anchor_14625.tsv"),
        (edir / "route_B_non_anchor_14625.fasta", args.output_root / "stage1_B_common/02_non_anchor/route_B_uniprot_without_canonical_nucleus_non_anchor_14625.fasta"),
        (edir / "route_B_anchor_excluded_122.tsv", args.output_root / "stage1_B_common/03_audit/route_B_working_anchor_excluded_122.tsv"),
        (b69_path, args.output_root / "stage1_B_common/03_audit/route_B_isoform_only_nucleus_non_anchor_69.tsv"),
    ]
    for source, destination in interface_files:
        atomic_bytes(destination, source.read_bytes())
        outputs.append(destination)

    b_common_counts = [
        {"metric": "reviewed_canonical_reference", "value": 20416, "note": "frozen UniProt human reviewed canonical universe"},
        {"metric": "route_A_canonical_nucleus", "value": 5669, "note": "Nucleus annotation applies to canonical/displayed sequence"},
        {"metric": "route_B_without_canonical_nucleus", "value": 14747, "note": "exact canonical complement of Route A"},
        {"metric": "working_anchor_accessions", "value": 2368, "note": "current regulatory anchor interface"},
        {"metric": "working_anchor_atomic_genes", "value": 2372, "note": "semicolon/pipe split, uppercased"},
        {"metric": "route_B_anchor_excluded", "value": 122, "note": "accession or atomic-gene match"},
        {"metric": "route_B_common_non_anchor", "value": 14625, "note": "shared canonical input for B1 and B2"},
        {"metric": "route_B_isoform_only_nucleus_non_anchor", "value": 69, "note": "canonical records assigned to B; alternative isoforms handled by B3"},
    ]
    b_counts_path = args.output_root / "stage1_B_common/03_audit/route_B_common_input_counts.tsv"
    write_tsv(b_counts_path, ["metric", "value", "note"], b_common_counts)
    outputs.append(b_counts_path)

    # Independent cross-checks.
    if 5669 + 14747 != 20416 or 3418 + 14625 != 18043 or 2251 + 122 != 2373:
        raise RuntimeError("partition arithmetic failed")
    if {str(row["Entry"]) for row in a_nonanchor} & {str(row["Entry"]) for row in b_nonanchor}:
        raise RuntimeError("non-anchor A/B inputs overlap")

    manifest_rows = []
    for role, path, records in [
        ("INPUT", args.canonical_annotations, 20416),
        ("INPUT", args.alternative_products, 20416),
        ("INPUT", args.anchor_table, 2368),
        ("INPUT", args.classification_master, 3487),
        *[("OUTPUT", path, max(0, sum(1 for _ in path.open(encoding="utf-8", errors="ignore")) - 1) if path.suffix == ".tsv" else "") for path in outputs],
    ]:
        manifest_rows.append(
            {
                "role": role,
                "path": str(path),
                "records": records,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "run_id": RUN_ID,
            }
        )
    manifest_path = adir / "stage0_manifest.tsv"
    write_tsv(manifest_path, ["role", "path", "records", "bytes", "sha256", "run_id"], manifest_rows)

    l1 = Counter(row["regulatory_evidence_status"] for row in class_master if row["candidate_accession"] in {str(r["Entry"]) for r in a_nonanchor})
    print("STAGE0_BUILD_PASS")
    print(f"canonical_partition A={len(route_a)} B={len(route_b)}")
    print(f"anchor_exclusion A={len(a_excluded)} B={len(b_excluded)}")
    print(f"stage1_inputs A={len(a_nonanchor)} B={len(b_nonanchor)}")
    print(f"isoform_only_nucleus parents={len(isoform_only)} nonanchor={len(b_nonanchor_69)}")
    print(f"route_A_level1={dict(l1)}")


if __name__ == "__main__":
    main()
