#!/usr/bin/env python3
"""Freeze coordinate-resolved Pfam/InterPro architectures for needed proteins."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from stage25_v3_common import atomic_write_tsv, load_coarse_annotations, read_tsv, sha256


API_ROOT = "https://www.ebi.ac.uk/interpro/api/entry"


def clean(value: str) -> str:
    return (value or "").replace("|", "/").replace("\t", " ").replace("\n", " ").strip()


def validate_json(path: Path) -> dict:
    obj = json.loads(path.read_text())
    if not isinstance(obj, dict) or "results" not in obj:
        raise ValueError(f"invalid InterPro response: {path}")
    return obj


def fetch_one(acc: str, source: str, cache_dir: Path) -> tuple[str, str, Path, str]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / f"{acc}.{source}.json"
    if target.exists():
        validate_json(target)
        return acc, source, target, "reused_valid_cache"
    url = f"{API_ROOT}/{source}/protein/uniprot/{acc}/?page_size=200"
    last_error = ""
    for attempt in range(1, 5):
        partial = cache_dir / f"{acc}.{source}.json.partial.attempt{attempt}"
        proc = subprocess.run(
            [
                "curl", "-L", "--fail", "--silent", "--show-error",
                "--connect-timeout", "20", "--max-time", "120",
                "--retry", "2", "--retry-delay", "2",
                "-A", "Stage25HomologyReview/2.0", "-o", str(partial), url,
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode == 0:
            try:
                validate_json(partial)
                partial.replace(target)
                return acc, source, target, "downloaded"
            except Exception as exc:
                last_error = repr(exc)
        else:
            last_error = proc.stderr.strip()
        time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"InterPro fetch failed for {acc}/{source}: {last_error}")


def parse_architecture(acc: str, source: str, path: Path) -> tuple[dict, list[dict]]:
    obj = validate_json(path)
    features = []
    protein_length = 0
    for result in obj.get("results", []):
        metadata = result.get("metadata", {})
        feature_id = metadata.get("accession", "")
        feature_name = clean(metadata.get("name", ""))
        feature_type = metadata.get("type", "")
        if source == "interpro" and feature_type not in {"domain", "repeat", "family", "homologous_superfamily"}:
            continue
        for protein in result.get("proteins", []) or []:
            if str(protein.get("accession", "")).upper() != acc.upper():
                continue
            protein_length = int(protein.get("protein_length") or protein_length or 0)
            for location_index, location in enumerate(protein.get("entry_protein_locations", []) or [], start=1):
                fragments = location.get("fragments", []) or []
                if not fragments:
                    continue
                start = min(int(x["start"]) for x in fragments)
                end = max(int(x["end"]) for x in fragments)
                features.append(
                    {
                        "accession": acc,
                        "architecture_source": source.upper(),
                        "feature_id": feature_id,
                        "feature_name": feature_name,
                        "feature_type": feature_type,
                        "start": start,
                        "end": end,
                        "location_index": location_index,
                        "fragment_count": len(fragments),
                        "response_sha256": sha256(path),
                    }
                )
    features.sort(key=lambda x: (x["start"], x["end"], x["feature_id"], x["location_index"]))
    status = "available" if features and protein_length else "no_coordinate_resolved_features"
    row = {
        "accession": acc,
        "architecture_source": source.upper(),
        "architecture_status": status,
        "protein_length": protein_length,
        "feature_count": len(features),
        "ordered_feature_ids": "|".join(x["feature_id"] for x in features),
        "ordered_feature_names": "|".join(x["feature_name"] for x in features),
        "ordered_intervals": "|".join(f"{x['start']}-{x['end']}" for x in features),
        "response_path": str(path),
        "response_sha256": sha256(path),
    }
    return row, features


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accessions", required=True, type=Path)
    parser.add_argument("--coarse-annotations", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    out = args.output_dir.resolve()
    cache = out / "resources" / "domain_api_cache"
    coarse = load_coarse_annotations(args.coarse_annotations.resolve())
    accessions = [r["accession"] for r in read_tsv(args.accessions.resolve())]
    jobs = []
    for acc in accessions:
        rec = coarse.get(acc, {})
        source = "pfam" if rec.get("pfam_ids") else "interpro"
        jobs.append((acc, source))

    fetch_rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(fetch_one, acc, source, cache): (acc, source) for acc, source in jobs}
        for future in as_completed(futures):
            acc, source, path, status = future.result()
            fetch_rows.append(
                {
                    "accession": acc,
                    "requested_source": source.upper(),
                    "fetch_status": status,
                    "response_path": str(path),
                    "response_bytes": path.stat().st_size,
                    "response_sha256": sha256(path),
                }
            )
            if len(fetch_rows) % 50 == 0:
                print(f"domain responses validated: {len(fetch_rows)}/{len(jobs)}", flush=True)

    architecture_rows = []
    feature_rows = []
    for acc, source in sorted(jobs):
        row, features = parse_architecture(acc, source, cache / f"{acc}.{source}.json")
        architecture_rows.append(row)
        feature_rows.extend(features)
    if len(architecture_rows) != len(accessions):
        raise RuntimeError("domain architecture row count mismatch")
    atomic_write_tsv(out / "resources" / "domain_fetch_manifest.tsv", sorted(fetch_rows, key=lambda r: r["accession"]))
    atomic_write_tsv(out / "resources" / "domain_architecture_v3.tsv", architecture_rows)
    atomic_write_tsv(
        out / "resources" / "domain_features_long_v3.tsv",
        feature_rows,
        [
            "accession", "architecture_source", "feature_id", "feature_name", "feature_type",
            "start", "end", "location_index", "fragment_count", "response_sha256",
        ],
    )
    print(
        json.dumps(
            {
                "requested_accessions": len(accessions),
                "architecture_available": sum(r["architecture_status"] == "available" for r in architecture_rows),
                "coordinate_features": len(feature_rows),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
