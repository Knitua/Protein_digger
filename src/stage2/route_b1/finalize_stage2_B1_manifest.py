#!/usr/bin/env python3
"""Create a relocatable manifest for the published Stage2-B1 directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import os
from datetime import datetime, timezone
from pathlib import Path


RUN_ID = "stage2_B1_current1350_hi_union_rf2_official_precision_20260904"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def row_count(path: Path) -> int:
    if path.suffix == ".tsv":
        with path.open(errors="replace") as handle:
            return max(0, sum(1 for _ in handle) - 1)
    return 1


def render(rows: list[dict[str, object]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=[
            "role", "path", "records", "size_bytes", "sha256",
            "description", "run_id", "generated_at_utc",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def safe_write(path: Path, content: bytes) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    partial = path.with_name(path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"Partial exists: {partial}")
    with partial.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage2-b1", type=Path, required=True)
    parser.add_argument("--stage1-b1", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--hi-edges", type=Path, required=True)
    parser.add_argument("--rf2-final80", type=Path, required=True)
    parser.add_argument("--rf2-final90", type=Path, required=True)
    args = parser.parse_args()

    audit_dir = args.stage2_b1 / "04_audit"
    manifest = audit_dir / "stage2_B1_ppi_screen_manifest.tsv"
    archived = audit_dir / "stage2_B1_ppi_screen_manifest.prepublish_paths.tsv"
    if not manifest.exists() or archived.exists():
        raise FileExistsError("Unexpected manifest publication state")
    manifest.rename(archived)

    timestamp = datetime.now(timezone.utc).isoformat()
    external = [
        ("INPUT", args.stage1_b1, "formal corrected Stage1-B1 candidates"),
        ("INPUT", args.anchors, "current 2,368 PPI anchors"),
        ("INPUT", args.hi_edges, "current-anchor HI-union edges"),
        ("INPUT", args.rf2_final80, "current-anchor RF2 final80 edges"),
        ("INPUT", args.rf2_final90, "current-anchor RF2 final90 edges"),
    ]
    rows: list[dict[str, object]] = []
    for role, path, description in external:
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append(
            {
                "role": role,
                "path": str(path),
                "records": row_count(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "description": description,
                "run_id": RUN_ID,
                "generated_at_utc": timestamp,
            }
        )

    excluded = {manifest.resolve(), archived.resolve()}
    for path in sorted(args.stage2_b1.rglob("*")):
        if not path.is_file() or path.resolve() in excluded or path.name.endswith(".partial"):
            continue
        relative = path.relative_to(args.stage2_b1)
        role = "CODE" if relative.parts[0] == "scripts" else "OUTPUT"
        if relative.parts[0] == "04_audit":
            role = "AUDIT"
        elif relative.name == "README.md":
            role = "DOCUMENTATION"
        rows.append(
            {
                "role": role,
                "path": str(relative),
                "records": row_count(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "description": "published Stage2-B1 artifact",
                "run_id": RUN_ID,
                "generated_at_utc": timestamp,
            }
        )
    rows.sort(key=lambda row: (str(row["role"]), str(row["path"])))
    safe_write(manifest, render(rows))
    print(f"manifest_records\t{len(rows)}")


if __name__ == "__main__":
    main()
