#!/usr/bin/env python3
"""Write immutable hashes and execution-order audit after the independent case audit."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from stage25_v3_common import atomic_write_tsv, json_dump, sha256


def stamp(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()

    core = out / "audit" / "STAGE2_5_V3_CORE_FREEZE.json"
    sensitivity = out / "validation" / "positive_control_threshold_sensitivity_summary.json"
    case = out / "postfreeze_case_audit" / "RAS_postfreeze_audit_v3.tsv"
    if not all(path.exists() for path in (core, sensitivity, case)):
        raise RuntimeError("core freeze, sensitivity output, and case audit must all exist")

    ordered = core.stat().st_mtime_ns < sensitivity.stat().st_mtime_ns < case.stat().st_mtime_ns
    chronology = {
        "run_id": "stage2_5_homology_function_review_v3_20261006",
        "core_freeze_utc": stamp(core),
        "positive_control_sensitivity_utc": stamp(sensitivity),
        "case_audit_utc": stamp(case),
        "execution_order_verified": ordered,
        "order": ["generic_core_freeze", "non_case_positive_control_sensitivity", "case_audit"],
        "core_freeze_sha256": sha256(core),
        "sensitivity_summary_sha256": sha256(sensitivity),
        "case_audit_sha256": sha256(case),
    }
    if not ordered:
        raise RuntimeError("execution-order timestamps are inconsistent")
    json_dump(out / "audit" / "stage2_5_v3_execution_order_audit.json", chronology)

    rows = []
    folders = ("validation", "postfreeze_case_audit")
    for folder in folders:
        for path in sorted((out / folder).rglob("*")):
            if path.is_file() and ".partial." not in path.name:
                rows.append(
                    {
                        "phase": folder,
                        "path": str(path.resolve()),
                        "bytes": path.stat().st_size,
                        "sha256": sha256(path),
                    }
                )
    audit_path = out / "audit" / "stage2_5_v3_postfreeze_manifest.tsv"
    atomic_write_tsv(audit_path, rows)
    print(json.dumps(chronology, indent=2))


if __name__ == "__main__":
    main()
