#!/usr/bin/env python3
"""Build evidence-coded GOA and direct-role Reactome resources for Stage2.5 v3."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from stage25_v3_common import atomic_write_tsv, json_dump, read_tsv, sha256


GO_PRIMARY_CODES = {"EXP", "IDA", "IMP", "IGI"}
GO_DISALLOWED_QUALIFIERS = {"NOT", "contributes_to", "colocalizes_with"}
REACTOME_API = "https://reactome.org/ContentService/data/query/ids"


def build_goa(goa_path: Path, wanted: set[str], out: Path) -> dict:
    rows = []
    included = 0
    with gzip.open(goa_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("!"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                continue
            acc = cols[1].split("-")[0]
            if acc not in wanted or cols[8] not in {"P", "F"}:
                continue
            qualifiers = {x for x in cols[3].split("|") if x}
            evidence = cols[6]
            positive = not bool(qualifiers & GO_DISALLOWED_QUALIFIERS)
            allowed = evidence in GO_PRIMARY_CODES and positive
            if allowed:
                reason = "experimental_function_code"
                included += 1
            elif not positive:
                reason = "disallowed_qualifier"
            else:
                reason = "non_primary_evidence_code"
            rows.append(
                {
                    "accession": acc,
                    "gene": cols[2],
                    "go_id": cols[4],
                    "aspect": cols[8],
                    "qualifier": cols[3],
                    "reference": cols[5],
                    "evidence_code": evidence,
                    "with_from": cols[7],
                    "assigned_by": cols[14],
                    "annotation_date": cols[13],
                    "included_in_primary": "yes" if allowed else "no",
                    "inclusion_reason": reason,
                }
            )
    rows.sort(key=lambda r: (r["accession"], r["aspect"], r["go_id"], r["evidence_code"], r["reference"]))
    atomic_write_tsv(out / "resources" / "goa_stage25_evidence_long_v3.tsv", rows)
    atomic_write_tsv(
        out / "resources" / "goa_direct_function_v3.tsv",
        [r for r in rows if r["included_in_primary"] == "yes"],
        list(rows[0]) if rows else [],
    )
    return {"goa_rows_for_stage2_universe": len(rows), "goa_primary_direct_rows": included}


def load_jsonl_by_key(path: Path, key_name: str) -> dict[str, dict]:
    result = {}
    if not path.exists():
        return result
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            obj = json.loads(line)
            key = obj.get(key_name)
            if key is not None:
                result[str(key)] = obj
    return result


def response_key(obj: dict, numeric: bool) -> str:
    return str(obj.get("dbId")) if numeric else str(obj.get("stId"))


def fetch_batches(ids: list[str], numeric: bool, namespace: str, cache_dir: Path) -> dict[str, dict]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    result = {}
    for offset in range(0, len(ids), 20):
        batch = ids[offset : offset + 20]
        payload = ",".join(batch)
        key = hashlib.sha256(payload.encode()).hexdigest()[:20]
        target = cache_dir / f"batch_{key}.json"
        if target.exists():
            envelope = json.loads(target.read_text())
            if envelope.get("requested") != batch:
                raise RuntimeError(f"cached Reactome request mismatch: {target}")
            objects = envelope.get("objects", [])
        else:
            last_error = ""
            for attempt in range(1, 5):
                partial = cache_dir / f"batch_{key}.json.partial.attempt{attempt}"
                proc = subprocess.run(
                    [
                        "curl", "-L", "--fail", "--silent", "--show-error",
                        "--connect-timeout", "20", "--max-time", "120",
                        "--retry", "2", "--retry-delay", "2",
                        "-A", "Stage25HomologyReview/2.0",
                        "-H", "Accept: application/json",
                        "-H", "Content-Type: text/plain",
                        "--data-binary", payload,
                        REACTOME_API,
                    ],
                    capture_output=True,
                    text=True,
                )
                if proc.returncode == 0:
                    try:
                        objects = json.loads(proc.stdout)
                        if not isinstance(objects, list):
                            raise ValueError("Reactome batch response is not a list")
                        envelope = {"requested": batch, "objects": objects}
                        partial.write_text(json.dumps(envelope, separators=(",", ":")) + "\n")
                        partial.replace(target)
                        break
                    except Exception as exc:
                        last_error = repr(exc)
                else:
                    last_error = proc.stderr.strip()
                time.sleep(2 ** (attempt - 1))
            else:
                raise RuntimeError(f"Reactome {namespace} batch failed: {last_error}")
        returned = {response_key(obj, numeric) for obj in objects if isinstance(obj, dict)}
        missing = set(batch) - returned
        if missing:
            raise RuntimeError(f"Reactome omitted {namespace} IDs: {sorted(missing)}")
        for obj in objects:
            result[response_key(obj, numeric)] = obj
        print(f"Reactome {namespace}: {min(offset + 20, len(ids))}/{len(ids)}", flush=True)
    return result


def merge_reused_and_fetch(
    ids: set[str], numeric: bool, namespace: str, reused: dict[str, dict], cache_root: Path
) -> dict[str, dict]:
    found = {key: reused[key] for key in ids if key in reused}
    missing = sorted(ids - found.keys(), key=(lambda x: int(x)) if numeric else None)
    if missing:
        found.update(fetch_batches(missing, numeric, namespace, cache_root / f"batches_{namespace}"))
    if set(found) != ids:
        raise RuntimeError(f"Reactome {namespace} coverage mismatch")
    return found


def direct_accession(entity: dict) -> str:
    if entity.get("schemaClass", entity.get("className")) != "EntityWithAccessionedSequence":
        return ""
    ref = entity.get("referenceEntity") or {}
    acc = ref.get("identifier", "") if isinstance(ref, dict) else ""
    return acc.split("-")[0]


def build_reactome(
    mapping_path: Path,
    target_accessions: set[str],
    out: Path,
    reused_cache_dir: Path,
) -> dict:
    raw_events: dict[str, set[str]] = defaultdict(set)
    raw_rows = []
    with mapping_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6 or cols[0] not in target_accessions or cols[5] != "Homo sapiens" or not cols[1].startswith("R-HSA-"):
                continue
            raw_events[cols[0]].add(cols[1])
            raw_rows.append(
                {
                    "accession": cols[0],
                    "event_id": cols[1],
                    "event_name": cols[3],
                    "mapping_evidence": cols[4],
                    "species": cols[5],
                }
            )
    atomic_write_tsv(out / "resources" / "reactome_raw_mapping_target_accessions_v3.tsv", raw_rows)
    event_ids = {event for values in raw_events.values() for event in values}

    reused_events = load_jsonl_by_key(reused_cache_dir / "A1_relevant_events.jsonl", "stId")
    events = merge_reused_and_fetch(
        event_ids,
        False,
        "events",
        reused_events,
        out / "resources" / "reactome_api_cache",
    )
    role_ids = set()
    for event in events.values():
        for stub in (event.get("catalystActivity", []) or []) + (event.get("regulatedBy", []) or []):
            if isinstance(stub, dict) and stub.get("dbId") is not None:
                role_ids.add(str(stub["dbId"]))
    reused_roles = load_jsonl_by_key(reused_cache_dir / "A1_role_objects.jsonl", "dbId")
    roles = merge_reused_and_fetch(
        role_ids,
        True,
        "roles",
        reused_roles,
        out / "resources" / "reactome_api_cache",
    )

    entity_ids = set()
    for role in roles.values():
        for unit in role.get("activeUnit", []) or []:
            if isinstance(unit, dict) and unit.get("dbId") is not None:
                entity_ids.add(str(unit["dbId"]))
        for key in ("physicalEntity", "regulator"):
            holder = role.get(key)
            if isinstance(holder, dict) and holder.get("dbId") is not None:
                entity_ids.add(str(holder["dbId"]))
    reused_entities = load_jsonl_by_key(reused_cache_dir / "A1_entity_closure.jsonl", "dbId")
    entities = merge_reused_and_fetch(
        entity_ids,
        True,
        "entities",
        reused_entities,
        out / "resources" / "reactome_api_cache",
    )

    evidence = []
    for event_id, event in sorted(events.items()):
        if event.get("isInferred") or event.get("isInDisease") or event.get("speciesName") != "Homo sapiens":
            continue
        event_name = event.get("displayName", "")
        event_direct = []
        for role_name, stubs in (("CATALYST", event.get("catalystActivity", [])), ("REGULATOR", event.get("regulatedBy", []))):
            for stub in stubs or []:
                if not isinstance(stub, dict) or stub.get("dbId") is None:
                    continue
                role = roles[str(stub["dbId"])]
                holders = []
                if role_name == "CATALYST":
                    active = role.get("activeUnit", []) or []
                    if active:
                        holders.extend(active)
                    else:
                        holder = role.get("physicalEntity")
                        if isinstance(holder, dict):
                            holders.append(holder)
                else:
                    holder = role.get("regulator")
                    if isinstance(holder, dict):
                        holders.append(holder)
                for holder in holders:
                    entity = entities.get(str(holder.get("dbId")), holder)
                    acc = direct_accession(entity)
                    if acc and acc in target_accessions:
                        event_direct.append((acc, role_name, str(role.get("dbId", "")), str(entity.get("dbId", ""))))
        for acc, role_name, role_id, entity_id in sorted(set(event_direct)):
            evidence.append(
                {
                    "accession": acc,
                    "event_id": event_id,
                    "event_name": event_name,
                    "event_role": role_name,
                    "participant_directness": "DIRECT_ENTITY",
                    "inferred": "false",
                    "in_disease": "false",
                    "role_object_id": role_id,
                    "physical_entity_id": entity_id,
                }
            )
    atomic_write_tsv(
        out / "resources" / "reactome_direct_function_v3.tsv",
        evidence,
        [
            "accession", "event_id", "event_name", "event_role", "participant_directness",
            "inferred", "in_disease", "role_object_id", "physical_entity_id",
        ],
    )
    return {
        "reactome_target_accessions": len(target_accessions),
        "reactome_raw_mapping_rows": len(raw_rows),
        "reactome_unique_events_queried": len(event_ids),
        "reactome_role_objects": len(role_ids),
        "reactome_entity_objects": len(entity_ids),
        "reactome_direct_role_rows": len(evidence),
        "reactome_direct_role_accessions": len({r["accession"] for r in evidence}),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--goa", required=True, type=Path)
    parser.add_argument("--reactome-map", required=True, type=Path)
    parser.add_argument("--reused-reactome-cache", required=True, type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve()

    wanted = set()
    for name, column in (
        ("stage2_5_v3_final80_positive_seeds_1862.tsv", "seed_accession"),
        ("stage2_5_v3_final80_failed_review_pool_2631.tsv", "candidate_accession"),
    ):
        wanted.update(row[column] for row in read_tsv(out / "inputs" / name))
    target_accessions = {r["accession"] for r in read_tsv(out / "resources" / "domain_accessions_needed.tsv")}

    summary = {}
    summary.update(build_goa(args.goa.resolve(), wanted, out))
    summary.update(
        build_reactome(
            args.reactome_map.resolve(),
            target_accessions,
            out,
            args.reused_reactome_cache.resolve(),
        )
    )
    summary["goa_sha256"] = sha256(args.goa.resolve())
    summary["reactome_map_sha256"] = sha256(args.reactome_map.resolve())
    json_dump(out / "audit" / "stage2_5_v3_function_resource_summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
