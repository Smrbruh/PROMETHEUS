"""Read-only Phase 1.1 evidence audit; writes one new report exclusively.

No mechanics, geometry reconstruction, field projection or matrix assembly is
performed. Labels are recomputed from archived scalar evidence. A failed audit
is written before exiting nonzero, and an existing report is never replaced.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from prometheus.reference_labels import reference_label

ANGLES = (-20., 0., 20.)
REQUIRED_DOCUMENTS = (
    "REFERENCE_LABEL_PROTOCOL.md", "NUMERICAL_TRUTH_AUDIT.md",
    "BINARY_GEOMETRY_AUDIT.md", "KINEMATIC_AUDIT.md",
    "PHASE1_1_DECISION.md", "PHASE1_1_REPRODUCE.md",
)
RESIDUAL_FAILURE = "P2 free-DOF equilibrium residual exceeds 1e-8"


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def binary_field_hash(path):
    """Reproduce the evaluator's field identity, not a repair provenance ID.

    The repaired catalogue candidate_hash identifies the construction rule and
    source history. Mechanics records instead hash the frozen little-endian
    float64 field plus its representation. Neither identity replaces the other.
    This loads only the small source mask and performs no mechanics operation.
    """
    import numpy as np
    values = np.asarray(np.load(path, allow_pickle=False), dtype="<f8")
    if values.ndim != 1 or not np.isin(values, [0., 1.]).all():
        raise ValueError("Expected a one-dimensional binary source field")
    return hashlib.sha256(values.tobytes()+b"binary").hexdigest()


class Audit:
    def __init__(self):
        self.checks = []
        self.inputs = {}

    def check(self, name, passed, detail=None):
        row = {"name": name, "passed": bool(passed)}
        if detail is not None:
            row["detail"] = detail
        self.checks.append(row)

    def read(self, path):
        path = Path(path).resolve()
        self.inputs[str(path.relative_to(ROOT))] = digest(path)
        return json.loads(path.read_text(encoding="utf-8"))

    def hash(self, path, expected, name):
        path = Path(path).resolve()
        actual = digest(path) if path.is_file() else None
        self.check(name, bool(expected) and actual == expected,
                   {"path": str(path.relative_to(ROOT)), "expected": expected, "actual": actual})

    def section(self, name, function):
        try:
            function()
        except Exception as exc:
            self.check(name+": completed", False,
                       {"exception": type(exc).__name__, "message": str(exc)})


def key(row):
    return row["catalog_key"], float(row["angle_degrees"])


def check_csv(audit, path, rows):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        actual = list(csv.DictReader(stream))
    columns = list(dict.fromkeys(column for row in rows for column in row))
    expected = [{column: "" if row.get(column) is None else str(row[column])
                 for column in columns} for row in rows]
    audit.check("CSV equals JSON: "+Path(path).name, actual == expected)


def run(output, refinement, tables):
    output, refinement, tables = map(lambda p: Path(p).resolve(), (output, refinement, tables))
    if output.exists():
        raise FileExistsError(output)
    audit = Audit()
    baseline_path = ROOT/"results/phase11/preservation_before.json"
    catalogue_path = ROOT/"results/phase11/geometry/catalog.json"
    context = {}

    def preservation():
        baseline = audit.read(baseline_path)
        context["baseline"] = baseline
        audit.check("Original inventory has 341 files", len(baseline["files"]) == 341)
        for relative, item in baseline["files"].items():
            audit.hash(ROOT/relative, item["sha256"], "Original unchanged: "+relative)
            path = ROOT/relative
            audit.check("Original byte count: "+relative, path.is_file() and path.stat().st_size == item["bytes"])

    def provenance():
        catalogue = audit.read(catalogue_path)["entries"]
        context["entries"] = {entry["catalog_key"]: entry for entry in catalogue}
        audit.check("Catalogue contains 13 distinct candidates", len(catalogue) == len(context["entries"]) == 13)
        for entry in catalogue:
            name = entry["catalog_key"]
            source = entry.get("provenance") or {}
            expected = source.get("array_file_sha256", entry["source_input_sha256"])
            audit.hash(ROOT/entry["array_path"], expected, "Analyzed catalogue array: "+name)
            audit.check("Catalogue geometry identity: "+name,
                        entry["geometry_sha256"] == entry["geometry"]["geometry_sha256"])
            if source:
                audit.hash(ROOT/source["original_array_path"], entry["source_input_sha256"], "Repair parent array: "+name)
                audit.check("Repair provenance matches catalogue geometry: "+name,
                            source["geometry_sha256"] == entry["geometry_sha256"])
            else:
                for path_key, hash_key in (("source_mesh_path", "source_mesh_file_sha256"),
                                           ("history_path", "historical_history_sha256")):
                    audit.hash(ROOT/entry[path_key], entry[hash_key], name+": "+path_key)
                audit.hash(ROOT/entry["history_path"].replace("evaluations.json", "verification.json"),
                           entry["historical_label_sha256"], "Historical label file: "+name)
        thresholds = audit.read(ROOT/"results/phase11/geometry/thresholds.json")["entries"]
        context["thresholds"] = thresholds
        audit.check("21 threshold reconstructions", len(thresholds) == 21)
        for entry in thresholds:
            audit.hash(ROOT/entry["array_path"], entry["array_file_sha256"], "Threshold array: "+entry["candidate_id"])
            audit.hash(ROOT/entry["density_source_path"], entry["source_input_sha256"], "Threshold density source: "+entry["candidate_id"])
        generation = audit.read(ROOT/"results/phase11/geometry/generation_cost.json")
        for relative, expected in generation["input_files_sha256"].items():
            audit.hash(ROOT/relative, expected, "Geometry generation input: "+relative)
        for relative, expected in generation["artifact_files_sha256"].items():
            audit.hash(ROOT/relative, expected, "Geometry archived artifact: "+relative)

    def labels():
        entries = context["entries"]
        aggregate = audit.read(refinement/"labels.json")
        archived = {key(row): row["label"] for row in aggregate}
        expected_keys = {(name, angle) for name in entries for angle in ANGLES}
        audit.check("39 unique aggregate labels", len(aggregate) == 39 and set(archived) == expected_keys)
        recomputed = {}
        for name, entry in entries.items():
            folder = refinement/"candidates"/name
            histories = audit.read(folder/"responses.json")
            saved = audit.read(folder/"labels.json")
            audit.check("Three response angles: "+name, set(histories) == set(saved) == {str(a) for a in ANGLES})
            for angle in ANGLES:
                label = reference_label(histories[str(angle)], historical_label=entry.get("historical_label"), geometry_audit=entry["geometry"])
                recomputed[name, angle] = label
                audit.check(f"Recomputed label equals both archives: {name}, {angle:g}",
                            label == saved[str(angle)] == archived.get((name, angle)))
                audit.check(f"Historical label unchanged: {name}, {angle:g}",
                            label["historical_label"] == entry.get("historical_label") and label["historical_label_modified"] is False)
                eligible = (label["passed"] and label["numerical_status"] == "numerically_accepted"
                            and label["geometry_status"] == "geometry_valid" and label["constraint_status"] == "within_budget")
                audit.check(f"Eligibility conjunction: {name}, {angle:g}", label["eligible_reference"] == eligible)
                audit.check(f"No unresolved numerical reference: {name}, {angle:g}",
                            label["numerical_status"] != "unresolved" or label["reference_estimate"] is None)
                audit.check(f"No physical or exact-truth label: {name}, {angle:g}",
                            label["physical_validation"] == "not_assessed" and label["finest_is_exact_truth"] is False)
        context["labels"] = recomputed

    def table_checks():
        rows = audit.read(tables/"candidate_results.json")
        indexed = {key(row): row for row in rows}
        audit.check("Candidate table has 39 unique rows", len(rows) == 39 and set(indexed) == set(context["labels"]))
        for identity, label in context["labels"].items():
            row = indexed[identity]
            for field in ("numerical_status", "geometry_status", "constraint_status", "classification", "eligible_reference"):
                audit.check(f"Table label {field}: {identity}", row[field] == label[field])
            estimate = label["reference_estimate"] or {}
            audit.check(f"Table computed reference only: {identity}",
                        row["reference_mesh_h_mm"] == estimate.get("mesh_h_mm")
                        and row["reference_compliance_n_mm"] == estimate.get("compliance_n_mm")
                        and row["reference_pad_mm"] == estimate.get("max_pad_displacement_mm"))
            audit.check(f"Shared-angle cost scope explicit: {identity}",
                        "shared across three angles" in row["cost_scope"] and "do not sum" in row["cost_scope"])
        imports = audit.read(refinement/"imported_evidence.json")
        costs = audit.read(refinement/"costs.json")
        context["costs"] = costs
        for name in context["entries"]:
            historical = sum(item.get("historical_wall_seconds") or 0 for item in imports if item["catalog_key"] == name)
            new = sum(item["wall_seconds"] for item in costs["attempts"] if item["catalog_key"] == name)
            audit.check("Candidate cost counted once per basis: "+name,
                        all(indexed[name, angle]["historical_evaluation_wall_seconds"] == historical
                            and indexed[name, angle]["new_attempt_wall_seconds"] == new for angle in ANGLES))
        thresholds = audit.read(tables/"threshold_results.json")
        audit.check("21 threshold table rows", len(thresholds) == 21)
        audit.check("Thresholds have no invented mechanics reference", all(row["numerical_status"] == "not_evaluated" and row["eligible_reference"] is False for row in thresholds))
        audit.check("Threshold table identities", {r["candidate_id"] for r in thresholds} == {r["candidate_id"] for r in context["thresholds"]})
        check_csv(audit, tables/"candidate_results.csv", rows)
        check_csv(audit, tables/"threshold_results.csv", thresholds)

    def requests_and_costs():
        plan = audit.read(refinement/"plan.json")
        audit.hash(catalogue_path, plan["catalogue_sha256"], "Execution-plan catalogue")
        audit.hash(ROOT/"research/REFERENCE_LABEL_PROTOCOL.md", plan["protocol_sha256"], "Execution-plan protocol")
        costs = context.get("costs") or audit.read(refinement/"costs.json")
        ledger = audit.read(tables/"evaluation_costs.json")
        derived = {(row["catalog_key"], row["mesh_h_mm"]): row for row in ledger if row["category"] == "new_refinement"}
        attempts = {(row["catalog_key"], row["mesh_h_mm"]): row for row in costs["attempts"]}
        audit.check("One derived ledger row per attempt", len(derived) == len(costs["attempts"]) and set(derived) == set(attempts))
        requests = list((refinement/"attempts").glob("*/*/request.json"))
        audit.check("Every attempt has one immutable request", len(requests) == len(attempts))
        seen = set()
        for path in requests:
            request = audit.read(path)
            entry = request["entry"]
            identity = entry["catalog_key"], request["mesh_h_mm"]
            seen.add(identity)
            audit.check("Request catalogue entry: "+str(identity), entry == context["entries"][identity[0]])
            source_hash = (entry.get("provenance") or {}).get("array_file_sha256", entry["source_input_sha256"])
            audit.hash(ROOT/entry["array_path"], source_hash, "Request analyzed array: "+str(identity))
            audit.hash(ROOT/"research/REFERENCE_LABEL_PROTOCOL.md", request["protocol_sha256"], "Request protocol: "+str(identity))
            for relative, expected in request["code_sha256"].items():
                audit.hash(ROOT/relative, expected, "Request code: "+str(identity)+": "+relative)
            cost, row = attempts[identity], derived[identity]
            # A preflight deferral has no worker and may omit record_path;
            # its skeleton record is still archived beside request.json.
            record_path = ROOT/cost["record_path"] if cost.get("record_path") else path.parent/"record.json"
            record = audit.read(record_path)
            expected_field_hash = binary_field_hash(ROOT/entry["array_path"])
            audit.check("Attempt frozen-field identity: "+str(identity),
                        record["hash"] == expected_field_hash and record["mesh_h_mm"] == identity[1])
            if not entry.get("is_new_geometry"):
                audit.check("Original catalogue field identity: "+str(identity),
                            expected_field_hash == entry["candidate_hash"])
            if "analyzed_array_sha256" in record:
                audit.check("Worker analyzed digest: "+str(identity), record["analyzed_array_sha256"] == source_hash)
                audit.check("Worker parent-input digest: "+str(identity), record["source_input_sha256"] == entry["source_input_sha256"])
            audit.check("Attempt wall time is not double charged: "+str(identity), row["wall_seconds"] == cost["wall_seconds"])
            exc = record.get("exception", {})
            residual = exc.get("type") == "ArithmeticError" and exc.get("message") == RESIDUAL_FAILURE
            # This archived traceback is after splu returned and before the
            # P2 factor.solve call. Charge the completed P2 LU only; no P2 RHS
            # completion is implied by allocation failure during pivot access.
            pivot_memory = (exc.get("type") == "MemoryError"
                            and "factor.U.diagonal()" in exc.get("traceback", "")
                            and "verification.py" in exc.get("traceback", ""))
            if residual:
                audit.check("Rejected P2 completed-stage charge: "+str(identity),
                            cost["factorizations"] == 1 and cost["rhs_solves"] == 2
                            and row["known_extra_p2_factorizations"] == 1 and row["known_extra_p2_rhs_solves"] == 2
                            and row["charged_factorizations"] == 2 and row["charged_rhs_solves"] == 4
                            and record["status"] != "evaluated")
            elif pivot_memory:
                audit.check("P2 pivot-allocation failure completed-stage charge: "+str(identity),
                            cost["factorizations"] == 1 and cost["rhs_solves"] == 2
                            and row["known_extra_p2_factorizations"] == 1 and row["known_extra_p2_rhs_solves"] == 0
                            and row["charged_factorizations"] == 2 and row["charged_rhs_solves"] == 2
                            and record["status"] != "evaluated")
            elif exc:
                audit.check("Other exception counts remain unknown: "+str(identity),
                            row["charged_factorizations"] is None and row["charged_rhs_solves"] is None)
            else:
                audit.check("Unmodified known stage counts: "+str(identity),
                            row["charged_factorizations"] == cost["factorizations"] and row["charged_rhs_solves"] == cost["rhs_solves"])
            if cost["status"] == "resource_preflight_deferred":
                audit.check("Deferred attempt has exactly zero solves: "+str(identity),
                            cost["factorizations"] == cost["rhs_solves"] == row["charged_factorizations"] == row["charged_rhs_solves"] == 0
                            and request["preflight"]["allowed"] is False)
        audit.check("Request identities cover all attempts", seen == set(attempts))
        validation = audit.read(ROOT/"results/phase11/validation_replay.json")["measured_cost"]
        validation_rows = [row for row in ledger if row["category"] == "validation_replay"]
        audit.check("Measured validation costs retained", len(validation_rows) == 1 and
                    validation_rows[0]["charged_factorizations"] == validation["factorizations"] and
                    validation_rows[0]["charged_rhs_solves"] == validation["rhs_solves"] and
                    validation_rows[0]["wall_seconds"] == validation["wall_seconds"])
        failed_path = ROOT/"results/phase11/kinematics/failure.json"
        if failed_path.exists():
            failure = audit.read(failed_path)
            failures = [row for row in ledger if row["category"] == "kinematics_csv_failure"]
            audit.check("Failed-exporter measured/unknown timing preserved", len(failures) == 1 and failures[0]["wall_seconds"] == failure["wall_seconds"])
        check_csv(audit, tables/"evaluation_costs.csv", ledger)

    def links():
        baseline = context.get("baseline") or audit.read(baseline_path)
        required = {ROOT/"research"/name for name in REQUIRED_DOCUMENTS}
        documents = required | {path for path in (ROOT/"research").rglob("*.md") if path.relative_to(ROOT).as_posix() not in baseline["files"]}
        count = 0
        for document in sorted(documents):
            audit.check("Required/new report exists: "+document.name, document.is_file())
            if not document.is_file():
                continue
            audit.inputs[str(document.relative_to(ROOT))] = digest(document)
            text = document.read_text(encoding="utf-8")
            text = re.sub(r"(?ms)^```[^\n]*\n.*?^```\s*$", "", text)
            targets = re.findall(r"!?\[[^\]\n]*\]\((<[^>]+>|[^)\s]+)(?:\s+[^)]*)?\)", text)
            targets += re.findall(r"(?m)^\s*\[[^\]]+\]:\s*(<[^>]+>|\S+)", text)
            for target in targets:
                target = target.strip("<>")
                if target.startswith(("#", "http://", "https://", "mailto:", "app:", "codex:", "data:")):
                    continue
                local = unquote(target.split("#", 1)[0].split("?", 1)[0])
                if not local:
                    continue
                path = (document.parent/local).resolve()
                # A report may link this audit artifact before its exclusive
                # write below. No other absent artifact receives an exemption.
                exists = path.exists() or path == output
                audit.check("Local Markdown link: "+document.name+" -> "+target, exists)
                count += 1
        audit.check("New reports contain local evidence links", count > 0, {"links_checked": count})

    for name, function in (("preservation", preservation), ("provenance", provenance),
                           ("labels", labels), ("tables", table_checks),
                           ("requests_and_costs", requests_and_costs), ("markdown_links", links)):
        audit.section(name, function)
    passed = all(row["passed"] for row in audit.checks)
    report = {"schema_version": "prometheus.phase11.evidence-audit.v1",
              "checked_utc": datetime.now(timezone.utc).isoformat(), "passed": passed,
              "scope": "Read-only file integrity, scalar label recomputation, table/cost consistency and local report-link existence; no mechanics or physical validation",
              "factorizations": 0, "rhs_solves": 0, "new_mechanics_performed": False,
              "script_sha256": digest(Path(__file__)), "inputs_sha256": audit.inputs,
              "checks_total": len(audit.checks),
              "checks_failed": sum(not row["passed"] for row in audit.checks), "checks": audit.checks}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"passed": passed, "checks": len(audit.checks), "failed": report["checks_failed"], "output": str(output)}))
    if not passed:
        raise AssertionError("Phase 1.1 evidence audit failed; report retained at "+str(output))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/phase11/evidence_audit.json")
    parser.add_argument("--refinement", default="results/phase11/refinement")
    parser.add_argument("--tables", default="results/phase11/tables")
    args = parser.parse_args()
    run(args.output, args.refinement, args.tables)
