"""Generate the frozen Phase 1.1 geometry catalogue without mechanics solves.

Run: python scripts/phase11_geometry.py --output results/phase11/geometry
Existing nonempty output directories are refused. Historical artifacts are
opened read-only, hashed before/after, and referenced under their original IDs.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from prometheus.binary_audit import (THRESHOLDS, binary_diagnostics,
    threshold_sweep, retain_pad_component)
from prometheus.geometry import Mesh
from prometheus.models import Material


RUNS = ((.35, "results/phase1_v035"), (.40, "results/phase1_v040_verified"),
        (.45, "results/phase1_v045"))
ORIGINAL_NAMES = ("B0", "B1", "B3_binary", "B3_threshold")


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode("utf-8")).hexdigest()


def relative(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def write_json_new(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as output:
        json.dump(value, output, indent=2, allow_nan=False)
        output.write("\n")


def save_array_new(path, value):
    with Path(path).open("xb") as output:
        np.save(output, value, allow_pickle=False)


def build_catalogue(output):
    started, cpu = perf_counter(), process_time()
    output = Path(output).resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("Refusing to overwrite existing geometry evidence: "+str(output))
    output.mkdir(parents=True, exist_ok=True)
    arrays_dir = output/"arrays"
    arrays_dir.mkdir()
    inputs = {}

    def track(path):
        path = Path(path)
        key = relative(path)
        digest = file_hash(path)
        if key in inputs and inputs[key] != digest:
            raise RuntimeError("Historical input changed during generation: "+key)
        inputs[key] = digest
        return path, digest

    def load_json(path):
        tracked, digest = track(path)
        return json.loads(tracked.read_text(encoding="utf-8")), digest

    def load_array(path):
        tracked, digest = track(path)
        return np.load(tracked, allow_pickle=False), digest

    code_protocol = {relative(path): file_hash(path) for path in (
        Path(__file__), ROOT/"prometheus/binary_audit.py",
        ROOT/"prometheus/geometry.py", ROOT/"prometheus/models.py",
        ROOT/"research/REFERENCE_LABEL_PROTOCOL.md")}
    catalogue, families, aliases = [], {}, {}
    for fraction, run_relative in RUNS:
        run = ROOT/run_relative
        manifest, manifest_hash = load_json(run/"manifest.json")
        if manifest["benchmark"]["volume_fraction"] != fraction:
            raise ValueError("Archived volume setting disagrees with frozen catalogue")
        material = Material(**manifest["material"])
        mesh_path, mesh_hash = track(run/"source_mesh.npz")
        with np.load(mesh_path, allow_pickle=False) as data:
            mesh = Mesh(data["nodes"], data["cells"], data["cell_indices"], float(data["h_mm"]))
        if mesh.h != 2.:
            raise ValueError("Phase 1.1 reconstruction catalogue requires the frozen h2 source")
        family_key = f"v{round(fraction*100):03d}"
        families[family_key] = {"run": run, "mesh": mesh, "material": material,
                                "volume_setting": fraction, "manifest_sha256": manifest_hash}
        for name in ORIGINAL_NAMES:
            generation_start, generation_cpu = perf_counter(), process_time()
            directory = run/"candidates"/name
            array_path = directory/"physical_source.npy"
            binary, input_hash = load_array(array_path)
            history_path = directory/"evaluations.json"
            history, history_hash = load_json(history_path)
            historical, label_hash = load_json(directory/"verification.json")
            if not history or any(record["candidate"] != name for record in history):
                raise ValueError("Historical candidate identity mismatch")
            expected_hash = hashlib.sha256(binary.astype("<f8").tobytes()+b"binary").hexdigest()
            if any(record["hash"] != expected_hash for record in history):
                raise ValueError("Historical binary field hash mismatch")
            diagnostics = binary_diagnostics(mesh, binary, material, fraction)
            key = family_key+"_"+name
            identity = {"catalog_key": key, "candidate_name": name, "candidate_hash": expected_hash,
                        "volume_setting": fraction, "source_input_sha256": input_hash,
                        "source_mesh_file_sha256": mesh_hash, "manifest_sha256": manifest_hash,
                        "historical_history_sha256": history_hash, "historical_label_sha256": label_hash,
                        "geometry_sha256": diagnostics["geometry_sha256"]}
            row = {**identity, "catalog_identity_sha256": canonical_hash(identity),
                   "array_path": relative(array_path), "history_path": relative(history_path),
                   "source_mesh_path": relative(mesh_path), "original_run": run_relative,
                   "historical_label": historical, "geometry": diagnostics,
                   "is_new_geometry": False, "material": asdict(material),
                   "generation_wall_seconds": perf_counter()-generation_start,
                   "generation_cpu_seconds": process_time()-generation_cpu,
                   "factorizations": 0, "rhs_solves": 0}
            aliases.setdefault(diagnostics["geometry_sha256"], "catalog-"+row["catalog_identity_sha256"])
            catalogue.append(row)

    threshold_rows = []
    for family_key, family in families.items():
        source_path = family["run"]/"candidates"/"B3"/"physical_source.npy"
        physical, source_hash = load_array(source_path)
        generated = threshold_sweep(family["mesh"], physical, source_hash,
                                    family["material"], family["volume_setting"], aliases=aliases)
        for binary, record in generated:
            array_path = arrays_dir/(family_key+"_threshold_"+record["candidate_id"][7:23]+".npy")
            save_array_new(array_path, binary)
            record.update(volume_setting=family["volume_setting"], original_run=relative(family["run"]),
                          density_source_path=relative(source_path), array_path=relative(array_path),
                          array_file_sha256=file_hash(array_path))
            threshold_rows.append(record)
            aliases.setdefault(record["geometry_sha256"], record["candidate_id"])

    family = families["v035"]
    original_path = family["run"]/"candidates"/"B3_binary"/"physical_source.npy"
    original, original_hash = load_array(original_path)
    repaired, repair = retain_pad_component(family["mesh"], original, original_hash,
                                            family["material"], family["volume_setting"])
    repair_path = arrays_dir/("v035_repair_"+repair["candidate_id"][7:23]+".npy")
    save_array_new(repair_path, repaired)
    repair.update(catalog_key="v035_repair", volume_setting=.35, original_run=relative(family["run"]),
                  original_array_path=relative(original_path), array_path=relative(repair_path),
                  array_file_sha256=file_hash(repair_path))
    identity = {"catalog_key": "v035_repair", "candidate_name": repair["candidate_id"],
                "candidate_hash": repair["candidate_id"], "source_input_sha256": original_hash,
                "geometry_sha256": repair["geometry_sha256"], "volume_setting": .35}
    catalogue.append({**identity, "catalog_identity_sha256": canonical_hash(identity),
                      "array_path": relative(repair_path), "history_path": None,
                      "source_mesh_path": relative(family["run"]/"source_mesh.npz"),
                      "original_run": relative(family["run"]), "historical_label": None,
                      "geometry": repair["diagnostics"], "is_new_geometry": True,
                      "provenance": repair, "material": asdict(family["material"]),
                      "generation_wall_seconds": repair["generation_wall_seconds"],
                      "generation_cpu_seconds": repair["generation_cpu_seconds"],
                      "factorizations": 0, "rhs_solves": 0})
    identities = [row["catalog_identity_sha256"] for row in catalogue]
    catalogue_hash = canonical_hash(identities)
    catalog_document = {"schema_version": "prometheus.phase11.geometry.v1", "entries": catalogue,
                        "catalogue_identity_sha256": catalogue_hash,
                        "input_files_sha256": inputs, "code_protocol_sha256": code_protocol}
    write_json_new(output/"catalog.json", catalog_document)
    write_json_new(output/"thresholds.json", {"thresholds": list(THRESHOLDS), "entries": threshold_rows,
                                             "mechanics_evaluated": False,
                                             "code_protocol_sha256": code_protocol})
    write_json_new(output/"repair.json", repair)
    for name, expected in inputs.items():
        if file_hash(ROOT/name) != expected:
            raise RuntimeError("Historical input changed during generation: "+name)
    artifact_hashes = {relative(path): file_hash(path) for path in sorted(output.rglob("*")) if path.is_file()}
    cost = {"schema_version": "prometheus.phase11.geometry-cost.v1",
            "original_geometry_audits": 12, "threshold_reconstructions": len(threshold_rows),
            "explicit_repair_operations": 1, "new_archived_masks": len(threshold_rows)+1,
            "factorizations": 0, "rhs_solves": 0,
            "generation_wall_seconds": perf_counter()-started,
            "generation_cpu_seconds": process_time()-cpu,
            "timing_scope": "input loading through preceding artifact writes and preservation checks; excludes startup/final cost-file write",
            "catalogue_identity_sha256": catalogue_hash, "historical_inputs_unchanged": True,
            "input_files_sha256": inputs, "code_protocol_sha256": code_protocol,
            "artifact_files_sha256": artifact_hashes}
    write_json_new(output/"generation_cost.json", cost)
    return {"output": relative(output), "catalogue_entries": len(catalogue),
            "threshold_masks": len(threshold_rows), "repair_removed_area_mm2": repair["removed_area_mm2"],
            "generation_wall_seconds": cost["generation_wall_seconds"],
            "factorizations": 0, "rhs_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(ROOT/"results/phase11/geometry"))
    args = parser.parse_args()
    print(json.dumps(build_catalogue(args.output), indent=2))


if __name__ == "__main__":
    main()
