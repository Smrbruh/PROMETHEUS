"""Refinement orchestration tests with mocked mechanics and subprocesses.

These checks exercise archival/cost handling only; they do not run a FEM
factorization or allocate any fine mesh.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from prometheus.geometry import make_rect_mesh


@pytest.fixture
def runner(tmp_path, monkeypatch):
    source = Path(__file__).resolve().parents[1]/"scripts"/"phase11_refine.py"
    spec = importlib.util.spec_from_file_location("phase11_refine_test_module", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    script = tmp_path/"scripts"/"phase11_refine.py"
    script.parent.mkdir()
    script.write_text("# mocked runner identity\n", encoding="utf-8")
    monkeypatch.setattr(module, "__file__", str(script))
    (tmp_path/"research").mkdir()
    (tmp_path/"research"/"REFERENCE_LABEL_PROTOCOL.md").write_text("Mock frozen protocol\n", encoding="utf-8")
    (tmp_path/"prometheus").mkdir()
    # The source topology is one quad; no source or fine benchmark is solved.
    monkeypatch.setattr(module, "make_l_mesh", lambda h: make_rect_mesh(2., 2., 2.))
    monkeypatch.setattr(module, "peak_memory_bytes", lambda: 123456)
    return module


def entry_for(runner, *, repaired=False):
    path = runner.ROOT/"mask.npy"
    np.save(path, np.ones(1))
    actual_digest = runner.digest(path)
    entry = {"candidate_name": "mock_binary", "catalog_key": "mock_v040",
             "array_path": "mask.npy", "volume_setting": .4,
             "source_input_sha256": actual_digest,
             "geometry": {"pads_preserved": True}}
    if repaired:
        entry["source_input_sha256"] = "parent-density-array-digest"
        entry["provenance"] = {"array_file_sha256": actual_digest}
    return entry, actual_digest


def write_request(runner, entry, h=.25):
    out = runner.ROOT/"attempt"
    out.mkdir()
    request = out/"request.json"
    request.write_text(json.dumps({"entry": entry, "mesh_h_mm": h}), encoding="utf-8")
    return request


def test_preflight_rejection_archives_zero_solves_without_launching_child(runner, monkeypatch):
    entry, _ = entry_for(runner)
    estimate = {"allowed": False, "reasons": ["Declared assembly cap"],
                "p2_dofs": 2000000, "q4_cells": 250000,
                "estimated_p2_coo_bytes": 5*1024**3, "matrix_assembled": False}
    monkeypatch.setattr(runner, "resource_estimate", lambda *args: estimate)

    def forbidden(*args, **kwargs):
        pytest.fail("Rejected preflight must not launch a worker or invoke mechanics")

    monkeypatch.setattr(runner.subprocess, "Popen", forbidden)
    monkeypatch.setattr(runner, "evaluate_candidate", forbidden)
    out = runner.ROOT/"deferred"
    record, arrays, cost = runner.attempt(entry, .125, out)
    assert record["status"] == cost["status"] == "resource_preflight_deferred"
    assert record["factorizations"] == record["rhs_solves"] == 0
    assert cost["factorizations"] == cost["rhs_solves"] == 0
    assert arrays == {} and cost["peak_working_set_bytes"] == 0
    assert "compliance_n_mm" not in record
    assert record["resource_preflight"] == estimate
    archived = json.loads((out/"request.json").read_text())
    assert archived["preflight"] == estimate
    assert archived["policy"]["max_p2_dofs"] == 1500000
    assert json.loads((out/"record.json").read_text()) == record
    assert not (out/"worker.log").exists()
    assert not (out/"fields.npz").exists()


def test_worker_preserves_completed_primary_evidence_after_independent_exception(runner, monkeypatch):
    entry, expected_digest = entry_for(runner, repaired=True)
    partial = runner.skeleton(entry, .25)
    partial.update(status="independent_verification_failed", solver_passed=True,
                   factorizations=1, rhs_solves=2, compliance_n_mm=42.,
                   diagnostics={"passed": True},
                   independent_failure={"type": "RuntimeError", "message": "mock P2 failure"})
    fields = {"nodes": np.array([[0., 80.], [120., 20.]]),
              "u_basis": np.arange(8., dtype=float).reshape(4, 2)}

    def fail_independent(*args):
        assert args[0].name == entry["candidate_name"]
        assert args[-3:] == (.25, True, False)
        exc = RuntimeError("mock P2 failure")
        exc.partial_evaluation = partial
        exc.partial_arrays = fields
        raise exc

    monkeypatch.setattr(runner, "evaluate_candidate", fail_independent)
    request = write_request(runner, entry)
    runner.worker(request)
    saved = json.loads((request.parent/"record.json").read_text())
    assert saved["status"] == "independent_verification_failed"
    assert saved["solver_passed"] is True
    assert saved["factorizations"] == 1 and saved["rhs_solves"] == 2
    assert saved["compliance_n_mm"] == 42.
    assert saved["exception"]["type"] == "RuntimeError"
    assert "mock P2 failure" in saved["exception"]["traceback"]
    assert saved["analyzed_array_sha256"] == expected_digest
    assert saved["source_input_sha256"] == "parent-density-array-digest"
    with np.load(request.parent/"fields.npz") as restored:
        assert set(restored.files) == set(fields)
        for key, value in fields.items():
            np.testing.assert_array_equal(restored[key], value)
    completion = json.loads((request.parent/"worker_completion.json").read_text())
    assert completion["fields_saved"] is True
    assert completion["record_sha256"] == runner.digest(request.parent/"record.json")
    assert completion["peak_memory_bytes"] == 123456


def test_worker_exception_before_record_creation_does_not_invent_responses(runner, monkeypatch):
    entry, _ = entry_for(runner)

    def fail_early(*args):
        raise ValueError("mock early failure")

    monkeypatch.setattr(runner, "evaluate_candidate", fail_early)
    request = write_request(runner, entry)
    runner.worker(request)
    saved = json.loads((request.parent/"record.json").read_text())
    assert saved["status"] == "solver_exception" and saved["solver_passed"] is False
    assert saved["factorizations"] == saved["rhs_solves"] == 0
    assert "compliance_n_mm" not in saved and "independent" not in saved
    assert not (request.parent/"fields.npz").exists()
    assert json.loads((request.parent/"worker_completion.json").read_text())["fields_saved"] is False


def test_changed_input_digest_prevents_mechanics(runner, monkeypatch):
    entry, _ = entry_for(runner)
    entry["source_input_sha256"] = "unexpected"

    def forbidden(*args, **kwargs):
        pytest.fail("Input identity failure must be detected before mechanics")

    monkeypatch.setattr(runner, "evaluate_candidate", forbidden)
    request = write_request(runner, entry)
    with pytest.raises(ValueError, match="digest"):
        runner.worker(request)
    assert not (request.parent/"record.json").exists()
