# PROMETHEUS: numerical research baseline

PROMETHEUS-01 is a small Python plane-stress elasticity and SIMP/OC benchmark. It contains no AI selector, service or product interface. Material, fixtures and loads are synthetic research assumptions; numerical verification does not establish physical safety or manufacturability.

Start with [Phase 1 decision](research/PHASE1_DECISION.md), [measured results](research/BASELINE_RESULTS.md) and [numerical validation](research/NUMERICAL_VALIDATION.md). The [Phase 0 research package](research/README.md) remains a historical record of the specification and prior-art review.

Phase 1 completed on 2026-10-11 local time: **104 automated tests passed**. The canonical optimization reproduces exactly. Density models pass the declared nominal numerical audits; reconstructed binary designs have geometry or resolution failures. Phase 2 remains deferred. See the reports for the setting-specific outcomes.

## Reproduce

Run from this directory with Python 3.14.4 and the versions in [requirements-tested.txt](requirements-tested.txt). NumPy, SciPy, pytest and matplotlib were already available for this run; no external solver is required. On a fresh environment:

```powershell
python -m pip install -r requirements-tested.txt
python -m pytest -q
python scripts/validation_evidence.py --output results/validation_replay.json
python -m prometheus.run --output results/replay_v040 --volume 0.40 --repeats 2 --audit-derivatives
python -m prometheus.run --output results/replay_v035 --volume 0.35 --repeats 1
python -m prometheus.run --output results/replay_v045 --volume 0.45 --repeats 1
python scripts/mesh_dependence.py --output results/replay_h1
python scripts/report_evidence.py --runs results/replay_v040 results/replay_v035 results/replay_v045 --output results/replay_figures
```

Run the experiments sequentially. Evidence directories must be new or empty; existing measurements are preserved. Timing and bitwise reproducibility are machine-specific. Tests use analytical solutions, manufactured solutions and independent assembly; they do not reproduce a physical experiment.

## Layout

| Path | Purpose |
|---|---|
| `prometheus/` | Geometry, Q4 elasticity, filtering/SIMP/OC, exact angular envelope, P2 verification, evaluation and runner |
| `tests/` | Analytical, derivative, convergence, invalid-input, acceptance and runner checks |
| `scripts/` | Recorded validation evidence, optimization mesh-dependence study and result plotting |
| `results/` | Immutable run inputs/snapshots, numerical fields, histories, diagnostics and figures |
| `research/` | Original Phase 0 references and specification plus Phase 1 reports |

Use [IMPLEMENTATION_REPORT.md](research/IMPLEMENTATION_REPORT.md) for interfaces and assumptions and [EXPERIMENT_PROTOCOL.md](research/EXPERIMENT_PROTOCOL.md) for comparison rules. No project license has been selected; the implementation does not vendor top88 source.

The shipped evidence can be checked without rerunning mechanics using `python scripts/check_phase1_evidence.py`. This checks the named archived run directories, source/specification hashes, arrays, costs, status consistency and report links; it does not validate physical behavior. See [evidence_audit.json](results/evidence_audit.json).
