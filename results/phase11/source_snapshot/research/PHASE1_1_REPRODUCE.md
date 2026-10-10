# Phase 1.1 reproduction and evidence map

Run from the repository root using the Phase 1 Python/version pins in [requirements-tested.txt](../requirements-tested.txt). No new dependencies or modified Phase 1 kernels are needed. Every destination below must be new. The shipped `results/phase11` evidence must be retained; substitute a fresh replay name when repeating a command.

## Baseline and preservation

The initial [inventory](../results/phase11/preservation_before.json) contains 341 preexisting source, test, specification, result and failed-attempt files. `scripts/phase11_preservation.py check` compares their bytes without rewriting them. Do not recapture this inventory after adding files. Runtime caches are excluded.

```powershell
python -B -m pytest -p no:cacheprovider -q tests/test_envelope.py tests/test_evaluation.py tests/test_fem.py tests/test_optimization.py tests/test_run.py tests/test_verification.py --junitxml=results/phase11/tests_baseline_replay.xml
python -B scripts/phase11_reproduce.py --output results/phase11/baseline_replay_again
python -B scripts/report_evidence.py --runs results/phase11/baseline_replay_again/replay_v040 results/phase11/baseline_replay_again/replay_v035 results/phase11/baseline_replay_again/replay_v045 --output results/phase11/figures_replay_again
python -B scripts/phase11_preservation.py check --output results/phase11/preservation_replay_check.json
```

The replay driver executes all four documented numerical jobs, including the h1 reoptimization, and the original validation script. It compares optimization histories, evaluation responses and diagnostics, stopping causes, failure records, labels and saved `.npy` arrays. Timing/memory fields are excluded from numerical equality. It resumes only already-completed replay jobs; an incomplete job is preserved and requires a new output identity.

## Added audit

```powershell
python -B -m pytest -p no:cacheprovider -q --junitxml=results/phase11/tests_all_replay.xml
python -B scripts/phase11_geometry.py --output results/phase11/geometry_replay
python -B scripts/phase11_refine.py --output results/phase11/refinement_replay
python -B scripts/phase11_kinematics.py --output results/phase11/kinematics_replay
python -B scripts/phase11_repair_kinematics.py --output results/phase11/repair_kinematics_replay
python -B scripts/phase11_tables.py --refinement results/phase11/refinement_replay --output results/phase11/tables_replay
python -B scripts/phase11_figures.py --output results/phase11/figures_replay
```

The refinement and table scripts deliberately use the shipped immutable `results/phase11/geometry/catalog.json` as the authoritative candidate catalogue. The independent geometry replay demonstrates reconstruction reproducibility without replacing that catalogue. The kinematic script uses the archived original h=.5 fields. The table script combines the selected refinement directory with the shipped baseline replay and geometry/kinematic cost records; it is not a fresh timing of every stage. A new end-to-end timing campaign must identify its chosen input bundles explicitly.

Numerical workers run sequentially with single-threaded backend environment variables. Each child has a 600-second wall-time limit and an observed 8-GiB working-set stop. Before launch, the unchanged pixel topology predicts P2 degrees of freedom and the original conservative assembly allocation; the 1.5-million DOF and 4-GiB assembly limits are applied without changing tolerances. A resource deferral means no solve, not convergence. Parent monitoring is polling, not an operating-system hard memory guarantee.

The [protocol](REFERENCE_LABEL_PROTOCOL.md) was written before the fine calculations. The [committed execution plan](../results/phase11/refinement/plan.json), also written before those calculations, clarifies one scope ambiguity: **any unresolved declared angle can trigger refinement of a connected candidate, including B0/B1**. Historical nominal acceptance never shortcuts a new-angle check. No threshold or acceptance tolerance changes. Stop after all three angle audits pass, the first failed/resource-limited attempt, or exhaustion of the declared sequence. Every terminal failure remains in its candidate's label history.

## Reading the outputs

- [Catalogue](../results/phase11/geometry/catalog.json), [threshold masks](../results/phase11/geometry/thresholds.json) and [repair](../results/phase11/geometry/repair.json): immutable identities, source byte hashes, mask paths, geometry diagnostics and generation costs.
- [All mesh responses](../results/phase11/refinement/responses.csv): original h2/h1/h.5 evidence and every new attempt, including failures. Full projected records and labels are under each candidate directory.
- [Candidate summary](../results/phase11/tables/candidate_results.csv): separate numerical, geometry, budget and eligibility predicates, with empirical reference estimates only where justified. Cost columns repeat across angles because the same two solved basis fields are reused; do not sum those repeated columns.
- [Threshold summary](../results/phase11/tables/threshold_results.csv): geometry-only results; no invented mechanics labels.
- [Evaluation costs](../results/phase11/tables/evaluation_costs.csv): completed runs, failed attempts, stage-supported solver accounting and explicitly unknown measurements.
- [Kinematic fields and summary](../results/phase11/kinematics/complete/manifest.json): archived-field diagnostics, separate density-region context, zero new FEM solves.

For a P2 residual exception, the original evaluator returns only the completed Q4 counters even though the traceback shows P2 LU and both RHS solves completed before rejection. The derived cost ledger adds that demonstrable work and preserves the original counters. Rejection never becomes acceptance. Unmeasured failed exporter time remains null rather than being fabricated as zero. The validation replay separately records 14 factorizations and 14 RHS solves.

For the recorded allocation failure at P2 `factor.U` inspection, the traceback establishes that P2 factorization returned but its RHS solve had not started. The ledger adds one factorization and zero RHS columns. A monitored child termination without a returned record has unknown completed-stage counts. In every case the complete measured attempt elapsed time is charged once; stage corrections add counts, not duplicate time.

All results concern the specified synthetic, small-deformation, isotropic plane-stress model. A numerical reference is an empirically resolved response for that model, not exact continuum truth or physical validation.
