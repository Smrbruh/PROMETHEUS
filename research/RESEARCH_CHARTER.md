# PROMETHEUS Phase 0 research charter
Date: 2026-10-10 (Asia/Qyzylorda). Status: evidence review and proposed experiments, not a validated engineering system.

## Mission and decision
Determine whether a small, reproducible mechanics research program can contribute beyond mature topology optimization, reliability analysis and adaptive experimentation. The current decision is **NARROW**. Build confidence in the benchmark and measure whether there is a problem worth a new method before attempting one.

PROMETHEUS-01, “The Impossible Bracket”, is a benchmark nickname, not a claim that the part is impossible, novel, safe or manufacturable. Its first version is an in-plane L-shaped plate. A physical bracket with joints, out-of-plane loads and manufacturing variability is a later, separately specified problem.

## Scope
Included: primary-source literature and product documentation; explicit prior-art threats; reproducible linear-elastic baseline specification; PC feasibility; mathematical checks; independent verification strategy; candidate hypotheses and rejection gates; architecture recommendations.

Excluded: full application, UI, autonomous procurement or testing, trained foundation model, GPU/distributed infrastructure, production CAD workflow, certified material allowables, legal patent conclusions, and any assertion of structural safety. No production software exists or was added.

The research engine must not depend on paid software. Existing solvers and optimizers should be reused where their interfaces, licenses and validation are suitable. Source availability does not establish permission to redistribute it.

## Scientific rules
1. Mark statements as established method, source-reported finding, proposed mechanism, local algebra check, or demonstrated bracket result.
2. Preserve failures, unconverged solves and invalid geometry in denominators and logs.
3. Separate optimization physics from independently verified geometry and observations.
4. Calibrate an algorithm using development cases; freeze it before held-out evaluation.
5. Charge all compute, including training, meshes, failed calls, preprocessing and tuning.
6. Establish constraints and acceptance thresholds before running a comparison.
7. Report uncertainty and counterevidence. Search absence is not novelty.
8. Physical experiments require a competent supervised laboratory and a test-specific risk assessment; passing an elastic calculation does not authorize a load test.

## Working assumptions and scope limits
The user has a Windows 11 PC with Intel N100, four cores, approximately 15.75 GiB installed RAM. This is evidence for a small CPU-first experiment, not a measured runtime forecast. Python 3.14.4, NumPy 2.4.6 and SciPy 1.17.1 are available. Candidate geometry and verification packages still require installation/compatibility evaluation.

The synthetic material, dimensions, prescribed forces and ideal clamp in [BASELINE_EXPERIMENT.md](BASELINE_EXPERIMENT.md) are numerical inputs. They do not describe a purchased alloy, process or service environment. Material/load/cost requirements for a real application remain unspecified. A cost model must eventually include stock, setup, cutting/printing, supports, finishing, inspection and rejects; mass alone is not cost.

## Measurable outcomes and gates
| Gate | Evidence required | Failure response |
|---|---|---|
| Numerical credibility | Patch/manufactured tests, force/energy checks, gradient check, mesh study and independently assembled solution | Fix model/solver; no optimization claim |
| Baseline credibility | Fixed-volume references and established SIMP method; same loads, material, domain, pads and compute accounting | Revise comparison, do not claim improvement |
| Useful discrepancy | Frozen candidates change feasibility or meaningful rank after separately measured discretization, reconstruction or physics changes | If below numerical tolerance, drop adaptive audit hypothesis for this family |
| Useful allocation method | Better verified feasible selection or lower total cost at equal selection quality against strongest applicable methods | Report null/negative result; retain simple verification |
| Transfer | Benefit on held-out load/geometry families, with calibration and failures reported | Restrict claims or revise representation |
| Physical validity, later | Independent specimens and held-out measured load–displacement/strain behavior with uncertainty | Revisit joint, material and process assumptions |

Primary metrics: verified feasible candidate rate, verified objective, selection regret relative to an exhaustively audited finite catalogue, number and cost of expensive evaluations, wall time, memory, false-feasibility fraction, invalid-run fraction, and consistency across runs. Secondary later metrics: actual manufacturing yield, cost and simulation-to-measurement error.

Do not set a promised percentage improvement. A development pilot estimates solver noise, timing variance, paired outcome variance and effect sizes that matter in practice; freeze a minimum useful difference and sample size before confirmatory trials. Numerical tolerances in the benchmark are provisional engineering acceptance rules, not empirical improvement targets.

## Completion standard for Phase 0
A bounded, auditable review; defensible baseline; explicitly falsifiable candidate; independent critical review; and a conditional next action. A novel method, a successful bracket or an exhaustive global prior-art search is not claimed. [RESEARCH_STATUS.md](RESEARCH_STATUS.md) records what was actually done and what remains open.

