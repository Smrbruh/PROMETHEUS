# Decision memo — PROMETHEUS Phase 0
**Decision: NARROW.** Date: 2026-10-10.

Investigate a reproducible **independent-verification and discrepancy benchmark for lightweight structures** before building an invention system. The evidence supports a practical, falsifiable research experiment. It does **not** yet establish an original optimization algorithm, a superior bracket, or a validated physical design.

## Why this decision
Classical topology optimization is mature. Robust load/stress optimization already exists; simulation-informed autonomous mechanical experimentation already exists. Even the narrower idea of selecting costly evaluations by expected information value has direct predecessors in multi-information-source optimization and constrained ranking/selection. [Robust TO](https://doi.org/10.1016/j.cma.2017.02.015), [simulation-assisted autonomous experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC8010472/), [information-source optimization](https://papers.nips.cc/paper_files/paper/2017/hash/df1f1d20ee86704251795841e6a9405a-Abstract.html).

The relevant question is narrower: **when coarse density optimization produces several plausible designs, do geometry reconstruction and independent verification change the engineering decision enough that selective auditing saves total computation?** That is testable and useful even if the answer is no. A new algorithm is warranted only after a measured gap and a defensible prior-art difference emerge.

## What to build first
A small CPU experiment, using Python/NumPy/SciPy, that reproduces established density-filtered SIMP/OC and independently checks the solution with a different element/assembly path. The initial reference is a 120×80 mm L-shaped plane-stress plate, 5 mm thick, with explicit preserved pads, a distributed in-plane load and a 40% total-volume budget. It is a numerical benchmark with synthetic material parameters, not a certified metal bracket.

Compare a solid reference, an exact-volume simple L-path, a grey-density diagnostic and the established optimizer. First verify identical geometry/material fields across meshes and solvers; then measure the effect of a declared binary reconstruction. Do not hide that those are different comparisons.

## Strongest candidate and required competitors
**H2: cost-aware allocation of independent verification**, conditional on the discrepancy audit. Its value-of-information principle is established. The first credible contribution may be a reusable benchmark and negative/positive empirical findings, not a new mathematical method.

Required competitors: direct exhaustive verification with honest cost accounting; coarse-objective ordering; random/space-filling ordering; expected verified improvement per cost; and the closest finite-selection/knowledge-gradient method. Add active reliability/U-learning only after an explicit response-feasibility threshold exists. For fixed linear load-direction uncertainty, use the exact two-RHS oracle; AI cannot earn a speedup by competing against redundant FE solves.

All final recommendations must have paid for a successful independent verification. Hidden audit labels and runtimes cannot guide query choices.

## What evidence justifies proceeding
1. The baseline passes element, derivative, equilibrium, energy, fixed-design mesh and independent-assembly checks.
2. Reconstruction/model changes produce repeatable, decision-relevant discrepancy exceeding numerical uncertainty.
3. There is a real cost tradeoff on this PC; exhaustive verification is not already cheaper than fitting/selecting.
4. A frozen method improves verified selection quality per **total** budget against strong comparators on held-out families.
5. The exact contribution survives a full-text prior-art review and independent domain critique.

No percentage speedup is promised. Pilot data must set effect sizes, confirmatory sample count and budgets before test results are opened.

## Change direction or stop when
- The solver/gradient checks fail: repair the mechanics, do not optimize around the error.
- Discrepancy is negligible: use conventional optimization and routine verification.
- Geometry conversion dominates: improve geometry/constraints before adding learning.
- Standard allocation methods match the proposed one: adopt them; stop claiming algorithmic novelty.
- Benefit vanishes after charging meshing, training and final verification: reject the efficiency hypothesis.
- A real use case lacks material, load, attachment or manufacturing definitions: remain a numerical benchmark.
- Physical behavior disagrees: revise the model rather than presenting optimized simulation as safety evidence.

## Highest-information next action
**Run one established SIMP baseline with a paired independent verification and reconstruction audit.** This includes prerequisite numerical checks. It can determine whether the proposed research problem exists before any AI model or application is built.

## What Phase 0 actually delivered
Twelve requested research documents, a deduplicated bibliography and source/search records, an independent adversarial review, and a small algebra-only utility. The utility checks the two-RHS load envelope and mass bookkeeping. **No bracket FE analysis, optimization campaign, learned model, commercial benchmark, manufacturing run or physical test was performed.** Access limitations and unresolved questions are listed in [RESEARCH_STATUS.md](RESEARCH_STATUS.md).

