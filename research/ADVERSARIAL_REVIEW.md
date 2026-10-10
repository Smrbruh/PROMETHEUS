# Independent adversarial review

Review date: 2026-10-10. Reviewer: a separately delegated mechanics/review agent, after the primary proposal was drafted. This is an independent reasoning pass within the same research team, not an external peer review or experimental replication.

Scope: [RESEARCH_HYPOTHESES.md](RESEARCH_HYPOTHESES.md), [ARCHITECTURE.md](ARCHITECTURE.md), the three mechanics/baseline/verification documents, and `exploratory/linear_load_envelope.py`. The review checked mathematical consistency, scientific distinguishability, comparison fairness, PC feasibility and missing gates. It did not run finite-element analysis, prove novelty, inspect every cited full text or validate a physical bracket.

## Review decision

**NARROW is justified.** The strongest present contribution is a proposed reproducible discrepancy audit, not an original optimization algorithm. Authorize the smallest checked numerical experiment: solve the declared plate, reproduce classical optimization, convert geometries by a frozen rule, and independently measure what changes. The evidence does not justify a general invention engine, neural surrogate training, or physical claims.

The proposal becomes scientifically useful if it can return a negative answer: reconstruction discrepancy may be small; simple coarse ranking may be adequate; independent verification may be cheaper than acquiring/training a model. The updated proposal now allows these outcomes instead of assuming that adaptive allocation must win.

## Findings requiring action or explicit limits

| Finding | Severity | Response and remaining gate |
|---|---|---|
| H2 uses familiar value-of-information/knowledge-gradient ideas | Major novelty risk | The hypothesis now states that the acquisition formula is conventional. A new algorithm claim requires a precise difference and closest-method reproduction; neither exists yet. A benchmark/application result remains possible. |
| Angle-only linear load adversaries are unnecessary | Major if ignored; resolved in scope | Two basis solves and a small analytic maximization resolve the whole interval. The exact method is mandatory. H1 was downgraded rather than marketed as novel AI. |
| v0 has no displacement/strength feasibility limit | Major protocol ambiguity; corrected | The initial audit now studies equal-mass compliance ranking and geometric validity. Response-based false-feasibility and U-learning are deferred until a defensible limit is frozen. No threshold should be invented to create a favourable classifier experiment. |
| Acquisition originally admitted unverified terminal recommendations | Major policy mismatch; corrected | Terminal actions now consist of H-audited feasible candidates plus abstain. A fantasy H query can enlarge that set; a cheaper query cannot. Initialization and any final verification consume the shared budget. |
| One-step selection can collapse to expected improvement per cost | Major novelty/comparator issue; acknowledged | With deterministic audited objectives and a feasible incumbent, the risk formulation can have no advantage over a standard EI-style policy. The simpler rule is now an explicit comparator. More elaborate posterior sampling is not an inherent contribution. |
| Twelve pilot candidates mean only four choices at each of three mass budgets | Major statistical limit; acknowledged | Appropriate for detecting discrepancy, insufficient to substantiate a complex GP policy. Enlarge the within-volume catalog only after measured variation/cost justify it; freeze the new catalog and confirmatory protocol first. |
| Deterministic FEA does not provide stochastic bandit samples | Major theorem-transfer risk; acknowledged | Repeated solves of an identical deterministic problem cannot create independent observation noise. Borrowed stochastic ranking/fidelity guarantees need their actual assumptions, not just analogous terminology. |
| Artificial void stiffness can hide an unusable topology | Major engineering risk; addressed by gate | The plan requires edge connectivity, a frozen-design Emin study and a material-only binary verifier. A numerically finite relaxed solve alone is not a valid part. No failures may be repaired without a new candidate record. |
| Optimized design changes during a refinement study | Major potential confounding; avoided | Freeze densities/geometry for solution convergence. Reoptimization at a fixed physical filter radius is a separate study. Refitering or thresholding again would conflate representation and numerical errors. |
| Peak stress at the L corner is not a convergence target | Major strength-claim risk; bounded | No stress safety claim is permitted. Physical-distance stress exclusions are diagnostics only, and all unrestricted peaks/trends remain visible. A later strength benchmark requires smooth geometry, material data and a different verification task. |
| PC resources may make full enumeration preferable | Major feasibility question; open | The N100/16 GiB environment supports a small prototype in principle, not an assumed large research campaign. Measure actual factorization fill, time and memory before expanding. The 8 GiB process cap is policy, not a demonstrated bound. |
| No verifier, optimization or convergence study has yet run | Blocking for results, not Phase 0 | Implementation is a future milestone. No claims of feasible optimized brackets, measured speedups or manufacturing success follow from the current algebra-only utility. |

## Mathematical audit

The declared unit system is consistent: E in MPa is N/mm²; the 5 mm thickness multiplies stiffness integrals; the boundary line load is force divided by 20 mm and must not acquire a second thickness factor. The total area is 6,400 mm², preserved area 600 mm², 40% volume budget 12,800 mm³, and mass at the synthetic density is 34.56 g. At h=2 mm, there are 1,600 active cells, including 150 fixed-solid cells and a 490-cell variable budget for the binary comparator.

The lower raw-density bound and the ersatz modulus are separate. The volume uses physical densities, includes pads and excludes the geometric cutout. The pad-aware filter makes a naïve raw initialization at 0.337931 incorrect; the specified scalar bisection on actual volume fixes this. The compliance and volume gradients correctly apply the transpose filter and exclude overridden pad rows. The OC multiplier search evaluates the same physical volume as the constraint.

The nominal objective is a standard nonconvex compliance problem. Compliance is a load-weighted displacement/energy quantity, not a guarantee about every nodal displacement or failure mode. Full-solid compliance is a useful monotonicity check for the same load and softer density field; it is not an equal-mass competitor.

The exact angle envelope formula was checked. The hypothesis/utility uses an x/right and y/**down** load basis, so its coefficient vector is `[sin θ, cos θ]` and cross term is positive. The mathematical document uses a y/**up** basis, so its vector is `[sin θ, −cos θ]` and cross term is negative. These conventions are equivalent after changing the sign of the second basis column and corresponding off-diagonal Q entries. There is no disagreement in physics. The derivative in the utility is `(Q00−Q11)sin(2θ)+2Q01 cos(2θ)`; its `atan2` construction enumerates stationary angles correctly.

The specified sinusoidal manufactured-solution body force was independently differentiated from `−div σ(u)`. The mixed derivatives and signs are consistent with the plane-stress Lamé parameter. Full boundary zero displacement, thickness weighting, physical units, three test meshes and quadrature sensitivity are now explicit. The fixture remains unexecuted; specifying it is not passing it.

## Utility inspection and additional checks

The existing utility's saved result reports 40 synthetic positive-definite systems, with maximum direct-versus-two-RHS relative error about 8.52×10⁻¹⁶ and maximum analytic-versus-dense-grid relative gap about 8.59×10⁻¹⁰. Its scope field explicitly excludes bracket FEA and optimization.

The reviewer separately imported its `envelope` function and compared four deterministic Q matrices with a 100,001-angle grid over [−20°,20°]:

| Case | Q | Analytic maximum | Maximizing angle |
|---|---|---:|---:|
| Constant response | identity | 1 | −20° returned, all angles equivalent |
| Endpoint maximum | diag(2,1) | 1.116977778441 | −20° returned, both endpoints equivalent |
| Interior maximum | diag(1,2) | 2 | 0° |
| Nonzero cross term | [[2,0.3],[0.3,3]] | 3.083095189485 | 15.48187827° |

All four comparisons passed the test `analytic ≥ grid maximum −10⁻¹²` and absolute gap <10⁻⁹. These tests add constant/endpoint/interior coverage; they do not establish solver correctness for mechanics. The function assumes a symmetric finite 2×2 Q and an ordered angle interval. Those assumptions hold in this utility's construction; a reusable production API would need input checks rather than silently accepting arbitrary matrices.

## Fairness and falsification audit

The same immutable geometry transform, pads, boundary mapping, volume accounting, mesh acceptance and final-verification rules must apply to all selectors. Full solid must not be ranked as an equal-mass design. A lower-density solution must not be compared to a binary part without disclosing their different constitutive and geometry definitions.

The offline exhaustive catalog supplies labels for replay; constructing it is a real cost and must be reported separately. A replay's virtual query count is not a measured live speedup. Charge surrogate training, acquisition computation, failures and terminal audits. On very small catalogs, the absolute seconds saved and the break-even number of future uses matter more than a percentage of solver calls.

The updated common development seed set is 0–9. Multiple seeds on one catalog measure policy variability, not geometry generalization. Near-duplicate iterates may not cross train/test boundaries. Confirmation needs independent geometry or parameter families, a preregistered budget, an effect size meaningful relative to numerical error, and enough independent cases; a ten-seed pilot alone cannot provide that.

Reject or revise H2 if any of these occurs:

1. Coarse-versus-verified ranking differences are within numerical uncertainty or absent after cheap screening.
2. A simple deterministic discrepancy indicator or coarse-rank audit achieves the same decision quality at lower total cost.
3. Cost-aware learning adds overhead without a useful verification saving, or any gain disappears against EI/knowledge-gradient/appropriate active-learning baselines.
4. A promising result depends on leaked hidden labels, favourable post-hoc thresholds, excluding failed cases or granting the proposed method more computation.
5. Predictive confidence is treated as permission to skip the required final H audit.
6. The chosen geometry is materially sensitive to arbitrary reconstruction, filter or attachment assumptions without being evaluated across them.

## Additional independent source and protocol reviews

Two other delegated researchers spot-checked close-prior-art bibliographic records and challenged the allocation protocol. No mismatch was found in the selected primary-source checks (R05–R07, R09, R12, L17 and L28); this is not an audit of every bibliography entry. Their substantive corrections were incorporated in the hypothesis document:

- Query selection uses only a prospective cost estimate from allowed data. Actual hidden-oracle runtimes are charged after a query, never leaked into its selection.
- Exhaustive audit is a full-budget reference. Its sequential version has an explicit fixed order and abstention rule; the catalogue optimum is not an available partial-budget algorithm.
- Four candidates per mass budget are a discrepancy pilot only. Within-volume model training is the default for an expanded study; any pooling policy must be declared and matched across methods.
- R11-style knowledge gradient and simple expected verified improvement replace a less directly matched multifidelity comparator as mandatory first-stage comparisons.

These corrections improve falsifiability; they do not supply experimental evidence of an advantage.

## What remains uncertain

No current result establishes convergence within 300 designs, h=0.5 mm solver cost, the stability of sharp pixel geometry conversion, or enough learnable discrepancy to motivate selection. The initial 1%/2% mesh-comparison gates and stopping settings are provisional; changes must be justified on development tests and applied to every comparator before confirmation. These limitations are reasons to stage the experiment, not to claim a benefit or to build the product first.

The next action with highest information value is a verified B0/B1/B3 run and reconstruction discrepancy audit. If it shows no actionable uncertainty or no cost advantage for selection, retain a transparent classical design/verification workflow and stop the AI allocation branch. If it reveals a structured discrepancy and a meaningful cost gap, design the larger H2 comparison before implementing its advanced selector.
