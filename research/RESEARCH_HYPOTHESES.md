# Candidate research hypotheses and novelty assessment
Date: 2026-10-10. Every hypothesis below is **untested**. Ranking is a research decision, not a result.

## What prior art rules out
Worst-case compliance and stress under uncertain loads are established [R01](https://doi.org/10.1016/j.cma.2017.02.015). Active learning near failure boundaries is established [R02](https://doi.org/10.1016/j.strusafe.2011.01.002). Robust optimization with adaptive multifidelity surrogates exists [R04](https://discovery.ucl.ac.uk/id/eprint/10134527/); nonhierarchical sample value per cost exists [R05](https://arxiv.org/abs/2310.03298); MFMC has already been placed inside robust continuum topology optimization [R06](https://doi.org/10.1016/j.cma.2025.118714). Thus “physics + AI + uncertainty + manufacturing + experiment selection” is not a defensible novelty claim.

[The literature review](LITERATURE_REVIEW.md) also identifies autonomous fabrication/testing and simulation-assisted discrepancy learning. A new acquisition name or an assembly of familiar modules is insufficient.

## H1 — Adaptive load counterexamples
1. **Problem:** reduce worst compliance over an uncertain load envelope at fixed material budget.
2. **Existing approaches:** robust topology optimization, scenario exchange, ellipsoidal worst-case formulations, reliability and fail-safe optimization.
3. **Strongest competitor:** exact reduced load-space optimization wherever linearity permits; robust SIMP/level-set methods otherwise.
4. **Unresolved question:** whether model changes such as contact or uncertain support position make a useful separation oracle expensive.
5. **Mechanism:** alternate a restricted design solve with an adversary finding a violated load/model condition.
6. **Formulation:** minimize t over x subject to C(x,z_i)≤t; add z*=argmax_z C(x,z). An approximate adversary supplies a found counterexample, not proof of its absence.
7. **Assumptions:** specified uncertainty set, reproducible solver and stable design sensitivities.
8. **Expected advantage:** fewer explicit scenarios than an unnecessarily dense grid.
9. **Failure modes:** adversary misses modes; added scenarios merely spend more compute; search learns solver artifacts.
10. **Feasibility:** easy for 2D linear mechanics, but that ease removes the need for AI.
11. **Validation:** compare to the exact load-envelope answer, then independent altered-support/contact cases if justified.
12. **Falsification:** cannot improve cost/quality over the exact oracle or standard exchange at equal budget.
13. **Novelty evidence needed:** full comparison to robust TO and separation-oracle literature under the actual nontrivial uncertainty model.
14. **Value/ranking:** useful baseline, **rank 4; reject as the initial original contribution**.

For fixed stiffness K and two assembled load bases F=[f_x,f_down], solve KU=F once per design. For a=[sinθ,cosθ], C=aᵀ(FᵀU)a. The interval maximum comes from endpoints and stationary angles, not a neural model or thousands of FE solves. Uniform force scaling changes C quadratically; uniform stiffness scaling changes it inversely if every stiffness term scales. These identities cease to cover moving load patches, changing supports or nonlinear response.

## H2 — Cost-aware independent verification of a fixed design catalogue
1. **Problem:** with limited verification budget, select a truly feasible good design from candidates whose coarse density results can misrank their reconstructed geometries.
2. **Existing approaches:** active reliability, multifidelity optimization, surrogate model management, discrepancy calibration and post-optimization validation.
3. **Strongest competitors:** exact cheap physics reductions; coarse-rank audit; space-filling audit; AK-MCS-style U learning; cost-aware multifidelity acquisition using the same descriptors and data. Multi-information-source knowledge gradient [R11](https://papers.nips.cc/paper_files/paper/2017/hash/df1f1d20ee86704251795841e6a9405a-Abstract.html), multifidelity best-arm identification [R12](https://papers.neurips.cc/paper_files/paper/2022/hash/71c31ebf577ffdad5f4a74156daad518-Abstract-Conference.html), and constrained ranking/selection [R13](https://doi.org/10.1016/j.automatica.2018.12.008) are especially close threats. R05/R06 are additional mechanism threats, not automatically drop-in baselines for finite selection.
4. **Unresolved question:** is reconstruction/verification discrepancy large, structured and learnable enough that adaptive allocation helps on a small CPU budget? No source search proves this unanswered universally.
5. **Mechanism:** learn separate response discrepancies for frozen candidate geometries; prioritize independent audits that reduce uncertainty in the final feasible choice; abstain until the chosen design is directly verified.
6. **Formulation:** the explicit finite decision experiment below.
7. **Assumptions:** an immutable catalogue, fixed reconstruction rule, enough cross-candidate correlation, meaningful verification fidelity, and stable call-cost measurement.
8. **Expected advantage:** fewer expensive audits for the same verified selection quality; not a guarantee of global optimum or structural safety.
9. **Failure modes:** discontinuous reconstruction; rare missed failure islands; descriptors omit topology; confidently wrong GP; all candidates invalid; independent verification cheap enough that learning overhead dominates.
10. **Feasibility:** pilot with 12–24 modest 2D candidates and CPU Gaussian processes or a simple regression baseline. No GPU. Scale only after measured timing.
11. **Validation:** exhaustive offline audits as hidden truth, equal-budget replay, paired seeds, unseen geometry/load families, complete failure accounting and component ablations.
12. **Falsification:** negligible baseline discrepancy; no benefit over U-learning/cost-aware methods; apparent benefit disappears when overhead is charged; or calibration/false-feasibility is worse.
13. **Novelty evidence needed:** closest-method reproductions and an explicit mathematical difference that survives the comparison. If the acquisition is merely standard value-of-information per cost, claim only a reproducible engineering benchmark/application.
14. **Value/ranking:** **rank 1, conditional recommendation**. Useful negative results can identify when full verification is cheaper and more trustworthy than AI.

### Proposed H2 experiment, not an implemented algorithm
Let catalogue D={d₁,…,d_N} be fixed before allocation. Each d has density field, immutable reconstruction T(d), material, loads and transforms. Coarse responses y_L and verified responses y_H must refer to explicitly paired quantities. Keep three effects separate:
- **Numerical error:** same material field/geometry, different mesh or element/assembler.
- **Representation discrepancy:** density model versus binary/reconstructed geometry.
- **Model discrepancy:** plane stress versus 3D/joints/contact or physical observations, deferred initially.

Do not call one an error estimate of another. Fine FEA remains a model. Save mass before and after T; a mass-violating reconstruction is infeasible, never silently rescaled.

For each applicable response j, use normalized feasibility margin g_j(d)=response_j/limit_j−1, with feasible iff all g_j≤0. Geometry predicates and solver success are additional requirements. Limits are frozen from application requirements, or declared development-calibrated experimental thresholds, never retrospectively chosen on held-out results. No material strength limit exists in v0; stress is a diagnostic.

**v0 has no application-derived response limit.** Its first audit therefore measures verified equal-mass compliance ranking and geometric validity only. It cannot report stress/displacement false-feasibility or run a meaningful response-boundary classifier until a separate thresholded experiment has frozen a defensible limit. Do not introduce a threshold merely to make an AI method look useful.

If that thresholded extension is justified, a prototype may fit δ_j(d)=g_H,j(d)−g_L,j(d) with a GP on a small descriptor vector: mass, nominal compliance, exact angular worst compliance, grey fraction, minimum measured feature, perimeter/area and connected components. Mass and connectivity can often be checked without FEA; do those first. Penalize invalid geometry explicitly, not through invented numeric response labels. Pure ranking can instead model verified compliance discrepancy, without fabricated feasibility labels.

Define the admissible terminal action set A(data) as **already H-audited feasible candidates plus abstain**. No unevaluated candidate may be recommended. For feasible r use loss L(r,Y)=[C_H(r)−C_H*]_+/s_C, where C_H* is the best feasible catalogue compliance. Abstain has declared loss λ_abstain. If no catalogue candidate is feasible, define a separate no-feasible-candidate state: abstention is then correct rather than inventing a regret value. Scaling and abstention losses are frozen on development cases. Predicted false-feasibility remains a diagnostic on hidden labels, never permission to recommend an unaudited design.

A candidate acquisition to investigate is expected reduction in posterior decision risk divided by a **prospective cost estimate** available before the query:
a(d,l) = {min_(r in A(data)) E[L(r,Y)|data] − E_y[min_(r in A(data+(d,l,y))) E[L(r,Y)|data,(d,l,y)]]}/estimated_cost(d,l).

Use a constant cost per audit type initially, or a frozen estimator trained only on development/previously queried timings and cheap features. Hidden exhaustive-oracle runtimes must not enter query choice. Charge actual measured cost only after querying; an oracle-known-cost policy may appear solely as a labeled unattainable upper bound. This prevents timing leakage in offline replay.
The fantasy terminal set gains d only after a successful H audit; a lower-fidelity query cannot make it admissible. All H calls, including initialization and any final audit, consume the common budget. Approximate using modest posterior samples; log Monte Carlo seed and acquisition cost. **This is a familiar value-of-information construction, not a claimed new formula.** With deterministic audited objectives and a known feasible incumbent, the one-H-query objective can reduce to ordinary constrained expected improvement per cost. That simpler rule is a mandatory comparator. R11 is direct prior art for the acquisition principle. R12 requires fidelity bias information this benchmark does not yet provide; R13 assumes stochastic simulation, while repeated deterministic FE solves do not supply independent noisy observations. Do not import their theorems without meeting assumptions. A simpler mismatch-probability-per-cost rule is an ablation, not an assumed superior method.

First stage has one expensive audit type H and all L results already available. Do not label that stage multifidelity optimization: it is sequential verification. Add a medium mesh only after its cost/correlation benefit is measured. Every final recommendation needs direct H verification; predicted feasibility is never acceptance. H2 does not require a safety certificate from posterior variance.

### Evaluation specification
A 12-candidate development catalogue is a proposed workload cap: 3 volume settings × 4 deterministic generator settings, keeping each volume family compared at its own mass budget. This supplies only four candidates per equal-mass selection task: sufficient to screen for discrepancy, insufficient to substantiate a complex GP allocator. Freeze the actual generator settings only after v0 converges. Expand the within-volume catalogue for method comparison only when the measured cost and outcome variation justify it, and preregister its size. Do not mix objectives across volume families or count intermediate optimization iterates as independent specimens. Reserve separate geometry/load families for confirmation.

The small pilot fits no confirmatory cross-volume model. A later study must freeze whether training pools volume families; the default is separate within-volume selection models. If pooling is tested, include volume explicitly, normalize responses with a declared reference, give every comparator the same training data, and keep held-out family labels inaccessible. Pooling is an ablation, not an invisible advantage.

Exhaustively audit development and later held-out catalogues offline; hide unqueried H outputs **and timings** from allocation methods. Charge queried H costs in replay and disclose total offline oracle-building cost separately. Reusing an oracle for replay does not turn those measurements into live speedups. The catalogue optimum is a full-budget evaluation reference, not a clairvoyant algorithm available at partial budgets. A sequential exhaustive policy uses a frozen candidate-ID order and the same audited-incumbent/abstention rule; it reaches the optimum only after its actual audit costs are paid.

Compare: exhaustive audit; ascending coarse objective with hard cheap screens; seeded random/space-filling; simple expected verified improvement per cost (constrained EI only after limits exist); U-learning applied to paired margins only in the thresholded extension (label adaptation from reliability estimation); cost-normalized uncertainty/mismatch learning; and an R11-inspired finite-selection knowledge-gradient/decision-risk comparator with identical descriptors, observations and terminal H-audit rule. R05-like multifidelity selection is additional only if multiple fidelity query types are actually introduced. Document any comparator that cannot be reproduced; do not claim state-of-the-art superiority when a strong comparator is missing.

For stochastic selectors use seeds 0–9 as a development design; determine confirmatory sample count from a paired variance/power pilot and freeze it. Deterministic SIMP has no meaningful ten-seed result unless its initialization actually differs. Pair catalogues and seeds across methods. Budget checkpoints are equal measured CPU/wall costs, with solver-evaluation counts also reported. Cap to the measured PC envelope; prefer sequential execution to competing processes.

Report feasible final recommendation/abstention rates, selected verified compliance, catalogue regret, false-feasibility among predicted positives (denominator explicit), failures, calibration and total time. Bootstrap/paired intervals across independent benchmark families; repeated seeds on one catalogue do not establish geometry generalization. Compare numerical/engineering relevance as well as intervals. If no candidate is feasible, report that fact without manufacturing a regret score.

## H3 — Manufacturing-tolerance and support-model robustness
1. **Problem:** designs optimized with ideal geometry/clamps can degrade under tolerances and realistic attachment models.
2. **Existing:** eroded/intermediate/dilated TO, machining accessibility, stress-constrained manufacturing tolerance, fail-safe methods.
3. **Competitors:** established robust projection and process-constrained optimization with independent checking.
4. **Unresolved:** which joint/tolerance parameters dominate this particular manufactured application; no application data yet.
5. **Mechanism:** coupled uncertain offsets and support stiffness, then nested robust optimization.
6. **Formulation:** min_x max_(δ,k_s) C(T_δ(x),k_s), with volume and explicit process predicates across the same uncertainty set.
7. **Assumptions:** calibrated offset/support ranges, stable geometry, no unmodeled failure mode dominating.
8. **Advantage:** more transferable performance than nominal optimization.
9. **Failures:** arbitrary uncertainty bounds; expensive nested mesh changes; process geometry cannot be represented by uniform offset.
10. **Compute:** small 2D sweeps possible; full 3D process optimization premature.
11. **Validation:** matched nominal/robust methods, measured parts and fixtures, held-out tolerances.
12. **Reject if:** standard robust projection already matches performance or estimated uncertainty is smaller than discretization error.
13. **Novelty:** must identify a nontrivial coupling missing in nearest process/joint-aware methods.
14. **Value/ranking:** practical engineering value, **rank 2**, application definition first.

## H4 — Representation-independent knowledge transfer
1. **Problem:** reuse expensive data across density, truss and CAD families.
2. **Existing:** mesh GNNs, neural operators, learned latent inverse design, differentiable geometry.
3. **Competitors:** from-scratch adjoint optimization, simple normalized descriptors, nearest-neighbor/warm-start transfer.
4. **Unresolved:** whether invariant features transfer constraint knowledge across genuinely different domains.
5. **Mechanism:** encode geometry/load/support graphs and predict correction or warm start.
6. **Formulation:** train f_θ(G,BC,material) to predict response residual; minimize held-out response/constraint loss, then independently solve candidates.
7. **Assumptions:** representative cross-family labelled data and faithful encoding of BCs/materials.
8. **Advantage:** amortized evaluation across many tasks.
9. **Failures:** data leakage; mesh dependence; training cost exceeds savings; out-of-distribution confident errors.
10. **Compute:** dataset-generation cost presently unjustified for N100; defer.
11. **Validation:** entire geometry families held out; count offline generation and training.
12. **Reject if:** simple warm starts match total-cost performance or feasibility degrades.
13. **Novelty:** explicit comparison to GNN/operator/inverse-design prior art under the same transfer definition.
14. **Value/ranking:** potentially broad, **rank 5** due to poor initial feasibility.

## H5 — Physical experiment choice for discrepancy identification
1. **Problem:** choose tests that separate material, joint and model discrepancy before redesign.
2. **Existing:** Bayesian calibration [R08](https://doi.org/10.1111/1467-9868.00294), discrepancy-aware experimental design [R09](https://doi.org/10.1016/j.jcp.2025.114469), autonomous mechanical experiments in the literature review.
3. **Competitors:** fixed coupon/fixture test plan and conventional information-gain design.
4. **Unresolved:** identifiability and practical value for a specified real bracket.
5. **Mechanism:** choose load/measurement to reduce uncertainty relevant to the next design decision.
6. **Formulation:** argmax_e expected reduction in decision loss divided by test cost, given y=h(x,θ)+δ(x)+ε.
7. **Assumptions:** safe repeatable apparatus, calibrated sensors, identifiable model, independent specimens.
8. **Advantage:** fewer physical tests to make a reliable redesign decision.
9. **Failures:** δ absorbs parameter error; sensor bias; unmodeled plasticity; adaptive tests overfit same specimen.
10. **Compute/resources:** feasible inference at small dimension, laboratory access is the binding unknown.
11. **Validation:** randomized fixed-plan comparator, test count/time, withheld specimens and loads.
12. **Reject if:** no benefit over planned calibration tests or parameter/discrepancy remains unidentifiable.
13. **Novelty:** distinctive mechanism and replication against prior autonomous mechanics and BED work.
14. **Value/ranking:** long-term practical value, **rank 3**, defer until numerical and laboratory gates pass.

## Patent and novelty boundary
The Autodesk disclosure [R10](https://patents.google.com/patent/US11468634B2/en) already describes learned fine-design updates driven by coarse structural analyses. Additional Chinese load-uncertainty and machining-constraint patents surfaced, but their complete claims/families were not reviewed. This is technical prior-art screening, not legal advice or freedom-to-operate clearance.

Before a publication novelty claim: read closest papers in full, trace backward/forward citations, search alternative terminology (model management, multifidelity best-arm identification, active reliability, verification-driven design, reconstruction uncertainty, decision-focused experimental design), reproduce strong baselines, and have a domain expert challenge the exact proposed delta. If nothing survives, deliver the benchmark and negative finding; do not rename existing algorithms.
