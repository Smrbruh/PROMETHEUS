# Phase 1 numerical verification and validation record

Research date: 2026-10-10. **Code checks support the implemented linear plane-stress equations; physical validation has not been performed.** For the canonical 40% setting, B0, B1 and the frozen B3 density field pass the declared nominal compliance and pad-displacement resolution gates. Both binary B3 reconstructions remain `resolution_unverified`. Passing algebraic checks alone does not change that result.

This report uses the measured [analytical/manufactured evidence](../results/validation.json), [canonical manifest](../results/phase1_v040_verified/manifest.json), [canonical summary](../results/phase1_v040_verified/summary.json), and [sensitivity evidence](../results/phase1_v040_verified/sensitivity.json). The final automated-test result and count belong to [tests_final.xml](../results/tests_final.xml); they are not inferred from the presence of test code. Results for the other volume settings and optimization-grid study are reported separately in [BASELINE_RESULTS.md](BASELINE_RESULTS.md).

## 1. Definitions and acceptance scope

Final automated command: `python -m pytest -q --junitxml=results/tests_final.xml`. **104 tests passed in 18.87 s**, with no failures or skips. This is separate from the timed benchmark bundles.

The separate [artifact integrity audit](../results/evidence_audit.json) passes for all four completed run bundles: source/specification snapshot hashes, raw/physical-array hashes, histories, coefficient volumes, saved status consistency and solve counts agree. The checker recalculates acceptance from archived metrics; it does not rerun mechanics or create a stronger reference label.

Code verification tests whether the discretized equations and mathematical components are implemented correctly. Solution verification assesses resolution for a particular frozen model. Neither establishes agreement with a manufactured specimen, actual material, attachment, loading fixture or failure mechanism. The governing fixture and thresholds are in [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md); implementation choices are in [IMPLEMENTATION_REPORT.md](IMPLEMENTATION_REPORT.md).

The `independent_verified` flag concerns **nominal downward-load compliance and nominal maximum pad-displacement magnitude**, for the recorded geometry, physical-density/coefficient field, material and boundary conditions. It is not a statement about all stresses, every displacement, buckling, manufacturing feasibility or load-bearing safety. The angular envelope is computed exactly for each frozen discrete linear model, but its maximum has no separate mesh-resolution acceptance gate in this release. It must not inherit a stronger interpretation from the nominal flag.

The full solid B0 can pass numerical checks while violating the 40% mass constraint. A density field such as B3 is a SIMP relaxation; it is not a verified binary part. `optimization_feasible`, `independent_verified`, binary connectivity and physical validation are distinct predicates.

## 2. Element, assembly and analytical references

The primary path uses full-integration Q4 elements. The second path uses straight-sided six-node quadratic triangles, splitting every parent quadrilateral along corner 0–2 and keeping its modulus on both children. Its engineering-strain material law is independently constructed from the plane-stress parameters lambda = E nu/(1-nu²), mu = E/[2(1+nu)]. The three-point triangle rule integrates the quadratic stiffness integrand exactly. Boundary traction is reconstructed from exterior edges; P2 endpoint/midpoint/endpoint weights are 1/6, 2/3, 1/6. Thickness is not applied twice to line loads.

The tests cover shape interpolation and partition of unity; stiffness symmetry; three rigid modes and the expected five Q4 or nine P2 positive deformation modes; explicit element-to-global scattering; affine displacement recovery on irregular conforming patches; essential constraints and nonzero prescribed-displacement work; reactions, energy, superposition and force/modulus scaling; missing boundary patches; free bodies, disconnected islands, invalid geometry and a point-connected P2 hinge. The binary geometry tests separately reject connections sharing only one vertex. No weak springs or artificial supports make these fixtures solvable.

For the analytical tensile strip, the rectangle is 100 by 10 mm, thickness 5 mm, with 1 MPa right-edge traction and resultant 50 N. Only ux is fixed on the left edge; uy is fixed at the origin. Thus lateral Poisson contraction remains free. The exact field is ux=x/E, uy=-nu y/E, the stress is (1,0,0) MPa, and compliance is 50×100/70,000 = **0.0714285714285714 N·mm**. The archived comparison uses h=2 mm in both element families.

| Measured quantity | Q4 | P2 |
|---|---:|---:|
| Compliance [N·mm] | 0.0714285714285732 | 0.0714285714285408 |
| Relative compliance error | 2.53e-14 | 4.29e-13 |
| Relative nodal displacement error | 5.10e-12 | 4.30e-11 |
| Maximum absolute stress-component error [MPa] | 9.61e-13 | 8.58e-12 |
| Reaction x resultant [N] | -50.0000000000024 | -49.9999999999680 |

Both satisfy the declared relative 1e-8 analytical gates. These are exact solutions of the same continuum equations and constraints, not experimental observations. See [test_fem.py](../tests/test_fem.py) and [test_verification.py](../tests/test_verification.py).

## 3. Smooth manufactured solution and quadrature stability

The implemented fixture is the prescribed 100 by 40 mm rectangle with ux=0.01 sin(pi x/100) sin(pi y/40) and uy=0.02 sin(pi x/100) sin(pi y/40), in mm. Exact zero displacement is imposed on every exterior boundary node. Body force is -div(sigma(exact u)), including the mixed derivatives, in N/mm³. Its weak-form integral includes the 5 mm thickness.

The displacement error is sqrt(integral ||u_h-u||² dA), with units mm². The strain-energy error is sqrt(integral t epsilon_errorᵀ D epsilon_error dA), with units sqrt(N·mm). The table uses order-six force/error quadrature; Q4 stiffness remains 2×2 and P2 stiffness remains its exact degree-two rule.

| h [mm] | Q4 displacement L² error | Q4 energy error | P2 displacement L² error | P2 energy error |
|---:|---:|---:|---:|---:|
| 5 | 9.44324196e-3 | 3.65944327 | 3.24092238e-4 | 0.251359222 |
| 2.5 | 2.36734662e-3 | 1.83271525 | 4.06377370e-5 | 0.0630902331 |
| 1.25 | 5.92250258e-4 | 0.916733710 | 5.08356528e-6 | 0.0157887900 |

Measured log2(error(h)/error(h/2)) orders are:

| Norm and element | 5→2.5 mm | 2.5→1.25 mm | Expected smooth-solution order |
|---|---:|---:|---:|
| Q4 displacement L² | 1.99601 | 1.99899 | 2 |
| Q4 energy | 0.997642 | 0.999408 | 1 |
| P2 displacement L² | 2.99551 | 2.99891 | 3 |
| P2 energy | 1.99426 | 1.99852 | 2 |

Both families follow the expected trends. Increasing body-force and error quadrature from order four to six changes the Q4 L² error by at most 1.26e-8 relative and the P2 L² error by at most 3.33e-5 relative; energy differences are smaller. The archived quadrature gates pass. P2 uses a tensor Gauss rule through a Duffy map, so its order labels do not denote the same integration-point arrangement as Q4.

These regularity-based rates are **not imposed on the sharp L corner or the corners of pixel geometry**. No Richardson extrapolation or certified error bound is claimed for those bracket results. The evidence utility records 14 factorizations and 14 right-hand-side solves for its analytical/MMS calculations; these are separate from benchmark costs.

## 4. Density-filter and compliance sensitivities

The filter tests independently calculate boundary, pad-adjacent and cutout rows of the area-weighted normalized density filter. They check exact solid-pad override, initialization at the actual filtered volume, the nonconstant volume gradient, and the transpose chain rule. Derivatives of overridden physical-pad rows are zero. OC checks cover its feasible bisection endpoint, move limits, invalid sensitivities, and the rule that returned candidates must already have been evaluated.

The canonical audit stores complete interior raw/physical arrays, seed 20261010, directions, step sizes, analytical directional derivatives, finite differences and absolute/relative residuals in [sensitivity_well_scaled.npz](../results/phase1_v040_verified/sensitivity_well_scaled.npz) and [sensitivity_signed.npz](../results/phase1_v040_verified/sensitivity_signed.npz). No perturbation is clipped. Steps are 1e-2, 1e-3, 1e-4, 1e-5 and 1e-6. Relative error is |dFD-g·v|/max(|dFD|,|g·v|,1e-10), with absolute residuals retained.

The first sweep uses ten unit directions with nonnegative free entries, including three localized near the support, load pad and re-entrant corner. This avoids cancellation of the negative compliance sensitivities and exposes a truncation-error region. The second sweep reuses those three local directions and adds seven global sign-changing directions; it is not twenty distinct directions.

| Audit measurement | Nonnegative sweep | Signed sweep |
|---|---:|---:|
| Directions passing relative compliance error <1e-4 at two adjacent steps | 10/10 | 10/10 |
| Largest relative compliance error at step 1e-2 | 5.31e-6 | 5.31e-6 |
| Largest relative compliance error at step 1e-3 | 8.00e-8 | 5.28e-6 |
| Largest relative volume error over steps 1e-2…1e-4 | 1.32e-10 | 3.45e-8 |
| Charged factorizations per sweep | 101 | 101 |

Every direction in the nonnegative sweep improves from step 1e-2 to 1e-3, establishing the required central-difference convergence region. The compliance tolerance remains 1e-4 at two adjacent steps, and the volume tolerance remains 1e-7 in the stable step region.

Smaller steps are not uniformly better. The largest signed compliance relative errors at 1e-5 and 1e-6 are 7.46e-4 and 5.22e-4; corresponding maximum volume relative errors are 8.38e-7 and 8.87e-6. These values are archived, not discarded. Some signed directional derivatives nearly cancel; subtraction of nearly equal objective/volume values then amplifies finite-precision error. The linear volume formula has no central-difference truncation error to remove by shrinking the step. Larger-step agreement and the nonnegative sweep support the derivatives; the signed sweep supplies additional accuracy evidence, not a claim of a monotonic truncation trend below roundoff. See [test_optimization.py](../tests/test_optimization.py).

## 5. Load superposition and exact envelope

The original exploratory mathematics is retained in [envelope.py](../prometheus/envelope.py). Basis loads are RIGHT and DOWN, each with the full 1,000 N resultant, and coefficients are [sin(theta),cos(theta)]. Endpoints and every stationary angle in [-20°,20°] are checked. Tests include constant, diagonal and coupled analytic matrices; degenerate intervals; unrestricted eigenvalue limits; invalid inputs; forty deterministic positive-definite systems compared with direct solves and dense angular samples; and a direct 13° finite-element load solve against two-RHS superposition.

The algebraic envelope interface still requires symmetry to 1e-12. Mechanics-generated Gram matrices have separately measured reciprocity error ||Q-Qᵀ||/||Q||. A mechanics adapter rejects this error above 1e-8, then passes (Q+Qᵀ)/2 to the algebraic routine while preserving raw Q. For any real coefficient vector a, aᵀQa equals aᵀ(Q+Qᵀ)a/2 because the skew part contributes zero. This is an explicit finite-precision interface choice, not stiffness averaging or a change to the load model. Canonical raw reciprocity error reaches 1.47e-10, which explains why a discrete mechanical solve can exceed the algebra utility's stricter input tolerance while satisfying the mechanical checks.

Exact angular maximization removes an angular sampling error for the frozen discrete model. It does not remove finite-element error, validate the continuum assumptions, or justify an adaptive selector that competes against repeated solves of this two-dimensional load space.

## 6. Per-solve and fixed-design gates

The bracket uses finite input/solution checks, positive Jacobians and material parameters, preserved boundary extents and pads, and no unresolved rigid motions. Sparse factorization failures propagate. An LU minimum/maximum pivot ratio below 1e-14 is rejected as singular or an unresolved near mechanism; this ordering-dependent ratio is **not a condition-number estimate or forward-error bound**.

| Check | Acceptance |
|---|---:|
| Free-DOF equilibrium residual relative to the loaded system | <=1e-8 |
| Assembled stiffness relative asymmetry, before any averaging | <=1e-12 |
| Applied-plus-reaction force and moment balance | <=1e-7 |
| Work versus twice strain energy | <=1e-8 |
| Q4 separate element-energy sum versus assembled energy | <=1e-8 |
| Physical volume excess | <=1e-6 times full volume |
| Final h=1→0.5 compliance change, both solvers | <1% |
| Final h=1→0.5 maximum pad-displacement change, both solvers | <2% |
| Finest Q4/P2 compliance / pad agreement | <1% / <2% |

Nonzero prescribed displacement contributes reaction work to the energy identity. The P2 energy check uses its assembled stiffness, whereas Q4 also resums element energies. For canonical candidate audits, the largest free residuals are 1.01e-9 for Q4 and 3.65e-9 for P2; largest work-energy errors are 9.08e-11 and 2.20e-9. The largest Q4 element-sum discrepancy is 2.38e-9. Thus all completed canonical candidate solves satisfy their algebraic gates, including the binary reconstructions that fail resolution gates.

Full acceptance also requires the complete h=2,1,0.5 mm sequence with matching candidate/representation, physical-field hash, source geometry, material/problem signature and source spacing. Child elements inherit frozen parent coefficients; binary refinement subdivides the same retained squares. No refiltering or reranking occurs. Relative solver gaps are recomputed from metrics, not trusted from cached display fields. Compliance refinement changes must improve, and compliance must be nondecreasing within a relative 1e-8 solve-roundoff allowance. The latter follows from the nested conforming spaces and unchanged load functional. Pointwise pad maxima have no corresponding variational monotonicity requirement.

## 7. Canonical 40% solution-verification results

Every percentage below is a relative difference, not a certified uncertainty interval. Full raw values and both refinement steps are in the [candidate evaluations](../results/phase1_v040_verified/summary.json).

| Candidate | Q4 compliance change 2→1 / 1→0.5 [%] | P2 compliance change 2→1 / 1→0.5 [%] | Finest compliance / pad gap [%] | Result |
|---|---:|---:|---:|---|
| B0 solid | 0.5220 / 0.2260 | 0.2182 / 0.1026 | 0.0969 / 0.0920 | Numerical gates pass; over mass budget |
| B1 binary L path | 1.8444 / 0.7018 | 0.4785 / 0.2252 | 0.3082 / 0.3097 | Numerical gates pass |
| B3 density | 2.5777 / 0.8471 | 0.4370 / 0.1782 | 0.3289 / 0.3148 | Numerical gates pass for density model |
| B3_binary | 5.9987 / 2.9017 | 2.5137 / 1.1573 | 1.4274 / 1.3296 | Resolution unverified |
| B3_threshold | 6.1523 / 2.9752 | 2.5713 / 1.1826 | 1.4637 / 1.3693 | Resolution unverified |

B2 has only its declared coarse diagnostic solve and remains `resolution_unverified`; no missing refinement is inferred. Both binary reconstructions are connected and pass their algebraic solves, but each fails the compliance refinement and independent-compliance thresholds. Their Q4 final pad changes also exceed 2%: 2.7080% and 2.7871%. A smaller finest pad gap does not override failed compliance gates. No binary result is promoted by loosening a threshold.

At h=0.5 mm, the Q4/P2 nominal compliance pairs [N·mm] are B0 449.893310/450.329736, B1 12987.600594/13027.747319, B3 density 1865.950061/1872.107811, B3_binary 1230.659572/1248.480615, and B3_threshold 1249.351534/1267.910339. These are approximations, not exact continuum reference values. For the frozen B3 field at h=2, changing Emin/E from 1e-6 to 1e-8 shifts compliance by +0.002165%; changing it to 1e-4 shifts it by -0.180106%. This is below the provisional 1% ersatz-dependence flag for that specific test; it does not validate density-to-solid conversion.

The canonical optimizer converges in 114 evaluated designs in each of two repetitions. Raw/physical arrays and numerical history, excluding timings, are identical; final coarse compliance is 1802.453026 N·mm. Final raw change is 0.00898742 and the last-eleven compliance spread is 0.000204500. This demonstrates reproducibility on the recorded machine, not a proof of optimality or bitwise portability. The canonical evidence contains no failed audit attempts and records 463 completed factorizations and 496 RHS solves, including both optimization repetitions and the two sensitivity sweeps.

## 8. Failures, review fixes and evidence integrity

The first benchmark launch stopped in Windows memory telemetry after its first completed numerical evaluation. The preserved [failure record](../results/phase1_v040/optimization_1/failure.json) shows a ctypes integer-conversion error. Explicit Windows handle/pointer argument types fix the telemetry call. That aborted launch remains separate from the successful canonical run; it is not presented as a mechanics failure or omitted evidence. The successful run records a process peak working set of 2,132,164,608 bytes. Checkpoints and the conservative P2 allocation guard are not an OS-enforced 8 GiB hard limit.

Review also separated numerical Gram reciprocity from the exact algebra utility's symmetry precondition; added missing-provenance, duplicate/mixed-identity, cached-gap and compliance-monotonicity rejection checks; and made per-mesh verification flags consistent with the complete audit. Completed primary and independent solves are charged before later postprocessing can fail; failures preserve known partial fields, diagnostics and timing. Unknown factorization stages inside unsuccessful solves are not invented. Failure work remains in elapsed/CPU totals, and conditional audit costs must not be described as complete selector costs.

The code and specification snapshots in the canonical run identify the exact implementation used for its successful calculations. Subsequent defensive bookkeeping changes do not retroactively alter those files. Final regression tests, including injected failures and false-acceptance fixtures, are recorded in [tests_final.xml](../results/tests_final.xml). The analytical evidence has its own source hashes and uses the same unchanged Q4/P2 kernels. See [tests/test_evaluation.py](../tests/test_evaluation.py) for status, identity, reciprocity and failure-preservation tests.

## 9. What remains unverified

Q4 and P2 share supplied geometry, frozen material coefficients and NumPy/SciPy sparse linear algebra. Independently written shapes, constitutive construction, assembly and traction integration reduce shared implementation risk; they do not exclude a shared modeling or geometry-interpretation error and are not a third-party certification.

Sharp corners prevent treating raw peak stress as a mesh-converged strength result. The regular-region percentile and unrestricted maps are diagnostics with no yield allowable or failure criterion. Several model responses trigger the displacement/120 mm >1% review flag: for example, B1's finest P2 pad maximum is 13.7412 mm. Numerical agreement for a geometrically linear system does not establish that small-deformation assumptions remain adequate at that load. No nonlinear correction, physical load rating or safety conclusion is inferred.

No physical specimens, measured material data, joints/contact, plasticity, fatigue, out-of-plane instability or manufacturing uncertainty were validated. The two unresolved binary candidates need further solution verification or an explicitly changed geometry before receiving numerical acceptance. No Phase 2 selective-verification advantage follows from these results.
