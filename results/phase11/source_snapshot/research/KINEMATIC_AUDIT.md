# Phase 1.1 kinematic audit of archived linear solutions

Date: 2026-10-11, Asia/Qyzylorda. **The simple B1 structures have large omitted finite-rotation strain terms at the stated load. Numerical refinement of the linear equations does not resolve that model limitation.** This audit computes kinematics from saved linear displacements; it does not solve nonlinear equilibrium, estimate nonlinear response error, establish strength, or identify a safe load.

The complete batch is [kinematics/complete](../results/phase11/kinematics/complete/manifest.json), with [full metric records](../results/phase11/kinematics/complete/records.json), [CSV](../results/phase11/kinematics/complete/summary.csv), [separate density context](../results/phase11/kinematics/complete/supplemental_context.json), and [measured cost](../results/phase11/kinematics/complete/cost.json). Existing Phase 1 files and statuses are preserved. This audit implements the kinematic portion of the [Phase 1.1 protocol](REFERENCE_LABEL_PROTOCOL.md).

## Scope, inputs and definitions

The twelve original catalogue entries are B0, B1, B3_binary and B3_threshold at each of the three volume settings. For each available field, the source is its original **Q4 h=0.5 mm material-only binary mesh**, displacement basis and record. The manifest hashes source NPZ files, source arrays, historical evaluations, meshes, original manifests, catalogue, code and protocol. Before/after input hashes establish read-only preservation for this batch.

Right/down basis columns already include the 1,000 N resultant. A displacement at angle theta is `sin(theta)*u_right + cos(theta)*u_down`. The audit uses -20°, 0° and +20° without another factorization. These are three reported scenarios, not the proven maximum of every kinematic statistic over the continuous interval.

There are **30 computed original candidate-angle records**. The two disconnected .35 reconstructions have no displacement field: their six angle records explicitly say `not_evaluated_geometry_failed_no_field` and contain null metrics. They are not assigned zero displacement or zero strain. The repaired candidate had no selected completed field for this batch and is listed separately in [omitted.json](../results/phase11/kinematics/complete/omitted.json); a later audit must use its own matching geometry hash, never its parent's displacement.

At Q4 Gauss points, derivatives are with respect to original coordinates X:

```text
H_ij = du_i/dX_j
e = (H + H.T)/2
Q = H.T @ H / 2
E_GL = e + Q = ((I+H).T @ (I+H) - I)/2
```

These are **in-plane 2x2 tensors**. No out-of-plane stretch or nonlinear plane-stress constitutive closure is inferred. Tensor shear e_xy is used, not engineering shear gamma_xy=2e_xy. For a tensor T, the reported area L2 norm is `sqrt(sum(w*||T||_F²))`; area RMS divides it by `sqrt(sum(w))`. The P95 is the inverse weighted empirical distribution of pointwise Frobenius norms. Area weights are reference mm²; dimensionless-tensor L2 norms consequently have units mm, while RMS and percentiles are dimensionless.

The ratio reported below is `||Q||_L2 / ||e||_L2`. It is a ratio of regional norms, not an average or percentile of pointwise ratios. A zero denominator produces null. The regular region retains points at least 6 mm from the five original singular locations. It **does not exclude all newly introduced pixel corners**. Both unrestricted and regular results remain available.

The audit also records `||H||_F`, the skew-gradient norm, infinitesimal rotation `(H_yx-H_xy)/2`, the proper polar rotation angle of `F=I+H` where det(F)>0, and determinant diagnostics. Infinitesimal skew rotation is not the exact finite rotation angle. Rigid rotation must not produce material strain: the Green–Lagrange tensor has this property, while the infinitesimal tensor does not at finite angles. [COMSOL deformation measures](https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.009.html) documents these distinctions. Its [Lagrangian formulation](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.008.html) defines the material-coordinate deformation gradient and positive-determinant requirement.

## Measured nominal results

All values in this table use the **archived linear Q4 h=0.5 solution**, not P2. `u_max` is the maximum on retained material nodes; the full JSON separately records pad maximum. For these B1/reconstruction cases the global and pad maxima coincide. B0 has global maximum .495609 mm and pad maximum .476145 mm. B0 is repeated in the catalogue for provenance but represents the same full-solid geometry and displacement at every budget.

| Candidate | u_max [mm] | u_max/120 | Global Q/e L2 | Regular Q/e L2 | Global absolute polar-angle P95 [degrees] |
|---|---:|---:|---:|---:|---:|
| B0, each budget | .495609 | .4130% | .026564 | .027580 | .3058 |
| .35 B1 | 22.365135 | 18.6376% | 4.787066 | 4.624617 | 13.6895 |
| .40 B1 | 13.698592 | 11.4155% | 2.465202 | 2.397034 | 8.3706 |
| .40 B3_binary | 1.248259 | 1.0402% | .125605 | .132868 | 1.4047 |
| .40 B3_threshold | 1.265662 | 1.0547% | .133025 | .140688 | 1.4565 |
| .45 B1 | 7.091327 | 5.9094% | 1.002871 | .973243 | 4.6333 |
| .45 B3_binary | 1.078765 | .8990% | .086756 | .091216 | 1.0517 |
| .45 B3_threshold | 1.059567 | .8830% | .082956 | .087188 | 1.0215 |

The original .35 B3_binary and B3_threshold have no corresponding rows of numbers because their geometry failed before analysis. No favorable kinematic value can be inferred for them.

For canonical B1, the global quadratic correction's L2 norm is approximately **2.47 times** that of the linear strain; the regular-region ratio remains **2.40**. The regular quadratic P95 is .0152825, compared with linear-strain P95 .00932184. Its global displacement-gradient P95 is .208091. The omitted term is therefore not confined to the five excluded singular neighborhoods, and the linear kinematic approximation warrants a direct nonlinear-model comparison. This finding is separate from B1's successful spatial-resolution checks for the linear equations.

The canonical ranked reconstruction's regular ratio is .132868 and its regular quadratic P95 is .000444965, compared with linear P95 .00195218. These are substantially smaller than B1's values but do not meet a defined physical-validity criterion: no such kinematic acceptance threshold has been specified. Its prior geometry/domain and numerical-resolution questions remain separate. A displacement/length ratio below or above 1% is a screening observation, not a safety classification.

## Load orientations and density context

| Candidate | Regular Q/e at -20° | At 0° | At +20° |
|---|---:|---:|---:|
| .35 B1 | 4.685470 | 4.624617 | 3.986284 |
| .40 B1 | 2.426562 | 2.397034 | 2.069892 |
| .45 B1 | .968917 | .973243 | .856360 |
| .40 B3_binary | .102208 | .132868 | .146336 |

The largest reported .35 B1 displacement among these angles is 22.907956 mm at -20°. Canonical B1 reaches 13.957292 mm at -20°. Different statistics need not reach their largest sampled value at the same angle; none of these three-angle tables proves a continuous-angle kinematic maximum.

The canonical original B3 **density** field is supplemental context outside the twelve-entry binary catalogue. To avoid interpreting displacement gradients in ersatz void as part deformation, its diagnostic includes only cells with rho>=0.9. At nominal load the maximum over that region is 1.940226 mm (1.6169% of 120 mm); global/regular Q/e ratios on that selected material region are .111039/.111006, and polar-angle P95 is 1.1755°. This region and coefficient model differ from the binary candidates; the numbers are not a like-for-like comparison of complete manufactured geometries.

## Verification, cost and retained exporter failure

The new kinematic component's **14 tests pass**. They cover exact translation, affine extension/shear, irregular Q4 coordinate mapping, exact finite rigid rotations, superposed-rotation objectivity, linear versus quadratic load-field scaling, area-weighted percentiles, finite physical exclusion masks, empty regions, orientation reversal and invalid inputs. Exact finite rigid rotation produces nonzero infinitesimal strain but Green–Lagrange strain at roundoff. These tests verify tensor calculations; they do not validate a nonlinear equilibrium solver.

The successful postprocessing batch took **17.4555 s wall time and 8.171875 s CPU time**, including recorded input reading, hashing, computation and preceding JSON/CSV writes. It used **zero factorizations, zero RHS solves and zero new mechanics solves**. Earlier mechanics costs belong to their original evidence bundles and are not made free by reusing archived fields.

The first batch computed and wrote its JSON results, then its CSV exporter failed on an explicit null kinematic record for rejected geometry. That [failure record](../results/phase11/kinematics/failure.json), original JSON and original source snapshot are preserved. Only the new CSV exporter was corrected. The complete batch was rerun into a new subdirectory; no existing artifact was overwritten and no numeric definition was altered. The failed batch did not retain its complete timing, which is reported as unknown rather than zero or an invented charge. The measured successful-run time therefore is not the total time spent on this task.

Reproduce into a **new directory**:

```powershell
python -m pytest tests/test_kinematic_audit.py -q
python scripts/phase11_kinematics.py --output results/phase11/kinematics_replay
```

## Suitable next model, not a safe-load extrapolation

The next numerical reference should solve geometrically nonlinear, objective **plane-stress equilibrium** for the same frozen geometry with incremental loading, a consistent tangent, equilibrium/energy checks and step-size/refinement studies. A total-Lagrangian formulation using material-coordinate kinematics and work-conjugate stresses is a reasonable research path; a verified corotational formulation is another option when small material strain with finite rotation is justified. The constitutive law and plane-stress thickness/stretch elimination must be stated explicitly. Large-deformation kinematics do not automatically make a small-strain elastic material law appropriate at arbitrary stretch. [Bathe's MIT lecture on elastic constitutive relations](https://ocw.mit.edu/courses/res-2-002-finite-element-procedures-for-solids-and-structures-spring-2010/resources/lecture-16/) discusses the required distinctions between total/updated Lagrangian transformations and large-displacement, small-strain cases.

This recommendation is an engineering inference from the diagnostics and those formulations, not an already tested solution. Freeze whether the load remains a global dead traction/resultant on reference material points or follows the deformed boundary: that choice changes the nonlinear problem. Establish rigid-rotation and appropriate nonlinear benchmark checks before using the new response as a reference. Contact, actual joint compliance, constitutive/material data and out-of-plane effects remain additional model questions.

Adding Q to the saved linear strain does not rebalance internal and external forces, update the tangent, or yield a new displacement solution. Likewise, scaling the current load until a selected diagnostic looks small does not establish a safe operating load. No candidate is approved for manufacture or physical loading by this audit, and no learned selective-verification claim follows.

## Late addendum: successful Q4 fields for the .35 repaired geometry

The earlier complete catalogue batch and its tables remain unchanged. After that batch, the new repaired candidate produced a successful **Q4 h=0.25 mm** solution, while its independent **P2 solve failed the 1e-8 free-residual gate**. The retained [attempt record](../results/phase11/refinement/attempts/v035_repair/h0.25/record.json) and [partial fields](../results/phase11/refinement/attempts/v035_repair/h0.25/fields.npz) permit displacement-only postprocessing. They do not supply an accepted P2 response or a numerically accepted reference label.

The repair retains the pad-connected component of the original .35 ranked reconstruction and removes floating material. Its actual volume fraction is .344375, area 2,204 mm² and mass 29.754 g; it is a changed geometry with a separate identity and is not an exact-equal-mass .35 competitor. The dedicated addendum verifies the actual source-array byte hash and Candidate content hash, catalogue identity, geometry hash, source-record status and Q4 diagnostics. It retains the P2 exception explicitly. The [repair manifest](../results/phase11/kinematics/repair/manifest.json), [full records](../results/phase11/kinematics/repair/records.json) and [CSV](../results/phase11/kinematics/repair/summary.csv) are new artifacts; no prior evidence or acceptance status was overwritten.

| Load angle | Q4 global and pad maximum [mm] | Maximum/120 | Global Q/e L2 | Regular Q/e L2 | Global absolute polar-angle P95 [degrees] |
|---:|---:|---:|---:|---:|---:|
| -20° | 33.147149 | 27.6226% | 45.287166 | 39.806754 | 44.9704 |
| 0° | 43.951818 | 36.6265% | 60.691943 | 53.278786 | 52.8593 |
| +20° | 49.460157 | 41.2168% | 68.422279 | 60.036719 | 55.9966 |

At nominal load, the archived Q4 compliance is 32,269.454 N·mm. The global omitted quadratic strain norm exceeds the linear strain norm by approximately **60.7 times**, and the regular-region ratio is approximately **53.3**. This is strong evidence that treating this linear displacement field as a small-rotation deformation is inappropriate; refining the same linear equations cannot repair that kinematic approximation. It is **not** a measurement of a nonlinear equilibrium displacement, nonlinear compliance error, collapse load, or physical safety. The failed independent solve and any separate domain/geometry interpretation issues remain unresolved. These late calculations neither override the original disconnected geometry rejection nor grant acceptance to its repaired descendant.

The three-angle addendum used **1.83084 s wall time, 1.59375 s CPU time and zero mechanics solves/factorizations/RHS solves**. Its [cost record](../results/phase11/kinematics/repair/cost.json) separately identifies the archived Q4/P2 attempt's 13.2931 s wall cost. Reproduce this addendum alone into another new directory:

```powershell
python scripts/phase11_repair_kinematics.py --output results/phase11/kinematics_repair_replay
```
