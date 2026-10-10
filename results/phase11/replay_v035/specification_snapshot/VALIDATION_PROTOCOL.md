# Verification and validation protocol

Status: **specified checks; not evidence that they have been run**. Applies to `PROMETHEUS-01-v0`. Numerical acceptance checks below are engineering research choices and must be calibrated on analytic tests and the pilot. They are not certification limits or values prescribed by a standard.

## 1. Distinguish three questions

Code verification asks whether the equations were implemented correctly. Solution verification asks how accurately a particular discretized calculation approximates those equations. Physical validation asks whether the model predicts measured behaviour for its intended use. ASME's public verification module distinguishes code and solution verification and describes reference-solution and mesh-refinement approaches [M03, [ASME educational module](https://www.asme.org/getmedia/e446ce87-543d-4e2d-ba2c-67a18b64bc5b/vvuq-intro-module-part-2-verification.pdf)]. The full paid standard was not reviewed; no ASME compliance is claimed.

Independent solver agreement is valuable evidence for numerical correctness. It does not validate shared assumptions such as perfect clamping, isotropy or plane stress. A cross-check that imports the same erroneous load vector cannot independently verify that load definition.

## 2. Gate V0: problem and units

Before solving, assert all geometry dimensions, region membership, area, thickness and volumes from [MATHEMATICAL_FOUNDATION.md](MATHEMATICAL_FOUNDATION.md). In h=2: 1,600 active cells, 150 preserved and 1,450 variable cells; 6,400 mm² domain area; 12,800 mm³ total material budget; 3,000 mm³ preserved volume. Reject units missing from the manifest. E in MPa means N/mm²; density in kg/mm³ is only used to convert volume to mass here.

Independently integrate force and moment from boundary segments. Nominal resultant is (0,−1,000) N and moment about the origin is −120,000 N·mm. The equivalent resultant location is (120,20) mm. Verify the fixed boundary has zero displacement and leaves no rigid-body motion of the full solid domain. Plotting a boundary is useful for inspection but does not replace checking tagged coordinates.

## 3. Gate V1: element, assembly and derivative tests

| Test | Procedure and proposed acceptance |
|---|---|
| Element geometry | Positive Jacobian at every quadrature point; counterclockwise node order; exact area for regular squares. |
| Rigid-body modes | An isolated Q4 has three 2D rigid modes and five positive deformational modes, within relative eigenvalue tolerance 10⁻¹⁰. Verify translations and in-plane rotation produce zero strain. |
| Constant-strain patch | Prescribe an affine displacement field on a small irregular conforming patch; recover constant strain and affine nodal displacements to relative 10⁻⁹ on a well-scaled test. |
| Analytic tensile strip | Ω=[0,100]×[0,10] mm², t=5 mm; right traction 1 MPa, left ux=0, only (0,0) uy=0, free top/bottom. Exact ux=x/E, uy=−νy/E mm for E in MPa; F=50 N. Check displacement, uniform stress, reaction and energy to relative 10⁻⁸. Do not clamp the left edge in y. |
| Convergence test | Use the smooth manufactured field defined below on h=5,2.5,1.25 mm grids. Q4 displacement L² error should approach order 2 and strain-energy norm error order 1; do not demand this smooth-solution rate at the singular L corner. |
| Volume derivative | Compare hj to central directional differences, including directions near each pad and the cutout. Relative error <10⁻⁷ for the linear filter case. |
| Compliance derivative | Interior densities and ten archived unit-norm directions; sweep step sizes 10⁻²…10⁻⁶. Require a central-difference convergence region and relative error <10⁻⁴ at at least two adjacent steps; use an absolute criterion for near-zero directional derivatives. |
| Pad and chain rule | Perturb only free variables; pads remain exactly solid. Check gradient entries affected by padding and normalized filter boundaries explicitly. |
| Load superposition | Compare a directly assembled 13° load solve with two-RHS reconstructed displacements/compliance to relative 10⁻⁸. Verify Q symmetry. |
| Limiting cases | Doubling load doubles u and quadruples C; scaling all E (including Emin) by two halves u and C; full solid C does not exceed any softer same-domain density field under the same load. |

For derivative error use |dFD−g·v|/max(|dFD|,|g·v|,10⁻¹⁰ in the corresponding derivative units), document the absolute residual as well, and avoid clipped perturbations at variable bounds. Nonlinear projection experiments need new derivative checks for both volume and compliance. A single lucky finite-difference step is not adequate evidence.

For the manufactured solution use Ωtest=[0,L]×[0,H] with L=100 mm, H=40 mm, the same E, ν and t, and impose the exact zero displacement on **all** boundary nodes:

```text
ux = a sin(kx x) sin(ky y),   uy = b sin(kx x) sin(ky y),
a=0.01 mm, b=0.02 mm,       kx=π/L, ky=π/H,
φ=sin(kx x)sin(ky y),       χ=cos(kx x)cos(ky y),
λps=Eν/(1−ν²),             μ=E/[2(1+ν)],
body_force_x = a[(λps+2μ)kx²+μky²]φ − (λps+μ)b kx ky χ,
body_force_y = b[μkx²+(λps+2μ)ky²]φ − (λps+μ)a kx ky χ,
body force units: N/mm³.
```

The forcing follows from −div σ(u), including the mixed derivative term; its weak-form load is ∫Ωtest t Nᵀ body_force dA. Use the exact analytical field and derivatives in error integrals, start with 4×4 quadrature for the sinusoidal body force and error norms, and confirm stability when quadrature is increased. Q4 stiffness remains full 2×2 integration. Report measured orders on both refinements, flag substantial departure from the expected asymptotic trend, and inspect implementation before relaxing a threshold. If the independent P2 triangle route uses this same smooth fixture, its expected asymptotic energy/L² orders are 2/3 under the usual regularity and mesh assumptions. This is a specified fixture, not an already executed test. No optimization result is accepted before the independent test suite runs.

## 4. Gate V2: every-solve checks

| Quantity | Proposed acceptance / action |
|---|---|
| Free-DOF residual | ||Kff uf−ff||₂/max(||ff||₂,1 N) ≤10⁻⁸ |
| Symmetry | ||K−Kᵀ||F/max(||K||F,1 N/mm) ≤10⁻¹² before any deliberate symmetrization; do not hide assembly errors by averaging K. |
| Equilibrium | Reaction resultant plus applied force ≤10⁻⁷ times max(total force magnitude,1 N); moment error ≤10⁻⁷ times max(F×120 mm,1 N·mm). |
| Energy | Positive C; |fᵀu−uᵀKu|/max(|C|,10⁻¹² N·mm) ≤10⁻⁸; compare assembled element strain energy as well. |
| Numerical health | Finite coordinates, matrices, u, C and gradients; successful factorization; no invalid Jacobians. Record condition estimates and warn on strong sensitivity to ordering/tolerance. |
| Optimization feasibility | Pads exactly ρ=1; raw bounds respected within 10⁻¹²; volume ≤ budget +10⁻⁶ VΩ. |
| Support/load transfer | Correct tagged entities remain present after every representation change; do not attach load to nearest surviving node without an explicit new model. |

A small algebraic residual alone does not control forward error in an ill-conditioned model. Inspect near-mechanisms, soft void paths and topology; repeat the final density analysis at Emin/E0∈{10⁻⁸,10⁻⁶,10⁻⁴} with the design frozen. A >1% compliance shift from the primary setting is a provisional flag for ersatz dependence, not an invitation to select the most favourable result. A binary material-only verifier is necessary to expose connections carried only by ersatz material.

## 5. Gate V3: fixed-design mesh and independent-solver comparison

Analyse frozen B0, B1 and B3 at h=2,1,0.5 mm. Hold geometry, thickness, loads, boundary extents, physical density and any stress-exclusion distances fixed. For grey fields, assign each child the parent's density; do not refilter during the convergence sequence. For binary geometry, subdivide the same union of squares; do not select a new set of solid pixels. Separately document optimization mesh dependence by reoptimizing at fixed physical r=6 mm.

Report C, load-patch mean displacement, maximum pad displacement, material-region maximum displacement, reactions and energy for each mesh. Proposed initial resolution gates are:

- |C0.5−C1|/|C0.5| <1% and |δpad,0.5−δpad,1|/max(|δpad,0.5|,10⁻¹² mm) <2%.
- The h=2→1 and h=1→0.5 changes must be reported; a single small last difference does not establish an asymptotic regime. Oscillation or a worsening trend requires refinement/diagnosis.
- Converged primary and independently assembled quadratic-triangle solutions should agree within 1% in compliance and 2% in pad displacement. These are investigation thresholds, not guaranteed accuracy estimates.

The independent path should integrate the same constitutive law directly with thickness and independently reconstruct the boundary traction. Use vector P2 triangles in scikit-fem, or another audited independent solver. For P2, the correct **plane-stress** Lamé parameter is λps=Eν/(1−ν²), with μ=E/[2(1+ν)]; importing the 3D λ formula without reduction would instead change the problem. Quadrature must integrate the chosen polynomial stiffness and boundary load accurately. Preserve coefficient discontinuities by aligning mesh edges to the frozen parent-cell boundaries.

Record all differences in element order, quadrature, mesh, DOF constraints, coefficient mapping and solver. If disagreement exceeds the gate, isolate load mapping, stiffness law, geometry transfer and discretization before changing tolerance. The use of scikit-fem is a proposed independent implementation path, not a claim that the package automatically validates the benchmark [M05, [software paper](https://joss.theoj.org/papers/10.21105/joss.02369)].

If successive differences support an asymptotic error model, report an estimated convergence rate and extrapolation with its assumptions. Do not force Richardson extrapolation onto oscillating results or imply that an estimated error bar is a mathematical certificate. If refinement exceeds the PC budget, label the metric `resolution_unverified` and report the resource boundary.

## 6. Geometry and stress checks

On a binary candidate construct **edge-connected** solid-cell components (sharing only a vertex does not count). Both pads and every retained solid component must belong to the support-connected component; otherwise fail connectivity and report the offending area. No silent removal, weak spring stabilization or invented support is allowed. Exported geometry must retain the loaded 20 mm edge, have no self-intersections, and agree in area with the recorded candidate. Record minimum ligament and corner radii as measured diagnostics; the density filter does not guarantee them.

The v0 L corner at (40,40), clamp endpoints (0,80),(40,80), and load-patch endpoints (120,10),(120,30) can generate singular or strongly localized stresses. Define the diagnostic regular region as material integration points at least **6 mm** from those five points and, for a density field, in cells with ρ≥0.9. Report its measure and use an area-weighted 95th percentile von Mises stress there. If that region is empty, report `not_defined`; do not fabricate a zero. Also report the unrestricted maximum over material cells and identify new geometric corners created by binarization.

These diagnostics do not cover every possible singularity in a pixel geometry and cannot become a strength pass criterion. Inspect stress convergence locally and retain the complete maps. A rising peak near a sharp re-entrant corner is not cured by refining indefinitely; a physical fillet or an appropriate fracture/local model changes the engineering question [M04, [official singularity explanation](https://www.comsol.com/blogs/how-identify-resolve-singularities-model-meshing)]. A later strength experiment must use declared smooth geometry, actual material data, applicable failure criteria and process-informed allowables before making a stress feasibility claim.

## 7. Uncertainty, discrepancy and fair method evaluation

Keep aleatory variability (a distribution justified by repeated observations), epistemic parameter uncertainty, numerical error, model-form discrepancy and manufacturing deviation in separate fields. Until data exist, θ∈[−20°,20°] is a bounded scenario set, not a statistical distribution. The exact two-RHS angular envelope already resolves this set for each frozen linear model.

For the adaptive-verification hypothesis, predefine a finite catalog and a trustworthy reference label for each case. Reference labels must themselves pass mesh/solver checks; unresolved cases remain unresolved. Freeze the discrepancy measure and decision threshold before test evaluation, after calibrating on a pilot. Compare random selection, simple deterministic priority rules, established uncertainty/active-learning rules and any proposed method using identical initial observations and cumulative cost accounting. The hypothesis document specifies the competing methods and decision metrics.

Do not claim generalization from a random split of strongly correlated optimization iterates. Split by geometry family, optimization run or parameter regime; keep all near-duplicates together. Keep final verification labels hidden from selection and hyperparameter tuning. Include failures, timeouts, invalid geometry, repairs and training in the cost ledger. Analyse paired differences across repeated seeds, report uncertainty and effect size, and investigate whether an apparent gain disappears after including setup cost or a stronger non-AI selector.

No calibrated failure probability, material safety factor, or surrogate confidence guarantee exists at Phase 0. A predictive interval that is calibrated on the observed catalog is not a certificate over arbitrary future designs.

## 8. Gate V4: physical validation is a later, separate study

Physical testing is justified only after numerical and manufacturability gates and a competent laboratory review. It is not required to complete Phase 0. The future test must use the **actual** manufactured geometry, measured thickness and mass, characterized material, loading fixture and in-plane or three-dimensional model appropriate to that fixture. A real L plate gripped along an edge is not automatically equivalent to a bolted angle bracket.

Before fabrication, choose a process and confirm feature sizes, achievable fillets, tolerance, surface condition, tool access and material/process variation with its operator. Use a fixture and calibrated load frame with a rated load path, shielding as appropriate, and no person in the failure path. A qualified supervisor sets the working range and abort limits based on specimen, machine and fixture, not the nominal 1,000 N research setting. Do not infer that applying 1,000 N to an unvalidated prototype is safe.

The proposed measurement sequence is: inspect dimensions and defects; characterize machine/fixture compliance; measure a solid reference specimen; perform low-load reversible load–unload cycles; measure actual force and displacement using a load cell plus independent displacement measurement or DIC; record hysteresis, slip, residual displacement and temperature; increase only within the approved elastic test range. Register a separate supervised destructive protocol if collapse testing becomes necessary.

Compare predicted and measured force–displacement slopes and full-field displacement/strain where available, including instrument uncertainty, repeated specimens, build orientation/batch and fixture variability. Three specimens per configuration can be a logistics pilot, not a statistical power justification. Choose confirmatory sample size and an engineering-relevant error tolerance after estimating variability, before final testing. Keep model calibration specimens separate from validation specimens. Report prediction error and combined uncertainty rather than retrospectively selecting only well-matched locations.

Material yielding, unmodelled contact/slip, out-of-plane motion, significant residual displacement or a systematic slope discrepancy is evidence to revise the model. It is not evidence that more topology optimization will solve the discrepancy. No physical testing or material validation has occurred in this phase.

## 9. Evidence bundle and release gate

Each candidate needs input and geometry hashes; units and all boundary tags; raw and physical density arrays; material assumptions; filter/optimizer/solver settings; mesh and quadrature; software versions and numerical backends; iteration history; random generator and actual initial arrays; timings and peak memory; conversion/repair history; diagnostics; independent verifier input and output; and a clear state: `proposed`, `evaluated`, `numerically_verified`, `geometry_failed`, `resolution_unverified`, or `physically_validated_for_declared_use`.

Publish negative results and failed cases with the successful ones. An acceptance status must identify the manifest, candidate, model and metrics to which it applies. Reproducibility requires replaying the recorded computation without consulting an LLM. Any assistant-generated requirement or load interpretation remains a proposal until schema, units, feasibility and human engineering intent are checked; the solver output remains authoritative only for the model it actually solved.
