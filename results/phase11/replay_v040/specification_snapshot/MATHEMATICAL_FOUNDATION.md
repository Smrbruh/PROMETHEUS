# PROMETHEUS-01 mathematical foundation

Status: **proposed specification, not an executed simulation or a qualified structural design**. Research date: 2026-10-10. Numerical constants below define a reproducible synthetic problem; they are not measured material properties, certified allowables, or demonstrated performance.

## 1. Scope and canonical problem

PROMETHEUS-01-v0 is an **in-plane, plane-stress L-shaped plate**. It is a tractable precursor to a three-dimensional bracket with joints, contact, out-of-plane instability, manufacturing defects, and measured material properties. Passing v0 establishes performance only within its mathematical model.

Use millimetres, newtons, and megapascals throughout the mechanics calculation. Coordinates have their origin at the lower left; positive x is right and positive y is up.

| Quantity | Definition |
|---|---|
| Domain | Ω = ([0,120] × [0,40]) ∪ ([0,40] × [40,80]) mm² |
| Thickness | t = 5 mm, uniform |
| Domain area and full volume | AΩ = 6,400 mm²; VΩ = 32,000 mm³ |
| Preserved support pad | Ps = [0,40] × [70,80] mm² |
| Preserved load pad | Pl = [110,120] × [10,30] mm² |
| Preserved area and volume | AP = 600 mm²; VP = 3,000 mm³ |
| Designable domain | D = Ω \ (Ps ∪ Pl), area 5,800 mm² |
| Fixed boundary ΓD | y = 80, 0 ≤ x ≤ 40; ux = uy = 0 |
| Loaded boundary ΓL | x = 120, 10 ≤ y ≤ 30, length L = 20 mm |
| Nominal load | F = (0, −1,000) N, distributed uniformly on ΓL |
| Other boundaries | Traction free; no gravity or other body force |
| Synthetic solid material | E0 = 70,000 N/mm², ν = 0.30, ρmass = 2.70 × 10⁻⁶ kg/mm³ |
| Material budget | V ≤ 0.40 VΩ = 12,800 mm³, including preserved pads |
| Equivalent budget mass | m ≤ 0.03456 kg = 34.56 g; full solid mass 86.4 g |

The cutout above the horizontal leg is excluded from the mesh, stiffness matrix, filter, and volume denominator. It is not an ersatz low-density region. Only the **top edge** of the support pad is fixed; fixing its entire area would define a different problem. There are no bolt holes or inferred fastening assumptions.

The available designable material volume is 9,800 mm³, equivalent to an unfiltered mean density 1,960/5,800 = 0.33793103448 in D. This number is **not** the initialization value when pad values enter a density filter; initialization must satisfy the actual filtered-volume equation.

## 2. Linear elasticity and discretization

Let u=(ux,uy), engineering strain ε(u)=(ux,x, uy,y, ux,y+uy,x)ᵀ, and

```text
D(E,ν) = E/(1−ν²) [[1, ν, 0], [ν, 1, 0], [0, 0, (1−ν)/2]].
σ = D ε(u).
```

Plane stress imposes σzz=τxz=τyz=0; it is distinct from plane strain. For any admissible virtual displacement v vanishing on ΓD,

```text
∫Ω t ε(v)ᵀ D ε(u) dA = ∫ΓL vᵀ q ds,
q = F/L   [N/mm].
```

Equivalently the physical surface traction is τ=F/(tL)=[0,−10] MPa for the nominal case. Do not multiply a line load by thickness a second time. Assemble boundary forces by integrating shape functions over each loaded edge; endpoint nodes receive half the share of an interior node on a uniform first-order edge mesh. The summed force and moment must equal those of a uniform load acting at (120,20) mm.

Use conforming, displacement-based, four-node bilinear quadrilaterals (Q4), two displacement degrees of freedom per node, and full 2×2 Gauss integration. On an element e,

```text
Ke = t ∫Ae Beᵀ D(Ee,ν) Be dA = Ee K̂e,
K = assemble(Ke);  Kff uf = ff;  uD = 0.
```

K̂e includes thickness and unit Young's modulus. The optimization mesh is a square grid with h=2 mm and 1,600 active elements. Delete nodes unused after removal of cutout elements. Full integration avoids the additional hourglass modes of an unstabilized reduced-integration Q4. Mesh refinement assesses in-plane discretization error; it cannot test omitted three-dimensional physics.

## 3. Variables, filter, and SIMP interpolation

Use raw design variables xj∈[10⁻³,1] for j∈D; set xj=1 for preserved elements j∈P. The positive lower bound is an optimizer convention, separate from the ersatz stiffness below. Let cj be element centroids and Aj their areas. For every active element,

```text
wij = max(0, r − ||ci−cj||₂) Aj,     r = 6 mm,
Hij = wij / Σk wik,
x̃i = Σj Hij xj.
```

Neighbours are active cells of Ω only. Preserved cells participate in the weighted average with x=1. Set physical density ρi=x̃i for i∈D, and **override ρi=1 for i∈P after filtering**. This fully defines filter behaviour at the external boundary, cutout, and pads. This normalized filter need not preserve volume near boundaries and pads; always constrain the resulting physical density. r remains 6 mm under mesh refinement, rather than a fixed number of elements.

The first baseline has **no projection** (identity mapping). Use

```text
Ee = Emin + ρeᵖ (E0−Emin),    p=3,    Emin/E0=10⁻⁶.
V(x) = Σe t Ae ρe,           m(x)=ρmass V(x).
```

Emin is a numerical device, not a claim that empty space supports load. SIMP intermediate densities are an optimization relaxation, not a prescription for an isotropic porous manufactured material. Density filtering and an optimality-criteria optimizer form an established educational baseline; the authors' top88 reference supplies a useful independent algorithmic reference, not evidence of novelty [M01, [DTU paper and code](https://www.topopt.mek.dtu.dk/apps-and-software/efficient-topology-optimization-in-matlab)].

The problem is

```text
minimize_x C(x) = fᵀu(x)
subject to K(x)u=f; u|ΓD=0;
           V(x)/VΩ ≤ 0.40;
           10⁻³ ≤ xj ≤ 1 (j∈D); ρi=1 (i∈P).
```

This nonconvex finite-dimensional problem has no asserted globally optimal solution. With a fixed load, compliance equals twice strain energy: C=fᵀu=uᵀKu=2U, in N·mm (1 N·mm=0.001 J). Lower compliance means lower load-weighted displacement and greater stiffness **for that loading**. It does not minimize every displacement component or certify strength.

## 4. Derivatives and the OC baseline

For fixed, design-independent loads and homogeneous displacement conditions,

```text
∂C/∂ρe = −p ρe^(p−1) (E0−Emin) ueᵀ K̂e ue.
∂V/∂ρe = t Ae.
```

Rows corresponding to preserved physical densities have zero derivative with respect to any raw design variable. For j∈D the correct derivatives are

```text
gj = ∂C/∂xj = Σi∈D (∂C/∂ρi) Hij,
hj = ∂V/∂xj = Σi∈D t Ai Hij.
```

Do not replace hj by a constant, omit the transpose of H, differentiate overridden pad densities, or reuse sensitivity-filter equations for this density-filter problem. Central directional differences must verify both derivatives before optimizing.

With move limit mOC=0.20, use the proposed update

```text
xj,new(λ) = clip[xj sqrt(max(0,−gj)/(λ hj)),
                 max(10⁻³,xj−mOC), min(1,xj+mOC)].
```

Choose λ≥0 by bisection using **V computed after the filter and pad override**. Expand the positive upper bracket until its candidate is feasible; consider λ→0 as the upper-move-limit candidate rather than dividing by zero. Use the feasible endpoint of the bisection. Since H has nonnegative entries, increasing λ cannot increase volume. Nonfinite derivatives, hj≤0, or substantial positive compliance sensitivities are diagnostic failures. The zero-sensitivity case is allowed. The update is a standard OC heuristic for this problem, not a proof of convergence or a universal method for arbitrary stress constraints.

For a later, separately versioned projection experiment, define η∈(0,1) and

```text
Pi(s;β,η) = [tanh(βη)+tanh(β(s−η))]/
            [tanh(βη)+tanh(β(1−η))],
Pi'(s) = β sech²(β(s−η)) /
         [tanh(βη)+tanh(β(1−η))].
```

Then ρi=Pi(x̃i) for i∈D and the chain rule includes Pi'(x̃i) in **both** gj and hj. At β=0 use the limit Pi(s)=s, Pi'=1. Projection continuation, erosion/intermediate/dilation formulations, and minimum feature controls are outside v0; all stages and evaluations must count toward a comparison budget. A single projection or a radius setting alone does not establish a guaranteed minimum feature size or manufacturing robustness [M02, [Wang, Lazarov and Sigmund](https://doi.org/10.1007/s00158-010-0602-y)].

## 5. Multiple loads and an exact low-dimensional robustness check

For a declared set of loads fk and nonnegative fixed weights αk summing to one, weighted compliance is Σk αk fkᵀuk and its sensitivity is the same weighted sum of single-load sensitivities. A minimax objective instead minimizes z subject to Ck(x)≤z. These are different objectives; do not compare their nominal compliance without also reporting their intended worst-case metric.

For this benchmark fix the same load patch and let

```text
F(θ)=1,000 [sin θ, −cos θ] N,    θ∈[−20°,20°].
```

Let B=[bx,by] contain dimensionless nodal force distributions whose resultant in their respective coordinate direction is one. Solve KU=B once with two right-hand sides and form Q=BᵀU. The derivation from superposition gives

```text
C(θ)=1,000² [Qxx sin²θ − 2Qxy sinθ cosθ + Qyy cos²θ].
dC/dθ=1,000² [(Qxx−Qyy) sin(2θ) − 2Qxy cos(2θ)].
```

Evaluate interval endpoints and all roots of the derivative within the interval; handle constant C explicitly. This gives the exact angular maximum **for the frozen discretized linear model**, within floating-point solve error. For unrestricted load directions the maximum is F² λmax(Q); that bound need not be attained inside the restricted interval. Degrees are input/output only; trigonometric calculations use radians.

An AI counterexample selector cannot justify an efficiency claim by comparing against many repeated solves of this two-dimensional load space. The exact two-RHS method is the required competitor. Load amplitude and uniform E changes are similarly cheap: C scales with F²/E and u with F/E when the full stiffness, including Emin, scales with E. Patch relocation, shape changes, contact, plasticity, or buckling can invalidate this shortcut and require a separately specified model. Exactness of the angle search does not establish exactness of the underlying physical model.

## 6. Alternative engineering objectives and their status

| Objective or constraint | Mathematical form | v0 status / missing evidence |
|---|---|---|
| Maximum displacement | max over a declared material region of ||u||₂≤δlim | Report on both pads and on cells with ρ≥0.9; no application-derived δlim is known. Displacements in ersatz void are not component displacements. |
| Mass minimization | min ρmass V with stiffness/displacement limits | Later alternative once a meaningful limit exists; do not silently substitute for fixed-volume compliance. |
| Local stress | σVM=√(σxx²−σxxσyy+σyy²+3τxy²), σVM≤σallow | Diagnostic only: no yield allowable, stress-relaxation method, or filleted part is established. |
| Buckling | Positive load multiplier from tangent/geometric stiffness eigenproblem, followed by imperfection-sensitive nonlinear assessment | Deferred. The plane-stress model omits out-of-plane modes. |
| Fatigue | Life or damage constraint under a specified cycle spectrum and material/process data | Deferred; static elastic stress cannot establish fatigue life. |
| Manufacturing | Minimum ligament, tool access, radius, tolerance, process-specific defect and cost constraints | Filter is only regularization; geometry inspection and process-qualified constraints needed later. |
| Robustness | max over a declared uncertainty set C(x,ξ), or probability of a declared limit-state violation | Angular set is available; a probability of failure requires a justified distribution and threshold. |
| Multiobjective design | Pareto trade-off among mass, stiffness, process cost, etc. | Use separate mass budgets initially; objective weights do not establish a complete Pareto front. |

Optimization feasibility means satisfying the volume/design/equilibrium restrictions actually imposed. Independent verification additionally checks model transfer, mesh error, connectivity, and numerical reliability. A density solution can be optimization-feasible and fail geometric verification.

## 7. Stress singularities and physical limits

The domain has a sharp re-entrant corner at (40,40), abrupt clamping endpoints, and traction-patch endpoints. A diverging pointwise peak stress under refinement is consistent with an idealized singularity; it must not be presented as a converged physical strength prediction. Global compliance can converge while peak stress does not [M04, [COMSOL singularity guidance](https://www.comsol.com/blogs/how-identify-resolve-singularities-model-meshing)].

Report stress maps and a clearly labelled diagnostic statistic away from fixed **physical** exclusion regions defined in the validation protocol. Always retain the unrestricted peak and its mesh trend. Do not shrink exclusions to a fixed element count and do not use exclusions to assert structural safety. A later stress-constrained benchmark needs declared fillets, a material allowable, suitable stress constraints, and an independent local convergence study.

Small strain, small rotation, linear isotropic elasticity, uniform thickness, perfect geometry, perfect fixation, and static in-plane loading are assumptions. No assertion follows about plastic collapse, fracture, fatigue, impact, corrosion, contact, buckling, or arbitrary load directions. A maximum pad displacement divided by 120 mm exceeding 1% is a provisional **model-review trigger**, not a universal small-deformation theorem or strength limit. Independent calculations and physical tests serve different purposes; see [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md).
