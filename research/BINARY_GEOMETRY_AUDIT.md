# Phase 1.1 binary geometry audit

Date: 2026-10-11, Asia/Qyzylorda. **This is a measured geometry audit, with zero mechanics solves.** It preserves the original Phase 1 candidates and their historical labels. New reconstructions and the single declared repair have no inherited numerical-verification or physical-validation status.

The evidence is the [13-entry catalogue](../results/phase11/geometry/catalog.json), [21 threshold reconstructions](../results/phase11/geometry/thresholds.json), [repair record](../results/phase11/geometry/repair.json), and [generation cost and hashes](../results/phase11/geometry/generation_cost.json). Generation follows the separately frozen [reference-label protocol](REFERENCE_LABEL_PROTOCOL.md). The implementation is [binary_audit.py](../prometheus/binary_audit.py), invoked by [phase11_geometry.py](../scripts/phase11_geometry.py).

## What is measured

Every mask is a union of the original h=2 mm squares. The domain has 1,600 cells and area 6,400 mm². One solid cell has area 4 mm², volume 20 mm³ at the declared 5 mm thickness, and synthetic mass 0.054 g at 2.70×10⁻⁶ kg/mm³. The .35/.40/.45 material budgets therefore correspond to 30.240/34.560/38.880 g. All 12 original catalogue masks and the repair retain both pads.

Area is counted directly from retained cells, without contour approximation. In the records, `solid_area_mm2`, `material_volume_mm3`, and `mass_kg` give those measures; `mass_minus_budget_kg` is the signed mass error before applying the existing feasibility tolerance. A negative error means less material than the budget, not an exact-mass match. The allowable area excess remains 10⁻⁶×6,400=0.0064 mm², corresponding to 0.0000864 g with these synthetic constants. No tabulated excess of practical size is hidden by that tolerance.

`edge_component_count` uses shared full cell edges; touching only at a vertex does not join components. `offending_area_mm2` is the material outside the component containing the support pad. Local point contacts are checked separately: at a grid vertex, exactly two diagonally opposite cells are solid and the other two are absent, including cells outside the domain. Their coordinates and participating cell IDs are retained in `local_point_contacts`, with `local_point_contact_detection_checked=true`. A globally connected geometry can still contain such a local contact; this case has a dedicated analytical regression test.

`axis_runs.x` and `axis_runs.y` contain the complete histograms, minima, medians, and maxima of maximal contiguous horizontal and vertical solid runs. Each run is counted once, irrespective of its length. The reported minimum run can occur at a tip, step, island, or pad transition. It is **not a Euclidean minimum ligament, a bound on every load-bearing connection, a guaranteed minimum feature, or a manufacturing limit**. The implementation explicitly leaves `minimum_ligament_mm=null` and `minimum_ligament_status="not_measured"`. Pixel corner counts likewise do not establish a smooth fillet radius or a stress value.

## Original catalogue and the separately identified repair

In the following table, component areas are listed in their deterministic catalogue order. The support component is the first listed component in these measured cases. “Contacts” means local point contacts, not the number of disconnected components. X/Y are the minimum horizontal/vertical maximal-run lengths, in mm.

Δmρ is the signed difference from the archived optimized B3 density-model mass in the same volume family. This is a comparison baseline for B0/B1, which were not reconstructed from B3. For the original rows it is calculated here from the archived [35%](../results/phase1_v035/candidates/B3/physical_source.npy), [40%](../results/phase1_v040_verified/candidates/B3/physical_source.npy), and [45%](../results/phase1_v045/candidates/B3/physical_source.npy) density arrays using 4×5×2.70×10⁻⁶×Σρ kg. The original catalogue's optional `reference_density_area_mm2` field is null; this report does not attribute the derived comparison to that field. The repair record's own reference-density argument is its **parent binary mask**, so its recorded `mass_change_from_density_kg` describes mass removed from that parent. The table below instead uses the same optimized-density reference as the other rows.

| Catalogue key | Components: areas [mm²] | Contacts | Actual volume fraction | Mass [g] | Δm budget [g] | Δmρ [g] | Minimum X/Y runs [mm] |
|---|---|---:|---:|---:|---:|---:|---:|
| v035_B0 | 1: 6400 | 0 | 1.000000 | 86.400 | +56.160 | +56.160000143 | 40 / 40 |
| v035_B1 | 1: 2240 | 0 | .350000 | 30.240 | 0 | +1.433×10⁻⁷ | 8 / 10 |
| v035_B3_binary | 3: 2204, 24, 12 | 2 | .350000 | 30.240 | 0 | +1.433×10⁻⁷ | 4 / 2 |
| v035_B3_threshold | 4: 2128, 8, 8, 8 | 1 | .336250 | 29.052 | −1.188 | −1.187999857 | 4 / 2 |
| v040_B0 | 1: 6400 | 0 | 1.000000 | 86.400 | +51.840 | +51.840000079 | 40 / 40 |
| v040_B1 | 1: 2560 | 0 | .400000 | 34.560 | 0 | +7.911×10⁻⁸ | 10 / 8 |
| v040_B3_binary | 1: 2560 | 0 | .400000 | 34.560 | 0 | +7.911×10⁻⁸ | 4 / 2 |
| v040_B3_threshold | 1: 2556 | 0 | .399375 | 34.506 | −.054 | −.053999921 | 4 / 2 |
| v045_B0 | 1: 6400 | 0 | 1.000000 | 86.400 | +47.520 | +47.520000482 | 40 / 40 |
| v045_B1 | 1: 2880 | 0 | .450000 | 38.880 | 0 | +4.822×10⁻⁷ | 10 / 10 |
| v045_B3_binary | 1: 2880 | 0 | .450000 | 38.880 | 0 | +4.822×10⁻⁷ | 6 / 4 |
| v045_B3_threshold | 1: 2896 | 0 | .452500 | 39.096 | +.216 | +.216000482 | 6 / 4 |
| v035_repair | 1: 2204 | 0 | .344375 | 29.754 | −.486 | −.485999857 | 4 / 2 |

The optimized density areas used for this comparison are 2239.9999893821223, 2559.9999941400115, and 2879.9999642800176 mm². Their tiny deficits from the exact budgets explain the nonzero Δmρ for an exact-budget binary. They are numerical volume-control differences, not measured material tolerances.

The .35 ranked binary retains two non-support components totaling **36 mm²**. The .35 threshold binary retains three non-support components totaling **24 mm²**. Both were already geometry failures in Phase 1, and remain unchanged. Their newly recorded local point contacts are:

| Original mask | Contact vertex [mm] | Participating h2 cell indices (ix, iy) | Original cell IDs |
|---|---|---|---|
| v035_B3_binary | (96, 32) | (48,15), (47,16) | 948, 1007 |
| v035_B3_binary | (90, 34) | (45,16), (44,17) | 1005, 1064 |
| v035_B3_threshold | (102, 30) | (51,14), (50,15) | 891, 950 |

These configurations have zero shared area. A shared-node finite-element representation would impose a local point constraint; the protocol withholds an unambiguous continuum-reference label for such geometry pending a declared interpretation. This audit does not split nodes or modify the old masks. **No local point contacts were found in the original .40 or .45 reconstructions**, so this diagnostic does not explain their previously reported resolution failures. Their remaining mesh-resolution question must be addressed separately.

## Predefined threshold sweep

For each frozen optimized density field, the exact generation rule is `physical_density_ge_threshold_then_solid_pad_override`: retain cells with ρ≥τ, then force both pads solid. No rank adjustment, smoothing, redistribution, or island removal is applied to these 21 masks. The .50 entries reproduce the original threshold masks; their records retain new generation provenance and explicit geometry aliases. All masks have immutable provenance IDs and their own saved arrays in the [threshold archive](../results/phase11/geometry/thresholds.json).

Δmρ below is taken from each threshold record's `diagnostics.mass_change_from_density_kg` and converted to grams. Minimum X/Y run lengths retain the restricted interpretation above.

| Source volume | Threshold τ | Actual fraction | Mass [g] | Δm budget [g] | Δmρ [g] | Components | Contacts | Minimum X/Y [mm] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| .35 | .35 | .435000 | 37.584 | +7.344 | +7.344000143 | 1 | 0 | 4 / 4 |
| .35 | .40 | .396250 | 34.236 | +3.996 | +3.996000143 | 1 | 0 | 4 / 4 |
| .35 | .45 | .369375 | 31.914 | +1.674 | +1.674000143 | 1 | 0 | 4 / 2 |
| .35 | .50 | .336250 | 29.052 | −1.188 | −1.187999857 | 4 | 1 | 4 / 2 |
| .35 | .55 | .289375 | 25.002 | −5.238 | −5.237999857 | 12 | 9 | 2 / 2 |
| .35 | .60 | .252500 | 21.816 | −8.424 | −8.423999857 | 8 | 4 | 2 / 2 |
| .35 | .65 | .200000 | 17.280 | −12.960 | −12.959999857 | 6 | 0 | 2 / 2 |
| .40 | .35 | .488125 | 42.174 | +7.614 | +7.614000079 | 3 | 0 | 2 / 4 |
| .40 | .40 | .457500 | 39.528 | +4.968 | +4.968000079 | 1 | 0 | 6 / 4 |
| .40 | .45 | .426875 | 36.882 | +2.322 | +2.322000079 | 1 | 0 | 6 / 4 |
| .40 | .50 | .399375 | 34.506 | −.054 | −.053999921 | 1 | 0 | 4 / 2 |
| .40 | .55 | .366250 | 31.644 | −2.916 | −2.915999921 | 2 | 2 | 4 / 2 |
| .40 | .60 | .326250 | 28.188 | −6.372 | −6.371999921 | 6 | 4 | 2 / 2 |
| .40 | .65 | .289375 | 25.002 | −9.558 | −9.557999921 | 8 | 7 | 2 / 2 |
| .45 | .35 | .534375 | 46.170 | +7.290 | +7.290000482 | 1 | 0 | 4 / 4 |
| .45 | .40 | .504375 | 43.578 | +4.698 | +4.698000482 | 1 | 0 | 2 / 4 |
| .45 | .45 | .477500 | 41.256 | +2.376 | +2.376000482 | 3 | 1 | 2 / 4 |
| .45 | .50 | .452500 | 39.096 | +.216 | +.216000482 | 1 | 0 | 6 / 4 |
| .45 | .55 | .435625 | 37.638 | −1.242 | −1.241999518 | 1 | 0 | 6 / 4 |
| .45 | .60 | .410625 | 35.478 | −3.402 | −3.401999518 | 1 | 0 | 4 / 2 |
| .45 | .65 | .366875 | 31.698 | −7.182 | −7.181999518 | 4 | 2 | 4 / 2 |

Mass is nonincreasing with increasing threshold, as required by nested retained-cell sets and the fixed pad override. Component counts need not be monotone: an isolated component can disappear entirely as the threshold rises. For example, the .45 family has three components at τ=.45 and one at τ=.50; this does not mean thresholding created a connecting ligament.

Of the 21 masks, 11 are edge-connected and 13 have no local point contacts. Only **three** simultaneously satisfy the declared geometry/pad checks, have no local contacts, and lie within their own mass budget: (.40, τ=.50), (.45, τ=.55), and (.45, τ=.60). None of the seven .35 threshold masks meets that combined geometry-and-mass screen. These are deterministic screening results, not mechanics acceptance. The sweep runs no FEM analyses and supplies no new compliance, displacement, strength, or reference label for the unevaluated masks.

## Why representation and identity remain separate

The optimized B3 field is a SIMP coefficient model with E=Emin+ρ³(E0−Emin) on every active domain cell. Intermediate densities are a numerical relaxation; the positive void modulus is an ersatz device. A binary analysis instead removes void cells and assigns the solid modulus to retained material. Equal modeled mass therefore does not imply equal stiffness fields, equal topology, or transferable verification.

The archived `B3_binary` rule preserves pads and selects the highest-density free cells to reach an exact cell-count budget, using deterministic y/x tie breaking. `B3_threshold` uses ρ≥.5 plus preserved pads and does not enforce an exact count. At .45, that threshold produces an explicit **0.216 g excess**; at .35 it produces less than budget mass but disconnected material. Neither a favorable coarse response nor a lower mass would repair those geometry/constraint issues.

Candidate provenance IDs include source archive digest, source field values, source mesh, threshold or repair rule, and output geometry hash. A different rule or source remains separately recorded even when the output geometry matches an existing geometry alias. An alias identifies matching pixel geometry; it does not inherit a numerical or physical acceptance claim. Original names and historical hashes are retained in their catalogue rows. The repair's full provenance ID is its `candidate_name`; its new numerical evaluator must compute the ordinary field hash from the saved mask rather than treat that provenance ID as the Phase 1 field-hash formula.

## The single declared .35 repair

The rule `retain_pad_connected_component` retains the existing component containing both complete pads and removes the other two edge-disconnected components. It removes nine cells, **36 mm²**, **180 mm³**, and **0.486 g**; it adds no material and makes no other modification. The original ranked binary remains unchanged.

The new mask has 551 cells, area 2,204 mm², volume 11,020 mm³, mass **29.754 g**, and volume fraction **.344375**. Its only component contains both pads, and the local-contact count is zero. Removing the islands changes mass, so this is a new under-budget candidate, not an exact-.35-mass substitute for the original. Its record says `status="generated_not_analyzed"` and `independent_verified=false`; subsequent mechanics evidence must be assigned to this new identity.

Immutable repair ID:

`binary-d51ed6f3d6179ed716f11b90d5e49f2f018436296d4e413422614fe9b0d152bd`

The [repair record](../results/phase11/geometry/repair.json) contains all nine removed (ix,iy) indices and the full pre/post diagnostics. The [saved new mask](../results/phase11/geometry/arrays/v035_repair_d51ed6f3d6179ed7.npy) is separate from the [archived original mask](../results/phase1_v035/candidates/B3_binary/physical_source.npy).

| Provenance item | SHA-256 |
|---|---|
| Original NPY file bytes (`source_input_sha256`) | `b0753c4d637c2d9cf06bb8ab6a2cded0fd6048a534aedcd3af7a09717f2cae35` |
| Original binary values as little-endian float64 (`source_field_sha256`) | `4e439d941c8e3b9c012a58b100af03988706a78d115e59d173b7f36ede46c68e` |
| Original pixel geometry | `9548c7a7d9c574a00c38d435be8bd78bb4d2fb3a89fe215ddce6a6fdb5f8bc2d` |
| New pixel geometry | `eda41e1e69baf7b8964ecb40ef05f4073d9af0c48e5ec4f15fbdb94fc82bb969` |
| New NPY file bytes (`array_file_sha256`) | `1a2febe9572c3a6a5ef6ddd27716017fd8817fcdc45f6d8d9ac3bb5b2e55c47a` |

These hashes cover different objects and are not interchangeable. In particular, a file-byte digest includes the NPY container, while the source-field digest hashes only its canonical numerical values.

## Generation cost, preservation, and verification boundary

The [generation ledger](../results/phase11/geometry/generation_cost.json) records **1.774375 s wall time**, **1.140625 s CPU time**, **zero factorizations**, and **zero RHS solves** for 12 original geometry audits, 21 threshold reconstructions, and one explicit repair. It archives 22 new masks. The timing covers input loading through the preceding artifact writes and input-preservation checks; it excludes process startup and the final cost-file write. Individual records also retain their measured generation times. This one local measurement is not a runtime guarantee or an algorithmic speedup result.

Every loaded historical input was rehashed after generation and remained unchanged. Existing nonempty output directories are refused; the refusal was exercised and the generated artifact bytes remained unchanged. The ledger contains input hashes, generated-artifact hashes, and the exact script/module/protocol hashes. The catalogue identity digest is `25761d24c206bb1f80f1a803d6a08ca3f37c0bd290b53e49e5a50716351da17a`.

The [23 automated geometry tests](../tests/test_binary_audit.py) passed in the focused run. They cover analytically counted areas and runs, disconnected and locally touching cells, a globally connected loop with a local point contact, threshold nesting and pad preservation, deterministic IDs/aliases, unchanged inputs, exact declared island removal, invalid inputs, and absent rather than fabricated minimum-ligament values. They do not constitute mechanics or physical validation.

To reproduce into a new destination:

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
python -m pytest tests/test_binary_audit.py -q -p no:cacheprovider
python scripts/phase11_geometry.py --output results/phase11/geometry_replay
```

This evidence identifies representation, mass, connectivity, and local-contact changes explicitly. It does not grant a finer numerical reference label, resolve the .40/.45 reconstruction convergence failures, or establish safety, strength, manufacturing feasibility, or a learned selective-verification advantage.
