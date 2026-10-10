# Mechanics evidence notes

Search date: 2026-10-10. Sources M01–M06 are recorded in `mechanics.json`; document citations use these stable IDs. Searches were targeted at baseline implementation, projection/mesh issues, verification terminology, independent FEM tooling and stress singularities. This is a focused mechanics evidence search, not an exhaustive systematic review.

Queries included:

- `site.topopt.mek.dtu.dk 88 line topology optimization code efficient topology optimization in MATLAB 2011 Andreassen`
- `site.topopt.mek.dtu.dk robust formulation manufacturing tolerant topology optimization Wang Lazarov Sigmund 2011 projection`
- `site.comsol.com mesh refinement stress singularities reentrant corner finite element`
- `site.osti.gov verification validation computational solid mechanics ASME V V 10`
- `scikit-fem paper Gustafsson McBain 2020`
- `site.asme.org "V&V 10" "2019"`
- `"Efficient topology optimization in MATLAB using 88 lines" "10.1007"`
- `"Material interpolation schemes in topology optimization" "s001580050070"`
- `site.calculix.de "Four-node plane stress element"`

The direct DOI request for top88 returned an error; its DOI, year, authors, pages and venue were verified from the DTU institutional publication record. Its author-hosted uncorrected PDF was accessible. The material-interpolation DOI open also returned an error and no independent source record for that paper was added to this subset. This does not imply that the article itself is unavailable or invalid.

The publisher page for M02 exposes a 2022 correction. Because its details were not inspected, only the narrow abstract-level caution about simple projection is used; no claimed guarantee from a corrected formula is relied upon.

A direct old ASME standards URL failed, but the current official VVUQ hub and its open 2025 educational module were accessible. The numerical tolerances in the protocol are proposed by this research specification, not attributed to the standard. The paid standard was not purchased or read.

CalculiX's official site was accessible. A search exposed a historical MIT mirror of its CPS4 documentation noting expansion into 3D elements. That old mirror was not used to pin current CalculiX behaviour or to promise identical Q4 results. The initial independent verifier is therefore specified at the weak-form level with quadratic triangles and must be implemented and tested, not assumed equivalent from an element name.

The documents contain no claimed FE result. Domain area, pad area, material volume, unit conversion and element counts are arithmetic consequences of the declared geometry. Equations and the two-RHS angle envelope are derivations of the specified linear model; they are not novel research claims. Material constants are deliberately synthetic, so no alloy grade or certified allowable is invented.

Open issues before implementation: pin and licence-review exact software versions; execute the specified manufactured-solution fixture; profile a prototype solve; verify derivative tolerances under realistic stiffness contrast; determine whether the 300-evaluation budget reaches the proposed stop criterion; calibrate convergence tolerances; audit grey-to-binary geometry discrepancy; establish a physical material/fixture model only if physical testing is later authorized and justified.
