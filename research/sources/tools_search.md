# Tool/competitor research search and access audit

Date: 2026-10-10 (client date). Reviewer scope: commercial capabilities, open engineering tools, candidate numerical/CAD libraries, licensing and implementation implications.

## Method

Searches used the internet search tool, followed by primary vendor documentation, official project repositories, actual license texts and source downloads. Search results from Wikipedia, Reddit and third-party tutorials were not used to substantiate technical claims. Vendor material establishes documented capability, not independent performance.

This is a bounded desk review. It is not an exhaustive feature certification, systematic software census, procurement exercise, legal opinion or executed product benchmark. No credentials, trials, paid APIs, license changes or installed applications were required. No production code was created.

## Representative executed searches

| Area | Query strings used | Follow-up and disposition |
|---|---|---|
| Autodesk | `site.autodesk.com Fusion generative design objectives manufacturing constraints API simulation`; `site.autodesk.com support Fusion simulation generative design API not supported`; `site.help.autodesk.com Fusion simulation study types static stress nonlinear` | Opened generative overview, design criteria, study types and official product-manager API announcement. Updated older “no Simulation API” assumption to limited Insider access. |
| Altair | `site.altair.com inspire topology optimization manufacturing constraints OptiStruct API`; `site.help.altair.com hwsolvers os topology optimization manufacturing constraints multi objective`; `site.help.altair.com licensing Altair Units Inspire OptiStruct 2026`; `site.help.altair.com optistruct multi objective weighted compliance minmax` | Followed current feature/API/response pages and 2025.1 manufacturing documentation. Distinguish Inspire interface from underlying OptiStruct capabilities. |
| Ansys | `site.ansys.com structural optimization topology level set density adjoint manufacturing constraints`; `site.ansyshelp.ansys.com v261 structural optimization density mixable level set`; `site.mechanical.docs.pyansys.com license Mechanical installation licensed`; `site.ansys.com academic students ansys student limits structural 128k`; `site.ansyshelp.ansys.com v261 "Objective" "Weight" "Structural Optimization"` | Search initially surfaced v242/v251. Opened v261 equivalents before finalizing classic-density limitations, level-set controls and objective aggregation. Read student/academic and Python licensing docs. |
| nTop | `site.ntop.com nTop automate implicit modeling design optimization`; `site.support.ntop.com topology optimization constraints static structural material linear`; `site.support.ntop.com design optimization multi objective`; `site.ntop.com legal software license agreement nTop`; `site.support.ntop.com "System Requirements Guide"` | Read official knowledge-base pages, EULA, headless workflow and desktop requirements. Followed October 2026 ONNX surrogate page. Avoided unrelated ntop.org network-monitoring project. |
| FreeCAD/CalculiX | `site.freecad.org license LGPL FEM CalculiX scripting`; `site.dhondt.de CalculiX license GPL features` | Opened official FreeCAD repository/raw LICENSE; CalculiX site, capability list and linked GPL text. No Windows binary installed. |
| Numerical libraries | `site.docs.scipy.org scipy license BSD optimize sparse`; `site.pymoo.org license Apache multi objective NSGA II`; `site.botorch.org license MIT multi objective constrained Bayesian optimization` | Opened current official docs plus raw SciPy, pymoo and BoTorch license files. Inspected scikit-fem and DOLFINx README/LICENSE/COPYING directly. |
| DTU | `site.topopt.mek.dtu.dk 88 line code copyright educational purposes`; `site.topopt.mek.dtu.dk "top88.m" "rights"` | Opened official software page. Web parser could not read source attachment; retrieved official file using read-only HTTP and inspected footer. Reserved-rights/education text, not an inferred MIT license. |
| Geometry | Direct opens of Gmsh official manual and Open Cascade licensing URL/repository exception | Gmsh manual identifies GPL v2-or-later with exception. OCCT old licensing URL redirected; current official page says LGPL 2.1 plus exception. |
| AI engineering | `site.ansys.com products structures ansys simai training data`; `site.neuralconcept.com platform deep learning geometry optimization engineering`; `site.physicsx.ai physics AI engineering simulation platform` | Opened vendor platform statements, SimAI data documentation and SimAI Pro bracket tutorial. No public comparable accuracy benchmark inferred. |

Queries above are representative actual strings, not a claim of exhaustive search coverage. The retained source ledger contains 50 entries with claim-level support and limitations.

## Access issues and resolution

1. FreeCAD GitHub LICENSE initially loaded; a subsequent find call timed out. Raw official LICENSE fetched successfully. The feature webpage returned a fetch error; use official repository description instead.
2. Fusion API Introduction URL returned an error. The official Autodesk product-manager August 2026 announcement was accessible and established limited early access. Private Insider functionality was not accessed.
3. DTU top88 attachment was rejected by the web parser because it is `application/octet-stream`. Read-only PowerShell `Invoke-WebRequest` retrieved the official page link and source; UTF-8 text decoding exposed the author/license footer. No copy was vendored or executed.
4. Open Cascade's historical licensing URL redirected to the current official OCCT page. That page and the official exception text were accessible.
5. Search often ranked older Ansys documentation. Current 2026 R1 HTML pages were opened explicitly; the final analysis uses those pages for method-specific limitations.
6. Search snippets and crawl dates were not used as publication years. Undated project pages are recorded as `n.d.`; dated product pages use displayed dates. Floating repositories were not pinned to commit hashes, so the later implementation must repeat exact-version license/dependency verification.
7. No private enterprise license/API contracts for PhysicsX or Neural Concept were accessible. Marked unresolved, not assumed open, free or unavailable.
8. Commercial software performance, solver reproducibility, exact job concurrency and account entitlement were not tested. Product-specific material/contact/manufacturing combinations remain a later verification task.

## Exclusions and follow-up

- No numerical performance comparison or price estimate.
- No license compatibility opinion for a future bundled product.
- No claim that unadvertised capabilities are absent.
- No audit of every FreeCAD plugin, alternative FE library, optimizer or CAD wrapper.
- Matplotlib/PyVista, testing frameworks, specific experiment trackers and all transitive dependencies still need exact-version selection/license checks at implementation time.
- Validate installed-machine compatibility and a small numerical test before adopting any proposed dependency.

