# Literature search record

Search/access date and inclusion cutoff: 2026-10-10 (user time context, Asia/Qyzylorda).
Scope owner: scientific literature stream of PROMETHEUS Phase 0.
Output: 33 source records L01–L33 in literature.json and the structured LITERATURE_REVIEW.md.

## Method

This was a targeted, iterative web search with backward/adjacent tracing, not a preregistered systematic review. No subscription database export or complete citation graph was available. The queries below are representative exact search strings from the session; repeated URL opens and minor query variants are omitted. Results were selected for primary-source authority, relevance to structural brackets and ability to constrain proposed novelty. Search rankings and citation counts were not treated as quality evidence.

Institutional repositories and author project sites were used to cross-check metadata and retrieve openly accessible material. Method existence, bibliography, abstract-level evidence and complete methodological inspection are distinguished in each record's access field. A URL returning a search excerpt while direct fetch failed was recorded as such.

## Query families and outcomes

| Area | Representative exact queries | Selected evidence |
|---|---|---|
| Density methods | "Material interpolation schemes in topology optimization"; site.orbit.dtu.dk efficient topology optimization 88 lines code 2011 Andreassen | L01, L02 |
| Foundational topology | "Bendsøe" "Kikuchi" "1988" "Optimal topologies" | Found 1988 DOI metadata; direct publisher 403. Not retained as a fully reviewed item. Seminal density foundation represented by L01. |
| Level set | site.sciencedirect.com level set structural optimization Allaire Jouve Toader 2004 | Author original PDF L03; unrelated similar-title papers excluded. |
| ESO/BESO | "A simple evolutionary procedure for structural optimization" Xie Steven; "Huang" "Xie" "Convergent and mesh-independent solutions" | L04/L05; author-uploaded manuscript checked against publisher record. |
| Multiobjective/gradient | "A fast and elitist multiobjective genetic algorithm" Deb 2002 pdf; "The method of moving asymptotes" Svanberg 1987 | L06/L30 |
| Fair classical baselines | "On the usefulness of non-gradient approaches in topology optimization"; "A new generation 99 line" Ferrari Sigmund 2020 | L31/L32 |
| Geometric robustness | "On projection methods, convergence and robust formulations in topology optimization" | L07 and correction notice; correction contents not audited. |
| Worst-load optimization | topology optimization worst case uncertain loads de Gournay Allaire Jouve; "A general framework for robust topology optimization" Thore Holmberg Klarbring 2017 | L08/L09: strong counterexample-generation prior art. |
| ROM/reanalysis | "reduced basis" "topology optimization" "Amir" | L10; thesis used only for discovery, original publication selected. |
| Differentiable simulation | "JAX-FEM" paper 2023; JAX-FEM differentiable GPU accelerated 2023 Xue Liao | L11 and publisher program summary. |
| Neural representations | "TOuNN" "Topology" neural networks Chandrasekhar Suresh; site.github.com UW-ERSL TOuNN license | L12 and GPL-3.0 author repository. |
| PINNs | "Physics-informed neural networks: A deep learning framework" Raissi Perdikaris Karniadakis 2019; "Characterizing possible failure modes in physics-informed neural networks" | L13/L14; limits kept within studied equations. |
| GNNs | "Learning Mesh-Based Simulation with Graph Networks" | L15, author project and official repository README. |
| Generative frontier | site.arxiv.org topology optimization bracket benchmark diffusion 2024 2025 2026 | L16/L17; generated-topology performance distinguished from verified manufactured performance. |
| Bracket datasets | site.arxiv.org SIMJEB simulation ready dataset 3d mechanical brackets; "DeepJEB" "041703" DOI; "DeepJEB" "10.1115" | L18/L19, datasets not downloaded. |
| Sequential constrained choice | "Bayesian optimization with inequality constraints" Gardner 2014 | L20; original PDF corrected the malformed third-author HTML metadata. |
| Reliability | "AK-MCS" "Echard" "2011"; "AK-MCS" site.sciencedirect.com | L21; later adaptive/multifidelity neighbors forwarded to hypothesis-review stream. |
| Calibration | "Bayesian calibration of computer models" Kennedy O Hagan 2001 | L22; publication citation year distinguished from a conflicting online metadata date. |
| DfAM | "An additive manufacturing filter for topology optimization of print-ready designs"; "topology optimization" "experimental validation" bracket additive manufacturing 2024 | L23/L26; process surrogate PATO found by related result and opened as L27. |
| Physical autonomous experiments | "Autonomous" "mechanical" "design" "Bayesian" robot experimentation 2024; "A Bayesian experimental autonomous researcher for mechanical design"; "Using simulation to accelerate autonomous experimentation" | L24/L25; direct prior art against broad closed-loop discrepancy claims. |
| RL | "reinforcement learning" "topology optimization" Hayashi Ohsaki; site.arxiv.org topology optimization reinforcement learning benchmark 2024 2025; "Structural design through reinforcement learning" "2026" DOI | L28/L33; 2025 final journal metadata cross-checked with 2024 preprint. |
| Nonlinear frontier | "Material interpolation schemes in topology optimization" yielded an explicit 2026 primary result | L29: publisher page explicitly dated 29 April 2026. |

## Access and provenance limitations

- Publisher DOI redirects for the 1988 homogenization paper, Huang–Xie BESO and AK-MCS failed in direct opening. BESO/AK-MCS publisher article URLs and indexed primary excerpts were available; their access fields describe the narrower evidence.
- Some ScienceDirect/Wiley articles exposed only abstracts, excerpts or program summaries. No full-text-only numerical assertion is presented as independently checked.
- SOgym's publisher direct page intermittently returned 403; its publisher search record supplied the final 2025 venue/DOI, and the original 2024 arXiv record and author code were accessible.
- A guessed TopoDiff GitHub URL returned 404. It is not cited as an implementation; the accessible author project page supplies the authoritative code/data links.
- Direct browsing of the MeshGraphNets directory was restricted, while its official README and author project page were accessible. Repository/data licenses were therefore not assumed.
- A SimJEB paper click failed, but an original Eurographics PDF and publisher metadata were located independently. Its dataset landing page returned insufficient license text. The website repository's MIT license does not establish dataset terms.
- DeepJEB's author page states a license and sample count, but the downloadable archive/manifest was not audited. Its displayed cluster/sample descriptions should not be treated as a reconciled inventory.
- The 2011 projection article displays a correction published in 2022. Retrieve and review the correction before implementing or reproducing its formulas.
- Search result phrases such as “published last year” and institutional PDF “downloaded” dates were discarded when bibliographic dates were explicitly available.
- No author-reported numerical improvement was reproduced. No source-code benchmark, training run, specimen production or physical test occurred in this review stream.
- Independent replication searches were opportunistic, not comprehensive. “Not identified” is a coverage statement, not a claim that replication does not exist.
- The paper list is selective across disciplines. Fatigue, contact, nonlinear buckling, metal-AM process qualification, testing standards, CAD reconstruction and formal verification each warrant deeper work before entering the implementation scope.
- Patent claims, legal status and freedom-to-operate are not resolved by this literature review. Any patents found elsewhere in Phase 0 are discovery leads, not legal conclusions.

## Minimum follow-up before a novelty claim

Read complete original methods and supplements for the selected nearest comparators; review corrections; follow backward/forward citations and alternate terminology such as multifidelity best-arm identification, constrained ranking-and-selection, multi-information-source optimization and value of information. Audit code/data licenses at pinned versions. Reproduce at least one strong baseline and one nearest proposed-method comparator with identical budgets. Register the finite benchmark catalog and held-out discrepancy tests before method tuning.

A comparison showing no benefit, or showing that a standard allocation method already works, is an informative Phase 1 result.

