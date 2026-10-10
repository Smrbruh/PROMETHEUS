# Focused prior-art and decision audit
Access date/cutoff: 2026-10-10. Read alongside [literature_search.md](literature_search.md) and [tools_search.md](tools_search.md).

## Search sequence and representative exact queries
1. “topology optimization uncertain loads worst case semidefinite eigenvalue adaptive exchange method”
2. “topology optimization adaptive sampling multifidelity Bayesian optimization model discrepancy robust design”
3. “adversarial topology optimization counterexample guided robust structural design”
4. “AK MCS active learning reliability method combining Kriging Monte Carlo simulation Echard Gayton Lemaire 2011”
5. “multifidelity optimization trust region model management Alexandrov Lewis Gumbert Green Newman 2001”
6. “topology optimization validation manufacturing postprocessing uncertainty experimental rank design performance”
7. “site:patents.google.com/patent topology optimization uncertainty”
8. “site:patents.google.com/patent generative design manufacturing constraints”
9. Exact titles of MuFASa and Hyun et al. multifidelity robust TO.
10. “Bayesian calibration of computer models Kennedy O Hagan 2001”
11. “Bayesian experimental design under model discrepancy Huchen Yang Dong Wu 2025”
12. “multi information source optimization knowledge gradient Poloczek Wang Frazier 2017”
13. “constrained ranking selection feasibility Bayesian value information engineering design simulation”
14. “multi fidelity best arm identification cost aware verification engineering optimization”

Primary evidence only was used for scientific claims: author/institution repositories, conference proceedings, publisher previews and public patent text. Secondary hits and generic summaries were leads, not authority.

## Decision changes caused by evidence
- Worst-case uncertainty formulations rule out claiming load-counterexample generation itself as novel.
- Active reliability and multifidelity robust TO rule out the broad sampling/uncertainty combination.
- Physical autonomous experimentation and discrepancy modeling rule out claiming the whole simulation-to-test loop.
- MISO/knowledge-gradient, best-arm identification and constrained ranking/selection directly challenge the finite-catalogue proposal. Its acquisition is explicitly conventional.
- Mathematical linearity supplies a much stronger angle-only baseline than repeated scenario solves.
- Final recommendation: empirical discrepancy audit first; a new method only if measured need and a specific prior-art delta survive.

## Access issues and follow-up
Direct ScienceDirect access to several URLs failed with 403/internal errors; indexed publisher abstracts/previews remained available. The UCL record for R04 and arXiv/NeurIPS records were accessible. Access labels in prior_art.json deliberately distinguish preview evidence from full-paper derivation review.

Trust-region model-management work was found, but not added to the root bibliography because complete primary author metadata was not confirmed in the retrieved page; avoid guessing authors from third-party aggregators. The exact title and DOI 10.1016/S0377-0427(00)00424-6 remain a follow-up lead.

The 2024 “Optimal Multi-Fidelity Best-Arm Identification” author preprint, https://arxiv.org/abs/2406.03033, was opened; its revision is dated 26 May 2025. It is a follow-up to the verified 2022 paper, not evidence that its guarantees apply to deterministic reconstruction bias.

Patent US11468634B2 was read for technical disclosure. CN111339645A and CN106650147B surfaced for uncertain-load topology optimization; the former direct open failed and full claims were not reviewed. Machining-constraint and lattice-design patent families also surfaced. They are search leads, not verified enforceability or clearance results.

No negative search result is treated as proof of novelty. Needed before a publication claim: full-text closest-method review, forward/backward citations, same-budget reproductions and specialist review.

