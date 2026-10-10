# PROMETHEUS Phase 0 research package
Research date: 2026-10-10. **Decision: NARROW.**

**Later update, 2026-10-11:** Phase 1 implementation and measured evidence are available in [PHASE1_DECISION.md](PHASE1_DECISION.md), [BASELINE_RESULTS.md](BASELINE_RESULTS.md) and [NUMERICAL_VALIDATION.md](NUMERICAL_VALIDATION.md). The original Phase 0 statements below describe the earlier research stage and are preserved as history. Reproduction commands are in the [root README](../README.md).

Start with [DECISION_MEMO.md](DECISION_MEMO.md), then [RESEARCH_STATUS.md](RESEARCH_STATUS.md) for the evidence boundary. No bracket performance, novel algorithm or physical validation is claimed.

## Required deliverables
| Document | Purpose |
|---|---|
| [RESEARCH_CHARTER.md](RESEARCH_CHARTER.md) | Scope, assumptions, scientific rules and success gates |
| [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md) | Methods, primary research, failure modes and competing work |
| [COMPETITOR_ANALYSIS.md](COMPETITOR_ANALYSIS.md) | Commercial/open tools, APIs, licenses and reproducibility |
| [MATHEMATICAL_FOUNDATION.md](MATHEMATICAL_FOUNDATION.md) | Governing equations, units, design/filter/gradient definitions |
| [RESEARCH_HYPOTHESES.md](RESEARCH_HYPOTHESES.md) | Five falsifiable hypotheses and ranked conditional choice |
| [BASELINE_EXPERIMENT.md](BASELINE_EXPERIMENT.md) | Frozen L-plate benchmark and established comparators |
| [ARCHITECTURE.md](ARCHITECTURE.md) | CPU Python modules, interfaces, provenance and AI strategy |
| [VALIDATION_PROTOCOL.md](VALIDATION_PROTOCOL.md) | Numerical verification, fairness and later physical testing |
| [IMPLEMENTATION_ROADMAP.md](IMPLEMENTATION_ROADMAP.md) | Small milestones and proceed/stop gates |
| [REFERENCES.md](REFERENCES.md) | Deduplicated bibliography with supported claims and access limits |
| [DECISION_MEMO.md](DECISION_MEMO.md) | Recommendation and single highest-information next action |
| [RESEARCH_STATUS.md](RESEARCH_STATUS.md) | Completed work, actual checks and unresolved evidence |

Additional [adversarial review](ADVERSARIAL_REVIEW.md), [source/search log](sources/prior_art_search.md), source ledgers and exploratory utilities support the package.

## Reproduce only the Phase 0 checks
From the workspace root:

~~~powershell
python research/exploratory/linear_load_envelope.py
python research/exploratory/build_bibliography.py
python research/exploratory/check_research_package.py
~~~

The first script checks synthetic linear algebra and benchmark arithmetic. It is not an FEA solver. The second renders the bibliography from existing source annotations; it does not re-check online sources. The third checks required documents, local links, JSON and citation identifiers, not scientific validity.
