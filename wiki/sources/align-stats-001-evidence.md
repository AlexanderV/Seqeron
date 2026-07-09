---
type: source
title: "Evidence: ALIGN-STATS-001 (Alignment statistics — identity/similarity/gaps)"
tags: [validation, alignment]
doc_path: docs/Evidence/ALIGN-STATS-001-Evidence.md
sources:
  - docs/Evidence/ALIGN-STATS-001-Evidence.md
source_commit: 9e4eb79fe421228e018c91fb093940f3feff24a0
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ALIGN-STATS-001

The validation-evidence artifact for test unit **ALIGN-STATS-001** (Pairwise Alignment
Statistics — Identity / Similarity / Gaps — and alignment formatting). One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
metrics themselves are summarized in [[alignment-statistics]]. See [[test-unit-registry]] for
how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13) — EMBOSS `needle` docs (release 6.6) for the
  identity/similarity/gaps/score definitions and the HBA_HUMAN vs HBB_HUMAN worked example;
  EMBOSS AlignFormats for the srspair markup legend; NCBI BLAST (NBK1734) for
  Identities/Positives/Gaps; a percent-identity denominator corroboration; and `pseqsid`
  (github.com/amaurypm/pseqsid) as a reference implementation.
- **Metric spec** — count / Length × 100 with the denominator **including gap columns**;
  Similarity = identical *or* positively-scoring columns (superset of Identity); the
  "positive substitution score ⇒ similar" rule confirmed independently by EMBOSS, BLAST, and
  pseqsid.
- **Datasets** — the EMBOSS 149-column example (65/90/9 → 43.6% / 60.4% / 6.0%) as a
  formula-level cross-check; a hand-built 9-column DNA alignment (SimpleDna ⇒ Similarity =
  Identity); and a positive-scoring-mismatch case (Mismatch = +1 ⇒ Similarity > Identity).
- **Corner cases** — gap columns counted in the denominator; Similarity ≥ Identity always;
  DNA simple-model ⇒ Similarity = Identity; dashes are neither identities nor positives.
- **Recommended coverage** — MUST tests for the formula/denominator, DNA equality, positive
  mismatch, gap counting, and srspair markup; SHOULD tests for empty/null/`lineWidth`
  contracts; COULD tests for line-wrapping and the partition invariant.

## Assumption (rendering-only)

The srspair `.` tier (small positive score) collapses away: the DNA scoring model exposes a
single integer Match/Mismatch scalar with no graded positive scores, so only `|` (identical)
or `:` (a configured positive non-identical score) is emitted. Documented as
non-correctness-affecting for the counted statistics. No contradictions with other sources.
