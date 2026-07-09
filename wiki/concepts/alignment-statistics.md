---
type: concept
title: "Alignment statistics (identity / similarity / gaps)"
tags: [alignment, algorithm]
sources:
  - docs/Evidence/ALIGN-STATS-001-Evidence.md
source_commit: 9e4eb79fe421228e018c91fb093940f3feff24a0
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: align-stats-001-evidence
      evidence: "Test Unit ID: ALIGN-STATS-001 ... Algorithm: Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Alignment Formatting"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:global-alignment-needleman-wunsch
      source: align-stats-001-evidence
      evidence: "Statistics are computed over the reported aligned region (columns of an existing pairwise alignment such as the EMBOSS needle global alignment); denominator = alignment Length including gap columns."
      confidence: high
      status: current
---

# Alignment statistics (identity / similarity / gaps)

The **post-alignment metric layer**: given an existing pairwise alignment (from
[[global-alignment-needleman-wunsch|global]], [[semi-global-alignment-fitting|semi-global]],
or local alignment), it summarizes how similar the two sequences are and renders the
familiar EMBOSS/BLAST three-line display. Validated as [[align-stats-001-evidence|ALIGN-STATS-001]]
against EMBOSS `needle`, NCBI BLAST, and `pseqsid`.

## The three metrics (EMBOSS / BLAST convention)

All three share the **same denominator: the full alignment Length, including gap columns**.
Percentage = count / Length × 100.

- **Identity** — count of columns that are an *exact* character match. (BLAST "Identities".)
- **Similarity** — count of columns that are identical *or* score positively under the
  substitution matrix. A superset of Identity, so **Similarity% ≥ Identity%** always. (BLAST
  "Positives".) The operative rule across EMBOSS, BLAST, and pseqsid: **positive substitution
  score ⇒ similar**.
- **Gaps** — count of columns where either side is a gap (`-`). Gap columns are neither
  identities nor mismatches.

Because the denominator is fixed, the identity / similar-only / mismatch-only / gap fractions
**partition to 100%**.

### Canonical worked example (EMBOSS needle, HBA_HUMAN vs HBB_HUMAN, EBLOSUM62)

```
# Length: 149
# Identity:   65/149 (43.6%)
# Similarity: 90/149 (60.4%)
# Gaps:        9/149 ( 6.0%)
# Score: 292.5
```

The 25-column gap between Identity (65) and Similarity (90) is exactly the non-identical
columns that still score positively. This protein example validates the **formula and
denominator**; the Seqeron implementation runs on DNA with a simple Match/Mismatch model and
does not ship EBLOSUM62, so it is a formula-level cross-check.

## DNA vs. graded scoring

- **SimpleDna model (Match +1, Mismatch −1):** no non-identical column ever scores positively,
  so **Similarity equals Identity**.
- **Positive-scoring mismatch (Mismatch = +1):** every non-identical column is counted similar,
  so **Similarity > Identity**.

## srspair markup (alignment formatting)

The three-line display puts a markup line between the two sequence lines. Legend:
`|` identical · `:` similar (substitution score > 1.0) · `.` small positive score · space =
mismatch or gap. Seqeron's DNA scoring uses a single integer Match/Mismatch scalar, so it
emits only `|` (identical) or `:` (a configured positive non-identical score); the graded `.`
tier is unreachable and never emitted — a rendering-only choice that affects no counted
statistic.

## API contract (outside the algorithm spec)

Empty alignment ⇒ `AlignmentStatistics.Empty` / `""` (denominator undefined); null alignment
⇒ `ArgumentNullException`; non-positive `lineWidth` ⇒ `ArgumentOutOfRangeException`.

This is one [[test-unit-registry|test unit]] (ALIGN-STATS-001); see
[[algorithm-validation-evidence]] for the evidence-artifact pattern behind its validation.
