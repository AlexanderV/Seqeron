---
type: concept
title: "Coding potential (CPAT hexamer usage-bias score)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/ANNOT-CODING-001-Evidence.md
  - docs/algorithms/Extended_Annotation/Coding_Potential_Calculation.md
  - docs/Validation/reports/ANNOT-CODING-001.md
source_commit: e748206486a14ab05fe3c14e312e74cd77874af2
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-coding-001-evidence
      evidence: "Test Unit ID: ANNOT-CODING-001 ... Algorithm: Coding Potential Calculation (CPAT hexamer usage-bias score)"
      confidence: high
      status: current
---

# Coding potential (CPAT hexamer usage-bias score)

A sequence-level **coding-potential** measure: how strongly a DNA sequence's in-frame
hexamer usage resembles known coding vs. noncoding sequence. It is the hexamer feature of
**CPAT** (Wang et al. 2013), one of the standard alignment-free discriminators between
protein-coding transcripts and lncRNAs/noncoding RNA. Validated as
[[annot-coding-001-evidence|ANNOT-CODING-001]] against the CPAT paper and the reference
`cpmodule/FrameKmer.py`. See [[test-unit-registry]] for how the unit is tracked.
The two-stage validation run ([[annot-coding-001-report]]) found and fixed a real code
defect in exactly the both-zero branch below (it was being counted, not skipped) — see that
report for the verdict and cross-check.

## The score

For a sequence, slide a window of `wordSize = 6` (a hexamer) with `stepSize = 3` starting at
frame 0, keeping only full-length words — i.e. hexamers on codon boundaries. For each scored
hexamer `k`, take the **natural log of the ratio of its coding to noncoding frequency**,
`ln(coding[k] / noncoding[k])`, and the sequence score is the **mean** of those per-hexamer
log-ratios:

    score = mean over in-frame hexamers of ln( coding[k] / noncoding[k] )

**Sign convention (a testable invariant):** a **positive** score indicates coding, a
**negative** score indicates noncoding. Coding-biased tables give positive scores; noncoding-
biased tables give negative scores.

## Branch / edge behaviour

The averaging follows `FrameKmer.py.kmer_ratio` exactly:

- Both frequencies > 0 ⇒ add `ln(coding[k]/noncoding[k])`, count the hexamer.
- coding > 0, noncoding == 0 ⇒ add **+1** (pseudo-score, avoids log of infinity), count it.
- coding == 0, noncoding > 0 ⇒ add **−1** (pseudo-score), count it.
- coding == 0 **and** noncoding == 0 ⇒ `continue`: the hexamer is **skipped and not counted**
  (not scored as 0).
- Hexamer missing from *either* table (including any `N`-containing hexamer) ⇒ skipped, not
  counted.
- Sequence shorter than `wordSize` ⇒ score 0.

The score is **unit-agnostic**: coding and noncoding tables may be raw counts (as built by
`kmer_freq_file`) or normalized proportions (as in CPAT's distributed tables), as long as both
use the same units — differing units add only a constant `ln(Σcoding/Σnoncoding)`.

## Relation to other coding measures

CPAT's hexamer score is one member of a family of "coding potential" statistics. A distinct,
authoritative alternative is the **Fickett (1982) TESTCODE** statistic (EMBOSS `tcode`), which
combines position-asymmetry and composition values rather than hexamer usage bias. TESTCODE is
recorded in the evidence as a related, **not-implemented** alternative — the registry's
"hexamer frequency bias" title selects the CPAT hexamer score.
