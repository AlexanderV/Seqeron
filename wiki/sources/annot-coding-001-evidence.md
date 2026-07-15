---
type: source
title: "Evidence: ANNOT-CODING-001 (Coding potential — CPAT hexamer usage-bias score)"
tags: [validation, annotation]
doc_path: docs/Evidence/ANNOT-CODING-001-Evidence.md
sources:
  - docs/Evidence/ANNOT-CODING-001-Evidence.md
source_commit: 9c179d10f250f0459ff34603c92634ceb60cf68c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ANNOT-CODING-001

The validation-evidence artifact for test unit **ANNOT-CODING-001** (Coding Potential
Calculation — the CPAT hexamer usage-bias score). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is
summarized in [[coding-potential-hexamer-score]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13) — the CPAT paper (Wang et al. 2013, *NAR* 41(6):e74)
  for the log-likelihood hexamer-usage-bias definition and sign convention; the CPAT/lncScore
  reference implementation `cpmodule/FrameKmer.py` (`kmer_ratio`, `word_generator`,
  `kmer_freq_file`) for the exact stepping, averaging, and pseudo-score branches; and the
  EMBOSS `tcode` manual (Fickett 1982 TESTCODE) as a related, not-implemented alternative
  coding measure.
- **Algorithm spec** — score = **mean of `ln(coding[k] / noncoding[k])`** over in-frame
  hexamers; window `wordSize = 6`, `stepSize = 3` (codon-boundary hexamers via
  `word_generator`); natural-log base; positive ⇒ coding, negative ⇒ noncoding.
- **Datasets** — hand-derived worked example `ATGAAA…` → **0.34657359027997264** (mean of
  ln 4 and ln 0.5); pseudo-score cases (coding-only hexamer ⇒ +1, noncoding-only ⇒ −1).
- **Corner cases** — seq shorter than wordSize ⇒ 0; hexamer missing from either table skipped
  (not counted); both-zero-in-both-tables hexamer `continue`d and **does not** increment the
  scored count (verified verbatim against canonical `liguowang/cpat` and `WGLab/lncScore`,
  2026-06-15); `N`-containing hexamers absent from tables so skipped.
- **Recommended coverage** — MUST tests for the mean log-ratio value, ±1 pseudo-scores, sign
  invariant, short-sequence guard, and missing-key skip; SHOULD tests for in-frame stepping
  and the null/argument-validation contract; COULD test for table-unit invariance.

## Assumptions (from the artifact)

1. **No-scorable-hexamer sentinel.** The reference `kmer_ratio` returns `-1` on a caught
   division-by-zero (frame-0 helper), but the two reference branches disagree (the public path
   returns `sum/count`). The C# port returns **0**, matching the documented
   `len(seq) < wordSize → 0` "no information" semantics. Only affects inputs with zero scorable
   hexamers, where the score is meaningless either way.
2. **Table units (counts vs proportions).** `kmer_freq_file` stores raw integer counts while
   CPAT's distributed tables are normalized proportions; `kmer_ratio` divides stored values
   directly, so it is unit-agnostic as long as both tables share units (the difference is an
   additive constant `ln(Σcoding/Σnoncoding)`). Tests use small explicit tables for exactness.

No contradictions between sources; the only deviation is the `-1 → 0` sentinel choice above,
which does not affect any real (nonempty) score.
