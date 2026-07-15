---
type: source
title: "Evidence: CODON-USAGE-001 (Codon Usage Analysis — CalculateCodonUsage + CompareCodonUsage)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-USAGE-001-Evidence.md
sources:
  - docs/Evidence/CODON-USAGE-001-Evidence.md
source_commit: dda05cc3a1a2b473ced5663ee7beeeb062b11907
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-USAGE-001

The validation-evidence artifact for test unit **CODON-USAGE-001** (Codon Usage Analysis —
`CalculateCodonUsage(string)` + `CompareCodonUsage(string, string)` on `CodonOptimizer`). One
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. This unit is the **raw codon-usage table + distribution-comparison** view of the
codon-usage family: it produces an unnormalized codon-count table and compares two such tables
with a Total Variation Distance similarity — the distinct measure captured in
[[codon-usage-comparison]]. See [[test-unit-registry]] for how units are tracked.

## The two methods

- **`CalculateCodonUsage`** — splits a coding sequence into non-overlapping triplets, counts each
  codon, and returns `Dictionary<string, int>` of **raw counts** (not frequencies, not RSCU). It
  converts `T→U` internally (RNA representation), uppercases, and ignores an incomplete trailing
  codon. This is the same non-overlapping-triplet tally that underlies
  [[relative-synonymous-codon-usage|RSCU]]'s `CountCodons` — but left **unnormalized**, so it is
  the raw base of the family, not a codon-usage-bias measure.
- **`CompareCodonUsage`** — normalizes both count tables to frequency distributions (count / total
  codons) and returns a **Total Variation Distance (TVD) similarity** in `[0, 1]`:

      Similarity = 1 − ( Σ_c |f₁(c) − f₂(c)| ) / 2

  where `f_i(c)` is the frequency of codon `c` in sequence `i`, summed over all codons. Empty input
  → similarity `0` (no data), not NaN or exception. See [[codon-usage-comparison]] for the measure.

## Datasets (oracles)

Raw-count oracles:

- `AUGGCUGCU` (M-A-A) → `{AUG:1, GCU:2}`.
- Concatenation of all 64 standard RNA codons → 64 distinct keys, each count 1.
- Invariant: Σ counts = total codons in sequence.

TVD-similarity oracles (all analytically derivable from the formula):

- **Identical** — `AUGGCUGCACUG` vs itself → `1.0` (TVD 0).
- **Disjoint** — all-`UUU` vs all-`GGG` → `0.0` (Σ|f₁−f₂| = 1+1 = 2).
- **2/3 shared** — `AUGGCUAUG` vs `AUGUUUAUG` → `2/3` (Σ = 0 + 1/3 + 1/3).
- **sim = 0.5** — `CUGCUGCUGCUA` vs `CUACUACUACUG` → 0.5 (Σ = 1/2 + 1/2 = 1).
- **sim = 0.75** (symmetry) — `AUGAUGCCCUUU` vs `AUGUUUUUUCCC` → 0.75, and equal both ways.
- **sim = 0.75** (low diff) — `AUGAUGAUGCCC` vs `AUGAUGCCCCCC` → 0.75.
- **sim = 0.25** (high diff) — `AUGAUGAUGAUG` vs `AUGCCCCCCCCC` → 0.25 (Σ = 3/4 + 3/4 = 3/2).
- **Empty** — `""` vs `""` → `0`.

## Reference-table (Kazusa) verification

The unit's predefined organism codon tables were verified against the **Kazusa Codon Usage
Database** (March 2026): E. coli K12 W3110 (species 316407, 4,332 CDS / 1,372,057 codons),
S. cerevisiae (4932, 14,411 CDS / 6,534,504 codons), H. sapiens (9606, 93,487 CDS / 40,662,582
codons) — all 64 relative fractions match. Method: Kazusa per-thousand frequencies converted to
relative fractions per amino acid, compared with implementation values to 2 decimal places.

## Assumptions and deviations (from the artifact)

1. **Comparison uses TVD, not cosine / correlation.** Wikipedia lists cosine similarity and
   correlation coefficients among comparison metrics; the implementation uses TVD-based similarity
   (`1 − Σ|f₁−f₂|/2`) instead. Every expected value in the test suite is analytically derivable
   from the TVD formula, and the proven properties (identity 1.0, symmetry, range [0,1], disjoint
   → 0) follow from TVD theory. A deliberate metric choice, not a departure from a fixed spec.
2. **Empty sequence → similarity 0** (convention: no data → 0), rather than NaN or an exception.
3. **`CalculateCodonUsage` returns counts, not frequencies** — normalization happens inside
   `CompareCodonUsage`. Incomplete trailing codons are ignored (Kazusa / EMBOSS convention);
   `T`/`U` treated as biologically equivalent; mixed case handled; 1–2-nt input → empty table.

**Contradictions:** none — Wikipedia (codon usage bias / degeneracy), the Kazusa table format,
and Sharp & Li 1987 (per-amino-acid normalization) agree on the underlying codon-usage biology;
the TVD similarity is an implementation metric choice the sources do not contradict.

## Related units

Same codon-usage family: [[seq-codon-freq-001-evidence|SEQ-CODON-FREQ-001]] (the normalized,
frame-aware frequency analog of this raw-count table), [[codon-cai-001-evidence|CODON-CAI-001]]
(uses codon frequencies for CAI), [[codon-rare-001-evidence|CODON-RARE-001]] (rare codons from usage
tables),
[[codon-opt-001-evidence|CODON-OPT-001]] (optimization consumes usage), and the RSCU units
[[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]] / [[codon-rscu-001-evidence|CODON-RSCU-001]].
