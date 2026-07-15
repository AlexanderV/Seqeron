---
type: source
title: "Evidence: SEQ-CODON-FREQ-001 (Codon Frequencies — CalculateCodonFrequencies with reading frame)"
tags: [validation, annotation]
doc_path: docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md
sources:
  - docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md
source_commit: ae4a6ae53a125dceb08b5ff21c344ab447afe335
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-CODON-FREQ-001

The validation-evidence artifact for test unit **SEQ-CODON-FREQ-001** (Codon Frequencies —
`SequenceStatistics.CalculateCodonFrequencies(string dnaSequence, int readingFrame = 0)` in the
Analysis assembly). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. This unit is the **normalized,
frame-aware frequency-table** view of the codon-usage family: it returns `count / total` fractions
directly, keyed by DNA codon, with a selectable reading frame. See [[test-unit-registry]] for how
units are tracked, and [[codon-usage-comparison]] for the raw-count sibling it normalizes.

## The method (and how it differs from the raw-count sibling)

`CalculateCodonFrequencies` splits a coding sequence into non-overlapping triplets from the chosen
frame offset, counts each codon, and returns `IReadOnlyDictionary<string, double>` of **frequencies**
(`count / total counted codons`), not raw counts:

    freq(c) = count(c) / Σ_c' count(c')

It is deliberately distinct from `CodonOptimizer.CalculateCodonUsage`
([[codon-usage-001-evidence|CODON-USAGE-001]]), the raw-count end of the family, on four points:

- **Output is normalized** (`double` frequency, count/total) — not `Dictionary<string,int>` counts.
- **Reading-frame parameter** (`readingFrame`, 0/1/2, default 0): counting starts at that offset, so
  the *same* sequence yields a different codon multiset per frame. `CalculateCodonUsage` is frame-0
  only. This is the distinctive new semantic.
- **Ambiguous (non-ACGT) triplets are excluded** from both the per-codon count and the total
  (Kazusa CUTG convention) — a triplet with any `N`/degenerate base contributes nothing.
- **DNA-native keys** (T retained); it does not apply the `T→U` RNA rewrite that
  `CalculateCodonUsage` does.

Shared with the rest of the family: case-insensitive (upper-cased internally), an incomplete
trailing 1–2-nt remainder is not a codon and is dropped, and frequencies over the counted codons
sum to 1.0 (invariant INV-02).

## Sources (accessed 2026-06-14)

- **Kazusa Codon Usage Database (CUTG) README** (authority 5, the canonical convention): codon usage
  is a per-codon count aggregated over all CDS, reported as raw count and **per-thousand** frequency;
  "codons containing ambiguous base were excluded from count." The metric under test = per-thousand
  frequency ÷ 1000 = count/total.
- **EMBOSS `cusp` documentation** (authority 3, reference implementation): five output columns
  (Codon, AA, Fraction, Frequency, Number). Its **Fraction** is a *per-amino-acid* proportion — a
  **different** metric (that is the [[relative-synonymous-codon-usage|RSCU]]/`Fraction` view, not
  this one). Its **Frequency** = expected codons per 1000 bases. The verbatim 386-codon sample output
  cross-checks the count/total identity: CGC 22/386×1000 = 56.995, GGC 23/386×1000 = 59.585 exactly.
- **Wikipedia — Codon usage bias** (authority 4): establishes codon frequency as the unit of
  measurement; cites Sharp & Li 1987 (CAI) and Ikemura 1981 as downstream indices, not needed for
  the raw frequency definition.
- **Nakamura, Gojobori, Ikemura (2000)**, *Nucleic Acids Research* 28(1):292
  (doi:10.1093/nar/28.1.292; authority 1): the paper underlying the Kazusa CUTG database.

## Datasets (oracles)

Cross-check (EMBOSS cusp): total codons Σ = 386; CGC count 22 → 22/386 = 0.0569948… = 56.995/1000;
GGC count 23 → 23/386 = 0.0595854… = 59.585/1000.

Hand-derived exact rationals (direct application of the count/total definition):

| Input | Frame | Codons read | Expected frequencies |
|-------|-------|-------------|----------------------|
| `ATGATGAAA` | 0 | ATG, ATG, AAA | ATG = 2/3, AAA = 1/3 |
| `ATGATGAAA` | 1 | TGA, TGA (AA leftover) | TGA = 1.0 |
| `ATGNNNAAA` | 0 | ATG, NNN *excluded*, AAA | ATG = 1/2, AAA = 1/2 |
| `ATGAA` | 0 | ATG (AA leftover) | ATG = 1.0 |
| `atgaaa` | 0 | ATG, AAA (case-normalized) | ATG = 1/2, AAA = 1/2 |

## Corner cases and the single assumption

- **Trailing partial codon** (length not a multiple of 3 from the frame start): the 1–2 leftover
  bases are not a codon and are ignored (documented non-overlapping-triplet remainder rule).
- **Reading-frame offset** yields a different multiset for the same sequence — the MUST-test that
  distinguishes this unit (`ATGATGAAA` frame 1 → TGA = 1.0).
- **No valid codon** (all triplets ambiguous, or sequence shorter than 3): total = 0, so an **empty
  frequency table** is the only well-defined result — the guard avoids division by zero.

**ASSUMPTION (the single open one):** empty result for zero valid codons. Kazusa defines frequency
as count/total and excludes ambiguous codons but does not define the result when total = 0; returning
an empty table is the only value consistent with count/total and matches the implementation guard.
Non-correctness-affecting for any input with at least one valid codon.

**Contradictions:** none — Kazusa, EMBOSS `cusp`, Wikipedia and Nakamura 2000 agree on the
count/total (per-thousand) convention and ambiguous-base exclusion. The one nuance recorded is that
`cusp`'s *Fraction* column is the per-amino-acid RSCU-style metric, not the frequency under test.

## Related units

Same codon-usage family: the raw-count + TVD-comparison sibling
[[codon-usage-001-evidence|CODON-USAGE-001]] (concept [[codon-usage-comparison]]);
[[codon-rscu-001-evidence|CODON-RSCU-001]] / [[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]]
(the per-amino-acid `Fraction`/RSCU metric); [[codon-cai-001-evidence|CODON-CAI-001]] (CAI over
frequencies); [[codon-stats-001-evidence|CODON-STATS-001]] (the aggregation/reporting view);
[[codon-rare-001-evidence|CODON-RARE-001]] and [[codon-opt-001-evidence|CODON-OPT-001]].
