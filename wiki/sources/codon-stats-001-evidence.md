---
type: source
title: "Evidence: CODON-STATS-001 (Codon Usage Statistics — GetStatistics + CAI)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-STATS-001-Evidence.md
sources:
  - docs/Evidence/CODON-STATS-001-Evidence.md
source_commit: 1256bb39beca90dab01e639a9bb8bbc4229010b1
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-STATS-001

The validation-evidence artifact for test unit **CODON-STATS-001** (Codon Usage Statistics —
`GetStatistics` + `CalculateCai`). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. This unit is the **aggregation /
reporting view** of the codon-usage family: `GetStatistics` bundles the family's individual
measures over one input plus a positional-GC composition block, and `CalculateCai` is the same
CAI computation validated on its own as [[codon-cai-001-evidence|CODON-CAI-001]]. It therefore
reuses rather than redefines the family concepts — see [[relative-synonymous-codon-usage]],
[[codon-adaptation-index]], and [[effective-number-of-codons]]; the CDS-rewriting sibling is
[[codon-optimization]] and the diagnostic sibling is [[rare-codon-analysis]]. See
[[test-unit-registry]] for how units are tracked.

## What GetStatistics aggregates

- **Codon counts** and **total codons** — the non-overlapping-triplet tally (same counting step
  behind [[relative-synonymous-codon-usage|RSCU]]'s `CountCodons`).
- **RSCU** per codon — [[relative-synonymous-codon-usage]] (`RSCU = n_i·x/Σx`, family-mean
  normalization, sense codons only).
- **ENC / Nc** — [[effective-number-of-codons]] (Wright 1990 reference-free bias measure in
  [20, 61]).
- **CAI** (`CalculateCai`) — [[codon-adaptation-index]] (geometric mean of `w = f/max_synonym_f`
  against an E. coli or human reference table).
- **Positional GC composition** — `GC1 / GC2 / GC3`, `GC3s`, and `OverallGc` (documented below).
- **Reference tables** — `EColiOptimalCodons` (Sharp & Li 1987 `w` values via Biopython
  `SharpEcoliIndex`) and `HumanOptimalCodons` (RSCU derived from Kazusa *H. sapiens* [gbpri]).

## Positional GC composition (the distinctive content)

This is the one measure not covered by an existing family concept, so it is recorded here.

- **GC1 / GC2 / GC3** — the fraction of G or C at codon positions 1, 2, 3 respectively, over
  **all** codons (EMBOSS `cusp` "1st/2nd/3rd letter GC"). `OverallGc` = (GC1+GC2+GC3)/3 — a
  derived convenience field.
- **GC3s** — GC content of the **third position of *synonymous* codons only**, i.e. **excluding
  Met (ATG), Trp (TGG) and the three stop codons** (Peden 1999 thesis §1.8.2.1.3; confirmed by
  the PMC7596632 "59 synonymous codons" restatement). GC3s ≠ GC3 precisely because GC3s drops
  the non-degenerate Met/Trp/stop third positions — the worked case `ATGGCA` gives **GC3s = 0**
  (Met excluded, Ala GCA third base = A) while **GC3 = 50** (ATG third = G, GCA third = A).
- **Empty synonymous denominator** — a sequence with only Met/Trp/stop has no synonymous third
  positions → GC3s reported as **0** (Peden empty-denominator corner case).

## Datasets (oracles)

- **E. coli `w` table** (Sharp & Li 1987 via Biopython `SharpEcoliIndex`): CTG=1.000, CTC=0.037,
  GCT=1.000, GCC=0.122, GCA=0.586, CGT=1.000, TTC=1.000, TTT=0.296 — reproduced by
  `EColiOptimalCodons`.
- **Human RSCU** (Kazusa *H. sapiens*, `RSCU = n·x/Σx`): CTG≈2.3713, GCC≈1.5988, ATG=1.0000,
  TGG=1.0000 — reproduced by `HumanOptimalCodons`.
- **Worked CAI / GC cases**: `CTGATCGTTGCTCGTAAA`→CAI 1.0 (all `w=1`); `GCTGCC`→CAI
  0.34928…=√(1×0.122); `CTAATAGTC`→CAI 0.01114…; `ATGTGGTAA`→CAI **0** (Met+Trp+stop only, no
  scorable codon); `GCCGCA`→GC3s 50.0; `ATGGCA`→GC3s 0.0 / GC3 50.0; `CTGGTTAAA`→GC1/GC2/GC3 =
  66.67/0.0/33.33.

## Assumptions and deviations (from the artifact)

1. **GC3s reported as a percentage (×100).** CodonW reports GC3s as a fraction in [0,1]; this
   implementation reports it as a percentage for consistency with the GC1/GC2/GC3 fields (EMBOSS
   `cusp` percentage style). A unit/labeling choice only — the synonymous-codon subset in the
   numerator/denominator is exactly per Peden.
2. **Zero-`w` codons are skipped, not floored to 0.01.** seqinr/EMBOSS substitute a small value
   (0.01, Bulmer 1988) to avoid `ln(0)`; this implementation instead skips codons whose relative
   adaptiveness is 0 (an entirely-zero gene → CAI 0). For the supplied reference tables no
   synonymous codon has `w=0`, so real-CDS CAI is unaffected; only a gene using a codon entirely
   absent from the reference differs. (Note: this is the CAI zero-handling recorded for this
   *stats* unit; the standalone [[codon-adaptation-index|CODON-CAI-001]] documents a `1e-6` clamp
   for the "codon absent but family present" case — the two artifacts describe the same guard from
   different angles.)
3. **Single-codon / stop exclusion** for CAI and GC3s follows seqinr / CodonW: Met, Trp and stop
   codons are excluded. A Met/Trp/stop-only sequence has no scorable codon → CAI 0, GC3s 0.

**Contradictions:** none between sources — Sharp & Li 1987 (+ Biopython reproduction), Wikipedia,
seqinr, CodonW/Peden, EMBOSS `cusp`, and Kazusa agree on the formulae and the synonymous-codon
exclusion set. Deviations are unit-labeling and zero-frequency-handling choices, both documented.
