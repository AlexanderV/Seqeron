---
type: source
title: "Evidence: SEQ-DINUC-001 (Dinucleotide analysis — frequencies, relative abundance, codon frequencies)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-DINUC-001-Evidence.md
sources:
  - docs/Evidence/SEQ-DINUC-001-Evidence.md
source_commit: 4a7f3b50df393c2ccf0fe505da489d087d4f22f4
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-DINUC-001

The validation-evidence artifact for test unit **SEQ-DINUC-001** — **Dinucleotide Analysis**
(dinucleotide frequencies, observed/expected relative abundance, and codon frequencies), method
`CalculateDinucleotideRatios` and its frequency companion in the Analysis assembly. One instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; see
[[test-unit-registry]] for how units are tracked. The synthesized write-up — definitions, the
Karlin odds ratio, the CpG-O/E relationship, oracles, and edge cases — lives on
[[dinucleotide-relative-abundance]]; this file records the source trace and datasets.

## What this file records

- **Online sources:**
  - **Karlin S., "Pervasive properties of the genomic signature" (PMC126251)** (rank 1) — the
    odds-ratio / relative-abundance definition `ρ_XY = f_XY / (f_X · f_Y)`, where `f_x` is the
    normalized single-base frequency and `f_xy` the normalized dinucleotide frequency on the leading
    strand. `r_xy = 1.0` ⇒ no bias (independence baseline); `>1` over-represented, `<1`
    under-represented. Frequencies are **normalized** (counts / respective totals), not raw counts.
  - **Tsirigos & Rigoutsos / Karlin & Burge — MBE 19(6):964** (rank 1) — restates
    `ρ*_XY = f*_XY / (f*_X f*_Y)` and gives the **interpretation thresholds** attributed to Karlin &
    Burge (1995): under-represented if `ρ ≤ 0.78`, over-represented if `ρ ≥ 1.23`. Documentation-only —
    the method returns the raw ratio, not the classification.
  - **Gardiner-Garden & Frommer (1987) CpG O/E** (rank 1/4) — `O/E = (#CpG/N)/((#C/N)·(#G/N))`, the
    same odds-ratio shape as Karlin's ρ but normalized by `N` rather than `N − 1`. Recorded as the
    `N/(N−1)` modeling-choice difference behind [[cpg-island-detection]], not an error.
  - **Kazusa Codon Usage Database (CUTG) readme** (rank 5/2) — codon frequency = count / total
    counted codons over **non-overlapping** in-frame triplets; ambiguous (non-ACGT) triplets excluded.

- **Datasets (hand-derived exact rationals, sequence `ATGCGCGT`, length 8):**
  mono counts A=1,T=2,G=3,C=2; dinucleotide counts (7 positions) AT=1,TG=1,GC=2,CG=2,GT=1.

  | Base step | `f_XY` | `ρ_XY` |
  |-----------|--------|--------|
  | GC | 2/7 | 64/21 = 3.047619047619048 |
  | CG | 2/7 | 64/21 = 3.047619047619048 |
  | AT | 1/7 | 32/7 = 4.571428571428571 |
  | TG | 1/7 | 32/21 = 1.523809523809524 |
  | GT | 1/7 | 32/21 = 1.523809523809524 |

  Codon frequencies: `ATGATGAAA` frame 0 → ATG=2/3, AAA=1/3; frame 1 → TGA=1.0.

- **Corner cases:** independence baseline `f_X·f_Y` is undefined when a constituent base is absent
  (`f_X = 0` ⇒ expected 0 ⇒ division by zero) — the method **returns ratio 0** (expected-0 guard);
  codon counting starts at the frame offset and steps by 3, trailing 1–2 bases ignored, non-ACGT
  triplets excluded.

- **Recommended coverage:** the ρ oracles above; dinucleotide-frequency values/sum semantics; codon
  frequencies for frame 0 and 1; ρ ≈ 1.0 for a uniform/independent sequence (no-bias baseline);
  null/empty/length-<2 (ratios, freqs) and length-<3 (codons) → empty; division-by-zero guard →
  ratio 0; lowercase and RNA (`U`) inputs.

## Assumptions (both documented modeling choices, not correctness gaps)

1. **Karlin `(N−1)` normalization vs Gardiner-Garden `(N)`** — dinucleotide frequency is
   `count/(N−1)` (Karlin odds-ratio convention) rather than `count/N` (the CpG convention). Both are
   published and give the same odds-ratio shape, differing only by `N/(N−1) → 1` for long sequences.
2. **`U` treated as a fifth base** — `CalculateDinucleotideRatios`' single-base denominator includes
   `U` for RNA inputs, matching `CalculateNucleotideComposition`; RNA dinucleotide signatures are
   defined the same way with `U` replacing `T`.

## Contradictions

**None.** Karlin (PMC126251), Karlin & Burge 1995 (via MBE 19(6):964), Gardiner-Garden & Frommer
1987, and the Kazusa CUTG readme agree on the odds-ratio definition, the independence baseline, and
the codon count/total convention. The two items above are documented normalization/alphabet choices,
each authoritative in its own literature.

## Related units

- **[[dinucleotide-relative-abundance]]** — the synthesizing concept (dinucleotide frequency + Karlin
  genomic signature; CpG O/E as the `CG` special case).
- **[[cpg-island-detection]]** (EPIGEN-CPG-001) — the `CG`-specialized odds ratio + island calling.
- **[[base-composition]]** (SEQ-COMPOSITION-001) — the single-base `f_X` layer this builds on.
- **[[seq-codon-freq-001-evidence]]** (SEQ-CODON-FREQ-001) / [[codon-usage-comparison]] — the codon
  frequency metric this source's codon output mirrors.
- **[[nucleotide-composition-skew]]** — the single-base compositional-asymmetry cousin.
