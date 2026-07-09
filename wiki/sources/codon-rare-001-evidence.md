---
type: source
title: "Evidence: CODON-RARE-001 (Rare codon detection — FindRareCodons + cluster methods)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-RARE-001-Evidence.md
sources:
  - docs/Evidence/CODON-RARE-001-Evidence.md
source_commit: 496eef96c5b3099b418a04b72a333aa2a248be2d
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-RARE-001

The validation-evidence artifact for test unit **CODON-RARE-001** (Rare Codon Detection —
`CodonOptimizer.FindRareCodons(string, CodonUsageTable, double)`, plus the 2026-06-24 addendum's
opt-in cluster methods `CalculateMinMaxProfile` and `FindRareCodonClusters`). One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm
itself is summarized in [[rare-codon-analysis]], the codon-usage family's thresholded-frequency +
cluster-detection unit whose sibling measures are [[relative-synonymous-codon-usage]] /
[[codon-adaptation-index]] / [[effective-number-of-codons]] and whose actuator is
[[codon-optimization]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Base method (`FindRareCodons`)** — a codon is **rare** when its usage frequency (within its
  synonymous family, from the target organism's `CodonUsageTable`) falls **below a threshold**;
  default **0.15**. Reports each rare codon's nucleotide position (0-indexed, multiple of 3), the
  translated amino acid, and the actual frequency.
- **Online sources (base)** — Wikipedia "Codon usage bias" (rare codons = low-abundance tRNAs,
  slower elongation, consecutive rare codons inhibit translation, folding/regulation significance);
  **GenScript GenRCA** (Fan et al. 2024, BMC Bioinformatics 25:309); **Kazusa Codon Usage Database**
  (35,799 organisms, the reference-frequency source); **Shu et al. 2006** (PMC6032470 — five
  consecutive rare CUA Leu codons cause ~3-fold translation inhibition in E. coli; 5′-end effect
  stronger than internal); **Sharp & Li 1987** (PMC340524 — rare codons have low relative
  adaptiveness, validating the frequency-threshold approach). Standard E. coli K12 rare codons
  (Kazusa MG1655, freq < 0.10): AGA 0.04, AGG 0.02, CGA 0.06, CUA 0.04, AUA 0.07, UAG(stop) 0.07.
- **Base edge cases** — empty seq → empty; length not ÷3 → trailing incomplete codons ignored; no
  rare / all rare handled; T→U converted internally; threshold 0 → all rare except freq-0 codons;
  threshold 1 → all rare; unknown codon → frequency 0, always reported.
- **Base invariants** — reported positions are multiples of 3 in `[0, len−3]`; codon length always 3;
  frequency ∈ [0, 1]; every reported codon has freq < threshold; amino acid matches the standard
  genetic code; deterministic.
- **Base limitations (documented)** — no codon **context** (consecutive rare codons), no positional
  weighting (5′ vs internal), single threshold for all amino acids, reports frequency not tRNA
  availability. The addendum below closes the "consecutive rare codons" limitation.

## Addendum (2026-06-24): cluster / run detection

Two published, citable cluster-detection methods, added as **opt-in** additions; the per-codon
default is unchanged.

- **%MinMax — Clarke & Clark 2008** ("Rare Codons Cluster", PLoS ONE 3(10):e3412; summed form
  cross-checked against Rodriguez et al. 2018, %MinMax, Protein Science). Sliding window over
  **per-amino-acid synonymous** quantities Xij (actual codon freq), Xmax,i / Xmin,i (most / least
  common synonymous codon), Xavg,i (family mean = Σfreq / n). Per window:
  `%Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i) × 100` when ΣXij > ΣXavg,i, else
  `%Min = Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i) × 100`. Each window yields one signed value
  (−100% = all rarest synonymous codons, +100% = all common, 0% = mean). Default window **18 codons**
  (Clarke & Clark 2008; the 2018 paper cites z = 17). Rare-codon **clusters appear as negative
  (%Min) peaks**.
- **Sherlocc — Chartier, Gaudreault & Najmanovich 2012** (Bioinformatics 28(11):1438–1445,
  DOI 10.1093/bioinformatics/bts149; reference impl `mtthchrtr/sherlocc`). A **7-codon-wide** window;
  a **'slow'/pause position** is one whose codon-usage frequency is **≤ threshold** (same per-codon
  rare criterion as `FindRareCodons`, default 0.15); a window with **≥ 4 of 7** slow positions is a
  **rare-codon cluster (RCC)**.

## Documented corner cases (cluster methods)

- **Single-codon amino acids (Met AUG, Trp UGG)** — Xmax = Xmin = Xavg = Xij, contributing 0 to both
  %MinMax numerator and denominator; the window value stays defined (no divide-by-zero / NaN).
- **Window longer than the sequence** — no windows / no clusters.
- **Isolated rare codons** — flagged by `FindRareCodons` but **not** a Sherlocc cluster unless ≥ 4
  fall in a 7-codon window (precisely the capability the addendum adds).
- **Overlapping qualifying windows** — merged into one maximal cluster so a long rare run is reported
  once (implementation choice; Sherlocc reports cluster regions, not raw windows).

## Worked oracles (E. coli K12 Arg family: CGU 0.38, CGC 0.40, CGA 0.06, CGG 0.10, AGA 0.04, AGG 0.02; Xavg = 1/6 = 0.16667, Xmax = 0.40, Xmin = 0.02)

- `AGA·AGA·AGA` (window 3, %MinMax) → **−86.3636…%** = −(0.16667−0.04)/(0.16667−0.02)×100.
- `CGC·CGC·CGC` (window 3, %MinMax) → **+100%** = (0.40−0.16667)/(0.40−0.16667)×100.
- `CUG·AGA` (window 2, %MinMax) → **+36.4706…%** = (0.54−0.33333)/((0.50−0.16667)+(0.40−0.16667))×100.
- `7× AGA` (window 7, ≥4) → **1 Sherlocc cluster**, codons 0–6, 7 rare.
- `3 AGA + 4 CGC` (window 7, ≥4) → **no cluster** (3 rare < 4).
- `4 AGA + 3 CGC` (window 7, ≥4) → **1 cluster**, exactly 4 rare (≥ 4).

Base-method test cases: `AUGAGA` → AGA at pos 3; `AUGAGAAGG` → AGA@3, AGG@6; `AUGCUG` → none;
`AGAAGGCGA` → all positions; `CUGCUA` → CUA@3 (CUG common, CUA rare).

## Contradictions / deviations

**None** — Wikipedia, GenScript/GenRCA, Kazusa, Shu 2006 and Sharp & Li 1987 agree that a
frequency threshold identifies rare codons; the two cluster methods (Clarke & Clark %MinMax,
Chartier Sherlocc) are complementary published algorithms retrieved and cross-checked this session,
each cited to peer-reviewed sources plus a reference implementation. The addendum's overlapping-window
merge is an explicitly flagged implementation choice (Sherlocc reports regions), not a departure from
any source rule. Additional references: Plotkin & Kudla 2011 (Nat. Rev. Genet. 12(1):32–42).
