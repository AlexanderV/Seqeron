---
type: source
title: "Evidence: CODON-ENC-001 (Effective Number of Codons — ENC / Nc)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-ENC-001-Evidence.md
sources:
  - docs/Evidence/CODON-ENC-001-Evidence.md
source_commit: 5bc4ea5003342f5c5c657d68183edfae40fba29a
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-ENC-001

The validation-evidence artifact for test unit **CODON-ENC-001** (Effective Number of Codons —
ENC / Nc, Wright 1990). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is summarized in
[[effective-number-of-codons]], a reference-free sibling of the
[[codon-adaptation-index]] in the [[relative-synonymous-codon-usage]] codon-usage family. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** — **Fuglsang (2004)** *"The 'effective number of codons' revisited"* (Biochem.
  Biophys. Res. Commun. 317:957–964, retrieved PDF, rank 1) which reproduces Wright's equations
  **verbatim**: codon homozygosity Eq. 1 `F̂ = (n·Σp_i² − 1)/(n − 1)`, per-aa `N̂c = 1/F̂` (Eq. 2),
  gene aggregation `N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆` (Eq. 3), within-class averaging for
  missing amino acids (Eq. 4), the isoleucine `F̂₃ = (F̂₂+F̂₄)/2` fallback (Eq. 5a), the re-adjust-
  down-to-61 rule, and the 20–61 range. **Fuglsang (2006)** *"Estimating the ENC: the Wright way…"*
  (Genetics 172(2):1301–1307, rank 1) confirms the sampling-without-replacement `F` and the
  constraints. Standard-genetic-code degeneracy partition (NCBI, rank 2): 5 quartets / 9 doublets /
  3 sextets / 1 triplet (Ile) / 2 singlets (Met, Trp) — exactly the constants in Eq. 3.
  Primary reference: **Wright, F. (1990)** *Gene* 87(1):23–29.
- **Algorithm spec** — reciprocal-homozygosity effective-codon count aggregated by degeneracy class;
  reference-free; ranges 20 (max bias, one codon per aa) to 61 (no bias, uniform usage). Calculable
  only when each represented amino acid has ≥ 2 codons (the `n−1` denominator).
- **Datasets** — (1) fully unbiased gene → asymptotic `F̂₂=0.5, F̂₃=1/3, F̂₄=0.25, F̂₆=1/6` summing to
  **61**; (2) Fuglsang (2004) "no-bias-discrepancy" simulation → **40.5** (per-aa Nc 1.5/2.0/2.5/3.5
  by class + 1.0 singlets); (3) hand-derived 2-fold Phe TTT×3,TTC×1 ⇒ `F̂=0.5` ⇒ `Nc(Phe)=2`, and
  even split ⇒ `Nc(Phe)=3` (per-aa overshoot).
- **Corner cases / failure modes** — amino acid with `n ≤ 1` (F̂ undefined → within-class average
  substitute); empty 3-fold class (Ile absent → `F̂₃=(F̂₂+F̂₄)/2`); Eq. 3 overshoot past 61
  (re-adjust to 61); very low F̂ on small counts (per-aa Nc exceeds degeneracy, still capped at 61).

## The one assumption (from the artifact)

1. **Lower clamp at 20.** Wright/Fuglsang prescribe re-adjusting **down to 61** at the top but give
   no explicit hard clamp at the bottom; 20 is the structural minimum. Retaining `Math.Max(20, …)`
   is consistent with the stated range and cannot lower a legitimately-computed value — a defensive
   bound, not an algorithmic parameter.

## Recommended coverage (from the artifact)

Fully-unbiased → 61; single-aa two-fold with hand-derived `F̂`; maximally-biased → 20; invariant
`20 ≤ Nc ≤ 61`; isoleucine-absent uses the `(F̂₂+F̂₄)/2` fallback; null → ArgumentNullException,
empty/whitespace → 0; lowercase-normalized / invalid-codon-skip; overshoot re-adjust to 61.

**Contradictions:** none — Fuglsang (2004) and (2006) reproduce the same Wright equations, and the
NCBI degeneracy partition matches the Eq. 3 constants exactly. The only recorded deviation is the
non-source-prescribed lower clamp at 20 (a harmless defensive bound). Flagged on
[[effective-number-of-codons]].
