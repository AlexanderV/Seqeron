---
type: source
title: "Evidence: MIRNA-PRECURSOR-001 (pre-miRNA precursor hairpin detection)"
tags: [validation, mirna]
doc_path: docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md
sources:
  - docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md
source_commit: e0541d580467f016d02636dab852b866b6e05940
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MIRNA-PRECURSOR-001

The validation-evidence artifact for test unit **MIRNA-PRECURSOR-001** — **Pre-miRNA Hairpin
Detection** (`MiRnaAnalyzer.FindPreMiRnaHairpins` and its opt-in MFE-fold / cleavage-ruler / trained-
classifier extensions). This is the **second ingested unit of the MiRNA family** and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
synthesizing concept is [[pre-mirna-hairpin-detection]]. [[test-unit-registry]] tracks the unit.

The unit builds on the shared [[rna-base-pairing]] primitive (the stem is scored as consecutive
{A-U, G-C} + G-U wobble pairs).

## What this file records

- **Online sources (rank 1 unless noted):** Bartel 2004 (*Cell* review — ~70 nt stem-loop, ~33 bp
  imperfect stem, 4-15 nt loop, Drosha cuts ~11 nt from base), Ambros 2003 (*RNA* uniform-annotation
  criteria — mature ~22 nt in one hairpin arm; MFE-path stem ≥16 complementary bases), Bartel 2009,
  Krol 2004 (stem 22-33 bp critical, G:U wobbles/bulges common, ≥18-22 bp effective stem),
  Wikipedia MicroRNA (rank 4 — biogenesis, 5p/3p arms, 2-nt overhang, strand selection), miRBase
  (rank 5 — 55-120 nt precursor range, hsa-mir-21 72 nt / hsa-let-7a-1 80 nt sequences).
  MFE-path & later sections add Bonnet 2004 (pre-miRNAs fold lower than shuffled), Zhang 2006
  (AMFE = 100·|MFE|/len, MFEI = AMFE/(G+C)%, pre-miRNA MFEI > 0.85), Meyers 2008 (single-arm duplex,
  minimal bulges), Han 2006 (Drosha ~11 bp from basal junction), Park 2011 (Dicer ~22 nt 5'-counting
  rule), Lee 2003 (RNase III 2-nt 3' overhang), Auyeung 2013 (CNNC 16-18 nt from Drosha cut),
  Altschul-Erickson 1985 (dinucleotide-preserving shuffle), and Turner 2004 (NNDB parameters).

- **Documented corner cases:** sequences < ~55 nt cannot form a valid precursor; no
  self-complementarity → no stem; loop < 3 nt or > ~25 nt; stem < 18 bp; DNA input needs T→U;
  G:U wobbles count as valid stem pairs; long sequences may yield overlapping candidates.

- **Datasets / oracles:**
  - **Dataset 1** — Watson-Crick hairpin (57 nt, 25+7+25); heuristic effective stem 23 bp (capped
    `n/2 − 5`); MFE fold ΔG° −48.48, 27 stem bp, loop 3, MFEI 1.9392 → ACCEPT.
  - **Dataset 2** — same hairpin with 4 G:U wobble pairs (Krol 2004) — wobbles counted in the stem.
  - **Dataset 3** — hsa-mir-21 (72 nt): heuristic **NOT detected** (only 16 consecutive end-pairs);
    MFE fold **ACCEPT** (ΔG° −35.13, 32 bp, MFEI 1.0037).
  - **Dataset 4** — hsa-let-7a-1 (80 nt): heuristic **NOT detected** (5 consecutive pairs); MFE fold
    **ACCEPT** (ΔG° −34.31, MFEI 1.0091).
  - **Dataset 5** — short-stem hairpin (55 nt, stem 15 bp): **rejected** (15 < 18), isolating the
    stem-length gate from the n < 55 early exit.
  - **Structure-over-ΔG guard** — a 120-nt multibranch `5S-rRNA-like` fold (ΔG° −47.04, stronger than
    Dataset 1) is **REJECTED** because it is not a single dominant hairpin.
  - **Cleavage cross-check** — hsa-miR-21-5p: basal junction + 11 bp Drosha cut + 22-nt Dicer span
    reproduces `UAGCUUAUCAGACUGAUGUUGA` (22 nt) exactly against miRBase.
  - **Classifier spot-check** — hsa-mir-21 P(natural) ≈ 0.99999, hsa-let-7a-1 ≈ 0.99989, a
    di-shuffled hsa-mir-21 ≈ 0.00052; held-out accuracy = AUC = 1.0.

- **Test-coverage recommendations:** MUST — valid hairpin detection, short-sequence/null-empty
  returns empty, stem ≥18 bp, loop 3-25 nt, T→U conversion, mature/star arm extraction, dot-bracket
  correctness, free-energy relative ordering, real miRBase NOT-detected (heuristic limitation).
  SHOULD — multiple hairpins, G:U wobble accepted, `maxHairpinLength` respected. COULD — O(n²)
  performance.

## Deviations and assumptions

Two recorded items on the **default** path, both accepted and mitigated by the opt-in MFE fold:

1. **ASM-03 / Deviation 1** — mature strand is extracted from the fixed 5' arm and star from the
   mirrored 3' arm; precursors whose dominant mature product comes from the opposite arm are
   represented asymmetrically.
2. **ASM-01 / Deviation 2** — the default requires **uninterrupted end-to-end stem pairing**, which
   is stricter than natural pre-miRNA folding, so real miRBase precursors with internal bulges/
   mismatches are frequently rejected (documented via hsa-mir-21 / hsa-let-7a-1, tests M18/M19). The
   opt-in `FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe` remove this by folding with the
   RNA-STRUCT-001 Zuker–Stiegler engine and reading the hairpin from the real MFE structure.

**No source contradictions.** The consecutive-pairing default and the opt-in MFE/cleavage/classifier
paths are consistent with all cited literature; the trained classifier uses only public-domain
miRBase positives and di-shuffled negatives (no GPL miRDeep2 code/weights). The only residual gap is
the read-stacking (small-RNA-seq) miRDeep2 signal, which is data-blocked (needs the caller's reads).
