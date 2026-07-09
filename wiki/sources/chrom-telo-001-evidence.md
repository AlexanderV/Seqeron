---
type: source
title: "Evidence: CHROM-TELO-001 (Telomere analysis — TTAGGG repeat detection + T/S ratio length)"
tags: [validation, chromosome]
doc_path: docs/Evidence/CHROM-TELO-001-Evidence.md
sources:
  - docs/Evidence/CHROM-TELO-001-Evidence.md
source_commit: 3b01c634a9332ecb839b0da26ecc154feebc1d56
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CHROM-TELO-001

The validation-evidence artifact for test unit **CHROM-TELO-001** — telomere analysis: detecting the
canonical vertebrate telomere repeat at each chromosome end and estimating telomere length (directly
from the detected repeat run, and from a qPCR T/S ratio). This is the fifth **Chromosome-analysis**
family Evidence file (after [[chrom-aneu-001-evidence]], [[chrom-cent-001-evidence]],
[[chrom-karyo-001-evidence]], [[chrom-synt-001-evidence]]) and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the repeat-detection rule,
purity model, T/S linearity and invariants are synthesized in [[telomere-analysis]], the anchor for the
chromosome telomere family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online / literature sources:**
  - **Wikipedia — "Telomere"** (encyclopedia) — the vertebrate `TTAGGG` repeat, per-division shortening
    (~50–100 bp), critically short telomeres → DNA-damage response + senescence, the 75–300 base 3′
    single-stranded `TTAGGG` overhang, and the cross-species repeat table.
  - **Meyne, Ratliff & Moyzis (1989)**, *PNAS* 86(18):7049 (PMID 2780561) — conservation of the human
    telomere sequence `(TTAGGG)n` across all vertebrates (the repeat-sequence authority).
  - **Cawthon (2002)**, *Nucleic Acids Res* 30(10):e47 (PMID 12000852) — telomere measurement by
    quantitative PCR; the **T/S ratio** (telomere copy number T over single-copy gene S) is proportional
    to average telomere length (r² = 0.677 vs Southern-blot TRF).
  - **Blackburn & Gall (1978)**, *J Mol Biol* (PMID 642006) — Nobel-recognised foundational telomere
    structure work.

- **Algorithm behaviour (from the implementation):**
  - **Repeat orientation** — 3′ end (chromosome terminus) carries forward `TTAGGG` repeats extending
    toward the end; 5′ end carries `CCCTAA` (the reverse complement of `TTAGGG`). Detection searches for
    `CCCTAA` from the start and `TTAGGG` from the end; each repeat unit is 6 bp.
  - **Per-window purity** — biological telomeres diverge from perfect repeats, so detection uses a
    configurable **70 % per-window similarity** threshold: for the 6 bp unit, 5/6 bases must match (1
    mismatch allowed); for a 7 bp unit (e.g. Arabidopsis `TTTAGGG`), 5/7 must match. Purity ∈ [0,1] is
    tracked alongside length; higher purity = younger/healthier telomere.
  - **Length by T/S ratio** — `EstimatedLength = referenceLength × (tsRatio / referenceRatio)` (linear,
    Cawthon 2002); a T/S ratio of 1.0 at referenceRatio 1.0 returns the reference length.
  - **Critical-length flag** — `IsCriticallyShort = (hasTelomere && length < criticalLength) OR empty`;
    the default critical threshold **3,000 bp** is a configurable implementation default, *not* a fixed
    biological constant.
  - **Configurable parameters (implementation defaults, not biological constants):** `criticalLength`
    3,000 bp, `minTelomereLength` (detection threshold) 500 bp, `searchLength` window (truncates
    reported length), and `referenceLength` 7,000 bp for T/S calculations.

- **Datasets / oracles:** hand-built repeat runs — `[1000 A's]+[200×TTAGGG]` → Has3PrimeTelomere,
  length 1200, purity 1.0; `[200×CCCTAA]+[1000 A's]` → Has5PrimeTelomere; both-ends case; no-telomere /
  empty (critically short); short-telomere gated by `minTelomereLength`; divergent `TTAGGA×200` →
  purity 5/6 ≈ 0.833; long `TTAGGG×2000` → length 12000; `searchLength`-limited truncation. T/S table:
  {1.0,1.5,0.5,2.0}@ref 7000 → {7000,10500,3500,14000}; referenceRatio 2.0 → 3500; T/S 0.0 → 0.

## Assumptions / limitations (from the artifact)

- **Deviations and Assumptions: None** — all repeat sequences (Wikipedia table / Meyne 1989), the T/S
  proportionality (Cawthon 2002), and 5′/3′ orientation (Wikipedia chromosome structure) are verified
  against the cited sources.
- The four configurable parameters (`criticalLength`, `minTelomereLength`, `searchLength`,
  `referenceLength`) are **implementation defaults**, explicitly flagged as engineering choices rather
  than biological constants.
- Species variation is documented (human/mouse/Xenopus `TTAGGG`, Arabidopsis `TTTAGGG`, Tetrahymena
  `TTGGGG`, S. cerevisiae `TGTGGGTGTGGTG`, Bombyx `TTAGG`), but the default detector targets the
  vertebrate `TTAGGG`/`CCCTAA` pair; non-vertebrate units require a different repeat with its own
  per-window match count.

No contradictions among the sources — the encyclopedic article, Meyne 1989 repeat conservation, and
Cawthon 2002 T/S proportionality agree; the `TTAGGG` unit and T/S linearity recur consistently.
