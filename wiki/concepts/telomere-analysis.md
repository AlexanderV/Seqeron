---
type: concept
title: "Telomere analysis (TTAGGG repeat detection + T/S ratio length)"
tags: [chromosome, algorithm]
mcp_tools:
  - analyze_telomeres
  - estimate_telomere_length_from_ts_ratio
sources:
  - docs/Evidence/CHROM-TELO-001-Evidence.md
  - docs/Validation/reports/CHROM-TELO-001.md
  - docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md
source_commit: 9dfe8fee4470a739dd91e9192efd5d7319ec5c50
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-telo-001-evidence
      evidence: "ID: CHROM-TELO-001 ... Title: Telomere Analysis ... Area: Chromosome"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-telo-001-report
      evidence: "Validation report CHROM-TELO-001 (Stage A/B PASS, CLEAN, 2026-06-24) validates the telomere-analysis unit; Area: Chromosome"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:aneuploidy-detection
      source: chrom-telo-001-evidence
      evidence: "Area: Chromosome — telomere analysis is a sibling of the chromosome-analysis copy-number/ploidy family anchored by aneuploidy-detection"
      confidence: medium
      status: current
---

# Telomere analysis (TTAGGG repeat detection + T/S ratio length)

Detecting the **telomere** — the repetitive cap at each chromosome end — and estimating its length. This
is the fifth ingested unit of the **Chromosome-analysis** family; sibling units get their own concepts —
see [[aneuploidy-detection]] (copy-number/ploidy anchor), [[centromere-analysis]] (centromere /
alpha-satellite), [[karyotype-analysis]] (descriptor karyotyping + ploidy), and
[[synteny-and-rearrangement-detection]] (synteny / rearrangement). Validated under test unit
**CHROM-TELO-001**; the literature-source record is [[chrom-telo-001-evidence]] and the two-stage
verdict (Stage A/B PASS, ✅ CLEAN, no code changed) is [[chrom-telo-001-report]], with
[[test-unit-registry]] tracking the unit. See [[algorithm-validation-evidence]] for the artifact
pattern.

A telomere is a run of the conserved vertebrate hexamer `TTAGGG` (Meyne et al. 1989; conserved across
all vertebrates). It shortens ~50–100 bp per cell division; once **critically short** it triggers a DNA
damage response and cellular senescence (Wikipedia). The unit has two independent parts: **repeat
detection** at each chromosome end, and **length estimation** — directly from the detected run or from a
qPCR **T/S ratio**.

## Repeat detection and orientation

The two chromosome ends carry the repeat in opposite orientation:

| End | Repeat searched | Direction |
|-----|-----------------|-----------|
| 3′ (terminus) | `TTAGGG` | forward, from the sequence end inward |
| 5′ (start) | `CCCTAA` (reverse complement of `TTAGGG`) | from the sequence start inward |

Each unit is 6 bp. Detection reports, per end, whether a telomere is present, its length, and its
**purity**.

## Purity (imperfect-repeat tolerance)

Biological telomeres diverge from a perfect repeat, so each window is accepted under a configurable
**70 % similarity** threshold rather than exact match:

- 6 bp unit (`TTAGGG`): **5/6** bases must match — 1 mismatch allowed.
- 7 bp unit (e.g. Arabidopsis `TTTAGGG`): **5/7** must match — 2 mismatches allowed.

Purity ∈ [0,1] is tracked alongside length; e.g. a divergent `TTAGGA` run scores 5/6 ≈ **0.833**. Higher
purity indicates a younger/healthier telomere.

## Length estimation — direct and by T/S ratio

The detected run gives length directly (200 × `TTAGGG` on 1000 A's → length 1200, purity 1.0). The qPCR
method (Cawthon 2002) estimates average length from the **T/S ratio** — telomere repeat copy number (T)
over a single-copy reference gene (S), proportional to length:

```
EstimatedLength = referenceLength × (tsRatio / referenceRatio)
```

| T/S ratio | referenceRatio | referenceLength | Estimated length |
|-----------|----------------|-----------------|------------------|
| 1.0 | 1.0 | 7000 | 7000 |
| 1.5 | 1.0 | 7000 | 10500 |
| 0.5 | 1.0 | 7000 | 3500 |
| 2.0 | 1.0 | 7000 | 14000 |
| 1.0 | 2.0 | 7000 | 3500 |
| 0.0 | 1.0 | 7000 | 0 |

## Invariants

- **Length** ≥ 0; **purity** ∈ [0,1].
- `Has5Prime/3Prime = true` ⇒ length ≥ `minTelomereLength` (threshold consistency).
- `IsCriticallyShort = (hasTelomere && length < criticalLength) OR empty`.
- **T/S linearity:** `EstimatedLength = referenceLength × (tsRatio / referenceRatio)`.
- **Orientation:** 5′ expects the reverse complement, 3′ the forward repeat.

## Configurable parameters (defaults, not biological constants)

- `criticalLength` = 3,000 bp — below this a detected telomere is flagged critically short.
- `minTelomereLength` = 500 bp — detection sensitivity threshold.
- `searchLength` — search-window cap; a shorter window **truncates** the reported length (search window
  600 on a 1200 bp run → length 600).
- `referenceLength` = 7,000 bp — reference for T/S length calculations.

All four are engineering defaults, explicitly **not** fixed biological constants (Cawthon's clinical
values differ by assay). Empty sequence → no telomere, marked critically short.

## Scope

The default detector targets the vertebrate `TTAGGG`/`CCCTAA` pair. The artifact documents cross-species
variation (Arabidopsis `TTTAGGG`, Tetrahymena `TTGGGG`, S. cerevisiae `TGTGGGTGTGGTG`, Bombyx `TTAGG`);
a non-vertebrate telomere requires its own repeat unit and per-window match count. A
[[research-grade-limitations|research-grade]] simplification: repeat-run detection, not the full
shelterin/t-loop biology.
