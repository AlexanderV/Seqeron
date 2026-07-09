---
type: concept
title: "Aneuploidy detection (copy number from read depth)"
tags: [chromosome, algorithm]
sources:
  - docs/Evidence/CHROM-ANEU-001-Evidence.md
  - docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-aneu-001-evidence
      evidence: "Test Unit ID: CHROM-ANEU-001 ... Algorithm Group: Chromosome Analysis ... Aneuploidy Detection"
      confidence: high
      status: current
---

# Aneuploidy detection (copy number from read depth)

Detecting an **abnormal chromosome copy number** — the canonical aneuploidies of human cytogenetics
(monosomy, trisomy) — from sequencing **read depth**. This is the first ingested unit of the
**Chromosome-analysis** family and its copy-number/ploidy anchor; sibling chromosome-family units get
their own concepts — see [[centromere-analysis]] for the centromere / alpha-satellite anchor and
[[karyotype-analysis]] for descriptor-based karyotyping + ploidy detection and
[[synteny-and-rearrangement-detection]] for the synteny-block / rearrangement anchor and
[[telomere-analysis]] for the telomere repeat-detection / T/S-length anchor; arm-ratio follows.
Validated under test unit
**CHROM-ANEU-001**; the validation record is [[chrom-aneu-001-evidence]], and [[test-unit-registry]]
tracks the unit. See [[algorithm-validation-evidence]] for the artifact pattern.

Aneuploidy is "an abnormal number of chromosomes in a cell" (Wikipedia) — a human somatic cell with 45
or 47 instead of the usual 46. The algorithm has two stages: a per-bin **copy-number estimate** from
depth, then a **whole-chromosome classification** against the normal disomic (CN=2) state.

## Copy number from read depth

The estimate is a log-ratio of observed depth to the sample median:

```
logRatio   = log2(observedDepth / medianDepth)
copyNumber = round(2^logRatio × 2)          # clamped to [0, 10]
```

The `× 2` re-scales the ratio onto the diploid baseline, so a depth **ratio** of 1.0 maps to CN 2. The
anchor points:

| Depth ratio | log2 ratio | Estimated CN | Term |
|-------------|-----------|--------------|------|
| 0.0 | −∞ | 0 | Nullisomy |
| 0.5 | −1.0 | 1 | Monosomy |
| 1.0 | 0.0 | 2 | Disomy (normal) |
| 1.5 | 0.58 | 3 | Trisomy |
| 2.0 | 1.0 | 4 | Tetrasomy |

Depth is aggregated into **bins** by position / `binSize`; multiple chromosomes are grouped by name;
output bins are ordered.

## Whole-chromosome classification

A chromosome is called aneuploid only when a **dominant CN spans ≥ `minFraction` of its bins** (default
80%). The CN→term map follows standard cytogenetic terminology: 0→Nullisomy, 1→Monosomy, 2→Normal
(**no call** — only CN ≠ 2 is returned), 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, and CN > 5 → the generic
"Copy number = N". The `minFraction` gate is also how **mosaicism** (variable CN across cells) is
tolerated: a chromosome whose CN is inconsistent across bins fails the fraction threshold.

## Confidence

Each call carries a confidence ∈ [0, 1]:

```
confidence = 1 − min(1, |expected − observed|)
             expected = copyNumber / 2
             observed = 2^logRatio
```

At every integer-CN depth ratio the expected and observed values coincide, so confidence = 1.0 — the
`S1` boundary test verifies this for ratios 0.0/0.5/1.0/1.5/2.0.

## Documented clinical oracles

The classic autosomal and sex-chromosome aneuploidies are the source's worked examples: **Down**
syndrome (chr21 trisomy), **Edwards** (chr18 trisomy), **Patau** (chr13 trisomy), **Turner** (chrX
monosomy, 45,X), and **Klinefelter** (chrX, XXY / CN 3).

## Limitations

- **Sex chromosomes are not special-cased.** X/Y are scored against the CN=2 baseline like autosomes,
  so a normal male (single X) would be flagged monosomic. Documented as a limitation; it does not
  affect the autosome-focused detection the unit targets. A [[research-grade-limitations|research-grade]]
  simplification.
- **Partial aneuploidy** (sub-chromosomal CN changes from translocations) is detected at the per-bin
  level but does not trigger a whole-chromosome call, which needs consistent CN across ≥80% of bins.
- **Guards:** empty input → empty result; zero or negative median depth → empty (prevents division by
  zero).
