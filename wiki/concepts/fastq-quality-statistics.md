---
type: concept
title: "FASTQ quality statistics (Q20/Q30, mean, median, min/max, std dev)"
tags: [file-io, algorithm, validation]
mcp_tools:
  - fastq_filter
  - fastq_statistics
sources:
  - docs/Evidence/QUALITY-STATS-001-Evidence.md
  - docs/algorithms/Quality/Quality_Statistics.md
source_commit: b0582a69a05e98839d57529bf7d538af4b125f29
created: 2026-07-10
updated: 2026-07-16
graph:
  relationships:
    - predicate: depends_on
      object: concept:phred-quality-encoding
      source: quality-stats-001-evidence
      evidence: "Cock et al. (2010) Phred+33 decode (ord−33) ... Used here only to convert the quality string to scores prior to statistics; the decode itself is the contract of QUALITY-PHRED-001 and is not re-tested."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: quality-stats-001-evidence
      evidence: "Test Unit ID: QUALITY-STATS-001 ... Algorithm: Quality Statistics (Phred score summary statistics; Q20/Q30 percentages)"
      confidence: high
      status: current
---

# FASTQ quality statistics (Q20/Q30, mean, median, min/max, std dev)

The **quality-statistics** surface summarizes the per-base Phred scores of a FASTQ read into
QC scalars: the **Q20 / Q30 fractions**, **mean**, **median**, **min / max**, and **population
standard deviation** (with its variance). It is the run-quality report over a decoded quality
string — distinct from, and downstream of, the [[phred-quality-encoding]] concept that turns the
ASCII quality line into scores. That encoding concept explicitly flagged the statistics as a
**separate surface**; this page is that surface. The dedicated unit
[[quality-stats-001-evidence]] anchors it; [[test-unit-registry]] tracks the owning unit.

## Inputs: decoded Phred scores, not characters

The statistics are computed over the **decoded Phred scores**, so the ASCII **offset is
irrelevant** to every result — a Phred+33 string and the Phred+64 string with the same underlying
scores yield identical statistics. The decode step (`ord(c) − 33` for Phred+33, `ord(c) − 64` for
Phred+64) is the contract of the sibling **Phred-handling** unit [[quality-phred-001-evidence]] and
the shared [[phred-quality-encoding]] concept; it is **consumed here, not re-tested**. The full
quality string is treated as a **complete population** of observed scores, which fixes the
variance/σ convention below.

## The statistics

| Statistic | Definition | Source |
|-----------|-----------|--------|
| **Mean** | arithmetic mean of the decoded scores `μ = (Σ Qᵢ)/N` | descriptive stats |
| **Median** | middle order statistic; **even N** averages the two central values | Math is Fun |
| **Min / Max** | smallest / largest decoded score | — |
| **Population variance** | `σ² = (1/N) Σ (Qᵢ − μ)²` — divide by **N** (not N−1) | Newcastle Univ. ASK |
| **Population std dev** | `σ = √σ²` | Newcastle Univ. ASK |
| **% ≥ Q20** | fraction of bases with score **≥ 20**, ×100 | Illumina |
| **% ≥ Q30** | fraction of bases with score **≥ 30**, ×100 | Illumina |

Two conventions are load-bearing and easy to get wrong:

- **Population, not sample:** variance/σ divide by **N** (Bessel's correction is *not* applied),
  because the quality string is the complete set of observed scores, not a sample of a larger set.
- **Inclusive Q thresholds:** `% ≥ Q30` counts a base **exactly at Q30** (1-in-1000 error);
  same for Q20 (`≥ 20`). The thresholds are `≥`, not `>`.

**Mean is arithmetic over log-scaled scores.** The mean here averages the Phred *scores*
directly; it is **not** an error-probability-averaged quality (which would average `p = 10^(−Q/10)`
and re-encode). Both are "mean quality" in the literature — this unit's contract is the arithmetic
score mean.

## Q30 as the NGS benchmark

Per the Illumina spec, **`% ≥ Q30`** is the standard run-quality metric — "Q30 is considered a
benchmark for quality in next-generation sequencing"; at Q30 "virtually all of the reads will be
perfect." The convenience method `CalculateQ30Percentage(...)` must equal
`CalculateStatistics(...).PercentAboveQ30` for the same input (one is a shortcut for the other).

## Worked oracle

Phred+33 string `5?I` decodes to scores **20 / 30 / 40**:

- Mean `(20+30+40)/3 = 30.0`; median (odd, n=3) `30`; min/max `20 / 40`.
- Population variance `((−10)² + 0² + 10²)/3 = 200/3 ≈ 66.6667`; σ `= √(200/3) ≈ 8.16497`.
- `% ≥ Q20 = 3/3 = 100.0`; `% ≥ Q30 = 2/3 ≈ 66.667` (the base at Q20 is excluded from Q30).

Even-length `5II?` → scores `20,40,40,30`, sorted `20,30,40,40` → **median `(30+40)/2 = 35.0`**,
mean `32.5` (the even-count median averages the two central order statistics). Single base `I`
(Q40) → mean = median = min = max = `40`, **σ = 0** (zero spread), `% ≥ Q20 = % ≥ Q30 = 100.0`.

## Corner cases

- **Even vs odd count:** the median branch differs — average the two central order statistics for
  even N, take the single middle for odd N.
- **Single element:** mean = median = min = max = that value; population σ = 0.
- **Empty / null input:** the cited sources do not define statistics over zero observations; the
  repository contract returns a **zeroed `QualityStatistics` (TotalBases = 0)** / 0.0 percentage
  rather than throwing. This is a documented API decision, not a sourced value.

## Implementation surface

`QualityScoreAnalyzer` (in `Seqeron.Genomics.IO`) exposes three entry points; the algorithm spec
(QUALITY-STATS-001) pins these to the theory above:

- `CalculateStatistics(string, QualityEncoding)` — full summary for one quality string, returning the
  `QualityStatistics` record: `MeanQuality`/`MedianQuality`/`StandardDeviation` (double),
  `MinQuality`/`MaxQuality`/`TotalBases`/`BasesAboveQ20`/`BasesAboveQ30` (int),
  `PercentAboveQ20`/`PercentAboveQ30` (double 0–100), and `PerPositionMeanQuality`.
- `CalculateStatistics(IEnumerable<string>, QualityEncoding)` — **delegate variant** that aggregates
  scores across a *set* of reads and reports **`PerPositionMeanQuality`**: the mean decoded score at
  each base position across reads (the FastQC "per-base sequence quality" surface). For a single
  string this per-position vector is just that string's per-base scores.
- `CalculateQ30Percentage(string, QualityEncoding)` — the canonical `% ≥ Q30` shortcut, guaranteed
  equal to `CalculateStatistics(...).PercentAboveQ30` (INV-04).

**Complexity:** `CalculateStatistics` is **O(n log n)** — the sort for the median dominates; mean, σ,
and the Q20/Q30 counts are each O(n). `CalculateQ30Percentage` is a single O(n) pass. This unit does
no substring/pattern search, so the repository suffix tree is **not applicable**.

**Probability-averaged quality lives elsewhere.** Because `MeanQuality` is the arithmetic mean of
log-scaled Q (the samtools/FastQC "average quality"), it *understates* a read's true error rate — one
very-low-Q base dominates the error probability far more than the Q-mean suggests. The
error-probability summaries (`P̄ = (1/N) Σ 10^(−Qᵢ/10)`, expected-error sum `Σ 10^(−Qᵢ/10)`) are
provided by the separate `CalculateExpectedErrors` / `PhredToErrorProbability` methods, **not** by
`MeanQuality`.
