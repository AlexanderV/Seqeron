---
type: source
title: "Evidence: QUALITY-STATS-001 (FASTQ quality statistics — Q20/Q30, mean, median, min/max, std dev)"
tags: [validation, file-io]
doc_path: docs/Evidence/QUALITY-STATS-001-Evidence.md
sources:
  - docs/Evidence/QUALITY-STATS-001-Evidence.md
source_commit: 8f4f606b371e3b0654cce99998263706c7d8c185
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: QUALITY-STATS-001

The validation-evidence artifact for test unit **QUALITY-STATS-001** — **Quality Statistics**:
summary statistics over the decoded per-base Phred scores of a FASTQ quality string
(mean, median, min/max, population variance/std dev, and the **% ≥ Q20 / ≥ Q30** fractions),
plus the `CalculateQ30Percentage` convenience shortcut. This is a **QUALITY**-family Evidence file
and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern.

Its subject is a **distinct surface** from Phred encoding: the artifact reuses the Phred+33/+64
decode (contract of the sibling unit [[quality-phred-001-evidence]], factored into
[[phred-quality-encoding]]) only to obtain scores, then computes descriptive statistics over them.
That distinctness is called out in the encoding concept itself, so this ingest **creates a
dedicated concept** — [[fastq-quality-statistics]] — which **depends on** (consumes decoded scores
from) [[phred-quality-encoding]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Authoritative sources:**
  - **Illumina — Sequencing Quality Scores** (rank 2, platform spec): `Q = −10·log₁₀(e)`; Q-score
    table (Q20 → 99%, Q30 → 99.9%); **`% ≥ Q30` is the standard NGS run-quality benchmark**
    (inclusive `≥`).
  - **Newcastle University ASK** (rank 1, teaching material): population variance
    `σ² = (1/N) Σ (xᵢ − μ)²` and population std dev `σ = √σ²` — divide by **N** (not N−1) because
    the quality string is the complete population, not a sample.
  - **Math is Fun** (rank 4): median = middle of the sorted list; **even count** averages the two
    central order statistics.
  - **Wikipedia / Ewing & Green 1998** (rank 4/1): the Phred `Q = −10 log₁₀ P` formula and Q→accuracy
    ladder (provenance for the score definition).
  - **Cock et al. 2010** (rank 1, cited only): Phred+33 decode `ord(c) − 33` (Phred+64 `ord(c) − 64`)
    — used to convert the quality string to scores; the decode itself is QUALITY-PHRED-001's contract
    and is **not re-tested here**.
- **Statistics contract:** mean = arithmetic mean of decoded scores; median (odd = middle, even =
  mean of the two central values); min/max; **population** variance/σ (÷N); `% ≥ Q20` and `% ≥ Q30`
  with **inclusive** thresholds; `CalculateQ30Percentage(...)` must equal
  `CalculateStatistics(...).PercentAboveQ30`.
- **Worked oracles (test fixtures):**
  - Phred+33 `5?I` → scores 20/30/40 → mean 30.0, median 30, min/max 20/40, variance 200/3 ≈ 66.6667,
    σ = √(200/3) ≈ 8.16496580927726, `% ≥ Q20` = 100.0, `% ≥ Q30` ≈ 66.6667 (the Q20 base excluded).
  - Even-length `5II?` → sorted 20/30/40/40 → **median (30+40)/2 = 35.0**, mean 32.5.
  - Single base `I` (Q40) → mean = median = min = max = 40, **σ = 0.0**, `% ≥ Q20 = % ≥ Q30 = 100.0`.
- **Encoding-independence:** the statistics operate on decoded scores, so a Phred+64 input decodes to
  the same scores and yields identical statistics (a `COULD`-test).
- **Must-test surface:** mean/median(odd)/min/max/population σ on a derived score set; `% ≥ Q30` and
  `% ≥ Q20` with inclusive thresholds on a below/at/above-Q30 mix; even-count median averaging;
  `CalculateQ30Percentage` == `CalculateStatistics(...).PercentAboveQ30`; single-base σ=0 boundary;
  all-≥Q30 → 100% and none-≥Q30 → 0% boundaries.

## Deviations and assumptions

The artifact records **no source contradictions**. One API-shape assumption:

- **Empty/null quality string** → the canonical methods return an all-zero `QualityStatistics`
  (TotalBases = 0) / 0.0 percentage rather than throwing. The cited sources do not define summary
  statistics over zero observations; the zeroed-result contract is a repository decision
  (non-correctness-affecting for any non-empty input) — no numeric value invented.
