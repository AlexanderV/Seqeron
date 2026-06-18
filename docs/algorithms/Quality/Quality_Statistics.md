# Quality Statistics (Phred Summary Statistics, Q20/Q30 Percentages)

| Field | Value |
|-------|-------|
| Algorithm Group | Quality |
| Test Unit ID | QUALITY-STATS-001 |
| Related Projects | Seqeron.Genomics.IO |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Computes descriptive summary statistics over the per-base Phred quality scores of a FASTQ
read (or a set of reads): mean, median, minimum, maximum, population standard deviation, and
the percentage of bases at or above the Q20 (≥20) and Q30 (≥30) quality thresholds. These are
the standard run/read quality-control metrics in next-generation sequencing; %≥Q30 in
particular is the de-facto benchmark for run quality [3]. The computation is exact (closed-form
descriptive statistics over decoded integer scores), not heuristic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A FASTQ quality character encodes a Phred quality score `Q`, which is tied to the estimated
base-calling error probability `e` by `Q = -10 log10 e` [1][3]. Thus Q20 ⇒ e = 1-in-100
(99% accuracy) and Q30 ⇒ e = 1-in-1000 (99.9% accuracy) [2][3]. Summary statistics over the
decoded scores describe the overall confidence of a read or dataset; the fraction of bases
≥ Q30 is the headline quality metric reported for sequencing runs [3].

### 2.2 Core Model

Let `q₁..q_N` be the decoded Phred scores (N = number of bases). Decoding uses the FASTQ
encoding: Phred+33 score = ord(c) − 33, Phred+64 score = ord(c) − 64 [4].

- Mean: `μ = (1/N) Σ qᵢ`.
- Median: middle value of the sorted scores for odd N; the arithmetic mean of the two central
  order statistics for even N [5].
- Population variance / standard deviation: `σ² = (1/N) Σ (qᵢ − μ)²`, `σ = √σ²` [6]. The
  population divisor (÷N) is used because the quality string is the complete observed set of
  scores, not a sample [6].
- %≥Q20: `100 · |{i : qᵢ ≥ 20}| / N` (inclusive) [2][3].
- %≥Q30: `100 · |{i : qᵢ ≥ 30}| / N` (inclusive) [2][3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | min ≤ mean ≤ max for non-empty input | mean is a convex combination of the scores |
| INV-02 | σ ≥ 0, and σ = 0 iff all scores are equal | σ = √(mean of squared deviations) [6] |
| INV-03 | 0 ≤ %≥Q30 ≤ %≥Q20 ≤ 100 | {q ≥ 30} ⊆ {q ≥ 20} [2] |
| INV-04 | `CalculateQ30Percentage(s)` = `CalculateStatistics(s).PercentAboveQ30` | both = 100·|{q ≥ 30}|/N [3] |
| INV-05 | statistics depend only on decoded scores, not the encoding used | score is encoding-invariant [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| qualityString | string | required | FASTQ quality line | null/empty → zeroed result |
| encoding | QualityEncoding | Phred33 | Phred+33 / Phred+64 / Auto | decode offset 33 or 64 [4] |

### 3.2 Output / Return Value

`QualityStatistics` record: `MeanQuality`, `MedianQuality` (double), `MinQuality`, `MaxQuality`
(int), `StandardDeviation` (double, population), `TotalBases` (int), `BasesAboveQ20`,
`BasesAboveQ30` (int counts), `PercentAboveQ20`, `PercentAboveQ30` (double, 0–100),
`PerPositionMeanQuality` (per-base means; for a single string this is the per-base score).
`CalculateQ30Percentage` returns the double `PercentAboveQ30`.

### 3.3 Preconditions and Validation

0-based scoring. Null or empty `qualityString` returns an all-zero `QualityStatistics`
(TotalBases = 0) / `0.0` for `CalculateQ30Percentage` rather than throwing. Q20/Q30 thresholds
are inclusive (≥). Encoding decoding is the contract of QUALITY-PHRED-001 [4].

## 4. Algorithm

### 4.1 High-Level Steps

1. Decode the quality string to integer Phred scores using the chosen encoding [4].
2. If no scores, return the zeroed result.
3. Compute mean, sorted-order median (odd/even rule), min, max.
4. Compute population variance (÷N) and its square root [6].
5. Count bases ≥ Q20 and ≥ Q30; convert each to a percentage of N [2][3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateStatistics | O(n log n) | O(n) | sort dominates the median; counts/mean/σ are O(n) |
| CalculateQ30Percentage | O(n) | O(n) | single pass over decoded scores |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [QualityScoreAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs)

- `QualityScoreAnalyzer.CalculateStatistics(string, QualityEncoding)`: full summary statistics for one quality string.
- `QualityScoreAnalyzer.CalculateStatistics(IEnumerable<string>, QualityEncoding)`: aggregates scores across reads and reports per-position means (delegate variant).
- `QualityScoreAnalyzer.CalculateQ30Percentage(string, QualityEncoding)`: %bases ≥ Q30 (canonical Q30 metric).

### 5.2 Current Behavior

Median for even-length input averages the two central order statistics. Standard deviation
uses the population divisor (÷N). Q20/Q30 use inclusive `>=`. This unit performs no
substring/pattern search, so the repository suffix tree is **not applicable**.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Population standard deviation `σ = √((1/N) Σ(qᵢ−μ)²)` [6].
- Median: middle for odd N, average of the two central values for even N [5].
- %≥Q30 and %≥Q20 with inclusive thresholds; Q30 = 99.9% accuracy benchmark [2][3].
- Phred+33 / Phred+64 decode prior to statistics [4].

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none — the canonical and delegate methods cover the QUALITY-STATS-001 scope)

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty string | zeroed `QualityStatistics` / `0.0` Q30% | no observations; documented contract |
| single base | mean=median=min=max=score, σ=0 | zero spread [6] |
| base exactly Q30 | counted in %≥Q30 | inclusive ≥ [2][3] |
| all bases ≥ Q30 | %≥Q30 = 100 | INV-03 |

### 6.2 Limitations

`MeanQuality` is the **arithmetic mean of the decoded Phred scores** (`(1/N) Σ qᵢ`), matching
the "average quality" reported by `samtools stats` (Σ base qualities ÷ length) and FastQC's
per-base mean. Because Phred Q is logarithmic (`Q = -10 log10 e`), this arithmetic mean is a
crude indicator of a read's *overall accuracy*: a single very-low-Q base contributes far more to
the true error rate than the mean of Q suggests. A read's expected error rate is more faithfully
summarised by averaging error probabilities (`P̄ = (1/N) Σ 10^(-qᵢ/10)`, then `Q̄ = -10 log10 P̄`)
or by the expected-error sum `Σ 10^(-qᵢ/10)`. This unit reports the conventional arithmetic mean
of Q (the samtools/FastQC metric); the probability-based summaries are provided separately by
`CalculateExpectedErrors` and `PhredToErrorProbability`, not by `MeanQuality`. See *Averaging
basecall quality scores the right way* (gigabaseorgigabyte.wordpress.com, 2017) and the USEARCH
"average Q is a bad idea" note (drive5.com/usearch/manual/avgq.html).

Statistics summarize quality only; they do not assess sequence content, adapter contamination,
or per-tile artefacts. For an empty input the percentages are defined as 0 by contract rather
than mathematically (0/0 is undefined).

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var stats = QualityScoreAnalyzer.CalculateStatistics("5?I"); // Phred+33 → 20,30,40
// stats.MeanQuality == 30.0; stats.StandardDeviation == 8.16496580927726
// stats.PercentAboveQ30 == 66.66666666666667 (2 of 3 ≥ Q30)
double q30 = QualityScoreAnalyzer.CalculateQ30Percentage("5?I"); // 66.66666666666667
```

Scores 20,30,40: μ = 30; deviations −10,0,10; σ² = 200/3 ≈ 66.667; σ ≈ 8.16497. Two of three
bases (30,40) are ≥ Q30 → 66.667%.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [QualityScoreAnalyzer_CalculateStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_CalculateStatistics_Tests.cs) — covers INV-01..INV-05
- Evidence: [QUALITY-STATS-001-Evidence.md](../../../docs/Evidence/QUALITY-STATS-001-Evidence.md)
- Related algorithms: [Phred Score Handling](../Quality/Phred_Score_Handling.md)

## 8. References

1. Ewing B, Hillier L, Wendl MC, Green P. 1998. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research 8(3):175–185. https://doi.org/10.1101/gr.8.3.175
2. Ewing B, Green P. 1998. Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Research 8(3):186–194. https://doi.org/10.1101/gr.8.3.186
3. Illumina, Inc. Sequencing Quality Scores. https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html
4. Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM. 2010. The Sanger FASTQ file format for sequences with quality scores. Nucleic Acids Research 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137
5. Math is Fun. How to Find the Median Value. https://www.mathsisfun.com/median.html
6. Newcastle University, ASK Academic Skills Kit. Variance and Standard Deviation. https://www.mas.ncl.ac.uk/ask/numeracy-maths-statistics/statistics/descriptive-statistics/variance-and-standard-deviation.html
