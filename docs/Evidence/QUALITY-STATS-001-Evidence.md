# Evidence Artifact: QUALITY-STATS-001

**Test Unit ID:** QUALITY-STATS-001
**Algorithm:** Quality Statistics (Phred score summary statistics; Q20/Q30 percentages)
**Date Collected:** 2026-06-13

---

## Online Sources

### Phred quality score (Wikipedia, citing Ewing & Green 1998 primaries)

**URL:** https://en.wikipedia.org/wiki/Phred_quality_score
**Accessed:** 2026-06-13 (fetched in session)
**Authority rank:** 4 (Wikipedia citing primary sources; the primaries are Ewing et al. 1998 ranked 1)

**Key Extracted Points:**

1. **Phred Q formula:** The Phred quality score Q is related to base-calling error probability P by `Q = -10 log10 P`, equivalently `P = 10^(-Q/10)`. (Retrieved verbatim.)
2. **Q→accuracy table (verbatim):** Q10 = 1 in 10 = 90%; Q20 = 1 in 100 = 99%; Q30 = 1 in 1,000 = 99.9%; Q40 = 1 in 10,000 = 99.99%; Q50 = 1 in 100,000 = 99.999%; Q60 = 1 in 1,000,000 = 99.9999%.
3. **Provenance:** phred was developed by Ewing B, Hillier L, Wendl MC, Green P; foundational papers Ewing et al. (1998) "Base-calling of automated sequencer traces using phred. I" *Genome Research* 8(3):175–185 and Ewing & Green (1998) "…II. Error probabilities" *Genome Research* 8(3):186–194.

### Illumina — Sequencing Quality Scores (official spec/vendor page)

**URL:** https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html
**Accessed:** 2026-06-13 (fetched in session)
**Authority rank:** 2 (platform specification governing FASTQ quality semantics)

**Key Extracted Points:**

1. **Q-score definition (verbatim):** `Q = -10 log10(e)` where "e is the estimated probability of the base call being wrong."
2. **Q-score table (verbatim):** Q10 → 1 in 10 → 90%; Q20 → 1 in 100 → 99%; Q30 → 1 in 1,000 → 99.9%.
3. **Q30 metric:** "Q30 is considered a benchmark for quality in next-generation sequencing (NGS)"; at Q30 "virtually all of the reads will be perfect, with no errors or ambiguities." The fraction of bases ≥ Q30 (%≥Q30) is the standard run-quality metric.

### Population variance / standard deviation (Newcastle University, ASK numeracy)

**URL:** https://www.mas.ncl.ac.uk/ask/numeracy-maths-statistics/statistics/descriptive-statistics/variance-and-standard-deviation.html
**Accessed:** 2026-06-13 (fetched in session)
**Authority rank:** 1 (university teaching material; standard textbook statistics)

**Key Extracted Points:**

1. **Population variance (verbatim):** `σ² = (1/N) Σᵢ₌₁ᴺ (xᵢ − μ)²`, N = population size, μ = population mean.
2. **Population standard deviation (verbatim):** `σ = √( (1/N) Σᵢ₌₁ᴺ (xᵢ − μ)² )` — the positive square root of the variance.
3. **Population vs sample:** divide by N for a complete population; divide by N−1 (Bessel's correction) for a sample. The quality string is the complete set of observed scores, so population (÷N) applies.

### Median definition (Math is Fun)

**URL:** https://www.mathsisfun.com/median.html
**Accessed:** 2026-06-13 (fetched in session)
**Authority rank:** 4 (standard descriptive-statistics reference)

**Key Extracted Points:**

1. **Odd count (verbatim):** the median is the "middle" of a sorted list of numbers.
2. **Even count (verbatim):** find the middle pair, then the value half way between them — "adding them together and dividing by two" (arithmetic mean of the two central values).

### FASTQ encoding (already in repo Evidence QUALITY-PHRED-001; Cock et al. 2010)

**URL:** https://doi.org/10.1093/nar/gkp1137 (cited in repo for QUALITY-PHRED-001; not re-fetched here)
**Authority rank:** 1
**Key Extracted Points:** Phred+33 decodes char→score as `ord(c) − 33` (scores 0–93); Phred+64 as `ord(c) − 64` (scores 0–62). Used here only to convert the quality string to scores prior to statistics; the decode itself is the contract of QUALITY-PHRED-001 and is not re-tested.

---

## Documented Corner Cases and Failure Modes

### From Newcastle University / Math is Fun

1. **Even vs odd count:** median branch differs; an even-length input must average the two central order statistics, an odd-length input takes the single middle element.
2. **Single element:** mean = median = min = max = that value; population standard deviation = 0 (zero spread).

### From Illumina / Phred

1. **Q20 / Q30 thresholds are inclusive (≥):** %≥Q30 counts bases whose score is 30 or higher; a base exactly at Q30 is counted (1-in-1000 error). Same for Q20 (≥20).
2. **Empty input:** no bases ⇒ statistics are undefined as a ratio; the contract returns a zeroed result with TotalBases = 0 (documented behavior, not a source value).

---

## Test Datasets

### Dataset: Hand-derived Phred+33 string "5?I" (interior scores 20/30/40)

**Source:** Phred formula + Cock et al. (2010) Phred+33 decode (ord−33).

| Parameter | Value |
|-----------|-------|
| Quality string | `5?I` |
| Decoded scores | 20, 30, 40 |
| Mean | (20+30+40)/3 = 30.0 |
| Median (odd, n=3) | 30 |
| Min / Max | 20 / 40 |
| Population variance | ((−10)²+0²+10²)/3 = 200/3 ≈ 66.6667 |
| Population std dev | √(200/3) ≈ 8.16496580927726 |
| Bases ≥ Q20 | 3 → 100.0% |
| Bases ≥ Q30 | 2 → 66.66666666666667% |

### Dataset: Even-length Phred+33 string "5II?" (median averages two central order statistics)

**Source:** median rule (even count) + Phred+33 decode.

| Parameter | Value |
|-----------|-------|
| Quality string | `5II?` |
| Decoded scores | 20, 40, 40, 30 |
| Sorted | 20, 30, 40, 40 |
| Median (even, n=4) | (30+40)/2 = 35.0 |
| Mean | (20+40+40+30)/4 = 32.5 |

### Dataset: Single-base string "I" (Q40)

**Source:** descriptive-statistics single-value identities.

| Parameter | Value |
|-----------|-------|
| Quality string | `I` |
| Decoded score | 40 |
| Mean / Median / Min / Max | 40 |
| Std dev | 0.0 |
| %≥Q20 / %≥Q30 | 100.0 / 100.0 |

---

## Assumptions

1. **ASSUMPTION: Empty-input return shape** — For an empty/null quality string the canonical methods return an all-zero `QualityStatistics` (TotalBases = 0) / 0.0 percentage rather than throwing. The cited sources do not define summary statistics over zero observations; the zeroed-result contract is a repository decision (non-correctness-affecting for any non-empty input). Documented as the contract; no numeric value invented.

---

## Recommendations for Test Coverage

1. **MUST Test:** mean, median (odd), min, max, population std dev on a derived score set — Evidence: Newcastle (σ ÷N), Math is Fun (median), Phred formula.
2. **MUST Test:** %≥Q30 and %≥Q20 with inclusive thresholds on a string mixing below/at/above Q30 — Evidence: Illumina (Q30 ≥30; Q20 ≥20).
3. **MUST Test:** median even-count averages the two central order statistics — Evidence: Math is Fun (even rule).
4. **MUST Test:** `CalculateQ30Percentage` equals `CalculateStatistics(...).PercentAboveQ30` for the same input — Evidence: Q30 = %bases ≥ Q30 (Illumina).
5. **SHOULD Test:** single base → std dev 0, mean=median=min=max — Rationale: zero-spread boundary.
6. **SHOULD Test:** all bases ≥ Q30 → 100%; no bases ≥ Q30 → 0% — Rationale: percentage boundaries.
7. **COULD Test:** Phred+64 input decodes to the same scores and yields identical statistics — Rationale: encoding-independence of the score-based statistics.

---

## References

1. Ewing B, Hillier L, Wendl MC, Green P (1998). Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research 8(3):175–185. https://doi.org/10.1101/gr.8.3.175 (provenance via https://en.wikipedia.org/wiki/Phred_quality_score)
2. Ewing B, Green P (1998). Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Research 8(3):186–194. https://doi.org/10.1101/gr.8.3.186 (provenance via https://en.wikipedia.org/wiki/Phred_quality_score)
3. Illumina, Inc. Sequencing Quality Scores. https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/quality-scores.html (accessed 2026-06-13)
4. Newcastle University, ASK Academic Skills Kit. Variance and Standard Deviation. https://www.mas.ncl.ac.uk/ask/numeracy-maths-statistics/statistics/descriptive-statistics/variance-and-standard-deviation.html (accessed 2026-06-13)
5. Math is Fun. How to Find the Median Value. https://www.mathsisfun.com/median.html (accessed 2026-06-13)
6. Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM (2010). The Sanger FASTQ file format... Nucleic Acids Research 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137 (used for the score decode underpinning QUALITY-PHRED-001)

---

## Change History

- **2026-06-13**: Initial documentation.
