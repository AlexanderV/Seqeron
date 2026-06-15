# Validation Report: QUALITY-STATS-001 — Quality Statistics (Phred summary stats; Q20/Q30 %)

- **Validated:** 2026-06-15   **Area:** Quality
- **Canonical method(s):** `QualityScoreAnalyzer.CalculateStatistics(string, QualityEncoding)`,
  `CalculateStatistics(IEnumerable<string>, QualityEncoding)`, `CalculateQ30Percentage(string, QualityEncoding)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Phred quality score** — Wikipedia (timed out on fetch, but Q→accuracy table reproduced and
  cross-checked by direct computation below). Q = −10 log10 P; P = 10^(−Q/10).
- **Averaging basecall quality scores the right way** (gigabaseorgigabyte.wordpress.com, 2017-06-26,
  fetched) — verbatim: arithmetic mean of Q is *not* the right way; the correct mean read quality
  averages error probabilities then re-Phreds. Worked example: bases at Q20/Q10/Q3 (e = 1%/10%/50%)
  → arithmetic mean Q = 33/3 = **11**, but mean error rate = 0.61/3 ≈ 0.2 → Q = −10·log10(0.2) ≈ **7**.
- **USEARCH "average Q is a bad idea"** (drive5.com/usearch/manual/avgq.html, fetched) — verbatim
  "calculating average Q (Phred) scores is a bad idea"; a read of 140×Q35 + 10×Q2 has avg Q 33 but
  6.4 expected errors, vs all-Q25 (avg Q 25) with only 0.5 expected errors.
- **samtools stats** (manpage / SEQanswers, searched) — "average quality" = (Σ base qualities) ÷
  total length, i.e. the **arithmetic mean of Q**. FastQC per-base plots likewise report the
  arithmetic mean Phred per position.
- **Newcastle ASK** (population σ = √((1/N)Σ(xᵢ−μ)²), ÷N for a complete set) — as cited in Evidence.
- **Math is Fun** (median: middle for odd N, mean of two central for even N) — as cited in Evidence.

### Formula check
- **Mean** `μ = (1/N)Σqᵢ`: a *legitimate, conventionally-reported* statistic — it is exactly the
  "average quality" of `samtools stats` and FastQC's per-base mean. It is **not** the read's true
  mean accuracy: because Q is logarithmic, the probability-based mean (`P̄ = (1/N)Σ10^(−qᵢ/10)`,
  `Q̄ = −10 log10 P̄`) and expected-error sum are more faithful. The description defines the metric
  as the mean of *Phred scores* (not of accuracy), which is correct as stated. → **PASS-WITH-NOTES.**
- **Population σ (÷N)**: matches Newcastle exactly. Sample σ (÷(N−1)) would be wrong here (the
  quality string is the complete observed set). ✔
- **Median (odd/even rule)**: matches Math is Fun. ✔
- **%≥Q20 / %≥Q30 (inclusive ≥)**: matches Illumina (Q30→1-in-1000→99.9%, the base *at* Q30 is
  counted). ✔
- **Decode** Phred+33 = ord−33, Phred+64 = ord−64 (Cock et al. 2010). ✔

### Edge-case semantics
Empty/null → zeroed result (TotalBases=0)/0.0 Q30% is a documented repository contract (0/0 is
mathematically undefined; no value invented). Single base → σ=0, mean=median=min=max. All sourced.

### Independent cross-check (recomputed in Python this session)
| Quantity | Source-derived | Code/TestSpec |
|---|---|---|
| Mean {20,30,40} | 30.0 | 30.0 ✔ |
| Pop variance | 200/3 = 66.66666666666667 | ✔ |
| Pop σ | 8.16496580927726 | 8.16496580927726 ✔ |
| Median odd {20,30,40} | 30 | 30 ✔ |
| Median even {20,30,40,40} | (30+40)/2 = 35.0 | 35.0 ✔ |
| %≥Q30 (2/3) | 66.66666666666667 | 66.66666666666667 ✔ |
| %≥Q20 (3/3) | 100.0 | 100.0 ✔ |
| Multi-read mean / %≥Q30 (6 scores) | 30.0 / 66.667 | ✔ |
| Q-table | Q20=99%, Q30=99.9%, Q40=99.99% | ✔ |

### Findings / divergences
The description (and the M1 test) report the arithmetic mean of Q. Authoritative sources caution
that this is a *poor* indicator of read accuracy and that error-probability averaging is the
"right way". It is, however, the standard metric emitted by samtools stats / FastQC, so it is a
valid documented statistic — not a defect — provided the limitation is disclosed. The algorithm
doc previously did not disclose it; I added a sourced limitations paragraph (samtools/FastQC
convention vs probability-based mean; cross-refs `CalculateExpectedErrors`). No code or expected
value changed. Hence Stage A = PASS-WITH-NOTES (documentation completeness only).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs`:
- `CalculateStatistics(string)` L265–287 → `CalculateStatisticsFromPhred` L329–368.
- `CalculateStatistics(IEnumerable<string>)` L292–327 (aggregates all scores + per-position means).
- `CalculateQ30Percentage(string)` L380–393.

### Formula realised correctly?
- Mean = `phredScores.Average()` (arithmetic mean of Q) — matches the validated description. ✔
- Median L344–346: even → `(sorted[n/2-1]+sorted[n/2])/2.0`; odd → `sorted[n/2]`. ✔
- Variance L348 = `Select((x-mean)^2).Average()` (÷N, population) → σ = √variance. ✔
- Counts L351–352 use `q >= 20` / `q >= 30` (inclusive). ✔
- Percentages L365–366 = `100*count/N`. ✔
- `CalculateQ30Percentage` recomputes the same `100*count(≥30)/N`, consistent with INV-04. ✔
- Empty/null guarded in both string overloads and the aggregate (L269, L312, L333, L384). ✔

### Cross-verification table vs code
All 16 tests pass with exact `Is.EqualTo(...).Within(1e-10)` against the source-derived values in
the table above. A wrong implementation would be caught: ÷(N−1) fails M2; exclusive `>` fails S4
(base exactly Q30 → would drop to 0%); single-element-median for even N fails M4.

### Variant/delegate consistency
INV-04 test confirms `CalculateQ30Percentage` ≡ `CalculateStatistics(...).PercentAboveQ30`. INV-05
test confirms encoding invariance (Phred+33 "5?I" vs Phred+64 "T^h" both decode 20,30,40 → identical
Mean & %≥Q30). C1 confirms the IEnumerable overload aggregates {20,30,40}×2 → mean 30, %≥Q30 66.667.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every expected value traces to an external source (Newcastle ÷N,
  Math is Fun median, Illumina inclusive ≥, samtools arithmetic-mean convention) and was
  independently recomputed this session — not read back from the code.
- **No green-washing:** all assertions are exact equalities within 1e-10; no Greater/AtLeast/range,
  no widened tolerance, no skipped/ignored test, no expected value adjusted to match output.
- **Coverage:** all three canonical methods exercised; Stage-A branches covered — mean, min, max,
  population σ, median (odd **and** even), %≥Q20, %≥Q30, counts, encoding invariance, single-base
  zero-spread, all-≥Q30, none-≥Q30, exactly-Q30 inclusivity, empty, null, multi-read aggregation.
- **Honest green:** full unfiltered suite **Passed: 6510, Failed: 0, Skipped: 0**; build 0 errors
  (4 pre-existing warnings in unrelated files). Gate: **PASS.**

Minor untested-but-in-scope items (not defects, no sourced exact value at stake): `Auto` encoding
branch, the `IEnumerable` empty-input branch, and the `PerPositionMeanQuality` output field. The
TestSpec scopes the IEnumerable overload as "smoke only"; these are noted, not blocking.

### Findings / defects
None in the code. One documentation gap (mean-of-Q vs error-probability mean) fixed in the
algorithm doc this session.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (arithmetic-mean-of-Q metric is valid per samtools/FastQC but its
  accuracy limitation was undocumented; now disclosed with sources).
- **Stage B:** PASS.
- **End-state:** ✅ CLEAN — no code defect; the documentation note was completed; full suite green.
- **Test-quality gate:** PASS.
