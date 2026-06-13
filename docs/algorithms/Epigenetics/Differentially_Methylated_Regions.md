# Differentially Methylated Regions (DMR)

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-DMR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

A Differentially Methylated Region (DMR) is a genomic interval of adjacent CpG sites whose DNA-methylation level differs between two samples or conditions [1]. This algorithm compares two single-sample methylation profiles, tiles the genome into fixed-size windows, and for each window computes the mean per-site methylation difference and a two-sided Fisher's exact test p-value over the pooled methylated/unmethylated read counts, following the methylKit tiling-window model [1][2]. It is specification-driven (the methylKit defaults and the exact Fisher test are reproduced), not heuristic. Use it to locate hyper- and hypo-methylated regions between, for example, tumour and normal samples.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Bisulfite sequencing reports, at each cytosine, the number of methylated (C) and unmethylated (T) reads. The per-site methylation level is the ratio C/(C+T) [1]. Differential methylation analysis asks whether this proportion differs between two groups at a base (DMC) or over a region (DMR) [1]. methylKit groups bases into tiling windows of a fixed width and tests the summed counts per window [2].

### 2.2 Core Model

For a window containing covered cytosines, let the pooled counts be:

|          | methylated | unmethylated |
|----------|------------|--------------|
| sample1  | numC1 (a)  | numT1 (b)    |
| sample2  | numC2 (c)  | numT2 (d)    |

The probability of one such 2×2 table with fixed margins is the hypergeometric probability [3]:

`p = (a+b)! (c+d)! (a+c)! (b+d)! / (a! b! c! d! n!)`, with `n = a+b+c+d`.

The two-sided Fisher's exact p-value is the sum of `p` over every table with the same margins whose probability is ≤ that of the observed table [3].

The region methylation difference is the mean over the window of the per-site differences (treatment − control), in fraction units: `meanDiff = mean(level2 − level1)` [1]. A region is reported when `|meanDiff| > minDifference` (strict, default 0.25 = 25%) and is labelled "Hypermethylated" when `meanDiff > 0` (treatment higher than control) or "Hypomethylated" when `meanDiff < 0` [1][2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | One sample per group (no replicates), so Fisher's exact test is the appropriate test [4]. | With biological replicates, logistic regression is required instead; Fisher's exact (on pooled reads) ignores between-replicate variance and overstates significance. |
| ASM-02 | Per-site methylated/unmethylated counts can be reconstructed from `MethylationLevel × Coverage` (rounded). | If the stored level/coverage do not encode integer read counts, the reconstructed 2×2 table — and hence the p-value — is approximate. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported DMR has `|MeanDifference| > minDifference` (strict). | getMethylDiff cutoff is `meth.diff > difference` [3]. |
| INV-02 | Annotation is "Hypermethylated" iff MeanDifference > 0, "Hypomethylated" iff < 0. | hyper = higher than control; hypo = lower [1][2]. |
| INV-03 | `0 ≤ PValue ≤ 1`. | It is a sum of hypergeometric probabilities clamped to 1 [3]. |
| INV-04 | Positions ≥ windowSize apart fall in different windows. | Fixed-width tiling (win.size) [2]. |
| INV-05 | A reported DMR contains ≥ minCpGCount covered sites. | A DMR is a region of adjacent CpGs [1]. |
| INV-06 | Output is ordered by Start ascending and deterministic. | Positions are sorted before tiling; no randomness [2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sample1 | `IEnumerable<MethylationSite>` | required | Control profile | non-null |
| sample2 | `IEnumerable<MethylationSite>` | required | Treatment profile | non-null |
| windowSize | `int` | 1000 | Tiling window width in bp [2] | > 0 |
| minDifference | `double` | 0.25 | Min absolute mean methylation difference, fraction [0,1]; strict `>` [3] | [0,1] |
| minCpGCount | `int` | 3 | Min covered sites per region [1] | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Start | int | 0-based start position of the window (first position of the window). |
| End | int | Position of the last covered site in the window. |
| MeanDifference | double | Mean(level2 − level1) over the window, fraction in [−1,1]. |
| PValue | double | Two-sided Fisher's exact p-value of the pooled 2×2 table [3]. |
| CpGCount | int | Number of covered sites in the region. |
| Annotation | string | "Hypermethylated" / "Hypomethylated" (or a feature label after `AnnotateDMRs`). |

### 3.3 Preconditions and Validation

`sample1`/`sample2` must be non-null (`ArgumentNullException` otherwise). Empty input yields no DMRs. Positions are 0-based and in the same coordinate space across both samples and any annotations. A position present in only one sample is compared against an implicit methylation level of 0 in the other; a position with `Coverage == 0` contributes no reads to the 2×2 table.

## 4. Algorithm

### 4.1 High-Level Steps

1. Index both samples by position; take the sorted union of positions. Empty union → no DMRs.
2. Tile positions into windows: open a new window when a position is ≥ `windowSize` from the current window start.
3. Per window: require ≥ `minCpGCount` covered sites; compute `meanDiff = mean(level2 − level1)`; pool methylated/unmethylated read counts into a 2×2 table.
4. If `|meanDiff| > minDifference` (strict), emit a DMR with a two-sided Fisher's exact p-value and a hyper/hypo annotation by the sign of `meanDiff`.
5. `AnnotateDMRs` relabels each DMR with the first overlapping genomic feature, if any.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Tiling defaults: `win.size = step.size = 1000` bp, `cov.bases = 0` [2].
- Reporting cutoff: `|meth.diff| > 25%` (0.25 fraction) and (in methylKit) `q-value < 0.01`; this unit reports the raw Fisher p-value and applies the strict difference cutoff [1][3].
- Hyper/hypo by the sign of the difference [1][3].
- Fisher's exact p-value computed in log-factorial space to avoid overflow.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindDMRs | O(n log n) | O(n) | n = total positions; dominated by the sort of the position union. Per-window Fisher enumeration is O(m) tables in the smaller margin m (read counts). |
| FisherExactProbability | O(k) | O(1) | k = total reads (log-factorial sum). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.FindDMRs(...)`: tiling-window DMR detection with two-sided Fisher's exact significance.
- `EpigeneticsAnalyzer.AnnotateDMRs(...)`: relabels DMRs with overlapping genomic features.
- `EpigeneticsAnalyzer.FisherExactProbability(a,b,c,d)`: single-table hypergeometric probability (public, used by the two-sided test and directly testable against the published worked example).

### 5.2 Current Behavior

- The two-sided Fisher's exact p-value enumerates all tables sharing the observed margins and sums those with probability ≤ the observed (with a 1e-7 tie tolerance). Degenerate margins (a zero row/column total) return p = 1.0.
- Read counts are reconstructed from `MethylationLevel × Coverage` rounded to the nearest integer (ASM-02), because `MethylationSite` stores a fraction plus coverage rather than raw C/T counts.
- This is not a search/matching unit; the repository suffix tree is not applicable (no substring or occurrence search is performed).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-site methylation level C/(C+T) and the meth.diff = treatment − control definition [1].
- Strict reporting cutoff `|meth.diff| > difference` and sign-based hyper/hypo labelling [1][3].
- Fixed-width tiling (default 1000 bp) [2].
- Two-sided Fisher's exact test (hypergeometric, sum of tables with probability ≤ observed) [3][4].

**Intentionally simplified:**

- q-value / multiple-testing correction: the unit returns the raw Fisher p-value; **consequence:** users wanting methylKit's SLIM q-value < 0.01 gate must apply correction downstream.
- Single sample per group only; **consequence:** with replicates the statistically correct test is logistic regression, not implemented here (ASM-01).

**Not implemented:**

- Logistic-regression differential methylation for replicated designs; **users should rely on:** methylKit `calculateDiffMeth` (R) for replicated data — no in-repo alternative.
- CpG/CHG/CHH context stratification of DMRs; **users should rely on:** `GetMethylationContext` upstream to pre-filter sites by context.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Pooled-window Fisher vs per-base test | Assumption | p-value reflects the window's pooled counts, matching tileMethylCounts→Fisher | accepted | ASM-01; Evidence Assumption 1 |
| 2 | Count reconstruction from level×coverage | Assumption | rounding can shift small counts by 1 | accepted | ASM-02; Evidence Assumption 2 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | No DMRs | No tiles to test [2] |
| Window with < minCpGCount sites | Not reported | A DMR is a region of adjacent CpGs [1] |
| `|meanDiff|` == minDifference exactly | Not reported | Strict cutoff `>` [3] |
| Zero coverage in one group within a window | Fisher p = 1.0 | Degenerate fixed margin [3] |
| Null sample | `ArgumentNullException` | Input-validation contract |

### 6.2 Limitations

Single-sample-per-group only (no replicate modelling); no q-value correction; reconstructs read counts from fractional levels; does not stratify by sequence context; the simple tiling does not merge adjacent significant windows into larger DMRs (each window is reported independently).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var control = Enumerable.Range(0, 3)
    .Select(i => new EpigeneticsAnalyzer.MethylationSite(i, EpigeneticsAnalyzer.MethylationType.CpG, "CG", 0.0, 20));
var treatment = Enumerable.Range(0, 3)
    .Select(i => new EpigeneticsAnalyzer.MethylationSite(i, EpigeneticsAnalyzer.MethylationType.CpG, "CG", 1.0, 20));

var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();
// dmrs[0].MeanDifference == 1.0, dmrs[0].Annotation == "Hypermethylated"
```

**Numerical walk-through (Fisher single-table probability):** for the table a=1, b=9, c=11, d=3 (n=24), `p = 10! · 14! · 12! · 12! / (1! · 9! · 11! · 3! · 24!) ≈ 0.001346076` [3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_DMR_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_DMR_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [EPIGEN-DMR-001-Evidence.md](../../../docs/Evidence/EPIGEN-DMR-001-Evidence.md)
- Related algorithms: [Methylation_Analysis](../Epigenetics/Methylation_Analysis.md)

## 8. References

1. Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE. 2012. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biology 13:R87. https://doi.org/10.1186/gb-2012-13-10-r87
2. methylKit `tileMethylCounts` reference manual (win.size=1000, step.size=1000, cov.bases=0). al2na/methylKit, GitHub. https://github.com/al2na/methylKit/blob/master/man/tileMethylCounts-methods.Rd
3. Fisher's exact test (hypergeometric probability of a 2×2 table; two-sided p-value; worked example). Wikipedia, citing Fisher RA (1922, 1935). https://en.wikipedia.org/wiki/Fisher's_exact_test
4. methylKit `calculateDiffMeth` reference manual (Fisher's exact test applied when one sample per group). al2na/methylKit, GitHub. https://github.com/al2na/methylKit/blob/master/man/calculateDiffMeth-methods.Rd
