# Expression Quantification (TPM, FPKM, Quantile Normalization)

| Field | Value |
|-------|-------|
| Algorithm Group | Transcriptome |
| Test Unit ID | TRANS-EXPR-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Expression quantification converts raw RNA-seq read counts into normalized abundance measures that are comparable within and (with care) across samples. This unit implements three established transforms: **TPM** (transcripts per million), **FPKM/RPKM** (fragments/reads per kilobase per million mapped reads), and **quantile normalization**. TPM and FPKM are deterministic per-gene formulas driven by read count, transcript length, and library depth; quantile normalization is a deterministic rank-based transform that forces every sample to share the same value distribution [1][2][4]. All three are exact (closed-form) rather than heuristic or probabilistic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A read count `X_i` mapped to transcript `i` confounds true abundance with transcript length and sequencing depth: longer transcripts and deeper libraries accrue more reads independent of molar concentration. FPKM corrects for both length and depth; TPM corrects the same factors but, unlike FPKM, normalizes so that the per-sample total is constant, restoring the relative-molar-concentration interpretation that FPKM violates [1][2]. Quantile normalization addresses a different problem — removing technical distribution differences between samples by making their value distributions identical [4].

### 2.2 Core Model

**TPM** [2 (verbatim), 3]:

```
TPM_i = (X_i / l_i) / Σ_j (X_j / l_j) · 10^6
```

where `X_i` = reads mapped to transcript i, `l_i` = transcript length. Equivalently `TPM = 10^6 · RPKM / Σ RPKM` [2].

**FPKM / RPKM** [2 (verbatim `RPKM = 10^9 · reads / (total · length)`), 3]:

```
FPKM_i = X_i · 10^9 / (l_i · N)
```

where `N` = total mapped reads/fragments in the sample. The `10^9` factor combines per-kilobase (`10^3`) and per-million-reads (`10^6`) scaling [3].

**Quantile normalization** [4, citing Bolstad et al. 2003]: sort each sample (column); the value at rank `r` across all samples is replaced by the arithmetic mean of the rank-`r` values; the rank means are placed back at each value's original position. Values tied within a sample receive the mean of the rank means they would otherwise span (tie-average rule) [4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Effective length ≈ annotated length (no fragment-length correction) | Absolute FPKM/TPM shift slightly for very short transcripts; within-sample ratios unaffected [3] |
| ASM-02 | Quantile normalization presumes samples share a common underlying distribution | Forcing identical distributions can mask genuine global expression differences [4] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Σ_i TPM_i = 10^6 for any sample with Σ(X/l) > 0 | TPM normalizes by Σ(X/l) then scales by 10^6 [1][2] |
| INV-02 | TPM_i ≥ 0; equal X/l ⇒ equal TPM | Non-negative counts/lengths; TPM is monotone in X/l [2] |
| INV-03 | FPKM_i ≥ 0; FPKM = 0 when l ≤ 0 or N ≤ 0 | Formula undefined for non-positive length/depth; convention returns 0 |
| INV-04 | Quantile normalization preserves within-column rank order | Rank → rank-mean is non-decreasing [4] |
| INV-05 | Untied columns map to a permutation of the same rank-mean multiset | Each rank is assigned exactly one rank mean [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| geneCounts | `IEnumerable<(string GeneId, double RawCount, int Length)>` | required | per-gene count and length | Length > 0 for a meaningful rate; RawCount ≥ 0 |
| rawCount (FPKM) | `double` | required | reads mapped to the transcript | ≥ 0 |
| length (FPKM) | `int` | required | transcript length | > 0 (else returns 0) |
| totalReads (FPKM) | `double` | required | total mapped reads N | > 0 (else returns 0) |
| samples (QN) | `IEnumerable<IEnumerable<double>>` | required | one inner sequence per sample | all samples same length |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CalculateTPM | `IEnumerable<GeneExpression>` | per-gene record with TPM and FPKM populated; TPM sums to 10^6 |
| CalculateFPKM | `double` | FPKM for the single gene |
| QuantileNormalize | `IEnumerable<IReadOnlyList<double>>` | one normalized vector per input sample, original positions preserved |

### 3.3 Preconditions and Validation

Empty `geneCounts` / empty `samples` / zero genes yield an empty sequence. TPM with an all-zero count denominator (Σ(X/l) = 0) returns TPM = 0 and FPKM = 0 for every gene (0/0 undefined). FPKM returns 0 when length ≤ 0 or totalReads ≤ 0. Quantile normalization uses the first sample's length as the gene count. No exceptions are raised for these degenerate inputs.

## 4. Algorithm

### 4.1 High-Level Steps

1. **TPM:** compute rate `X_i/l_i` per gene; sum the rates; divide each rate by the sum and multiply by 10^6.
2. **FPKM:** return `X_i · 10^9 / (l_i · N)`, guarding non-positive length/depth.
3. **Quantile normalization:** compute the across-sample mean at each rank; for each sample, walk its values in sorted order, group equal-valued (tied) runs, and assign each run the average of the rank means it spans, placed back at original positions.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- TPM scaling constant `10^6` and FPKM scaling constant `10^9` are named constants citing Sources [1][2][3].
- Tie handling in QN follows the Bolstad tie-average rule [4]; the suffix tree is not applicable (no substring search) — N/A for search reuse.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateTPM | O(n) | O(n) | n genes (one pass for rates, one for output) |
| CalculateFPKM | O(1) | O(1) | scalar |
| QuantileNormalize | O(s · n log n) | O(s · n) | s samples of n values; per-sample sort dominates |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [TranscriptomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs)

- `TranscriptomeAnalyzer.CalculateTPM(...)`: per-gene TPM (and FPKM) from raw counts.
- `TranscriptomeAnalyzer.CalculateFPKM(rawCount, length, totalReads)`: single-gene FPKM.
- `TranscriptomeAnalyzer.QuantileNormalize(samples)`: rank-based cross-sample normalization with tie averaging.

### 5.2 Current Behavior

TPM also fills the `FPKM` field of each `GeneExpression` using the sample's total raw count as N. QuantileNormalize re-sorts each sample twice (once for rank means, once for placement) for clarity rather than caching; inputs are materialized to lists. No substring/pattern search is involved, so the repository suffix tree is not used (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- TPM `= (X_i/l_i)/Σ(X_j/l_j)·10^6` exactly as in [2][3].
- FPKM `= X_i·10^9/(l_i·N)` exactly as in [2][3].
- Quantile normalization rank-mean transform with the tie-average rule, per [4].

**Intentionally simplified:**

- Effective length: uses annotated length (l̃_i = l_i); **consequence:** no fragment-length/bias correction — within-sample ratios match the source, very-short-transcript absolutes differ slightly (ASM-01).

**Not implemented:**

- Library-size scaling factors (TMM, DESeq median-of-ratios); **users should rely on:** dedicated differential-expression tooling (out of scope for this unit).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Pre-fix QN ignored tied ranks | Deviation | tied values got wrong rank means | fixed | corrected to Bolstad tie-average rule [4] |
| 2 | All-zero TPM denominator → 0 | Assumption | degenerate 0/0 input | accepted | ASSUMPTION-01; no source specifies; convention emits 0 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty geneCounts / samples | empty sequence | nothing to quantify |
| All-zero counts (TPM) | all TPM = 0 | Σ(X/l)=0 ⇒ 0/0 undefined; convention 0 |
| length ≤ 0 or N ≤ 0 (FPKM) | 0 | formula undefined |
| Tied values within a sample (QN) | averaged rank means | Bolstad tie rule [4] |
| Identical columns (QN) | output equals input | mean of equal values is the value |

### 6.2 Limitations

TPM and FPKM are within-sample relative measures and are misused when compared directly across samples or protocols [1][2]. No effective-length or GC/bias correction is applied. Quantile normalization assumes a shared underlying distribution and can mask genuine global shifts (ASM-02).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var tpm = TranscriptomeAnalyzer.CalculateTPM(new[]
{
    ("A", 10.0, 2000), ("B", 20.0, 4000), ("C", 30.0, 1000)
}).ToList();
// TPM: A=125000, B=125000, C=750000 (sum = 1,000,000)

double fpkm = TranscriptomeAnalyzer.CalculateFPKM(1000, 2000, 1_000_000); // 500
```

**Numerical walk-through (TPM):** RPK = (0.005, 0.005, 0.030), ΣRPK = 0.04; TPM = RPK/0.04·10^6 = (125000, 125000, 750000); total 10^6 (INV-01).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs) — covers INV-01..INV-05
- Evidence: [TRANS-EXPR-001-Evidence.md](../../../docs/Evidence/TRANS-EXPR-001-Evidence.md)

## 8. References

1. Wagner GP, Kin K, Lynch VJ. 2012. Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Theory in Biosciences 131(4):281–285. https://doi.org/10.1007/s12064-012-0162-3
2. Zhao S, Ye Z, Stanton R. 2020. Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. RNA 26(8):903–909. https://pmc.ncbi.nlm.nih.gov/articles/PMC7373998/
3. Pimentel H. 2014. What the FPKM? A review of RNA-Seq expression units. https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
4. Bolstad BM, Irizarry RA, Astrand M, Speed TP. 2003. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics 19(2):185–193. (worked example via) https://en.wikipedia.org/wiki/Quantile_normalization
5. Mortazavi A, Williams BA, McCue K, Schaeffer L, Wold B. 2008. Mapping and quantifying mammalian transcriptomes by RNA-Seq. Nature Methods 5(7):621–628. https://doi.org/10.1038/nmeth.1226
