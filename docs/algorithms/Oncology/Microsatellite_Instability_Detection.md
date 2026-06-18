# Microsatellite Instability (MSI) Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-MSI-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Microsatellite instability (MSI) is a hypermutable phenotype caused by defective DNA mismatch repair, used clinically to predict response to immune-checkpoint blockade and to screen for Lynch syndrome. This unit implements the **scoring-and-classification layer** of MSI detection: given per-locus stability calls, it computes the MSI score as the fraction of unstable microsatellite loci and classifies the sample as MSI-High or microsatellite-stable using the MSIsensor2 20% cutoff [1][2], and provides the categorical NCI/Bethesda marker-count classification (MSS / MSI-L / MSI-H) [3]. It is specification/threshold-driven and deterministic. The upstream per-locus instability call (chi-square comparison of tumor-vs-normal repeat-length distributions) is out of scope.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Microsatellites are short tandem repeats (e.g. mononucleotide runs such as BAT-25, BAT-26) whose length is destabilised when mismatch repair fails. MSIsensor-family tools build expected (normal) and observed (tumor) repeat-length distributions per locus, test each locus, and summarise the sample by the fraction of loci that are unstable [1]. The classic clinical assay (NCI/Bethesda) scores instability across a validated 5-marker panel [3].

### 2.2 Core Model

**Continuous MSI score.** For a sample with `u` unstable loci among `n` valid evaluated loci [1][2]:

```
MSI score = u / n          (reported as a fraction in [0,1]; ×100 for percent)
```

A site is called unstable in MSIsensor when the chi-square comparison of tumor vs normal length distributions is significant under a default FDR of 0.05 [1] (upstream, out of scope here).

**Continuous classification.** MSIsensor2 [2]: a sample is **MSI-High** when `MSI score ≥ 20%` (inclusive); otherwise stable.

**Categorical (Bethesda) classification.** Over a fixed reference panel of markers, the NCI workshop [3] defines: **MSI-H** if ≥ 2 of 5 markers are unstable; **MSI-L** if exactly 1 of 5 is unstable; **MSS** if none is unstable.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ MSI score ≤ 1 | score = u/n with 0 ≤ u ≤ n [1][2] |
| INV-02 | score = u/n exactly | direct definition [2] |
| INV-03 | continuous status = MSI-High iff score ≥ 0.20 | MSIsensor2 cutoff [2] |
| INV-04 | Bethesda: 0→MSS, 1→MSI-L, ≥2→MSI-H | Boland et al. 1998 [3] |

### 2.5 Comparison with Related Methods

| Aspect | Computational (MSIsensor2) | Clinical (Bethesda panel) |
|--------|---------------------------|---------------------------|
| Input | fraction of unstable loci (NGS, many loci) | discrete count over 5 validated markers |
| MSI-H rule | score ≥ 20% [2] | ≥ 2/5 markers unstable [3] |
| MSI-L band | not defined on the continuous score | exactly 1/5 marker [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| unstableLoci | int | required | unstable loci count | 0 ≤ unstableLoci ≤ totalLoci |
| totalLoci | int | required | valid evaluated loci | > 0 |
| score | double | required | MSI score (fraction) | finite, in [0,1] |
| unstableMarkers / totalMarkers | int | required | Bethesda marker counts | 0 ≤ unstable ≤ total; total > 0 |
| locusUnstableFlags | IEnumerable&lt;bool&gt; | required | per-locus stability (true=unstable) | non-null, non-empty |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (CalculateMSIScore) | double | MSI score in [0,1] |
| (ClassifyMSIStatus / ClassifyBethesdaPanel) | MsiStatus | MSS / MSI_Low / MSI_High |
| MsiResult | record struct | UnstableLoci, TotalLoci, Score, Status |

### 3.3 Preconditions and Validation

`CalculateMSIScore` throws `ArgumentOutOfRangeException` if `totalLoci ≤ 0`, `unstableLoci < 0`, or `unstableLoci > totalLoci`. `ClassifyMSIStatus` throws if `score` is non-finite or outside [0,1]. `ClassifyBethesdaPanel` throws if `totalMarkers ≤ 0`, `unstableMarkers < 0`, or `unstableMarkers > totalMarkers`. `DetectMSI` throws `ArgumentNullException` on null and `ArgumentOutOfRangeException` on an empty sequence (no valid loci). The 20% continuous cutoff is inclusive.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count unstable loci `u` and valid loci `n` (or accept them directly).
2. Compute MSI score = `u / n` [1][2].
3. Continuous classification: MSI-High iff score ≥ 0.20 [2].
4. Categorical classification (Bethesda): ≥2→MSI-H, 1→MSI-L, 0→MSS [3].

### 4.2 Decision Rules, Scoring, Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| MsiHighScoreThreshold | 0.20 (≥, inclusive) | niu-lab/msisensor2 README [2] |
| BethesdaMsiHighMarkerCount | 2 | Boland et al. 1998 [3] |
| BethesdaMsiLowMarkerCount | 1 | Boland et al. 1998 [3] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateMSIScore / Classify* | O(1) | O(1) | arithmetic / comparison |
| DetectMSI | O(n) | O(1) | single pass over n loci flags |

This is not a search/matching unit, so the repository suffix tree is **not applicable** (no substring search, pattern matching, or occurrence enumeration is performed).

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateMSIScore(int, int)`: fraction unstable/total.
- `OncologyAnalyzer.ClassifyMSIStatus(double)`: continuous MSI-H (≥20%) vs MSS.
- `OncologyAnalyzer.ClassifyBethesdaPanel(int, int)`: categorical MSS/MSI-L/MSI-H.
- `OncologyAnalyzer.DetectMSI(IEnumerable<bool>)`: end-to-end count → score → status.

### 5.2 Current Behavior

The continuous `ClassifyMSIStatus` returns only MSI_High or MSS (no MSI_Low band), because MSIsensor2 defines only a binary MSI-H cutoff on the continuous score [2]; MSI_Low is produced solely by `ClassifyBethesdaPanel` from a discrete marker count [3]. `DetectMSI` uses the MSIsensor2 continuous cutoff. No suffix tree is used (Section 4.3).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- MSI score = unstable loci / valid loci [1][2].
- Continuous MSI-H at score ≥ 20% (inclusive) [2].
- Bethesda categorical: ≥2/5→MSI-H, 1/5→MSI-L, 0→MSS [3].

**Intentionally simplified:**

- Per-locus instability calling (chi-square tumor-vs-normal length distributions, FDR 0.05): **approximation:** accepted as upstream input (per-locus boolean flags); **consequence:** the user must supply already-called locus stability, not raw BAM/read data [1].
- No MSI_Low band on the continuous score: **consequence:** sub-20% samples report MSS, matching MSIsensor2's binary cutoff [2].

**Not implemented:**

- Tumor-vs-normal distribution modelling / chi-square testing from reads; **users should rely on:** an upstream caller (MSIsensor/MSIsensor2) to produce per-locus flags.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Continuous score has no MSI-L band | Assumption | sub-threshold samples → MSS | accepted | MSIsensor2 defines only binary MSI-H [2]; MSI-L via Bethesda count [3] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| 0 valid loci | throws | score u/n undefined [2] |
| unstable > valid | throws | 0 ≤ u ≤ n [2] |
| score = 0.20 | MSI-H | inclusive cutoff [2] |
| score = 0.0 | MSS | < 20% [2] |
| 1 of 5 markers | MSI-L | Boland 1998 [3] |
| empty flags / null | throws | no valid loci / guard |

### 6.2 Limitations

MSS vs MSI-L is unreliable on small panels [3]. The unit does not perform per-locus calling from sequencing reads; it requires upstream stability calls. Only the MSIsensor2 20% cutoff and the Bethesda 5-marker rule are implemented; tumor-type-specific or alternative cut-points (e.g. MSIsensor's 3.5% cohort boundary [1]) are not.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
// 6 unstable of 20 valid loci → score 0.30 ≥ 0.20 → MSI-High (MSIsensor2).
var r = OncologyAnalyzer.DetectMSI(new[] {
    true,true,true,true,true,true,                 // 6 unstable
    false,false,false,false,false,false,false,
    false,false,false,false,false,false,false });  // 14 stable
// r.Score == 0.30, r.Status == MsiStatus.MSI_High

// Bethesda: 2 of 5 markers unstable → MSI-H.
var s = OncologyAnalyzer.ClassifyBethesdaPanel(2, 5); // MsiStatus.MSI_High
```

### 7.2 Applications and Use Cases

- **Immunotherapy selection:** MSI-H solid tumors are eligible for checkpoint-inhibitor therapy.
- **Lynch syndrome screening:** the Bethesda panel identifies tumors warranting germline MMR testing [3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectMSI_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectMSI_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-MSI-001-Evidence.md](../../../docs/Evidence/ONCO-MSI-001-Evidence.md)

## 8. References

1. Niu B, Ye K, Zhang Q, Lu C, Xie M, McLellan MD, Wendl MC, Ding L. 2014. MSIsensor: microsatellite instability detection using paired tumor-normal sequence data. Bioinformatics 30(7):1015–1016. https://doi.org/10.1093/bioinformatics/btt755
2. niu-lab. msisensor2 — Microsatellite instability (MSI) detection for tumor only data. README (accessed 2026-06-14). https://github.com/niu-lab/msisensor2
3. Boland CR, Thibodeau SN, Hamilton SR, Sidransky D, Eshleman JR, Burt RW, et al. 1998. A National Cancer Institute Workshop on Microsatellite Instability for cancer detection and familial predisposition: development of international criteria for the determination of microsatellite instability in colorectal cancer. Cancer Res 58(22):5248–5257. https://pubmed.ncbi.nlm.nih.gov/9823339/
