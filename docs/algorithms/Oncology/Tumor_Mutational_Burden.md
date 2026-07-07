# Tumor Mutational Burden (TMB)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-TMB-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Tumor mutational burden (TMB) is the number of somatic coding mutations per megabase of sequenced
genomic region, reported in mutations/megabase (mut/Mb) [1]. It is a predictive biomarker for response
to immune-checkpoint inhibitors: the FDA approved pembrolizumab for solid tumors with TMB ≥ 10 mut/Mb [2].
This unit computes TMB as a simple, exact ratio (`mutationCount / targetRegionMb`) and classifies the
result against the FDA TMB-High cutoff. The computation is deterministic and specification-driven; mutation
filtering (germline / known-driver removal) is performed upstream by the somatic caller (ONCO-SOMATIC-001).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A high somatic mutation count generates more neoantigens and correlates with immunotherapy response [1][3].
To make counts comparable across assays of different panel sizes, the raw mutation count is normalized to
the megabases of coding sequence interrogated, giving mut/Mb [1][3]. Foundation Medicine's comprehensive
genomic profiling counts "all base substitutions and indels in the coding region of targeted genes,
including synonymous alterations" before filtering out known driver and germline alterations [1].

### 2.2 Core Model

TMB is defined as the number of somatic, coding mutations per megabase of genome examined [1]:

```
TMB = mutationCount / targetRegionMb        (mut/Mb)
```

The FoundationOne 315-gene panel uses a denominator of **1.1 Mb** of coding genome [1]. Classification uses
the FDA-approved cutoff [2]:

```
TMB-High  ⇔  TMB ≥ 10 mut/Mb   (inclusive)
TMB-Low   ⇔  TMB < 10 mut/Mb
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `CalculateTMB(n, r) = n / r` for r > 0 | direct definition [1] |
| INV-02 | TMB ≥ 0 for n ≥ 0, r > 0 | quotient of non-negatives |
| INV-03 | TMB is non-decreasing in count (fixed r) and non-increasing in r (fixed count) | division monotonicity |
| INV-04 | `ClassifyTMB(tmb) = High ⇔ tmb ≥ 10` | FDA inclusive cutoff [2] |

### 2.5 Comparison with Related Methods

| Aspect | TMB (this unit) | MSI (ONCO-MSI-001) |
|--------|-----------------|--------------------|
| Quantity | somatic mutations per Mb | fraction of unstable microsatellite loci |
| Cutoff | FDA TMB-High ≥ 10 mut/Mb [2] | MSS / MSI-L / MSI-H |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| mutationCount | int | required | counted somatic coding mutations | ≥ 0 |
| targetRegionMb | double | required | sequenced coding region in megabases | finite, > 0 |
| calls | IEnumerable&lt;SomaticCall&gt; | required (overload) | classified somatic calls; only Somatic counted | non-null |
| tmb | double | required | TMB value to classify | finite, ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (CalculateTMB) | double | TMB in mut/Mb (≥ 0) |
| (ClassifyTMB) | TmbStatus | `High` (tmb ≥ 10) or `Low` (tmb < 10) |

### 3.3 Preconditions and Validation

`CalculateTMB(int, double)` throws `ArgumentOutOfRangeException` when `mutationCount < 0` or when
`targetRegionMb` is NaN, infinite, or ≤ 0 (TMB is undefined when the megabase denominator is 0). The
`SomaticCall` overload throws `ArgumentNullException` on a null collection, then delegates. `ClassifyTMB`
throws `ArgumentOutOfRangeException` for negative or non-finite TMB. Counts are integers; the region is a
positive real in megabases.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `mutationCount ≥ 0` and `targetRegionMb` finite and > 0.
2. Return `mutationCount / targetRegionMb`.
3. (Overload) Count `Somatic`-status calls, then apply steps 1–2.
4. `ClassifyTMB`: return `High` if `tmb ≥ 10`, else `Low`.

### 4.2 Decision Rules, Scoring, Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| TMB-High cutoff | 10.0 mut/Mb (inclusive ≥) | FDA pembrolizumab approval / FoundationOne CDx [2] |
| Reference panel denominator (example) | 1.1 Mb (315-gene panel) | Chalmers 2017 [1] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateTMB(int, double) | O(1) | O(1) | single division |
| CalculateTMB(calls, double) | O(n) | O(1) | one pass to count Somatic calls |
| ClassifyTMB | O(1) | O(1) | single comparison |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CalculateTMB(int, double)`: TMB = count / Mb.
- `OncologyAnalyzer.CalculateTMB(IEnumerable<SomaticCall>, double)`: counts Somatic calls, delegates.
- `OncologyAnalyzer.ClassifyTMB(double)`: TMB-High when tmb ≥ 10.

### 5.2 Current Behavior

This is not a search/matching unit, so the repository suffix tree is not applicable (N/A) — the algorithm
is an O(1) arithmetic ratio and an O(n) count, with no substring/occurrence search. Mutation filtering
(germline, known-driver) is done upstream and not duplicated here; the count passed in is the already-
filtered somatic count, consistent with Chalmers 2017's pre-count filtering [1].

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- TMB = somatic mutations per megabase of examined region (`count / Mb`) [1].
- TMB-High classification at the FDA-approved inclusive cutoff of ≥ 10 mut/Mb [2].

**Intentionally simplified:**

- Two-tier (High/Low) classification only; **consequence:** no Low/Intermediate/High sub-banding — the
  Registry's 6/20 boundaries have no retrieved authoritative source, so they are not implemented rather than
  fabricated.

**Not implemented:**

- Mutation filtering (COSMIC drivers, TSG truncations, dbSNP/ExAC germline, somatic-germline-zygosity) [1];
  **users should rely on:** the upstream somatic caller (ONCO-SOMATIC-001) to produce the filtered count.
- Panel-size / coding-region computation from a BED/target file; **users should rely on:** supplying the
  Mb denominator directly.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Registry "Low <6 / Intermediate 6–20 / High >20" vs FDA ≥10 | Deviation | classification boundary | accepted | No source retrieved for 6/20; FDA ≥10 used (evidence-first). Checklist updated. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| targetRegionMb = 0 | throws ArgumentOutOfRangeException | division by zero; TMB undefined [1] |
| mutationCount = 0 | TMB = 0 | quotient of 0 |
| tmb = 10.0 | High | inclusive cutoff [2] |
| Panel < 0.5 Mb | value still computed (no throw) | mathematically defined; instability documented, not an error [1] |

### 6.2 Limitations

Below ~0.5 Mb of sequenced region, panel-based TMB deviates substantially from whole-exome TMB and is
unreliable as a clinical estimate [1] — the value is returned but should be interpreted with caution. The
unit does not validate that the supplied count was correctly filtered; it trusts the upstream caller.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double tmb = OncologyAnalyzer.CalculateTMB(mutationCount: 11, targetRegionMb: 1.1); // 10.0 mut/Mb
OncologyAnalyzer.TmbStatus status = OncologyAnalyzer.ClassifyTMB(tmb);              // High (>= 10)
```

**Numerical walk-through:** 11 somatic mutations over the 315-gene 1.1 Mb panel [1] gives
11 / 1.1 = 10.0 mut/Mb, which is exactly the FDA TMB-High cutoff [2], so the tumor is TMB-High.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CalculateTMB_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_CalculateTMB_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-TMB-001-Evidence.md](../../../docs/Evidence/ONCO-TMB-001-Evidence.md)

## 8. References

1. Chalmers ZR, Connelly CF, Fabrizio D, et al. 2017. Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Genome Medicine 9:34. https://doi.org/10.1186/s13073-017-0424-2 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719)
2. Marcus L, Fashoyin-Aje LA, Donoghue M, et al. 2021. FDA Approval Summary: Pembrolizumab for the Treatment of Tumor Mutational Burden–High Solid Tumors. Clinical Cancer Research 27(17):4685–4689. https://doi.org/10.1158/1078-0432.CCR-21-0327 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8416776/)
3. Fancello L, Gandini S, Pelicci PG, Mazzarella L; Friends of Cancer Research TMB Harmonization Project (Merino DM, McShane LM, et al. 2020, J Immunother Cancer 8:e000147). Tumor Mutational Burden as a Predictive Biomarker in Solid Tumors. https://pmc.ncbi.nlm.nih.gov/articles/PMC7710563/
