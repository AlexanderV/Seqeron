# Loss of Heterozygosity (HRD-LOH)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-LOH-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Loss of heterozygosity (LOH) occurs when a heterozygous germline locus retains only one parental allele in the tumour: the minor-allele copy number drops to zero while the major allele is retained. This unit detects LOH regions from allele-specific copy-number segments and computes the **HRD-LOH genomic-scar score** — the number of LOH regions longer than 15 Mb that do not span a whole chromosome [1][2]. That count is one of the three components of the composite HRD score (with TAI and LST). The algorithm is specification-driven and exact: it reproduces the segment-counting rule of the scarHRD reference implementation [2], itself based on Abkevich et al. (2012) [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In a normal diploid genome a heterozygous SNP has both parental alleles (allele-specific copy numbers 1|1). Homologous-recombination-deficient (HRD) tumours accumulate genomic scars; one scar type is intermediate-size LOH. Abkevich et al. (2012) found that LOH regions of "intermediate size" (longer than 15 Mb but shorter than a whole chromosome) are strongly associated with BRCA1/2 deficiency, whereas whole-chromosome LOH is not — it is thought to arise through a competing mechanism not involving double-strand breaks [1].

### 2.2 Core Model

Given allele-specific copy-number segments, each with major-allele copy number `cn_major` and minor-allele copy number `cn_minor`:

- A segment is **LOH** iff `cn_minor == 0 AND cn_major != 0` (the minor allele is lost while the major allele is retained; a homozygous deletion `cn_major == 0` is excluded) [2].
- A chromosome is **whole-chromosome LOH** iff every one of its segments has `cn_minor == 0`; its regions are excluded [1][2].
- The **HRD-LOH score** is the number of LOH regions whose length (`end − start`) is strictly greater than `15,000,000` bp, after excluding whole-chromosome-LOH chromosomes [1][2][3].

scarHRD `calc.hrd` realises this as: `segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,]`; `segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,]` with `sizelimit1 = sizelimitLOH = 15e6`; `segLOH <- segLOH[!segLOH[,2] %in% chrDel,]`; `output[i] <- nrow(segLOH)` [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | HRD-LOH score ≥ 0 | it is a count of regions (`nrow(segLOH) ≥ 0`) [2] |
| INV-02 | 0 ≤ LOH_fraction ≤ 1 per chromosome | LOH length is a subset of total covered length on the chromosome |
| INV-03 | A counted LOH region has `cn_minor == 0` and `cn_major != 0` | LOH definition [2] |
| INV-04 | A counted LOH region has length strictly > 15,000,000 bp | size filter `> sizelimit1`, `sizelimitLOH = 15e6` [2] |
| INV-05 | No counted region lies on a whole-chromosome-LOH chromosome | `chrDel` exclusion [2]; Abkevich "< whole chromosome" [1] |
| INV-06 | The score is independent of input segment order | per-chromosome aggregation is set-based |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments | `IEnumerable<AlleleSpecificSegment>` | required | Allele-specific CN segments (chr, start, end, major CN, minor CN) | not null; each `End > Start`; copy numbers ≥ 0 |
| chromosome | `string` | required (fraction only) | Chromosome to score for `CalculateLOHFraction` | not null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `LohResult.Regions` | `IReadOnlyList<LohRegion>` | The qualifying LOH regions (> 15 Mb, not whole-chromosome) |
| `LohResult.Score` | `int` | HRD-LOH score = number of qualifying regions |
| `CalculateHrdLohScore` | `int` | The HRD-LOH score directly |
| `CalculateLOHFraction` | `double` | Length-weighted LOH fraction of one chromosome ∈ [0,1] |

### 3.3 Preconditions and Validation

Coordinates are in base pairs; segment length = `End − Start` (per scarHRD `seg[,4]-seg[,3]`, no `+1`) [2]. Null `segments` (or null `chromosome`) → `ArgumentNullException`. A segment with `End ≤ Start` or a negative copy number → `ArgumentException`. `CalculateLOHFraction` for a chromosome absent from the input returns `0.0` (no covered length). Chromosome identifiers are matched with ordinal (case-sensitive) string comparison.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate and group segments by chromosome (first-seen order).
2. Skip chromosomes where every segment is LOH (`cn_minor == 0`) — whole-chromosome LOH.
3. Within each remaining chromosome, sort by start and merge adjacent same-LOH-state segments (gap ≤ 1 bp).
4. Keep merged LOH regions with length strictly > 15 Mb; the count is the HRD-LOH score.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Rule | Value | Source |
|------|-------|--------|
| LOH segment | `cn_minor == 0 && cn_major != 0` | scarHRD `calc.hrd` [2] |
| Min region length | `> 15,000,000 bp` (strict) | Abkevich 2012 [1]; scarHRD `sizelimitLOH=15e6` [2] |
| Whole-chromosome exclusion | all segments `cn_minor == 0` → excluded | scarHRD `chrDel` [2]; Abkevich [1] |
| Adjacent-segment merge gap | ≤ 1 bp, same LOH state | oncoscanR `score_loh` [3] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DetectLOH` | O(n log n) | O(n) | n = segment count; dominated by per-chromosome sort for merging |
| `CalculateLOHFraction` | O(n) | O(1) | single pass over segments |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectLOH(IEnumerable<AlleleSpecificSegment>)`: returns the qualifying LOH regions and the HRD-LOH score.
- `OncologyAnalyzer.CalculateHrdLohScore(IEnumerable<AlleleSpecificSegment>)`: returns the HRD-LOH score directly.
- `OncologyAnalyzer.CalculateLOHFraction(IEnumerable<AlleleSpecificSegment>, string)`: per-chromosome length-weighted LOH fraction ∈ [0,1].

### 5.2 Current Behavior

Detection operates on already-segmented allele-specific copy number (the scarHRD `seg`-table shape); upstream segmentation / B-allele-frequency modelling is out of scope (cf. ONCO-CNA-001). No substring/pattern search is involved, so the repository suffix tree is **not applicable** (the input is numeric segments, not sequence text). `CalculateLOHFraction` deliberately applies neither the 15 Mb size filter nor the whole-chromosome exclusion — it reports a raw per-chromosome LOH burden satisfying INV-02, distinct from the HRD-LOH count.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- LOH segment condition `cn_minor == 0 && cn_major != 0` [2].
- Strict 15 Mb length filter (`> 15,000,000 bp`) [1][2].
- Whole-chromosome-LOH exclusion [1][2].
- Adjacent same-state merge (≤ 1 bp gap) [3].

**Intentionally simplified:**

- Segment merging uses the LOH/non-LOH state and a 1 bp adjacency gap (oncoscanR rule [3]); scarHRD's `shrink.seg.ai.wrapper` caps major CN to 1 before merging on full allelic state — **consequence:** identical LOH-region counts for typical inputs, since only the LOH state drives whether a run is counted.

**Not implemented:**

- Computing allele-specific copy-number segments from raw reads / BAF / VCF — **users should rely on:** an upstream segmenter (e.g. ASCAT/Sequenza/FACETS) feeding `AlleleSpecificSegment`; tracked as ONCO-CNA-001.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Length exactly 15,000,000 bp | not counted | strict `>` filter [2] |
| Homozygous deletion (minor=0, major=0) | not counted | `cn_major != 0` clause [2] |
| Heterozygous retained (minor≠0) | not counted | `cn_minor == 0` clause [2] |
| Whole-chromosome LOH | not counted | `chrDel` exclusion [1][2] |
| Empty input | score 0; fraction 0 | empty domain |
| Null input | `ArgumentNullException` | input validation |
| `End ≤ Start` or negative CN | `ArgumentException` | invalid segment |

### 6.2 Limitations

Requires allele-specific copy-number segments as input; tumour purity / ploidy correction and the BAF→copy-number inference are assumed already applied upstream. Whole-chromosome detection is segment-coverage based (all provided segments on the chromosome are LOH); incomplete chromosome coverage in the input could under-call the whole-chromosome condition.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var segments = new[]
{
    new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 20_000_000, 1, 0),       // 20 Mb LOH → counted
    new OncologyAnalyzer.AlleleSpecificSegment("1", 20_000_000, 60_000_000, 1, 1), // het, not LOH
    new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 10_000_000, 2, 0),       // LOH but ≤ 15 Mb
};
int score = OncologyAnalyzer.CalculateHrdLohScore(segments); // 1
double frac = OncologyAnalyzer.CalculateLOHFraction(segments, "1"); // 20M / 60M = 0.3333...
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectLOH_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_DetectLOH_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [ONCO-LOH-001-Evidence.md](../../../docs/Evidence/ONCO-LOH-001-Evidence.md)
- Related algorithms: [HRD composite score](../Oncology/HRD_Score.md)

## 8. References

1. Abkevich V, Timms KM, Hennessy BT, et al. 2012. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107(10):1776–1782. https://doi.org/10.1038/bjc.2012.451 (PMID 23047548)
2. Sztupinszki Z, et al. scarHRD R package — `calc.hrd` / `scar_score`. https://github.com/sztup/scarHRD/blob/master/R/calc.hrd.R
3. Christinat Y. oncoscanR `score_loh` documentation. https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html
