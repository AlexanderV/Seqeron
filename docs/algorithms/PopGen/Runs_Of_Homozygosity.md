# Runs of Homozygosity (ROH) and Genomic Inbreeding (F_ROH)

| Field | Value |
|-------|-------|
| Algorithm Group | PopGen |
| Test Unit ID | POP-ROH-001 |
| Related Projects | Seqeron.Genomics.Population |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

A run of homozygosity (ROH) is a contiguous stretch of the genome over which an individual is homozygous at every (or nearly every) marker, arising when two copies of an identical-by-descent haplotype are inherited from a common ancestor [1]. `FindROH` detects such runs from per-SNP genotypes using the window-free *consecutive-runs* method [3], and `CalculateInbreedingFromROH` converts the detected runs into the genomic inbreeding coefficient F_ROH, the proportion of the autosomal genome lying in ROH [1]. The detection is a deterministic linear scan; F_ROH is an exact ratio. The method is used to estimate individual autozygosity / inbreeding directly from SNP genotypes without a pedigree.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Offspring of related parents inherit long homozygous tracts because the same ancestral haplotype reaches them through both lineages. The length distribution of ROH reflects the depth of relatedness (recent consanguinity → long ROH; ancient drift → many short ROH) [1]. Summed ROH length, normalized by the analyzable autosomal genome, estimates the inbreeding coefficient F [1].

### 2.2 Core Model

**ROH detection (consecutive-runs method).** SNPs are processed in ascending position order. A candidate run grows while it satisfies all of: at most `maxOppRun` opposite (heterozygous) genotypes, at most `maxMissRun` missing genotypes, and no inter-SNP gap exceeding `maxGap`; crossing any threshold terminates the run [3]. A terminated run is retained only if it contains at least `minSNP` homozygous SNPs and spans at least `minLengthBps` base pairs [2][3].

**Genomic inbreeding coefficient.** McQuillan et al. define [1]:

> F_ROH = ΣL_ROH / L_AUTO

where ΣL_ROH is the total length of an individual's runs of homozygosity above a chosen minimum length and L_AUTO is the length of the autosomal genome covered by SNPs (excluding centromeres) [1]. McQuillan used L_AUTO = 2,673,768 kb [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Homozygosity reflects identity-by-descent (autozygosity) | ROH from LD/population homozygosity inflate F_ROH; choose a length threshold that excludes short non-IBD tracts [1] |
| ASM-02 | Tolerated opposite genotypes are genotyping errors, not true heterozygosity | Too-permissive `maxOppRun` merges distinct runs; too-strict fragments true runs [3] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported run has SnpCount ≥ minSnps and End − Start ≥ minLength | both thresholds are checked before emission [2][3] |
| INV-02 | A reported run has ≤ maxHeterozygotes opposite genotypes and no inter-SNP gap > maxGap | run is closed as soon as either limit is crossed [3] |
| INV-03 | Runs are emitted in ascending Start order with Start ≤ End | SNPs are sorted ascending and a run ends at its last homozygous SNP |
| INV-04 | 0 ≤ F_ROH ≤ 1 when all ROH lie within [0, genomeLength] | F_ROH = ΣL_ROH / L_AUTO with ΣL_ROH ≤ L_AUTO [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| genotypes | `IEnumerable<(int Position, int Genotype)>` | required | per-SNP position (bp) and genotype (0=hom-ref, 1=het, 2=hom-alt) | not null; ordered internally |
| minSnps | `int` | 100 | minimum SNPs in a retained run (PLINK --homozyg-snp) [2] | ≥ 1 |
| minLength | `int` | 1,000,000 | minimum run length in bp (PLINK --homozyg-kb 1000) [2] | ≥ 0 |
| maxHeterozygotes | `int` | 1 | max opposite genotypes inside a run (Marras maxOppRun; PLINK --homozyg-window-het) [2][3] | ≥ 0 |
| maxGap | `int` | 1,000,000 | max bp gap between consecutive SNPs in a run (PLINK --homozyg-gap 1000) [2] | ≥ 0 |
| rohSegments | `IEnumerable<(int Start, int End)>` | required | half-open ROH intervals for F_ROH | not null |
| genomeLength | `int` | required | L_AUTO in bp | returns 0 if ≤ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| FindROH result | `IEnumerable<(int Start, int End, int SnpCount)>` | each run: first/last SNP positions (inclusive, bp) and SNP count, ascending by Start |
| CalculateInbreedingFromROH result | `double` | F_ROH = ΣL_ROH / L_AUTO ∈ [0, 1] for in-range input |

### 3.3 Preconditions and Validation

Positions are 0-based bp; runs are reported on the inclusive `[Start, End]` SNP positions. Genotype `1` is the opposite (heterozygous) call; any other integer is treated as homozygous (no missing-data sentinel). `FindROH` throws `ArgumentNullException` for null genotypes and `ArgumentOutOfRangeException` for `minSnps < 1` or negative `minLength`/`maxHeterozygotes`/`maxGap`. `CalculateInbreedingFromROH` throws `ArgumentNullException` for null segments and returns 0 for `genomeLength ≤ 0` (no defined denominator). F_ROH segment lengths are computed as `End − Start` (half-open).

## 4. Algorithm

### 4.1 High-Level Steps

1. Sort genotypes by position ascending.
2. Scan SNP by SNP, growing the current run; track opposite-genotype count and the last homozygous SNP.
3. On a gap > maxGap or one het beyond `maxHeterozygotes`, close the run at its last homozygous SNP and restart at the breaking SNP.
4. Emit a closed run only if SnpCount ≥ minSnps and (End − Start) ≥ minLength.
5. F_ROH: sum retained-run lengths and divide by L_AUTO.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Default thresholds are PLINK 1.9 `--homozyg` defaults: `--homozyg-snp 100`, `--homozyg-kb 1000`, `--homozyg-window-het 1`, `--homozyg-gap 1000 kb` [2]. F_ROH uses L_AUTO supplied by the caller; McQuillan's value was 2,673,768 kb [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindROH | O(n log n) | O(n) | dominated by the position sort of n SNPs; the scan itself is O(n) (registry O(n) assumes pre-sorted input) |
| CalculateInbreedingFromROH | O(m) | O(1) | m = number of ROH segments |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.FindROH(...)`: window-free consecutive-runs ROH detection.
- `PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(...)`: F_ROH = ΣL_ROH / L_AUTO.

### 5.2 Current Behavior

A run is reported on the inclusive interval `[first homozygous SNP .. last homozygous SNP]`: trailing tolerated heterozygotes are not part of the emitted run, and a heterozygous SNP never seeds a new run. Argument validation is performed eagerly (before iteration begins) by wrapping the deferred `yield` iterator. No suffix tree or other substring-search structure is used: ROH detection is a single linear pass over already-positional SNP data, not a "find occurrences of X in Y" search, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- F_ROH = ΣL_ROH / L_AUTO exactly as McQuillan et al. (2008) [1].
- Consecutive-runs SNP-by-SNP scan with opposite-genotype tolerance and gap break (Marras et al. 2015) [3].
- PLINK 1.9 default thresholds for minimum SNP count, minimum length, max heterozygotes, and max gap [2].

**Intentionally simplified:**

- PLINK's two-phase sliding-window scan (`--homozyg-window-snp`, `--homozyg-window-threshold`, `--homozyg-window-missing`, `--homozyg-density`); **consequence:** ROH boundaries may differ slightly from PLINK when genotyping noise is window-dependent. The implemented consecutive-runs method is itself an established alternative [3].

**Not implemented:**

- Missing-genotype tolerance (`maxMissRun` / `--homozyg-window-missing`); **users should rely on:** pre-filtering missing calls upstream — the input model carries no missing sentinel.
- SNP-density filter (`--homozyg-density`); **users should rely on:** the `maxGap` constraint, which bounds the maximum inter-SNP spacing.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Window-free vs PLINK sliding window | Deviation | minor ROH-boundary differences vs PLINK | accepted | Marras 2015 reference method; ASM-02 |
| 2 | No missing-genotype model | Assumption | missing calls not tolerated as a separate class | accepted | input has no missing sentinel |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty genotype input | no runs | nothing to scan |
| Unsorted input | sorted internally; same runs | positional algorithm [3] |
| Leading heterozygotes | run starts at first homozygous SNP | a het cannot seed a homozygous run [3] |
| Run passing count but not length (or vice-versa) | discarded | both thresholds required [2] |
| genomeLength ≤ 0 | F_ROH = 0 | no defined denominator |

### 6.2 Limitations

Missing genotypes and PLINK's window/density refinements are not modeled (§5.3). The choice of `minLength` controls the IBD vs background-homozygosity trade-off (ASM-01); the default 1 Mb follows PLINK but a 1.5 Mb threshold was more discriminating for endogamy in McQuillan et al. [1]. The complexity is O(n log n) due to the internal sort.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var snps = Enumerable.Range(0, 100).Select(i => (i * 20_000, 0)).ToList();
var runs = PopulationGeneticsAnalyzer.FindROH(snps).ToList();
// runs[0] == (Start: 0, End: 1_980_000, SnpCount: 100)

double fRoh = PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(
    new[] { (0, 10_000_000), (50_000_000, 60_000_000) }, 100_000_000);
// fRoh == 0.20  (ΣL_roh = 20 Mb, L_auto = 100 Mb)  [McQuillan 2008]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PopulationGeneticsAnalyzer_FindROH_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_FindROH_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [POP-ROH-001-Evidence.md](../../../docs/Evidence/POP-ROH-001-Evidence.md)

## 8. References

1. McQuillan R, Leutenegger A-L, Abdel-Rahman R, et al. 2008. Runs of Homozygosity in European Populations. American Journal of Human Genetics 83(3):359-372. https://pmc.ncbi.nlm.nih.gov/articles/PMC2556426/ (DOI: 10.1016/j.ajhg.2008.08.007)
2. Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. 2015. Second-generation PLINK. GigaScience 4:7. PLINK 1.9 "Runs of homozygosity": https://www.cog-genomics.org/plink/1.9/ibd (DOI: 10.1186/s13742-015-0047-8)
3. Marras G, Gaspa G, Sorbolini S, et al. 2015. Analysis of runs of homozygosity and their relationship with inbreeding in five cattle breeds. Animal Genetics 46(2):110-121. detectRUNS vignette: https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html (DOI: 10.1111/age.12259)
