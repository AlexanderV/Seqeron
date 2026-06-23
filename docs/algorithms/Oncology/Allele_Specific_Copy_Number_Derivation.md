# Allele-Specific Copy-Number Derivation (ASCAT-style purity/ploidy fit + multiplicity)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ASCAT-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-23 |

## 1. Overview

This unit derives the upstream quantities that the tumour purity/ploidy/CCF/clonality methods consume, so the
pipeline can run from per-locus signal rather than caller-supplied pre-made allele-specific segments and
multiplicity. From per-locus log-R ratio (logR) and B-allele frequency (BAF) at germline-heterozygous SNPs
(observed measurements), it (1) **segments** the genome into (logR, BAF) summaries, (2) **jointly fits** tumour
purity ρ and ploidy ψ by grid search using the ASCAT equations and goodness-of-fit objective, emitting
allele-specific **integer** copy-number segments, and (3) **derives** somatic mutation multiplicity from VAF,
purity and copy number. It is a faithful but simplified single-sample realisation of ASCAT (Van Loo et al. 2010)
[1][2] plus the McGranahan multiplicity convention [3][4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A bulk tumour sample is a mixture of aberrant (cancer) cells at fraction ρ (purity) and normal cells at 1 − ρ.
SNP/sequencing assays yield two tracks per locus: **logR** (total signal intensity) and **BAF** (allelic
contrast at germline-het SNPs) [1]. To recover allele-specific integer copy numbers one must jointly estimate ρ
and the tumour ploidy ψ, because the same observed (logR, BAF) is consistent with multiple (ρ, ψ) solutions [1].

### 2.2 Core Model

For a segment with fitted logR r and BAF b, ASCAT maps to allele-specific copy numbers (nA, nB) given ρ, ψ and
the platform parameter γ (verbatim from ascat.runAscat.R [2]):

```
nA = (ρ − 1 − (b − 1)·2^(r/γ) · ((1−ρ)·2 + ρ·ψ)) / ρ
nB = (ρ − 1 +  b   ·2^(r/γ) · ((1−ρ)·2 + ρ·ψ)) / ρ
```

For massively parallel sequencing data **γ = 1** [2]. The joint fit grid-searches (ρ, ψ) to minimise the
segment-length-weighted squared distance of the minor allele to the nearest non-negative integer [1][2]:

```
d(ρ,ψ) = Σ_segments  (n_minor − round(n_minor))² · length · w_b      where w_b = 0.05 if b = 0.5 else 1
TheoretMaxdist = Σ_segments 0.25 · length · w_b
goodnessOfFit  = (1 − d / TheoretMaxdist) · 100   [%]
```

(0.25 = (½)² is the worst-case squared distance to an integer [2].)

Mutation multiplicity (number of mutated copies per cancer cell) is the rounded observed mutation copy number
[3][4]:

```
n_mut = VAF · (1/ρ) · [ρ·N_T + 2(1−ρ)]
m     = clamp( round(n_mut), 1, majorCopyNumber )
```

which is the inversion of the PICTograph generative model VAF = m·CCF·ρ / (N_T·ρ + 2(1−ρ)) at clonal CCF = 1 [4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A single tumour clone / single (ρ, ψ) explains the genome | Subclonal copy number is not modelled; segments forced to nearest integer CN |
| ASM-02 | γ matches the platform (γ = 1 for sequencing) | Wrong γ rescales logR → biased copy number [2] |
| ASM-03 | BAF is measured at germline-heterozygous SNPs | BAF at homozygous sites is uninformative; folding assumes het loci |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | At the true (ρ₀, ψ₀) of an integer-CN genome, d ≈ 0 and GoF ≈ 100% | nA, nB equal exact integers by construction [1][2] |
| INV-02 | Derived multiplicity ∈ [1, majorCopyNumber] | explicit clamp; a variant sits on ≥ 1 and ≤ major-allele copies [3][4] |
| INV-03 | GoF ≤ 100% | d ≥ 0 and d ≤ TheoretMaxdist by the 0.25 worst-case bound [2] |
| INV-04 | major ≥ minor in every emitted segment | nA, nB sorted before rounding |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loci | IEnumerable\<AlleleSpecificLocus\> | required | per-locus (chrom, pos, logR, BAF) measurements | non-null; chrom non-null |
| logRChangeThreshold | double | required | mean-shift split threshold (logR units) | > 0 |
| minLociPerSegment | int | 1 | min loci before a split | ≥ 1 |
| segments | IReadOnlyList\<AlleleSpecificSegmentSummary\> | required | segment summaries for the fit | non-empty |
| purityMin/Max/Step | double | 0.05 / 1.0 / 0.01 | purity grid | (0,1]; max ≥ min; step > 0 |
| ploidyMin/Max/Step | double | 1.5 / 5.0 / 0.05 | ploidy grid | > 0; max ≥ min; step > 0 |
| gamma | double | 1.0 | platform γ | > 0 |
| vaf, purity, totalCopyNumber, majorCopyNumber | double/int | required | multiplicity inputs | vaf∈[0,1]; ρ∈(0,1]; N_T≥1; major∈[1,N_T] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| AlleleSpecificSegmentSummary | record | per-segment mean logR, mirrored mean BAF, locus count |
| PurityPloidyFit | record | recovered ρ, ψ, GoF %, and the implied integer `AlleleSpecificSegment`s |
| DeriveMultiplicity | int | integer multiplicity m ∈ [1, majorCopyNumber] |

### 3.3 Preconditions and Validation

Positions are 0-based. Null `loci`/`segments` → `ArgumentNullException`; empty `segments` → `ArgumentException`;
out-of-range thresholds, grid bounds, or multiplicity arguments → `ArgumentOutOfRangeException`. BAF is mirrored
about 0.5 (b' = 0.5 + |b − 0.5|) during segmentation so the two symmetric het clusters reinforce.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Segment:** scan loci in order; start a new segment on chromosome change, on a logR mean-shift, OR on a
   (mirrored) BAF mean-shift (after `minLociPerSegment`). Summarise each run by mean logR and mirrored mean BAF.
   Segmenting on BAF as well as logR is essential because copy-neutral LOH (e.g. 2:0) shares a balanced region's
   logR but not its BAF.
2. **Fit:** for each (ρ, ψ) on the grid, map every segment to (nA, nB) via the ASCAT equations, accumulate the
   length-weighted squared minor-allele integer distance (the reported GoF objective). For *selection* the
   major-allele integer distance is added too, and exact ties prefer the lower ploidy ψ — the ASCAT parsimony
   convention that resolves the 2n vs 4n degeneracy. Emit the rounded, clamped integer `AlleleSpecificSegment`s
   and the percentage GoF (computed from the minor-allele distance, per ascat.runAscat.R).
3. **Multiplicity:** m = clamp(round(VAF·[ρ·N_T + 2(1−ρ)]/ρ), 1, major). Feed (VAF, ρ, N_T, m) into `EstimateCcf`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- γ = 1 (sequencing) [2]; balanced (BAF = 0.5) segments weighted ×0.05 [2]; worst-case integer distance 0.25 [2].
- Segments with End == Start are emitted with a 1 bp span so `AlleleSpecificSegment.Length > 0`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Segmentation | O(L) | O(L) | L = loci |
| Purity/ploidy fit | O(P·Q·S) | O(S) | P,Q = grid sizes, S = segments |
| Multiplicity | O(1) | O(1) | closed form |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.SegmentAlleleSpecific(...)`: mean-shift segmentation of per-locus logR/BAF.
- `OncologyAnalyzer.FitPurityPloidy(...)`: ASCAT grid fit → ρ, ψ, GoF, integer segments.
- `OncologyAnalyzer.DeriveMultiplicity(...)`: McGranahan multiplicity (rounded, clamped).

### 5.2 Current Behavior

Single-sample, single-clone fit on a fixed (ρ, ψ) grid; minor-allele integer-distance objective; mirrored-BAF
summaries. Not a search/matching task, so the repository suffix tree is **not used** (no occurrence enumeration).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ASCAT nA/nB equations and the minor-allele integer-distance goodness-of-fit, including the 0.05 balanced
  down-weight and the 0.25 worst-case term [2].
- McGranahan observed mutation copy number n_mut and the [1, major] multiplicity clamp [3][4].

**Intentionally simplified:**

- Segmentation: a joint logR + mirrored-BAF mean-shift (PCF/CBS-style), **not** the full ASPCF penalised
  least-squares multi-track segmentation; **consequence:** breakpoints are detected greedily left-to-right rather
  than by a global penalised fit.
- Fit: fixed-grid minimum with a lower-ploidy tie-break, **not** the ASCAT refit/sub-clonal heuristics;
  **consequence:** sub-clonal copy number is rounded to the nearest integer state.

**Not implemented:**

- Multi-sample / sub-clonal copy number and the full ASCAT refitting; **users should rely on:** ASCAT/FACETS for
  production allele-specific calling; this unit covers the canonical single-sample derivation feeding the
  downstream `EstimatePurity`/`EstimatePloidy`/`EstimateCcf`/`ClassifyClonality`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | logR mean-shift segmentation (not ASPCF) | Deviation | breakpoint sensitivity | accepted | ASM-01; use ASCAT for production segmentation |
| 2 | Single (ρ, ψ) per genome | Assumption | no subclonal CN | accepted | ASM-01 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Balanced-only genome (all b=0.5) | fit completes; segments ×0.05 weighted | source [2] |
| Single locus per chromosome | one segment per chromosome, LocusCount=1 | segmentation contract |
| VAF rounds to 0 | multiplicity clamped to 1 | INV-02 |
| VAF rounds above major CN | multiplicity clamped to major CN | INV-02 |

### 6.2 Limitations

logR and BAF are observed measurements and are always a caller input — this is inherent, not a limitation of the
derivation. The derivation does not model sub-clonal copy number, whole-genome-doubling disambiguation beyond the
grid, or contamination, and assumes germline-het BAF loci.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var loci = /* per-locus (chrom, pos, logR, BAF) measurements */;
var summaries = OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2);
var fit = OncologyAnalyzer.FitPurityPloidy(summaries);          // → ρ, ψ, integer segments
double ploidy = OncologyAnalyzer.EstimatePloidy(fit.Segments);  // downstream consumer
var seg = fit.Segments[0];
int m = OncologyAnalyzer.DeriveMultiplicity(vaf: 0.40, purity: fit.Purity,
            totalCopyNumber: seg.MajorCopyNumber + seg.MinorCopyNumber,
            majorCopyNumber: seg.MajorCopyNumber);
var ccf = OncologyAnalyzer.EstimateCcf(0.40, fit.Purity,
            seg.MajorCopyNumber + seg.MinorCopyNumber, m);       // ≈ 1.0 for a clonal mutation
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_AscatDerivation_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AscatDerivation_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-ASCAT-001-Evidence.md](../../../docs/Evidence/ONCO-ASCAT-001-Evidence.md)
- Related algorithms: [Tumor_Ploidy_Estimation](./Tumor_Ploidy_Estimation.md), [Cancer_Cell_Fraction_Estimation](./Cancer_Cell_Fraction_Estimation.md), [Tumor_Purity_Estimation](./Tumor_Purity_Estimation.md)

## 8. References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R`. https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ. 2021. DeCiFering the elusive cancer cell fraction. PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
