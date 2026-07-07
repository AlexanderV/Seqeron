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
[1][2] plus the McGranahan multiplicity convention [3][4]. In addition it provides the **ASPCF** penalised
least-squares segmentation (Nilsen et al. 2012 [6]; Ross et al. 2021 [7]) — the global-optimum joint logR/BAF
changepoint method ASCAT uses — and **sub-clonal copy number** modelling (Battenberg two-population model,
Nik-Zainal et al. 2012 [8]), which expresses a segment that does not fit a single integer state as a mixture of
two adjacent integer states with a sub-clonal cellular fraction.

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

**ASPCF segmentation (penalised least squares).** ASCAT segments via Piecewise Constant Fitting [6][7]: minimise

```
L(S | y, γ) = Σ_{I∈S} Σ_{j∈I} (y_j − ȳ_I)² + γ·|S|
```

over all segmentations S, where the first term is the within-segment sum of squared errors (SSE), ȳ_I is the
segment mean, |S| the number of segments, and γ > 0 a penalty per breakpoint [6]. The global optimum is found by
the dynamic-program recurrence (O(n²)) [6]:

```
e_k = min_{j ∈ {1,…,k}} ( d_{jk} + e_{j−1} + γ ),    e_0 = 0
```

with d_{jk} the within-segment SSE of loci j..k. For allele-specific (joint) segmentation the logR (y₁) and
mirrored-BAF (y₂) tracks are segmented with **common breakpoints** but separate per-track means [6][7]:

```
L(S | y₁, y₂, γ) = L(S | y₁, γ) + L(S | y₂, γ)
```

so the per-segment data cost is (logR-SSE + mirroredBAF-SSE) and γ is charged once per segment. BAF is mirrored
to a single allelic-imbalance track (b' = 0.5 + |b − 0.5|) before segmentation [7].

**Sub-clonal copy number (Battenberg two-population model).** A segment has either one integer state (clonal, all
tumour cells) or two integer states (sub-clonal, two cell populations whose fractions sum to 1) [8]. The
real-valued ASCAT allele-specific copy numbers (nA, nB) at the fitted (ρ, ψ) are decomposed, when not
(near-)integer, as a single shared fraction f mixing the two bracketing integers:

```
major_obs = f·a_hi + (1−f)·a_lo,   minor_obs = f·b_hi + (1−f)·b_lo,   f ∈ [0,1]
```

where (a_hi, b_hi)/(a_lo, b_lo) are the two integer states; f is solved by least squares over both candidate
allele pairings and the lower-residual pairing is kept. A near-integer segment collapses to a single state (f≈0/1)
[8].

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
| INV-05 | ASPCF penalised cost is the global minimum (≤ any greedy segmentation) | DP recurrence over all S [6] |
| INV-06 | A segment with no logR/BAF change is one ASPCF segment; γ→∞ ⇒ \|S\|=1 | penalty dominates SSE gains [6] |
| INV-07 | Sub-clonal state fractions sum to 1; a clonal segment has f=1 and no second state | two-population model [8] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loci | IEnumerable\<AlleleSpecificLocus\> | required | per-locus (chrom, pos, logR, BAF) measurements | non-null; chrom non-null |
| logRChangeThreshold | double | required | mean-shift split threshold (logR units) | > 0 |
| minLociPerSegment | int | 1 | min loci before a split | ≥ 1 |
| penalty (ASPCF γ) | double | 40.0 | per-segment penalty in the PCF cost | > 0 |
| purity, ploidy (sub-clonal) | double | required | fitted ρ, ψ for the sub-clonal decomposition | ρ∈(0,1]; ψ>0 |
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
| SubclonalSegmentFit | record | per-segment primary/secondary `SubclonalCopyNumberState` (major, minor, cellFraction) + IsSubclonal flag |

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
4. **ASPCF (alternative to step 1):** per chromosome, run the PCF dynamic program over the joint (logR, mirrored
   BAF) SSE with penalty γ; backtrack the global-optimum breakpoint set and emit segment summaries. This replaces
   the greedy mean-shift with the penalised-least-squares optimum [6][7].
5. **Sub-clonal fit:** for each segment compute (nA, nB) at the fitted (ρ, ψ); if both alleles snap to integers
   within tolerance the segment is clonal (one state, f=1), else decompose into two adjacent integer states with a
   shared least-squares fraction f (Battenberg) [8].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- γ = 1 (sequencing) [2]; balanced (BAF = 0.5) segments weighted ×0.05 [2]; worst-case integer distance 0.25 [2].
- Segments with End == Start are emitted with a 1 bp span so `AlleleSpecificSegment.Length > 0`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Segmentation (greedy) | O(L) | O(L) | L = loci |
| ASPCF segmentation | O(L²) per chromosome | O(L) | PCF dynamic program (Nilsen 2012) [6] |
| Purity/ploidy fit | O(P·Q·S) | O(S) | P,Q = grid sizes, S = segments |
| Multiplicity | O(1) | O(1) | closed form |
| Sub-clonal fit | O(S) | O(S) | closed-form decomposition per segment |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.SegmentAlleleSpecific(...)`: greedy mean-shift segmentation of per-locus logR/BAF.
- `OncologyAnalyzer.SegmentAlleleSpecificAspcf(...)`: ASPCF penalised-least-squares (PCF DP) joint segmentation.
- `OncologyAnalyzer.FitPurityPloidy(...)`: ASCAT grid fit → ρ, ψ, GoF, integer segments.
- `OncologyAnalyzer.DeriveMultiplicity(...)`: McGranahan multiplicity (rounded, clamped).
- `OncologyAnalyzer.FitSubclonalCopyNumber(...)`: Battenberg two-state sub-clonal copy-number decomposition.

### 5.2 Current Behavior

Single-sample fit on a fixed (ρ, ψ) grid; minor-allele integer-distance objective; mirrored-BAF summaries.
Segmentation is available both as the original greedy mean-shift and as the global-optimum ASPCF (PCF DP); the
sub-clonal fit adds the Battenberg two-population decomposition. Not a search/matching task, so the repository
suffix tree is **not used** (no occurrence enumeration).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ASCAT nA/nB equations and the minor-allele integer-distance goodness-of-fit, including the 0.05 balanced
  down-weight and the 0.25 worst-case term [2].
- McGranahan observed mutation copy number n_mut and the [1, major] multiplicity clamp [3][4].
- ASPCF penalised-least-squares segmentation: the PCF cost `Σ SSE + γ·|S|`, the O(n²) DP recurrence, and the
  joint logR + mirrored-BAF cost with common breakpoints [6][7].
- Sub-clonal copy number: the Battenberg two-population model — one or two integer states whose cellular fractions
  sum to 1 — with the two states being the bracketing integers and a shared fraction f [8].

**Intentionally simplified:**

- Purity/ploidy fit: fixed-grid minimum with a lower-ploidy tie-break, **not** the ASCAT iterative refit;
  **consequence:** the (ρ, ψ) optimum is grid-resolution limited.
- Sub-clonal fit: two adjacent integer states with one shared fraction, **not** an arbitrary multi-state mixture;
  **consequence:** three-or-more cell populations per segment are not modelled.

**Not implemented:**

- Multi-sample (asmultipcf) segmentation and whole-genome-doubling refit search; **users should rely on:**
  ASCAT/Battenberg/FACETS for those; this unit covers the single-sample ASPCF + two-state sub-clonal derivation
  feeding the downstream `EstimatePurity`/`EstimatePloidy`/`EstimateCcf`/`ClassifyClonality`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Greedy mean-shift segmentation retained alongside ASPCF | Deviation | breakpoint sensitivity of the greedy path | accepted | ASPCF (`SegmentAlleleSpecificAspcf`) is the global-optimum path [6] |
| 2 | Two-state sub-clonal mixture (adjacent integers) | Assumption | no 3+-population segments | accepted | Battenberg single-fraction segment [8] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Balanced-only genome (all b=0.5) | fit completes; segments ×0.05 weighted | source [2] |
| Single locus per chromosome | one segment per chromosome, LocusCount=1 | segmentation contract |
| VAF rounds to 0 | multiplicity clamped to 1 | INV-02 |
| VAF rounds above major CN | multiplicity clamped to major CN | INV-02 |
| ASPCF flat track (no change) | single segment | INV-06 |
| ASPCF γ → ∞ | single segment per chromosome | INV-06 [6] |
| Sub-clonal: integer (nA,nB) | single clonal state, f=1 | INV-07 [8] |

### 6.2 Limitations

logR and BAF are observed measurements and are always a caller input — this is inherent, not a limitation of the
derivation. The derivation models sub-clonal copy number only as a two-population (two adjacent integer states)
mixture per segment; it does not model 3+ populations per segment, multi-sample (asmultipcf) segmentation, or a
whole-genome-doubling refit search, and assumes germline-het BAF loci.

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

- Tests: [OncologyAnalyzer_AscatDerivation_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatDerivation_Tests.cs) — covers `INV-01`–`INV-07`
- Evidence: [ONCO-ASCAT-001-Evidence.md](../../../docs/Evidence/ONCO-ASCAT-001-Evidence.md)
- Related algorithms: [Tumor_Ploidy_Estimation](./Tumor_Ploidy_Estimation.md), [Cancer_Cell_Fraction_Estimation](./Cancer_Cell_Fraction_Estimation.md), [Tumor_Purity_Estimation](./Tumor_Purity_Estimation.md)

## 8. References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R`. https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ. 2021. DeCiFering the elusive cancer cell fraction. PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
6. Nilsen G, Liestøl K, Van Loo P, et al. 2012. Copynumber: Efficient algorithms for single- and multi-track copy number segmentation. BMC Genomics 13:591. https://doi.org/10.1186/1471-2164-13-591
7. Ross EM, Haase K, Van Loo P, Markowetz F. 2021. Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37(13):1909–1911. https://doi.org/10.1093/bioinformatics/btaa538
8. Nik-Zainal S, Van Loo P, Wedge DC, et al. 2012. The Life History of 21 Breast Cancers. Cell 149(5):994–1007. https://doi.org/10.1016/j.cell.2012.04.023 ; Battenberg, https://github.com/Wedge-lab/battenberg
