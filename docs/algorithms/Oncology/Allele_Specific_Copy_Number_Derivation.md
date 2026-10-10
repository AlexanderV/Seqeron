# Allele-Specific Copy-Number Derivation (ASCAT-style purity/ploidy fit + multiplicity)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ASCAT-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Complete (ASCAT runASCAT / ascat.aspcf, Battenberg determine_copynumber ports; `SegmentAlleleSpecific` = ASPCF at penalty 70 since B24 F35) |
| Last Reviewed | 2026-09-28 |

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

For massively parallel sequencing data **γ = 1** [2]. ASCAT (`runASCAT`, ascat.runAscat.R [2]) computes a
distance matrix over ψ ∈ seq(min_ploidy − 0.5, max_ploidy + 0.5, 0.05) × ρ ∈ seq(min_purity, max_purity, 0.01)
(defaults 1.5 / 5.5 / 0.1 / 1.05) from the **autosomal** segments, `length` being the number of heterozygous probes
(`make_segments`) and the minor allele chosen genome-wide (`sum(nA) < sum(nB)`):

```
d(ρ,ψ) = Σ_segments  (n_minor − max(round(n_minor), 0))² · length · w_b      w_b = 0.05 if b = 0.5 else 1
TheoretMaxdist = Σ_segments 0.25 · length · w_b
goodnessOfFit  = (1 − d / TheoretMaxdist) · 100   [%]
```

(0.25 = (½)² is the worst-case squared distance to an integer [2].) Candidate solutions are the cells that are the
strict minimum of their 7 × 7 neighbourhood; a four-pass filter cascade keeps candidates with recomputed ploidy
Σ(nA+nB)·length/Σlength in (min_ploidy, max_ploidy), ρ ≥ 0.2, GoF > 80 % and percentzero > 0.02 (pass 1), then
relaxes to 1.7 < ploidy < 2.3 with perczeroAbb > 0.1 (pass 2), ρ > 1 columns masked with percentzero / perczeroAbb /
percOddEven alternatives (pass 3) and strict ploidy alone (pass 4). The smallest-distance candidate wins (ρ > 1
reported as 1); with no candidate ASCAT returns rho = NA. Integer segments (`seg_raw`) fold a negative allele into
the other, round half-to-even and, for BAF = 0.5, apply the odd-total rule (`limitround = 0.5`).

Mutation multiplicity (number of mutated copies per cancer cell) is the rounded observed mutation copy number
[3][4]:

```
n_mut = VAF · (1/ρ) · [ρ·N_T + 2(1−ρ)]
m     = clamp( round(n_mut), 1, majorCopyNumber )     round = R round(), half-to-even (2.5 → 2, 1.5 → 2)
```

The rounding is that of facets-suite `expected_mutant_copies` (mskcc/facets-suite `R/ccf-annotate-maf.R`, "Based on
PMID 28270531"): `mu < 1 → 1`, then R `round` (IEC 60559 half-to-even). The floor at 1 is facets-suite's; the cap at
the major-allele copy number is a Seqeron extra (facets-suite leaves `expected_alt_copies` uncapped). Exact `.5`
ties: `(VAF, ρ, N_T) = (0.625, 1, 4)` → n_mut = 2.5 → m = 2; `(0.375, 0.5, 2)` → 1.5 → 2; `(0.5625, 1, 8)` → 4.5 → 4
(R 4.x `expected_mutant_copies` output; FIN-B24 F26).

which is the inversion of the PICTograph generative model VAF = m·CCF·ρ / (N_T·ρ + 2(1−ρ)) at clonal CCF = 1 [4].

**ASPCF segmentation (penalised least squares).** ASCAT segments via Piecewise Constant Fitting [6][7] (the
implementation is a port of `ascat.aspcf` / `fastAspcf` / `aspcfpart` [2], see §4.1 step 4): minimise

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

where in ASCAT each track's SSE is divided by its MAD-based variance (`getMad`: MAD of the residuals from a running
median), so the per-segment data cost is (logR-SSE/sd₁² + BAF-SSE/sd₂²), γ is charged per breakpoint and every
segment has at least kmin = 6 loci [2]. BAF is mirrored to a single allelic-imbalance track before segmentation [7];
a segment's BAF is 0.5 + mean|b − 0.5|, shrunk to 0.5 when sqrt(sd₂² + μ²) < 2·sd₂ [2].

**Sub-clonal copy number (Battenberg two-population model).** A segment has either one integer state (clonal, all
tumour cells) or two integer states (sub-clonal, two cell populations whose fractions sum to 1) [8]. Battenberg's
`determine_copynumber` (R/fitcopynumber.R, R/orderEdges.R) takes the corners of the **nearest edge** of the
copy-number square around (nMajor, nMinor) — two states differing in one allele — chosen from the segment BAF l
relative to the corner BAF levels (1 − ρ + ρ·M)/(2 − 2ρ + ρ·(M + m)) and the priority `ntot < x + y + 1`, and solves the
BAF mixture for the fraction of state 1:

```
τ = (1 − ρ + ρ·M₂ − 2l(1 − ρ) − lρ(m₂ + M₂)) / (lρ(m₁ + M₁) − lρ(m₂ + M₂) − ρM₁ + ρM₂)
```

A segment is clonal when |l − BAF of the closest corner| < maxdist = 0.01 (or the per-SNP t-test is not significant;
a segment summary has constant SNP BAF, for which Battenberg sets pval = 0) [8].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A single (ρ, ψ) explains the genome | `FitPurityPloidy` rounds to integers; sub-clonal segments are modelled separately by `FitSubclonalCopyNumber` |
| ASM-02 | γ matches the platform (γ = 1 for sequencing) | Wrong γ rescales logR → biased copy number [2] |
| ASM-03 | BAF is measured at germline-heterozygous SNPs | BAF at homozygous sites is uninformative; folding assumes het loci. With germline genotypes (`SegmentAlleleSpecificAspcf(loci, germlineHeterozygous)`, B24 F36) homozygous loci contribute logR only, as in `ascat.aspcf` |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | At the true (ρ₀, ψ₀) of an integer-CN genome, d ≈ 0 and GoF ≈ 100% (`EvaluatePurityPloidy`) | nA, nB equal exact integers by construction [1][2] |
| INV-02 | Derived multiplicity ∈ [1, majorCopyNumber] | explicit clamp; a variant sits on ≥ 1 and ≤ major-allele copies [3][4] |
| INV-03 | GoF ≤ 100% | d ≥ 0 and d ≤ TheoretMaxdist by the 0.25 worst-case bound [2] |
| INV-04 | major ≥ minor ≥ 0 in every emitted segment | ASCAT BAF ≤ 0.5 orientation ⇒ nA ≥ nB; negative-value correction [2] |
| INV-05 | ASPCF output equals `ascat.aspcf` (per window the standardised cost is minimised exactly over segmentations with ≥ kmin loci) | aspcfpart DP [2][6] |
| INV-06 | A segment with no logR/BAF change is one ASPCF segment; γ→∞ ⇒ \|S\|=1 | penalty dominates SSE gains [6] |
| INV-07 | Sub-clonal state fractions sum to 1 (τ unclamped, as Battenberg); a clonal segment has f=1 and no second state | two-population model [8] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| loci | IEnumerable\<AlleleSpecificLocus\> | required | per-locus (chrom, pos, logR, BAF) measurements | non-null; chrom non-null; finite logR; BAF ∈ [0,1] |
| logRChangeThreshold, bafChangeThreshold, minLociPerSegment | double, double, int | required, 0.1, 1 | `SegmentAlleleSpecific` legacy parameters — **ignored** since B24 F35 (the method runs ASPCF at γ = 70; no ASPCF equivalent), still range-checked for compatibility | > 0, > 0, ≥ 1 |
| germlineHeterozygous | IReadOnlyList\<bool\> | — (overload) | germline genotype per locus for the germline-aware ASPCF (true = heterozygous; ASCAT `germlinegenotypes == FALSE`) (B24 F36) | non-null; one per locus; BAF ∈ [0,1] required at heterozygous loci only (homozygous BAF ignored, may be NaN) |
| penalty (ASPCF γ) | double | 70.0 | per-breakpoint penalty on the standardised cost (`ascat.aspcf`) | > 0, finite |
| purity, ploidy (sub-clonal) | double | required | fitted ρ, ψ for the sub-clonal decomposition | ρ∈(0,1]; ψ>0 |
| segments | IReadOnlyList\<AlleleSpecificSegmentSummary\> | required | segment summaries for the fit | non-empty; finite logR; BAF ∈ [0,1]; LocusCount ≥ 1; ≥ 1 autosomal |
| purityMin/Max/Step | double | 0.1 / 1.05 / 0.01 | purity grid (ASCAT min/max_purity) | min ∈ (0,1]; max ≥ min, finite; step > 0 |
| ploidyMin/Max/Step | double | 1.5 / 5.5 / 0.05 | ASCAT min/max_ploidy filter; ψ grid spans ±0.5 beyond | > 0; max ≥ min; step > 0 |
| gamma | double | 1.0 | platform γ | > 0 |
| vaf, purity, totalCopyNumber, majorCopyNumber | double/int | required | multiplicity inputs | vaf∈[0,1]; ρ∈(0,1]; N_T≥1; major∈[1,N_T] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| AlleleSpecificSegmentSummary | record | per-segment mean logR, mirrored mean BAF, locus count |
| AspcfSegmentation | class | germline-aware ASPCF (F36): per-locus `SegmentedLogR` (ASCAT `Tumor_LogR_segmented`, every locus), `SegmentedBaf` (mirrored, NaN at homozygous loci) and `Segments` (`AspcfSegment`: runASCAT logR runs within a chromosome, BAF of the first heterozygous locus or NaN, `LocusCount`, `HeterozygousLocusCount`) |
| PurityPloidyFit | record | ρ, ASCAT output ploidy (probe-weighted mean integer total CN), GoF %, integer `AlleleSpecificSegment`s (one per summary, ASCAT `seg_raw`), `Psi` (ψ), `IsNonAberrant` |
| DeriveMultiplicity | int | integer multiplicity m ∈ [1, majorCopyNumber] |
| SubclonalSegmentFit | record | per-segment primary/secondary `SubclonalCopyNumberState` (major, minor, cellFraction) + IsSubclonal flag |

### 3.3 Preconditions and Validation

Positions are 0-based. Null `loci`/`segments` → `ArgumentNullException`; empty/malformed `segments` (or no autosomal
segment) → `ArgumentException`; out-of-range thresholds, grid bounds, or multiplicity arguments →
`ArgumentOutOfRangeException`; no acceptable ASCAT optimum → `InvalidOperationException` from `FitPurityPloidy`
(`TryFitPurityPloidy` returns false). Both segmentation entry points (`SegmentAlleleSpecific` delegates to
`SegmentAlleleSpecificAspcf`) reject a locus with a non-finite logR or a BAF outside [0, 1] (`ArgumentException`). BAF is
mirrored about 0.5 (ascat.aspcf `ifelse(b > 0.5, b, 1 − b)`) during segmentation so the two symmetric het clusters reinforce.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Segment (ASPCF, `ascat.aspcf`):** per contiguous chromosome run, minimise the MAD-standardised joint logR +
   mirrored-BAF squared error plus γ per breakpoint (kmin 6; γ = 70 for `SegmentAlleleSpecific`, B24 F35). Summarise
   each segment by mean raw logR and the ASPCF mirrored BAF. Segmenting on BAF as well as logR is essential because copy-neutral LOH (e.g. 2:0) shares a balanced region's
   logR but not its BAF.
2. **Fit (`runASCAT`):** build the distance matrix over the autosomal segments, collect strict 7 × 7 local minima
   through the four-pass filter cascade, keep the smallest distance (ρ > 1 ⇒ 1), emit the `seg_raw` integer
   segments for every summary (sex chromosomes with the diploid model), the ASCAT ploidy and GoF. No candidate ⇒
   rho = NA. `EvaluatePurityPloidy` is the rho_manual/psi_manual path.
3. **Multiplicity:** m = clamp(round_half_even(VAF·[ρ·N_T + 2(1−ρ)]/ρ), 1, major). Feed (VAF, ρ, N_T, m) into `EstimateCcf`.
4. **ASPCF (`ascat.aspcf`, alternative to step 1):** per chromosome, MAD-winsorise logR and mirrored BAF; < 6 loci ⇒
   one segment; else `fastAspcf`: 1000-locus windows (overlap 100), per window MAD sd of both tracks (a window with
   sd 0 is skipped), `aspcfpart` exact DP with kmin 6 on the standardised joint cost; segment logR = mean raw logR,
   BAF = 0.5 + μ (shrunk to 0.5 when sqrt(sd₂² + μ²) < 2·sd₂); repeat with the next penalty of the ladder while ≥ 800
   distinct logR levels remain [2][6][7].
   **With germline genotypes (B24 F36, `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous, penalty)`)** the full
   `ascat.aspcf` runs: homozygous stretches = runs of ≥ ceiling(log(0.001, h)) homozygous loci (h = homozygous fraction,
   `predictGermlineHomozygousStretches`); logR of all loci is winsorised and averaged onto each heterozygous locus
   (midpoint windows); ASPCF runs on the heterozygous loci; logR breakpoints between heterozygous loci of different
   levels are placed in the homozygous gap at the minimum absolute deviation ("find best breakpoint", R's `1:0`
   quirk kept); a chromosome without heterozygous loci is one logR segment without BAF; each homozygous stretch is
   re-segmented by `exactPcf` (kmin 6, penalty floor(γ/4)) over ± 100 loci, replacing ± 5 loci that differ by > 0.3 when
   > 5 differ; genome-wide `fillNA(zeroIsNA = TRUE)` and level re-estimation (mean raw logR per run).
5. **Sub-clonal fit (Battenberg `determine_copynumber`):** nMajor/nMinor at (ρ, ψ) from l = max(b, 1 − b); nearest
   edge (`orderEdges` option 1); clonal if the closest corner is within 0.01 BAF, else two states with fraction τ [8].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- γ = 1 (sequencing) [2]; balanced (BAF = 0.5, exact equality as in R) segments weighted ×0.05 [2]; worst-case
  integer distance 0.25 [2]; ASCAT constants MINRHO 0.2, MINGOODNESSOFFIT 80, MINPERCZERO 0.02, MINPERCZEROABB 0.1,
  MINPERCODDEVEN 0.05, MINPLOIDYSTRICT 1.7, MAXPLOIDYSTRICT 2.3, MINABB 0.03, MINABBREGION 0.005 [2].
- Battenberg constants maxdist 0.01, cn_upper_limit 1000, minimum minor CN 0.01 [8].
- Segments with End == Start are emitted with a 1 bp span so `AlleleSpecificSegment.Length > 0`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ASPCF segmentation | O(L·W) per chromosome (W = 1000-locus window) | O(W) | windowed PCF DP [2][6] |
| Purity/ploidy fit | O(P·Q·S) | O(P·Q + S) | P,Q = grid sizes (≤ 4·10⁶ cells), S = segments |
| Multiplicity | O(1) | O(1) | closed form |
| Sub-clonal fit | O(S) | O(S) | closed-form decomposition per segment |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.CopyNumberPloidy.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.CopyNumberPloidy.cs)

- `OncologyAnalyzer.SegmentAlleleSpecific(...)`: `SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty = 70)` (B24 F35). It formerly ran an unsourced greedy mean-shift heuristic; that heuristic was removed and its threshold parameters are now ignored (kept for source compatibility, still range-checked).
- `OncologyAnalyzer.SegmentAlleleSpecificAspcf(...)`: ASCAT ASPCF (`ascat.aspcf` port) on heterozygous loci.
- `OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, germlineHeterozygous, penalty)`: complete `ascat.aspcf` with germline genotypes (homozygous probes, homozygous-stretch resegmentation) → `AspcfSegmentation` (B24 F36).
- `OncologyAnalyzer.FitPurityPloidy(...)` / `TryFitPurityPloidy(...)`: ASCAT `runASCAT` fit → ρ, ψ, ploidy, GoF, integer segments.
- `OncologyAnalyzer.EvaluatePurityPloidy(...)`: ASCAT rho_manual/psi_manual path.
- `OncologyAnalyzer.DeriveMultiplicity(...)`: McGranahan multiplicity (rounded, clamped).
- `OncologyAnalyzer.FitSubclonalCopyNumber(...)`: Battenberg `determine_copynumber` (nearest edge, τ, maxdist).

### 5.2 Current Behavior

Single-sample ASCAT fit (runASCAT port, R-verified 150/150); ASPCF = `ascat.aspcf` port (R-verified 60/60);
sub-clonal = Battenberg `determine_copynumber` port (R-verified 402/402). `SegmentAlleleSpecific` runs the same ASPCF
at ASCAT's default penalty 70 (the former greedy heuristic was removed, B24 F35). Not a search/matching task, so the repository
suffix tree is **not used** (no occurrence enumeration).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ASCAT `runASCAT`: nA/nB equations, distance matrix (probe-count weights, autosomes, genome-wide minor allele),
  local-minimum scan, filter cascade, ρ clamp, `seg_raw` rounding, GoF, non-aberrant flag, manual (ρ, ψ) path [2].
- McGranahan observed mutation copy number n_mut and the [1, major] multiplicity clamp [3][4].
- ASCAT ASPCF: winsorisation, MAD-standardised joint cost, kmin 6, windows, BAF shrinkage, penalty ladder [2][6][7].
- ASCAT `ascat.aspcf` with germline genotypes (B24 F36): heterozygous-only BAF, logR from all probes averaged onto the
  heterozygous probes, gap breakpoint search, `predictGermlineHomozygousStretches` + `exactPcf` resegmentation,
  `fillNA`, genome-wide level re-estimation — R-verified on 3 genomes (26 segments, ≤ 5.6e-16 per locus) [2].
- Battenberg `determine_copynumber`: nearest edge (`orderEdges` option 1), τ, maxdist clonality, negative-minor
  adjustment [8].

**Intentionally simplified:**

- ~~ASCAT inputs are segment summaries: germline-homozygous probes (their logR averaging and the homozygous-stretch
  resegmentation of `ascat.aspcf`) … are not available~~ — **resolved by F36** (germline-aware
  `SegmentAlleleSpecificAspcf` overload). The gender-specific haploid X/Y model is not available; sex-chromosome
  segments are excluded from the fit and emitted with the diploid (gender "XX") model.
- Sub-clonal fit: a summary carries no per-SNP BAF spread, so Battenberg's t-test cannot be run (pval = 0 as for a
  constant-BAF segment; only maxdist decides clonality); bootstrap CIs and alternative solutions B–F are not produced.

**Not implemented:**

- Multi-sample (asmultipcf) segmentation and whole-genome-doubling refit search; **users should rely on:**
  ASCAT/Battenberg/FACETS for those; this unit covers the single-sample ASPCF + two-state sub-clonal derivation
  feeding the downstream `EstimatePurity`/`EstimatePloidy`/`EstimateCcf`/`ClassifyClonality`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ~~Greedy mean-shift segmentation retained alongside ASPCF~~ | Deviation | — | **resolved (B24 F35)** | `SegmentAlleleSpecific` now delegates to ASPCF (`ascat.aspcf`, penalty 70) [2]; legacy thresholds ignored |
| 2 | Segment summaries instead of probes | Assumption | ~~no homozygous-probe logR~~ (resolved F36), no haploid X/Y model | partly resolved | see §5.3 |
| 3 | No per-SNP t-test in the sub-clonal fit | Assumption | clonality by maxdist only | accepted | constant-BAF branch of Battenberg [8] |
| 4 | `PurityPloidyFit.Ploidy` = probe-weighted mean integer CN over heterozygous probes | Assumption | ASCAT averages over all probes | accepted | `Psi` carries ψ |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Balanced-only genome (all b=0.5) | non-aberrant flag; fit may fail (rho = NA) | source [2] |
| Single-segment genome | no strict local minimum passes the filters ⇒ rho = NA | source [2] (R run) |
| Noise-free logR/BAF | ASPCF places no breakpoint (MAD sd = 0) | source [2] (R run) |
| Single locus per chromosome | one segment per chromosome, LocusCount=1 | segmentation contract |
| VAF rounds to 0 | multiplicity clamped to 1 | INV-02 |
| VAF rounds above major CN | multiplicity clamped to major CN | INV-02 |
| ASPCF flat track (no change) | single segment | INV-06 |
| ASPCF γ → ∞ | single segment per chromosome | INV-06 [6] |
| ASPCF chromosome with < 6 loci | one segment, mean winsorised mirrored BAF | source [2] |
| Germline-aware ASPCF: chromosome without heterozygous loci | one logR segment (mean raw logR), no BAF (`AspcfSegment.HasBaf` false) | source [2] (R run, F36) |
| Germline-aware ASPCF: copy-number change inside a homozygous stretch | found by the `exactPcf` resegmentation when > 5 loci differ by > 0.3 | source [2] (R run, F36) |
| Sub-clonal: BAF within 0.01 of a corner | single clonal state, f=1 | INV-07 [8] |

### 6.2 Limitations

logR and BAF are observed measurements and are always a caller input — this is inherent, not a limitation of the
derivation. The unit works on heterozygous-locus segment summaries (germline-homozygous probes: use the germline-aware
ASPCF overload, F36); the haploid X/Y
(male) model, multi-sample (asmultipcf) segmentation and Battenberg's haplotype phasing / per-SNP t-test / bootstrap
are out of scope (phased BAFs need an imputation reference panel). `FitPurityPloidy` fails (like ASCAT) when no local
minimum passes the filters — use `TryFitPurityPloidy`, or `EvaluatePurityPloidy` with externally chosen (ρ, ψ).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var loci = /* per-locus (chrom, pos, logR, BAF) measurements */;
var summaries2 = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci); // ASCAT ASPCF segmentation (penalty 70)
// SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2) returns the same (legacy name; thresholds ignored, F35)
var fit = OncologyAnalyzer.FitPurityPloidy(summaries2);         // → ρ, ψ, integer segments (throws if ASCAT finds none)
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
- Tests: [OncologyAnalyzer_AscatGermlineHomozygous_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatGermlineHomozygous_Tests.cs) — germline-homozygous probes (F36/F37), R-locked
- Evidence: [ONCO-ASCAT-001-Evidence.md](../../../docs/Evidence/ONCO-ASCAT-001-Evidence.md)
- Related algorithms: [Tumor_Ploidy_Estimation](./Tumor_Ploidy_Estimation.md), [Cancer_Cell_Fraction_Estimation](./Cancer_Cell_Fraction_Estimation.md), [Tumor_Purity_Estimation](./Tumor_Purity_Estimation.md)

## 8. References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R`, `ASCAT/R/ascat.aspcf.R` (master, read 2026-09-28). https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ. 2021. DeCiFering the elusive cancer cell fraction. PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
6. Nilsen G, Liestøl K, Van Loo P, et al. 2012. Copynumber: Efficient algorithms for single- and multi-track copy number segmentation. BMC Genomics 13:591. https://doi.org/10.1186/1471-2164-13-591
7. Ross EM, Haase K, Van Loo P, Markowetz F. 2021. Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37(13):1909–1911. https://doi.org/10.1093/bioinformatics/btaa538
8. Nik-Zainal S, Van Loo P, Wedge DC, et al. 2012. The Life History of 21 Breast Cancers. Cell 149(5):994–1007. https://doi.org/10.1016/j.cell.2012.04.023 ; Battenberg `R/fitcopynumber.R` (`determine_copynumber`), `R/orderEdges.R` (master, read 2026-09-28), https://github.com/Wedge-lab/battenberg
