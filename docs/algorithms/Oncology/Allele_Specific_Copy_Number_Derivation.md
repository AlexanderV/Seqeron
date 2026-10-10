# Allele-Specific Copy-Number Derivation (ASCAT-style purity/ploidy fit + multiplicity)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-ASCAT-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Complete (ASCAT runASCAT / ascat.aspcf / ascat.asmultipcf, Battenberg determine_copynumber / segment.baf.phased / merge_segments / mask_high_cn_segments / callSubclones ports, R-locked; `SegmentAlleleSpecific` = ASPCF at penalty 70 since B24 F35; only Battenberg's external haplotype imputation is not built in) |
| Last Reviewed | 2026-10-10 |

## 1. Overview

This unit derives the upstream quantities that the tumour purity/ploidy/CCF/clonality methods consume, so the
pipeline can run from per-locus signal rather than caller-supplied pre-made allele-specific segments and
multiplicity. From per-locus log-R ratio (logR) and B-allele frequency (BAF) at germline-heterozygous SNPs
(observed measurements), it (1) **segments** the genome into (logR, BAF) summaries, (2) **jointly fits** tumour
purity ρ and ploidy ψ by grid search using the ASCAT equations and goodness-of-fit objective, emitting
allele-specific **integer** copy-number segments, and (3) **derives** somatic mutation multiplicity from VAF,
purity and copy number. It is a port of ASCAT (Van Loo et al. 2010) [1][2] — `runASCAT` (incl. the germline-aware,
homozygous-segment and haploid X/Y male paths), `ascat.aspcf` and multi-sample `ascat.asmultipcf`, R-locked (§5.2) — plus
the McGranahan multiplicity convention [3][4]. In addition it provides the **ASPCF** penalised
least-squares segmentation (Nilsen et al. 2012 [6]; Ross et al. 2021 [7]) — the global-optimum joint logR/BAF
changepoint method ASCAT uses — and **sub-clonal copy number** modelling (Battenberg two-population model,
Nik-Zainal et al. 2012 [8]), which expresses a segment that does not fit a single integer state as a mixture of
two adjacent integer states with a sub-clonal cellular fraction; the Battenberg chain (per-SNP t-test, phased-BAF
segmentation, bootstrap CIs, `merge_segments`, `mask_high_cn_segments`, `callSubclones`) is ported and R-locked. The
only Battenberg stage not built in is haplotype imputation (IMPUTE2/Beagle5 against the 1000 Genomes panel; callers
supply phased BAFs, §5.3).

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

A segment is clonal when |l − BAF of the closest corner| < maxdist = 0.01 (pval = 1) or when the two-sided one-sample
t-test of its phased SNP BAFs against that corner's BAF is not significant (pval > siglevel = 0.05); fewer than two SNPs
or constant SNP BAFs give pval = 0 [8]. `FitSubclonalCopyNumber` (segment summaries, constant BAF) therefore decides by
maxdist alone; `FitSubclonalCopyNumberWithSnpTest` (B24 F39) runs the t-test on the supplied SNP BAFs.

**Alternative solutions and confidence intervals (Battenberg `determine_copynumber`, B24 F41).** For a sub-clonal
segment Battenberg reports all six `orderEdges` solutions A–F (A = nearest edge; rows needing a negative copy number
are NA and moved last), each with τ, the delta-method `SDfrac = |τ(l + s) − τ|/2 + |τ(l − s) − τ|/2` (s = sd(BAFke)/√n)
and a non-parametric bootstrap: noperms = 1000 resamples of the segment's SNP BAFs with replacement, τ at each resample
mean; `SDfrac_BS = sd(τ*)`, `frac1_0.025 = sort(τ*)[25]`, `frac1_0.975 = sort(τ*)[975]` (fixed order statistics — the
2.5 %/97.5 % points only for noperms = 1000; NA below). The RNG is R's (`set.seed` → Mersenne-Twister, rejection
`sample`), so with the same seed the resamples are Battenberg's [8].

**Merging and masking (Battenberg `merge_segments` / `mask_high_cn_segments`, B24 F60).** After the first fit,
`callSubclones` merges adjacent segments "if there is not enough evidence for them to be separate": per chromosome the
smallest segment (width end − start + 1) with an unchecked neighbour is compared with its neighbours, closest first —
never across > 3 Mb (`GenomicRanges::distance`); merged when both are clonal with the same (nMaj1_A, nMin1_A); otherwise
merged when they lie in the same square (round(nmin) or round(nmaj) equal, with
nmin = (ρ − 1 − (BAF − 1)·2^(LogR/γ)·(2(1 − ρ) + ρψ))/ρ, nmaj = (ρ − 1 + BAF·2^(LogR/γ)·(2(1 − ρ) + ρψ))/ρ and ψ the
psi of all cells — Battenberg's formula verbatim) and neither the Welch two-sample t-test on their logR nor on their
`BAFphased` is significant (p < 0.05; > 10 values each). The surviving neighbour is extended, its BAF recomputed
(`calc_seg_baf_option`) and written into `BAFseg`, its LogR re-averaged, and both sides re-checked. After the second fit,
segments with nMaj1_A or nMin1_A > `max_allowed_state` (250) are masked (A's copy numbers NA; `BAFseg` NA for
startpos < Position ≤ endpos) [8]. `CallBattenbergSubclones` (B24 F61) runs `callSubclones`' copy-number sequence
`set.seed(seed)` → `determine_copynumber` → `merge_segments` → `determine_copynumber` → `mask_high_cn_segments` with
**one** RNG stream, so the second fit's bootstrap resamples continue where the first fit's stopped, exactly as in R.

**Phased-BAF segmentation (Battenberg `segment.baf.phased`, B24 F40).** Battenberg phases germline-heterozygous SNPs
with IMPUTE2/Beagle5 against the 1000 Genomes reference panel (external executables + multi-GB reference bundle; not
run here — the caller supplies phased BAFs), then per chromosome pre-segments at prior (SV) breakpoints and SNP gaps
≥ 3 Mb and, per presegment, estimates sd = max(0.09, getMad(min(BAF, 1 − BAF), 25)); a first fast PCF
(`selectFastPcf`, kmin = phasekmin = 3, penalty phasegamma·sd = 3·sd) on the raw phased BAF locates switched haplotype
blocks, `BAFphased = BAF if BAFsegm > 0.5 else 1 − BAF`; a second PCF (kmin 3, penalty gamma·sd = 10·sd) on BAFphased
gives the segments, whose BAF is the median of BAFphased (`calc_seg_baf_option` 3: the mean when the median is 0 or 1;
1: median; 2: PCF mean). Presegments of < 50 SNPs are not segmented (mean). `selectFastPcf` is the copynumber-package
fast PCF (`filterMarkS4` candidate breakpoints from 6L/6L2 high-pass and kmin filters with type-7 quantile limits,
then exact Potts filtering on the compacted array; windowed `runPcfSubset` above 15 000 SNPs) [8].

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
| germlineHeterozygous | IReadOnlyList\<bool\> | — (overload) | germline genotype per locus for the germline-aware ASPCF (true = heterozygous; ASCAT `germlinegenotypes == FALSE`) (B24 F36) | non-null; one per locus; BAF ∈ [0,1] or NaN (= R NA, B24 F59) at heterozygous loci (homozygous BAF ignored); logR finite or NaN (NA) |
| penalty (ASPCF γ) | double | 70.0 | per-breakpoint penalty on the standardised cost (`ascat.aspcf`); for `SegmentAlleleSpecificAsMultiPcf` on the raw joint cost (`ascat.asmultipcf`, ladder 25/50/100/200/400/800) | > 0, finite |
| samples (multi-sample) | IReadOnlyList\<IReadOnlyList\<AlleleSpecificLocus\>\> | — (`SegmentAlleleSpecificAsMultiPcf`) | ≥ 2 samples on one probe set (same chromosome/position per index), one germline (`germlineHeterozygous`, null = all heterozygous) (B24 F53) | logR finite or NaN (= R NA); heterozygous BAF ∈ [0,1] or NaN (B24 F59); ±∞ rejected; every chromosome part ≥ 2 probes (R fails on 1 sample / 1-probe part) |
| options (multi-sample) | AsMultiPcfOptions | ASCAT defaults | `SampleWeights` (`wsample`, S or 2·S positive values, null = 1), `Algorithm` (`Exact` / `Fast` = `selectAlg`), `Refine` (true) | weights > 0, finite |
| purity, ploidy (sub-clonal) | double | required | fitted ρ, ψ for the sub-clonal decomposition | ρ∈(0,1]; ψ>0 |
| segments (t-test) | IReadOnlyList\<SubclonalSegmentSnpBafs\> | — (`FitSubclonalCopyNumberWithSnpTest`) | segment summary + its phased SNP BAFs (Battenberg `BAFphased`, used unmirrored) (B24 F39) | non-null; SNP BAFs ∈ [0,1] (list may be empty) |
| significanceLevel, maxBafDistance | double | 0.05, 0.01 | Battenberg `siglevel`, `maxdist` | siglevel ∈ [0,1]; maxdist ≥ 0, finite |
| seed, permutations | int, int | required, 1000 | `FitSubclonalCopyNumberWithBootstrap`: R `set.seed` value and Battenberg `noperms` (B24 F41) | any int; permutations ≥ 1 |
| calls, segmentedSnps, logR, bafOption (merge) | IReadOnlyList\<BattenbergSegmentCall\>, rows, IReadOnlyList\<LogRProbe\>, BattenbergSegmentBafOption | — (`MergeBattenbergSegments`), option 3 | fitted profile + its `BAFsegmented` rows + raw logR (B24 F60) | BAFs ∈ [0,1]; logR finite or NaN (= NA); every chromosome of the rows needs ≥ 1 call and ≥ 1 logR probe (R `stopifnot`) |
| maxAllowedState | int | 250 | `MaskHighCopyNumberSegments` / `CallBattenbergSubclones`: Battenberg `max_allowed_state` (B24 F60/F61) | any int |
| segmentedSnps, logR, seed (driver) | rows, IReadOnlyList\<LogRProbe\>, int | — (`CallBattenbergSubclones`) | `BAFsegmented` + raw logR + `set.seed` value; also γ, siglevel, maxdist, permutations, maxAllowedState, bafOption (B24 F61) | as `BuildBattenbergSegments` / `MergeBattenbergSegments`; permutations ≥ 1 |
| snps (phased segmentation) | IReadOnlyList\<PhasedBafSnp\> | — (`SegmentPhasedBaf`) | caller-phased SNP BAFs (chromosome, position, BAF of haplotype 1) (B24 F40) | non-null; BAF ∈ [0,1] or NaN (= missing, dropped) |
| options (phased segmentation) | BattenbergPhasedSegmentationOptions | Battenberg defaults | gamma 10, phasegamma 3, kmin 3, phasekmin 3, no_segmentation false, calc_seg_baf_option 3, prior breakpoints none | gammas finite ≥ 0; kmins ∈ [1, 14] |
| segmentedSnps, logR | IReadOnlyList\<PhasedBafSegmentedSnp\>, IReadOnlyList\<LogRProbe\> | — (`BuildBattenbergSegments`) | `segment.baf.phased` rows + raw logR probes → `SubclonalSegmentSnpBafs` for the t-test fit (B24 F40) | BAFs ∈ [0,1]; non-finite logR ignored |
| segments | IReadOnlyList\<AlleleSpecificSegmentSummary\> | required | segment summaries for the fit | non-empty; finite logR; BAF ∈ [0,1]; LocusCount ≥ 1; ≥ 1 autosomal |
| purityMin/Max/Step | double | 0.1 / 1.05 / 0.01 | purity grid (ASCAT min/max_purity) | min ∈ (0,1]; max ≥ min, finite; step > 0 |
| ploidyMin/Max/Step | double | 1.5 / 5.5 / 0.05 | ASCAT min/max_ploidy filter; ψ grid spans ±0.5 beyond | > 0; max ≥ min; step > 0 |
| gamma | double | 1.0 | platform γ | > 0 |
| sexModel | AscatSexModel | — (overload; gender-less = `Female`) | runASCAT `gender` / `X_nonPAR` for the emitted X/Y segments (B24 F38): `Female` (XX, ASCAT default), `Male` (XY, whole X haploid = `genomeVersion = NULL`), `MaleWithXNonPar(GRCh37 \| GRCh38)` (ASCAT hg19/hg38 constants; `XNonParChm13` also exposed) or `new(gender, (start, end))` | non-null; Start ≤ End |
| vaf, purity, totalCopyNumber, majorCopyNumber | double/int | required | multiplicity inputs | vaf∈[0,1]; ρ∈(0,1]; N_T≥1; major∈[1,N_T] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| AlleleSpecificSegmentSummary | record | per-segment mean logR, mirrored mean BAF, locus count |
| AspcfSegmentation | class | germline-aware ASPCF (F36): per-locus `SegmentedLogR` (ASCAT `Tumor_LogR_segmented`, every locus), `SegmentedBaf` (mirrored, NaN at homozygous loci and at heterozygous loci without a segmented BAF — missing data, F59) and `Segments` (`AspcfSegment`: runASCAT logR runs within a chromosome, BAF of the first heterozygous locus or NaN, `LocusCount`, `HeterozygousLocusCount`) |
| IReadOnlyList\<AspcfSegmentation\> | list | `SegmentAlleleSpecificAsMultiPcf` (F53): one germline-aware segmentation per sample (`Tumor_LogR_segmented[, s]`, mirrored `Tumor_BAF_segmented[[s]]`), breakpoints common to the samples except where `refine` removed them |
| PurityPloidyFit | record | ρ, ASCAT output ploidy (probe-weighted mean integer total CN; over all probes from `FitPurityPloidyFromAspcf`, F37), GoF %, integer `AlleleSpecificSegment`s (one per summary / per `AspcfSegment`, ASCAT `seg_raw`), `Psi` (ψ), `IsNonAberrant` |
| DeriveMultiplicity | int | integer multiplicity m ∈ [1, majorCopyNumber] |
| SubclonalSegmentFit | record | per-segment primary/secondary `SubclonalCopyNumberState` (major, minor, cellFraction) + IsSubclonal flag |

### 3.3 Preconditions and Validation

Positions are 0-based. Null `loci`/`segments` → `ArgumentNullException`; empty/malformed `segments` (or no autosomal
segment) → `ArgumentException`; out-of-range thresholds, grid bounds, or multiplicity arguments →
`ArgumentOutOfRangeException`; no acceptable ASCAT optimum → `InvalidOperationException` from `FitPurityPloidy`
(`TryFitPurityPloidy` returns false). `SegmentAlleleSpecific` / `SegmentAlleleSpecificAspcf(loci, penalty)` (the
all-heterozygous fastAspcf path) reject a locus with a non-finite logR or a BAF outside [0, 1] (`ArgumentException`). The
`ascat.aspcf` / `ascat.asmultipcf` drivers (`SegmentAlleleSpecificAspcf(loci, germlineHeterozygous, …)`,
`SegmentAlleleSpecificAsMultiPcf`) accept NaN as R's `NA` in logR and heterozygous BAF and follow R's NA path (B24 F59);
±∞ logR, and a heterozygous BAF that is ±∞ or outside [0, 1], are rejected. BAF is
mirrored about 0.5 (ascat.aspcf `ifelse(b > 0.5, b, 1 − b)`) during segmentation so the two symmetric het clusters reinforce.

## 4. Algorithm

### 4.1 High-Level Steps

1. **Segment (ASPCF, `ascat.aspcf`):** per contiguous chromosome run, minimise the MAD-standardised joint logR +
   mirrored-BAF squared error plus γ per breakpoint (kmin 6; γ = 70 for `SegmentAlleleSpecific`, B24 F35). Summarise
   each segment by mean raw logR and the ASPCF mirrored BAF. Segmenting on BAF as well as logR is essential because copy-neutral LOH (e.g. 2:0) shares a balanced region's
   logR but not its BAF.
2. **Fit (`runASCAT`):** build the distance matrix over the autosomal segments, collect strict 7 × 7 local minima
   through the four-pass filter cascade, keep the smallest distance (ρ > 1 ⇒ 1), emit the `seg_raw` integer
   segments for every summary (sex chromosomes with the diploid model unless an `AscatSexModel` says XY, F38), the
   ASCAT ploidy and GoF. No candidate ⇒
   rho = NA. `EvaluatePurityPloidy` is the rho_manual/psi_manual path.
   **With homozygous probes (B24 F37, `FitPurityPloidyFromAspcf` / `TryFitPurityPloidyFromAspcf` /
   `EvaluatePurityPloidyFromAspcf` on an `AspcfSegmentation`):** the fit uses the heterozygous autosomal probes only
   (`r = lrrsegmented[names(bafsegmented)]`) grouped by `make_segments` (runs of identical (r, b), length = het probes);
   `seg_raw` has one row per logR segment with `bafke` = BAF of its first heterozygous probe, and `bafke = 0` when the
   segment has none (NA), which after the negative-value correction puts the whole total on nA (only nA + nB matters);
   ploidy = `mean(nA + nB)` over all probes (heterozygous and homozygous, ascat.runAscat l. 98).
   **Sex chromosomes (B24 F38, `AscatSexModel` overloads of all six fit entry points):** the fit always excludes X/Y
   (`autoprobes = !(chr %in% sexchromosomes)`, any gender). For XX (default) `haploidchrs` is empty and X/Y use the
   diploid equations. For XY, X and Y are haploid: `nAraw = (ρ − 1 + (2(1 − ρ) + ρψ)·2^(r/γ))/ρ`, `nBraw = 0` (ψ = the
   selected grid ψ / `psi_manual`), then the usual negative-value correction and R rounding. With an X non-PAR interval
   (`diploidprobes_fixnonPAR`) an X segment is haploid only if its closed overlap with the interval exceeds 50 % of
   its closed [first probe, last probe] span (IRanges widths); otherwise it is diploid (PAR). Y stays haploid.
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
   edge (`orderEdges` option 1); clonal if the closest corner is within 0.01 BAF or (with SNP BAFs) the t-test against
   its level gives p > 0.05, else two states with fraction τ [8].
6. **Phased-BAF segmentation (Battenberg `segment.baf.phased`, B24 F40):** `SegmentPhasedBaf` (presegments, MAD sd,
   two `selectFastPcf` passes, `calc_seg_baf_option`) → `BuildBattenbergSegments` (`determine_copynumber`
   switchpoints: new segment where BAFseg or chromosome changes; LogR = mean finite logR within [startpos, endpos],
   0 when none) → `FitSubclonalCopyNumberWithSnpTest` [8].
7. **Alternative solutions + bootstrap CIs (B24 F41, `FitSubclonalCopyNumberWithBootstrap`):** after step 5, for a
   sub-clonal segment the six `orderEdges` options (NA rows last), τ, SDfrac, and per option `noperms` bootstrap
   resamples (R RNG stream: segments in order, options A–F, NA options included) → SDfrac_BS and the 25th/975th order
   statistics [8].
   7a. **Merge + mask (B24 F60, `MergeBattenbergSegments` / `MaskHighCopyNumberSegments`):** the `merge_segments` loop of
   §2.2 on the fitted rows (output chromosomes in `GenomeInfoDb` seqlevel order — the order `makeGRangesFromDataFrame`
   gives, e.g. "10", "2", "X" → 2, 10, X — which reorders `BAFsegmented` for unsorted input); `BuildBattenbergSegments`
   of the merged rows feeds the second fit; masking flags rows with a state-1 copy number above `max_allowed_state` [8].
   7b. **`callSubclones` driver (B24 F61, `CallBattenbergSubclones`):** one `RMersenneTwister(seed)` → fit (step 7) →
   merge → fit of the merged rows on the same stream → mask; returns the initial calls (`_1.txt`), the merge and the
   masked final calls (plots, files, PGA-clonal and the ploidy recalculation are not part of it) [8].
8. **Multi-sample segmentation (ASCAT `ascat.asmultipcf`, B24 F53, `SegmentAlleleSpecificAsMultiPcf`):** per chromosome
   part, each sample's winsorised logR (every probe) and mirrored winsorised BAF (heterozygous probes; weight 0 at
   homozygous probes) form 2·S tracks; `ASmultiPCFcompact` minimises Σ_segments Σ_tracks −(Σw·y)²/Σw + γ·#breaks on the
   raw values (fewer than 6 probes ⇒ one segment); with `refine`, each sample is re-segmented on the compacted joint
   segments with γ/S; BAF per run = 0.5 + |b − 0.5| (shrunk to 0.5 by the 2·sd rule), logR levels → mean raw logR;
   ladder (penalty, 25 … 800) while any sample has ≥ 800 levels [2][7].
   **Missing data (R `NA`, B24 F59; NaN in the input):** `Select_sites` = probes with a (heterozygous) BAF or a logR in
   at least one sample — only these enter the joint segmentation, each with weight 0 for every missing value;
   `Select_sites2` = probes with both in at least one sample get a segmented BAF (a sample's missing BAF there takes the
   joint segment's value). Winsorisation (`madWinsMatrixWithNA`) and the per-site logR averages skip NA (with a single
   selected site R winsorises its BAF across the samples); a part with no selected site gets `mean(lr, na.rm = TRUE)`;
   NaN/0 levels take the closest level of the part, then the genome-wide re-estimation uses `na.rm = TRUE` and carries the
   previous level over an all-missing run (`prevlevel`). `ascat.aspcf` has its own NA path, ported in the same way:
   `Select_het = !homo & !is.na(baf) & !is.na(lr)`, NA skipped in the winsorisation, gap breakpoint (`na.rm`), level means
   and the homozygous-stretch `exactPcf` (`pcfed` = 0 at NA, `!anyNA(dif)`), then `fillNA` and `prevlevel`. runASCAT then
   fits on the probes with a segmented BAF (`SNPposhet = SNPpos[names(bafsegmented), ]`, `bafke` = first of them) and
   `ploidy = mean(nA + nB, na.rm = TRUE)` (NA at missing raw logR, and at a probe with a segmented but missing raw BAF).
   Every R `mean()` of the ASPCF ports (`mean(logRaveraged)`, `mean(bafselwinsmirrored)`, fastAspcf `yhat1 =
   mean(logR[frst:last])` and `mu = mean(abs(yi2 − 0.5))`, exactPcf `mean(y)`, and the `na.rm = TRUE` level / averaging
   means of `ascat.aspcf` / `ascat.asmultipcf`) is R's `real_mean` reproduced bit for bit — long-double Σx/n plus the
   residual-refinement pass, emulated exactly by `StatisticsHelper.ExtendedPrecisionMean` (B24 F66); the R-locked segment
   levels are therefore bit-identical (formerly ≤ 7.8e−16 from a plain double Σx/n).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- γ = 1 (sequencing) [2]; balanced (BAF = 0.5, exact equality as in R) segments weighted ×0.05 [2]; worst-case
  integer distance 0.25 [2]; ASCAT constants MINRHO 0.2, MINGOODNESSOFFIT 80, MINPERCZERO 0.02, MINPERCZEROABB 0.1,
  MINPERCODDEVEN 0.05, MINPLOIDYSTRICT 1.7, MAXPLOIDYSTRICT 2.3, MINABB 0.03, MINABBREGION 0.005 [2].
- Battenberg constants maxdist 0.01, siglevel 0.05, cn_upper_limit 1000, minimum minor CN 0.01 [8].
- Segments with End == Start are emitted with a 1 bp span so `AlleleSpecificSegment.Length > 0`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ASPCF segmentation | O(L·W) per chromosome (W = 1000-locus window) | O(W) | windowed PCF DP [2][6] |
| Multi-sample ASPCF (exact) | O(N²·2S) per chromosome part (N probes, S samples) | O(N·2S) | `ASmultiPCFcompact` DP; `Fast`: 5001-probe windows + compacted DP [2][7] |
| Purity/ploidy fit | O(P·Q·S) | O(P·Q + S) | P,Q = grid sizes (≤ 4·10⁶ cells), S = segments |
| Multiplicity | O(1) | O(1) | closed form |
| Sub-clonal fit | O(S) | O(S) | closed-form decomposition per segment |
| Bootstrap CIs | O(6·noperms·Σnₛ) | O(noperms + n) | nₛ = SNPs of sub-clonal segment s [8] |
| merge_segments | O(S²·n_c) per chromosome (S segments, n_c SNPs/probes of the chromosome; ranges found by linear scan) | O(n_c) | each iteration merges or closes one side [8] |
| Phased-BAF segmentation | O(C²) per presegment, C = filterMarkS4 candidates (≈ 15–30 % of SNPs); 5000-SNP windows above 15 000 | O(n) | Potts DP on the compacted array [8] |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.CopyNumberPloidy.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.CopyNumberPloidy.cs)

- `OncologyAnalyzer.SegmentAlleleSpecific(...)`: `SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty = 70)` (B24 F35). It formerly ran an unsourced greedy mean-shift heuristic; that heuristic was removed and its threshold parameters are now ignored (kept for source compatibility, still range-checked).
- `OncologyAnalyzer.SegmentAlleleSpecificAspcf(...)`: ASCAT ASPCF (`ascat.aspcf` port) on heterozygous loci.
- `OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, germlineHeterozygous, penalty)`: complete `ascat.aspcf` with germline genotypes (homozygous probes, homozygous-stretch resegmentation) → `AspcfSegmentation` (B24 F36).
- `OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples, germlineHeterozygous, penalty = 70, AsMultiPcfOptions?)`: joint multi-sample ASPCF (`ascat.asmultipcf` port: `ASmultiPCFcompact`, `compactASMulti`, `runFastASMultiPCF`, `expandMulti`, per-sample refinement) → one `AspcfSegmentation` per sample, each usable by `FitPurityPloidyFromAspcf` (B24 F53).
- `OncologyAnalyzer.FitPurityPloidy(...)` / `TryFitPurityPloidy(...)`: ASCAT `runASCAT` fit → ρ, ψ, ploidy, GoF, integer segments.
- `OncologyAnalyzer.EvaluatePurityPloidy(...)`: ASCAT rho_manual/psi_manual path.
- `OncologyAnalyzer.FitPurityPloidyFromAspcf(...)` / `TryFitPurityPloidyFromAspcf(...)` / `EvaluatePurityPloidyFromAspcf(...)`: runASCAT on a germline-aware `AspcfSegmentation` (homozygous segments, all-probe ploidy; B24 F37).
- `AscatSexModel` / `AscatGender` + an `AscatSexModel` overload of each of the six fit entry points above: runASCAT `gender` / `X_nonPAR` (B24 F38).
- `AscatMaleXGenotyping(sexModel, seed?, germlineBaf?)` with `SegmentAlleleSpecificAspcf(loci, germlineHeterozygous, maleX, penalty)` and `AsMultiPcfOptions.MaleXGenotyping`: the male `X_nonPAR` germline re-genotyping of `ascat.aspcf` / `ascat.asmultipcf` (`set.seed(seed)`; `rank(DIST, ties.method = "random")` with germline BAF, `sample()` without), reusing the private `RMersenneTwister` (B24 F58). `SegmentAlleleSpecificAspcfSamples(samples, seed?, penalty)` (`AspcfSampleInput`: loci, genotypes, per-sample `SexModel` / `GermlineBaf`): `ascat.aspcf` over several samples after one `set.seed`, the re-genotyping draws continuing on one R stream from male sample to male sample (B24 F65).
- `OncologyAnalyzer.DeriveMultiplicity(...)`: McGranahan multiplicity (rounded, clamped).
- `OncologyAnalyzer.FitSubclonalCopyNumber(...)`: Battenberg `determine_copynumber` (nearest edge, τ, maxdist).
- `OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(SubclonalSegmentSnpBafs[], ρ, ψ, γ, siglevel = 0.05, maxdist = 0.01)` →
  `SubclonalSegmentTestedFit` (fit + Battenberg pval): the per-SNP t-test (B24 F39; t-test =
  `StatisticsHelper.OneSampleTTestPValue`).
- `OncologyAnalyzer.SegmentPhasedBaf(PhasedBafSnp[], BattenbergPhasedSegmentationOptions?)` → `PhasedBafSegmentedSnp[]`
  (Battenberg `segment.baf.phased`; private port of `selectFastPcf`/`filterMarkS4`/`PottsCompact`/`runPcfSubset`, shared
  `GetMad`; type-7 quantile = `StatisticsHelper.SampleQuantileType7`) and
  `OncologyAnalyzer.BuildBattenbergSegments(rows, LogRProbe[])` → `SubclonalSegmentSnpBafs[]` (B24 F40).
- `OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(SubclonalSegmentSnpBafs[], ρ, ψ, seed, γ, siglevel, maxdist,
  permutations = 1000)` → `BattenbergSegmentCall` (Fit = solution A, pval, BAF l, ntot, `BattenbergSubclonalSolution` A–F
  with SDfrac / SDfrac_BS / frac1_0.025 / frac1_0.975); private `RMersenneTwister` = R's default RNG (B24 F41).
- `OncologyAnalyzer.MergeBattenbergSegments(calls, rows, LogRProbe[], ρ, ψ, γ, bafOption)` → `BattenbergSegmentMerge`
  (`BattenbergMergedSegment` rows + updated `BAFsegmented`) and `OncologyAnalyzer.MaskHighCopyNumberSegments(calls, rows,
  maxAllowedState = 250)` → `BattenbergHighCopyNumberMask` (`BattenbergSegmentCall.IsMasked`, masked rows,
  `masked_count`, `masked_size`); Welch test = `StatisticsHelper.WelchTTestPValue` (R `t.test(x, y)`); seqlevel order =
  internal `BattenbergSeqlevelOrder` (GenomeInfoDb `rankSeqlevels`) (B24 F60).
- `OncologyAnalyzer.CallBattenbergSubclones(rows, LogRProbe[], ρ, psit, seed, γ, siglevel, maxdist, permutations,
  maxAllowedState, bafOption)` → `BattenbergSubcloneCalls(InitialCalls, Merge, Mask)` (`Calls` = final table); shares the
  private `BattenbergDetermineCopyNumber` (RNG passed in) with `FitSubclonalCopyNumberWithBootstrap` (B24 F61).

### 5.2 Current Behavior

ASCAT fit (runASCAT port, R-verified 150/150; germline-aware / homozygous-segment, male X/Y and NA paths, F36–F38/F58/F59);
ASPCF = `ascat.aspcf` port (R-verified 60/60), multi-sample `ascat.asmultipcf` (F53); sub-clonal = Battenberg
`determine_copynumber` port (R-verified 402/402) plus the per-SNP t-test, phased-BAF segmentation, bootstrap,
`merge_segments` / `mask_high_cn_segments` and `callSubclones` driver (F39–F41, F60, F61). `SegmentAlleleSpecific` runs the same ASPCF
at ASCAT's default penalty 70 (the former greedy heuristic was removed, B24 F35). Not a search/matching task, so the repository
suffix tree is **not used** (no occurrence enumeration).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ASCAT `runASCAT`: nA/nB equations, distance matrix (probe-count weights, autosomes, genome-wide minor allele),
  local-minimum scan, filter cascade, ρ clamp, `seg_raw` rounding, GoF, non-aberrant flag, manual (ρ, ψ) path [2].
- runASCAT `gender` / `X_nonPAR` (B24 F38): haploid male X/Y `seg_raw` equations, `diploidprobes_fixnonPAR`, XX
  default — R-verified on 3 male genomes (germline-aware and het-only paths, XY / XY without non-PAR / XX) [2].
- runASCAT with germline-homozygous probes (B24 F37): het-only `make_segments` fit, `bafke` NA ⇒ 0 for homozygous
  segments, `ploidy = mean(nA + nB)` over all probes — R-verified end to end (ascat.aspcf → runASCAT) on 3 genomes [2].
- McGranahan observed mutation copy number n_mut and the [1, major] multiplicity clamp [3][4].
- ASCAT ASPCF: winsorisation, MAD-standardised joint cost, kmin 6, windows, BAF shrinkage, penalty ladder [2][6][7].
- ASCAT `ascat.aspcf` with germline genotypes (B24 F36): heterozygous-only BAF, logR from all probes averaged onto the
  heterozygous probes, gap breakpoint search, `predictGermlineHomozygousStretches` + `exactPcf` resegmentation,
  `fillNA`, genome-wide level re-estimation — R-verified on 3 genomes (26 segments, ≤ 5.6e-16 per locus) [2].
- ASCAT `ascat.asmultipcf` (B24 F53): joint exact (`ASmultiPCFcompact`) or fast (`runFastASMultiPCF`) segmentation of
  S logR + S mirrored-BAF tracks on the raw cost, `wsample` weights, per-sample `refine` (γ/S), BAF shrinkage, closest-probe
  fill of 0 levels, genome-wide level re-estimation, penalty ladder — R-verified on 5 cohorts (2–4 samples, 800–6500
  probes; private breakpoints, homozygous stretches) in 14 runs, every per-probe logR/BAF ≤ 7.8e-16, breakpoints
  identical; R quirks kept (NaN `bestCost[1]` for a zero-weight first column, fast-path unweighted window sums) [2][7].
- ASCAT NA path of `ascat.asmultipcf` / `ascat.aspcf` + runASCAT (B24 F59): R-verified on 3 NA cohorts (scattered NA in
  one sample, whole-chromosome NA runs, NA BAF at heterozygous probes; 2–3 samples, 1032–1505 probes) — 14 segmented
  tracks with identical breakpoints and NA positions and every per-probe logR/BAF within 1e-12; runASCAT purity, ψ,
  ploidy, goodness of fit and seg_raw identical on all 12 fitted tracks (R finds no solution for the 13th/14th) [2].
- Battenberg `determine_copynumber`: nearest edge (`orderEdges` option 1), τ, maxdist clonality, negative-minor
  adjustment [8]; per-SNP t-test clonality (`t.test(BAFke, mu = test.level)`, siglevel) via
  `FitSubclonalCopyNumberWithSnpTest` — R-verified 180/180 segments on 4 genomes (B24 F39) [8].

**Intentionally simplified:**

- ~~ASCAT inputs are segment summaries: germline-homozygous probes (their logR averaging and the homozygous-stretch
  resegmentation of `ascat.aspcf`) … are not available~~ — **resolved by F36** (germline-aware
  `SegmentAlleleSpecificAspcf` overload). ~~The gender-specific haploid X/Y model is not available; sex-chromosome
  segments are … emitted with the diploid (gender "XX") model~~ — **resolved by F38** (`AscatSexModel` overloads:
  runASCAT `gender = "XY"` haploid X/Y, `X_nonPAR` > 50 % overlap rule, ASCAT hg19/hg38/CHM13 constants; R-verified on
  3 male genomes, 14 runs). Sex-chromosome segments are still excluded from the fit, exactly as in ASCAT. ~~Not ported:
  the male-only `X_nonPAR` germline re-genotyping inside `ascat.aspcf` (random draw)~~ — **resolved by F58**
  (`AscatMaleXGenotyping`: all non-PAR X probes homozygous, then `round(m · h_auto)` re-marked heterozygous — the
  smallest germline-BAF distance to 0/1 with R's runif tie-break, or `sample()` without germline BAF — from R's seeded
  Mersenne-Twister; R-verified on 4 aspcf tracks and 3 asmultipcf cohorts: selected probes identical, segments identical,
  levels ≤ 1e−12; bit-identical since F66).
- Sub-clonal fit: ~~a summary carries no per-SNP BAF spread, so Battenberg's t-test cannot be run~~ — **resolved by
  F39** (`FitSubclonalCopyNumberWithSnpTest` takes the phased SNP BAFs; the summary overload still decides by maxdist,
  exactly Battenberg's result for constant BAF). Where R's `t.test` would stop ("data are essentially constant") and
  abort Battenberg, pval 0 is used. ~~Haplotype phasing (`BAFphased` is the caller's input)~~ — **phased-BAF path resolved
  by F40**: `SegmentPhasedBaf` ports `segment.baf.phased` (haplotype-block correction + segmentation) on caller-phased
  BAFs, R-verified on 6 tracks (300–16 000 SNPs, all three `selectFastPcf` regimes, prior breakpoints, ≥ 3 Mb gaps,
  options 1/2/3, no_segmentation): every segment extent identical and BAFseg ≤ 1e−12; end-to-end with
  `determine_copynumber` 9/9 segments identical. The phasing itself (IMPUTE2/Beagle5 + 1000 Genomes panel) is not
  run (BLOCKED: external executables + multi-GB reference bundle). ~~Bootstrap CIs and alternative solutions B–F are not
  produced~~ — **resolved by F41** (`FitSubclonalCopyNumberWithBootstrap`): R-verified on 7 single segments (NA solution,
  noperms 500/20 ⇒ NA bounds, 1 SNP ⇒ SDfrac NA, seeds 0/−12345/2³¹−1) and the 2 multi-segment F40 tracks — every
  column ≤ 1e−12 including the bootstrap ones (R's RNG reproduced; resamples identical); 100 other seeds agree with R's
  400-seed Monte-Carlo means within 5 standard errors. ~~`merge_segments` / `mask_high_cn_segments` not ported~~ —
  **resolved by F60** (`MergeBattenbergSegments`, `MaskHighCopyNumberSegments`): R-verified (GenomicRanges 1.54.1) on 4
  genomes covering every merge branch (same clonal state, same square with non-significant / significant Welch tests,
  too few values, different squares, > 3 Mb), options 1/2/3 (incl. a median of exactly 1 ⇒ mean), chromosome reordering
  and two masked segments (max_allowed_state 5 and the BAF = 1 ⇒ cn_upper_limit segment at 250): merged extents and
  states identical, BAF/LogR ≤ 1e−12, `BAFsegmented` runs identical, masked count/size/rows identical. ~~`callSubclones`
  calls `determine_copynumber` twice after one `set.seed`; a seeded run matches only the second call's state~~ —
  **resolved by F61** (`CallBattenbergSubclones`): on the same 4 genomes both tables (initial and final masked, all 67
  columns incl. SDfrac_BS and the bootstrap bounds of the second fit) agree with R ≤ 1e−12 (pval 1e−9 rel).

**Not implemented:**

- ~~Multi-sample (asmultipcf) segmentation~~ — **resolved by F53** (`SegmentAlleleSpecificAsMultiPcf`); its missing-data
  (NA) path is ported since F59 and the male-only `X_nonPAR` random re-genotyping since F58
  (`AsMultiPcfOptions.MaleXGenotyping`). The former "whole-genome-doubling refit search" item was removed: ASCAT has no such
  procedure (`grep -rniE "wgd|whole.?genome.?doubl|refit" ASCAT/R` finds only the `ascat.metrics` WGD status, ported in
  F30), so it had no reference counterpart.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | ~~Greedy mean-shift segmentation retained alongside ASPCF~~ | Deviation | — | **resolved (B24 F35)** | `SegmentAlleleSpecific` now delegates to ASPCF (`ascat.aspcf`, penalty 70) [2]; legacy thresholds ignored |
| 2 | Segment summaries instead of probes | Assumption | ~~no homozygous-probe logR~~ (resolved F36), ~~no haploid X/Y model~~ (resolved F38) | **resolved (B24 F36/F38)** | see §5.3; the male non-PAR re-genotyping is ported too (F58) |
| 3 | ~~No per-SNP t-test in the sub-clonal fit~~ | Assumption | ~~clonality by maxdist only~~ | **resolved (B24 F39)** | `FitSubclonalCopyNumberWithSnpTest` runs Battenberg's `t.test(BAFke, mu = test.level)`; the summary overload is Battenberg's constant-BAF branch [8] |
| 4 | ~~`PurityPloidyFit.Ploidy` = probe-weighted mean integer CN over heterozygous probes~~ | Assumption | ASCAT averages over all probes | **resolved (B24 F37)** | `FitPurityPloidyFromAspcf` averages over all probes; summary-based `FitPurityPloidy` has only heterozygous loci, where both coincide |

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
| Multi-sample ASPCF: one sample, or a 1-probe chromosome part | `ArgumentException` (R stops: "incorrect number of subscripts on matrix") | source [2] (R run, F53) |
| Multi-sample ASPCF: breakpoint supported by one sample only | joint breakpoint kept in that sample, removed from the others by `refine` (γ/S) | source [2] (R run, F53) |
| Multi-sample ASPCF with refine: chromosome part without heterozygous probes | its joint breakpoint is always merged away (R: `bestCost[1]` = NaN for the zero-BAF-weight first block) | source [2] (R run, F53) |
| Missing logR/BAF (NaN = R NA) in `ascat.asmultipcf` / `ascat.aspcf` | R's NA path: weight 0 / `na.rm`, NA-free probes selected, no segmented BAF where a probe has no BAF (or no logR) in any sample (aspcf: in that sample); whole-chromosome NA → `prevlevel` (asmultipcf) or `fillNA` (aspcf); ±∞ still `ArgumentException` | source [2] (R run, F59) |
| Last averaging window of a part of 2, 32, 128, 16384 … probes | `lr[n:(n−1)]` = probes n and n − 1 (R descending range; IEEE `floor(n + 0.01 − 0.01)`) | source [2] (R run, F53; also fixed in the F36 path) |
| Male (`AscatSexModel` XY): X/Y segment | haploid: nB = 0, nA = round((ρ − 1 + (2(1−ρ)+ρψ)·2^(r/γ))/ρ), negative ⇒ 0:0 | source [2] (R run, F38) |
| Male with X non-PAR: X segment overlapping non-PAR by exactly 50 % | diploid (rule is strictly > 0.5) | source [2] (`diploidprobes_fixnonPAR`) |
| Sub-clonal: BAF within 0.01 of a corner | single clonal state, f=1 | INV-07 [8] |
| merge_segments: Welch t-test on essentially constant data | `InvalidOperationException` (R `t.test` stops) | source [8] (F60) |
| merge_segments: ±∞ logR; chromosome without calls or logR probes | `ArgumentException` (R fails: non-finite means/tests; `stopifnot`) | source [8] (F60) |
| mask_high_cn_segments: first SNP of a masked segment | keeps its BAFseg (`startpos < Position`) | source [8] (R run, F60) |

### 6.2 Limitations

logR and BAF are observed measurements and are always a caller input — this is inherent, not a limitation of the
derivation. The unit works on heterozygous-locus segment summaries (germline-homozygous probes: use the germline-aware
ASPCF overload, F36); the haploid X/Y (male) model is available through `AscatSexModel` (F38). Battenberg's per-SNP t-test is available through
`FitSubclonalCopyNumberWithSnpTest` (F39; phased SNP BAFs supplied by the caller). Multi-sample segmentation is available
through `SegmentAlleleSpecificAsMultiPcf` (`ascat.asmultipcf`, F53; missing logR/BAF as R's NA since F59). ASCAT's male `X_nonPAR`
germline re-genotyping is available through `AscatMaleXGenotyping` (F58); it is random, so it reproduces R only for an
explicit seed (ASCAT's default seed is `as.integer(Sys.time())`, exposed as `AscatMaleXGenotyping.DefaultSeed()` and used when `seed` is null); a multi-sample `ascat.aspcf` run (one seed, one stream across samples) is reproduced by `SegmentAlleleSpecificAspcfSamples` (F65). Battenberg's built-in haplotype imputation (IMPUTE2/Beagle5 against the 1000 Genomes reference panel —
external executables and a multi-GB reference bundle) is out of scope; the downstream phased path is available
(`SegmentPhasedBaf` → `BuildBattenbergSegments` → `FitSubclonalCopyNumberWithSnpTest`, F40) on caller-phased BAFs; alternative
solutions B–F and the seeded bootstrap CIs are available through `FitSubclonalCopyNumberWithBootstrap` (F41). Battenberg's
`merge_segments` / `mask_high_cn_segments` post-processing is available through `MergeBattenbergSegments` /
`MaskHighCopyNumberSegments` (F60), and the whole copy-number sequence of `callSubclones` (one seeded RNG stream) through
`CallBattenbergSubclones` (F61). `FitPurityPloidy` fails (like ASCAT) when no local
minimum passes the filters — use `TryFitPurityPloidy`, or `EvaluatePurityPloidy` with externally chosen (ρ, ψ).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var loci = /* per-locus (chrom, pos, logR, BAF) measurements */;
var summaries2 = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci); // ASCAT ASPCF segmentation (penalty 70)
// SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2) returns the same (legacy name; thresholds ignored, F35)
var fit = OncologyAnalyzer.FitPurityPloidy(summaries2);         // → ρ, ψ, integer segments (throws if ASCAT finds none)
// With germline genotypes (true = heterozygous): full ascat.aspcf + runASCAT incl. homozygous probes (F36/F37)
var aspcf = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, germlineHeterozygous);
var fitAll = OncologyAnalyzer.FitPurityPloidyFromAspcf(aspcf);  // ploidy = mean(nA + nB) over all probes
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
- Tests: [StatisticsHelper_StudentT_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/StatisticsHelper_StudentT_Tests.cs) — Student t / incomplete beta / one-sample t-test vs R (F39); Battenberg t-test rows in `OncologyAnalyzer_AscatDerivation_Tests.cs` (F39), R-locked
- Tests: [OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_BattenbergPhasedSegmentation_Tests.cs) — `segment.baf.phased` on 6 tracks + end-to-end `determine_copynumber` (F40), R-locked; [StatisticsHelper_QuantileType7_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/StatisticsHelper_QuantileType7_Tests.cs) (F40)
- Tests: [OncologyAnalyzer_BattenbergBootstrap_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_BattenbergBootstrap_Tests.cs) — solutions A–F, SDfrac, bootstrap CIs (F41), R-locked + Monte-Carlo agreement
- Tests: [OncologyAnalyzer_BattenbergMergeMask_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_BattenbergMergeMask_Tests.cs) — `merge_segments` / `mask_high_cn_segments` (F60), R-locked on 4 genomes; seqlevel order vs GenomeInfoDb; [OncologyAnalyzer_BattenbergCallSubclones_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_BattenbergCallSubclones_Tests.cs) — `callSubclones` driver (F61), R-locked; [StatisticsHelper_WelchTTest_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/StatisticsHelper_WelchTTest_Tests.cs) (F60)
- Tests: [OncologyAnalyzer_AscatAsMultiPcf_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatAsMultiPcf_Tests.cs) — multi-sample `ascat.asmultipcf` (F53), R-locked on 5 cohorts / 14 runs
- Tests: [OncologyAnalyzer_AscatMissingData_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatMissingData_Tests.cs) — NA path of `ascat.asmultipcf` / `ascat.aspcf` + runASCAT (F59), R-locked on 3 NA cohorts (6 segmentations, 12 fits)
- Tests: [OncologyAnalyzer_AscatMaleXNonPar_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatMaleXNonPar_Tests.cs) — male X non-PAR re-genotyping in aspcf / asmultipcf (F58), R-locked on 7 tracks
- Tests: [OncologyAnalyzer_AscatSexChromosome_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_AscatSexChromosome_Tests.cs) — male haploid X/Y, X non-PAR, XX default (F38), R-locked
- Evidence: [ONCO-ASCAT-001-Evidence.md](../../../docs/Evidence/ONCO-ASCAT-001-Evidence.md)
- Related algorithms: [Tumor_Ploidy_Estimation](./Tumor_Ploidy_Estimation.md), [Cancer_Cell_Fraction_Estimation](./Cancer_Cell_Fraction_Estimation.md), [Tumor_Purity_Estimation](./Tumor_Purity_Estimation.md)

## 8. References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. 2010. Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R`, `ASCAT/R/ascat.aspcf.R` (master, read 2026-09-28), `ASCAT/R/ascat.loadData.R` (`gender`, `X_nonPAR` constants; read 2026-10-10). https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. 2022. Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ. 2021. DeCiFering the elusive cancer cell fraction. PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
6. Nilsen G, Liestøl K, Van Loo P, et al. 2012. Copynumber: Efficient algorithms for single- and multi-track copy number segmentation. BMC Genomics 13:591. https://doi.org/10.1186/1471-2164-13-591
7. Ross EM, Haase K, Van Loo P, Markowetz F. 2021. Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37(13):1909–1911. https://doi.org/10.1093/bioinformatics/btaa538
8. Nik-Zainal S, Van Loo P, Wedge DC, et al. 2012. The Life History of 21 Breast Cancers. Cell 149(5):994–1007. https://doi.org/10.1016/j.cell.2012.04.023 ; Battenberg `R/fitcopynumber.R` (`determine_copynumber`; `callSubclones`, `merge_segments`, `mask_high_cn_segments`, read 2026-10-10 at 57a8f7e), `R/orderEdges.R` (master, read 2026-09-28), https://github.com/Wedge-lab/battenberg
