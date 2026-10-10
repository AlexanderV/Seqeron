# Evidence Artifact: ONCO-ASCAT-001

**Test Unit ID:** ONCO-ASCAT-001
**Algorithm:** Upstream derivation of allele-specific copy-number segments, joint purity/ploidy fit (ASCAT), and mutation multiplicity
**Date Collected:** 2026-06-23

---

## Online Sources

### Van Loo et al. (2010) — ASCAT, "Allele-specific copy number analysis of tumors" (PNAS)

**URL:** https://www.pnas.org/doi/10.1073/pnas.1009843107 (full text PDF mirror: https://unclineberger.org/peroulab/wp-content/uploads/sites/1008/2019/06/Aug-23-ASCAT-aCGH-VanLoo-PNAS-2010.pdf)
**Accessed:** 2026-06-23 (PDF retrieved and read; pages 2–4 examined directly)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points (from the retrieved PDF, Fig. 1 caption and main text):**

1. **Two input tracks:** SNP arrays (and, with γ=1, sequencing) deliver **Log R** (total signal intensity, "r") and **B-allele frequency** (BAF, "b", allelic contrast). For each segment a single fitted logR value is obtained; BAF gives one or two values per segment.
2. **Grid search over ψ and ρ:** "ASCAT first determines the ploidy of the tumor cells ψ_t and the fraction of aberrant cells ρ. This procedure evaluates the goodness of fit for a grid of possible values for both parameters" (Fig. 1 caption). The blue/red "sunrise" plot shows goodness of fit over (ploidy, aberrant cell fraction); the optimal solution (green cross) minimises distance of the allele-specific copy numbers to non-negative integers.
3. **Integer-closeness objective:** "ASCAT evaluates a plurality of possible combinations of tumor ploidy and tumor fractions, based on the assumption that the associated allele-specific copy number calls should be as close as possible to nonnegative whole numbers for germline heterozygous SNPs."
4. **Ploidy scale:** ploidy is "the amount of DNA relative to a haploid genome"; a pure diploid genome → ψ = 2; ">2.7n" marks aneuploidy (Fig. 2 text).

### ASCAT R reference implementation — `ascat.runAscat.R` (VanLoo-lab/ascat, master)

**URL:** https://raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.runAscat.R
**Accessed:** 2026-06-23 (raw source fetched; core fitting lines quoted)
**Authority rank:** 3 (reference implementation of the primary paper, by the original authors)

**Key Extracted Points (verbatim R lines from the source):**

1. **Allele-specific copy number from (r, b, ρ, ψ, γ):**
   ```r
   nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
   nB = (rho-1 +  b   *2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
   ```
   where `r` = segment logR, `b` = segment BAF, `rho` = aberrant-cell fraction (purity ρ), `psi` = tumour ploidy ψ, `gamma` = platform parameter.
2. **Goodness-of-fit distance** (minimised over the (ρ, ψ) grid):
   ```r
   d[i,j] = sum( abs(nMinor - pmax(round(nMinor),0))^2 * length * ifelse(b==0.5, 0.05, 1), na.rm=TRUE )
   ```
   i.e. segment-length-weighted **squared distance of the minor allele to the nearest non-negative integer**, with BAF=0.5 (balanced) segments down-weighted ×0.05.
3. **Theoretical maximum distance and percentage GoF:**
   ```r
   TheoretMaxdist = sum( rep(0.25, n) * length * ifelse(b==0.5, 0.05, 1), na.rm=TRUE )
   goodnessOfFit  = (1 - m/TheoretMaxdist) * 100
   ```
   (0.25 = (½)² is the worst-case distance to an integer; GoF reported as a percentage.)
4. **Integer assignment:** per-probe plotting values use `nA = pmax(round(nAfull),0)`; the **segment** output (`seg`/`seg_raw`) first corrects negative values (`nA+nB<0 ⇒ 0,0`; a negative allele is added to the other), then rounds with R `round` (half-to-even) and, for BAF = 0.5 segments, applies `limitround = 0.5` (odd total ⇒ nA+1 or nB−1). With ASCAT's segmented BAF ≤ 0.5, nA is the major allele. *(Corrected 2026-09; the earlier reading omitted the segment rules.)*

### ASCAT README — γ (gamma) platform parameter

**URL:** https://github.com/VanLoo-lab/ascat/blob/master/README.md (and Crick / MD Anderson Van Loo lab software pages)
**Accessed:** 2026-06-23
**Authority rank:** 3

**Key Extracted Points:**

1. **γ for sequencing = 1:** "For massively parallel sequencing data, gamma should always be set to 1"; "for HTS data (WGS, WES and TS), gamma must be set to 1 in ascat.runASCAT." (Default γ=0.55 is for SNP arrays only.)

### Nilsen et al. (2012) — Copynumber / PCF + ASPCF (BMC Genomics 13:591) [ASPCF segmentation half]

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3582591/
**Accessed:** 2026-06-23 (PMC HTML full text fetched)
**Authority rank:** 1 (peer-reviewed methods paper)

**Key Extracted Points (verbatim from the retrieved text):**

1. **PCF penalised-least-squares criterion:**
   `L(S | y, γ) = Σ_{I∈S} Σ_{j∈I} (y_j − ȳ_I)² + γ·|S|` — within-segment SSE plus a penalty `γ > 0` per
   segment (`|S|` = number of segments); γ trades goodness-of-fit against parsimony.
2. **Dynamic-programming recurrence (global optimum, O(n²)):**
   `e_k = min_{j ∈ {1,…,k}} ( d_{jk} + e_{j−1} + γ )`, `e_0 = 0`, where `d_{jk}` is the within-segment SSE of the
   run `j..k`.
3. **ASPCF / multi-track joint cost:** `L(S | y₁, y₂, γ) = L(S | y₁, γ) + L(S | y₂, γ)` — a single segmentation
   (common breakpoints) with **separate per-track segment means**; the per-segment cost is the sum of the two
   tracks' SSE and γ is charged once per segment.
4. **Default penalty:** "A fairly conservative penalty of γ = 40 is the default in the copynumber package."

### Ross et al. (2021) — Allele-specific multi-sample segmentation in ASCAT (Bioinformatics 37:1909) [ASPCF half]

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8317109/
**Accessed:** 2026-06-23 (PMC HTML full text fetched)
**Authority rank:** 1

**Key Extracted Points:**

1. **Joint allele-specific objective (verbatim):**
   `L(S|Y,W,γ) = Σ_i Σ_{I∈S} Σ_{j∈I} w_{ij}(y_{ij} − ȳ_{i,I})² + γ|S|` — logR and BAF jointly segmented with
   **common change points**, each track keeping its own segment mean.
2. **BAF mirroring (verbatim):** preprocessing includes "mirroring BAFs to obtain a single track in regions of
   allelic imbalance" — confirms folding BAF to its distance from 0.5 before joint segmentation.

### Nik-Zainal et al. (2012) — Battenberg sub-clonal copy number (Cell 149:994) [sub-clonal half]

**URL:** https://www.cell.com/cell/fulltext/S0092-8674(12)00527-2 ;
https://github.com/Wedge-lab/battenberg/blob/master/README.md
**Accessed:** 2026-06-23 (Battenberg README fetched; Cell paper via search-result confirmation)
**Authority rank:** 1 / 3

**Key Extracted Points (Battenberg README, verbatim):**

1. **Two-state model:** "Each segment in the tumour genome will have either one or two copy number states. If there
   is one state it represents the clonal copy number (i.e. all tumour cells have this state); if there are two
   states it represents subclonal copy number (i.e. there are two populations of cells, each with a different
   state). A copy number state consists of a major and a minor allele and their frequencies, which together add
   give the total copy number for that segment and an estimate fraction of tumour cells that carry each allele."
2. **Output columns:** `nMaj1_A, nMin1_A, frac1_A` (state 1 CN + tumour-cell fraction), `nMaj2_A, nMin2_A,
   frac2_A` (state 2; NA for clonal). `frac1 + frac2 = 1`.
3. **Decomposition** *(corrected 2026-09 from the Battenberg source, R/fitcopynumber.R `determine_copynumber` +
   R/orderEdges.R)*: the two states are the corners of the **nearest edge** of the copy-number square around
   (nMajor, nMinor) — they differ in one allele — chosen by the segment BAF l relative to the corner BAF levels and
   the total-copy-number priority `ntot < x + y + 1`; the fraction of state 1 is the BAF-mixture solution
   `τ = (1−ρ+ρM₂−2l(1−ρ)−lρ(m₂+M₂)) / (lρ(m₁+M₁)−lρ(m₂+M₂)−ρM₁+ρM₂)` (unclamped). A segment is clonal when
   |l − closest corner level| < maxdist = 0.01 or the per-SNP t-test is not significant (constant SNP BAF ⇒ pval 0).
   The earlier "both alleles bracketed with one shared least-squares fraction" reading was not Battenberg's model.

### McGranahan et al. (2016) — clonal neoantigens (Science) [already cited in repo]

**URL:** https://www.science.org/doi/10.1126/science.aaf1490 (search-result snippet retrieved 2026-06-23; full formula already cited in OncologyAnalyzer.cs line 6965)
**Accessed:** 2026-06-23
**Authority rank:** 1

**Key Extracted Points:**

1. **Expected allele frequency:** AF_expected = p·M / (p·C_t + C_n·(1−p)), where p = purity, M = mutated-allele copy number (multiplicity), C_t = tumour copy number, C_n = normal copy number (= 2). Equivalently the observed mutation copy number n_mut = VAF·(1/p)·[p·C_t + C_n·(1−p)], and CCF = n_mut / M.

### Zheng et al. (2022) — PICTograph (Bioinformatics) [already cited in repo]

**URL:** https://academic.oup.com/bioinformatics/article/38/15/3677/6596597
**Accessed:** 2026-06-23 (fetched; generative-model VAF equation quoted)
**Authority rank:** 1

**Key Extracted Points:**

1. **VAF generative model (verbatim):** `VAF = (m × CCF × p) / (c × p + 2 × (1 − p))`, where m = multiplicity, CCF = cancer cell fraction, p = purity, c = tumour total copy number.
2. **Multiplicity by inversion:** at clonal CCF = 1, m = VAF·(c·p + 2(1−p)) / p; rounding to the nearest integer and clamping to [1, major-allele CN] gives the integer multiplicity (this is the McGranahan n_mut rounded — n_mut = m when CCF = 1).

### DeCiFering (2021, PMC8542635) — CCF closed form

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
**Accessed:** 2026-06-23 (fetched)
**Authority rank:** 1

**Key Extracted Points:**

1. **CCF closed form (verbatim):** c = (F·v)/(ρ·M) with F = ρ·N_tot + 2(1−ρ); v = VAF, M = multiplicity. Confirms the CCF formula already implemented in `EstimateCcf`.

---

## Documented Corner Cases and Failure Modes

### From ASCAT (Van Loo 2010 / source)

1. **Balanced (BAF = 0.5) segments** carry little allele-specific information and are down-weighted ×0.05 in the goodness-of-fit (they cannot distinguish e.g. 1+1 from 2+2 except via logR).
2. **Non-identifiability / multiple optima:** the sunrise plot can show several local minima (e.g. a 2n vs 4n solution). *(Corrected 2026-09.)* ASCAT does **not** take the global grid minimum: `runASCAT` keeps strict minima of a 7×7 window that pass a 4-pass filter cascade (ploidy in (min_ploidy, max_ploidy), ρ ≥ 0.2, GoF > 80 %, percentzero > 0.02; fallbacks with strict ploidy 1.7–2.3, perczeroAbb > 0.1, percOddEven > 0.05 and ρ>1 columns masked) and selects the smallest distance; when none passes it returns rho = NA.
3. **γ must match the platform:** γ=1 for sequencing, ≈0.55 for arrays; a wrong γ rescales logR and biases copy number.

### From McGranahan / PICTograph

1. **Multiplicity clamp:** rounded multiplicity must be clamped to [1, major-allele CN]; a raw value < 0.5 would round to 0 (no mutated copy) which is non-physical for an observed variant → clamp to ≥ 1.

---

## Test Datasets

### Dataset: Planted-truth synthetic genome (deterministic)

**Source:** Synthesised by inverting the ASCAT forward model (Van Loo 2010 equations) — for KNOWN ρ₀, ψ₀ and integer (nA, nB) per segment, compute the per-locus logR r and BAF b that ASCAT would observe, then run the derivation and assert recovery.

Forward model (algebraic inverse of the two nA/nB equations, γ=1):
- total tumour copy number n = nA + nB
- BAF b = (ρ·nB + (1−ρ)·1) / (ρ·n + (1−ρ)·2)   (B-allele fraction at a germline-het SNP: tumour contributes nB of n, normal contributes 1 of 2)
- logR r = log2( (ρ·n + (1−ρ)·2) / ψ_overall ),  where the average total ploidy used to normalise is ρ·ψ_tumour + 2(1−ρ) so that a balanced reference segment has r = 0.

| Parameter | Value |
|-----------|-------|
| ρ₀ (purity) | 0.80 |
| ψ₀ (tumour ploidy, nA+nB length-weighted) | 2.0 (diploid planted) and 3.0 (triploid planted) |
| γ | 1 (sequencing) |
| Segment A | nA=1, nB=1 (balanced diploid), b=0.5, r=0 |
| Segment B | nA=2, nB=0 (copy-neutral LOH), b→ extreme, r computed |
| Segment C | nA=2, nB=1 (gain), b, r computed |
| Planted clonal mutation | VAF synthesised from m=1 on a CN=2 (1+1) segment at CCF=1 |

---

### Dataset: Planted two-level logR track (ASPCF breakpoint recovery)

> **Superseded 2026-09:** ASCAT's ASPCF places no breakpoint on noise-free data (MAD sd = 0); the unit tests now use the R-verified noisy tracks of §"2026-09 review".

**Source:** synthesised per Nilsen et al. 2012 PCF objective (deterministic).

| Parameter | Value |
|-----------|-------|
| Track | logR = 0.0 for loci 0–9, logR = 1.0 for loci 10–19 (one true breakpoint at index 10) |
| BAF | 0.5 throughout (balanced) |
| γ | 0.5 (small enough that 1 breakpoint beats 0: ΔSSE between merged and split = 25 ≫ γ) |
| Expected | 2 segments; breakpoint between position 9 and 10; means 0.0 and 1.0 |

### Dataset: Planted sub-clonal segment (mixture recovery)

> **Superseded 2026-09:** Battenberg mixes the two corners of the nearest edge; for this input it returns (2,0)@0.142857 + (2,1)@0.857143 (§"2026-09 review").

**Source:** Battenberg two-state model (Nik-Zainal 2012).

| Parameter | Value |
|-----------|-------|
| ρ (purity) | 1.0 |
| ψ (ploidy) | 2.0 ; γ = 1 |
| True mixture | nA_obs = 0.4·2 + 0.6·1 = 1.4 ; nB_obs = 0.4·0 + 0.6·1 = 0.6 (total 2.0) |
| Expected fit | states (2,0) and (1,1); fraction f ≈ 0.4 (within 0.05) |
| Pure-clonal control | nA_obs = 2, nB_obs = 1 → single state, f ≈ 0 or 1 |

---

## Assumptions

1. **ASSUMPTION: Germline-heterozygous-SNP BAF forward model.** The per-locus BAF at a germline heterozygous SNP is b = (ρ·nB + (1−ρ)) / (ρ·(nA+nB) + 2(1−ρ)) — the B-allele copies (tumour nB + one normal copy) over total copies. This is the standard ASCAT allelic-contrast model implied by the nA/nB inversion; it is used only to **synthesise planted-truth test inputs**, not in the production derivation (the production code consumes measured logR/BAF). Justified because it is the exact algebraic inverse of the two cited ASCAT equations.
2. **ASSUMPTION: logR normalisation reference = average sample ploidy.** Planted logR uses r = log2( (ρ·n + 2(1−ρ)) / (ρ·ψ + 2(1−ρ)) ). This makes a segment at the genome-average copy number have r = 0, matching ASCAT's tumour-baseline-corrected logR. Used only for planted-truth synthesis.

---

## Recommendations for Test Coverage

1. **MUST Test:** Allele-specific segmentation recovers planted breakpoints from per-locus logR/BAF — Evidence: ASCAT segmentation (Van Loo 2010 §, integer-closeness over segments).
2. **MUST Test:** Joint (ρ, ψ) grid fit recovers planted ρ₀=0.80, ψ₀∈{2,3} within tolerance and the integer (nA,nB) per segment — Evidence: ASCAT nA/nB equations + goodness-of-fit (source lines).
3. **MUST Test:** Derived multiplicity equals planted m for a clonal mutation — Evidence: McGranahan n_mut rounding / PICTograph inversion.
4. **MUST Test:** End-to-end CCF on a planted clonal mutation ≈ 1.0 — Evidence: DeCiFering / McGranahan CCF closed form.
5. **SHOULD Test:** GoF percentage at the true (ρ,ψ) is higher (distance lower) than at a deliberately wrong (ρ,ψ) — Rationale: the objective must actually discriminate.
6. **SHOULD Test:** Null/empty inputs and invalid (ρ,ψ) grid bounds throw — Rationale: contract robustness.
7. **COULD Test:** A balanced-only genome (all 1+1) yields BAF≈0.5 segments down-weighted in GoF — Rationale: documented corner case.
8. **MUST Test (ASPCF):** ASPCF recovers the planted single breakpoint on a two-level logR track with a sourced γ — Evidence: Nilsen 2012 PCF objective + DP recurrence.
9. **MUST Test (ASPCF):** On a constructed noisy track, the ASPCF penalised cost ≤ the cost of any other segmentation (DP global optimum) — Evidence: Nilsen 2012 (DP returns the global minimum). *(B24 F35: the greedy `SegmentAlleleSpecific` comparator no longer exists — the method now delegates to ASPCF at penalty 70; the R-locked M-ASPCF-1/2 outputs pin the optimum.)*
14. **MUST Test (F35):** `SegmentAlleleSpecific(loci, …)` ≡ `SegmentAlleleSpecificAspcf(loci, 70)`: on the R-locked noisy step track it returns the ascat.aspcf output (2 segments, loci 1–40 / 41–80, logR −0.023564999999999999 / 0.60321000000000002, BAF 0.5 / 0.75129407874999998), independent of the ignored legacy thresholds; NaN/±∞ logR or BAF ∉ [0, 1] → `ArgumentException` (ASPCF input contract).
10. **MUST Test (ASPCF):** Mirrored-BAF joint cost separates a copy-neutral-LOH segment from a balanced segment that share logR — Evidence: ASCAT joint segmentation (Ross 2021).
11. **MUST Test (ASPCF):** Large γ collapses to a single segment; small γ recovers each level — Evidence: Nilsen 2012.
12. **MUST Test (sub-clonal):** Sub-clonal fit recovers planted f₀ = 0.4 with states (2,0)/(1,1) within tolerance — Evidence: Battenberg two-state model.
13. **MUST Test (sub-clonal):** A pure-clonal (integer) segment collapses to a single state (f ≈ 0 or 1) — Evidence: Battenberg (one state = all tumour cells).

### ASPCF / sub-clonal assumptions

A. **ASSUMPTION: γ exposed as a sourced parameter rather than hard-coded.** The penalty *form* (`+ γ·|S|`) and DP
   recurrence are sourced verbatim (Nilsen 2012). The numeric default (copynumber γ = 40; ASCAT later 70) is
   probe-scale-specific; the repository API operates on caller-supplied summary scales, so γ is a required exposed
   parameter and tests use a γ derived from each dataset's ΔSSE so the optimum is provable.
B. **ASSUMPTION: two-state mixture uses the two bracketing integers.** A single fractional value `n_obs` has a
   unique two-state mixture with `f ∈ [0,1]` using `⌊n_obs⌋, ⌈n_obs⌉`. Three-or-more-population mixtures and
   non-adjacent states are out of scope (documented limitation).

---

## 2026-09 review (B24 — F12, F13, F14)

### Sources opened (2026-09-28)

| Source | What it confirmed |
|---|---|
| `raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.runAscat.R` | `runASCAT`, `make_segments` (length = probe count), `create_distance_matrix` (grid seq(min_ploidy−0.5, max_ploidy+0.5, 0.05) × seq(min_purity, max_purity, 0.01); genome-wide minor-allele choice), local-minimum scan + filter cascade + constants, ρ>1 ⇒ 1, `seg_raw` rounding, autosome-only fit, rho_manual/psi_manual path |
| `.../ASCAT/R/ascat.aspcf.R` | `ascat.aspcf` (penalty 70, ladder 35/50/70/100/140 while ≥ 800 levels, MAD winsorisation, < 6 loci ⇒ one segment), `fastAspcf` (1000/100 windows, sd validity, BAF shrinkage), `aspcfpart` (kmin 6, standardised costs), `getMad`, `madWins`, `medianFilter` |
| `.../ASCAT/R/ascat.metrics.R` | ASCAT's WGD / ploidy metrics (not used by this unit) |
| `raw.githubusercontent.com/Wedge-lab/battenberg/master/R/fitcopynumber.R`, `R/orderEdges.R`, `R/clonal_ascat.R` | `callSubclones` / `determine_copynumber` (maxdist 0.01, siglevel 0.05, cn_upper_limit 1000, τ formula), `orderEdges` |
| R 4.3 `stats::runmed`, `stats::smoothEnds` (printed source) | running median with `endrule = "median"` (ends kept, then Tukey end rule) |

### ASCAT R cross-check (executed, not transcribed)

R 4.3 with the original files sourced; plotting functions stubbed. `runASCAT(lrr, baf, lrrsegmented, bafsegmented,
"XX", SNPpos, ch, chrs, c("X","Y"), FALSE, …, gamma = 1, min_ploidy 1.5, max_ploidy 5.5, min_purity 0.1,
max_purity 1.05)` on one probe per heterozygous locus (a segment = LocusCount probes with its r and 1 − mirrored b);
integer segments from the verbatim `seg_raw` block per segment. `ascat.aspcf(obj, penalty, out.dir = NA)` with
Germline_BAF = 0.5 (all heterozygous). Battenberg: the verbatim `determine_copynumber` per-segment body with
`orderEdges` sourced.

| Check | Cases | Result (C# vs R) | Pre-fix code vs R |
|---|---|---|---|
| `FitPurityPloidy` / `TryFitPurityPloidy` (ρ, ψ, GoF, nonaberrant, integer segments) | 150 random genomes (3–9 segments, ρ 0.25–1, doubled genomes, noise, chrX, 2 NA) | 150/150 identical (max |Δ| 4.3e-14) | 123/148 differ (ρ/ψ in 123, segments in 54) |
| `SegmentAlleleSpecificAspcf` (breakpoints, logR, BAF) | 60 genomes, 681 segments, up to 4 500 loci/chromosome, penalties 5–150 | 60/60 identical (max |Δ| 5.0e-16) | 52/60 breakpoint sets differ (83 vs 681 segments) |
| `FitSubclonalCopyNumber` | 402 segments (random + exact clonal, BAF 1, unmirrored BAF) | 402/402 identical | — (different model) |

Locked reference values (unit tests): runASCAT case A (chrX, 3:1…) ρ = 1, ψ = 2.7, GoF 99.781420571107006, 2:1 ×3;
case B ρ = 0.85, ψ = 2.2, GoF 99.999772627448223, 2:0 2:1 1:1 (pre-fix 0.72 / 4.45); two genomes with rho = NA;
ascat.aspcf noisy step: logR −0.023565 / 0.60321, BAF 0.5 / 0.75129407875; LOH track BAF 0.966975; penalty 1e6
single segment logR 0.2898225, BAF 0.637905789375; noise-free step ⇒ one segment; 5-locus chromosome logR 0.13,
BAF 0.774; Battenberg (1,1)@0.70000000000000051 + (2,1)@0.29999999999999949; (2,0)@0.14285714285714241 + (2,1).

### Notes

- ASCAT cannot segment noise-free data (every window has MAD sd = 0 and is skipped); planted-truth tests of the
  segmenter therefore use deterministic noisy tracks.
- A single-segment genome never yields an accepted ASCAT optimum (no strict 7×7 minimum passes the filters): runASCAT
  returns NA for all 21 single-segment edge inputs of the fuzz suite.
- Battenberg reproduces a floating-point artefact at exact-integer inputs (e.g. ρ = 1, BAF 0.66666666666666674,
  logR log2(1.5): ntot rounds to 3 ⇒ sub-clonal (3,0)/(3,1) with τ = −0.5); the port reproduces it bit-for-bit.

## 2026-10 FIN-B24 F26 — DeriveMultiplicity tie rule (facets-suite)

- Source opened: mskcc/facets-suite `R/ccf-annotate-maf.R` (master @ 7d54d0f6, raw.githubusercontent.com),
  `expected_mutant_copies(t_var_freq, total_copies, purity)` — "Based on PMID 28270531":
  `mu = t_var_freq·(1/purity)·(purity·total_copies + (1−purity)·2)`; `alt_copies = ifelse(mu < 1, 1, abs(mu))`;
  `round(alt_copies)` (R `round` = IEC 60559 half-to-even). `total_copies == 0 → 1` (Seqeron rejects N_T < 1 instead).
  facets-suite does not cap at the major copy number — Seqeron's cap at `majorCopyNumber` is a documented extra.
- Same quantity as `DeriveMultiplicity` (n_mut = `AdjustVAFForPurity`, identical formula); floor at 1 is equivalent
  (`round(mu) ≤ 1` for mu < 1). Only exact .5 ties differ from the former away-from-zero rounding.
- R 4.3.3 output of `expected_mutant_copies` (sourced from the fetched file; all inputs dyadic ⇒ mu exact):
  (0.625, 4, 1.0) mu 2.5 → 2; (0.375, 2, 0.5) mu 1.5 → 2; (0.75, 2, 1.0) 1.5 → 2; (0.875, 4, 1.0) 3.5 → 4;
  (0.5625, 8, 1.0) 4.5 → 4; (0.125, 4, 1.0) 0.5 → 1; (0.3, 2, 1.0) 0.59999999999999998 → 1;
  (0.55, 4, 0.8) 2.4750000000000001 → 2. C# `Math.Round(·, MidpointRounding.ToEven)` reproduces all.

## 2026-10 FIN-B24 F36 — `ascat.aspcf` with germline-homozygous probes

- Source opened: `raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.aspcf.R` — `ascat.aspcf`
  (`gg = ascat.gg$germlinegenotypes`, TRUE = homozygous; `Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)`;
  logR winsorised over all probes, BAF over the heterozygous probes; `averageIndices`/`startindices`/`endindices`
  midpoint averaging onto each heterozygous probe; `fastAspcf` (≥ 6 het probes) or the mean (< 6); "find best
  breakpoint" in homozygous gaps; no-het chromosome ⇒ `mean(lr)` logR segment; "correct wrong segments in germline
  homozygous stretches" (± 100 PCF context, ± 5 replacement margin, `exactPcf(winsed, 6, floor(segmentlength/4))`,
  `dif > 0.3`, `sum(dif > 0.3) > 5`); `fillNA(logRPCFed, zeroIsNA = TRUE)`; "adapt levels again"),
  `predictGermlineHomozygousStretches` (`homthres = ceiling(log(0.001, perchom))`), `exactPcf`, `fillNA`.
- R quirks kept verbatim: `(1:bp)` at bp = 0 is `c(1, 0)` (the gap search compares probes `at+1` and `at` with the
  left level); `fillNA` stretches of 2–3 NAs get `(end+1):end` (only the last element takes the next value);
  `endindices = floor(length(lr) + 0.01 − 0.01)` evaluated in IEEE doubles; genome-wide `rle` re-estimation crosses
  chromosome boundaries when two adjacent levels are bit-equal.
- R cross-check (executed): R 4.3.3, `ascat.aspcf.R` + `ascat.runAscat.R` sourced (plots/`png`/`dev.*` stubbed),
  `ascat.aspcf(obj, ascat.gg = list(germlinegenotypes = !het), penalty = 70, out.dir = NA)` with `chr = ch` = one part per
  chromosome. Input: three deterministic genomes (1080 / 1002 / 1100 probes; 608 / 600 / 525 heterozygous) generated by
  the MINSTD script below (bit-identical in C#, `OncologyAnalyzer_AscatGermlineHomozygous_Tests.Simulate`):

  ```r
  mk <- function(seed) { s <- seed; function() { s <<- (s * 48271) %% 2147483647; s / 2147483647 } }
  # per probe: ug, a1..a4, us, c1..c4 = 10 draws; hom <- forced || ug < pHom
  # logR = r + 0.2*((((a1+a2)+a3)+a4) - 2); het BAF = min(1, max(0, (us<0.5 ? 1-b : b) + 0.1*((((c1+c2)+c3)+c4) - 2)))
  # hom BAF = (us < 0.5 ? 0 : 1); position = 1000 * index
  ```
  Genomes (blocks: chromosome, probes, r, mirrored b, forced homozygous): G1 seed 11, pHom 0.3 — 1:150 −0.189/0.5,
  1:100 0.2439/0.6296, 1:30 hom, 1:50 hom −1.926 (homozygous deletion inside a 110-probe homozygous stretch), 1:30 hom,
  1:100, 2:150 −0.8105/0.7692, 2:150 0.5765/0.5, 3:120 hom −0.189 (no heterozygous probe), 4:100 0.5765/0.7059,
  4:100 −0.189/0.5. G2 seed 22, pHom 0.25 — CN breakpoint inside an 80-probe homozygous gap on chr1, a 12-probe chr3,
  homozygous-only chr4, chrX with a 30-probe homozygous tail. G3 seed 33, pHom 0.35 — 60-probe focal gain (0.2615)
  inside a 180-probe homozygous stretch on chr2 (−0.6645), homozygous-only chr4.
- Result: per-locus `Tumor_LogR_segmented` C# vs R max |Δ| 3.3e-16 / 1.7e-16 / 5.6e-16, segmented BAF ≤ 2.2e-16, NA
  pattern identical; runASCAT logR segments 11 / 7 / 8 identical (extents, levels, first-het BAF). The homozygous
  deletion (G1 probes 281–330, level −1.9130208073104829) and the focal gain (G3 481–540, 0.26558251092493651) are
  produced only by the homozygous-stretch resegmentation; G2's breakpoint lands at probe 240 (end of the planted
  block) by the gap search.
- Gap shown by the same data (pre-F36 API): dropping homozygous probes (het-only `SegmentAlleleSpecificAspcf`) loses the
  deletion/focal gain and shifts the logR levels (G3 fit ρ 0.91 / ψ 3.2 instead of ASCAT's 0.92 / 3.15); passing them
  as heterozygous (BAF 0/1) corrupts the BAF segmentation (G1 fit ρ 0.86 / ψ 3.25 instead of 0.7 / 2.45).
- All-heterozygous input: the new overload reproduces the het-only overload bit for bit on the 3 genomes' het subsets
  (the two R genome-wide steps it adds — `fillNA` of a level exactly 0, re-averaging of bit-equal adjacent levels across
  chromosomes — do not occur there).

## 2026-10 FIN-B24 F37 — runASCAT with homozygous segments

- Source opened: `raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.runAscat.R` — `runASCAT`:
  `b = bafsegmented; r = lrrsegmented[names(bafsegmented)]` (heterozygous probes only), `autoprobes`, `make_segments(r2, b2)`
  (runs of identical (r, b), no chromosome check; length = probes), distance matrix / TheoretMaxdist / nonaberrant /
  filter cascade on those segments only; `seg` loop over `rle(lrrsegmented)` ∪ chromosome ends with
  `bafke = bafsegmented[bafpos][1]` and "if bafke is NA … germline homozygous stretch … just their sum matters":
  `bafke = 0`; `ascat.runAscat`: `ploidy = mean(nA + nB, na.rm = TRUE)` over `n1all`/`n2all` = every non-NA probe
  (CN probes and homozygous probes get nMajor + nMinor).
- Consequences ported: homozygous probes never enter the distance matrix or the GoF (only via the logR levels); a
  segment without heterozygous probes gets nAraw = total + (1 − ρ)/ρ, nBraw = (ρ − 1)/ρ ≤ 0 ⇒ negative-value correction
  ⇒ (round(total), 0); ploidy is weighted by all probes.
- R cross-check (executed; same genomes and harness as F36, `ascat.runAscat(gamma = 1)`):

  | Genome | purity | psi | ploidy (all probes) | goodnessOfFit | seg_raw (nMajor, nMinor) | het-probe mean nA+nB |
  |---|---|---|---|---|---|---|
  | G1 | 0.7 | 2.45 | 2.4722222222222223 | 99.937918942624279 | 1:1 2:1 2:1 0:0 2:1 2:1 1:0 2:2 2:0 3:1 1:1 | 2.6134868421052633 |
  | G2 | 0.53 | 2.35 | 2.1816367265469063 | 99.634551649302338 | 1:1 3:0 1:0 2:1 2:2 1:0 1:1 | 2.2716666666666665 |
  | G3 | 0.92 | 3.15 | 3.3818181818181818 | 98.658353206884712 | 2:2 3:1 2:0 4:0 2:0 3:2 2:1 3:0 | 3.539047619047619 |

  C# `FitPurityPloidyFromAspcf`: ρ, ψ, ploidy and GoF bit-identical, every seg_raw row identical; nonaberrant FALSE.
  Homozygous-only segments: G1 chr3 (nAraw 2.0120388012516401, nBraw 0 ⇒ 2:0), G1 deletion 0:0, G3 focal 4:0, chr4 3:0.
- Gap shown (pre-F37): the summary fit on heterozygous loci reports ploidy over heterozygous probes only
  (G1 2.6134868421052633 vs ASCAT 2.4722222222222223) and cannot emit homozygous segments.

## 2026-10 FIN-B24 F38 — Sex-chromosome model (runASCAT `gender` / `X_nonPAR`)

- Sources opened: `raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.runAscat.R` — `runASCAT`:
  `autoprobes = !(SNPposhet[,1] %in% sexchromosomes)` (the fit excludes X/Y for **every** gender);
  `haploidchrs = unique(c(substring(gender,1,1), substring(gender,2,2)))`, minus the letter when both are equal
  (XX ⇒ none; XY ⇒ X, Y); `nullchrs = setdiff(sexchromosomes, …)` (XX ⇒ Y) but `nullprobes` is consulted only for
  non-diploid probes, so in XX Y uses the diploid equations; `seg_raw` for a non-diploid segment:
  `nAraw = (rho-1 + ((1-rho)*2+rho*psi)*2^(logR/gamma))/rho`, `nBraw = 0`, ψ = `psi_opt1` (grid value) or `psi_manual`,
  then the negative-value correction and rounding (`limitround` cannot fire with nBraw = 0);
  `diploidprobes_fixnonPAR` (only when `!is.null(X_nonPAR) && gender == "XY"`): X probes reset to diploid, X segments =
  `rle` runs of the segmentation, a segment is non-diploid when `width(pintersect(nonPAR, segment))/width(segment) > 0.5`.
  `ascat.loadData.R`: `gender = NULL ⇒ "XX"`; `X_nonPAR` = hg19 `c(2699521, 154931043)`, hg38 `c(2781480, 155701382)`,
  CHM13 `c(2394411, 153925834)`, `genomeVersion = NULL ⇒ X_nonPAR = NULL` (whole X haploid in a male).
  `ascat.aspcf.R`: for a male with `X_nonPAR` the non-PAR germline genotypes are re-drawn (all homozygous, then a random
  autosome-matched fraction heterozygous) — **not ported** (caller's genotypes are used).
- R cross-check (executed; R 4.3.3, sourced ascat.aspcf.R + ascat.runAscat.R, GenomicRanges 1.54.1; harness of F36 with
  per-block positions; `ascat.aspcf(X_nonPAR = NULL)`, then `gender` / `X_nonPAR` set, `ascat.runAscat(gamma = 1)`):

  | Genome / model | purity | psi | ploidy | goodnessOfFit | X / Y seg_raw (nMajor:nMinor, nAraw) |
  |---|---|---|---|---|---|
  | M1 XY (X_nonPAR NULL) | 0.7 | 2.45 | 2.181159420289855 | 99.937918942624279 | X 1:0 (1.03218649423988) 2:0 (1.8962760923968001) 2:0 (2.0480422819682); Y 1:0 1:0 0:0 (0.017976122283112) |
  | M1 XX | 0.7 | 2.45 | 2.1775362318840581 | 99.937918942624279 | X 1:0 1:0 1:1; Y 1:0 1:0 0:0 |
  | M2 XY non-PAR [890000, 1090000] | 0.53 | 2.35 | 2.0981228668941978 | 99.634551649302338 | X 873000–907000 3:0 (overlap 17001/34001), 908000–912000 3:0, 1:0, 2:0, 1078000–1122000 2:1 (12001/45001 ⇒ diploid); Y 1:0 |
  | M2 XX | 0.53 | 2.35 | 1.8805460750853242 | 99.634551649302338 | X 1:1 1:1 0:0 1:0 2:1; Y 0:0 |
  | M2 XY (X_nonPAR NULL) | 0.53 | 2.35 | 2.1365187713310578 | 99.634551649302338 | X 3:0 3:0 1:0 2:0 4:0 (4.0928384899644898); Y 1:0 |
  | M3 XY hg19 | 0.92 | 3.15 | 2.9482758620689653 | 98.658353206884712 | X PAR1 1:1 1:1, 3–152 Mb 1:0, 153–155.02 Mb 3:0 (2.75988943219894), PAR2 1:1; Y 2:0 0:0 |
  | M3 XX | 0.92 | 3.15 | 2.9482758620689653 | 98.658353206884712 | X 1:1 1:1 1:0 2:1 1:1; Y 2:0 0:0 |
  | M3 XY (X_nonPAR NULL) | 0.92 | 3.15 | 2.9482758620689653 | 98.658353206884712 | X 2:0 2:0 1:0 3:0 2:0; Y 2:0 0:0 |
  | het-only M1 XY / XX | 0.7 | 2.45 | 2.4579831932773111 (both) | 99.950135206703436 | X 2:0 vs 1:1; Y 1:0 vs 1:0 |
  | het-only M2 XY / XX | 0.53 | 2.35 | 2.3914529914529914 / 2.3384615384615386 | 99.661294751227899 | X 3:0 vs 1:1, PAR2 2:1 (both) |
  | het-only M3 XY / XX | 0.91 | 3.2 | 3.4351687388987568 (both) | 98.662279834887002 | one X segment 100000–155135000: 2:0 vs 1:1 |

  C# (`AscatSexModel` overloads of `FitPurityPloidyFromAspcf` / `FitPurityPloidy`): ρ, ψ, ploidy (≤ 1e-14), GoF and
  every seg_raw row identical in all 14 runs; ρ/ψ/GoF never depend on the sex model (fit excludes X/Y).
- Gap shown (pre-F38): X/Y were always emitted with the diploid model — e.g. M1 male X segment of 2 tumour copies came
  out 1:1 instead of ASCAT's 2:0 and the ploidy 2.1775362318840581 instead of 2.181159420289855.

---

## References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. (2010). Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R`, `ascat.aspcf.R`, `ascat.loadData.R` (master). https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. (2016). Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. (2022). PICTograph: estimation of cancer cell fractions and clone trees. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ (2021). DeCiFering the elusive cancer cell fraction. Cell Systems / PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
6. Nilsen G, Liestøl K, Van Loo P, et al. (2012). Copynumber: Efficient algorithms for single- and multi-track copy number segmentation. BMC Genomics 13:591. https://doi.org/10.1186/1471-2164-13-591 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3582591/)
7. Ross EM, Haase K, Van Loo P, Markowetz F (2021). Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37(13):1909–1911. https://doi.org/10.1093/bioinformatics/btaa538 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8317109/)
8. Nik-Zainal S, Van Loo P, Wedge DC, et al. (2012). The Life History of 21 Breast Cancers. Cell 149(5):994–1007. https://doi.org/10.1016/j.cell.2012.04.023 ; Battenberg, https://github.com/Wedge-lab/battenberg

---

## Change History

- **2026-10-10**: FIN-B24 F38 — runASCAT sex-chromosome model (male haploid X/Y, `X_nonPAR` rule, XX default) R cross-check section added.
- **2026-10-10**: FIN-B24 F37 — runASCAT with homozygous segments (`bafke` NA ⇒ 0, all-probe ploidy) R cross-check section added.
- **2026-10-10**: FIN-B24 F36 — germline-aware `ascat.aspcf` (homozygous probes, homozygous-stretch resegmentation) R cross-check section added.
- **2026-10-09**: FIN-B24 F35 — the unsourced greedy `SegmentAlleleSpecific` heuristic was removed; the public name now delegates to `SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty = 70)` (ascat.aspcf port, R-verified 60/60), with the ASPCF finite-logR / BAF ∈ [0, 1] validation; legacy thresholds ignored (range-checked for compatibility). Test point 14 added; point 9 reworded.

- **2026-06-23**: Initial documentation.
- **2026-06-23**: Added ASPCF penalised-least-squares segmentation (Nilsen 2012, Ross 2021) and sub-clonal copy-number two-state mixture (Nik-Zainal 2012 / Battenberg) evidence for the residual-closing fix.
- **2026-10-09**: FIN-B24 F26 — `DeriveMultiplicity` ties half-to-even per facets-suite `expected_mutant_copies`; F27 `SubclonalIntegerTolerance` `[Obsolete]`.
- **2026-09-28**: B24 review — FitPurityPloidy = runASCAT port, ASPCF = ascat.aspcf port, sub-clonal fit = Battenberg determine_copynumber port; corrected corner case 2, integer-assignment and Battenberg decomposition statements; R cross-check section added.
