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
  autosome-matched fraction heterozygous) — ~~not ported~~ ported by F58 (see § F58).
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

## 2026-10 FIN-B24 F39 — Battenberg per-SNP clonality t-test

- Source opened: `raw.githubusercontent.com/Wedge-lab/battenberg/master/R/fitcopynumber.R` — `determine_copynumber`
  (called twice by `callSubclones(..., siglevel = 0.05, maxdist = 0.01, ...)`) and `R/segmentation.R`:
  `BAFke = BAFphased[segment]` (per-SNP BAFs flipped to the segment side, `BAFphased = ifelse(BAFsegm > 0.5, BAF, 1 − BAF)`;
  **not** mirrored again), `l = max(BAFseg, 1 − BAFseg)`, `test.levels` = BAF of the two corners of the nearest edge,
  `whichclosestlevel.test = which.min(abs(test.levels − l))`;
  `if (is.na(sd(BAFke)) || sd(BAFke) == 0) pval = 0 else pval = t.test(BAFke, alternative = "two.sided", mu = test.levels[whichclosestlevel.test])$p.value`;
  `if (abs(l − test.level) < maxdist) pval = 1`; sub-clonal iff `pval <= siglevel`. < 2 SNPs ⇒ sd NA ⇒ pval 0; constant ⇒ 0;
  no tryCatch — R `t.test`'s "data are essentially constant" stop (stderr < 10·eps·|mean|) would abort Battenberg (port: pval 0).
- Student t: R `stats:::t.test.default` (stderr = sqrt(var/n), df = n − 1, p = 2·pt(−|t|, df)); R nmath `pt.c`
  (pbeta(1/(1 + t²/n), n/2, ½)), `lbeta.c`, `toms708.c` (`brcomp` large-parameter prefactor). New
  `StatisticsHelper.RegularizedIncompleteBeta` / `StudentTCdf` / `OneSampleTTestPValue` (continued fraction in
  double-double). Sweep vs R 4.3.3: `pt` 400 cases (ν 1–3·10⁵, |t| 1e−3–40): median 1.1e−16, max 1.7e−14 for P > 1e−30,
  max 1.7e−13 in the far tail (P ~ 1e−175); `pbeta` 300 cases (a, b 0.01–3·10⁴): 90 % ≤ 8.7e−15, max 1.5e−13. scipy 1.18.1
  `t.cdf` / `betainc` / `ttest_1samp` agree to the same digits (e.g. pt(−6, 2·10⁵) R 9.8827494396739874e−10, scipy
  9.882749439674027e−10, exact (mpmath) 9.8827494396740285e−10).
- R cross-check (executed; `determine_copynumber` sourced verbatim + `orderEdges.R`, R 4.3.3, noperms 1000): 4 genomes
  (ρ/ψ 0.8/2.5, 1/2, 0.55/3.1, 0.35/1.9), 45 segments each (40 random: integer or planted-mixture states, BAF offsets
  0/±0.004/0.012/−0.015, 1–200 SNPs, SNP sd 0–0.15; 5 designed): **180/180 identical** (pval rel ≤ 1.2e−13, states equal,
  fractions ≤ 1e−13). The maxdist-only `FitSubclonalCopyNumber` disagrees with Battenberg on 32/180.

  | Segment | SNP BAFs | R pval | Battenberg call | maxdist-only call |
  |---|---|---|---|---|
  | G1S41 (ρ 0.8, ψ 2.5, l = level(2,1) + 0.03) | 0.5529 0.7629 0.5729 0.7829 0.6929 | 0.56090425422510204 | clonal 2:1 | sub-clonal |
  | G1S42 (same, tight SNPs) | 0.6689 0.6759 0.6739 0.6709 0.6749 | 2.1020858614373501e−05 | (2,0)@0.15605024544591201 + (2,1) | same |
  | G1S43 / S44 (constant / 1 SNP) | 0.6729 ×4 / 0.6729 | 0 / 0 | (2,0)@0.15605… + (2,1) | same |
  | G1S45 (l − level 0.005) | 0.6379 0.6579 0.6484 | 1 (maxdist) | clonal 2:1 | same |
  | G1S28 | 0.8864 0.7955 0.8486 0.7765 0.617 | 0.067325006438243504 | clonal 2:0 | sub-clonal |
  | G1S40 | 8 SNPs | 0.0354791289104208 | (3,0)@0.69556883269516401 + (3,1) | same |
  | G3S40 (ρ 0.55, ψ 3.1) | 0.7938 0.7866 0.8057 | 0.036949670267513897 | (3,0)@0.84149606261457699 + (3,1) | same |
  | G4S27 (ρ 0.35, ψ 1.9) | 0.693 0.7098 0.7167 0.7037 0.6761 | 0.024885230920557801 | (2,0)@0.554567052328231 + (3,0) | same |

  Parameter variants (single-segment R runs): G1S28 with siglevel 0.1 ⇒ (2,0)@0.61666732876527997 + (2,1)@0.38333267123471998;
  G1S45 with maxdist 0.001 or 0 ⇒ pval 0.46228196278297101, clonal; BAFseg 1 − 0.672857 with SNPs 1 − G1S42 (unmirrored)
  ⇒ pval 1.7441626943676701e−09, sub-clonal (2,0)@0.15605024544591201.
- Gap shown (pre-F39): only maxdist decided clonality (a summary has no SNP spread ⇒ Battenberg pval 0), so noisy segments
  just beyond 0.01 BAF were always sub-clonal (G1S41, G1S28).

## 2026-10 FIN-B24 F40 — Battenberg phased-BAF segmentation (`segment.baf.phased`); built-in imputation BLOCKED

- Sources opened (Wedge-lab/battenberg master, commit 57a8f7e, 2025-12-08): `R/segmentation.R` (`segment.baf.phased`:
  `bkps_to_presegment_breakpoints` with `maxsnpdist = 3000000` and `addin_bigholes`; `run_pcf`: `sdev = getMad(ifelse(BAF < 0.5,
  BAF, 1 − BAF), k = 25)`, NA ⇒ 0, `< 0.09` ⇒ 0.09; `length(BAF) < 50` ⇒ mean; `selectFastPcf(BAF, phasekmin, phasegamma·sdev, T)`;
  `BAFphased = ifelse(BAFsegm > 0.5, BAF, 1 − BAF)`; `selectFastPcf(BAFphased, kmin, gamma·sdev, T)`; `calc_seg_baf_option`
  1/2/3 via `adjustSegmValues` (rle median); defaults gamma 10, phasegamma 3, kmin 3, phasekmin 3, option 3), `R/fastPCF.R`
  (`selectFastPcf`: < 1000 ⇒ `runFastPcf(…, 0.15, 0.15)`, < 15000 ⇒ `runFastPcf(…, 0.12, 0.05)`, else `runPcfSubset`;
  `filterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9)`, `compact`, `PottsCompact`, `findEst`, `markWithPotts`,
  `findMarks`, `getMad`, `medianFilter`). Battenberg's `getMad`/`medianFilter` are identical to ASCAT's (reused);
  `selectFastPcf` is **not** ASCAT's `exactPcf`/`fastAspcf` (ASCAT's ascat.aspcf.R has no `selectFastPcf`), so it is
  ported once (private) — no duplication.
- **Haplotype imputation BLOCKED (proof).** `R/impute.R::run_haplotyping` calls `run.impute` → `system("impute2 -m
  <genetic_map> -h <impute_hap> -l <impute_legend> -g <input> -int … -Ne 20000 -phase …")` per 5 Mb region, or
  `run.beagle5` → `system("java -Xmx10g -jar beagle.jar gt=… ref=<chrN.1kg.phase3.v5a.b37.bref3> map=<plink.chrN.GRCh37.map>
  impute=false")`. Both are external binaries driven by the README's "Required reference files": GRCh37 bundle
  (ora.ox.ac.uk uuid:2c1fec09…) `battenberg_1000genomesloci2012_v3.tar.gz`, `battenberg_impute_1000G_v3.tar.gz`,
  `probloci_270415.txt.gz`, GC/replication-timing correction tarballs; GRCh38 bundle (doi 10.48420/30406441)
  `1000G_loci_hg38.zip`, `imputation.zip`, `shapeit2.zip`, `beagle5.zip`, … . Sizes: a public mirror of the same
  bundle (bcgsc.ca morinlab/reference, via web search) lists `battenberg_impute_grch37.tar.gz` 3.7 GB and
  `battenberg_1000genomesloci_grch37.tar.gz` 240 MB; direct downloads from this sandbox are refused (proxy HTTP 403 for
  ora.ox.ac.uk / doi.org / figshare.manchester.ac.uk). A statistical phasing engine + 1000 Genomes haplotype panel
  cannot be shipped in or run by a pure C# library ⇒ the caller supplies phased BAFs (the `combine.baf.files` output).
- R cross-check (executed; `fastPCF.R` + `segmentation.R` sourced verbatim, plots stubbed, R 4.3.3). Inputs from a
  deterministic integer LCG reproduced bit-for-bit in the C# test (switched haplotype blocks, noise, clamping):

  | Track | SNPs / regime | Options | R segments (first–last pos, n, BAFseg) | Σ BAFphased / Σ BAFseg |
  |---|---|---|---|---|
  | t1 | 300, `runFastPcf` 0.15/0.15 | defaults | 2000–101000 100 0.49788217067718499; 102000–201000 100 0.71850197553634598; 202000–301000 100 0.58956708192825302 | 176.48445771217348 / 180.59512281417841 |
  | t2 | chr2 2000 (0.12/0.05) with 4 Mb gap + chr3 40 (< 50 ⇒ mean) | option 3 = option 1 | 6 segments (chr2 split at the gap; chr3 0.50574477136135099) | 1246.572572066784 / 1248.5744446992874 |
  | t2 | same | option 2 (mean) | 0.54146806168556205, 0.80023094666004202, 0.50063822877407005, 0.50009172968069704, 0.67059457948207901, 0.505747132062912 | Σ BAFseg 1246.572572066784 |
  | t3 | 16 000, `runPcfSubset` | defaults | 5 segments 3000/4000/2000/3500/3500 SNPs, BAFseg 0.50033228099346105 … 0.50094892680644998 | 10639.269587749242 / 10632.704436182976 |
  | t4 | 400, prior breakpoints 150500, 320000 | gamma 5, kmin 5 | 149/101/69/81 SNPs (0.59599266052246103, 0.60799378395080605, 0.80668817043304497, 0.80411475181579595) | 270.88505104064944 / 271.00505725383766 |
  | t4 | same | no_segmentation | one segment 0.63947307586669999 | Σ BAFseg 255.78923034668 |
  | t5 | 1500, sd 0.1, 7-SNP blocks | defaults | 336/367/396/401 SNPs (breakpoints at 100801/210901/329701 ≠ truth) | 869.7522670245171 / 876.50769936084771 |

  C#: every segment extent identical, BAFseg ≤ 1e−12, sums ≤ 1e−12 relative (R `cumsum`/`sum`/`mean` accumulate in
  80-bit long double; only Σ over 16 000 SNPs shows the 1e−13-relative difference).
- End-to-end (`set.seed(s); determine_copynumber(BAFvals, LogRvals, ρ, ψ_all, 1, ctrans, ctrans, 0.01, 0.05, 1000, 1000)`
  on the R `segment.baf.phased` output; logR = lv[seg] + ((i mod 7) − 3)·0.01 per chromosome-local SNP index, plus one `Inf`
  probe): e1 (t1, ρ 0.7, ψ 2.6): clonal 1:1 (pval 1, LogR −0.0003); (2,0)@0.47709449087559203 + (2,1) (pval 6.2072657065460898e−36,
  LogR 0.1501); (1,1)@0.37649750519839498 + (2,1) (pval 5.0025818343502595e−13). e2 (t2, ρ 0.85, ψ 3): 6 segments, e.g.
  (3,0)@0.39163557478371203 + (3,1) pval 2.2067308446755200e−101; clonal 1:1, 2:2, 1:1 (chr3). C# `BuildBattenbergSegments` +
  `FitSubclonalCopyNumberWithSnpTest`: 9/9 segments identical (startpos/endpos, LogR ≤ 1e−12, pval ≤ 1e−9 rel, states, fractions ≤ 1e−12).
- Harness note: R `read.table` turns numeric chromosome labels into integers, which `ctrans[...]` then indexes by
  position (NA startpos); the harness casts chromosomes to character as Battenberg's `read_bafsegmented` does.
- New shared helper: `StatisticsHelper.SampleQuantileType7` (R `quantile(type = 7)`, 1-based index arithmetic) —
  R 4.3.3 e.g. x = {0.31, 0.12, 0.97, 0.55, 0.12, 0.44, 0.08}: p 0.05 ⇒ 0.091999999999999998, p 0.88 ⇒ 0.66760000000000019.

## 2026-10 FIN-B24 F41 — Battenberg alternative solutions A–F, SDfrac and bootstrap CIs

- Sources opened: Wedge-lab/battenberg master 57a8f7e `R/fitcopynumber.R` (`callSubclones(..., noperms = 1000, seed =
  as.integer(Sys.time()))` → `set.seed(seed)` once, `determine_copynumber` called twice, before and after `merge_segments`;
  `determine_copynumber` sub-clonal branch: `all.edges = orderEdges(...)`, NA rows (`which(nMaj1<0|nMin1<0|...)`) moved
  last, `tau`, `sdl = sd(BAFke, na.rm = T)/sqrt(sum(!is.na(BAFke)))`, `sdtau = |τ(l+sdl) − τ|/2 + |τ(l−sdl) − τ|/2`;
  "Bootstrapping to obtain 95% confidence intervals": for each of the 6 options, `noperms` × `permBAFs = sample(BAFke,
  length(BAFke), replace = T)`, `permFraction[j] = τ(mean(permBAFs))`, `sdtaubootstrap = sd(permFraction)`, `tau25 =
  sort(permFraction)[25]`, `tau975 = sort(permFraction)[975]`; columns `nMaj1_X … frac1_X_0.975` for X = A–F, NA for
  clonal rows), `R/orderEdges.R` (six options per case). The bounds are fixed order statistics, **not** `quantile()`
  (type 7 is not used): exact 2.5 %/97.5 % points only for noperms = 1000; NA when fewer than 25 / 975 non-NA resamples.
  Alternative solutions B–F are the same `determine_copynumber` output row ⇒ in scope, ported.
- R RNG (R 4.3.3 `src/main/RNG.c`, `src/main/random.c`): `set.seed` → `RNG_Init` (50 × `seed = 69069·seed + 1`, then
  625 LCG values into `dummy[]`, `FixupSeeds` sets `mti = 624`), `MT_genrand` (MT19937 tempering, ·2.3283064365386963e−10,
  `fixup` into (0, 1)), `sample.kind = "Rejection"` (R ≥ 3.6): `R_unif_index(n)` = `rbits(ceil(log2 n))` (16-bit chunks of
  `floor(unif_rand()·65536)`) until < n. Ported as a private class; `set.seed(42); sample(5, 10, TRUE)` = 1 5 1 1 2 4 2 2 1 4.
- R cross-check (executed; `determine_copynumber` + `orderEdges` sourced verbatim, `set.seed(seed)` then one call):

  | Case | Input | Key R values (solution A unless noted) |
  |---|---|---|
  | s1 (ρ 0.8, ψ 2.5, seed 7) | 5 SNPs around (2,0)/(2,1) | τ 0.15605024544591201, SDfrac 0.0064798217571942303, SDfrac_BS 0.0058290159705588597, CI [0.145295959445356, 0.16716042067841799]; F = NA (x − 1 < 0) |
  | s1n500 / s1n20 | noperms 500 / 20 | SDfrac_BS 0.0060085874798433301 / 0.0084178392254478495; frac1_0.975 NA / both NA |
  | s2 (seed 3) | 1 SNP | SDfrac NA, SDfrac_BS 0, bounds 0.15626393223361601 (= τ(SNP BAF) ≠ τ(l)) |
  | s3 (ρ 0.55, ψ 3.1, seed −12345) | 3 SNPs | SDfrac 0.0335351207441906, SDfrac_BS 0.027099948479094901, CI [0.78233132237708902, 0.89740146907827201] |
  | s4 (ρ 0.35, ψ 1.9, seed 2³¹ − 1) | 5 SNPs | SDfrac 0.14517451315980001, SDfrac_BS 0.12696742413780401, CI [0.28061684037572598, 0.75502869611456302] |
  | s5 (ρ 1, ψ 2, seed 0) | 12 SNPs | SDfrac 0.0187275272879537, SDfrac_BS 0.018022906735212799, CI [0.66620363939436, 0.73567095563232898] |
  | e1 / e2 (F40 tracks, seeds 4711 / 99) | 3 / 6 segments, RNG stream across sub-clonal segments | all 60 columns of every row |

  C# `FitSubclonalCopyNumberWithBootstrap`: every column of every row (all six solutions, incl. NA placement and NA
  bounds) ≤ 1e−12 — the bootstrap columns too, i.e. the resamples are R's (R's long-double `mean`/`sd` change only
  the last bits).
- Statistical agreement (seed-independent): R, s4, seeds 1..400: mean (sd) SDfrac_A_BS 0.130353 (0.002746),
  frac1_A_0.025 0.280433 (0.007417), frac1_A_0.975 0.785209 (0.012439); C# seeds 100001..100100 agree within
  5·sd·√(1/100 + 1/400).
- Gap shown (pre-F41): Battenberg's `SDfrac_*`, `SDfrac_*_BS`, `frac1_*_0.025/0.975` and solutions B–F had no counterpart.

## 2026-10 FIN-B24 F53 — Multi-sample `ascat.asmultipcf`

- Source opened: `raw.githubusercontent.com/VanLoo-lab/ascat/master/ASCAT/R/ascat.asmultipcf.R` (master 61ddf3b, read in
  full): `ascat.asmultipcf(ASCATobj, ascat.gg, penalty = 70, out.dir, wsample = NULL, selectAlg = "exact", refine = TRUE,
  seed)` — one germline (`gg[, 1]`), `segmentlengths = unique(c(penalty, 25, 50, 100, 200, 400, 800))` (≥ penalty), per
  `chr` part: `bafna[homo, ] <- NA`; `useLogRonlySites = TRUE` ⇒ `Select_sites` = probes with any logR/BAF (every probe
  when logR is complete), `Select_sites2` = probes with a BAF; `mirrorBafMatrix`, `madWinsMatrixWithNA(·, 2.5, 25)`,
  `bafwins = mirror(madWins(mirror(baf)))`; `logRaveraged` over `startindices:endindices`; `nrow < 6` ⇒ column means;
  else `lrANDbaf = cbind(logRaveraged, bafwins)`, `w` = 0 at NA × `wsample` (length S ⇒ `c(wsample, wsample)`, length 2S,
  else `stop`), `ASmultiPCFcompact(nr = w, wSum = t(lrANDbaf)·w, gamma = segmentlength)` (exact) or `runFastASMultiPCF`
  (subsize 5000, step 4000); `refine`: per sample `compactASMulti` on the joint breakpoints + `ASmultiPCFcompact(gamma =
  segmentlength / S)`, taken when `nIntervals` drops; BAF correction (`rle`, `getMadwithNA`, `sqrt(sd² + μ²) < 2·sd`);
  NA/0 logR filled from the closest probe; genome-wide "adapt levels again"; ladder while any sample has ≥ 800 levels.
  No homozygous-stretch resegmentation ("included in the segmentation from the start"). Helpers `madWins`,
  `medianFilter` from `ascat.aspcf.R`.
- "WGD refit search" (algorithm doc §5.3, formerly "not implemented"): `grep -rniE "wgd|whole.?genome.?doubl|refit"
  ASCAT/R` on the same checkout finds only the `ascat.metrics` WGD *status* (`WGD` 0/1/"1+", GI score — ported in F30);
  ASCAT has no whole-genome-doubling refit search, so the phrase had no reference counterpart and was deleted.
- R quirks kept verbatim: (a) `ASmultiPCFcompact` `bestCost[1] = helper %*% (−Sum[, 1]·bestAver[, 1])` is NaN when the
  first column has a zero weight (homozygous first probe ⇒ BAF 0/0); `which.min` skips NaN, so no split after column 1 —
  with `refine`, a part without heterozygous probes therefore always merges its joint breakpoint away (cohort C1 chr4:
  S1's −1.5 deletion is found jointly, `refine = FALSE` keeps 150 + 150, `refine = TRUE` returns one level
  −0.7425429726108953); (b) `runFastASMultiPCF` passes the **unweighted** `t(x[…])` as `wSum` in the windows and marks
  segment *starts* (`mark[start0 + res$sta − 1]`) as compact-block ends, so its breakpoint can sit one probe right of the
  exact one (C4: 2501 vs 2500); (c) `wsample` is never used inside `ASmultiPCFcompact`; (d) the last averaging window is
  `lr[n:floor(n + 0.01 − 0.01)]`, which R evaluates as the descending range `n:(n − 1)` for n = 2, 32, 128, 16384, 65536
  (checked: `which(floor(n + 0.01 − 0.01) != n)` over 1..10⁵).
- R fails (verified): one sample, and any single-probe `chr` part — `Tumor_LogR[chr[[k]], ]` drops to a vector and
  `bafna[homo | is.na(homo), ] <- NA` stops with "incorrect number of subscripts on matrix". The port rejects both
  (`ArgumentException`; one sample ⇒ use `SegmentAlleleSpecificAspcf`).
- R cross-check (executed): R 4.3.3, `ascat.aspcf.R` + `ascat.asmultipcf.R` sourced (`print` stubbed), `chr = ch` one part
  per chromosome, `ascat.gg = list(germlinegenotypes = matrix(!het))`, `out.dir = NA`. Cohorts from the MINSTD generator
  of § F36 extended to S samples (per probe: 1 genotype draw shared, then 9 draws per sample — logR noise ×4, BAF
  orientation, BAF noise ×4; bit-identical in C#, `OncologyAnalyzer_AscatAsMultiPcf_Tests.Simulate`):

  | Cohort | S, seed, pHom | Probes (het) | Features |
  |---|---|---|---|
  | C1 | 2, 101, 0.3 | 1632 (858) | S2-only break at 400; forced 100-probe homozygous stretch; S1-only break at 900; 128-probe chr2; 4-probe chr3; homozygous-only chr4 with S1-only deletion |
  | C2 | 3, 202, 0.25 | 1437 (1016) | S1+S2 break at 500 (not S3); shared break at 800; 32-probe chr2; 5-probe chr3; S3-only break inside a 60-probe homozygous stretch |
  | C3 | 4, 303, 0.35 | 900 (522) | sample-specific subsets of breaks 200/400; S4-only deletion confined to an 80-probe homozygous stretch |
  | C4 | 2, 404, 0.3 | 6500 (4574) | one 6500-probe chromosome (fast windows) |
  | C5 | 2, 505, 0.3 | 800 (568) | S2-only 0.4 step below the joint penalty |

  Runs (14): C1, C1 `refine = FALSE`, C1 `wsample = c(1, 0.5, 2, 1)`; C2; C3, C3 `refine = FALSE`, C3 `wsample = c(1, 2,
  0.5, 1)`, C3 penalty 25, C3 penalty 0.001 (≥ 800 levels ⇒ ladder climbs to 25, output = penalty 25); C5, C5 `wsample =
  c(1, 4)` (S2 step called: 1–300 / 301–600); C4 exact, C4 `selectAlg = "fast"`, C4 fast + `wsample = c(1, 3)`.
- Result: every run, every sample: per-probe `Tumor_LogR_segmented` C# vs R max |Δ| ≤ 7.8e-16, `Tumor_BAF_segmented`
  ≤ 1.2e-16, NA pattern identical, breakpoints identical (0 mismatches); segments per sample C1 5/5 (no refine 7/7), C2
  6/6/6, C3 3/4/3/4 (no refine 6/6/6/6), C5 2/2 (weighted 2/3), C4 3/2.
- F36 defect found by (d): the single-sample germline-aware `ascat.aspcf` port averaged the empty range `[n, n − 1)` (NaN)
  for the last heterozygous probe when it and its predecessor are the last two probes of a part of 2/32/128/… probes. R
  (`ascat.aspcf`, all-heterozygous genome seed 7, chr1 300 + chr2 32 probes): segments 200 / 100 / 16 / 16, logR
  −0.22801451477129689 / 0.38581407161327735 / 0.32626688956574856 / −0.51100735978619538, BAF 0.5 / 0.70045650485928002
  / 0.61705520200405972 / 0.814577442655958; old C# split chr2 into 32 single-probe segments with BAF 0.7158163223300089.
  Fixed by the shared `RColonMean` (R `mean(x[a:b])`, descending when b < a).

---

## 2026-10 FIN-B24 F58 — Male X non-PAR germline re-genotyping (`ascat.aspcf` / `ascat.asmultipcf`, seeded)

- Sources opened: `ASCAT/R/ascat.aspcf.R` (master 61ddf3b) — `ascat.aspcf(…, seed = as.integer(Sys.time()))`:
  `set.seed(seed)`; `gg` from `ascat.gg` or `Germline_BAF < 0.3 | > 0.7`; `ghs = predictGermlineHomozygousStretches(chr,
  gg)` is computed **before** the per-sample block `if (!is.null(X_nonPAR) && gender[sample] == "XY")`:
  `nonPAR_index = which(SNPpos[, 1] %in% c("chrX", "X") & pos in [X_nonPAR[1], X_nonPAR[2]] & !is.na(gg[, sample]))`,
  `autosomes_info = table(gg[autosomes, sample])` (autosomes = `setdiff(chrs, sexchromosomes)`); if `length > 5`: all
  `TRUE` (homozygous), then with `Germline_BAF`: `DIST = 1 − max(x, 1 − x)`,
  `gg[nonPAR_index[which(rank(DIST, ties.method = "random") <= round(length(DIST) · (autosomes_info["FALSE"] /
  sum(autosomes_info))))]] = FALSE`; without it `gg[sample(nonPAR_index, round(…))] = FALSE`.
  `ASCAT/R/ascat.asmultipcf.R` ll. 30–63: identical block on the first germline column (`gender[1]`,
  `Germline_BAF[nonPAR_index, 1]`), no homozygous stretches. R 4.3.3 `base::rank`: `random =
  sort.list(order(x, stats::runif(sum(!nas))))` (printed from the R session; `set.seed(1); rank(c(3,1,3,2,3), "random")`
  = `sort.list(order(x, runif(5)))` = 4 1 5 2 3). `src/main/random.c` (R-4-3-branch, `do_sample`): without replacement
  `x[i] = i; for i < k: j = R_unif_index(n); y[i] = x[j] + 1; x[j] = x[--n]` (k < 2: `R_unif_index(n) + 1`, same first
  draw). R `round` is half-even (`nearbyint`). Reused: the F41 `RMersenneTwister` (set.seed scrambling, MT19937,
  `unif_rand` fixup, `R_unif_index` rejection sampling) + new `SampleWithoutReplacement`.
- Semantics noted: the probes ASCAT marks **heterozygous** are those *closest* to germline BAF 0/1 (smallest DIST), with
  exact 0/1 BAFs tied at DIST = 0 and ordered by the `runif` draw; the stretches use the caller's genotypes. ASCAT's default
  seed is the wall-clock second, so R itself is reproducible only with an explicit seed.
- R cross-check (executed; R 4.3.3, sourced ascat.aspcf.R + ascat.asmultipcf.R; harness: MINSTD tracks — per probe 3
  germline draws + 9 per sample — autosomes 1/2 plus X PAR1 (30 probes, 100 000–2 420 000), non-PAR (300 + 100 probes,
  3 000 000–152 700 000, logR step) and PAR2 (20 probes from 155 000 000), `X_nonPAR = c(2699521, 154931043)` (hg19),
  `gender = "XY"`, `Germline_BAF` (ascat.gg NULL) or `ascat.gg` (Germline_BAF NULL); `asmultipcf` with `penalty = 10`):

  | Track | Call | Branch | seed | m | h_auto | k (selected) | exact 0/1 ties | segments |
  |-------|------|--------|------|---|--------|--------------|----------------|----------|
  | M1 | aspcf | Germline_BAF / rank | 1 | 400 | 0.7033 | 281 | 159 | 9 |
  | M2 | aspcf | Germline_BAF / rank | 7777 | 400 | 0.3672 | 147 | 167 (tie-break decides) | 7 |
  | M3 | aspcf | ascat.gg / sample | 42 | 400 | 0.594 | 238 | — | 7 |
  | M4 | aspcf | Germline_BAF / rank | 1760090000 (time-type, as `as.integer(Sys.time())`) | 400 | 0.692 | 277 | 153 | 7 |
  | K1 | asmultipcf, 2 samples | rank | 1 | 400 | 0.7067 | 283 | 158 | 7 / 5 |
  | K2 | asmultipcf, 3 samples | rank | 2024 | 400 | 0.3073 | 123 | 146 (tie-break decides) | 5 / 6 / 4 |
  | K3 | asmultipcf, 2 samples | sample | 99 | 400 | 0.718 | 287 | — | 6 / 5 |

  C# (`SegmentAlleleSpecificAspcf(loci, het, AscatMaleXGenotyping)` / `AsMultiPcfOptions.MaleXGenotyping`): the
  re-genotyped probe set is identical to R's (`rownames(Tumor_BAF_segmented)`) on all 7 tracks; every segment extent,
  probe and heterozygous count identical; segment logR and BAF ≤ 1e−12 (observed max 5.6e−16 — R's long-double `mean`).
  A female model, a male without `X_nonPAR`, or ≤ 5 non-PAR probes leave the output identical to the default overload.
  Harness and generated data: scratchpad `wp28/harness.R`, `run.R`, `gen.py` (same MINSTD generator as `Simulate()` in
  `OncologyAnalyzer_AscatMaleXNonPar_Tests.cs`).

## 2026-10 FIN-B24 F59 — Missing data: the NA path of `ascat.asmultipcf` / `ascat.aspcf` (+ runASCAT)

- Sources opened (VanLoo-lab/ascat master 61ddf3b): `ASCAT/R/ascat.asmultipcf.R` ll. 79–129 (`bafna[homo | is.na(homo), ]
  <- NA`; `useLogRonlySites = TRUE`: `Select_sites = !(all-NA BAF row & all-NA logR row)`, `Select_sites2 = any BAF & any
  logR`, subset to the selected sites; `madWinsMatrixWithNA` — per column on `!is.na(x)`, **across the samples when
  `nrow(x) == 1`**; `logRaveraged[i, ] = mean(lrwins[start:end], na.rm = TRUE)` around each selected site), ll. 131–166
  (`length(indices) == 0` ⇒ `mean(lr, na.rm = TRUE)`; `< 6` sites ⇒ per-sample `mean(…, na.rm = TRUE)`; else `w = 0` at
  `is.na(t(lrANDbaf))`, `lrANDbaf[is.na] = 0`), ll. 206–253 (`bafASPCF[Select_sites2, ]`; `getMadwithNA` drops 0 and NA;
  `rle` per sample; `if (!all(is.na(yi)))`), ll. 274–291 (NA / 0 levels ← closest non-NA non-zero level of the part, `which.min`),
  ll. 304–324 (genome-wide `mean(Tumor_LogR[run], na.rm = TRUE)`, NaN ⇒ `prevlevel`, starting at 0), l. 342 (`Tumor_BAF_segmented`
  keeps the non-NA rows only). `ASCAT/R/ascat.aspcf.R` ll. 81–96 (`lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)])`,
  `Select_het = !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)`), ll. 131–158 (gap breakpoint `sum(…, na.rm = TRUE)`,
  `mean(lr[run], na.rm = TRUE)`), ll. 160–163 (no het ⇒ `mean(lr, na.rm = TRUE)`), ll. 190–203 (`pcfed[!is.na(towins)] =
  exactPcf(madWins(towins[!is.na(towins)]))`, 0 elsewhere; replace only if `!anyNA(dif)`), ll. 206–229 (`fillNA(zeroIsNA =
  TRUE)`, `prevlevel`). So the single-sample `ascat.aspcf` **also** accepts NA by default (its own path, not the
  multi-sample weights) — ported too (decision: R's default accepts NA, so the `ascat.aspcf` driver follows it; the
  all-heterozygous `SegmentAlleleSpecificAspcf(loci, penalty)` / `SegmentAlleleSpecific` fastAspcf path is not the
  `ascat.aspcf` driver — no `Select_het`, no NA handling in `fastAspcf` itself — and keeps rejecting NaN).
  `ASCAT/R/ascat.runAscat.R` ll. 227–231 (`r = lrrsegmented[names(bafsegmented)]`, `SNPposhet = SNPpos[names(bafsegmented), ]`),
  ll. 545–567 (`bafke = bafsegmented[bafpos][1]`, first probe of the segment WITH a segmented BAF), ll. 658–669 + l. 98
  (`n1all/n2all` NA at `is.na(lrr)`; `ifelse(baf[heteroprobes] <= 0.5, …)` is NA where the raw BAF is NA;
  `ploidy = mean(nA + nB, na.rm = TRUE)`).
- R behaviour checked in the session: `madWins(numeric(0), 2.5, 25)$ywin` = `numeric(0)`; `madWinsMatrixWithNA` of an all-NA
  column = all NA; `getMadwithNA(c(NA, NA))` = NA; `exactPcf(numeric(0), 6, 10)` = `numeric(0)`.
- R cross-check (executed; R 4.3.3, sourced ascat.aspcf.R + ascat.asmultipcf.R + ascat.runAscat.R; harness: the MINSTD
  `simulateMulti` of F53, then NA written by deterministic rules — `mod`: 1-based probe i with i·a ≡ 0 (mod b); `range`;
  `chr`; `hetevery`: every a-th heterozygous probe), `ascat.asmultipcf(obj, ascat.gg, penalty = 70)` and
  `ascat.aspcf(obj, ascat.gg, penalty = 70)` on the same object, `ascat.runAscat(gamma = 1)` on each:

  | Cohort | S | probes | NA rules (sample: track) | NA logR / NA het BAF per sample |
  |--------|---|--------|--------------------------|----------------------------------|
  | N1 (seed 611) | 2 | 1032 | S1 logR i·37 ≡ 0 (11); S1 BAF i·29 ≡ 0 (23); S2 BAF i·53 ≡ 0 (17) — scattered | 93 / 27; 0 / 33 |
  | N2 (seed 622) | 3 | 1505 | S2 logR + BAF all of chr2; S1 logR all of chr3 (5 probes); S3 logR 1340–1350 (in a forced homozygous stretch) — whole-chromosome runs | 5 / 0; 200 / 159; 11 / 0 |
  | N3 (seed 633) | 2 | 1021 | S1 BAF 301–450; S2 BAF every 3rd het probe; both: all of chr5, chr6 probes 1011–1012 (one selected site), logR 860–865; S1 logR all of chr7 — NA BAF at het probes | 26 / 114; 18 / 245 |

  Result: for all 6 segmentations (14 sample tracks) the per-probe `Tumor_LogR_segmented` and mirrored
  `Tumor_BAF_segmented` agree within 1e-12 (asserted per probe), the breakpoints (runs of equal level) are identical and
  the NA positions are identical (e.g. N1 asmultipcf: only probe 782 — BAF missing in both samples — has no segmented
  BAF, while ascat.aspcf drops 83 / 33 heterozygous probes; N2 asmultipcf S2 takes chr1's last level on chr2 by
  `prevlevel`). runASCAT on the 14 tracks: N1 S2 has no solution in R (both segmenters; not locked); the other 12 —
  purity, ψ, ploidy (1e-12), goodnessOfFit (1e-9) and seg_raw — are identical, e.g. N1 asmultipcf S1 ρ = 0.84,
  ψ = 2.35, ploidy 2.0164113785557989, GoF 98.219103302540461; N3 asmultipcf S1 ρ = 0.24, ψ = 2.15, ploidy
  4.7448648648648648. Before the runASCAT fixes the ploidy differed (e.g. N2 asmultipcf S2 2.2658 vs R 2.4598: probes with
  NA logR, and asmultipcf probes with a segmented BAF but a missing raw BAF, are NA in `n1all + n2all`).
- Harness: scratchpad `wp29/harness.R`, `cohorts.R`, `run.R` (segmentations → `expected.cs`), `fit.R` (runASCAT → `fits.cs`).

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

- **2026-10-10**: FIN-B24 F59 — NA path of `ascat.asmultipcf` / `ascat.aspcf` and runASCAT on NA segmentations R cross-check (3 NA cohorts, 14 tracks, 12 fits).
- **2026-10-10**: FIN-B24 F58 — male X non-PAR seeded germline re-genotyping (`AscatMaleXGenotyping`) in aspcf / asmultipcf R cross-check (7 tracks); the F38 "not ported" note is superseded.
- **2026-10-10**: FIN-B24 F53 — multi-sample `ascat.asmultipcf` port (`SegmentAlleleSpecificAsMultiPcf`) R cross-check (5 cohorts, 14 runs); "WGD refit search" phrase removed (no ASCAT counterpart); F36 last-window descending-range fix.
- **2026-10-10**: FIN-B24 F41 — Battenberg solutions A–F, SDfrac and seeded bootstrap CIs (`FitSubclonalCopyNumberWithBootstrap`, R RNG port) R cross-check + Monte-Carlo agreement section added.
- **2026-10-10**: FIN-B24 F40 — Battenberg `segment.baf.phased` port (`SegmentPhasedBaf`, `BuildBattenbergSegments`) R cross-check; built-in IMPUTE2/Beagle5 haplotype imputation BLOCKED (proof recorded).
- **2026-10-10**: FIN-B24 F38 — runASCAT sex-chromosome model (male haploid X/Y, `X_nonPAR` rule, XX default) R cross-check section added.
- **2026-10-10**: FIN-B24 F37 — runASCAT with homozygous segments (`bafke` NA ⇒ 0, all-probe ploidy) R cross-check section added.
- **2026-10-10**: FIN-B24 F39 — Battenberg per-SNP t-test (`FitSubclonalCopyNumberWithSnpTest`) and StatisticsHelper Student-t R cross-check section added.
- **2026-10-10**: FIN-B24 F36 — germline-aware `ascat.aspcf` (homozygous probes, homozygous-stretch resegmentation) R cross-check section added.
- **2026-10-09**: FIN-B24 F35 — the unsourced greedy `SegmentAlleleSpecific` heuristic was removed; the public name now delegates to `SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty = 70)` (ascat.aspcf port, R-verified 60/60), with the ASPCF finite-logR / BAF ∈ [0, 1] validation; legacy thresholds ignored (range-checked for compatibility). Test point 14 added; point 9 reworded.

- **2026-06-23**: Initial documentation.
- **2026-06-23**: Added ASPCF penalised-least-squares segmentation (Nilsen 2012, Ross 2021) and sub-clonal copy-number two-state mixture (Nik-Zainal 2012 / Battenberg) evidence for the residual-closing fix.
- **2026-10-09**: FIN-B24 F26 — `DeriveMultiplicity` ties half-to-even per facets-suite `expected_mutant_copies`; F27 `SubclonalIntegerTolerance` `[Obsolete]`.
- **2026-09-28**: B24 review — FitPurityPloidy = runASCAT port, ASPCF = ascat.aspcf port, sub-clonal fit = Battenberg determine_copynumber port; corrected corner case 2, integer-assignment and Battenberg decomposition statements; R cross-check section added.
