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
4. **Integer assignment:** `nA = pmax(round(nAfull),0)`, `nB = pmax(round(nBfull),0)` — round to nearest, clamp at 0; the major allele is the larger of the rounded {nA,nB}.

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
3. **Decomposition (algebra of the two-population model):** the observed real-valued allele-specific copy number is
   the fraction-weighted average of the two integer states, `n_obs = f·n₁ + (1 − f)·n₂`, with the two states the
   integers bracketing `n_obs` (`n₂ = ⌊n_obs⌋`, `n₁ = ⌈n_obs⌉`) and `f = (n_obs − n₂)/(n₁ − n₂) ∈ [0,1]`.
   Integer `n_obs` ⇒ clonal single state (`f ≈ 0` or `1`).

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
2. **Non-identifiability / multiple optima:** the sunrise plot can show several local minima (e.g. a 2n vs 4n solution); ASCAT selects the global minimum over the grid. A planted single-solution genome must have a unique minimum within tolerance.
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

**Source:** synthesised per Nilsen et al. 2012 PCF objective (deterministic).

| Parameter | Value |
|-----------|-------|
| Track | logR = 0.0 for loci 0–9, logR = 1.0 for loci 10–19 (one true breakpoint at index 10) |
| BAF | 0.5 throughout (balanced) |
| γ | 0.5 (small enough that 1 breakpoint beats 0: ΔSSE between merged and split = 25 ≫ γ) |
| Expected | 2 segments; breakpoint between position 9 and 10; means 0.0 and 1.0 |

### Dataset: Planted sub-clonal segment (mixture recovery)

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
9. **MUST Test (ASPCF):** On a constructed noisy track, the ASPCF penalised cost ≤ the greedy `SegmentAlleleSpecific` cost (DP global optimum) — Evidence: Nilsen 2012 (DP returns the global minimum).
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

## References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. (2010). Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R` (master). https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. (2016). Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. (2022). PICTograph: estimation of cancer cell fractions and clone trees. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ (2021). DeCiFering the elusive cancer cell fraction. Cell Systems / PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/
6. Nilsen G, Liestøl K, Van Loo P, et al. (2012). Copynumber: Efficient algorithms for single- and multi-track copy number segmentation. BMC Genomics 13:591. https://doi.org/10.1186/1471-2164-13-591 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3582591/)
7. Ross EM, Haase K, Van Loo P, Markowetz F (2021). Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37(13):1909–1911. https://doi.org/10.1093/bioinformatics/btaa538 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8317109/)
8. Nik-Zainal S, Van Loo P, Wedge DC, et al. (2012). The Life History of 21 Breast Cancers. Cell 149(5):994–1007. https://doi.org/10.1016/j.cell.2012.04.023 ; Battenberg, https://github.com/Wedge-lab/battenberg

---

## Change History

- **2026-06-23**: Initial documentation.
- **2026-06-23**: Added ASPCF penalised-least-squares segmentation (Nilsen 2012, Ross 2021) and sub-clonal copy-number two-state mixture (Nik-Zainal 2012 / Battenberg) evidence for the residual-closing fix.
