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

---

## References

1. Van Loo P, Nordgard SH, Lingjærde OC, et al. (2010). Allele-specific copy number analysis of tumors. PNAS 107(39):16910–16915. https://doi.org/10.1073/pnas.1009843107
2. VanLoo-lab/ascat reference implementation, `ASCAT/R/ascat.runAscat.R` (master). https://github.com/VanLoo-lab/ascat
3. McGranahan N, Furness AJS, Rosenthal R, et al. (2016). Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. Science 351(6280):1463–1469. https://doi.org/10.1126/science.aaf1490
4. Zheng L, et al. (2022). PICTograph: estimation of cancer cell fractions and clone trees. Bioinformatics 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac440
5. Satas G, Zaccaria S, El-Kebir M, Raphael BJ (2021). DeCiFering the elusive cancer cell fraction. Cell Systems / PMC8542635. https://pmc.ncbi.nlm.nih.gov/articles/PMC8542635/

---

## Change History

- **2026-06-23**: Initial documentation.
