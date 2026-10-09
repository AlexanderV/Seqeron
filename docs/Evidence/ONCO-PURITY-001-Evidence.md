# Evidence Artifact: ONCO-PURITY-001

**Test Unit ID:** ONCO-PURITY-001
**Algorithm:** Tumor Purity Estimation (from somatic SNV VAF / allele-specific copy number)
**Date Collected:** 2026-06-14

---

## Online Sources

### CNAqc — Quality Control of allele-specific copy numbers, mutations and tumour purity (vignette)

**URL:** https://caravagnalab.github.io/CNAqc/articles/CNAqc.html
**Retrieved by:** WebFetch of the URL above (CNAqc package vignette), 2026-06-14, prompting for the expected-VAF formula and its variable definitions.
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation / official package documentation), backing the peer-reviewed Genome Biology 2024 paper (rank 1).

**Key Extracted Points:**

1. **Expected-VAF formula (verbatim):** `v_m(c) = mπc / [2(1-π) + π(n_A+n_B)]`.
2. **Variable definitions (verbatim phrasing):** `m` = "mutations present in m copies of the tumour genome" (multiplicity); `π` = "tumour purity"; `c` = "mutations present in a percentage 0<c<1 of tumour cells" (cancer cell fraction / clonality); a segment is written `n_A:n_B` for the allele-specific copy numbers `n_A` and `n_B`.
3. **Denominator interpretation (verbatim phrasing):** `2(1−π)` represents "a healthy diploid normal" contribution and `π(n_A+n_B)` is "the proportion of all reads from the tumour"; total denominator = mean number of allele copies per cell across the mixture.

### CNAqc — Computational validation of clonal and subclonal CNAs (Genome Biology 2024, paper)

**URL:** https://link.springer.com/article/10.1186/s13059-024-03170-5 (search-indexed full text; abstract/figure captions returned via web search snippet of this DOI)
**Retrieved by:** WebSearch queries `CNAqc tumor purity expected VAF formula …` and `CNAqc "expected VAF" diploid heterozygous "purity" peak 0.5 …`, 2026-06-14, surfacing the paper's worked numeric examples. (Full-text and bioRxiv PDF returned HTTP 403 to direct fetch; numeric examples taken from the indexed snippets of this DOI.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Genome Biology).

**Key Extracted Points:**

1. **General clonal/subclonal form (verbatim from snippet):** for mutations on (possibly two-state) CNAs, `v = (m₁ρ₁ + m₂ρ₂)π / { 2(1−π) + π[(n_{A,1}+n_{B,1})ρ₁ + (n_{A,2}+n_{B,2})ρ₂] }`. For a single clonal state this reduces to the vignette form with c = 1.
2. **Worked example — diploid heterozygous, purity 60%:** "real purity of 60%" for a heterozygous diploid mutation (m=1, n_A+n_B=2) corresponds to expected VAF 30%; the purity tolerance band 55–65% maps to "a VAF range of 27.50–32.5%". (Confirms v = π/2 ⇒ π = 2·VAF.)
3. **Worked example — 2:1 segment, purity = 1:** "For 2:1 segments, there are two peaks of clonal mutations at 33% and 66% VAF" — i.e. m=1 ⇒ 1/3 ≈ 0.333 and m=2 ⇒ 2/3 ≈ 0.667 with n_A+n_B = 3, π = 1.

### FACETS — allele-specific copy number and clonal heterogeneity (NAR 2016)

**URL:** https://academic.oup.com/nar/article/44/16/e131/2460163
**Retrieved by:** WebFetch of the URL above, 2026-06-14, prompting for the purity/copy-number mixing model and full citation.
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nucleic Acids Research).

**Key Extracted Points:**

1. **Mixing model (verbatim phrasing):** tumor parental copy numbers adjusted for cellular fraction Φ are `m* = mΦ + (1 − Φ)` and `p* = pΦ + (1 − Φ)`; the model mixes a "normal diploid (1,1)" genotype with an "aberrant (m,p)" genotype at mixing proportion Φ (cellular fraction). This independently confirms the `2(1−π) + π·n_tot` denominator structure used by CNAqc (normal contributes 2 copies weighted 1−π).

### ABSOLUTE — Absolute quantification of somatic DNA alterations (Carter et al. 2012)

**URL / retrieval:** Citation record retrieved via Europe PMC REST API
`https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:22544022&format=json&resultType=core` (WebFetch), 2026-06-14. (PubMed/Nature direct pages returned reCAPTCHA / 403; the EuropePMC core record supplied the verified bibliographic metadata.)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, Nature Biotechnology).

**Key Extracted Points:**

1. **Verified citation:** Carter SL, Cibulskis K, Helman E, McKenna A, Shen H, Zack T, Laird PW, Onofrio RC, Winckler W, Weir BA, Beroukhim R, Pellman D, Levine DA, Lander ES, Meyerson M, Getz G. "Absolute quantification of somatic DNA alterations in human cancer." *Nature Biotechnology* 30(5):413–421, 2012. DOI: 10.1038/nbt.2203.
2. **Context (from search overview):** ABSOLUTE converts allelic fractions of point mutations into per-cancer-cell allele counts (cellular multiplicity) "by correcting for sample purity and local copy-numbers" — the same purity/copy-number correction inverted here to estimate purity.

### CNAqc R source — `expected_vaf_fun` (reference implementation)

**URL:** https://github.com/caravagnalab/CNAqc — `R/equations.R` (shallow git clone of master, 2026-09-28, B24 review)
**Authority rank:** 3 (original tool code)

**Key Extracted Points:**

1. `expected_vaf_fun(m, M, mut.allele, p) = mut.allele·p / (2(1−p) + p(m+M))` — identical to the vignette formula (c = 1).
2. `expectations_generalised` enumerates multiplicities `1:m` and `1:M` only, i.e. 1 ≤ multiplicity ≤ Major ≤ n_tot.
3. Python port cross-check (IEEE double): forward VAFs 3:1 m=3 p=0.7 → 0.6176470588235293; 3:1 m=1 → 0.20588235294117646; 2:0 m=2 p=0.45 → 0.45; 1:0 p=0.8 → 0.6666666666666667; 2:2 m=2 p=0.35 → 0.25925925925925924; the inversion π = 2v/[m + v(2−n_tot)] recovers p to ≤ 2 ulp. At p = 1 the clonal peak v = m/n_tot inverts to 1 + k·ulp for 66 of the 210 (m ≤ n_tot ≤ 20) pairs (e.g. m=1, n_tot=5 → 1.0000000000000002); maximum excess over n_tot ≤ 2000 is 922 ulp ≈ 0.47·(n_tot+4)·ε.
4. FIN-B24 re-verification (2026-10-09): CNAqc `expected_vaf_fun` evaluated in R: (1,4,1,1) → 0.20000000000000001, (2,3,2,1) → 0.40000000000000002 (the doubles of literals 0.2 / 0.4); `sort(unique(c(1:m,1:M)))` for 1:4 → 1..4, for 2:3 → 1..3 (max multiplicity = Major). Tolerance (n_tot+4)·ε, ε = 2⁻⁵² (Python IEEE double, same operation order as C#): all 2 001 000 exact peaks v = m/n_tot, 1 ≤ m ≤ n_tot ≤ 2000, accepted (max excess 0.468·(n_tot+4)·ε) and clamp to 1.0; with v = (m/n_tot)(1+1e-9) 0 of 1 999 000 pairs (v ≤ 1) falsely accepted.

### R `stats::density.default`, `bw.nrd0`, `peakPick::peakpick` — CNAqc's peak detector (FIN-B24 F33)

**Sources (opened 2026-10-09):** R trunk (github.com/wch/r-source, sparse clone) `src/library/stats/R/density.R` (`density.default`, incl. `old.coords`), `src/library/stats/src/massdist.c` (`BinDist`), `src/library/stats/src/approx.c` (`approx1`), `src/main/seq.c` (`seq.int(from, to, length.out)`: `from + i·by`, last element = `to`), `src/nmath/fround.c` (R ≥ 4 `round(x, d)`), `src/library/graphics/R/hist.R` (break fuzz `1e-7·median(diff)`); R 4.3.3 `body(density.default)` / `bw.nrd0` (installed); cran/peakPick 0.11 `R/peakpicking.R` (`peakpick`, `smallpeaks`, `helperpeak`, `keepmax`); caravagnalab/CNAqc 1.1.5 (commit 4b7cea4a) `R/peak_algorithms.R` `smooth_data` = `density(VAF, kernel = 'gaussian', adjust = kernel_adjust, na.rm = T)` (defaults: `bw = "nrd0"`, `n = 512`, `cut = 3`) and `simple_peak_detector` = `peakPick::peakpick(cbind(x, y), neighlim = 1..5)[, 2]`.
**Authority rank:** 3 (reference implementation code).

**Key Extracted Points:**

1. `density.default`: `bw = adjust·bw.nrd0(x)`; `from = min − 3bw`, `to = max + 3bw`; mass 1/N linearly binned (`BinDist`) on `n = 512` points over `[from − 4bw, to + 4bw]`; kernel ordinates `seq.int(0, L, length.out = 2n)` mirrored, `L = (2n−1)/(n−1)·(up−lo)` (R ≥ 4.4, `old.coords = FALSE`) or `2·(up−lo)` (R ≤ 4.3 / `old.coords = TRUE`, PR#18337); FFT circular convolution, `pmax(0, ·)`, `approx` onto `seq(from, to, length.out = n)`.
2. `bw.nrd0(x) = 0.9·min(sd, IQR/1.34)·n^(−1/5)`, fall-backs sd → |x₁| → 1 (type-7 quantiles).
3. `peakpick(mat, neighlim, deriv.lim = 0.04, peak.min.sd = 0.5, peak.npos = 10)`: per column — centred derivative `diff(lag = 2)/2`; candidate = non-flat +→− sign change adjacent with `|der| < deriv.lim`; drop npos head/tail; keep iff `v > mean + 0.5·sd/√(2npos+1)` over ±npos (sd of one value is NA → kept); repeatedly drop the lower point of the closest pair while ≤ neighlim apart. CNAqc passes the 2-column (x, y) matrix and reads column 2 (y).
4. Port (`StatisticsHelper.BandwidthNrd0`, `GaussianKernelDensity`, `PeakPick`, Infrastructure, additive): the convolution is evaluated directly (binned mass occupies the first half of the zero-padded buffer, so circular = linear) instead of by FFT.
5. Reference numbers (R, `sprintf("%.17g")`, data = TestData/CNAqc seeded sets): `bw.nrd0` D1 1:1 (n = 400) 0.012856625638346407, D3 1:1 (n = 700) 0.018444382013876173, D1 2:1 (n = 250) 0.038758432750663083; `bw.nrd0(c(1,2,3,4,10))` 0.97358462285063574, `(5,5,5)` 3.6123370279210381, `(0,0)` 0.78349550696651171, `(1,1,1,1,9)` 2.3337454992375779. D1 1:1 density, new coords: y[1] 0.00095861337313726954, y[256] 7.8057017230652654 (x[256] 0.35777992026975591), y[512] 0.0017138350430965449; old coords (R 4.3.3 native): y[256] 7.8128770632778197. Full x/y at indices 1, 2, 100, 200, 256, 300, 400, 511, 512 for six (dataset, adjust, coords) cases are locked in `StatisticsHelper_GaussianKde_Tests`; C# agrees to ≤ 1e−12 relative (y) / 1e−14 (x, bw). `peakpick` on those densities (1-based, neighlim 1–5): D1 1:1 adj 1 → 243; adj 0.5 → 193,234,254,291,362,479 (old coords 253 for 254); D3 1:1 → 140,292; D1 2:1 → 181,368. Synthetic `v = 0.05 sin(i/4) + 0.03 sin(i/1.7) + 0.0004 i`, i = 0..199 (0-based): neighlim 0 → 25 peaks, 1/8 → 12,25,34,56,78,87,100,109,131,153,162,175,184, 12/20 → 12,34,56,78,109,131,153,184; npos = 3 & sd = 0, or npos = 0 → 19 peaks (4,12,…,195).

### CNAqc `analyze_peaks` — peak-based purity QC (FIN-B24 F34)

**Source:** caravagnalab/CNAqc 1.1.5, commit 4b7cea4a (git clone via the session proxy, 2026-10-09): `R/analyze_peaks.R` (defaults `karyotypes = c('1:0','1:1','2:0','2:1','2:2')`, `min_karyotype_size = 0`, `min_absolute_karyotype_mutations = 100`, `p_binsize_peaks = 0.005`, `purity_error = 0.05`, `VAF_tolerance = 0.015`, `n_bootstrap = 1`, `kernel_adjust = 1`, `matching_strategy = "closest"`, `min_VAF = 0`), `R/peak_algorithms.R` (`analyze_peaks_common`, `smooth_data`, `simple_peak_detector`, `mixture_peak_detector`, `combined_peak_detector`), `R/vaf_functions.R` (`expected_vaf_peak`, `ascat`), `R/equations.R` (`delta_vaf_karyo`, `purity_from_vaf`, `compute_delta_purity`), `R/peak_detector.R` (`overlap_bands`; legacy `peak_detector` = rightmost matching), `R/print.cnaqc.R`; caravagnalab/BMix `R/bmixfit.R`, `R/bmixfit_EM.R`.
**Authority rank:** 3 (reference implementation code).

**Key Extracted Points:**

1. `analyze_peaks` filters `VAF > min_VAF`, recomputes `n_karyotype` over all karyotypes and calls `analyze_peaks_common` **without** `matching_strategy` (and `p_binsize_peaks` is unused inside) — 1.1.5 always runs closest matching; "rightmost" survives only in the legacy `peak_detector`.
2. Analysed karyotypes: in `karyotypes` ∧ n ≥ `min_absolute_karyotype_mutations` ∧ n / Σ n_all ≥ `min_karyotype_size`. Expected peaks m ∈ `unique(c(1, Major))`: `m·π / (2(1−π) + π·ploidy)`; bands `delta_vaf = 2·m·ε / (2 + π(ploidy − 2))²`.
3. Data peaks = KDE peaks (`density` → `peakpick(neighlim = 1..5)` union, `round(x, 2)`, `round(y, 2)`, `distinct(x)`, [0, 1]; `counts_per_bin = hist(VAF, seq(0,1,0.01))$counts[round(100x)]`, first index 0 → 1; `discarded = y ≤ max(y)/20`) **plus** BMix peaks (`bmixfit(NV, DP, K.Binomials = 1:4, K.BetaBinomials = 0)`, ICL-selected; each `B.params` snapped to `den$x[which.min(|den$x − p|)]`, never discarded). BMix is stochastic (`kmeans(nstart = 100)`, `sample`, `runif` jitter) → not ported; its means are an input (`MixturePeaks`).
4. Matching: `which.min(|x − peak|)` over non-discarded peaks; `offset = compute_delta_purity(x, peak − x, ploidy, m) = 2·m·Δ / (m + x(2 − ploidy))²`; `weight = n_K / Σ n_analysed` (1.1.5 does not divide by the number of peaks — the legacy detector did); `matched = overlap_bands(x, VAF_tolerance, peak ± delta_vaf)`.
5. QC: per karyotype `arrange(desc(counts_per_bin))` first row's `matched` → PASS/FAIL; sample `group_by(QC) %>% summarise(prop = sum(weight)) %>% arrange(desc(prop))` first (tie → FAIL, alphabetical group order + stable arrange); `score = weight %*% offset`. **No |score| ≤ ε test**, no purity proposal — `print.cnaqc` prints "Purity correction: round(100·score)%".
6. Harness (scratch, not committed): CNAqc/peakPick/BMix R files sourced with `CNAqc:::`, `peakPick::`, `BMix::` prefixes stripped and `easypar::run` replaced by `lapply` + error filter; dplyr 1.2.1 (+ rlang, vctrs, tibble, pillar, cli, glue, … current cran mirrors) built from github.com/cran mirrors (CRAN host blocked); R ≥ 4.4 `density` = R-trunk `density.R` evaluated in the stats namespace. Rightmost reference = `analyze_peaks_common` with its two matching lines replaced by the legacy `peak_detector` matching block (verbatim).
7. Data generator (R, `set.seed(20261009/10/11)`): `DP = pmax(10, rpois(n, depth))`, `NV = rbinom(n, DP, m·π·ccf / (2(1−π) + π(M+m)))`, VAF = NV/DP; D1 π 0.7 depth 120 (1:1 ×400 m1; 2:1 ×250 m1 60 % / m2; 2:2 ×150 50/50; 1:0 ×120; 3:1 ×60 m1/m3); D3 π 0.45 depth 80 (1:1 ×500; 1:1 ×200 CCF 0.35; 2:0 ×150 m1 30 % / m2; 4:0 ×60 m4 relabelled 2:0); D4 π 0.3 depth 60 (1:1 ×150; 2:1 ×90; 1:0 ×110). Written with `%.17g` to `TestData/CNAqc`.
8. Reference results (λ, sample QC): D1 π 0.7 0.0023756289876209163 PASS (R 4.3.3 native density: same λ, y 8.64/7.95 vs 8.63/7.94); D1 π 0.5 −0.29645326597409133 FAIL; D1 adjust 0.5 0.023857243903529196 PASS; D1 ε 0.01 tol 0.005: 1:0 FAIL, sample PASS; D1 π 0.62 −0.11362481021020236 PASS (1:0 FAIL); D3 π 0.45 0.014615384615384629 PASS (min_VAF 0.12: 0.030000000000000027); D4 π 0.3 0.0039996443291489122 (size 0.35: 0; min 80: −0.010660215188956853; min 200: no analysis); BMix (set.seed 7) D1 0.0045358591466179345, D3 0.0037889790909820253, D3 rightmost −0.13824761798675933 (2:0 FAIL). Every match row (peak, δ, x, y, counts_per_bin, offset_VAF, offset, weight, matched, QC) and peak table is locked in `OncologyAnalyzer_AnalyzePurityPeaks_Tests`.

---

## Documented Corner Cases and Failure Modes

### From CNAqc

1. **Multiplicity ambiguity on amplified segments:** on a 2:1 (n_tot=3) segment clonal mutations form two VAF peaks (1/3 and 2/3) because multiplicity m may be 1 or 2; purity cannot be inferred from VAF alone without knowing m and the copy-number state. Copy-neutral diploid heterozygous (1:1) loci avoid this — m = 1, n_tot = 2 — which is why purity = 2·VAF is the robust closed form there.
2. **Subclonal mutations (c < 1):** for c < 1 the VAF is depressed (v ∝ c); treating a subclonal VAF as clonal underestimates purity. Purity must be estimated from clonal mutations.

### From the checklist by-area definition (ONCO-PURITY-001)

1. **Purity < 0.1 (below detection limit):** very low purity yields VAFs near sequencing noise.
2. **No heterozygous SNPs / no informative variants:** with no usable variant evidence purity is undefined.
3. **High stromal contamination:** equivalent to low purity.

---

## Test Datasets

### Dataset: CNAqc clonal heterozygous diploid worked example

**Source:** CNAqc (Genome Biology 2024), doi:10.1186/s13059-024-03170-5; CNAqc vignette.

| Parameter | Value |
|-----------|-------|
| Segment (n_A:n_B) | 1:1 (copy-neutral diploid, n_tot = 2) |
| Multiplicity m | 1 (heterozygous somatic SNV) |
| Clonality c | 1 (clonal) |
| Expected VAF at purity 0.60 | 0.30 (= π/2) |
| Purity from VAF 0.30 | 0.60 (= 2·VAF) |
| Purity tolerance band 0.55–0.65 | VAF band 0.275–0.325 |

### Dataset: CNAqc 2:1 amplified segment worked example (purity = 1)

**Source:** CNAqc (Genome Biology 2024), doi:10.1186/s13059-024-03170-5.

| Parameter | Value |
|-----------|-------|
| Segment (n_A:n_B) | 2:1 (n_tot = 3) |
| Purity π | 1.0 |
| Clonal VAF, m = 1 | 1/3 ≈ 0.3333… |
| Clonal VAF, m = 2 | 2/3 ≈ 0.6667… |

---

## Assumptions

1. **ASSUMPTION: VAF-only purity estimator uses the copy-neutral diploid heterozygous model.** `EstimatePurityFromVAF` assumes the supplied variants are clonal (c = 1) heterozygous (m = 1) somatic SNVs at copy-neutral diploid (n_tot = 2) loci, giving the closed form ρ = 2·VAF. This is the textbook special case and the band the CNAqc example uses; it is stated explicitly in the API contract. It is a modeling scope choice, not an invented numeric constant — the formula itself is fully source-derived.
2. **ASSUMPTION: aggregation across variants uses the median VAF.** When multiple clonal heterozygous SNVs are supplied, their per-variant purity estimates are combined by the median (robust to subclonal/outlier VAFs). The literature establishes the per-variant relation; the choice of a robust central estimator over the set is a documented, non-correctness-affecting aggregation policy (it does not change the single-variant formula). Recorded as an assumption for transparency.

---

## Recommendations for Test Coverage

1. **MUST Test:** `EstimatePurityFromVAF` on a clonal heterozygous diploid SNV with VAF 0.30 returns purity 0.60. — Evidence: CNAqc worked example (π/2 = 0.30 ⇔ ρ = 2·VAF).
2. **MUST Test:** `EstimatePurityFromVAF` boundary VAF 0.50 → purity 1.0; VAF 0.0 → purity 0.0. — Evidence: ρ = 2·VAF closed form, INV 0 ≤ ρ ≤ 1.
3. **MUST Test:** `EstimatePurity` (allele-specific) inverting v = mπ/[2(1−π)+π·n_tot] recovers π for a clonal SNV on a 2:1 (n_tot=3) segment: π=1, m=1, v=1/3 ⇒ π=1; m=2, v=2/3 ⇒ π=1. — Evidence: CNAqc 2:1 peaks 33%/66%.
4. **MUST Test:** `EstimatePurity` on diploid heterozygous segment (n_tot=2, m=1) with VAF 0.30 returns 0.60 (must agree with VAF-only estimator). — Evidence: CNAqc 60%/30% example.
5. **MUST Test:** invalid inputs (VAF outside [0,1], VAF > 0.5 for the diploid model implying purity > 1, empty variant list, non-positive copy number) throw / are rejected. — Evidence: invariant 0 ≤ ρ ≤ 1; formula domain.
6. **SHOULD Test:** median aggregation across several heterozygous SNVs with mixed VAFs returns the median-derived purity. — Rationale: robustness corner case.
7. **COULD Test:** purity below detection (VAF near 0) yields a small purity near 0 without error. — Rationale: low-purity edge case from checklist.

---

## References

1. Antonello A, Bergamin R, Calonaci N, Househam J, Milite S, Williams MJ, Anselmi F, d'Onofrio A, Sundaram V, Sosinsky A, Cross WCH, Caravagna G (2024). Computational validation of clonal and subclonal copy number alterations from bulk tumor sequencing using CNAqc. *Genome Biology* 25(1):38. https://doi.org/10.1186/s13059-024-03170-5
2. CNAqc package vignette (Caravagna lab). Quality control of allele-specific copy numbers, mutations and tumour purity. https://caravagnalab.github.io/CNAqc/articles/CNAqc.html (accessed 2026-06-14)
3. Carter SL, Cibulskis K, Helman E, et al. (2012). Absolute quantification of somatic DNA alterations in human cancer. *Nature Biotechnology* 30(5):413–421. https://doi.org/10.1038/nbt.2203
4. Shen R, Seshan VE (2016). FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. *Nucleic Acids Research* 44(16):e131. https://doi.org/10.1093/nar/gkw520

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: B24 review — CNAqc R source cross-check; boundary-rounding and multiplicity-bound defects recorded (cross-batch, B22 file).
- **2026-10-09**: FIN-B24 F24/F25 — m > n_tot rejected; π = 1 rounding tolerance (n_tot+4)·ε with clamp; R + Python re-verification recorded (point 4).
- **2026-10-09**: FIN-B24 F33 — R `density.default`/`bw.nrd0`/`peakPick::peakpick` (CNAqc peak detector) sources and reference numbers recorded.
- **2026-10-09**: FIN-B24 F34 — CNAqc `analyze_peaks` peak-based purity QC source review, R harness, seeded datasets and reference results recorded.
