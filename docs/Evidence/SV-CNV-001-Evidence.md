# Evidence Artifact: SV-CNV-001

**Test Unit ID:** SV-CNV-001
**Algorithm:** Read-Depth Copy Number Variation Detection (windowed read depth → log2 ratio → integer copy number)
**Date Collected:** 2026-06-13

---

## Online Sources

### Yoon S, Xuan Z, Makarov V, Ye K, Sebat J (2009) — "Sensitive and accurate detection of copy number variants using read depth of coverage", Genome Research 19(9):1586–1592

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2752127/ (open-access PMC mirror of DOI:10.1101/gr.092981.109)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper)
**Retrieved how:** WebSearch `Yoon 2009 read depth copy number variants Genome Research log2 ratio formula normalization`; the publisher page genome.cshlp.org/content/19/9/1586 redirected to an authentication wall, so the open-access PMC copy (PMCID PMC2752127) was located by WebSearch `Yoon 2009 read depth copy number PMC PMCID full text event-wise testing GC correction` and opened with WebFetch.

**Key Extracted Points:**

1. **Read depth ∝ copy number (core hypothesis):** "We use the GC-adjusted RD within 100-bp windows as a quantitative measurement of genome copy number." The authors observed "a linear relationship between coverage and copy number" across regions with copy numbers of 1, 2, and 3. → mean read depth of a region is proportional to its copy number.
2. **Windowed read counting:** "RD was measured by counting the number of mapped reads in 100-bp windows, assigning each read only once by its start position." → read depth is summarised over non-overlapping windows of a fixed size.
3. **GC-content correction formula (verbatim):** "a simple adjustment was made according to the equation [r_i' = r_i × m / m_GC], where r_i are read counts of the i_th window, m_GC is the median read counts of all windows that have the same G+C percentage as the i_th window, and m is the overall median of all the windows." → per-window read counts are normalised to the overall median.

### CNVkit `call` command source (`cnvlib/call.py`, etal/cnvkit, master) — reference implementation

**URL:** https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation; primary paper Talevich et al. 2016 PLoS Comput Biol 12(4):e1004873, DOI:10.1371/journal.pcbi.1004873)
**Retrieved how:** WebSearch `CNVkit log2 absolute copy number formula 2^(log2+1) round ploidy threshold method calling`; the raw source file was fetched with WebFetch and the function bodies `absolute_threshold`, `_log2_ratio_to_absolute_pure`, `_log2_ratio_to_absolute` were read verbatim.

**Key Extracted Points:**

1. **log2 → absolute copy number ("round" method), verbatim:** `_log2_ratio_to_absolute_pure(log2_ratio, ref_copies)` returns `ncopies = ref_copies * 2**log2_ratio` with docstring `.. math :: n = r*2^v`. The integer copy number is then `absolutes.round().astype("int")`. → for a diploid reference (ref_copies = 2), copy number = round(2 · 2^log2).
2. **log2 ratio definition, verbatim (docstring of `_log2_ratio_to_absolute`):** `log2_ratio = log2(ncopies / ploidy)` ⇒ `2^log2_ratio = ncopies / ploidy` ⇒ `ncopies = ploidy * 2^log2_ratio`. → the log2 ratio is log2(observed copies / reference ploidy); equivalently log2(observed depth / reference depth).
3. **Reference log2 anchors, verbatim:** "a single-copy gain in a perfectly pure, homogeneous sample has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585, and a single-copy loss is log2(1/2) = -1.0." (CNVkit calling docs, https://cnvkit.readthedocs.io/en/stable/calling.html, accessed 2026-06-13 via WebFetch). → CN 3 ⇒ log2 ≈ +0.585; CN 1 ⇒ log2 = −1.0; CN 2 ⇒ log2 = 0.
4. **"threshold" method, verbatim defaults and mapping:** the docstring of `absolute_threshold` gives, for a 50%-clonality tumour heuristic, `R> log2(2:6 / 4)` = `-1.0  -0.4150375  0.0  0.3219281  0.5849625` and cutoffs `DEL(0) < -1.1`, `LOSS(1) < -0.25`, `GAIN(3) >= +0.2`, `AMP(4) >= +0.7`; the code iterates `for cnum, thresh in enumerate(thresholds): if row.log2 <= thresh: ... break`, assigning the index of the first threshold the log2 value falls at or below, and `np.ceil(_log2_ratio_to_absolute_pure(...))` above the last threshold. Default `thresholds = (-1.1, -0.25, 0.2, 0.7)`.
5. **Copy number is non-negative (verbatim):** the impure-sample path clamps `ncopies = max(0.0, ncopies)` so the call command "never emits a nonsensical negative copy number". → integer copy number ≥ 0.

### CNVkit "Calling copy number gains and losses" documentation — reference implementation docs

**URL:** https://cnvkit.readthedocs.io/en/stable/calling.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (canonical project documentation)
**Retrieved how:** WebFetch of the stable docs URL; the diploid conversion sentence and the single-copy gain/loss anchors were read verbatim.

**Key Extracted Points:**

1. **Diploid conversion (verbatim):** "the absolute copy number is calculated as 2 * 2^(log2 value). For example, a log2 value of 0.38 corresponds to 2^(0.38) = 1.3 times the reference ploidy, which for a diploid genome would be 2 * 2^(0.38) = 2.6." → CN = ploidy · 2^log2 with ploidy = 2.

---

## Documented Corner Cases and Failure Modes

### From Yoon et al. (2009)

1. **Zero / very low coverage windows:** RD is "counting the number of mapped reads" in a window; a window with zero mapped reads has RD = 0, for which log2(0/ref) is undefined (−∞). Such windows cannot be assigned a finite log2 ratio and are treated as an unobserved (no-call) window rather than producing −∞.

### From CNVkit `call.py`

1. **NaN log2:** `absolute_threshold` replaces a `nan` log2 with the neutral reference copy number ("log2=nan found; replacing with neutral copy number"). → a window with no usable signal defaults to copy-number-neutral (2) rather than a call.
2. **Negative extrapolation clamped:** copy number is clamped at 0 (`max(0.0, ncopies)`); copy number is physically ≥ 0.

---

## Test Datasets

### Dataset: Worked read-depth windows (derived from the cited formulas)

**Source:** Yoon et al. 2009 (RD ∝ CN; windowed counting) + CNVkit `_log2_ratio_to_absolute_pure` (CN = round(2·2^log2)) and the verbatim anchors log2(1/2)=−1.0, log2(3/2)=0.585, log2(2/2)=0.

Reference (diploid baseline) RD per window = 100 (= overall median of windows). Window size = 4 positions; depth listed per position.

| Window (per-position depth) | Window mean RD | log2(RD/100) | CN = round(2·2^log2) | Call |
|------------------------------|----------------|--------------|----------------------|------|
| 100,100,100,100 | 100 | log2(1.0) = 0.000 | 2 | Neutral |
| 50,50,50,50 | 50 | log2(0.5) = −1.000 | 1 | Deletion (loss) |
| 0,0,0,0 | 0 | log2(0) = −∞ (undefined) | — | No-call (homozygous-deletion candidate) |
| 150,150,150,150 | 150 | log2(1.5) = 0.585 | 3 | Duplication (gain) |
| 200,200,200,200 | 200 | log2(2.0) = 1.000 | 4 | Duplication (amplification) |

Worked check (CN formula): round(2·2^(−1.0)) = round(2·0.5) = round(1.0) = 1; round(2·2^0.585) = round(2·1.5) = round(3.0) = 3; round(2·2^1.0) = round(4.0) = 4; round(2·2^0) = round(2.0) = 2.

---

## Assumptions

1. **ASSUMPTION: Reference (diploid baseline) depth is the overall median of the windowed read depths.** Yoon et al. normalise to the overall median `m` of all windows (GC-correction equation `r_i' = r_i·m/m_GC`); CNVkit's log2 ratio is taken against a reference profile. When no external reference is supplied to `DetectCNV`, the per-sample reference RD is taken as the overall median of the non-zero window means, mirroring the Yoon overall-median baseline and the CNVkit "ratio against reference" definition. This is correctness-affecting (it sets the log2=0 anchor) but is the source-supported self-reference choice; an explicit baseline may be supplied to override it.
2. **ASSUMPTION: Diploid ploidy (2) is the copy-number baseline.** The CNVkit conversion `CN = ploidy·2^log2` and the anchors log2(1/2)=−1, log2(3/2)=0.585 are all stated for a diploid genome (ploidy 2). `DetectCNV` uses ploidy 2 as the default reference copy number; this is the standard human autosomal baseline used by both cited sources.

---

## Recommendations for Test Coverage

1. **MUST Test:** A window whose mean RD equals the reference RD has log2 ratio 0 and copy number 2 (Neutral). — Evidence: CNVkit `_log2_ratio_to_absolute_pure` (CN = 2·2^0 = 2); log2(2/2)=0.
2. **MUST Test:** A window with half the reference RD has log2 ratio −1.0 and copy number 1 (Deletion/loss). — Evidence: CNVkit anchor "single-copy loss is log2(1/2) = -1.0"; round(2·2^−1)=1.
3. **MUST Test:** A window with 1.5× the reference RD has log2 ratio 0.585 and copy number 3 (Duplication/gain). — Evidence: CNVkit anchor "single-copy gain ... log2(3/2) = 0.585"; round(2·2^0.585)=3.
4. **MUST Test:** A window with 2× the reference RD has log2 ratio 1.0 and copy number 4. — Evidence: CNVkit diploid conversion CN = 2·2^log2; round(2·2^1)=4.
5. **MUST Test:** Windowing — depth data is summarised into non-overlapping windows of the given size; mean RD per window is computed. — Evidence: Yoon et al. 2009 ("counting the number of mapped reads in 100-bp windows").
6. **MUST Test:** A zero-depth window produces a no-call (not a −∞ log2), because log2(0) is undefined. — Evidence: Yoon et al. RD = read count (0 ⇒ undefined ratio); CNVkit treats unusable signal as neutral/no-call.
7. **SHOULD Test:** Reference RD defaults to the overall median of window means when no baseline is given. — Rationale: Yoon overall-median baseline (ASSUMPTION 1).
8. **SHOULD Test:** Copy number is non-negative and the log2/CN relationship is monotonic non-decreasing in RD. — Evidence: CNVkit `max(0.0, ncopies)`; CN = 2·2^log2 is strictly increasing in log2.
9. **COULD Test:** Empty input yields empty output; null input throws ArgumentNullException. — Rationale: defined trivial / input-validation behaviour.

---

## References

1. Yoon S, Xuan Z, Makarov V, Ye K, Sebat J. 2009. Sensitive and accurate detection of copy number variants using read depth of coverage. Genome Research 19(9):1586–1592. https://doi.org/10.1101/gr.092981.109 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC2752127/)
2. Talevich E, Shain AH, Botton T, Bastian BC. 2016. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol 12(4):e1004873. https://doi.org/10.1371/journal.pcbi.1004873 (source: https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py; docs: https://cnvkit.readthedocs.io/en/stable/calling.html)

---

## Change History

- **2026-06-13**: Initial documentation.
