# Evidence Artifact: ONCO-CNA-001

**Test Unit ID:** ONCO-CNA-001
**Algorithm:** Copy-Number Alteration Classification (log2 copy ratio → absolute copy number → CNA state)
**Date Collected:** 2026-06-14

---

## Online Sources

### CNVkit reference implementation — `cnvlib/call.py` (`absolute_threshold`, `_log2_ratio_to_absolute_pure`)

**URL:** https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py
**Accessed:** 2026-06-14 (fetched with WebFetch; retrieved the function bodies verbatim)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **log2 → absolute copy number (pure sample):** function `_log2_ratio_to_absolute_pure(log2_ratio, ref_copies)` returns `ncopies = ref_copies * 2**log2_ratio`. For a diploid reference `ref_copies = ploidy = 2`, this is `n = 2 · 2^log2`.
2. **Hard-threshold integer calling (`absolute_threshold`):** docstring — "Call integer copy number using hard thresholds for each level. Integer values are assigned for log2 ratio values less than each given threshold value in sequence, counting up from zero. Above the last threshold value, integer copy numbers are called assuming full purity, diploidy, and rounding up."
3. **Binning loop (verbatim logic):** `cnum = 0; for cnum, thresh in enumerate(thresholds): if row.log2 <= thresh: break; else: cnum = int(np.ceil(_log2_ratio_to_absolute_pure(row.log2, ref_copies)))`. So the called copy number is the index of the FIRST threshold the log2 value is `<=` (boundary inclusive); if the log2 value exceeds every threshold, the copy number is `ceil(2 · 2^log2)`.
4. **Default thresholds (verbatim from `do_call`):** `thresholds = (-1.1, -0.25, 0.2, 0.7)`, mapping to copy states `[0, 1, 2, 3, 4+]`.
5. **Verbatim heuristic (docstring):** for single-copy gains/losses assuming 50% tumor cell clonality, `R> log2(2:6 / 4) = -1.0  -0.4150375  0.0  0.3219281  0.5849625`; "Allowing for random noise of +/- 0.1, the cutoffs are: `DEL(0) < -1.1`, `LOSS(1) < -0.25`, `GAIN(3) >= +0.2`, `AMP(4) >= +0.7`."
6. **NaN handling:** "log2=nan found; replacing with neutral copy number" — a NaN log2 ratio is treated as the neutral reference copy number (no call → diploid).

### CNVkit user documentation — `call` command threshold method

**URL:** https://cnvkit.readthedocs.io/en/stable/pipeline.html
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 3 (official project documentation for the reference implementation)

**Key Extracted Points:**

1. **Mapping table (verbatim):** "If log2 value is up to | Copy number" — `-1.1 → 0`, `-0.4 → 1`, `0.3 → 2`, `0.7 → 3`, `… → …`. (NOTE: the docs page shows the germline-tuned variant `-0.4 / 0.3`; the source code `call.py` default for tumor samples is `-0.25 / 0.2`. The implementation follows the source-code defaults, which the same docstring labels the tumor-sample heuristic; the docs `-0.4 / -0.3` are documented there as the germline-precision alternative.)
2. **Purity caveat (verbatim):** "The default threshold values are reasonably 'safe' for a tumor sample with purity of at least 30%."
3. **Threshold derivation:** Python example `np.log2((copy_nums + .5) / 2)` shows the cutoffs derive from log-transforming integer copy numbers (+0.5 for rounding) over ploidy 2.

### GISTIC2.0 — Mermel et al. (2011), Genome Biology

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
**Accessed:** 2026-06-14 (found via WebSearch; PMC article fetched with WebFetch)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Low-amplitude (noise) threshold (verbatim):** "a low amplitude threshold (log2 ratio of ± 0.1) that only eliminates low-level artifactual segments." This corroborates the CNVkit ±0.1 noise band around the neutral state.
2. **High-amplitude thresholds (verbatim):** "a high amplitude threshold (log2 ratio of 0.848 and -0.737 for amplifications/deletions)." Establishes that high-level amplification/deletion calling uses amplitude cutoffs well outside the ±0.1 noise band, consistent with the CNVkit AMP (≥ +0.7) and DEL (≤ −1.1) extremes.

### GISTIC2 documentation — `-ta` / `-td` parameters

**URL:** https://broadinstitute.github.io/gistic2/
**Accessed:** 2026-06-14 (fetched with WebFetch)
**Authority rank:** 2 (official tool documentation / standard)

**Key Extracted Points:**

1. **Amplification threshold `-ta` (verbatim):** "Threshold for copy number amplifications. Regions with a copy number gain above this positive value are considered amplified. Regions with a copy number gain smaller than this value are considered noise and set to 0." Default = 0.1.
2. **Deletion threshold `-td` (verbatim):** "Threshold for copy number deletions. Regions with a copy number loss below the negative of this positive value are considered deletions." Default = 0.1.

### Seqeron in-repo — `StructuralVariantAnalyzer.DetectCNV` / `SegmentCopyNumber` (SV-CNV-001)

**URL:** in-repo file `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs` (read in session)
**Accessed:** 2026-06-14
**Authority rank:** 3 (sibling repository implementation; overlap check)

**Key Extracted Points:**

1. **Overlap:** SV-CNV-001 already converts read-depth → log2 → **integer** copy number via `round(2 · 2^log2)` (CNVkit `_log2_ratio_to_absolute_pure`) and merges segments. It does NOT classify into discrete CNA states (deep deletion / loss / neutral / gain / amplification). ONCO-CNA-001 is the **oncology classification layer** on top of a log2 ratio: it adds the hard-threshold state calling, which SV-CNV-001 does not provide. The two are complementary, not duplicates; ONCO-CNA-001 reuses the same `n = 2·2^log2` conversion formula (cited to the same CNVkit source) but uses CNVkit's `absolute_threshold` hard-threshold binning rather than simple rounding.

---

## Documented Corner Cases and Failure Modes

### From CNVkit `cnvlib/call.py`

1. **NaN log2 ratio:** treated as a no-call and replaced with the neutral reference copy number (diploid → CN 2, Neutral).
2. **Boundary inclusivity:** the comparison is `log2 <= thresh`, so a value exactly on a threshold is assigned the LOWER copy-number state of the bin (e.g. log2 = −1.1 → CN 0; log2 = 0.7 → CN 3).
3. **Above the last threshold:** copy number is `ceil(2 · 2^log2)`, NOT a fixed value, so high amplifications get progressively larger integer CN (the AMP class).

### From GISTIC2 (Mermel 2011)

1. **Noise band:** |log2| ≤ 0.1 segments are low-amplitude noise; CNVkit's neutral bin (−0.25, 0.2] subsumes this band.

---

## Test Datasets

### Dataset: CNVkit `absolute_threshold` worked values (default thresholds −1.1, −0.25, 0.2, 0.7)

**Source:** CNVkit `cnvlib/call.py` binning rule applied to chosen log2 ratios; absolute CN from `n = 2·2^log2`.

| log2 ratio | Bin rule | Called integer CN | CNA state |
|-----------|----------|-------------------|-----------|
| −2.0 | ≤ −1.1 | 0 | DeepDeletion |
| −1.1 | ≤ −1.1 (boundary) | 0 | DeepDeletion |
| −1.0 | (−1.1, −0.25] | 1 | Loss |
| −0.25 | ≤ −0.25 (boundary) | 1 | Loss |
| 0.0 | (−0.25, 0.2] | 2 | Neutral |
| 0.2 | ≤ 0.2 (boundary) | 2 | Neutral |
| 0.5849625 = log2(3/2) | (0.2, 0.7] | 3 | Gain |
| 0.7 | ≤ 0.7 (boundary) | 3 | Gain |
| 1.0 | > 0.7 → ceil(2·2^1.0)=ceil(4.0)=4 | 4 | Amplification |
| 2.0 | > 0.7 → ceil(2·2^2.0)=ceil(8.0)=8 | 8 | Amplification |

### Dataset: absolute copy-number formula n = 2·2^log2 (`_log2_ratio_to_absolute_pure`, ploidy 2)

**Source:** CNVkit `cnvlib/call.py` `_log2_ratio_to_absolute_pure`.

| log2 ratio | n = 2·2^log2 (continuous) |
|-----------|---------------------------|
| 0.0 | 2.0 |
| 1.0 | 4.0 |
| −1.0 | 1.0 |
| log2(3/2) | 3.0 |
| 0.8 | 2·2^0.8 = 3.4822022… → ceil = 4 |

---

## Assumptions

1. **ASSUMPTION: Diploid (autosomal) reference ploidy = 2.** The classification and the absolute-CN formula use ploidy = 2. This is the CNVkit default for autosomes (`ref_copies = ploidy = 2`) and is the documented diploid baseline; sex chromosomes / haploid references are out of scope for this unit. Source-backed default, but the caller cannot currently override it — recorded as an assumption because a non-diploid baseline would change every output.

---

## Recommendations for Test Coverage

1. **MUST Test:** each of the five CNA states is returned for a representative log2 in its bin (deep deletion, loss, neutral, gain, amplification) — Evidence: CNVkit `absolute_threshold` default thresholds (−1.1, −0.25, 0.2, 0.7).
2. **MUST Test:** boundary inclusivity (`log2 <= thresh`): log2 = −1.1 → CN 0; −0.25 → CN 1; 0.2 → CN 2; 0.7 → CN 3 — Evidence: CNVkit binning loop `if row.log2 <= thresh`.
3. **MUST Test:** absolute CN formula n = 2·2^log2: log2 0→2, 1→4, −1→1, log2(3/2)→3 — Evidence: `_log2_ratio_to_absolute_pure`.
4. **MUST Test:** amplification integer CN above the last threshold uses ceil(2·2^log2): log2 1.0 → 4, log2 2.0 → 8 — Evidence: `absolute_threshold` else-branch.
5. **MUST Test:** NaN log2 → Neutral / CN 2 (no-call) — Evidence: `absolute_threshold` NaN handling.
6. **SHOULD Test:** classifying a batch of log2 ratios preserves input order and length — Rationale: deterministic per-element mapping.
7. **COULD Test:** custom thresholds override default bins — Rationale: CNVkit exposes `-t/--thresholds`.

---

## References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12(4):R41. https://doi.org/10.1186/gb-2011-12-4-r41
2. Talevich E, Shain AH, Botton T, Bastian BC. (CNVkit) `cnvlib/call.py` — `absolute_threshold`, `_log2_ratio_to_absolute_pure`. https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py
3. CNVkit documentation — `call` command threshold method. https://cnvkit.readthedocs.io/en/stable/pipeline.html
4. GISTIC2 documentation — `-ta` / `-td` thresholds. https://broadinstitute.github.io/gistic2/

---

## Change History

- **2026-06-14**: Initial documentation (ONCO-CNA-001).
