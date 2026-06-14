# Evidence Artifact: ONCO-CNA-002

**Test Unit ID:** ONCO-CNA-002
**Algorithm:** Focal Amplification Detection (GISTIC2 length-based focal/broad split + oncogene mapping)
**Date Collected:** 2026-06-14

---

## Online Sources

### Mermel et al. (2011) — GISTIC2.0 (Genome Biology)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
**Accessed:** 2026-06-14 (fetched via WebFetch of the PMC full-text page)
**Authority rank:** 1 (peer-reviewed paper, Genome Biology 12:R41)

**Key Extracted Points:**

1. **Length-based focal/arm-level split:** GISTIC2.0 separates a copy-number profile into arm-level and focal components by length. The fetched text states focal SCNAs are those with *"length < 98% of a chromosome arm"*, while events *"occupying more than 98% of a chromosome arm"* are classified as arm-level. The procedure removes *"all SCNAs occupying more than 98% of a chromosome arm, leaving only the focal events."*
2. **Length is the natural classifier:** *"This reproducible distribution provides a natural basis for classifying events as 'arm-level' and 'focal' based purely on length."* — i.e. GISTIC2.0 moves away from amplitude-based filtering toward length-based filtering.
3. **Amplitude thresholds (historical/low-level):** A low-amplitude threshold of log2 ratio ±0.1 *"eliminates only low-level artifactual segments"*; older high-amplitude thresholds (log2 0.848 amp / −0.737 del) were used in prior versions to exclude arm-level events. GISTIC2.0 itself relies on length, not a single default amplitude cutoff.

### Broad Institute — GISTIC2 documentation (parameter reference)

**URL:** https://broadinstitute.github.io/gistic2/
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (canonical project / reference-implementation documentation)

**Key Extracted Points:**

1. **`broad_len_cutoff` (default 0.98):** verbatim — *"Threshold used to distinguish broad from focal events, given in units of fraction of chromosome arm."* Default value 0.98. This is the configurable form of the 98% figure in the paper.
2. **`t_amp` (default 0.1):** verbatim — *"Threshold for copy number amplifications. Regions with a copy number gain above this positive value are considered amplified."*
3. **`t_del` (default 0.1):** verbatim — *"Threshold for copy number deletions. Regions with a copy number loss below the negative of this positive value are considered deletions."* (out of scope here; recorded for completeness — ONCO-CNA-003 covers deletions.)

### CNVkit — Calling copy number gains and losses (documentation)

**URL:** https://cnvkit.readthedocs.io/en/stable/calling.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (established reference implementation documentation)

**Key Extracted Points:**

1. **Single-copy-gain log2 value:** verbatim — *"a single-copy gain in a perfectly pure, homogeneous sample has a copy ratio of 3/2. In log2 scale, this is log2(3/2) = 0.585."* Establishes that any gain (CN > 2) has log2 > 0; an amplitude threshold of 0.1 is well below a single-copy gain, so it admits all gains as amplitude-positive.
2. **Focal vs non-focal interpretation:** CNVkit's calling chapter and gainloss discussion describe that non-focal amplified segments show surrounding segments with similar copy ratios and are filtered out as artifacts; focal high-amplitude events are the therapeutically actionable ones. Confirms the length/extent-based notion of "focal".

### NCBI Gene — oncogene cytogenetic locations

**URLs (each fetched via WebFetch on 2026-06-14):**
- ERBB2: https://www.ncbi.nlm.nih.gov/gene/2064 → Location "17q12"
- MYC: https://www.ncbi.nlm.nih.gov/gene/4609 → Location "8q24.21"
- EGFR: https://www.ncbi.nlm.nih.gov/gene/1956 → Location "7p11.2"
- CCND1: https://www.ncbi.nlm.nih.gov/gene/595 → Location "11q13.3"
- MDM2: https://www.ncbi.nlm.nih.gov/gene/4193 → Location "12q15"
- CDK4: https://www.ncbi.nlm.nih.gov/gene/1019 → Location "12q14.1"

**Authority rank:** 5 (curated database, NCBI Gene).

**Key Extracted Points:**

1. **Cytogenetic locations (chromosome + arm) for the registry oncogene panel:** ERBB2 = chr17 q-arm (17q12); MYC = chr8 q-arm (8q24.21); EGFR = chr7 p-arm (7p11.2); CCND1 = chr11 q-arm (11q13.3); MDM2 = chr12 q-arm (12q15); CDK4 = chr12 q-arm (12q14.1). The chromosome+arm prefix (e.g. "17q", "7p", "8q") is what a focal-amplification segment's arm label is matched against for oncogene mapping.

---

## Documented Corner Cases and Failure Modes

### From Mermel et al. (2011) / GISTIC2 docs

1. **Whole-arm event excluded as focal:** a segment occupying ≥ 98% of its arm is arm-level/broad, not focal — it must NOT be reported as a focal amplification even if highly amplified.
2. **Boundary at exactly 98%:** the paper's wording is "more than 98%" → arm-level, "less than 98%" → focal. The 98% point itself is the cutoff; the doc parameter is "fraction of chromosome arm" = 0.98. We treat a fraction strictly less than the cutoff as focal (length/arm < 0.98 → focal).

### From CNVkit docs

1. **Amplitude below threshold:** a segment whose log2 gain does not exceed `t_amp` (0.1) is not amplified at all and is excluded regardless of length.

---

## Test Datasets

### Dataset: Synthetic GISTIC2-rule worked segments

**Source:** Derivation from Mermel et al. (2011) 98%-of-arm rule and GISTIC2 `t_amp` = 0.1; oncogene arms from NCBI Gene.

| Segment | Arm | ArmLength (bp) | Start–End (bp) | SegLen | SegLen/Arm | log2 | Amplified? (>0.1) | Focal? (<0.98) | Focal amplification? |
|---------|-----|---------------|----------------|--------|-----------|------|-------------------|----------------|----------------------|
| A (ERBB2) | 17q | 1,000,000 | 100,000–600,000 | 500,000 | 0.50 | 1.0 | yes | yes | **yes** |
| B (whole arm) | 8q | 1,000,000 | 0–990,000 | 990,000 | 0.99 | 1.5 | yes | no | no (arm-level) |
| C (low amp) | 7p | 1,000,000 | 0–300,000 | 300,000 | 0.30 | 0.05 | no | yes | no (not amplified) |
| D (boundary) | 11q | 1,000,000 | 0–980,000 | 980,000 | 0.98 | 1.0 | yes | no (= cutoff) | no (not < 0.98) |

---

## Assumptions

1. **ASSUMPTION: amplitude test for "amplified".** GISTIC2.0 itself classifies focal-vs-broad purely by length; "amplification" (gain direction with positive amplitude) is taken from the GISTIC2 `t_amp` = 0.1 parameter (gain above +0.1). This is source-backed (GISTIC2 docs) but combining the length rule (paper) with the `t_amp` amplitude rule (docs) into a single `DetectFocalAmplifications` predicate is the integration choice of this unit; it is documented rather than invented.
2. **ASSUMPTION: arm fraction provided as input.** GISTIC2 derives chromosome-arm boundaries from the genome assembly/cytoband file. This unit does not bundle a cytoband table; the caller supplies each segment's arm label and the arm's length, and the algorithm computes the segment-length / arm-length fraction. The 0.98 cutoff and the amplitude rule are unchanged.

---

## Recommendations for Test Coverage

1. **MUST Test:** A high-amplitude segment spanning < 98% of its arm is reported as a focal amplification; a ≥ 98%-of-arm segment is not. — Evidence: Mermel et al. (2011) 98%-of-arm rule; GISTIC2 `broad_len_cutoff` = 0.98.
2. **MUST Test:** A gain whose log2 does not exceed `t_amp` (0.1) is not reported even if focal in length. — Evidence: GISTIC2 `t_amp` = 0.1.
3. **MUST Test:** Boundary at exactly 0.98 fraction is NOT focal (treated as arm-level). — Evidence: paper "more than 98% ⇒ arm-level".
4. **MUST Test:** `IdentifyAmplifiedOncogenes` maps a focal amplification on arm "17q" to ERBB2, "8q" to MYC, "7p" to EGFR, "11q" to CCND1, "12q" to MDM2 and CDK4. — Evidence: NCBI Gene locations.
5. **SHOULD Test:** A non-amplified (neutral/loss) segment is never an oncogene amplification. — Rationale: only focal amplifications feed the mapper.
6. **COULD Test:** Null / empty input handling. — Rationale: documented failure modes.

---

## References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12:R41. https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
2. Broad Institute. GISTIC2 documentation (parameter reference: `broad_len_cutoff`, `t_amp`, `t_del`). https://broadinstitute.github.io/gistic2/
3. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit — Calling copy number gains and losses. https://cnvkit.readthedocs.io/en/stable/calling.html
4. NCBI Gene database: ERBB2 (2064), MYC (4609), EGFR (1956), CCND1 (595), MDM2 (4193), CDK4 (1019). https://www.ncbi.nlm.nih.gov/gene/

---

## Change History

- **2026-06-14**: Initial documentation.
