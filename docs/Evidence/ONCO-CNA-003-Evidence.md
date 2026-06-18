# Evidence Artifact: ONCO-CNA-003

**Test Unit ID:** ONCO-CNA-003
**Algorithm:** Homozygous (Deep) Deletion Detection
**Date Collected:** 2026-06-14

---

## Online Sources

### cBioPortal — File Formats: Discrete Copy Number Data

**URL:** https://docs.cbioportal.org/file-formats/
**Retrieved by:** WebFetch of the URL above (prompt: discrete copy number values and the meaning of −2/−1/0/1/2 and Deep Deletion / homozygous deletion).
**Accessed:** 2026-06-14
**Authority rank:** 5 (well-maintained bioinformatics database / reference platform; the de-facto discrete-CNA convention)

**Key Extracted Points:**

1. **Discrete CNA scale (verbatim):** `"-2" is a deep loss, possibly a homozygous deletion`; `"-1" is a single-copy loss (heterozygous deletion)`; `"0" is diploid`; `"1" indicates a low-level gain`; `"2" is a high-level amplification.`
2. **−2 ↔ homozygous deletion:** The value −2 ("Deep Deletion") is explicitly the call associated with a (possible) **homozygous deletion** — i.e. the deepest discrete loss state.

### cBioPortal — FAQ (DNA: Mutations, Copy Number & Fusions)

**URL:** https://docs.cbioportal.org/user-guide/faq/
**Retrieved by:** WebFetch of the URL above (prompt: FAQ entries defining −2 Deep Deletion / homozygous deletion, −1 Shallow, 0, 1, 2).
**Accessed:** 2026-06-14
**Authority rank:** 5

**Key Extracted Points:**

1. **Deep Deletion (verbatim):** `"-2 or Deep Deletion indicates a deep loss, possibly a homozygous deletion"`.
2. **Shallow Deletion (verbatim):** `"-1 or Shallow Deletion indicates a shallow loss, possibley a heterozygous deletion"`.
3. **Neutral / gain / amp (verbatim):** `"0 is diploid"`; `"1 or Gain indicates a low-level gain (a few additional copies, often broad)"`; `"2 or Amplification indicate a high-level amplification (more copies, often focal)"`.

### Cheng et al. (2017) — Pan-cancer analysis of homozygous deletions (Nat Commun)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/
**Retrieved by:** WebSearch ("homozygous deletion total copy number zero both alleles lost tumor suppressor definition") → WebFetch of the PMC article (prompt: definition of homozygous deletion in terms of copy number; distinction from hemizygous).
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary literature)

**Key Extracted Points:**

1. **Homozygous deletion = total copy number 0 (verbatim):** homozygous deletions are `regions having zero copies of both alleles in the tumour cells` — i.e. absolute/total copy number 0.
2. **Distinction from hemizygous:** homozygous deletion = complete absence of **both** gene copies, whereas a hemizygous (single-copy / heterozygous) deletion loses only **one** allele. Homozygous deletions "require two independent hits".
3. **Biological role:** homozygous deletions are rare and recurrently target **tumour suppressor** genes; the discovery of cancer-specific homozygous deletions delineated tumour suppressors including RB1, CDKN2A, PTEN.

### CNVkit — `cnvlib/call.py` `absolute_threshold` (integer CN, shared with ONCO-CNA-001)

**URL:** (established and retrieved in ONCO-CNA-001) https://cnvkit.readthedocs.io/ ; source `cnvlib/call.py`.
**Retrieved by:** Re-used from ONCO-CNA-001 Evidence (the integer-CN mapping is already in `OncologyAnalyzer`); the integer copy number 0 ⇒ `DeepDeletion` state is the existing `ClassifyCopyNumber` contract.
**Accessed:** (via ONCO-CNA-001) — see docs/Evidence/ONCO-CNA-001-Evidence.md
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Integer CN 0 ⇒ deep/homozygous deletion:** the existing `OncologyAnalyzer.CopyNumberState.DeepDeletion` corresponds to integer copy number 0 (log2 ≤ −1.1 by CNVkit `absolute_threshold` default). A homozygous deletion is therefore exactly the integer copy number 0 state.

### NCBI Gene — tumour-suppressor cytogenetic locations

**URL / Retrieved by (WebFetch each NCBI Gene record, prompt = cytogenetic Location string):**
- TP53  — https://www.ncbi.nlm.nih.gov/gene/7157 → **17p13.1**
- RB1   — https://www.ncbi.nlm.nih.gov/gene/5925 → **13q14.2**
- CDKN2A— https://www.ncbi.nlm.nih.gov/gene/1029 → **9p21.3**
- PTEN  — https://www.ncbi.nlm.nih.gov/gene/5728 → **10q23.31**
- BRCA1 — https://www.ncbi.nlm.nih.gov/gene/672  → **17q21.31**
- BRCA2 — https://www.ncbi.nlm.nih.gov/gene/675  → **13q13.1**

**Accessed:** 2026-06-14
**Authority rank:** 5 (curated database — NCBI Gene)

**Key Extracted Points:**

1. **Arm assignment** (chromosome + p/q from the band): TP53→17p, RB1→13q, CDKN2A→9p, PTEN→10q, BRCA1→17q, BRCA2→13q. These are the arms used to map a homozygous-deletion segment to a deleted tumour suppressor.

---

## Documented Corner Cases and Failure Modes

### From cBioPortal

1. **Shallow vs deep:** a single-copy loss (−1, heterozygous) is NOT a homozygous deletion; only −2 / total CN 0 qualifies. A loss/shallow-deletion segment must NOT be reported as a homozygous deletion.
2. **Putative calls:** discrete calls are putative and unreviewed; purity/ploidy differences cause false +/−. (Interpretation caveat; does not change the CN-0 definition.)

### From Cheng et al. (2017)

1. **Both alleles required:** a region with one remaining copy (CN ≥ 1) is hemizygous, not homozygous — the defining line is total copy number exactly 0.

---

## Test Datasets

### Dataset: Discrete CNA scale (cBioPortal)

**Source:** cBioPortal File-Formats / FAQ (verbatim above)

| Discrete value | Meaning | Homozygous deletion? |
|----------------|---------|----------------------|
| −2 | Deep Deletion (possible homozygous deletion) | Yes |
| −1 | Shallow / single-copy (heterozygous) loss | No |
| 0 | Diploid | No |
| 1 | Gain | No |
| 2 | Amplification | No |

### Dataset: Integer copy number ⇒ state (CNVkit `absolute_threshold`, via ONCO-CNA-001)

**Source:** CNVkit `cnvlib/call.py`; `OncologyAnalyzer.ClassifyCopyNumber`

| log2 (default thresholds −1.1,−0.25,0.2,0.7) | Integer CN | State | Homozygous deletion |
|---|---|---|---|
| −2.0 (≤ −1.1) | 0 | DeepDeletion | **Yes** |
| −0.5 (−1.1 < log2 ≤ −0.25) | 1 | Loss | No |
| 0.0 | 2 | Neutral | No |

### Dataset: Tumour-suppressor arms (NCBI Gene)

**Source:** NCBI Gene cytogenetic locations (verbatim above)

| Gene | Cytoband | Arm |
|------|----------|-----|
| TP53 | 17p13.1 | 17p |
| RB1 | 13q14.2 | 13q |
| CDKN2A | 9p21.3 | 9p |
| PTEN | 10q23.31 | 10q |
| BRCA1 | 17q21.31 | 17q |
| BRCA2 | 13q13.1 | 13q |

---

## Assumptions

1. **ASSUMPTION: Curated tumour-suppressor panel is caller-supplied / fixed.** The registry pins the panel to TP53, RB1, CDKN2A, PTEN, BRCA1, BRCA2. Arm membership is source-backed (NCBI Gene); the *choice* of panel is a registry-supplied curated list (analogous to the ONCO-CNA-002 oncogene panel), not a derived threshold. This is non-correctness-affecting for the deletion-detection logic itself (it only labels which arms map to which gene names).
2. **ASSUMPTION: Homozygous deletion is identified at the integer-CN level via the existing CN-0 (DeepDeletion) classification.** Per cBioPortal a homozygous deletion is the deepest discrete loss (−2) and per Cheng et al. it is total CN 0; the repository already realizes "integer CN 0 ⇒ DeepDeletion" (CNVkit, ONCO-CNA-001). A segment is therefore a homozygous deletion iff its classified integer copy number is 0. No new numeric threshold is invented.

---

## Recommendations for Test Coverage

1. **MUST Test:** A segment whose log2 classifies to integer CN 0 (e.g. log2 = −2.0 with default thresholds) is reported as a homozygous deletion. — Evidence: cBioPortal (−2 = Deep Deletion = homozygous); Cheng et al. (total CN 0).
2. **MUST Test:** A single-copy loss segment (integer CN 1, e.g. log2 = −0.5) is NOT reported as a homozygous deletion. — Evidence: cBioPortal (−1 = shallow/heterozygous, not homozygous); Cheng et al. (one allele remains).
3. **MUST Test:** Neutral / gain / amplification segments are NOT reported. — Evidence: cBioPortal scale.
4. **MUST Test:** `IdentifyDeletedTumorSuppressors` maps a homozygous-deletion segment on each arm to the right gene(s): 17p→TP53, 13q→RB1+BRCA2, 9p→CDKN2A, 10q→PTEN, 17q→BRCA1. — Evidence: NCBI Gene locations.
5. **MUST Test:** Custom thresholds change the CN-0 boundary (e.g. raising the deletion cutoff makes a previously-CN-1 segment CN 0). — Evidence: CNVkit `absolute_threshold` thresholds are parameters.
6. **SHOULD Test:** Boundary log2 exactly at the deletion cutoff (−1.1) is CN 0 (≤ cutoff), so homozygous. — Rationale: off-by-one / inclusive boundary per CNVkit "less than or equal to each threshold".
7. **SHOULD Test:** Detection preserves input order and reports each qualifying segment once. — Rationale: order-preserving filter (mirror DetectFocalAmplifications).
8. **COULD Test:** Empty input ⇒ empty result; null input ⇒ ArgumentNullException; invalid arm length / End≤Start ⇒ ArgumentException. — Rationale: mirror sibling validation in ONCO-CNA-002.

---

## References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12:R41. https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
2. Cheng J, Demeulemeester J, Wedge DC, et al. (2017). Pan-cancer analysis of homozygous deletions in primary tumours uncovers rare tumour suppressors. Nature Communications 8:1221. https://pmc.ncbi.nlm.nih.gov/articles/PMC5663922/
3. cBioPortal — Discrete Copy Number data file format. https://docs.cbioportal.org/file-formats/ (accessed 2026-06-14)
4. cBioPortal — FAQ: meaning of Amplification / Gain / Deep Deletion / Shallow Deletion / −2..2. https://docs.cbioportal.org/user-guide/faq/ (accessed 2026-06-14)
5. Talevich E, Shain AH, Botton T, Bastian BC. (2016). CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol 12(4):e1004873. `cnvlib/call.py` `absolute_threshold`. https://cnvkit.readthedocs.io/
6. NCBI Gene records (accessed 2026-06-14): TP53 https://www.ncbi.nlm.nih.gov/gene/7157 ; RB1 https://www.ncbi.nlm.nih.gov/gene/5925 ; CDKN2A https://www.ncbi.nlm.nih.gov/gene/1029 ; PTEN https://www.ncbi.nlm.nih.gov/gene/5728 ; BRCA1 https://www.ncbi.nlm.nih.gov/gene/672 ; BRCA2 https://www.ncbi.nlm.nih.gov/gene/675

---

## Change History

- **2026-06-14**: Initial documentation.
