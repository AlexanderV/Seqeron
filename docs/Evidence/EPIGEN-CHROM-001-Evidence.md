# Evidence Artifact: EPIGEN-CHROM-001

**Test Unit ID:** EPIGEN-CHROM-001
**Algorithm:** Chromatin State Prediction from histone modification marks
**Date Collected:** 2026-06-13

---

## Online Sources

### Ernst & Kellis — ChromHMM software documentation (binarization model)

**URL:** http://compbio.mit.edu/ChromHMM/
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation / canonical project documentation)
**Retrieved how:** WebSearch `Ernst Kellis ChromHMM discovery characterization chromatin states combinatorial histone modification Nature Methods 2012`, then WebFetch of `http://compbio.mit.edu/ChromHMM/` (after a 301 redirect from the https URL).

**Key Extracted Points:**

1. **Binary mark model:** ChromHMM "is based on a multivariate Hidden Markov Model that explicitly models the presence or absence of each chromatin mark." Each mark contributes to a state as a present/absent (1/0) call, not a continuous value.
2. **Binarization step:** Before learning states, raw signal is converted to present/absent calls per mark via the `BinarizeBed` / `BinarizeBam` commands; the `LearnModel` step then operates on these binary calls. This justifies treating chromatin-state prediction as a function of binarized (present/absent) marks rather than raw magnitudes.

### Ernst & Kellis (2012) — ChromHMM: automating chromatin-state discovery and characterization

**URL:** https://www.nature.com/articles/nmeth.1906 (DOI https://doi.org/10.1038/nmeth.1906)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed, Nature Methods)
**Retrieved how:** WebSearch `Ernst Kellis ChromHMM discovery characterization chromatin states combinatorial histone modification Nature Methods 2012` returning the Nature Methods article record (vol 9, pp. 215–216, 2012).

**Key Extracted Points:**

1. **Combinatorial patterns define states:** Chromatin states are the "major re-occurring combinatorial and spatial patterns of marks"; states "often capture known classes of genomic elements such as promoters, enhancers, transcribed, repressed, and repetitive regions."
2. **Marks model presence/absence:** the model "explicitly models the presence or absence of each chromatin mark" — confirming the binary treatment used here.

### Roadmap Epigenomics — Chromatin state learning (15-state and 18-state model definitions)

**URL:** https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
**Accessed:** 2026-06-13
**Authority rank:** 2 (consortium standard / official Roadmap Epigenomics reference)
**Retrieved how:** WebSearch result for the Roadmap 18-state model, then WebFetch of the chromatin state learning portal page.

**Key Extracted Points:**

1. **Core 15-state model marks:** five marks — H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3.
2. **Expanded 18-state model marks:** six marks — H3K4me3, H3K4me1, **H3K27ac**, H3K36me3, H3K27me3, H3K9me3. These are exactly the six marks taken by the implementation's `PredictChromatinState`.
3. **State → characteristic mark mapping (core model):**
   - TssA "Active TSS" → H3K4me3
   - Tx "Strong transcription" / TxWk "Weak transcription" → H3K36me3
   - Enh "Enhancers" / EnhG "Genic enhancers" → H3K4me1
   - Het "Heterochromatin" → H3K9me3
   - TssBiv "Bivalent/Poised TSS" → H3K4me3 + H3K27me3
   - EnhBiv "Bivalent Enhancer" → H3K4me1 + H3K27me3
   - ReprPC "Repressed PolyComb" / ReprPCWk → H3K27me3
   - Quies "Quiescent/Low" → no enrichment
4. **18-state addition:** adding H3K27ac subdivides enhancers/TSS into active vs weak (active enhancer = H3K4me1 with H3K27ac present).

### H3K4me3 — Wikipedia (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/H3K4me3
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary: Liang et al. 2004, PNAS)
**Retrieved how:** WebFetch of the Wikipedia H3K4me3 article.

**Key Extracted Points:**

1. **Active promoter mark:** H3K4me3 "is highly enriched at active promoters near transcription start sites (TSS) and positively correlated with transcription."
2. **Primary citation:** Liang et al. (2004), PNAS 101(19):7357–7362, "Distinct localization of histone H3 acetylation and H3-K4 methylation to the transcription start sites in the human genome."

### H3K4me1 — Wikipedia (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/H3K4me1
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary: Rada-Iglesias 2018, Nature Genetics)
**Retrieved how:** WebFetch of the Wikipedia H3K4me1 article.

**Key Extracted Points:**

1. **Enhancer mark:** "H3K4me1 is enriched at active and primed enhancers."
2. **Primary citation:** Rada-Iglesias A (2018), Nature Genetics 50(1):4–5, "Is H3K4me1 at enhancers correlative or causative?"

### H3K27ac — Wikipedia (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/H3K27ac
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary: Creyghton et al. 2010, PNAS)
**Retrieved how:** WebFetch of the Wikipedia H3K27ac article.

**Key Extracted Points:**

1. **Active enhancer mark:** H3K27ac is "an active enhancer mark" and "separates active from poised enhancers and predicts developmental state." H3K27ac present on an enhancer (H3K4me1) marks it as active; absent leaves it weak/poised.
2. **Primary citation:** Creyghton MP et al. (2010), PNAS 107(50):21931–21936, "Histone H3K27ac separates active from poised enhancers and predicts developmental state."

### H3K27me3 — Wikipedia (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/H3K27me3
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary: Ferrari et al. 2014, Molecular Cell)
**Retrieved how:** WebFetch of the Wikipedia H3K27me3 article.

**Key Extracted Points:**

1. **Polycomb repression mark:** H3K27me3 "is associated with the downregulation of nearby genes via the formation of heterochromatic regions"; deposited by the Polycomb complex (PRC2). Roadmap maps it to the ReprPC "Repressed PolyComb" state.
2. **Primary citation:** Ferrari KJ et al. (2014), Molecular Cell 53(1):49–62.

### H3K9me3 — Wikipedia (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/H3K9me3
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary: Nicetto et al. 2019, Science)
**Retrieved how:** WebFetch of the Wikipedia H3K9me3 article.

**Key Extracted Points:**

1. **Heterochromatin mark:** H3K9me3 "is often associated with heterochromatin." Roadmap maps it to the Het "Heterochromatin" state.
2. **Primary citation:** Nicetto D et al. (2019), Science 363(6424):294–297.

### H3K36me3 — Wikipedia + Journal of Human Genetics review

**URL:** https://en.wikipedia.org/wiki/H3K36me3 ; https://www.nature.com/articles/jhg201366
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia) / 1 (Kimura 2013, J Hum Genet review — record retrieved; full text gated)
**Retrieved how:** WebFetch of the Wikipedia H3K36me3 article; WebSearch surfaced the J Hum Genet review whose abstract states gene bodies of actively transcribed genes carry H3K36me3 (full text behind a login redirect, so only the search-surfaced statement is used; the Roadmap Tx/TxWk mapping is the load-bearing source).

**Key Extracted Points:**

1. **Transcribed gene body mark:** H3K36me3 is "often associated with gene bodies"; the J Hum Genet review states "Gene bodies of actively transcribed genes are associated with trimethylated H3K36 (H3K36me3)." Roadmap maps it to the Tx "Strong transcription" / TxWk states.

---

## Documented Corner Cases and Failure Modes

### From ChromHMM (Ernst & Kellis 2012) / Roadmap

1. **No mark present (Quiescent/Low):** When no mark exceeds the presence call, Roadmap assigns the Quies "Quiescent/Low — no enrichment" state. Our analogue is `LowSignal`.
2. **Co-occurring active + repressive marks (bivalent):** H3K4me3 together with H3K27me3 is the canonical *bivalent/poised* TSS signature (Roadmap TssBiv), not a contradiction — must be classified as bivalent rather than as either active or repressed alone.
3. **Combinatorial precedence:** A state is defined by a *combination* of marks; a promoter signature (H3K4me3) must take precedence over an enhancer signature (H3K4me1) when both are present at the same locus, matching the Roadmap TSS-vs-enhancer distinction (promoter marks dominate).

### From the binarization model

1. **Magnitude is not ordinal beyond the threshold:** Once a mark is "present" (above the call), increasing its magnitude does not change the state. State is a function of the *set* of present marks (ChromHMM binary model), so two inputs with the same present/absent pattern yield the same state.

---

## Test Datasets

### Dataset: Canonical Roadmap chromatin-state signatures (present/absent marks)

**Source:** Roadmap Epigenomics chromatin state learning (15/18-state models); per-mark primaries above.

| Present marks (above call) | Expected state |
|----------------------------|----------------|
| H3K4me3 (+H3K27ac) | ActivePromoter (TssA) |
| H3K4me3 only | ActivePromoter (TssA) |
| H3K4me1 + H3K27ac | ActiveEnhancer (active Enh) |
| H3K4me1, no H3K27ac | WeakEnhancer (poised/weak Enh) |
| H3K36me3 | Transcribed (Tx) |
| H3K27me3 (alone) | Repressed (ReprPC) |
| H3K9me3 (alone) | Heterochromatin (Het) |
| H3K4me3 + H3K27me3 | BivalentPromoter (TssBiv) |
| H3K4me1 + H3K27me3 | BivalentEnhancer (EnhBiv) |
| none | LowSignal (Quies) |

---

## Assumptions

1. **ASSUMPTION: presence-call threshold value** — ChromHMM performs binarization with a Poisson background model from raw read counts; a single fixed numeric threshold on an already-normalized [0,1] signal is not specified by the sources. The implementation exposes the presence threshold as a caller-supplied parameter (default 0.5 on a normalized [0,1] enrichment signal) and documents it; the *state-assignment logic given the present/absent pattern* is fully source-backed and is what the tests verify. Tests choose mark magnitudes unambiguously above/below the call so the result does not depend on the exact default.

2. **ASSUMPTION: precedence of promoter over enhancer when both H3K4me3 and H3K4me1 are present without repressive marks** — Roadmap separates TSS (H3K4me3) from enhancer (H3K4me1) states; when both active marks co-occur at one locus we classify as promoter (H3K4me3 dominates), consistent with TSS states ranking above enhancer states in the Roadmap mnemonic ordering. Marked as assumption because Roadmap derives this from spatial HMM context, not a single-locus rule.

---

## Recommendations for Test Coverage

1. **MUST Test:** Each canonical single/combination mark signature maps to its Roadmap state (table above). — Evidence: Roadmap state definitions + per-mark primaries.
2. **MUST Test:** Bivalent signature (H3K4me3 + H3K27me3) → BivalentPromoter, not ActivePromoter or Repressed. — Evidence: Roadmap TssBiv.
3. **MUST Test:** No mark present → LowSignal. — Evidence: Roadmap Quies.
4. **MUST Test:** Binary invariance — same present/absent pattern, different magnitudes → same state. — Evidence: ChromHMM binary model.
5. **SHOULD Test:** `AnnotateHistoneModifications` assigns each region the state of its single mark (delegation/per-mark mapping). — Rationale: it labels regions by mark identity.
6. **SHOULD Test:** `FindAccessibleRegions` merges contiguous above-threshold positions into one region and excludes sub-`minWidth` regions. — Rationale: peak-calling contract.
7. **COULD Test:** Negative/zero signals treated as absent. — Rationale: robustness.

---

## References

1. Ernst J, Kellis M (2012). ChromHMM: automating chromatin-state discovery and characterization. Nature Methods 9(3):215–216. https://doi.org/10.1038/nmeth.1906
2. Ernst J, Kellis M. ChromHMM software and manual (binarization: present/absent marks). http://compbio.mit.edu/ChromHMM/
3. Roadmap Epigenomics Consortium. Chromatin state learning (15-state core and 18-state expanded models). https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
4. Liang G et al. (2004). Distinct localization of histone H3 acetylation and H3-K4 methylation to the transcription start sites in the human genome. PNAS 101(19):7357–7362. https://doi.org/10.1073/pnas.0401866101
5. Rada-Iglesias A (2018). Is H3K4me1 at enhancers correlative or causative? Nature Genetics 50(1):4–5. https://doi.org/10.1038/s41588-017-0018-3
6. Creyghton MP et al. (2010). Histone H3K27ac separates active from poised enhancers and predicts developmental state. PNAS 107(50):21931–21936. https://doi.org/10.1073/pnas.1016071107
7. Ferrari KJ et al. (2014). Polycomb-dependent H3K27me1 and H3K27me2 regulate active transcription and enhancer fidelity. Molecular Cell 53(1):49–62. https://doi.org/10.1016/j.molcel.2013.10.030
8. Nicetto D et al. (2019). H3K9me3-heterochromatin loss at protein-coding genes enables developmental lineage specification. Science 363(6424):294–297. https://doi.org/10.1126/science.aau0583
9. Kimura H (2013). Histone modifications for human epigenome analysis. Journal of Human Genetics 58(7):439–445. https://doi.org/10.1038/jhg.2013.66

---

## Change History

- **2026-06-13**: Initial documentation.
