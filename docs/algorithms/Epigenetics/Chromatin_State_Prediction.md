# Chromatin State Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-CHROM-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Reference |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Chromatin state prediction labels a genomic locus as a regulatory element class
(active promoter, enhancer, transcribed gene body, Polycomb-repressed, heterochromatin,
bivalent, or quiescent) from the combination of histone-modification signals present at
that locus. This implementation reproduces the deterministic, source-defined core of the
ChromHMM/Roadmap Epigenomics approach: each mark is binarized to present/absent and the
present-mark *set* is mapped to its canonical chromatin state [1][3]. It is a rule-based,
specification-driven classifier — not the full unsupervised HMM, whose emission and
transition probabilities are learned per dataset.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Post-translational histone modifications mark functionally distinct chromatin. The
canonical marks and their meanings: H3K4me3 marks active promoters / transcription start
sites [4]; H3K4me1 marks enhancers [5]; H3K27ac marks active (vs poised) enhancers [6];
H3K36me3 marks transcribed gene bodies [9]; H3K27me3 is the Polycomb repressive mark [7];
H3K9me3 marks constitutive heterochromatin [8]. Recurring *combinations* of these marks
define chromatin states [1].

### 2.2 Core Model

ChromHMM "explicitly models the presence or absence of each chromatin mark" [1]; raw
signal is binarized before state learning (the `BinarizeBed`/`BinarizeBam` step) [2]. A
state is the combinatorial pattern of present marks. The Roadmap Epigenomics 15-state core
model uses H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3; the 18-state expanded model adds
H3K27ac [3]. The canonical mark-set → state mapping used here [3]:

| Present marks | State (Roadmap mnemonic) |
|---------------|--------------------------|
| H3K4me3 + H3K27me3 | BivalentPromoter (TssBiv) |
| H3K4me1 + H3K27me3 | BivalentEnhancer (EnhBiv) |
| H3K4me3 (± H3K4me1) | ActivePromoter (TssA) |
| H3K4me1 + H3K27ac | ActiveEnhancer (active Enh) |
| H3K4me1, no H3K27ac | WeakEnhancer (Enh) |
| H3K36me3 | Transcribed (Tx) |
| H3K27me3 (alone) | Repressed (ReprPC) |
| H3K9me3 (alone) | Heterochromatin (Het) |
| none | LowSignal (Quies) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | State depends only on the present/absent mark pattern; equal patterns ⇒ equal state | binary mark model [1][2] |
| INV-02 | No mark present ⇒ `LowSignal` | Roadmap Quies state [3] |
| INV-03 | H3K4me3 ∧ H3K27me3 ⇒ `BivalentPromoter` (overrides plain active/repressed) | Roadmap TssBiv [3] |
| INV-04 | Return value is always a defined `ChromatinState` (total function) | exhaustive rules end in `LowSignal` |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This (rule-based signature) | ChromHMM (full) |
|--------|----------------------------|-----------------|
| Parameters | fixed canonical mark→state map | emission/transition probabilities learned per dataset |
| Spatial context | none (per-locus) | HMM transitions across the genome |
| Determinism | deterministic, dataset-independent | depends on training data |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| h3k4me3 | double | required | H3K4me3 signal | normalized enrichment, typ. [0,1] |
| h3k4me1 | double | required | H3K4me1 signal | as above |
| h3k27ac | double | required | H3K27ac signal | as above |
| h3k36me3 | double | required | H3K36me3 signal | as above |
| h3k27me3 | double | required | H3K27me3 signal | as above |
| h3k9me3 | double | required | H3K9me3 signal | as above |
| presenceThreshold | double | 0.5 | inclusive present call (`signal >= threshold`) | normalized scale |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | ChromatinState | predicted state for the locus |

### 3.3 Preconditions and Validation

Total function: any six `double` inputs yield a defined `ChromatinState`. Signals below
`presenceThreshold` (including 0 and negatives) are treated as absent; an all-absent input
returns `LowSignal`. No nulls (value-type parameters). `AnnotateHistoneModifications`
matches mark names case-insensitively; unrecognised marks map to `LowSignal`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Binarize each of the six marks to present/absent at `presenceThreshold` [1][2].
2. If an active mark co-occurs with H3K27me3, return the matching bivalent state [3].
3. Otherwise apply the canonical signature priority: promoter (H3K4me3) → enhancer
   (H3K4me1, active if H3K27ac present) → transcribed (H3K36me3) → repressed (H3K27me3) →
   heterochromatin (H3K9me3) → low signal [3].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

The mark→state table in §2.2 is the reference table; every row is sourced from the Roadmap
state definitions [3] and per-mark primaries [4]–[9]. Promoter precedence over enhancer
when both H3K4me3 and H3K4me1 are present reflects the Roadmap TSS-above-Enh ordering
(recorded as an assumption — see §5.4).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PredictChromatinState | O(1) | O(1) | constant six-mark decision |
| AnnotateHistoneModifications | O(n) | O(1) streamed | n = number of regions |
| FindAccessibleRegions | O(n log n) | O(n) | dominated by sorting positions |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.PredictChromatinState(...)`: six-mark signature → `ChromatinState`.
- `EpigeneticsAnalyzer.AnnotateHistoneModifications(...)`: labels each region by its single mark's state.
- `EpigeneticsAnalyzer.FindAccessibleRegions(...)`: peak-calls an accessibility (ATAC-seq-like) track.

### 5.2 Current Behavior

`presenceThreshold` defaults to 0.5 on a normalized [0,1] signal and is an inclusive call
(`signal >= threshold`). The search/peak-calling in `FindAccessibleRegions` is a single
left-to-right scan over position-sorted samples merging gaps ≤ `maxGap`; the repository
suffix tree is **not** applicable (this is threshold-based numeric peak calling over a
signal track, not exact substring matching), so a linear scan is used. `PeakType` is a
cosmetic strength label and does not affect which regions are returned.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Present/absent (binary) treatment of each mark [1][2].
- Canonical mark-set → state mapping from the Roadmap 15/18-state models [3].
- Bivalent state detection for H3K4me3+H3K27me3 and H3K4me1+H3K27me3 [3].

**Intentionally simplified:**

- Per-locus rule-based classification instead of the learned multivariate HMM; **consequence:** no spatial/transition modelling and no dataset-specific emission probabilities — output is a deterministic function of the present marks only.
- Single fixed presence threshold instead of ChromHMM's Poisson background binarization; **consequence:** the caller must supply normalized signals and may tune the threshold.

**Not implemented:**

- Model learning (`LearnModel`), genome segmentation, and the full 15/18-state granularity (sub-states like TssAFlnk, EnhG1/G2); **users should rely on:** ChromHMM for de novo state discovery.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | presence-call threshold value | Assumption | controls binarization; state logic itself is sourced | accepted | default 0.5; tests use unambiguous magnitudes |
| 2 | H3K4me3 precedence over H3K4me1 at one locus | Assumption | classification when both active marks co-occur | accepted | Roadmap ranks TSS above Enh [3] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No mark ≥ threshold | `LowSignal` | Roadmap Quies [3] (INV-02) |
| H3K4me3 + H3K27me3 | `BivalentPromoter` | Roadmap TssBiv [3] (INV-03) |
| H3K4me3 + H3K4me1 (no repressive) | `ActivePromoter` | TSS ranks above Enh [3] |
| Negative / zero signals | treated as absent | below presence call |
| Mark exactly at threshold | counts as present | inclusive `>=` |

### 6.2 Limitations

Per-locus only — no spatial context or genome segmentation; assumes pre-normalized input
signals; does not reproduce dataset-specific ChromHMM probabilities or the finer-grained
sub-states. Not a substitute for training a model on real ChIP-seq data.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// H3K4me3 present with H3K27me3 present → bivalent/poised promoter (Roadmap TssBiv).
var state = EpigeneticsAnalyzer.PredictChromatinState(
    h3k4me3: 0.9, h3k4me1: 0.0, h3k27ac: 0.0,
    h3k36me3: 0.0, h3k27me3: 0.8, h3k9me3: 0.0);
// state == ChromatinState.BivalentPromoter
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_ChromatinState_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_ChromatinState_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [EPIGEN-CHROM-001-Evidence.md](../../../docs/Evidence/EPIGEN-CHROM-001-Evidence.md)
- Related algorithms: [CpG_Site_Detection](./CpG_Site_Detection.md)

## 8. References

1. Ernst J, Kellis M. 2012. ChromHMM: automating chromatin-state discovery and characterization. Nature Methods 9(3):215–216. https://doi.org/10.1038/nmeth.1906
2. Ernst J, Kellis M. ChromHMM software and manual (binarization into present/absent marks). http://compbio.mit.edu/ChromHMM/
3. Roadmap Epigenomics Consortium. Chromatin state learning (15-state core and 18-state expanded models). https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
4. Liang G et al. 2004. Distinct localization of histone H3 acetylation and H3-K4 methylation to the transcription start sites in the human genome. PNAS 101(19):7357–7362. https://doi.org/10.1073/pnas.0401866101
5. Rada-Iglesias A. 2018. Is H3K4me1 at enhancers correlative or causative? Nature Genetics 50(1):4–5. https://doi.org/10.1038/s41588-017-0018-3
6. Creyghton MP et al. 2010. Histone H3K27ac separates active from poised enhancers and predicts developmental state. PNAS 107(50):21931–21936. https://doi.org/10.1073/pnas.1016071107
7. Ferrari KJ et al. 2014. Polycomb-dependent H3K27me1 and H3K27me2 regulate active transcription and enhancer fidelity. Molecular Cell 53(1):49–62. https://doi.org/10.1016/j.molcel.2013.10.030
8. Nicetto D et al. 2019. H3K9me3-heterochromatin loss at protein-coding genes enables developmental lineage specification. Science 363(6424):294–297. https://doi.org/10.1126/science.aau0583
9. Kimura H. 2013. Histone modifications for human epigenome analysis. Journal of Human Genetics 58(7):439–445. https://doi.org/10.1038/jhg.2013.66
