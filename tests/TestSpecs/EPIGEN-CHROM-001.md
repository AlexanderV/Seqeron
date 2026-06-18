# Test Specification: EPIGEN-CHROM-001

**Test Unit ID:** EPIGEN-CHROM-001
**Area:** Epigenetics
**Algorithm:** Chromatin State Prediction from histone modification marks
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Ernst & Kellis (2012). ChromHMM. Nature Methods 9(3):215–216 | 1 | https://doi.org/10.1038/nmeth.1906 | 2026-06-13 |
| 2 | ChromHMM software/manual (present/absent binarization) | 3 | http://compbio.mit.edu/ChromHMM/ | 2026-06-13 |
| 3 | Roadmap Epigenomics chromatin state learning (15/18-state) | 2 | https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html | 2026-06-13 |
| 4 | Liang et al. (2004). PNAS 101(19):7357–7362 (H3K4me3 → active promoter) | 1 | https://doi.org/10.1073/pnas.0401866101 | 2026-06-13 |
| 5 | Rada-Iglesias (2018). Nat Genet 50(1):4–5 (H3K4me1 → enhancer) | 1 | https://doi.org/10.1038/s41588-017-0018-3 | 2026-06-13 |
| 6 | Creyghton et al. (2010). PNAS 107(50):21931–21936 (H3K27ac → active enhancer) | 1 | https://doi.org/10.1073/pnas.1016071107 | 2026-06-13 |
| 7 | Ferrari et al. (2014). Mol Cell 53(1):49–62 (H3K27me3 → Polycomb) | 1 | https://doi.org/10.1016/j.molcel.2013.10.030 | 2026-06-13 |
| 8 | Nicetto et al. (2019). Science 363(6424):294–297 (H3K9me3 → heterochromatin) | 1 | https://doi.org/10.1126/science.aau0583 | 2026-06-13 |

### 1.2 Key Evidence Points

1. ChromHMM models marks as **present/absent** (binary); state is a function of the *set* of present marks — Source 1, 2.
2. Canonical mark → state signatures (Roadmap): H3K4me3→TssA (active promoter); H3K4me1→Enh (enhancer); H3K4me1+H3K27ac→active enhancer; H3K36me3→Tx (transcribed); H3K27me3→ReprPC (Polycomb repressed); H3K9me3→Het (heterochromatin); H3K4me3+H3K27me3→TssBiv (bivalent promoter); H3K4me1+H3K27me3→EnhBiv (bivalent enhancer); none→Quies (low) — Source 3.
3. Per-mark function confirmed by primaries — Sources 4–8.

### 1.3 Documented Corner Cases

- No mark present → Quiescent/Low (Roadmap Quies). — Source 3.
- H3K4me3 + H3K27me3 co-occurrence → bivalent/poised TSS (TssBiv), not active or repressed alone. — Source 3.
- Promoter signature (H3K4me3) takes precedence over enhancer signature (H3K4me1) at the same locus. — Source 3 (TSS ranks above Enh).

### 1.4 Known Failure Modes / Pitfalls

1. Treating signal magnitude as ordinal beyond the presence call: ChromHMM binarizes, so once a mark is present its exact magnitude must not change the state. — Source 1, 2.
2. Misclassifying bivalent (H3K4me3+H3K27me3) as purely active or purely repressed. — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictChromatinState(h3k4me3, h3k4me1, h3k27ac, h3k36me3, h3k27me3, h3k9me3, presenceThreshold)` | EpigeneticsAnalyzer | Canonical | Binary mark signature → ChromatinState |
| `AnnotateHistoneModifications(modifications)` | EpigeneticsAnalyzer | Delegate | Per-region single-mark → state label |
| `FindAccessibleRegions(signal, threshold, minWidth, maxGap)` | EpigeneticsAnalyzer | Delegate | Peak calling over an accessibility track |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | State is a function of the present/absent mark set: two inputs with identical above/below-threshold patterns yield the same state | Yes | Source 1, 2 (binary model) |
| INV-2 | No mark present ⇒ `LowSignal` | Yes | Source 3 (Quies) |
| INV-3 | H3K4me3 present and H3K27me3 present ⇒ `BivalentPromoter` (overrides plain active/repressed) | Yes | Source 3 (TssBiv) |
| INV-4 | Return value is always a defined `ChromatinState` enum member (total function) | Yes | Implementation contract |
| INV-5 | `FindAccessibleRegions` returns regions with `End >= Start` and `AccessibilityScore >= threshold` | Yes | Peak-calling contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | ActivePromoter | H3K4me3 present (+H3K27ac) | `ActivePromoter` | Source 3 TssA; Source 4 |
| M2 | ActivePromoter H3K4me3 only | H3K4me3 present alone | `ActivePromoter` | Source 3 TssA; Source 4 |
| M3 | ActiveEnhancer | H3K4me1 + H3K27ac present | `ActiveEnhancer` | Source 3; Source 5,6 |
| M4 | WeakEnhancer | H3K4me1 present, H3K27ac absent | `WeakEnhancer` | Source 3 Enh; Source 5,6 |
| M5 | Transcribed | H3K36me3 present | `Transcribed` | Source 3 Tx |
| M6 | Repressed | H3K27me3 present alone | `Repressed` | Source 3 ReprPC; Source 7 |
| M7 | Heterochromatin | H3K9me3 present alone | `Heterochromatin` | Source 3 Het; Source 8 |
| M8 | BivalentPromoter | H3K4me3 + H3K27me3 present | `BivalentPromoter` | Source 3 TssBiv |
| M9 | BivalentEnhancer | H3K4me1 + H3K27me3 present | `BivalentEnhancer` | Source 3 EnhBiv |
| M10 | LowSignal | No mark present | `LowSignal` | Source 3 Quies (INV-2) |
| M11 | Binary invariance | Same present/absent pattern, magnitudes 0.6 vs 0.99 | identical state | Source 1,2 (INV-1) |
| M12 | Promoter precedence | H3K4me3 + H3K4me1 both present | `ActivePromoter` | Source 3 (TSS > Enh) |
| M13 | AnnotateHistoneModifications | One region per mark | each region labeled with its mark's state | Source 3 |
| M14 | FindAccessibleRegions merge | contiguous above-threshold positions | single region End>=Start, score>=threshold | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Threshold boundary | mark exactly at threshold counts as present (>=) | present | `>=` semantics |
| S2 | Negative/zero signal | all marks 0 or negative | `LowSignal` | absent |
| S3 | FindAccessibleRegions empty | empty signal | empty result | no regions |
| S4 | AnnotateHistoneModifications empty | empty input | empty result | delegation |
| S5 | FindAccessibleRegions sub-minWidth | region narrower than minWidth | excluded | width filter |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Custom threshold | high presenceThreshold suppresses a marginal mark | downgraded state | parameterization |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `PredictChromatinState`, `AnnotateHistoneModifications`, or `FindAccessibleRegions`. `EpigeneticsAnalyzer_Bisulfite_Tests.cs` and `EpigeneticsAnalyzer_CpGDetection_Tests.cs` cover other methods only.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14 | ❌ Missing | New unit, no prior tests |
| S1–S5 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_ChromatinState_Tests.cs` — all cases for the three methods.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| EpigeneticsAnalyzer_ChromatinState_Tests.cs | Canonical (this unit) | 20 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented (property/invariance) | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | S3 | ❌ Missing | Implemented | ✅ Done |
| 18 | S4 | ❌ Missing | Implemented | ✅ Done |
| 19 | S5 | ❌ Missing | Implemented | ✅ Done |
| 20 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | PredictChromatinState_ActivePromoter |
| M2 | ✅ | PredictChromatinState_H3K4me3Only_ActivePromoter |
| M3 | ✅ | PredictChromatinState_ActiveEnhancer |
| M4 | ✅ | PredictChromatinState_WeakEnhancer |
| M5 | ✅ | PredictChromatinState_Transcribed |
| M6 | ✅ | PredictChromatinState_Repressed |
| M7 | ✅ | PredictChromatinState_Heterochromatin |
| M8 | ✅ | PredictChromatinState_BivalentPromoter |
| M9 | ✅ | PredictChromatinState_BivalentEnhancer |
| M10 | ✅ | PredictChromatinState_NoMarks_LowSignal |
| M11 | ✅ | PredictChromatinState_BinaryInvariance |
| M12 | ✅ | PredictChromatinState_PromoterPrecedence |
| M13 | ✅ | AnnotateHistoneModifications_PerMark |
| M14 | ✅ | FindAccessibleRegions_Merge |
| S1 | ✅ | PredictChromatinState_ThresholdBoundary |
| S2 | ✅ | PredictChromatinState_NegativeSignal_LowSignal |
| S3 | ✅ | FindAccessibleRegions_Empty |
| S4 | ✅ | AnnotateHistoneModifications_Empty |
| S5 | ✅ | FindAccessibleRegions_SubMinWidth |
| C1 | ✅ | PredictChromatinState_HighThreshold |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Presence-call threshold value (default 0.5 on normalized [0,1]); state logic given present/absent pattern is source-backed | PredictChromatinState parameter |
| 2 | Promoter (H3K4me3) precedence over enhancer (H3K4me1) at one locus | M12 |

---

## 7. Open Questions / Decisions

1. The exact ChromHMM emission/transition probabilities are learned per-dataset and are not reproducible as a fixed pure function; this unit implements the canonical *signature → state* mapping over binarized marks (the deterministic, source-defined core), with the binarization threshold exposed as a parameter. Recorded as Assumption 1.
