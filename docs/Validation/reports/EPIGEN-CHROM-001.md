# Validation Report: EPIGEN-CHROM-001 — Chromatin State Prediction

- **Validated:** 2026-06-15   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.PredictChromatinState(...)`, `AnnotateHistoneModifications(...)`, `FindAccessibleRegions(...)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Roadmap Epigenomics — Chromatin state learning** (https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html). Confirmed the 15-state core model uses 5 marks (H3K4me3, H3K4me1, H3K36me3, H3K27me3, H3K9me3) and the 18-state expanded model adds **H3K27ac**. State→mark mapping confirmed: TssA "Active TSS" = H3K4me3; Tx/TxWk = H3K36me3; Enh/EnhG = H3K4me1; Het "Heterochromatin" = H3K9me3; TssBiv "Bivalent/Poised TSS" = H3K4me3 + H3K27me3; EnhBiv "Bivalent Enhancer" = H3K4me1 + H3K27me3; ReprPC/ReprPCWk "Repressed PolyComb" = H3K27me3; Quies "Quiescent/Low" = minimal enrichment. The H3K27ac addition (18-state EnhA vs EnhWk) separates active from weak/poised enhancers.
- **ChromHMM review, Ernst & Kellis, Nat Protoc 2017** (PMC5945550). Confirmed verbatim: *"ChromHMM focuses its modeling power on combinations of epigenomic marks, by using binary presence/absence input features"* and *"ChromHMM then determines the presence or absence of each mark based on the significance of observed count of sequencing reads relative to a Poisson background distribution."* This validates the binarize-then-classify model and INV-01 (state is a function of the present/absent mark set, not magnitude).
- Per-mark biology (carried from Evidence, all standard textbook/primary): H3K4me3→active promoter (Liang 2004); H3K4me1→enhancer (Rada-Iglesias 2018); H3K27ac→active vs poised enhancer (Creyghton 2010); H3K36me3→transcribed gene body (Kimura 2013); H3K27me3→Polycomb repression (Ferrari 2014); H3K9me3→heterochromatin (Nicetto 2019).

### Formula / mapping check
The implementation is a deterministic rule-based realisation of the ChromHMM/Roadmap "signature → state" core (not the learned HMM, which is correctly disclaimed as out of scope). Every row of the §2.2 mark-set→state table matches the Roadmap state definitions I retrieved:
- H3K4me3 → ActivePromoter (TssA) ✓
- H3K4me1+H3K27ac → ActiveEnhancer; H3K4me1 alone → WeakEnhancer (18-state H3K27ac split) ✓
- H3K36me3 → Transcribed (Tx) ✓
- H3K27me3 alone → Repressed (ReprPC) ✓
- H3K9me3 alone → Heterochromatin (Het) ✓
- H3K4me3+H3K27me3 → BivalentPromoter (TssBiv) ✓
- H3K4me1+H3K27me3 → BivalentEnhancer (EnhBiv) ✓
- none → LowSignal (Quies) ✓

### Edge-case semantics check
- No mark present ⇒ LowSignal (Quies) — sourced (INV-02). ✓
- Bivalent co-occurrence checked before plain active/repressed (INV-03) — the *combination* defines the state per Roadmap TssBiv/EnhBiv; correct ordering. ✓
- Promoter-over-enhancer precedence (H3K4me3+H3K4me1 → ActivePromoter): recorded as an explicit Assumption; consistent with Roadmap TSS-above-Enh mnemonic ordering. Acceptable as documented.
- Binarization threshold value (default 0.5 on normalized [0,1]): correctly recorded as an Assumption — ChromHMM uses a Poisson background, not a fixed numeric cut; tests use unambiguous magnitudes so the result is threshold-independent.

### Independent cross-check (numbers)
Hand-derived expected states from the retrieved Roadmap table, all reproduced by the code:
| Present marks | Expected (Roadmap) | Code |
|---|---|---|
| H3K4me3 | TssA / ActivePromoter | ActivePromoter ✓ |
| H3K4me1+H3K27ac | active Enh / ActiveEnhancer | ActiveEnhancer ✓ |
| H3K4me1 | Enh / WeakEnhancer | WeakEnhancer ✓ |
| H3K36me3 | Tx / Transcribed | Transcribed ✓ |
| H3K27me3 | ReprPC / Repressed | Repressed ✓ |
| H3K9me3 | Het / Heterochromatin | Heterochromatin ✓ |
| H3K4me3+H3K27me3 | TssBiv / BivalentPromoter | BivalentPromoter ✓ |
| H3K4me1+H3K27me3 | EnhBiv / BivalentEnhancer | BivalentEnhancer ✓ |
| none | Quies / LowSignal | LowSignal ✓ |

### Findings / divergences
None. Description and per-mark biology are authoritative and correctly cited. The two recorded Assumptions (threshold value; promoter precedence) are honest and sourced as far as the literature allows.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`:
- `PredictChromatinState` (lines 883–933) — binarize six marks at `presenceThreshold` (inclusive `>=`), bivalent checks first, then promoter→enhancer→transcribed→repressed→heterochromatin→LowSignal.
- `AnnotateHistoneModifications` + `InferStateFromMark` (947–982) — per-region single-mark → state, case-insensitive, unknown mark → LowSignal.
- `FindAccessibleRegions` (991–1066) — single position-sorted scan, gap-merge ≤ maxGap, minWidth filter, max-signal score, cosmetic PeakType.

### Formula realised correctly?
Yes. Binarization is `signal >= presenceThreshold` (matches the inclusive convention in the doc). Bivalent precedence is checked before plain active states, so the *combination* governs (INV-03). Promoter ranks above enhancer (INV via Roadmap ordering). The exhaustive fall-through to LowSignal makes it a total function (INV-04). FindAccessibleRegions enforces `End >= Start` and only emits regions whose samples meet the threshold (INV-05).

### Cross-verification vs code
All nine Roadmap signatures above reproduced exactly by the running tests. New tests added this session also confirm the two-active-marks-plus-repressive branches (see below).

### Variant/delegate consistency
`AnnotateHistoneModifications` single-mark mapping agrees with `PredictChromatinState` single-mark results (H3K4me3→ActivePromoter, H3K27me3→Repressed, etc.). H3K4me1 alone maps to WeakEnhancer in both. Consistent.

### Test quality audit (HARD gate)
Pre-existing file `EpigeneticsAnalyzer_ChromatinState_Tests.cs` (24 tests) uses **exact** `Is.EqualTo` assertions against sourced Roadmap states — no `Greater/AtLeast/Contains`, no widened tolerances, no skips. Expectations trace to the Roadmap state table, not to code output. Gate result: **PASS**, after closing three genuine coverage gaps I identified ("cover all the logic"):
- `PredictChromatinState_K4me3AndK4me1AndK27me3_BivalentPromoter` — exercises the BivalentPromoter-over-BivalentEnhancer branch when both active marks co-occur with H3K27me3 (Roadmap TssBiv > EnhBiv).
- `PredictChromatinState_K4me1AndK27acAndK27me3_BivalentEnhancer` — confirms the H3K27ac active call does NOT override the repressive co-occurrence (Roadmap EnhBiv).
- `PredictChromatinState_BinaryInvariance_MultiMarkSignature` — INV-01 for a two-mark signature (H3K4me1+H3K27ac at 0.51 vs 0.99 → identical ActiveEnhancer).

All three pass against the current (correct) implementation — they are coverage additions, not bug fixes. No assertion was weakened; no expected value was adjusted to match output.

### Findings / defects
No defect. The implementation faithfully realises the validated Stage-A description.

## Verdict & follow-ups
- **Stage A: PASS.** Description, formulae, edge cases, and citations independently confirmed against the Roadmap chromatin-state portal and the Ernst & Kellis ChromHMM review.
- **Stage B: PASS.** Code matches the validated description on all branches; tests assert exact sourced values; three coverage gaps closed.
- **End-state: ✅ CLEAN.** Full unfiltered suite: 6542 passed, 0 failed; build 0 warnings/0 errors.
- No findings logged in the register (no defect).
