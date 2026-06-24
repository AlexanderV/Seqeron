# Test Specification: MIRNA-PRECURSOR-001

**Test Unit ID:** MIRNA-PRECURSOR-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Hairpin Detection
**Status:** ☐ Pending re-validation (MFE-fold opt-in added 2026-06-24)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-24

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Bartel (2004), Cell 116:281-297 | 1 | doi:10.1016/S0092-8674(04)00045-5 | 2026-02-10 |
| 2 | Ambros et al. (2003), RNA 9:277-279 | 1 | doi:10.1261/rna.2183803 | 2026-02-10 |
| 3 | Bartel (2009), Cell 136:215-233 | 1 | doi:10.1016/j.cell.2009.01.002 | 2026-02-10 |
| 4 | Krol et al. (2004), J Biol Chem 279:42230-42239 | 1 | doi:10.1074/jbc.M404931200 | 2026-02-10 |
| 5 | Wikipedia: MicroRNA | 4 | https://en.wikipedia.org/wiki/MicroRNA | 2026-02-10 |
| 6 | miRBase (Griffiths-Jones, 2006) | 5 | https://mirbase.org/ | 2026-02-10 |
| 7 | Bonnet et al. (2004), Bioinformatics 20(17):2911 | 1 | doi:10.1093/bioinformatics/bth374 | 2026-06-24 |
| 8 | Zhang et al. (2006), Cell Mol Life Sci 63:246 (AMFE/MFEI) | 1 | Cell Mol Life Sci 63:246-254 | 2026-06-24 |
| 9 | Meyers/Bartel et al. (2008), Plant Cell 20:3186 | 1 | doi:10.1105/tpc.108.064311 (PMC2630443) | 2026-06-24 |
| 10 | Han et al. (2006), Cell 125:887-901 (Drosha ~11 bp ruler) | 1 | doi:10.1016/j.cell.2006.03.043 (PMID 16751099) | 2026-06-24 |
| 11 | Park et al. (2011), Nature 475:201-205 (Dicer 5' counting ~22 nt) | 1 | doi:10.1038/nature10198 (PMID 21753850) | 2026-06-24 |
| 12 | Auyeung et al. (2013), Cell 152:844-858 (UG/UGU/CNNC motifs) | 1 | doi:10.1016/j.cell.2013.01.031 (PMID 23415231) | 2026-06-24 |
| 13 | miRBase hsa-mir-21 (MI0000077; mature MIMAT0000076) | 5 | https://mirbase.org/hairpin/MI0000077 | 2026-06-24 |

### 1.2 Key Evidence Points

1. Pre-miRNAs are ~60–120 nt hairpin structures — Bartel (2004)
2. Stem length ≥18–22 bp required for Drosha/Dicer processing — Krol (2004)
3. Terminal loop size typically 3–15 nt (up to 25 nt in some cases) — Bartel (2004)
4. Mature miRNA is ~22 nt, from one arm of the hairpin — Bartel (2009)
5. Star (passenger) strand from opposite arm — Wikipedia/MicroRNA (citing Bartel 2004)
6. G:U wobble pairs are valid base pairs in RNA stems — Krol (2004)
7. Dot-bracket notation: '(' for 5' stem bases, ')' for 3' stem bases, '.' for loop — standard RNA notation

**MFE-structure-based (opt-in) evidence:**

8. Pre-miRNA precursors fold to a free energy considerably lower (more stable) than shuffled sequences, unlike tRNA/rRNA — Bonnet et al. (2004).
9. AMFE = 100·MFE/length; MFEI = AMFE/(G+C)% — Zhang et al. (2006).
10. Pre-miRNA MFEI is typically > 0.85, remarkably higher than other RNAs — Zhang et al. (2006). Default acceptance cutoff `minMfei = 0.85`.
11. The mature miRNA sits in one arm of a fold-back hairpin with ≥16 complementary bases to the opposite arm — Ambros et al. (2003); single-arm duplex with ≤4 mismatches / minimal bulges — Meyers/Bartel (2008). ⇒ acceptance requires a single dominant hairpin with stem bp ≥ 16.

**Drosha/Dicer cleavage-site (opt-in) evidence:**

12. **Drosha basal ruler:** "The cleavage site is determined mainly by the distance (approximately 11 bp) from the stem-ssRNA junction" — Han et al. (2006). ⇒ Drosha 5' cut = basal junction + 11 bp (`DroshaCutBpFromBasalJunction = 11`).
13. **Dicer 5' counting rule:** "the cleavage site determined mainly by the distance (∼22 nucleotides) from the 5' end (5' counting rule)" — Park et al. (2011). ⇒ mature length fixed at ~22 nt (`DicerCutNtFrom5PrimeEnd = 22`).
14. **2-nt 3' overhang:** "Cleavage by RNase III domains results in 2-nt 3'-overhang end" — Lee et al. (2003)/Han et al. (2006). ⇒ both Drosha and Dicer leave a 2-nt 3' overhang (`RNaseIII3PrimeOverhang = 2`).
15. **CNNC motif:** "positioned 16–18 nt from the Drosha cut" — Auyeung et al. (2013). ⇒ optional confidence flag only (`HasCnncMotif`).
16. **miRBase cross-check:** hsa-miR-21-5p (MIMAT0000076) = `UAGCUUAUCAGACUGAUGUUGA` (22 nt). With an 11-nt lower stem prepended so the +11 ruler lands at the 5p start, the predicted 5p mature equals this exactly.

### 1.3 Documented Corner Cases

1. **Sequence too short:** Input < minHairpinLength → no candidates possible.
2. **No complementarity:** Random sequence with no self-complementary regions → no hairpin found.
3. **Loop too small/large:** Loop < 3 or > 25 nt → rejected.
4. **Stem too short:** < 18 consecutive base pairs → rejected.
5. **DNA input:** T must be converted to U before analysis.
6. **Empty/null input:** Must return empty, not throw.

### 1.4 Known Failure Modes / Pitfalls

1. **Simplified consecutive-pairing model** misses real pre-miRNAs with internal bulges — implementation limitation (documented in algorithm doc).
2. **Overlapping candidates** — scanning may yield multiple overlapping hairpins from same region.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength)` | MiRnaAnalyzer | Canonical | Public API for pre-miRNA detection (default consecutive-pairing heuristic) |
| `AnalyzeHairpin(sequence, matureLength)` | MiRnaAnalyzer | Internal | Private helper, tested indirectly |
| `AssessHairpinByMfe(candidate, minMfei, minLoopSize)` | MiRnaAnalyzer | Canonical | Opt-in: assess one candidate from its real MFE structure (RNA-STRUCT-001) |
| `FindPreMiRnaHairpinsByMfe(sequence, minHairpinLength, maxHairpinLength, minMfei, minLoopSize)` | MiRnaAnalyzer | Canonical | Opt-in: window scan that folds each candidate with the MFE engine |
| `CalculateMfeIndex(freeEnergy, length, gcPercent)` | MiRnaAnalyzer | Canonical | MFEI = AMFE/(G+C)%, AMFE = 100·\|ΔG°\|/length (Zhang 2006) |
| `PredictDroshaDicerCleavage(sequence, basalJunction)` | MiRnaAnalyzer | Canonical | Opt-in: predicts Drosha (~11 bp from basal junction) + Dicer (~22 nt 5' counting) cuts, mature/star spans, 2-nt 3' overhang, optional CNNC flag (Han 2006 / Park 2011 / Auyeung 2013) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned PreMiRna has Sequence.Length ∈ [minHairpinLength, maxHairpinLength] | Yes | Definition of scanning window |
| INV-2 | Every returned PreMiRna has Start ≥ 0 and End < input.Length | Yes | Array bounds |
| INV-3 | MatureSequence.Length ≤ matureLength and > 0 | Yes | Bartel (2009): ~22 nt mature |
| INV-4 | StarSequence.Length == MatureSequence.Length | Yes | Duplex symmetry |
| INV-5 | Structure.Length == Sequence.Length | Yes | Dot-bracket notation definition |
| INV-6 | Structure contains only '(', ')', and '.' characters | Yes | Standard notation |
| INV-7 | Count of '(' == Count of ')' in Structure | Yes | Balanced base pairs |
| INV-8 | FreeEnergy < 0 for all valid hairpins (stabilizing) | Yes | Turner 2004 nearest-neighbor model (NNDB) |
| INV-9 | A hairpin with a longer effective stem (23 bp / 11 nt loop) has more negative FreeEnergy than one with a shorter stem (20 bp / 7 nt loop); the added paired stem bases dominate the energy despite the differing loop size (see M11) | Yes | Turner (2004) principles |
| INV-10 | Sequence is uppercase RNA (A, U, G, C only) | Yes | T→U conversion in implementation |
| INV-11 | (MFE path) `AssessHairpinByMfe(s).FreeEnergy == RnaSecondaryStructure.CalculateMinimumFreeEnergy(s)` | Yes | MfeStructure.FreeEnergy equals scalar MFE by construction (RNA-STRUCT-001) |
| INV-12 | (MFE path) An accepted `PreMiRnaMfe` has `StemBasePairs ≥ 16`, `TerminalLoopSize ∈ [3,25]`, `Mfei ≥ minMfei` | Yes | Ambros (2003), Bartel (2004), Zhang (2006) |
| INV-13 | (MFE path) `Mfei == CalculateMfeIndex(FreeEnergy, length, GC%)` = (100·\|ΔG°\|/length)/(G+C)% | Yes | Zhang (2006) |
| INV-14 | (Cleavage path) `DroshaCut5Prime == BasalJunction + 11` | Yes | Han et al. (2006) |
| INV-15 | (Cleavage path) `MatureSequence.Length == 22` and `MatureEnd − MatureStart + 1 == 22` | Yes | Park et al. (2011) |
| INV-16 | (Cleavage path) `DroshaCut3Prime − MatureEnd == 2` and `ThreePrimeOverhang == 2` | Yes | Lee (2003)/Han (2006) RNase III 2-nt 3' overhang |
| INV-17 | (Cleavage path) `StarEnd − StarStart + 1 == 22` (3p span length) | Yes | Park et al. (2011) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | NullInput_ReturnsEmpty | Null sequence returns empty enumerable | Empty | Defensive coding |
| M2 | EmptyInput_ReturnsEmpty | Empty string returns empty | Empty | Defensive coding |
| M3 | ShortSequence_ReturnsEmpty | Sequence shorter than minHairpinLength | Empty | Bartel (2004): pre-miRNA ≥ ~55 nt |
| M4 | ValidHairpin_DetectsCandidate | Sequence with known stem-loop (≥18 bp stem, 5-10 nt loop, ≥55 nt total) | ≥1 PreMiRna returned | Bartel (2004), Krol (2004) |
| M5 | Position_Correct | Start/End of returned hairpin match expected position | Start/End within input | Array bounds invariant |
| M6 | MatureSequence_From5PrimeArm | Mature extracted from first matureLength bases of candidate | First N bases of hairpin | Bartel (2009): 5' arm convention |
| M7 | StarSequence_From3PrimeArm | Star extracted from last matureLength bases of candidate | Last N bases of hairpin | Bartel (2009): 3' arm |
| M8 | DotBracket_Structure_Correct | Structure has '(' for stem-5', '.' for loop, ')' for stem-3' | Matches stem/loop/stem pattern | Standard RNA notation |
| M9 | DotBracket_Balanced | Count('(') == Count(')') in structure | Equal counts | Notation definition |
| M10 | FreeEnergy_Negative | Free energy < 0 for valid hairpins | Negative value | Thermodynamic stability |
| M11 | FreeEnergy_Ordering | Longer effective stem (23 bp) more negative than shorter (20 bp) | E(stem=23) < E(stem=20) | Turner (2004) nearest-neighbor model |
| M12 | StemTooShort_Rejected | 55 nt sequence with only 15 bp stem (< 18 required); tests stem rejection, not n<55 | No hairpin found | Krol (2004) |
| M13 | LoopTooLarge_Rejected | 66 nt candidate with 30 nt loop (> 25 max) | Rejected | Bartel (2004) |
| M14 | TtoU_Conversion | DNA input (with T) handled correctly | T converted to U in output | RNA biology standard |
| M15 | GU_WobblePairs_InStem | G-U pairs count as valid stem pairs | Hairpin accepted | Krol (2004) |
| M16 | SequenceLength_InRange | All returned PreMiRnas have length within [min, max] | INV-1 verified | Scanning window definition |
| M17 | Invariants_AllHold | All invariants verified on results | INV-1 through INV-10 | Multiple sources |
| M18 | RealMiRBase_HsaMir21_NotDetected | hsa-mir-21 (MI0000077, 72 nt) — real pre-miRNA not detected by the **heuristic** | Empty (heuristic limitation) | miRBase v22 |
| M19 | RealMiRBase_HsaLet7a1_NotDetected | hsa-let-7a-1 (MI0000060, 80 nt) — real pre-miRNA not detected by the **heuristic** | Empty (heuristic limitation) | miRBase v22 |
| MF1 | AssessHairpinByMfe_PerfectHairpin | `ValidHairpin57` folded by the engine: single hairpin, ΔG°=−48.48, 27 bp, loop 3, AMFE=85.052632, MFEI=1.939200; ΔG° equals `CalculateMinimumFreeEnergy` | Accepted, exact values | RNA-STRUCT-001 engine; Zhang (2006) |
| MF2 | AssessHairpinByMfe_HsaMir21_Detected | hsa-mir-21 folded: single hairpin, ΔG°=−35.13, 32 bp, loop 3, MFEI=1.003714 | **Accepted** (MFE fold detects it) | miRBase v22; RNA-STRUCT-001 |
| MF3 | AssessHairpinByMfe_HsaLet7a1_Detected | hsa-let-7a-1 (80 nt) folded: single hairpin, ΔG°=−34.31, 32 bp, loop 4, MFEI=1.009118 | **Accepted** (MFE fold detects it) | miRBase v22; RNA-STRUCT-001 |
| MF4 | HeuristicRejects_ButMfeDetects | Same hsa-mir-21: heuristic empty, MFE fold accepts | Heuristic ∅, MFE ✓ | Limitation removed |
| MF5 | AssessHairpinByMfe_NoComplementarity_Rejected | `NoComplementarity` (ΔG°=0, all dots) | Rejected (null) | Bonnet (2004) |
| MF6 | AssessHairpinByMfe_Multibranch_RejectedDespiteStrongEnergy | 5S-rRNA-like (120 nt): ΔG°=−47.04 but multibranch | Rejected on structure, not energy | Meyers/Bartel (2008) single-arm criterion |
| MF7 | FindPreMiRnaHairpinsByMfe_DesignedHairpin | Window scan over `ValidHairpin57`: yields candidate with ΔG°=−48.48, 27 bp, MFEI ≥ 0.85 | ≥1 accepted candidate | RNA-STRUCT-001 |
| MF8 | AssessHairpinByMfe_MfeiBelowThreshold_Rejected | hsa-let-7a-1 with `minMfei=1.5` (MFEI 1.009 < 1.5) | Rejected; accepted at 0.85 | Zhang (2006) cutoff |
| MF9 | CalculateMfeIndex_MatchesZhang2006 | MFEI(−48.48, 57, 25/57·100) = 1.939200; zero length / zero GC% ⇒ 0 | Exact + guards | Zhang (2006) |
| MF10 | MfeFoldMethods_NullOrEmpty | null/empty/too-short inputs to both MFE methods | null / empty | Defensive coding |
| DD1 | PredictDroshaDicerCleavage_DroshaCut_Is11BpFromBasalJunction | Synthetic pri-miRNA (11-nt lower stem + miR-21 stem), junction=0 | DroshaCut5Prime == 11 | Han et al. (2006): ~11 bp ruler |
| DD2 | PredictDroshaDicerCleavage_MatureLength_Is22Nt | Same pri-miRNA | MatureSequence.Length == 22; MatureStart == DroshaCut5Prime | Park et al. (2011): 5' counting rule |
| DD3 | PredictDroshaDicerCleavage_ThreePrimeOverhang_Is2Nt | Same pri-miRNA | ThreePrimeOverhang == 2; DroshaCut3Prime − MatureEnd == 2; star length 22 | Lee (2003)/Han (2006) RNase III |
| DD4 | PredictDroshaDicerCleavage_HsaMir21_MatchesMirBaseMature5p | miRBase cross-check | MatureSequence == `UAGCUUAUCAGACUGAUGUUGA` (MIMAT0000076) | miRBase MI0000077 |
| DD5 | PredictDroshaDicerCleavage_StarSpan_HasExpectedCoordinatesAndSequence | Same pri-miRNA | matureStart=11, matureEnd=32, starEnd=34, starStart=13; star = pri[13..34] | Park (2011) + 2-nt overhang |
| DD6 | PredictDroshaDicerCleavage_DnaInput_NormalisesToRna | DNA spelling of the pri-miRNA | Mature == RNA `UAGCUUAUCAGACUGAUGUUGA`, no 'T' | RNA biology standard |
| DD7 | PredictDroshaDicerCleavage_NullOrEmpty_ReturnsNull | null / empty | null | Defensive coding |
| DD8 | PredictDroshaDicerCleavage_JunctionOutOfRange_ReturnsNull | junction < 0 or ≥ length | null | Defensive coding |
| DD9 | PredictDroshaDicerCleavage_TooShortForCuts_ReturnsNull | 30 nt (< 11+22+2 from junction) | null | Geometry guard |
| DD10 | PredictDroshaDicerCleavage_CnncMotif_DetectedWhenPresentDownstream | CNNC 16 nt 3' of Drosha cut vs none | HasCnncMotif true vs false | Auyeung et al. (2013) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | MultipleHairpins_LongSequence | Long sequence with two distinct hairpin regions | ≥2 candidates | Biological: miRNA clusters |
| S2 | MaxHairpinLength_Respected | Candidate exceeding maxHairpinLength filtered | Not returned | Parameter contract |
| S3 | MinHairpinLength_CustomValue | Custom minimum applied correctly | Only candidates ≥ min | Parameter contract |
| S4 | NoComplementarity_ReturnsEmpty | Random non-complementary sequence | Empty | No stem → no hairpin |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | CaseInsensitive_Input | Mixed-case input handled | Same results as uppercase | Robustness |
| C2 | MatureLength_Parameter | Custom matureLength affects MatureSequence.Length | Correct length | Parameterization |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Original weak tests in `MiRnaAnalyzerTests.cs` lines 78–112 (#region Pre-miRNA Tests): 3 methods
- Consolidated into canonical file `MiRnaAnalyzer_PreMiRna_Tests.cs`: 25 tests

### 5.2 Coverage Classification

| Test ID | Status | Notes |
|---------|--------|-------|
| M1 | ✅ Covered | Null input returns empty |
| M2 | ✅ Covered | Empty input returns empty |
| M3 | ✅ Covered | Short sequence rejected |
| M4 | ✅ Covered | Valid hairpin detected |
| M5 | ✅ Covered | Exact Start=0/End=56 + sub-window Start=1/End=55; count=2 |
| M6 | ✅ Covered | Exact mature = first 22 nt ("GCAUAGCUAGCUAGCUAGCUAG") |
| M7 | ✅ Covered | Exact star = last 22 nt ("CUAGCUAGCUAGCUAGCUAUGC") |
| M8 | ✅ Covered | Exact structure: 23×'(' + 11×'.' + 23×')' |
| M9 | ✅ Covered | Balanced parentheses + dot count ≥ 3 |
| M10 | ✅ Covered | Exact Turner energy −4⁣⁣⁣3.50 kcal/mol (hand-calculated from NNDB) |
| M11 | ✅ Covered | Exact energies: stem-23 = −4⁣⁣⁣3.50, stem-20 = −36.02 (hand-calculated) |
| M12 | ✅ Covered | 15 bp stem rejected (stem < 18) |
| M13 | ✅ Covered | 30 nt loop rejected (> 25) |
| M14 | ✅ Covered | DNA T→U conversion verified |
| M15 | ✅ Covered | G:U wobble pairs accepted |
| M16 | ✅ Covered | Sequence length in [min, max] range |
| M17 | ✅ Covered | All 10 invariants verified in composite check |
| M18 | ✅ Covered | hsa-mir-21 not detected by heuristic |
| M19 | ✅ Covered | hsa-let-7a-1 not detected by heuristic |
| MF1 | ✅ Covered | Exact ΔG°=−48.48, 27 bp, loop 3, AMFE 85.052632, MFEI 1.939200; ΔG° == engine scalar MFE |
| MF2 | ✅ Covered | hsa-mir-21 ACCEPTED: ΔG°=−35.13, 32 bp, loop 3, MFEI 1.003714 |
| MF3 | ✅ Covered | hsa-let-7a-1 ACCEPTED: ΔG°=−34.31, 32 bp, loop 4, MFEI 1.009118 |
| MF4 | ✅ Covered | Heuristic ∅ vs MFE ✓ on the same real pre-miRNA |
| MF5 | ✅ Covered | No complementarity (ΔG°=0) rejected |
| MF6 | ✅ Covered | Multibranch (ΔG°=−47.04) rejected on structure |
| MF7 | ✅ Covered | Window scan yields accepted candidate |
| MF8 | ✅ Covered | minMfei gate: rejected at 1.5, accepted at 0.85 |
| MF9 | ✅ Covered | CalculateMfeIndex exact 1.939200 + zero-guards |
| MF10 | ✅ Covered | null/empty/too-short handled |
| DD1 | ✅ Covered | DroshaCut5Prime == 11 (exact) |
| DD2 | ✅ Covered | Mature length == 22 (exact) |
| DD3 | ✅ Covered | 2-nt 3' overhang; star length 22 |
| DD4 | ✅ Covered | Mature == miRBase hsa-miR-21-5p exactly |
| DD5 | ✅ Covered | Exact star coords 11/32/34/13 + sequence |
| DD6 | ✅ Covered | DNA T→U normalisation |
| DD7 | ✅ Covered | null/empty ⇒ null |
| DD8 | ✅ Covered | out-of-range junction ⇒ null |
| DD9 | ✅ Covered | too short ⇒ null |
| DD10 | ✅ Covered | CNNC flag true/false |
| S1 | ✅ Covered | Multiple hairpins in long sequence |
| S2 | ✅ Covered | maxHairpinLength upper bound enforced |
| S3 | ✅ Covered | Custom minHairpinLength applied |
| S4 | ✅ Covered | No complementarity returns empty |
| C1 | ✅ Covered | Mixed-case input produces same results as uppercase |
| C2 | ✅ Covered | Custom matureLength=18 yields exact 18-nt mature/star |

**Summary:** 0 missing, 0 weak, 0 duplicate. 25 heuristic + 10 MFE-fold + 10 cleavage-site = 45 cases covered.

### 5.3 Strengthening Log (2026-03-16)

| Test | Change | Rationale |
|------|--------|-----------|
| M5 | Permissive bounds → exact Start/End + exact candidate count | Removes ambiguity; values derived from scanning algorithm + stem pairing |
| M6 | "length ≤ 22" → exact mature sequence | "GCAUAGCUAGCUAGCUAGCUAG" derived from matureEnd=min(22,23)=22 |
| M7 | "not null" → exact star sequence | "CUAGCUAGCUAGCUAGCUAUGC" derived from starStart=57−22=35 |
| M8 | Regex pattern → exact structure string | 23×'(' + 11×'.' + 23×')' from stem=23, loop=57−46=11 |
| M10 | "< 0" → exact −4⁣⁣⁣3.50 kcal/mol | Hand-calculated from Turner 2004 NNDB: 22 stacking pairs + loop(11) + TM(CUAG) |
| M11 | Ordering only → ordering + exact magnitudes | −4⁣⁣⁣3.50 (stem-23) and −36.02 (stem-20) both hand-calculated |
| C1 | New | Mixed-case → same as uppercase (ToUpperInvariant) |
| C2 | New | matureLength=18 yields exact 18-nt mature/star |

### 5.4 Consolidation Plan

- **Canonical file:** `MiRnaAnalyzer_PreMiRna_Tests.cs` — all MIRNA-PRECURSOR-001 tests
- **Removed from existing file:** Pre-miRNA Tests region from `MiRnaAnalyzerTests.cs` (3 weak tests replaced)
- **Kept in existing file:** All other tests (reverse complement, base pairing, context, family, utility)

### 5.5 Final State

| File | Role | Test Count |
|------|------|------------|
| `MiRnaAnalyzer_PreMiRna_Tests.cs` | Canonical for MIRNA-PRECURSOR-001 | 25 |
| `MiRnaAnalyzerTests.cs` | Residual tests for other units | Same minus 3 |

---

## 6. Design Limitations

| # | Limitation | Impact | Tests |
|---|-----------|--------|-------|
| 1 | **Default heuristic** uses consecutive stem pairing from ends; no tolerance for internal mismatches or bulges. Real pre-miRNAs (e.g., hsa-mir-21) have asymmetric internal loops that offset pairing alignment. | Real miRBase pre-miRNAs are not detected by the **default** model. | M18, M19 |
| 2 | **Resolved (opt-in):** `AssessHairpinByMfe` / `FindPreMiRnaHairpinsByMfe` fold the candidate with the validated Zuker–Stiegler MFE engine (RNA-STRUCT-001) and read the hairpin from the real MFE structure, detecting natural miRBase precursors. | hsa-mir-21 / hsa-let-7a-1 now detected via the MFE path. | MF2, MF3, MF4 |
| 3 | **Resolved (opt-in):** `PredictDroshaDicerCleavage` predicts the Drosha cut (~11 bp from the basal junction, Han 2006), the Dicer cut / mature length (~22 nt, 5' counting rule, Park 2011), the mature (5p) and star (3p) coordinates, the 2-nt 3' overhang, and an optional CNNC confidence flag (Auyeung 2013). Cross-checked against miRBase hsa-miR-21-5p. | Cleavage coordinates now predictable from the published rules. | DD1–DD10 |
| 4 | **Residual:** a competitive **trained** natural-vs-background precursor classifier (e.g. miRDeep2: a fitted probabilistic model using read-stacking signatures) — requires a trained model + labelled data, which are data-blocked. | Decision-grade probabilistic precursor scoring needs a trained model. | — |

---

## 7. Open Questions / Decisions

None. Behaviour is testable via the public APIs `FindPreMiRnaHairpins` (default heuristic) and
`FindPreMiRnaHairpinsByMfe` / `AssessHairpinByMfe` / `CalculateMfeIndex` (opt-in MFE-fold path).
