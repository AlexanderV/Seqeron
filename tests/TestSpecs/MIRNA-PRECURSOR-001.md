# Test Specification: MIRNA-PRECURSOR-001

**Test Unit ID:** MIRNA-PRECURSOR-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Hairpin Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-16

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

### 1.2 Key Evidence Points

1. Pre-miRNAs are ~60–120 nt hairpin structures — Bartel (2004)
2. Stem length ≥18–22 bp required for Drosha/Dicer processing — Krol (2004)
3. Terminal loop size typically 3–15 nt (up to 25 nt in some cases) — Bartel (2004)
4. Mature miRNA is ~22 nt, from one arm of the hairpin — Bartel (2009)
5. Star (passenger) strand from opposite arm — Wikipedia/MicroRNA (citing Bartel 2004)
6. G:U wobble pairs are valid base pairs in RNA stems — Krol (2004)
7. Dot-bracket notation: '(' for 5' stem bases, ')' for 3' stem bases, '.' for loop — standard RNA notation

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
| `FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength)` | MiRnaAnalyzer | Canonical | Public API for pre-miRNA detection |
| `AnalyzeHairpin(sequence, matureLength)` | MiRnaAnalyzer | Internal | Private helper, tested indirectly |

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
| INV-9 | Longer stem → more negative FreeEnergy (all else equal) | Yes | Turner (2004) principles |
| INV-10 | Sequence is uppercase RNA (A, U, G, C only) | Yes | T→U conversion in implementation |

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
| M18 | RealMiRBase_HsaMir21_NotDetected | hsa-mir-21 (MI0000077, 72 nt) — real pre-miRNA not detected | Empty (known limitation) | miRBase v22 |
| M19 | RealMiRBase_HsaLet7a1_NotDetected | hsa-let-7a-1 (MI0000060, 80 nt) — real pre-miRNA not detected | Empty (known limitation) | miRBase v22 |

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
| M18 | ✅ Covered | hsa-mir-21 not detected (known limitation) |
| M19 | ✅ Covered | hsa-let-7a-1 not detected (known limitation) |
| S1 | ✅ Covered | Multiple hairpins in long sequence |
| S2 | ✅ Covered | maxHairpinLength upper bound enforced |
| S3 | ✅ Covered | Custom minHairpinLength applied |
| S4 | ✅ Covered | No complementarity returns empty |
| C1 | ✅ Covered | Mixed-case input produces same results as uppercase |
| C2 | ✅ Covered | Custom matureLength=18 yields exact 18-nt mature/star |

**Summary:** 0 missing, 0 weak, 0 duplicate. All 25 tests covered.

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
| 1 | Consecutive stem pairing from ends; no tolerance for internal mismatches or bulges. Real pre-miRNAs (e.g., hsa-mir-21) have asymmetric internal loops that offset pairing alignment. A full RNA folding algorithm (Zuker/Nussinov) would be required to detect them. | Real miRBase pre-miRNAs are not detected by this model. | M18, M19 |

---

## 7. Open Questions / Decisions

None. All behavior is testable via the public API `FindPreMiRnaHairpins`.
