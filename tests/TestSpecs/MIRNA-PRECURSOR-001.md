# Test Specification: MIRNA-PRECURSOR-001

**Test Unit ID:** MIRNA-PRECURSOR-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Hairpin Detection
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-10

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
| INV-8 | FreeEnergy < 0 for all valid hairpins (stabilizing) | Yes | Thermodynamic requirement — **ASSUMPTION** (simplified model) |
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
| M11 | FreeEnergy_Ordering | Longer stem produces more negative energy | E(stem=25) < E(stem=20) | **ASSUMPTION** (simplified model) |
| M12 | StemTooShort_Rejected | Sequence with < 18 bp complementary stem | No hairpin found | Krol (2004) |
| M13 | LoopTooSmall_Rejected | Candidate with < 3 nt loop | Rejected | Bartel (2004) |
| M14 | TtoU_Conversion | DNA input (with T) handled correctly | T converted to U in output | RNA biology standard |
| M15 | GU_WobblePairs_InStem | G-U pairs count as valid stem pairs | Hairpin accepted | Krol (2004) |
| M16 | SequenceLength_InRange | All returned PreMiRnas have length within [min, max] | INV-1 verified | Scanning window definition |
| M17 | Invariants_AllHold | All invariants verified on results | INV-1 through INV-10 | Multiple sources |

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

- Existing tests in `MiRnaAnalyzerTests.cs` lines 78–112 (#region Pre-miRNA Tests)
- 3 test methods found

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1: NullInput | ❌ Missing | Not tested |
| M2: EmptyInput | ❌ Missing | Not tested |
| M3: ShortSequence | ⚠ Weak | `FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty` exists but uses arbitrary input, no evidence citation |
| M4: ValidHairpin | ⚠ Weak | `FindPreMiRnaHairpins_ValidHairpin_FindsPreMiRNA` has assertion `Has.Count.GreaterThanOrEqualTo(0)` — always passes |
| M5-M17 | ❌ Missing | No position, structure, energy, or invariant tests |
| S1-S4 | ❌ Missing | No edge case tests |
| Structure test | ⚠ Weak | `FindPreMiRnaHairpins_ReturnsStructureInfo` checks `Is.Not.Empty` — no structural validation |

### 5.3 Consolidation Plan

- **Canonical file:** `MiRnaAnalyzer_PreMiRna_Tests.cs` — new file for all MIRNA-PRECURSOR-001 tests
- **Remove from existing file:** Pre-miRNA Tests region from `MiRnaAnalyzerTests.cs` (3 tests, all weak/duplicate)
- **Keep in existing file:** All other tests (reverse complement, base pairing, context, family, utility) — belong to other test units

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MiRnaAnalyzer_PreMiRna_Tests.cs` | Canonical for MIRNA-PRECURSOR-001 | ~17 |
| `MiRnaAnalyzerTests.cs` | Residual tests for other units | Same minus 3 |

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Simplified energy model; test relative ordering not absolute values | M10, M11, INV-8, INV-9 |
| 2 | Consecutive stem pairing; no bulge tolerance | M4, M12 |
| 3 | ValidateHairpin = AnalyzeHairpin (private); tested indirectly | Methods table |

---

## 7. Open Questions / Decisions

None. All behavior is testable via the public API `FindPreMiRnaHairpins`.
