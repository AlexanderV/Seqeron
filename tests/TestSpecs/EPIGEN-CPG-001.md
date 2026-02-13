# Test Specification: EPIGEN-CPG-001

**Test Unit ID:** EPIGEN-CPG-001
**Area:** Epigenetics
**Algorithm:** CpG Site Detection
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Gardiner-Garden M, Frommer M (1987). "CpG islands in vertebrate genomes." J Mol Biol. 196(2):261–282. | 1 | https://doi.org/10.1016/0022-2836(87)90689-9 | 2026-02-13 |
| 2 | Takai D, Jones PA (2002). "Comprehensive analysis of CpG islands in human chromosomes 21 and 22." PNAS. 99(6):3740–5. | 1 | https://doi.org/10.1073/pnas.052410099 | 2026-02-13 |
| 3 | Saxonov S, Berg P, Brutlag DL (2006). "A genome-wide analysis of CpG dinucleotides." PNAS. 103(5):1412–1417. | 1 | https://doi.org/10.1073/pnas.0510310103 | 2026-02-13 |
| 4 | Wikipedia. "CpG site." | 4 | https://en.wikipedia.org/wiki/CpG_site | 2026-02-13 |

### 1.2 Key Evidence Points

1. CpG = cytosine immediately followed by guanine in 5'→3' direction — Wikipedia citing CpG definition.
2. CpG O/E = CpG_count / (C_count × G_count / Length) — Gardiner-Garden & Frommer (1987).
3. CpG island: ≥200 bp, GC% > 50%, O/E > 0.6 — Gardiner-Garden & Frommer (1987).
4. Stricter criteria: ≥500 bp, GC% > 55%, O/E > 0.65 — Takai & Jones (2002).
5. CpG is NOT GpC — Wikipedia CpG definition.

### 1.3 Documented Corner Cases

1. Sequence with no C or G: O/E = 0 (expected = 0, guarded).
2. Sequence shorter than minimum island length: no CpG island possible.
3. Adjacent CpG: "CGCG" contains 2 CpG sites at positions 0 and 2.
4. Single nucleotide: no dinucleotide possible → 0 sites.

### 1.4 Known Failure Modes / Pitfalls

1. Confusing GpC with CpG — must scan only C-then-G in 5'→3' — CpG definition.
2. Case sensitivity — sequences may be lowercase; must normalize — standard convention.
3. Division by zero in O/E when C=0 or G=0 — Gardiner-Garden formula.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindCpGSites(sequence)` | EpigeneticsAnalyzer | Canonical | Returns 0-based positions of CpG dinucleotides |
| `FindCpGIslands(sequence, minLength, minGc, minCpGRatio)` | EpigeneticsAnalyzer | Canonical | Sliding window CpG island detection |
| `CalculateCpGObservedExpected(sequence)` | EpigeneticsAnalyzer | Canonical | Gardiner-Garden & Frommer O/E ratio |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned CpG position i satisfies sequence[i]='C' and sequence[i+1]='G' | Yes | CpG definition |
| INV-2 | O/E ratio = CpG_count / (C_count × G_count / Length) when expected > 0 | Yes | Gardiner-Garden & Frommer (1987) |
| INV-3 | O/E ratio = 0 when C=0 or G=0 or sequence is empty/null | Yes | Formula guard |
| INV-4 | Every CpG island has length ≥ minLength, GC% ≥ minGc, O/E ≥ minCpGRatio | Yes | Gardiner-Garden & Frommer (1987) |
| INV-5 | FindCpGSites count ≤ Length - 1 (maximum possible CpG count) | Yes | Dinucleotide definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FindCpGSites exact positions | "ACGTCGACG" → positions [1, 4, 7] | Exactly 3 CpG at stated positions | CpG definition |
| M2 | FindCpGSites adjacent CpG | "CGCGCG" → positions [0, 2, 4] | Exactly 3 CpG sites | CpG definition |
| M3 | FindCpGSites no CpG | "AATTAATT" → empty | 0 sites | CpG definition |
| M4 | FindCpGSites null/empty | null, "" → empty | 0 sites | Trivial correctness |
| M5 | FindCpGSites case insensitive | "acgtcg" → positions [1, 4] | Same as uppercase | Standard convention |
| M6 | FindCpGSites boundary start | "CGAAA" → position [0] | CpG at start of sequence | CpG definition |
| M7 | FindCpGSites boundary end | "AACG" → position [2] | CpG at end of sequence | CpG definition |
| M8 | FindCpGSites single char | "C" → empty | Cannot form dinucleotide | CpG definition |
| M9 | O/E CGCG repeat | "CGCGCGCGCGCGCGCGCGCG" → 2.0 | 10 CpG, expected=5.0, O/E=2.0 | Gardiner-Garden formula |
| M10 | O/E AT-only | "AATTAATTAATTAATTAATT" → 0.0 | No CpG, no C or G | Gardiner-Garden formula |
| M11 | O/E mixed | "ACGTCGACG" → 3.0 | 3 CpG, C=3, G=3, exp=1.0 | Gardiner-Garden formula |
| M12 | O/E minimal | "ACGT" → 4.0 | 1 CpG, C=1, G=1, exp=0.25 | Gardiner-Garden formula |
| M13 | O/E null/empty | null → 0.0 | Guard against invalid input | Trivial correctness |
| M14 | O/E single char | "A" → 0.0 | Length < 2 | Gardiner-Garden formula |
| M15 | FindCpGIslands positive | 400 bp CGCG repeat → 1 island | All criteria met: L≥200, GC=100%, O/E=2.0 | Gardiner-Garden criteria |
| M16 | FindCpGIslands no island | AT-only → empty | Fails GC% and O/E criteria | Gardiner-Garden criteria |
| M17 | FindCpGIslands too short | "CGCGCG" with minLength=200 → empty | Length < minLength | Gardiner-Garden criteria |
| M18 | GpC not counted as CpG | "GCGCGC" → CpG at [1, 3] only | GC at 0 is GpC, not CpG | CpG definition |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindCpGSites minimal CG | "CG" → position [0] | Minimal valid CpG | Edge case |
| S2 | O/E case insensitive | "cgcgcgcgcgcgcgcgcgcg" → 2.0 | Same result as uppercase | Convention |
| S3 | FindCpGIslands GC content check | Island must have GcContent and CpGRatio in result tuple | Verify returned values | Structural |
| S4 | FindCpGIslands custom params | Custom minLength=100, minGc=0.4 | Respects parameter overrides | API correctness |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | O/E with only C, no G | "CCCC" → 0.0 | Expected = 0, guarded | Edge case |
| C2 | FindCpGIslands null input | null → empty | Guard | Robustness |

---

## 5. Test Inventory

**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CpGDetection_Tests.cs`
**Total tests:** 25
**Coverage:** ✅ Full — 0 missing, 0 weak, 0 duplicates

| ID | Test Method | Assertion Type |
|----|-------------|---------------|
| M1 | `FindCpGSites_SimpleCpG_ReturnsExactPositions` | Exact count + exact positions with messages |
| M2 | `FindCpGSites_AdjacentCpG_ReturnsAllPositions` | Exact count + exact positions |
| M3 | `FindCpGSites_NoCpG_ReturnsEmpty` | `Is.Empty` |
| M4a | `FindCpGSites_NullSequence_ReturnsEmpty` | `Is.Empty` |
| M4b | `FindCpGSites_EmptySequence_ReturnsEmpty` | `Is.Empty` |
| M5 | `FindCpGSites_LowercaseSequence_ReturnsExactPositions` | Exact count + exact positions |
| M6 | `FindCpGSites_CpGAtStart_Detected` | Exact position 0 |
| M7 | `FindCpGSites_CpGAtEnd_Detected` | Exact position |
| M8 | `FindCpGSites_SingleNucleotide_ReturnsEmpty` | `Is.Empty` |
| M9 | `CalculateCpGObservedExpected_PureCGRepeat_ReturnsExact2` | `Is.EqualTo(2.0)` with formula derivation |
| M10 | `CalculateCpGObservedExpected_ATOnly_ReturnsZero` | `Is.EqualTo(0.0)` |
| M11 | `CalculateCpGObservedExpected_MixedSequence_ReturnsExact3` | `Is.EqualTo(3.0)` with formula derivation |
| M12 | `CalculateCpGObservedExpected_MinimalSequence_ReturnsExact4` | `Is.EqualTo(4.0)` with formula derivation |
| M13 | `CalculateCpGObservedExpected_NullSequence_ReturnsZero` | `Is.EqualTo(0.0)` |
| M14 | `CalculateCpGObservedExpected_SingleChar_ReturnsZero` | `Is.EqualTo(0.0)` |
| M15 | `FindCpGIslands_CpGRichRegion_DetectsIsland` | Exact GcContent + CpGRatio assertions |
| M16 | `FindCpGIslands_ATRichSequence_ReturnsEmpty` | `Is.Empty` |
| M17 | `FindCpGIslands_ShortSequence_ReturnsEmpty` | `Is.Empty` |
| M18 | `FindCpGSites_GpCNotCountedAsCpG` | Exact count + exact positions |
| S1 | `FindCpGSites_MinimalCG_ReturnsOnePosition` | Exact position 0 |
| S2 | `CalculateCpGObservedExpected_Lowercase_ReturnsSameAsUppercase` | `Is.EqualTo(2.0)` both cases |
| S3 | `FindCpGIslands_ResultContainsValidMetrics` | Start/End/GcContent/CpGRatio structural |
| S4 | `FindCpGIslands_CustomParameters_Respected` | Default finds, impossible threshold filters |
| C1 | `CalculateCpGObservedExpected_OnlyCNoG_ReturnsZero` | `Is.EqualTo(0.0)` |
| C2 | `FindCpGIslands_NullSequence_ReturnsEmpty` | `Is.Empty` |

---

## 6. Assumption Register

None.

---

## 7. Open Questions / Decisions

None.
