# Test Specification: MOTIF-REGULATORY-001

**Test Unit ID:** MOTIF-REGULATORY-001
**Area:** Matching
**Algorithm:** Regulatory Elements (scan a DNA sequence for known regulatory consensus motifs)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Bucher (1990) J Mol Biol — eukaryotic Pol II promoter weight matrices (TATA, CCAAT) | 1 | https://doi.org/10.1016/0022-2836(90)90223-9 | 2026-06-14 |
| 2 | Harley & Reynolds (1987) Nucleic Acids Res — E. coli -10/-35 hexamers | 1 | https://doi.org/10.1093/nar/15.5.2343 | 2026-06-14 |
| 3 | Lundin, Nehlin & Ronne (1994) Mol Cell Biol — GC box GGGCGG | 1 | https://doi.org/10.1128/mcb.14.3.1979 | 2026-06-14 |
| 4 | Kozak (1987) Nucleic Acids Res — Kozak gccRccATGG | 1 | https://doi.org/10.1093/nar/15.20.8125 | 2026-06-14 |
| 5 | Proudfoot & Brownlee (1976) Nature — poly(A) AATAAA | 1 | https://doi.org/10.1038/263211a0 | 2026-06-14 |
| 6 | Massari & Murre (2000) Mol Cell Biol — E-box CANNTG | 1 | https://doi.org/10.1128/MCB.20.2.429-440.2000 | 2026-06-14 |
| 7 | Lee, Mitchell & Tjian (1987) Cell — AP-1 TGACTCA | 1 | https://pubmed.ncbi.nlm.nih.gov/3034433/ | 2026-06-14 |
| 8 | Sen & Baltimore (1986) Cell — NF-κB κB site | 1 | https://doi.org/10.1016/0092-8674(86)90346-6 | 2026-06-14 |
| 9 | Montminy et al. (1986) PNAS — CREB CRE TGACGTCA | 1 | https://doi.org/10.1073/pnas.83.18.6682 | 2026-06-14 |
| 10 | Shine–Dalgarno (Wikipedia, citing primaries) — AGGAGG | 4 | https://en.wikipedia.org/wiki/Shine%E2%80%93Dalgarno_sequence | 2026-06-14 |

### 1.2 Key Evidence Points

1. TATA box consensus = `TATAAA` — Bucher (1990) / Wikipedia TATA box.
2. -10 (Pribnow) box = `TATAAT`, -35 box = `TTGACA` — Harley & Reynolds (1987).
3. CCAAT box pentanucleotide = `CCAAT` (~30% of promoters) — Bucher (1990).
4. GC box = `GGGCGG` — Lundin, Nehlin & Ronne (1994).
5. Kozak most-preferred string = `GCCGCCACCATGG`; -3 purine, +4 G — Kozak (1987).
6. Shine-Dalgarno = `AGGAGG`, complementary to 3' 16S rRNA — Wikipedia citing primaries.
7. Poly(A) signal = `AATAAA` — Proudfoot & Brownlee (1976).
8. E-box = `CANNTG` (IUPAC N), canonical palindrome `CACGTG` — Massari & Murre (2000).
9. AP-1 (TRE) motif = `TGACTCA` (NOT TGAGTCA) — Lee, Mitchell & Tjian (1987).
10. NF-κB reference κB site = `GGGACTTTCC` (consensus GGGRNWYYCC) — Sen & Baltimore (1986).
11. CREB CRE palindrome = `TGACGTCA` — Montminy et al. (1986).

### 1.3 Documented Corner Cases

- E-box `CANNTG` is degenerate (IUPAC N): any base matches the two central positions.
- Multiple occurrences of one element and occurrences of several distinct elements may co-exist; each is reported with a 0-based start position.
- Empty sequence → no occurrences. Null sequence → `ArgumentNullException`.

### 1.4 Known Failure Modes / Pitfalls

1. Using `TGAGTCA` for AP-1 (G at position 4) is wrong; the consensus is `TGACTCA` — Lee, Mitchell & Tjian (1987). (Defect corrected in this unit.)
2. Confusing eukaryotic TATA (`TATAAA`) with the prokaryotic -10 hexamer (`TATAAT`) — distinct elements, differ only in the last base — Harley & Reynolds (1987).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindRegulatoryElements(DnaSequence)` | MotifFinder | **Canonical** | Scans for the 12 known regulatory consensus motifs. |
| `KnownMotifs` (consensus constants) | MotifFinder.KnownMotifs | **Internal** | Consensus strings verified equal to cited sources. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported element's `Sequence` (the matched substring) has the same length as its `Pattern`, and occurs at `Position` (0-based) in the input. | Yes | Scanning-window definition |
| INV-2 | Each reported element's `Sequence` IUPAC-matches its `Pattern`. | Yes | IUPAC matching (FindDegenerateMotif) |
| INV-3 | A consensus string and its element name are exactly those in the cited primary sources (no fabricated constants). | Yes | Sources §1.1 |
| INV-4 | Scanning is exhaustive: all 0-based start positions `0 <= i <= n - m` that match are reported (no missed/extra occurrences). | Yes | Window definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | TATA box detected | `GGGTATAAAGGG` | one "TATA Box" at pos 3, Pattern `TATAAA`, Sequence `TATAAA` | Bucher (1990) |
| M2 | -10 box detected | `CCTATAATCC` | one "-10 Box" at pos 2, Pattern `TATAAT` | Harley & Reynolds (1987) |
| M3 | -35 box detected | `AATTGACAGG` | one "-35 Box" at pos 2, Pattern `TTGACA` | Harley & Reynolds (1987) |
| M4 | CAAT box detected | `GGCCAATGG` | one "CAAT Box" at pos 2, Pattern `CCAAT` | Bucher (1990) |
| M5 | GC box detected | `AAGGGCGGTT` | one "GC Box" at pos 2, Pattern `GGGCGG` | Lundin et al. (1994) |
| M6 | Kozak detected | `TTGCCGCCACCATGGAA` | one "Kozak" at pos 2, Pattern `GCCGCCACCATGG` | Kozak (1987) |
| M7 | Shine-Dalgarno detected | `TTAGGAGGTTT` | one "Shine-Dalgarno" at pos 2, Pattern `AGGAGG` | Wikipedia (primaries) |
| M8 | Poly(A) signal detected | `CCAATAAACC` | one "Poly(A) Signal" at pos 2, Pattern `AATAAA` | Proudfoot & Brownlee (1976) |
| M9 | E-box degenerate match | `GGCACGTGGG` | one "E-box" at pos 2, Pattern `CANNTG`, Sequence `CACGTG` | Massari & Murre (2000) |
| M10 | AP-1 detected (corrected consensus) | `AATGACTCAGG` | one "AP-1" at pos 2, Pattern `TGACTCA` | Lee, Mitchell & Tjian (1987) |
| M11 | AP-1 old wrong pattern not matched | `AATGAGTCAGG` | zero "AP-1" hits (TGAGTCA ≠ consensus) | Lee, Mitchell & Tjian (1987) |
| M12 | NF-κB detected | `AAGGGACTTTCCAA` | one "NF-κB" at pos 2, Pattern `GGGACTTTCC` | Sen & Baltimore (1986) |
| M13 | CREB detected | `CCTGACGTCAGG` | one "CREB" at pos 2, Pattern `TGACGTCA` | Montminy et al. (1986) |
| M14 | Null sequence | `null` | `ArgumentNullException` | Contract |
| M15 | Empty sequence | `""` | empty result | Window definition |
| M16 | Constants equal source consensus | `KnownMotifs.*` | each constant equals its cited consensus string | §1.1 sources / INV-3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Multiple occurrences of one element | `AATAAACGAATAAA` | two "Poly(A) Signal" hits at pos 0 and 8 | INV-4 exhaustiveness |
| S2 | Multiple distinct elements | `TATAAA` + `AGGAGG` in one sequence | both "TATA Box" and "Shine-Dalgarno" reported | corner case |
| S3 | No regulatory element present | `GGGGCCCCGGGG` (no consensus) | empty result | corner case |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Matched substring length invariant | any detected element | `Sequence.Length == Pattern.Length` | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinderTests.cs` contained a legacy "Regulatory Element Tests" region (`FindRegulatoryElements_FindsTataBox/FindsPolyASignal/FindsEBox/ReturnsDescription`, `KnownMotifs_ContainsExpectedPatterns`) plus `FindRegulatoryElements_NullSequence_ThrowsException` in the Edge Cases region.
- No canonical `MotifFinder_FindRegulatoryElements_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 TATA box | ⚠ Weak | legacy used permissive `.Any(e => e.Name=="TATA Box")`, no position/pattern check |
| M2 -10 box | ❌ Missing | element did not exist before this unit |
| M3 -35 box | ❌ Missing | element did not exist before this unit |
| M4 CAAT box | ❌ Missing | not tested |
| M5 GC box | ❌ Missing | not tested |
| M6 Kozak | ❌ Missing | not tested |
| M7 Shine-Dalgarno | ❌ Missing | only asserted via constant |
| M8 Poly(A) signal | ⚠ Weak | legacy `.Any()` only |
| M9 E-box degenerate | ⚠ Weak | legacy `.Any()` only |
| M10 AP-1 corrected | ❌ Missing | not tested; pattern was wrong |
| M11 AP-1 wrong-pattern regression | ❌ Missing | not tested |
| M12 NF-κB | ❌ Missing | not tested |
| M13 CREB | ❌ Missing | not tested |
| M14 null | ✅ Covered | legacy null test (moved to canonical file) |
| M15 empty | ❌ Missing | not tested |
| M16 constants | ⚠ Weak | legacy checked only 3 constants |
| S1 multiple occurrences | ❌ Missing | not tested |
| S2 multiple distinct elements | ❌ Missing | not tested |
| S3 no element | ❌ Missing | not tested |
| C1 length invariant | ❌ Missing | not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindRegulatoryElements_Tests.cs` — all MUST/SHOULD/COULD cases.
- **Remove:** the "Regulatory Element Tests" region and the `FindRegulatoryElements_NullSequence_ThrowsException` test from `MotifFinderTests.cs` (duplicates/weak), replaced by a NOTE pointer.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_FindRegulatoryElements_Tests.cs` | Canonical MOTIF-REGULATORY-001 | 20 |
| `MotifFinderTests.cs` | Other MotifFinder methods (regulatory region removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with exact pos/pattern/sequence | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ⚠ Weak | rewrote with exact pos | ✅ Done |
| 9 | M9 | ⚠ Weak | rewrote with exact pos/sequence | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | M11 | ❌ Missing | implemented (regression) | ✅ Done |
| 12 | M12 | ❌ Missing | implemented | ✅ Done |
| 13 | M13 | ❌ Missing | implemented | ✅ Done |
| 14 | M14 | ✅ Covered | moved to canonical file | ✅ Done |
| 15 | M15 | ❌ Missing | implemented | ✅ Done |
| 16 | M16 | ⚠ Weak | rewrote — all 12 constants vs sources | ✅ Done |
| 17 | S1 | ❌ Missing | implemented | ✅ Done |
| 18 | S2 | ❌ Missing | implemented | ✅ Done |
| 19 | S3 | ❌ Missing | implemented | ✅ Done |
| 20 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M16 | ✅ | canonical file, exact evidence-based values |
| S1–S3 | ✅ | canonical file |
| C1 | ✅ | canonical file (invariant) |

All 20 in-scope cases ✅.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | NF-κB scanned as the strong reference site `GGGACTTTCC` rather than the full degenerate consensus `GGGRNWYYCC`. The exact string is source-backed (Sen & Baltimore 1986; p50-binding reference). | M12 |
| 2 | Kozak scanned as the single most-preferred-base string `GCCGCCACCATGG` rather than expanding -3 purine / +4 G degeneracy. | M6 |

---

## 7. Open Questions / Decisions

1. **Decision:** AP-1 pattern corrected `TGAGTCA` → `TGACTCA` per Lee, Mitchell & Tjian (1987). Existing repository value was a defect.
2. **Decision:** Added prokaryotic -10 (`TATAAT`) and -35 (`TTGACA`) hexamers (Harley & Reynolds 1987), which the task note lists as expected regulatory elements and which were absent before.
3. **Decision:** Used the per-position IUPAC scan (`FindDegenerateMotif`); the suffix tree was not used because the IUPAC E-box pattern is degenerate (not an exact substring), so a single linear scan over a short pattern set is the correct algorithm.
