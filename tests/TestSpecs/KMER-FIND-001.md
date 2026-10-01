# Test Specification: KMER-FIND-001

**Test Unit ID:** KMER-FIND-001
**Area:** K-mer Analysis
**Title:** K-mer Search (Most Frequent, Unique, Clumps)
**Created:** 2026-01-23
**Status:** Complete

---

## Canonical Methods

| Method | Class | Type | Description |
|--------|-------|------|-------------|
| `FindMostFrequentKmers(string, k)` | KmerAnalyzer | Canonical | Find k-mers with maximum count |
| `FindUniqueKmers(string, k)` | KmerAnalyzer | Canonical | Find k-mers appearing exactly once |
| `FindClumps(string, k, windowSize, minOccurrences)` | KmerAnalyzer | Canonical | Find (L,t)-clumps |

---

## Evidence Summary

| Source | Type | Key Contributions |
|--------|------|-------------------|
| Wikipedia (K-mer) | Primary | K-mer definition, unique k-mers for genomic fingerprinting |
| Rosalind BA1B | Primary | Most frequent k-mers problem, sample dataset |
| Rosalind BA1E | Primary | (L,t)-clump definition, sample dataset with expected outputs |

---

## Must Tests (M) - Evidence-Backed

### M1: FindMostFrequentKmers - Rosalind BA1B Sample Dataset
**Source:** Rosalind BA1B
**Test:** Sequence "ACGTTGCATGTCGCATGATGCATGAGAGCT" with k=4
**Expected:** Contains "CATG" and "GCAT" (both appear 3 times)

### M2: FindMostFrequentKmers - Multiple k-mers with Same Max Count
**Source:** Rosalind BA1B definition: "All most frequent k-mers"
**Test:** When multiple k-mers share maximum count, all should be returned
**Cases:**
- "ACGTACGT" k=4: only "ACGT" appears twice, others once → returns {"ACGT"}
- Sequence where 2+ k-mers tie for max → returns all

### M3: FindMostFrequentKmers - Single Most Frequent
**Source:** Rosalind BA1B
**Test:** "AAACGT" with k=2 → "AA" appears twice, others once → returns {"AA"}

### M4: FindMostFrequentKmers - Empty Sequence
**Source:** Implementation contract (edge case)
**Test:** Empty sequence → returns empty collection

### M5: FindMostFrequentKmers - k > Sequence Length
**Source:** Wikipedia pseudocode (L - k + 1 k-mers; if L < k, no k-mers)
**Test:** Short sequence with large k → returns empty collection

### M6: FindUniqueKmers - Basic Uniqueness
**Source:** Wikipedia (K-mer): k-mers appearing exactly once
**Test:** "ACGTACGT" k=4:
- ACGT appears 2× → NOT unique
- CGTA, GTAC, TACG appear 1× → unique
**Expected:** Returns {CGTA, GTAC, TACG}

### M7: FindUniqueKmers - All Unique
**Source:** Mathematical invariant
**Test:** When no k-mer appears more than once, all k-mers are unique
**Case:** "ACGT" k=2 → 3 distinct k-mers (AC, CG, GT), each appears once

### M8: FindUniqueKmers - No Unique (Homopolymer)
**Source:** Edge case
**Test:** "AAAA" k=2 → only "AA" which appears 3× → no unique k-mers

### M9: FindClumps - Rosalind BA1E Sample Dataset
**Source:** Rosalind BA1E
**Test:**
```
Sequence: CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
k=5, L=75, t=4
```
**Expected:** Contains "CGACA", "GAAGA", "AATGT"

### M10: FindClumps - Simple Clump Detection
**Source:** Rosalind BA1E definition
**Test:** "AAAAA" k=3, windowSize=5, minOccurrences=3
**Expected:** "AAA" forms a clump (appears 3× in 5-bp window)

### M11: FindClumps - No Clump Found
**Source:** Rosalind BA1E definition
**Test:** "ACGT" k=2, windowSize=4, minOccurrences=3
**Expected:** Empty (no k-mer appears 3× in window)

### M12: FindClumps - Invalid Parameters
**Source:** Implementation contract
**Tests:**
- k > windowSize → empty
- Empty sequence → empty
- windowSize > sequence length → empty

---

## Should Tests (S) - Additional Coverage

### S1: FindMostFrequentKmers - Case Insensitivity
**Test:** Mixed case sequences return same results as uppercase

### S2: FindMostFrequentKmers - Single Character Sequence
**Test:** "A" with k=1 → returns {"A"}

### S3: FindUniqueKmers - Empty Sequence
**Test:** Empty sequence → returns empty collection

### S4: FindClumps - Clump at Boundary
**Test:** Clump exists only at start or end of sequence

### S5: FindClumps - Multiple Distinct Clumps
**Test:** Sequence with multiple different k-mers forming clumps

---

## Could Tests (C) - Extended Coverage

### C1: Performance with Large Sequence
**Test:** FindClumps on sequence > 10,000 bp completes in reasonable time

### C2: Overlapping K-mers in Clumps
**Test:** Verify clumps correctly count overlapping occurrences

---

### C3: FindClumps - Window boundary (2026-09 review)
**Source:** Rosalind BA1E statement example (TGCA forms a (25,3)-clump in `gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac`)
**Test:** TGCA at 18, 23, 36 → three occurrences span 22 bp: L=25 → {TGCA}; L=22 → {TGCA}; L=21 → ∅ (occurrence must lie wholly inside the window). Lowercase input gives the same result.

### C4: FindClumps - Last window, overlaps, t=1
**Tests:** `CGTACGTTTTT` k=3 L=5 t=3 → {TTT} (only the final window); `AAAAA` k=3 L=5 t=3 → exactly {AAA}; `ACGTAC` k=2 L=4 t=1 → {AC, CG, GT, TA}, each once.

### Out-of-suite cross-check (not a unit test — 4.6 MB input)
E. coli genome (Compeau & Pevzner textbook dataset), k=9, L=500, t=3 → 1904 distinct 9-mers; Seqeron result set-equal to an independent Python reference.

---

## Test Audit

### Existing Tests (KmerAnalyzerTests.cs)

| Test Name | Verdict | Action |
|-----------|---------|--------|
| `FindMostFrequentKmers_SingleMostFrequent_ReturnsIt` | Weak | Replace with M3 |
| `FindMostFrequentKmers_MultipleMostFrequent_ReturnsAll` | Weak | Replace with M1, M2 |
| `FindMostFrequentKmers_EmptySequence_ReturnsEmpty` | Keep | Maps to M4 |
| `FindUniqueKmers_ReturnsKmersAppearingOnce` | Keep | Maps to M6 |
| `FindClumps_SimpleClump_Found` | Keep | Maps to M10 |
| `FindClumps_NoClump_ReturnsEmpty` | Keep | Maps to M11 |
| `FindClumps_InvalidParameters_ReturnsEmpty` | Keep | Maps to M12 |
| `FindClumps_RepetitiveRegion` | Weak | Remove (not evidence-based) |

### Consolidation Plan

1. Create `KmerAnalyzer_Find_Tests.cs` with all KMER-FIND-001 tests
2. Remove FindMostFrequentKmers, FindUniqueKmers, FindClumps tests from KmerAnalyzerTests.cs
3. Keep auxiliary methods (KmerDistance, GenerateAllKmers, FindKmerPositions, AnalyzeKmers, FindKmersWithMinCount) in KmerAnalyzerTests.cs

### Audit round 1, WP3 (B06, F12) — `FindClumpWindows`

`FindClumpWindows(sequence, k, L, t)` returns each clump k-mer with the maximal runs of window starts i for which
`Genome[i..i+L−1]` holds ≥ t occurrences (Compeau & Pevzner ch. 1 definition; Rosalind BA1E window convention).
`FindClumps` and `FindClumpWindows` share one sliding pass. Tests (`KmerAnalyzer_HistogramClumpWindowsFilters_Tests`):

| ID | Test | Evidence |
|---|---|---|
| W1 | `FindClumpWindows_MatchesBruteForce` ×7 | Python brute force over every window: BA1E `CGACA:0-6 GAAGA:0-16 AATGT:16-21`; BA1B L=12 t=2 split runs `GCAT:4-5,11-12 CATG:5-6,12-13 ATGA:13-14`; homopolymers; lower case |
| W2 | `FindClumpWindows_SameKmerSetAsFindClumps` ×7 | the k-mer set equals `FindClumps` |
| W3 | `FindClumpWindows_FirstWindowStart_IsLeftmostQualifyingWindow` | AATGT at 21, 73, 81, 86 → windows [86+5−75, 21] = [16, 21] |
| W4 | `FindClumpWindows_SameLeavingAndEnteringKmer_DoesNotSplitRun` | A^10 k=2 L=4 t=3 → one run 0–6 |
| W5 | `FindClumpWindows_DegenerateInput_Empty` ×6 | same empty conditions as `FindClumps` |

Out-of-suite: 3000 random cases vs the brute force, 0 mismatches; 2 Mbp random (k=6, L=300, t=4) vs an independent
occurrence-interval method, 113 k-mers / 142 runs identical (K-mer_Search.md §7.3).

---

## Assumptions

None. All tests backed by Rosalind or Wikipedia sources.

---

*TestSpec generated: 2026-01-23*
