# Test Specification: COMPGEN-DOTPLOT-001

**Test Unit ID:** COMPGEN-DOTPLOT-001
**Area:** Comparative
**Algorithm:** Dot Plot Generation (word-match / k-tuple dot matrix)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Gibbs & McIntyre (1970), *Eur. J. Biochem.* 16:1–11 | 1 | https://doi.org/10.1111/j.1432-1033.1970.tb01046.x | 2026-06-14 |
| 2 | EMBOSS dottup manual (word-match dot plot) | 3 | https://www.bioinformatics.nl/cgi-bin/emboss/help/dottup | 2026-06-14 |
| 3 | EMBOSS dottup manpage (default wordsize 10) | 3 | https://manpages.ubuntu.com/manpages/xenial/man1/dottup.1e.html | 2026-06-14 |
| 4 | Wikipedia — Dot plot (bioinformatics) | 4 | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | 2026-06-14 |
| 5 | Huttley TIB — Dotplot (worked example) | 4 | https://gavinhuttley.github.io/tib/seqcomp/dotplot.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. A dot is placed at (i, j) iff the word starting at i in sequence1 exactly matches the word starting at j in sequence2 — "dottup looks for places where words (tuples) of a specified length have an exact match in both sequences" (Source 2); for k=1 "if X[i] == Y[j]" (Source 5).
2. Default word size = 10 (Source 3).
3. Self-comparison produces the full main diagonal — "The main diagonal represents the sequence's alignment with itself" (Source 4).
4. Longer words reduce noise/sensitivity, shorter words increase both (Source 2).
5. Worked example: X=`AGCGT`, Y=`AT`, k=1 ⇒ dots at (0,0) and (4,1) (Source 5).

### 1.3 Documented Corner Cases

- Sequence shorter than word size ⇒ no words ⇒ no dots (Source 2/3).
- Self-comparison ⇒ full main diagonal (Source 4).
- Short words ⇒ extra off-diagonal chance dots (documented noise, Source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing `dottup` (exact word match — this unit) with `dotmatcher` (substitution-matrix scored window) — Source 2.
2. Forgetting overlapping occurrences: all matching start positions in sequence2 are reported, not just the first (SuffixTree.FindAllOccurrences; Source 2 "all places").

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GenerateDotPlot(string, string, int wordSize=10, int stepSize=1)` | ComparativeGenomics | **Canonical** | Word-match dot plot; returns (x,y) match coordinates |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A pair (x,y) is returned iff `sequence1[x..x+w]` equals `sequence2[y..y+w]` (case-insensitive), w=wordSize | Yes | Sources 2, 5 |
| INV-2 | Self-comparison contains every (i,i) for 0≤i≤n−w (full main diagonal) | Yes | Source 4 |
| INV-3 | No dots when either sequence length < wordSize, is null, or is empty | Yes | Sources 2, 3 |
| INV-4 | Number of returned pairs ≤ (#word-starts sampled) × (sequence2 length); each x is a multiple of stepSize | Yes | O(n×m) word-match model (Source 2) |
| INV-5 | wordSize ≤ 0 or stepSize ≤ 0 throws ArgumentOutOfRangeException | Yes | Sibling validation convention; undefined window for dottup |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Huttley k=1 example | `AGCGT` vs `AT`, wordSize=1 | exactly {(0,0),(4,1)} | Source 5 worked example |
| M2 | Exact word, wordSize=4 | `ACGTACGT` self, wordSize=4 (every word start x=0..4) | exactly {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} | Source 2 exact word match + all overlapping occurrences |
| M3 | Self main diagonal | `ACGT` self, wordSize=1 | contains (0,0),(1,1),(2,2),(3,3) | Source 4 main diagonal |
| M4 | No similarity | `AAAA` vs `CCCC`, wordSize=1 | empty | dot-only-on-match rule (Sources 2,5) |
| M5 | Word longer than sequence | `ACG` vs `ACGT`, wordSize=4 | empty | Source 2/3 word formation |
| M6 | Null / empty inputs | null/`""` either arg | empty | INV-3 |
| M7 | Invalid wordSize | wordSize=0 and −1 | throws ArgumentOutOfRangeException | INV-5 |
| M8 | Invalid stepSize | stepSize=0 and −1 | throws ArgumentOutOfRangeException | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | stepSize sampling | `ACGTACGT` self, wordSize=4, stepSize=4 | only x∈{0,4}: {(0,0),(0,4),(4,0),(4,4)} | every x is a multiple of stepSize (INV-4) |
| S2 | Case-insensitivity | `acgtacgt` vs `ACGTACGT`, wordSize=4 | same as M2 set {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} | implementation upper-cases both |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-4 property | random-free constructed pair | all x multiples of stepSize; count ≤ bound | property check of O(n×m) bound |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior test file for `GenerateDotPlot`. Sibling fixtures exist for other `ComparativeGenomics` methods under `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_*_Tests.cs`. New canonical file created: `ComparativeGenomics_GenerateDotPlot_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| M8 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_GenerateDotPlot_Tests.cs` — all cases for this unit.
- **Remove:** nothing (no pre-existing tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| ComparativeGenomics_GenerateDotPlot_Tests.cs | Canonical | 11 |

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
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | GenerateDotPlot_HuttleyExample_ReturnsTwoDots |
| M2 | ✅ | GenerateDotPlot_ExactWordMatch_ReturnsAllOverlappingOccurrences |
| M3 | ✅ | GenerateDotPlot_SelfComparison_ContainsMainDiagonal |
| M4 | ✅ | GenerateDotPlot_NoSharedResidues_ReturnsEmpty |
| M5 | ✅ | GenerateDotPlot_WordLongerThanSequence_ReturnsEmpty |
| M6 | ✅ | GenerateDotPlot_NullOrEmptyInput_ReturnsEmpty |
| M7 | ✅ | GenerateDotPlot_NonPositiveWordSize_Throws |
| M8 | ✅ | GenerateDotPlot_NonPositiveStepSize_Throws |
| S1 | ✅ | GenerateDotPlot_StepSizeFour_SamplesEveryFourthPosition |
| S2 | ✅ | GenerateDotPlot_CaseInsensitive_MatchesMixedCase |
| C1 | ✅ | GenerateDotPlot_AllXCoordinatesAreMultiplesOfStep_Property |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | x = sequence1, y = sequence2 (axis orientation is presentation convention) | M1, M2, S1 |
| 2 | Case-insensitive comparison (both sequences upper-cased) | S2 |

---

## 7. Open Questions / Decisions

1. Decision: build on the repository SuffixTree (exact-match enumeration, many word queries against one text — the suffix tree's strong case) rather than a naive O(n×m) scan. Output is identical; recorded in the algorithm doc §5.2.
