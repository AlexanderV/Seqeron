# Test Specification: GENOMIC-COMMON-001

**Test Unit ID:** GENOMIC-COMMON-001
**Area:** Analysis
**Algorithm:** Longest Common Substring / Common Region Detection (generalized suffix tree)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-18 (B09 review)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Gusfield (1997), Algorithms on Strings, Trees and Sequences (ISBN 0-521-58519-8) | 1 | citation retrieved from https://en.wikipedia.org/wiki/Longest_common_substring | 2026-06-13 |
| 2 | Wikipedia — Longest common substring | 4 | https://en.wikipedia.org/wiki/Longest_common_substring | 2026-06-13 |
| 3 | GeeksforGeeks — Suffix Tree Application 5 (LCS) | 3 | https://www.geeksforgeeks.org/dsa/suffix-tree-application-5-longest-common-substring-2/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. LCS = "a longest string which is substring of both S and T" — Wikipedia (source 2).
2. The common region is **contiguous**, unlike a subsequence — Wikipedia (source 2).
3. Found via a generalized suffix tree: deepest node whose subtree has leaves from both strings; Θ(n+m) — Gusfield 1997 / Wikipedia (sources 1, 2); GeeksforGeeks O(M+N) (source 3).
4. Ties exist: "BADANAT"/"CANADAS" share two maximal substrings "ADA" and "ANA" — Wikipedia (source 2).
5. Worked length example: "xabxac"/"abcabxabcd" → "abxa", length 4 — GeeksforGeeks (source 3).

### 1.3 Documented Corner Cases

- Multiple maximal substrings (ties): a deterministic representative must be chosen — Wikipedia (source 2).
- No common character → only the empty string qualifies → length-0 LCS — derived from the Wikipedia definition (source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing longest common **substring** (contiguous) with longest common **subsequence** (gapped) — Wikipedia (source 2). The repository XML doc previously mislabeled the method "subsequence"; the implementation computes substring.
2. Nondeterministic tie-break — must return a documented, deterministic representative — Wikipedia (source 2) + SuffixTree XML doc.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindLongestCommonRegion(seq1, seq2)` | GenomicAnalyzer | Canonical | longest contiguous common substring with positions |
| `FindCommonRegions(seq1, seq2, minLength)` | GenomicAnalyzer | Canonical | for each start position in seq2, the single longest common substring of length ≥ minLength (right-maximal, deduplicated) — NOT every common substring |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | The returned substring is contiguous and occurs in **both** sequences at the reported 0-based positions | Yes | Wikipedia definition (source 2) |
| INV-2 | No common contiguous substring strictly longer than the returned one exists | Yes | Wikipedia definition (source 2) |
| INV-3 | Empty / no-shared-character input → `CommonRegion.None` (empty, length 0, positions −1) | Yes | Wikipedia definition (source 2) |
| INV-4 | Every region from `FindCommonRegions(minLength)` has length ≥ minLength, is contiguous, occurs in both sequences at the reported 0-based positions, and is the single longest match at its start position in sequence2 (right-maximal; shorter prefixes at the same start are NOT reported) | Yes | substring/contiguity definition (source 2); right-maximal contract derived from the implementation, hand-verified |
| INV-5 | On ties, the representative is "first found in sequence2" (deterministic) | Yes | **ASSUMPTION** (SuffixTree XML doc; see §6) |
| INV-6 | Reported positions are the **first (lowest) 0-based occurrence** in each sequence (not an arbitrary suffix-tree-order occurrence) | Yes | independent Python `str.find` oracle (B09 F15) |
| INV-7 | Each `FindCommonRegions` region length equals the matching statistic ms(i) of its start i (longest prefix of seq2[i..] occurring in seq1), for every i with ms(i) ≥ max(1, minLength) | Yes | Chang & Lawler 1994; brute-force `in` test (B09 F16) |
| INV-8 | `null` sequence → `ArgumentNullException` (eager, also for the lazy `FindCommonRegions`); lower case folded by `DnaSequence`; non-ACGT rejected by `DnaSequence` | Yes | contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Unique LCS | `ACGTACGT` vs `TTACGTGG` | `TACGT`, length 5, pos1=3, pos2=1 | Wikipedia def (source 2); brute-force-verified |
| M2 | Tie | `CACAGAG` vs `TACATAGAT` | `ACA` (first-in-other of {ACA,AGA}), length 3, pos1=1, pos2=1 | Wikipedia tie property (source 2); SuffixTree tie-break |
| M3 | Identical | `ACGT` vs `ACGT` | `ACGT`, length 4, pos1=0, pos2=0 | Wikipedia def (source 2) |
| M4 | No common substring | `AAAA` vs `GGGG` | `CommonRegion.None` (empty, length 0, IsEmpty, pos −1/−1) | Wikipedia def (source 2) |
| M5 | Right-maximal regions, minLength=4 | `FindCommonRegions(ACGTACGT, TTACGTGG, 4)` | `{(TACGT,3,1),(ACGT,0,2)}` (longest match per start in seq2) | contiguity def (source 2); right-maximal contract hand-verified |
| M6 | Right-maximal regions, minLength=3 | `FindCommonRegions(ACGTACGT, TTACGTGG, 3)` | `{(TACGT,3,1),(ACGT,0,2),(CGT,1,3)}` — note the prefixes TAC/TACG/ACG, which are ALSO common substrings ≥3, are NOT returned (only the longest per start) | contiguity def (source 2); right-maximal contract hand-verified |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty sequence1 | `FindLongestCommonRegion("", ACGT)` | `CommonRegion.None` | empty handling (INV-3) |
| S2 | Empty sequence2 | `FindLongestCommonRegion(ACGT, "")` | `CommonRegion.None` | empty handling (INV-3) |
| S3 | `FindCommonRegions` no match | `FindCommonRegions(AAAA, GGGG, 2)` | empty enumeration | INV-4 |
| S4 | Single-char overlap | `FindLongestCommonRegion("A", "TTTAT")` | `A`, length 1, pos1=0, pos2=3 | minimal non-empty LCS |

### 4.2b B09 review tests (GENOMIC-COMMON-001 F15/F16)

| ID | Test | Expected |
|----|------|----------|
| B1 | `FindCommonRegions(AGGAAGAAGGAGGCCGATTACT, CGTGGAAGGAGGCCGA, 2)` | 13 regions; last = (`GA`, 2, 14) — first occurrence, tree order gave 15 |
| B2 | `FindLongestCommonRegion(CACAGATAAAGGACTGAGGCTCTACG, GA)` | (`GA`, 4, 0) |
| B3 | `ACACACACAC` / `GGACACGG` | (`ACAC`, 0, 2) |
| B4 | null arguments, lower case, non-ACGT (N, U), `minLength` = int.MinValue/0/1/6/int.MaxValue | contract (INV-8) |
| B5 | embedded oracle `TestData/Analysis/common_region_oracle.tsv` (283 pairs, independent Python DP + suffix automaton + brute-force ms) | LCS (tie, positions) and all regions equal |
| B6 | 150 random pairs vs in-test brute force (ms by `Contains`, distinct strings, `IndexOf` positions) | equal |
| B7 | 1 Mbp random pair with a planted 5 kb shared block | linear time (< 60 s guard, ≈ 4 s), planted block found, positions are first occurrences |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-1 property | returned substring occurs at reported positions in both | substring matches both slices | property test |
| C2 | INV-4 property | all FindCommonRegions results length ≥ minLength and ⊆ seq1 | all satisfy | O(n·m) property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `FindLongestCommonRegion` / `FindCommonRegions`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit, no prior test |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindCommonRegion_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomicAnalyzer_FindCommonRegion_Tests.cs` | canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | S1 | ❌ Missing | implemented | ✅ Done |
| 8 | S2 | ❌ Missing | implemented | ✅ Done |
| 9 | S3 | ❌ Missing | implemented | ✅ Done |
| 10 | S4 | ❌ Missing | implemented | ✅ Done |
| 11 | C1 | ❌ Missing | implemented | ✅ Done |
| 12 | C2 | ❌ Missing | implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindLongestCommonRegion_UniqueLcs_ReturnsTacgt` |
| M2 | ✅ Covered | `FindLongestCommonRegion_Tie_ReturnsFirstInOther` |
| M3 | ✅ Covered | `FindLongestCommonRegion_Identical_ReturnsWholeSequence` |
| M4 | ✅ Covered | `FindLongestCommonRegion_NoCommon_ReturnsNone` |
| M5 | ✅ Covered | `FindCommonRegions_MinLengthFour_ReturnsTwoRegions` |
| M6 | ✅ Covered | `FindCommonRegions_MinLengthThree_ReturnsThreeRegions` |
| S1 | ✅ Covered | `FindLongestCommonRegion_EmptyFirst_ReturnsNone` |
| S2 | ✅ Covered | `FindLongestCommonRegion_EmptySecond_ReturnsNone` |
| S3 | ✅ Covered | `FindCommonRegions_NoMatch_ReturnsEmpty` |
| S4 | ✅ Covered | `FindLongestCommonRegion_SingleCharOverlap_ReturnsA` |
| C1 | ✅ Covered | `FindLongestCommonRegion_ReturnedSubstring_OccursInBothAtPositions` |
| C2 | ✅ Covered | `FindCommonRegions_AllRegions_SatisfyMinLengthAndContainment` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Tie-break = "first found in sequence2" (deterministic, documented in SuffixTree XML doc; does not change maximal length, only the representative) | INV-5, M2 |

---

## 7. Open Questions / Decisions

1. The repository method `FindLongestCommonRegion` computes the longest common **substring** (contiguous) but its XML doc previously said "subsequence". Corrected in Phase 5; behavior unchanged (the SuffixTree call already computed substring).

---

## 2026-11 audit round 2 (A2-40, B09 F86): cancellation

`FindCommonRegions(s1, s2, minLength, CancellationToken)` polls the matching-statistics scan (every 1024 positions) and the region enumeration; tests in `SequenceAligner_Cancellation_Tests` (`PreCancelled_GenomicAnalyzerScans_Throw`, `CancelledWhileRunning_GenomicAnalyzerScans_StopPromptly`).

## Round 2 additions (audit A2-30 / A2-35, F94)

| Id | Test | Requirement |
|----|------|-------------|
| M-A2-30a | `FindCommonRegions_ByMode_BothStrands_MatchesMummerAndBruteForce` | 4716 rows: `CommonRegionMode` x `StrandSelection.Both` equals mummer 3.23 `-maxmatch / -mumreference / -mum -b` (== brute force; N rows brute force) |
| M-A2-30b | `FindCommonRegions_StrandSelection_AndDnaSequenceOverload`, `RightMaximalMode_ForwardEqualsTheExistingMethod_AndReverseIsTheSameEngineOnTheReverseComplement` | strand selection, `DnaSequence` overload, reverse positions are forward-strand windows of revcomp(Sequence) |
| M-A2-30c | `FindLongestCommonRegion_Strands_AgreeWithTheMaxmatchOracle` | longest region per strand / both (tie -> forward) |
| M-A2-30d | `LongestCommonRegionWithMismatches_MatchesBruteForce_AllRows`, `..._ZeroMismatches_IsTheExactMethod`, `..._ReverseStrand` | k-mismatch LCS (1680 rows), k = 0 equals the exact method, reverse strand |
| M-A2-30e | `FindLongestCommonRegion_SeveralSequences_MatchesBruteForce_AllRows` | N-sequence LCS with support q (694 rows) |
| S-A2-35 | `StringOverloads_WithN_ReportAnOrdinaryRegionNextToAnNRunUnchanged`, `..._Validation` | string / IUPAC input: N never matches |

| A3-13a | `SeqMatchAll_EqualsSeqmatchall_AllSets`, `SeqMatchAll_OnNucleotides_IsTheMaximalExactMatchModeOfFindCommonRegions`, `SeqMatchAll_WorkedExample_AndArguments` | 330 seqmatchall 6.6.0 sets (`seqmatchall_oracle.tsv`: nucleotide / N / lower case / protein, k 2-12) | matches equal as sets; the N-free nucleotide case equals `FindCommonRegions` MaximalExactMatches | A3-13 |
| A3-13b | `MegaMerge_EqualsMegamerger_AllCases`, `WorkedValues_RunOnTheBinary`, `Validation` | 600 megamerger 6.6.0 runs (`megamerger_oracle.jsonl`): merged sequence and report text, minimal match set, every start / end / mismatch case, `-prefer`, the empty merge | as specified | A3-13b |

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

`EmbossGapCharacters_Oracle_Tests` (`emboss_gap_oracle.tsv`, `generate_emboss_gap_oracle.py`): megamerger (240) and seqmatchall (240) - the installed EMBOSS 6.6.0 `megamerger, seqmatchall` on inputs that contain `- . ~ ? # + _`, digits, blanks, punctuation, `*` and non-ASCII letters equal `GenomicAnalyzer.MegaMerge`, `GenomicAnalyzer.SeqMatchAll` on every row, and the result of the cleaned input equals the result of the input. The reader rule (one shared `EmbossSequenceReader.Clean` / `CleanGapped`, Infrastructure): only `A-Z a-z *` survive in a plain sequence (`ajseqread.c seqAppend` + `ajseqtype.c ajSeqTypeCheckIn` + `ajstr.c ajStrRemoveGap`); coordinates and lengths count the cleaned sequence. Also `Clean_KeepsLettersAndStar_AndDropsEveryOtherCharacter` (all 256 code points) and `IsNucleotide_IsTheAjSeqIsNucRule_OfEveryPort`.

## Audit A4-1 (B09 F157): cancellation

`SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `MegaMerge` (`MinimalMatchSet`) poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

## Audit A4-2 (B09 F158): resource guards

`SequenceAligner_ResourceLimits_Tests`: `SeqMatchAll` / `MegaMerge` (`WordMatchMaxMatches`) - documented limit constants (`SequenceAligner.ResourceLimits.cs`), the exception of the `SequenceAligner.ResourceLimit` style raised before the allocation in 64-bit arithmetic, and a low-allocation assertion per limit.

## Audit A4-9 (B09 F164): megamerger edge cases

| ID | Test | Locks |
|---|---|---|
| A4-9 | `GenomicAnalyzer_MegaMerge_Tests.MegaMerge_EqualsMegamerger_OnTheEdgeCaseRows` | 480 rows (`megamerger_edge_oracle.jsonl`): lower-case-only, one lower / one upper, N runs longer than the word, short / identical / nested sequences, status ok / nomerge / abort |


## Audit R14-4a (B09 F291): matching statistics by the suffix automaton

`GenomicAnalyzer_CommonRegionsAutomaton_Tests`: the right-maximal regions (DnaSequence, string and coordinates-only overloads, forward / reverse / both strands) equal an independent brute force (ms(i) by naive extension, distinct strings, lowest occurrence by `IndexOf`) on 400 random / low-complexity pairs and 120 reverse-strand pairs; the documented example `ACGTACGT` / `TTACGTGG` min 3 = {TACGT, ACGT, CGT}; poly-A and (AC)n of 20 000 run in linear time, allocate < 60 MB for the coordinates and build no suffix tree; strings equal spans; cancellation; null. MCP `find_common_regions` `includeSequence = false` = the same list without the strings.


## Audit R14-22 (B09 F295): argument contract

`GenomicAnalyzer_ArgumentContract_Tests` (minLength <= 0 = minLength 1 for every finder; tandem structural parameters throw) and MCP `find_common_regions` refuses `minLength < 1`; see `docs/algorithms/Sequence_Comparison/Genomic_Argument_Contract.md`.

## Audit R14-17 (B09 F298): MCP gaps

MCP `FindLongestCommonRegionMultiTests` (Rosalind LCSM GATTACA / TAGACCA / ATACA = AC, CA, TA with first occurrences; minSupport 2; N runs; no common string; validation), SuffixTree.Mcp.Core `FindLongestRepeatCommonRegionOptionsTests` (N / lower case, strands ReverseComplement / Both, k-mismatch LCS, both-strand longest repeat with reversePositions, cancellation, IUPAC similarity): expected values from an independent brute force.
