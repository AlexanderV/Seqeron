# Validation Report: GENOMIC-COMMON-001 — Common Region Detection (Longest Common Substring)

- **Validated:** 2026-06-15   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2)`, `GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength)` (underlying: `SuffixTree.LongestCommonSubstringInfo` → `SuffixTreeAlgorithms.FindAllLcs`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** CLEAN (defect found and completely fixed this session)

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS-WITH-NOTES. **State:** FIXED.
- **Fixes (B09.md, F-numbers):** F15 `PositionInFirst` = first (lowest) occurrence, F16 `FindCommonRegions` quadratic -> matching statistics (Chang & Lawler; 1 Mbp pair ~ 4 s), F17 MCP / doc claims.
- **Audit round 2 (A2-30 / A2-35, F94 / F99):** `CommonRegionMode` (right-maximal, MEM, MUM-reference, MUM), strands (`mummer -r` / `-b`), `FindLongestCommonRegionWithMismatches` (k-mismatch LCS), N-sequence LCS, string / IUPAC overloads; oracles: `mummer` 3.23 on 4 716 (reference, query, L, mode) rows, 1 680 k-mismatch and 694 N-sequence brute-force rows (`scripts/analysis/generate_common_region_strand_oracle.py`); MCP `find_common_regions` (mode / strands), `find_longest_common_region_mismatches`.
- **Oracle numbers:** 283-pair independent oracle (O(nm) DP + suffix automaton + brute-force matching statistics): every LCS and region equal after F15.
- **Open items:** registered, not doable here: a public per-position matching-statistics API on `ISuffixTreeAnalysis` (B05, QUEUE R17) would replace `GenomicAnalyzer.MatchingStatistics`; the present implementation is exact and non-quadratic.
- The unit's TestSpec (`tests/TestSpecs/GENOMIC-COMMON-001.md`), Evidence (`docs/Evidence/GENOMIC-COMMON-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Stage A — Description

### Sources opened & what they confirm (retrieved 2026-06-15)
- **Wikipedia, "Longest common substring"** (https://en.wikipedia.org/wiki/Longest_common_substring), fetched this session. Confirms verbatim:
  - Definition: *"Given two strings, S of length m and T of length n, find a longest string which is substring of both S and T."*
  - Contiguity: *"…the longest common substring problem seeks a contiguous substring shared by both texts."* (vs subsequence, which allows gaps).
  - Complexity: *"…in Θ(n + m) time with the help of a generalized suffix tree"*, attributed to Gusfield (1997).
  - Tie example: "BADANAT"/"CANADAS" share two maximal substrings **"ADA"** and **"ANA"**.
  - Three-string example: "ABABC","BABCA","ABCBA" → single LCS **"ABC"**, length 3.
- **Rosalind LCSM** (https://rosalind.info/problems/lcsm/), fetched this session. Confirms the shared-motif / longest-common-substring across a collection; sample {GATTACA, TAGACCA, ATACA} → **AC**; explicitly notes the LCS is **not necessarily unique** ("AA"/"CC" both LCS of "AACC"/"CCAA"), so any maximal answer is acceptable.

### Formula / definition check
The repository's `FindLongestCommonRegion` description (contiguous LCS via generalized suffix tree, Θ(n+m), deterministic first-in-`sequence2` tie-break) matches the sources exactly. Since sources state any maximal answer is acceptable, the documented deterministic tie-break is a valid representative.

### Edge-case semantics
Empty input / no shared character → length-0 LCS → `CommonRegion.None`: consistent with "a longest string which is substring of both" (only the empty string qualifies). Identical sequences → whole sequence at 0/0: a string is a substring of itself. All sourced and correct.

### Independent cross-check (my own brute force this session, independent of repo code)
All-substrings enumeration in Python (`s[i:j] in t`):

| Case | Inputs | Brute-force maximal | pos1 | pos2 (first) |
|------|--------|---------------------|------|------|
| M1 | ACGTACGT / TTACGTGG | TACGT (len 5) | 3 | 1 |
| M2 | CACAGAG / TACATAGAT | tie ACA, AGA (len 3); ACA first in t | ACA→1 | ACA→1 |
| M3 | ACGT / ACGT | ACGT (len 4) | 0 | 0 |
| M4 | AAAA / GGGG | none (len 0) | -1 | -1 |
| S4 | A / TTTAT | A (len 1) | 0 | 3 |
| C1 | GATTACACGT / CCTTACAGG | TTACA (len 5) | 2 | 2 |
- Rosalind LCSM cross-check: my brute force gives maximal {AC, CA, TA} length 2; Rosalind's published "AC" is among them → methodology validated.

### Findings / divergences (Stage A defect — fixed)
**DEFECT A1 (description, fixed):** The description of `FindCommonRegions` was internally inconsistent and partly **false**. Multiple locations claimed it returns *"all distinct common substrings of length ≥ minLength"* / *"All common substrings … enumerated"* and that the per-position binary search yields *"identical set of regions"* as a full deepest-common-node DFS enumeration:
- `GenomicAnalyzer.cs` XML doc ("every distinct common **substring**")
- `Common_Region_Detection.md` §5.1, §5.3 (two lines)
- TestSpec §2 method-table note; INV-4 ("definition extended to all common substrings")

This is mathematically untrue. For `ACGTACGT` vs `TTACGTGG`, minLength 3, the **full** set of common substrings ≥ 3 is `{TAC, TACG, TACGT, ACG, ACGT, CGT}` (6), but the method (and the *correct* parts of the description — INV-4 prose, §4.1) returns only the **right-maximal** per-start set `{TACGT, ACGT, CGT}` (3); the prefixes `TAC, TACG, ACG` are dropped. Verified by my brute force this session.

**Resolution:** corrected the description in `GenomicAnalyzer.cs`, `Common_Region_Detection.md` (§5.1, §5.3, §6.2-style note), and the TestSpec (method table, INV-4, M5/M6 rows) to state the true contract: *for each start position in `sequence2`, the single longest common substring of length ≥ max(1, minLength), deduplicated — right-maximal matches, NOT every common substring.* The `FindLongestCommonRegion` description needed no change.

## Stage B — Implementation

### Code path reviewed
- `GenomicAnalyzer.cs:246-257` `FindLongestCommonRegion` → `SuffixTree.LongestCommonSubstringInfo` (`SuffixTree.Algorithms.cs:69`) → `SuffixTreeAlgorithms.FindAllLcs` (`SuffixTreeAlgorithms.cs:28-126`).
- `GenomicAnalyzer.cs:271-314` `FindCommonRegions`: per-start binary search over `seq2` using `tree.Contains`, dedup via `HashSet`, position from `FindAllOccurrences(bestMatch)[0]`.

### Formula realised correctly?
Yes. `FindLongestCommonRegion` computes the contiguous LCS via suffix-link streaming (Gusfield matching-statistics style); `maxLen` is updated only on strict increase, so the first match reaching the maximal length wins → the "first-found-in-`sequence2`" tie-break, exactly as documented. Recomputed M1–M4, S4, C1 against the code (tests) — every value matches my external brute force above.

### Cross-verification table recomputed vs code
All MUST/SHOULD values reproduced by the passing tests and matched to the independent brute-force numbers above. `FindCommonRegions` right-maximal sets (M5/M6/M7/S5) reproduced by an independent Python mirror of the per-start-longest algorithm and the full-substring brute force; the returned set is exactly the right-maximal subset.

### Variant/delegate consistency
`LongestCommonSubstring(string)` / `(ReadOnlySpan<char>)` delegate to `LongestCommonSubstringInfo`; `FindAllLongestCommonSubstrings` shares `FindAllLcs(firstOnly:false)`. Consistent.

### Numerical robustness
Pure integer index arithmetic; empty-input guard at `FindAllLcs:34`; binary-search bounds in `FindCommonRegions` are correct (`hi = seq2.Length - i`). No overflow/precision concerns for DNA-scale inputs.

### Test quality audit (HARD gate)
**DEFECT B1 (tests, fixed):** M5/M6 carried comments claiming *"Brute-force-verified set"* for the **all-common-substrings** definition, while the asserted 3-/2-element sets are actually the **right-maximal** sets — a genuine brute force of "all common substrings ≥ minLength" yields 6/3 elements, not 3/2. The expected values were therefore code-echoes mis-attributed to a definition the code does not implement.
**Resolution:**
- Rewrote M5/M6 comments to truthfully describe the right-maximal contract and to call out the omitted prefixes.
- Added **M7** (`FindCommonRegions_DoesNotReturnPrefixesOfLongerMatches`): asserts TAC/TACG/ACG (genuine common substrings ≥ 3, hand-enumerated) are NOT returned, and the result equals the right-maximal `{TACGT, ACGT, CGT}` — locks the actual contract against an external hand computation.
- Added **S5** (`FindCommonRegions_MinLengthBelowOne_TreatedAsOne`): exercises the documented `minLength < 1 → 1` branch (previously untested), hand-verified to `{ACGT@0/0, CGT@1/1, GT@2/2, T@3/3}` with no empty regions.

Other tests (M1–M4, S1–S4, C1, C2) assert exact sourced values and cover the Stage-A branches (unique LCS, tie, identical, no-match, empty seq1, empty seq2, single-char, position invariants). No weakened assertions, no skips, no widened tolerances. C2 uses property-style `GreaterThanOrEqualTo`/`Contains` appropriately for a variable-result invariant test (INV-4), with exact substring/position checks per region.

**Gate result: PASS.** Full unfiltered suite `Failed: 0, Passed: 6573` (1 `[Explicit]` benchmark skipped, by design). Changed test file builds warning-free; the 4 build warnings are pre-existing in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects
- Code behavior is correct and self-consistent; **no code-behavior change** was warranted. Both defects were a mislabeled contract (docs + test comments) overclaiming "all common substrings"; the implemented right-maximal-per-start semantics is the sensible, useful one and is now described and locked truthfully.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (LCS biology/maths correct & sourced; one false "all common substrings" claim for `FindCommonRegions` corrected).
- **Stage B:** PASS-WITH-NOTES (code correct; mis-documented/code-echoing `FindCommonRegions` tests corrected, contract-lock + edge-case tests added).
- **End-state:** CLEAN — defects A1/B1 completely fixed this session; build 0 errors, full suite Failed: 0.
- Logged in FINDINGS_REGISTER as a description+test-quality defect (no behavioral bug).

---

## B09 re-review (2026-11-18, campaign 2026-09)

- **Stage A:** PASS-WITH-NOTES — definitions unchanged (Wikipedia/Gusfield/GeeksforGeeks as above; Rosalind/Wikipedia re-fetch blocked in this campaign, earlier fetches stand). Added the matching-statistics basis (Chang & Lawler 1994, search snippets) for `FindCommonRegions`: the per-start longest right-maximal match is ms(i). Tie-break and "first occurrence" positions are now explicit.
- **Stage B:** PASS-WITH-NOTES → FIXED. Independent Python oracle (DP + suffix automaton + brute-force ms, 283 pairs, `scripts/analysis/generate_common_region_oracle.py`) and a 150-pair in-test brute force:
  - **F15** `PositionInFirst` was an arbitrary suffix-tree-order occurrence (`FindAllOccurrences(...)[0]`, `LongestCommonSubstringInfo`), e.g. `FindCommonRegions(AGGAAGAAGGAGGCCGATTACT, CGTGGAAGGAGGCCGA, 2)` last region `GA` at 15 (first occurrence 2); `FindLongestCommonRegion(CACAGATAAAGGACTGAGGCTCTACG, GA)` 15 (first 4). Now the lowest occurrence (`FindAllLongestCommonSubstrings` ascending lists / minimum of the occurrence list).
  - **F16** `FindCommonRegions` was O(m²) time and memory (binary search with `Substring(i, mid)` probes up to m characters; 100 kbp pair 2.6 s, 1 Mbp ≈ 4 min and GBs of garbage). Now matching statistics: short matches by `Contains` extension from ms(i−1)−1, long matches from the canonical `FindMaximalExactMatches` stream; 1 Mbp pair ≈ 4 s.
  - `null` arguments → `ArgumentNullException` (was `NullReferenceException`); MCP `find_common_regions` description said "All common regions" (it is the right-maximal per-start set) — corrected; `find_longest_common_region` doc example (`ATGCATGC`/`CATGCAT` → `ATGC`) was wrong (LCS is `ATGCAT`, 0/1, length 6).
- **End-state:** CLEAN. Build 0 errors; unit/combinatorial/property/metamorphic/fuzz tests of the unit green (247).
- **Audit r3 (A3-13, F139):** `GenomicAnalyzer.SeqMatchAll` (EMBOSS seqmatchall, MCP `seqmatchall`), 330 / 330 oracle sets equal.
- **Audit r3 (A3-13b, F149):** EMBOSS `megamerger` ported (`GenomicAnalyzer.MegaMerge`, MCP `megamerger`): 600 / 600 oracle pairs (merged sequence and report); shared `PairWordMatcher` with `seqmatchall`.

- **Audit A4-3 (F153):** the EMBOSS reader drops gap characters and every other non-letter (except `*`) from plain input sequences; one shared `EmbossSequenceReader` (Infrastructure) serves megamerger, seqmatchall: 240 oracle rows each equal (`emboss_gap_oracle.tsv`).
