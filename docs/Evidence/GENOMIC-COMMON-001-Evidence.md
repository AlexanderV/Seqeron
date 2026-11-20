# Evidence Artifact: GENOMIC-COMMON-001

**Test Unit ID:** GENOMIC-COMMON-001
**Algorithm:** Longest Common Substring / Common Region Detection (generalized suffix tree)
**Date Collected:** 2026-06-13

---

## Online Sources

### Wikipedia — "Longest common substring"

**URL:** https://en.wikipedia.org/wiki/Longest_common_substring
**Retrieved:** WebFetch of the URL on 2026-06-13 (also re-fetched the same URL to read the References section).
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia article citing primaries; the primary is Gusfield 1997, reference 2 below).

**Key Extracted Points:**

1. **Formal definition (verbatim):** "Given two strings, S of length m and T of length n, find a longest string which is substring of both S and T."
2. **Distinction from LCSubsequence (verbatim):** "Unlike the longest common subsequence problem, which finds insertions or deletions within the common text, the longest common substring problem seeks a contiguous substring shared by both texts." — i.e. the common region MUST be contiguous.
3. **Suffix-tree solution & complexity (verbatim):** "One can find the lengths and starting positions of the longest common substrings of S and T in Θ(n + m) time with the help of a generalized suffix tree." This Θ(n+m) claim is attributed in the References section to reference 2 (Gusfield 1997).
4. **Worked example (ties):** For S = "BADANAT" and T = "CANADAS", the article states these "share the maximal-length substrings 'ADA' and 'ANA'." → multiple distinct LCS of equal maximal length can co-exist.
5. **Worked example (3 strings, unique):** For "ABABC", "BABCA", "ABCBA" the result is "only one longest common substring, viz. 'ABC' of length 3."

### GeeksforGeeks — "Suffix Tree Application 5 — Longest Common Substring"

**URL:** https://www.geeksforgeeks.org/dsa/suffix-tree-application-5-longest-common-substring-2/
**Retrieved:** WebFetch of the URL on 2026-06-13.
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference-implementation-style description of the generalized-suffix-tree algorithm).

**Key Extracted Points:**

1. **Worked example (verbatim output):** "Longest Common Substring in xabxac and abcabxabcd is: abxa, of length: 4".
2. **GST mechanism (verbatim):** "The path label from root to the deepest node marked as XY will give the LCS of X and Y. The deepest node is highlighted in above figure and path label 'abx' from root to that node is the LCS of X and Y." → the LCS is the deepest internal node whose subtree contains leaves from *both* strings.
3. **Complexity (verbatim):** "If two strings are of size M and N, then Generalized Suffix Tree construction takes O(M+N) and LCS finding is a DFS on tree which is again O(M+N). So overall complexity is linear in time and space."

---

## Documented Corner Cases and Failure Modes

### From Wikipedia "Longest common substring"

1. **Multiple maximal substrings (ties):** "BADANAT"/"CANADAS" share two distinct maximal substrings ("ADA", "ANA"). A correct implementation must pick a deterministic one (or report all); reporting an arbitrary nondeterministic answer is a defect.
2. **No common substring:** When the two strings share no common character, the LCS is the empty string of length 0 (definition: a longest string that is a substring of both; if only the empty string qualifies, length is 0).

### From GeeksforGeeks "Suffix Tree Application 5"

1. **Contiguity:** the returned region is a *contiguous* substring of both inputs (path label of a single deepest common node), not a gapped subsequence.

---

## Test Datasets

### Dataset: GeeksforGeeks worked example (structure-preserving DNA mapping)

**Source:** GeeksforGeeks "Suffix Tree Application 5" — "abxa, of length: 4". The DNA alphabet is A/C/G/T only, so the literal strings cannot be used; the *length and contiguity* property is the transferable fact. DNA cases below are derived from the contiguity definition (Point 1/Point 2 of Wikipedia) and verified by an independent O(n³) brute-force enumeration of all common substrings (not by the repository implementation).

| Parameter | Value |
|-----------|-------|
| Definition used | longest *contiguous* substring common to both strings |
| Verification | independent brute-force `for i,j: if s[i:j] in t` enumeration |

### Dataset: Tie example (Wikipedia "BADANAT"/"CANADAS" → "ADA","ANA"), DNA analogue

**Source:** Wikipedia tie property (two distinct maximal substrings). DNA analogue `CACAGAG` vs `TACATAGAT` shares two distinct maximal length-3 substrings `ACA` and `AGA`; tie-break = "first found in the second string" selects `ACA` (ends earlier in `TACATAGAT`).

| Parameter | Value |
|-----------|-------|
| sequence1 | `CACAGAG` |
| sequence2 | `TACATAGAT` |
| Maximal substrings | `ACA`, `AGA` (both length 3) |
| Selected (first-in-other tie-break) | `ACA` |

---

## B09 review additions (2026-11-18)

- **Chang & Lawler 1994** (Algorithmica 12:327-344), matching statistics: ms(i) = length of the longest prefix of T[i..] occurring in the pattern text, computed for all i in one left-to-right scan with suffix links (WebSearch snippets / Helsinki lecture notes 2013 / Aboelhoda-Kurtz-Ohlebusch 2004 listing; the journal paper is paywalled; the authors' technical report UCB/CSD-91-654 was opened 2026-11-19, see "Primary sources opened" below). Basis of `FindCommonRegions`, whose per-start longest right-maximal match IS ms(i); property used by the implementation: ms(i) ≥ ms(i−1) − 1.
- **MUMmer (Kurtz et al. 2004, Genome Biol 5:R12)** maximal exact matches, used through the canonical `SuffixTree.FindMaximalExactMatches` (B05 F20, MUMmer-identical): the longest match at i lies on a diagonal whose maximal run is a MEM ⇒ ms(i) = max over MEMs with start ≤ i of (end − i) wherever ms(i) ≥ the MEM length threshold.
- **Oracle** (`scripts/analysis/generate_common_region_oracle.py`, no Seqeron code): O(nm) DP over all end cells and a suffix automaton streamed with t (asserted equal) for the LCS (tie = earliest end in t; positions = Python `str.find`); brute-force `in` matching statistics for the regions. 283 pairs (260 small incl. binary alphabet/tandem repeats/empty, 16 medium 300-700 bp with planted (mutated) 3-250 bp shared blocks, 7 literals incl. the M1/M2 examples); the Seqeron code agreed on every LCS length/sequence and every region, and differed only in `PositionInFirst` (first occurrence) before fix F15 (in many random pairs with repeated substrings).
- Timing (Debug build, random 1 Mbp pair): LCS 2.8 s incl. tree construction, `FindCommonRegions(min = 1)` 666 618 regions 3.8 s, `min = 20` 1.7 s; 1 Mbp pair with 1 % substitutions, min 20: 861 267 regions / 131.6 Mchar 4.6 s. Former implementation: 100 kbp pair 2.6 s (quadratic: ≈ 4 min extrapolated at 1 Mbp).

## Assumptions

1. **ASSUMPTION: Tie-break rule** — When several distinct substrings share the maximal length, no authoritative source mandates *which* one to return (Wikipedia reports all; GeeksforGeeks returns one). The repository's `SuffixTree.LongestCommonSubstringInfo` documents and implements "the first one found in 'other' is returned" (`src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Algorithms.cs`, XML doc lines 44–45, 61–62). This is a deterministic, documented choice consistent with returning *a* correct LCS; it does not change which lengths are maximal, only the representative. Tests assert the documented deterministic choice.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindLongestCommonRegion` returns the maximal contiguous shared substring with correct length and 0-based positions in both sequences — Evidence: Wikipedia definition (Point 1), GeeksforGeeks length example (Point 1).
2. **MUST Test:** tie example returns the deterministic first-in-other maximal substring (`CACAGAG`/`TACATAGAT` → `ACA`) — Evidence: Wikipedia tie property; SuffixTree documented tie-break.
3. **MUST Test:** no common substring → `CommonRegion.None` (empty, length 0, positions −1) — Evidence: Wikipedia definition (empty qualifies → length 0).
4. **MUST Test:** identical sequences → whole sequence is the LCS at positions 0/0 — Evidence: definition (the string is a substring of itself).
5. **MUST Test:** `FindCommonRegions(minLength)` enumerates all distinct contiguous shared substrings of length ≥ minLength with their positions — Evidence: definition extended to all common substrings (contiguity).
6. **SHOULD Test:** empty sequence input → `None` / empty enumeration — Rationale: documented empty handling (length-0 LCS).
7. **COULD Test:** invariant — returned substring occurs in both sequences at the reported positions — Rationale: cross-check that positions are consistent with the substring.
8. **COULD Test:** O(n·m) `FindCommonRegions` property — every returned region length ≥ minLength and is contained in sequence1 — Rationale: O(n·m) algorithm property test.

---

## References

1. Charalampopoulos P., Kociumaka T., Pissis S.P., Radoszewski J. (2021). Faster Algorithms for Longest Common Substring. European Symposium on Algorithms (ESA 2021), LIPIcs vol. 204, Schloss Dagstuhl. https://doi.org/10.4230/LIPIcs.ESA.2021.30 (cited in the Wikipedia References section retrieved 2026-06-13).
2. Gusfield, Dan (1997). Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology. Cambridge University Press. ISBN 0-521-58519-8. (Reference 2 in the Wikipedia "Longest common substring" References section; backs the Θ(n+m) generalized-suffix-tree claim.) Citation details retrieved verbatim from https://en.wikipedia.org/wiki/Longest_common_substring on 2026-06-13.
3. Wikipedia contributors. Longest common substring. https://en.wikipedia.org/wiki/Longest_common_substring (accessed 2026-06-13).
4. GeeksforGeeks. Suffix Tree Application 5 — Longest Common Substring. https://www.geeksforgeeks.org/dsa/suffix-tree-application-5-longest-common-substring-2/ (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.

## Primary sources opened (audit A1-22, 2026-11-19)

- **Chang & Lawler** - the journal version (Algorithmica 12:327-344, 1994) is paywalled (Springer HTTP 303 -> login; the CSHL repository entry `repository.cshl.edu/id/eprint/31388` is metadata and abstract only). The authors' technical report of the same work, **UCB/CSD-91-654 "Sublinear Expected Time Approximate String Matching and Biological Applications"** (`www2.eecs.berkeley.edu/Pubs/TechRpts/1991/CSD-91-654.pdf`, 24 pages, image scan), was opened and pages 1, 2, 7, 9 and 10 were read: **section 2.2 "Matching Statistics"** defines `M(T,P)_i = l` if `l` is the length of the longest substring of `P$` matching exactly a prefix of `T[i..n]`, computed "in order of increasing i, in a single left-to-right scan of the text" on the suffix tree of `P$` ("matching letter-at-a-time ... subsequent longest matches are found by following suffix links"), O(n + m) time; the "Matching Statistics Algorithm" with the invariant "`T[i..k-1]` is a word; `M(T,P)_i = k - i`". This is exactly `ms(i)` of `GenomicAnalyzer.MatchingStatistics` (text = sequence 1, query = sequence 2) and the suffix-link stream of `SuffixTreeAlgorithms.AdvanceMatchingStatistics`. The inequality `ms(i) >= ms(i-1) - 1` used by the implementation follows from the suffix-link step (same section).

## 2026-11-19 audit round 2 (A2-30 / A2-35, B09 F94): modes, strands, k mismatches, N sequences, string overloads

- **Source definitions.** MUMmer 3 (Kurtz et al. 2004, Genome Biol 5:R12): `mummer -maxmatch` (all maximal exact matches), `-mumreference` (unique in the reference), `-mum` (unique in both; Delcher et al. 1999 NAR 27:2369), `-r` / `-b` (reverse complement / both) and `-c` ("report the query position of a reverse complement match relative to the forward strand of the query"); the k-mismatch LCS of Flouri, Giaquinta, Kobert & Ukkonen 2015 (IPL 115:643) and Thankachan, Apostolico & Aluru 2016; the k-common substring of Gusfield 1997 §7.6. Opened: the installed `mummer` 3.23 usage text and its output on the worked pair `ACGGATCCGTTTTACGATCCGTAAGT` / `GGATCCGAAAACGGATCCGTTTT`.
- **Oracle numbers.** `scripts/analysis/generate_common_region_strand_oracle.py mum`: 4716 (reference, query, L in 2/3/5, mode) rows; on every ACGT row `mummer -b` (forward and reverse lines, positions converted to forward 0-based starts; the `-c` run asserted to give `p - len` for the same match) equals an independent Python brute force (MEM definition; uniqueness by occurrence counting; reverse = the same on revcomp(query)); 675 rows with N runs use the brute force only. `kmism`: 1680 rows (k = 0..3), brute force over every (i, j, length); `multi`: 694 rows, brute force of the longest substring in >= q of 2-4 sequences. 0 mismatches.
- **Implementation.** `GenomicAnalyzer.CommonRegions.cs`: `FindCommonRegions(…, CommonRegionMode, StrandSelection)` calls the canonical `SuffixTree.FindMaximalExactMatches` / `FindMaximalUniqueMatches` and the existing matching-statistics engine; the reverse strand runs them on the reverse-complemented query and maps `PositionInSecond = m - q - length`. `FindLongestCommonRegionWithMismatches`: per-diagonal two pointers. `FindLongestCommonRegion(IReadOnlyList<…>, minSupport)`: `SuffixTree.FindLongestCommonSubstrings` (minSupport >= 2: support 1 would return the N-run separators). String overloads: non-A/C/G/T runs become unique separators of the reference tree and a non-matching character of the query.
- **Tests.** `GenomicAnalyzer_CommonRegionStrands_Oracle_Tests` (11), MCP `find_common_regions` (`mode`, `strands`), new `find_longest_common_region_mismatches`; N-sequence LCS is covered by the existing `suffix_tree_k_common_substrings` tool of the SuffixTree server.

## 2026-11 audit round 3 (A3-13, B09 F139): EMBOSS seqmatchall

Source: `emboss/seqmatchall.c`, `nucleus/embword.c` (`embWordGetTable`, `embWordBuildMatchTable`; EMBOSS 6.6.0, raw.githubusercontent.com/kimrutherford/EMBOSS) and the installed binary. `GenomicAnalyzer.SeqMatchAll(sequences, wordSize)` reports, for every pair i < j with length(i) > k and length(j) >= k, every maximal exact match of at least k symbols (forward strand): by the canonical suffix-tree engine (`FindCommonRegions` mode `MaximalExactMatches`, `FindMaximalExactMatches`) when the indexed sequence holds no unknown symbol; otherwise by the program's own table scan. **Program quirk reproduced:** `embWordGetTable` skips words with the unknown symbol of the sequence type (N when every character of the indexed sequence is a nucleotide code, else X) but inspects only the symbol entering a window on the right (and the whole first window): after a skipped N the window restarts right behind it and its interior may still contain an N, so such words ARE in the table and match words with the same N positions (rows with 15-25 % N: matches of length 14 and 15 across an N). The query side has no skipping. Oracle `generate_seqmatchall_oracle.py` (330 sets of 2-5 sequences: ACGT / AC / ACGTN / protein, related by point mutations and shared segments, mixed case, k = 2..12; 700+ pairs): matches equal as sets on **330 / 330 sets** (307 non-empty). Differential on N-free nucleotide sets: `SeqMatchAll` = `FindCommonRegions(..., MaximalExactMatches)` per pair on 200 sets.

## 2026-11 audit round 3 (A3-13, B09 F139): EMBOSS palindrome vs `GenomicAnalyzer.FindPalindromes`

The documents cite "EMBOSS einverted: only A/C/G/T pair" and B04 (F10 / F21 / F22) has the full comparison of `palindrome` / `einverted` with `RepeatFinder` (`FindInvertedRepeats`, `maxMismatches`, `FindInvertedRepeatsScored`); B09's facade `FindPalindromes` is the contiguous (no-loop) reverse-complement palindrome of a restriction-site style definition and is NOT the EMBOSS program. Exact relation, from `palindrome -gaplimit 0 -nummismatches 0 -overlap` (`scripts/analysis/generate_palindrome_gap0_oracle.py`, 330 random sequences with planted palindromes, minpallen 2-5, maxpallen 6 / 12 / 20; 1000+ hits): EMBOSS reports, per centre, only the **longest** palindrome of arm >= minpallen, capped to the inner 2 x maxpallen bases; `FindPalindromes(2 m, 2 M)` lists every even palindromic substring of that length range, nested ones included. The per-centre longest member of `FindPalindromes` equals the EMBOSS hit set (start, end, arm) on **330 / 330 sequences**, and `FindPalindromes` returned 50+ extra nested palindromes that EMBOSS never reports. Gap- and mismatch-tolerant inverted repeats (`-gaplimit`, `-nummismatches`, `einverted`) are `RepeatFinder.FindInvertedRepeats` / `FindInvertedRepeatsScored` (B04).

## Audit A3-13b (B09 F149): megamerger

Sources actually opened: EMBOSS 6.6.0 `emboss/megamerger.c` (`megamerger_Merge`), `nucleus/embword.c` (`embWordGetTable`, `embWordBuildMatchTable` with `orderit`, `wordMatchCmp`, `wordMatchCmpPos`, `embWordMatchMin`, `wordDeadZone`; the word comparison / hash are case-insensitive), `megamerger.acd`, installed `megamerger` 6.6.0 (raw.githubusercontent.com/kimrutherford/EMBOSS, User-Agent `Seqeron-review/1.0`).

`GenomicAnalyzer.MegaMerge` + MCP `megamerger`: the word matches are those of `SeqMatchAll` (the per-pair engine is now the shared `PairWordMatcher`; the seqmatchall oracle still equal), ordered by `wordMatchCmp` (longest, start1, start2), reduced by `embWordMatchMin` (pop the longest, `wordDeadZone`: live 0 / dead 1 / truncated 2, a truncated match shorter than the word is dead, re-sort after a truncation, final order by position), merged by `megamerger_Merge` with the program's unsigned arithmetic (`AJMIN(mid, len - mid - 1)`) and its `ajStrAssignSubS` clipping (`EmbossText.Substring`, also used by showalign: a reversed or out-of-range window yields one character). Measured on the binary: a second sequence shorter than the word only warns and writes nothing; a first sequence shorter than the word is "No match found" (exit 1); matched regions are lower case, regions chosen over a mismatch or an unmatched end upper case.

Oracle `megamerger_oracle.jsonl` (`generate_megamerger_oracle.py`): 600 pairs (a shared core of 60-300 bp with substitutions / indels / N / IUPAC letters between the copies, flanks of random length or none, lower-case stretches, unrelated pairs; word size 2-16; `-prefer`): **538 merged sequences + reports equal, 62 pairs on which the program writes nothing equal (empty merge), 0 aborts**; every report branch is exercised (mismatch regions, both "neither overlap is at the start / end" warnings, both initial and final sequence choices, both "closer to the ends" choices). Tests: `GenomicAnalyzer_MegaMerge_Tests` (3), MCP `megamerger`.

## Audit A4-3 (B09 F153): EMBOSS reader - gap characters and other non-letters

Sources actually opened: EMBOSS `ajax/core/ajseqread.c` (`seqAppend`: `ajStrKeepSetAlphaC(.., "*.~?#+-")`), `ajseqtype.c` (`ajSeqTypeCheckIn`, `seqTypeFix`: an ungapped type calls `ajStrRemoveGap`; `seqType[]` table: every ACD type of the programs is ungapped, `gapany` for consambig / showalign), `ajstr.c` (`ajStrRemoveGap` keeps `A-Za-z*`), `ajseq.c` (`ajSeqIsNuc`, `ajSeqGapStandard`) from raw.githubusercontent.com/kimrutherford/EMBOSS (User-Agent `Seqeron-review/1.0`), the installed ACD files (`/usr/share/EMBOSS/acd`: `type:` of every sequence input) and the installed programs. Measured with `scripts/alignment/measure_emboss_reader_chars.py` (each of the characters 1..255 except newline / carriage return / `>` inserted into the middle of a sequence, both sequences, nucleotide and protein): for megamerger, seqmatchall the output equals the output without the character for exactly 199 characters (all except the 52 ASCII letters and `*`); `*` is kept (a protein stop; an error in a nucleotide-typed input).

Oracle `emboss_gap_oracle.tsv` (`generate_emboss_gap_oracle.py`, seed 20261126): 18 programs x 240 rows, megamerger (240) and seqmatchall (240): **equal on every row**. Audit example: `needle AC-GTTAGC.TT* ACGTTAGCTT` prints `# Length: 11`, `# Score: 58.0` (the port scored 43 on the uncleaned text).

## Audit A4-1 (B09 F157): cancellation

Mechanism: `AlignmentCancellation` (SequenceAligner.Cancellation.cs, A2-40); tests `SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `MegaMerge` (`MinimalMatchSet`) poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

## Audit A4-2 (B09 F158): resource guards

`SequenceAligner_ResourceLimits_Tests`: `SeqMatchAll` / `MegaMerge` (`WordMatchMaxMatches`) - documented limit constants (`SequenceAligner.ResourceLimits.cs`), the exception of the `SequenceAligner.ResourceLimit` style raised before the allocation in 64-bit arithmetic, and a low-allocation assertion per limit.
The programs have no such limits (they die at the allocation: measured OutOfMemory / SIGKILL on the audit inputs); minimap2 drops frequent minimizers (`mid_occ`, `mm_seed_select`), ported in A4-7 (the anchor-count guard now only fires when the minimizers that survive that filter are still above `MinimapMaxAnchors`).

## Audit A4-9 (B09 F164): megamerger edge cases

`scripts/analysis/generate_megamerger_edge_oracle.py` -> `TestData/Analysis/megamerger_edge_oracle.jsonl`, **480 rows** of the installed megamerger 6.6.0: 120 sequences entirely in lower case, 120 with one lower- and one upper-case sequence, 120 with runs of N longer than the word size (the word matcher skips N runs: they never match each other) and 120 short / equal-to-the-word / identical / nested sequences (status ok 3xx, nomerge, abort = the program dies, the port throws its `ArgumentException`). `GenomicAnalyzer_MegaMerge_Tests.MegaMerge_EqualsMegamerger_OnTheEdgeCaseRows`: merged sequence and report equal on 480 / 480 (no code change needed).
