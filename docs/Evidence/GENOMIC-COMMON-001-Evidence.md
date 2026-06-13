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
