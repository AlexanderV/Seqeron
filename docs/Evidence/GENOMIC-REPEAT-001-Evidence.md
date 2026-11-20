# Evidence Artifact: GENOMIC-REPEAT-001

**Test Unit ID:** GENOMIC-REPEAT-001
**Algorithm:** Repeat Detection — Longest Repeated Substring (LRS) and all repeated substrings via suffix tree
**Date Collected:** 2026-06-13

---

## Online Sources

### Carnegie Mellon University 15-451/651 (Algorithm Design & Analysis), Lecture #10: Suffix Trees and Arrays

**URL:** https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-instruction university course notes, textbook-level)
**Retrieval:** WebSearch "longest repeated substring suffix tree deepest internal node algorithm Gusfield" → fetched the PDF URL; WebFetch returned only binary, so the saved PDF was converted with `pdftotext` and grepped.

**Key Extracted Points:**

1. **Longest-repeat definition + algorithm (verbatim, §2.1):** "Find the longest repeat in T . That is, find the longest string r such [that] r occurs at least twice in T : Find the deepest node that has ≥ 2 leaves under it." (lines 199–200 of the extracted text). "Deepest" is in the sense of the longest string spelled from the root.
2. **Suffix-tree node structure:** "each internal node now has degree at least 2, hence the total number of nodes in the tree is at most twice the number of leaves" and "each leaf corresponds to some suffix of T". This is why a substring that occurs ≥ 2 times corresponds to an internal node (≥ 2 leaves below), and a substring occurring once corresponds to a leaf.
3. **Linear construction:** The suffix tree representation "uses O(t) space" for text length t; the longest-repeat query is a single traversal, so the overall solution is linear.

### Wikipedia — Longest repeated substring problem

**URL:** https://en.wikipedia.org/wiki/Longest_repeated_substring_problem
**Accessed:** 2026-06-13 · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (encyclopedic; used for the worked example and to corroborate CMU)
**Retrieval:** WebSearch (same query as above) → WebFetch of the article URL.

**Key Extracted Points:**

1. **Problem statement:** "The problem of finding the longest substring of a string that occurs at least twice."
2. **Algorithm:** Build a suffix tree with an end-of-string marker (e.g. `$`) and find "the deepest internal node in the tree with more than one child"; the path from root to that node spells a longest repeated substring. Depth = number of characters traversed from the root.
3. **Complexity:** Linear time and space, Θ(n).
4. **Worked example (verbatim):** Input `ATCGATCGA$` → longest repeated substring `ATCGA`.

### GeeksforGeeks — Suffix Tree Application 3: Longest Repeated Substring

**URL:** https://www.geeksforgeeks.org/dsa/suffix-tree-application-3-longest-repeated-substring/
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation walk-through; used for additional worked examples)
**Retrieval:** WebSearch "longest repeated substring suffix tree deepest internal node algorithm Gusfield" → WebFetch of the article URL.

**Key Extracted Points:**

1. **Selection rule (verbatim):** "longest repeated substring will end at the internal node which is farthest from the root (i.e. deepest node in the tree), because length of substring is the path label length from root to that internal node." A repeated substring's path "terminates at an internal node (a node with multiple children), not a leaf."
2. **Worked examples (input → LRS):**
   - `GEEKSFORGEEKS` → `GEEKS`
   - `AAAAAAAAAA` → `AAAAAAAAA`
   - `ABCDEFG` → No repeated substring
   - `ABABABA` → `ABABA`
   - `banana` → `ana`
3. **Complexity:** Building the suffix tree (Ukkonen) is O(N); finding the deepest node is O(N); overall linear.

### Johns Hopkins University (Ben Langmead), Lecture: Suffix Trees

**URL:** https://www.cs.jhu.edu/~langmea/resources/lecture_notes/08_suffix_trees_v2.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (university bioinformatics course notes, explicitly tracking Gusfield)
**Retrieval:** WebSearch "maximal repeat suffix tree bioinformatics definition Gusfield algorithms on strings" → fetched PDF, converted with `pdftotext` and grepped.

**Key Extracted Points:**

1. **Gusfield grounding:** The notes' suffix-tree applications cite "Gusfield 5.4" (Gusfield, *Algorithms on Strings, Trees, and Sequences*, ch. 5–7) as the source for the suffix-tree application family that includes repeat finding.
2. **Internal-node ⇒ shared substring:** The notes use the same characterization for the related longest-common-substring application — every internal node's path label is a substring that recurs; leaves are unique-position suffixes — confirming the "internal node ⇔ repeated substring" mapping used by the LRS algorithm.
3. **Maximal repeats (search-result summary, corroborated by the notes' references to REPuter/Kurtz):** A maximal repeat's path label is an internal node, and there are at most n maximal repeats in a string of length n.

---

## Documented Corner Cases and Failure Modes

### From CMU 15-451 §2.1 and GeeksforGeeks

1. **No repeat exists:** When every substring occurs at most once (no internal node has ≥ 2 leaves), there is no longest repeated substring. GeeksforGeeks: `ABCDEFG` → "No repeated substring". Expected output for this implementation: `RepeatInfo.None` (empty sequence).
2. **Overlapping repeats are allowed:** The definition counts occurrences regardless of overlap. GeeksforGeeks `AAAAAAAAAA` (length 10) → `AAAAAAAAA` (length 9), whose two occurrences at positions 0 and 1 overlap; and `ABABABA` → `ABABA` (occurrences at 0 and 2 overlap). The suffix-tree definition (≥ 2 leaves under a node) does not require non-overlapping occurrences.
3. **Ties:** When several substrings share the maximum length, the problem statement asks for *a* longest repeated substring; any one of equal length is correct (Wikipedia: "a longest repeated substring").

---

## Test Datasets

### Dataset: Wikipedia LRS example

**Source:** Wikipedia, Longest repeated substring problem (accessed 2026-06-13).

| Parameter | Value |
|-----------|-------|
| Input | `ATCGATCGA` (Wikipedia appends terminal `$`; the repo builds the tree with its own terminal internally) |
| Longest repeated substring | `ATCGA` |
| Length | 5 |
| Occurrence count | 2 |
| Occurrence start positions (0-based) | 0, 4 |

### Dataset: GeeksforGeeks worked examples (DNA-alphabet-compatible subset and analogs)

**Source:** GeeksforGeeks, Suffix Tree Application 3 (accessed 2026-06-13).

| Input (as published) | DNA analog used in tests | LRS | Length | Count | Positions |
|----------------------|--------------------------|-----|--------|-------|-----------|
| `AAAAAAAAAA` | `AAAAAAAAAA` (valid DNA) | `AAAAAAAAA` | 9 | 2 | 0, 1 |
| `ABCDEFG` (no repeat) | `ACGT` (no repeat) | (none) | 0 | 0 | — |
| `ABABABA` | `ATATATA` | `ATATA` | 5 | 2 | 0, 2 |

> Rationale for DNA analogs: `GenomicAnalyzer.FindLongestRepeat` takes a `DnaSequence`, whose alphabet is restricted to A/C/G/T (`DnaSequence.ValidateSequence`). The non-DNA GeeksforGeeks strings (`GEEKSFORGEEKS`, `banana`, `ABCDEFG`, `ABABABA`) are mapped to A/C/G/T strings with the identical repeat structure; the structural property (deepest internal node) is alphabet-independent, so the expected LRS, length, count, and overlap pattern are preserved.

### Dataset: FindRepeats enumeration (all repeated substrings ≥ minLength)

**Source:** Definition (CMU §2.1: a repeat is *any* substring occurring ≥ 2 times). Ground truth is an independent brute-force enumeration of **every** distinct substring of `ACGTACGTTTTTACGT` of length ≥ minLength that occurs ≥ 2 times (all `s[i:j]`, keep those with ≥ 2 occurrences), cross-checked against the sorted-suffix-LCP prefix set. *(Updated 2026-06-15: the earlier 5-entry row only listed the maximal-length LCP of each adjacent suffix pair and was a code echo of a defective implementation; it omitted the shorter repeated prefixes `ACG`, `TAC`, `TACG`. See report GENOMIC-REPEAT-001.md / FINDINGS_REGISTER.)*

| minLength | Repeated substrings (occurrences ≥ 2), with positions |
|-----------|--------------------------------------------------------|
| 3 | `ACG` @ {0,4,12}; `ACGT` @ {0,4,12}; `CGT` @ {1,5,13}; `TAC` @ {3,11}; `TACG` @ {3,11}; `TACGT` @ {3,11}; `TTT` @ {7,8,9}; `TTTT` @ {7,8} |

---

## Assumptions

1. **ASSUMPTION: Tie-breaking among equal-length longest repeats** — When more than one substring shares the maximum repeated length, the authoritative sources require only "a longest repeated substring" (Wikipedia), not a specific one. The repository returns whichever the suffix tree's deepest-internal-node bookkeeping records. Tests therefore avoid asserting a unique winner where the worked example has multiple equal-length maximal repeats; all cited test inputs have a single longest repeat, so this assumption is not exercised by any MUST case.
2. **ASSUMPTION: Occurrence ordering** — `RepeatInfo.Positions` is sorted ascending by 0-based start index. This is an output-shape convention (the implementation calls `OrderBy(p => p)`); the *set* of positions is fixed by the definition, only their listing order is conventional.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindLongestRepeat("ATCGATCGA")` returns sequence `ATCGA`, length 5, count 2, positions {0,4}. — Evidence: Wikipedia worked example.
2. **MUST Test:** `FindLongestRepeat("AAAAAAAAAA")` returns `AAAAAAAAA` (overlapping, length 9, count 2, positions {0,1}). — Evidence: GeeksforGeeks.
3. **MUST Test:** `FindLongestRepeat("ATATATA")` returns `ATATA` (overlapping, length 5, positions {0,2}). — Evidence: GeeksforGeeks `ABABABA`→`ABABA` analog.
4. **MUST Test:** `FindLongestRepeat("ACGT")` (no repeat) returns `RepeatInfo.None`/IsEmpty. — Evidence: GeeksforGeeks `ABCDEFG`→no repeat.
5. **MUST Test:** `FindLongestRepeat("")` (empty) returns `RepeatInfo.None`. — Evidence: definition (no substring occurs twice in the empty string).
6. **MUST Test:** `FindRepeats("ACGTACGTTTTTACGT", 3)` yields exactly {ACG@{0,4,12}, ACGT@{0,4,12}, CGT@{1,5,13}, TAC@{3,11}, TACG@{3,11}, TACGT@{3,11}, TTT@{7,8,9}, TTTT@{7,8}} (8 substrings); every result has count ≥ 2 and length ≥ minLength. — Evidence: definition (every substring occurring ≥ 2 times); brute-force ground truth.
7. **SHOULD Test:** `FindRepeats` returns empty when minLength exceeds any repeat length, and only repeats meeting minLength are returned. — Rationale: minLength filter boundary.
8. **SHOULD Test:** Every `RepeatInfo` from both methods satisfies the invariants (count ≥ 2; length = sequence length; positions sorted; each position is a true occurrence). — Rationale: invariant guard (property test for the O(n²)-or-worse enumeration path).
9. **COULD Test:** `FindRepeats` with minLength ≤ 0 still returns only substrings occurring ≥ 2 times (no zero-length repeats). — Rationale: degenerate parameter.

---

## References

1. Sleator, D. et al. (course staff). 2017. *15-451/651 Algorithm Design and Analysis, Lecture #10: Suffix Trees and Arrays.* Carnegie Mellon University. https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf
2. Wikipedia contributors. 2026. *Longest repeated substring problem.* Wikipedia. https://en.wikipedia.org/wiki/Longest_repeated_substring_problem
3. GeeksforGeeks. *Suffix Tree Application 3 – Longest Repeated Substring.* https://www.geeksforgeeks.org/dsa/suffix-tree-application-3-longest-repeated-substring/
4. Langmead, B. *Suffix Trees (lecture notes).* Johns Hopkins University. https://www.cs.jhu.edu/~langmea/resources/lecture_notes/08_suffix_trees_v2.pdf
5. Gusfield, D. 1997. *Algorithms on Strings, Trees, and Sequences: Computer Science and Computational Biology* (ch. 5–7, suffix-tree applications; cited as "Gusfield 5.4" by ref. 4). Cambridge University Press. https://doi.org/10.1017/CBO9780511574931 *(not retrieved in full text this session; cited only as the upstream source that refs. 1 and 4 attribute the suffix-tree repeat applications to — no formula or expected value is taken from it directly).*

---

## B09 Review Additions (2026-11-18)

**Sources actually opened (2026-11-18; superseded by "Primary sources opened" at the end of this file, 2026-11-19):** WebSearch snippets (no page fetched) restating Gusfield 1997 §7.12 — UCF CAP 5937 "Applications of Suffix Trees" (http://www.cs.ucf.edu/courses/cap5937/fall2004/Applications%20of%20suffix%20trees.pdf) and Abouelhoda/Kurtz/Ohlebusch 2004 (J. Discrete Algorithms 2:53–86): left character of a position i is S(i−1); an internal node is *left-diverse* when its subtree has leaves with ≥ 2 distinct left characters; a string is a **maximal repeat iff its suffix-tree node is left-diverse** (internal nodes are right-maximal by construction); a *supermaximal* repeat is a maximal repeat not occurring inside another maximal repeat. The book text was not retrievable. Biopython 1.88 (`Seq.reverse_complement`) installed locally; Rosalind REVP sample reused from `RepeatFinder_Palindrome_Tests` (rosalind.info is blocked).

**Oracle numbers.**
- `scripts/analysis/generate_repeat_oracle.py`: 628 (sequence, minLength) rows over 157 sequences (the worked examples above, homopolymers, Fibonacci(6/8), Thue–Morse(40), `ACGT`×6, planted repeats, random ACGT/AC/ACG of length 3–45) × minLength {1,2,3,5}; repeated substrings, maximal repeats and all longest repeats computed by (a) brute force over all substrings and (b) suffix array + Kasai LCP + lcp-interval left characters — asserted equal on every row; the C# methods match all 628 rows (count + SHA-256 prefix of the ordered `seq:positions` list; explicit maximal and longest lists). Worked values: `ATCGATCGA` minLength 1 → 14 repeated substrings, maximal repeats `A@{0,4,8}`, `ATCGA@{0,4}`; `ACGTACGT` → 10 repeated substrings, 1 maximal (`ACGT@{0,4}`); `TTGGCACGTAACGTCTTGG` → ties `TTGG@{0,15}`, `ACGT@{5,10}`; `AAAAA` → maximal `A, AA, AAA, AAAA` (text ends are diverse).
- In-test brute force: 480 random strings (4 alphabets) × 4 minLengths incl. order; maximal repeats = strings of `RepeatFinder.FindDirectRepeats` pairs (no limits) on 240 random strings.
- 1 Mbp random DNA + planted 200-mer, minLength 20: `FindRepeats` (16 653 repeats) + `FindMaximalRepeats` (2) in 4.7 s including the suffix-tree build.
- `scripts/analysis/generate_palindrome_oracle.py` (Biopython `reverse_complement` over every window, == pure-str translate/reverse): 910 rows (130 sequences × 7 (min,max) ranges), 0 mismatches against `GenomicAnalyzer.FindPalindromes`; REVP sample `4 6, 5 4, 6 6, 7 4, 17 4, 18 4, 20 6, 21 4` (1-based) reproduced.

**Corrected claims.** The MCP page `find_longest_repeat.md` / `.mcp.json` example claimed `ATGATGATG` → `ATG@[0,3,6]`; the longest repeated substring is `ATGATG@[0,3]` (length 6; the unit test already asserted this). The `find_repeats` page claimed O(n²) adjacent-suffix enumeration; the implementation is a single tree traversal, O(n + output).

**Assumption 1 retired:** tie-break is now deterministic (leftmost first occurrence).

---

## Change History

- **2026-11-18 (B09 A1-21)**: string overload `FindPalindromes(string, …)` (canonical raw-string `RepeatFinder` contract: only A/C/G/T pair); oracle `generate_palindrome_oracle.py --ambiguous` (861 rows, Biopython `reverse_complement` on A/C/G/T-only windows).
- **2026-06-13**: Initial documentation.
- **2026-11-18**: B09 review additions (maximal repeats, tie-break, palindrome facade, oracles).

## Primary sources opened (audit A1-22, 2026-11-19)

- **Gusfield 1997 §7.12.1** - the book itself is not retrievable. Opened instead: the UCF CAP 5937 lecture slides "Applications of Suffix Trees" (`www.cs.ucf.edu/courses/cap5937/fall2004/Applications%20of%20suffix%20trees.pdf`, 14 pages, `pdftotext`), a chapter-by-chapter summary of the book: "An internal node v in T(S) is called **left-diverse** if at least two leaves in v's subtree have different left characters"; "Theorem [Gusfield, page 144]: the string α labeling the path to a node v of T(S) is a **maximal repeat if and only if v is left diverse**"; "a supermaximal repeat is a maximal repeat that never occurs as a substring of any other maximal repeat"; the linear-time bottom-up computation of left-diversity ("if any child of v is left diverse, then v is left diverse"). This is the definition and algorithm of `FindMaximalRepeats` (`RepeatLocusCollector`: a bit set of left characters per node, a sentinel for suffix 0). Still a secondary source for the page-144 theorem; no primary text was reachable (Cambridge University Press / Google Books not fetched).

## 2026-11-19 audit round 2 (A2-33 / A2-35, B09 F97 / F99): both-strand repeats, repeat pairs, string overloads

- **Source definitions.** MUMmer `repeat-match` (Kurtz et al. 2004; default both strands, `-f` forward only, `r` lines), Vmatch `-p` and REPuter palindromic repeats (Kurtz & Schleiermacher 1999), Gusfield 1997 §7.12 (maximal pairs / repeats) applied to `S # revcomp(S)`. The installed `repeat-match` 3.23 was run; its quirk (S and revcomp(S) share the final terminator, so a palindromic pair such as (0,0,2) of `ATAAT` is lost and an exact-palindrome genome aborts) is the one B04 F25 found — the oracle runs it on `S + 'N'`.
- **Oracle numbers.** `scripts/analysis/generate_strand_repeat_oracle.py`: 1900 (sequence, L = 1..4) rows over random / planted ACGT, AC, AT sequences (950 with N runs): the maximal and all-substring classes of an independent Python brute force over `S # revcomp(S)` (unique separators for non-ACGT), and on the ACGT rows the `repeat-match` pair lists (both strands and `-f`); the script asserts that the class set of the `repeat-match` pairs equals the brute-force maximal classes on every ACGT row (0 disagreements; one palindromic quirk class removed by the sentinel). The C# results equal all columns.
- **Implementation.** `GenomicAnalyzer.Repeats.cs`: one suffix tree of `S # revcomp(S)` per call, traversed by the same `EnumerateRepeatsCore` (generalised to any suffix tree; separator characters above `Z` never form a repeat and count as a distinct left context) that serves the forward search; `FindRepeatPairs` unites `RepeatFinder.FindDirectRepeats` / `FindReverseComplementRepeats` with the filters lifted; string overloads build the tree of the masked text (`MaskText`).
- **Tests.** `GenomicAnalyzer_StrandRepeats_Oracle_Tests` (13: 1900-row class oracle, all-substring oracle, `DnaSequence` parity, `repeat-match` pairs (both / `-f`), pair strings == classes, longest classes, palindromic site, reverse-complement copy, in-test brute force of the string overloads on 600 N / IUPAC strings, parity on ACGT, N runs, validation); MCP `find_repeats` (`bothStrands`, IUPAC), new `find_all_longest_repeats`.

## 2026-11 audit round 3 (A3-13, B09 F139): `FindPalindromes` vs EMBOSS palindrome

See the GENOMIC-COMMON-001 Evidence file section "EMBOSS palindrome vs `GenomicAnalyzer.FindPalindromes`": per-centre longest palindromes of `FindPalindromes(2m, 2M)` equal `palindrome -gaplimit 0 -nummismatches 0` on 330 / 330 sequences; the facade additionally lists nested palindromes; gap / mismatch-tolerant inverted repeats are the B04 `RepeatFinder` methods.
