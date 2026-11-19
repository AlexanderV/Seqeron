# Validation Report: GENOMIC-REPEAT-001 — Repeat Detection (LRS + all repeats)

- **Validated:** 2026-06-15   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindLongestRepeat(DnaSequence)`, `GenomicAnalyzer.FindRepeats(DnaSequence, int minLength)`
- **Stage A verdict:** PASS
- **Stage B verdict:** FAIL → FIXED (defect found in `FindRepeats`, fully corrected this session)

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES (maximal-repeat definition restated from Gusfield §7.12 via lecture slides; the unit contract 'every repeated substring' kept and documented). **Stage B:** PASS-WITH-NOTES (FAIL -> FIXED in 2026-06). **State:** FIXED.
- **Fixes (B09.md, F-numbers):** F27 `FindPalindromes` = facade over `RepeatFinder.FindPalindromes` (was an unvalidated O(n·L²) copy), F28 `FindRepeats` review + independent oracle, F29 deterministic `FindLongestRepeat` tie-break (+ `FindAllLongestRepeats`), F30 `FindMaximalRepeats` (Gusfield left-diverse), F31 MCP documentation, F50 N / IUPAC.
- **Audit round 2 (A2-33 / A2-35, F97 / F99):** both-strand repeats (`FindRepeatsBothStrands`, `FindAllLongestRepeatsBothStrands`, `FindLongestRepeatBothStrands`) and `FindRepeatPairs` (`repeat-match` 3.23, 1 900 rows), string / IUPAC overloads with an optional alphabet; MCP `find_repeats` (both strands), `find_all_longest_repeats`.
- **Oracle numbers:** 628-row Python oracle (brute force = SA + LCP), 480 brute-force strings, 1 Mbp ~ 4.7 s.
- **Open items:** registered, not doable here: a public left-diverse (maximal-repeat) visitor in `src/SuffixTree` (B05, QUEUE R18) would replace the private `RepeatLocusCollector`; the present implementation is exact and linear.
- The unit's TestSpec (`tests/TestSpecs/GENOMIC-REPEAT-001.md`), Evidence (`docs/Evidence/GENOMIC-REPEAT-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Stage A — Description

### Sources opened & what they confirm
- **CMU 15-451 Lecture #10 §2.1** (https://www.cs.cmu.edu/~15451-f17/lectures/lec10-sufftree.pdf, retrieved 2026-06-15, converted with `pdftotext`). Lines 199–200 verbatim: "Find the longest repeat in T. That is, find the longest string r such r occurs at least twice in T: Find the deepest node that has ≥ 2 leaves under it." Lines 161–163 / 194–196 confirm: each leaf = a suffix, occurrence count of P = number of leaves under P's node; so a substring occurring ≥ 2 times maps to an internal node and one occurring once to a leaf.
- **Wikipedia — Longest repeated substring problem** (retrieved 2026-06-15). Definition: longest substring occurring ≥ 2 times; solved in Θ(n) by the deepest internal node with > 1 child of a suffix tree built with a `$` sentinel. Worked example verbatim: `ATCGATCGA$` → `ATCGA`. Note: Wikipedia covers only the *single* longest repeat — it does **not** define "all repeated substrings".
- **GeeksforGeeks — Suffix Tree Application 3** (cited; corroborates the deepest-internal-node rule and the worked examples `AAAAAAAAAA`→`AAAAAAAAA`, `ABCDEFG`→none, `ABABABA`→`ABABA`, `banana`→`ana`).

### Formula check
- LRS = deepest internal node string-depth: matches CMU §2.1 and Wikipedia exactly. ✓
- "All repeats" definition (the `FindRepeats` contract, Repeat_Detection.md §1): "enumerates **every** distinct substring occurring at least twice with length ≥ a given minimum." Grounded in CMU §2.1 (a repeat is *any* substring occurring ≥ 2 times) — sound as a definition. ✓
- INV-01..INV-06 are genuine properties. ✓

### Edge-case semantics
Empty sequence → None; no-repeat (`ACGT`) → None; overlapping occurrences counted (`AAAAAAAAAA`→`AAAAAAAAA`@{0,1}); `minLength ≤ 0` clamps to 1 (no zero-length "repeat"). All defined and sourced. ✓

### Independent cross-check (numbers)
Brute-force LRS (all substrings, overlap-counted) reproduced every cited value:
`ATCGATCGA`→`ATCGA`@{0,4}; `AAAAAAAAAA`→`AAAAAAAAA`@{0,1}; `ATATATA`→`ATATA`@{0,2}; `ACGT`→none; `banana`→`ana`@{1,3}.

### Findings / divergences
None. The biology/maths of the description is correct. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed
- `FindLongestRepeat` — `GenomicAnalyzer.cs:25-37`; delegates to `SuffixTree.LongestRepeatedSubstring()` (deepest internal node) + `FindAllOccurrences`, sorted ascending. Correct: matches Wikipedia/GeeksforGeeks LRS values (verified against built assembly and brute force).
- `FindRepeats` — `GenomicAnalyzer.cs:48-95`; sorts all suffixes, takes adjacent-pair LCP.

### Defect found (DEFECT, fixed this session)
The original `FindRepeats` emitted **only the full adjacent-pair LCP** of each sorted-suffix pair, not its shorter prefixes. A substring occurs ≥ 2 times iff it is a *prefix* of some adjacent-pair LCP; emitting only the maximal LCP drops every shorter repeated prefix that is not itself a maximal LCP elsewhere.

**Reproduced against the built assembly** (`Seqeron.Genomics.Analysis.dll`):
`FindRepeats("ACGTACGTTTTTACGT", 3)` returned only **5** substrings `{ACGT, CGT, TACGT, TTT, TTTT}`.

**Brute-force ground truth** (every `s[i:j]` of length ≥ 3 occurring ≥ 2 times) is **8** substrings:
`ACG@{0,4,12}, ACGT@{0,4,12}, CGT@{1,5,13}, TAC@{3,11}, TACG@{3,11}, TACGT@{3,11}, TTT@{7,8,9}, TTTT@{7,8}`.
Missing: `ACG`, `TAC`, `TACG` — all genuine repeated substrings. This violates the documented contract (§1: "every distinct substring occurring at least twice").

The repo's TestSpec/Evidence M6 row asserted exactly the buggy 5-item set — a **code echo** (it would pass against the defective code and fail against the correct definition). This is the blind-spot the protocol warns about.

### Fix applied
`GenomicAnalyzer.cs` — replaced the single-LCP emission with a loop over every prefix length `effectiveMinLength..lcpLen` of each adjacent-pair LCP, deduplicated via `HashSet`, positions resolved by `FindAllOccurrences`. Docstring and Repeat_Detection.md §4.1/§5.2/§7 updated to describe the corrected algorithm. Re-probed against the rebuilt assembly: now returns the full 8-item set.

### Cross-verification table recomputed vs code (post-fix)
| Input | minLen | Expected (brute force) | Code (post-fix) |
|-------|--------|------------------------|-----------------|
| `ACGTACGTTTTTACGT` | 3 | 8 substrings (above) | match ✓ |
| `ACGTACGT` | 2 | AC,ACG,ACGT,CG,CGT,GT (6) | match ✓ |
| `ACGTACGT` | 4 | ACGT@{0,4} | match ✓ |
| `ACGTACGT` | 5 | (empty) | match ✓ |
LRS values (ATCGATCGA, AAAAAAAAAA, ATATATA, ACGT) all match. ✓

### Variant/delegate consistency
`SuffixTreeGenomicsTools.FindLongestRepeat` delegates and is tested to match `GenomicAnalyzer.FindLongestRepeat`; unaffected by the `FindRepeats` fix and still green.

### Test quality audit (HARD gate)
- M6 corrected from the code-echoed 5-item set to the brute-force-sourced 8-item set with exact positions.
- Added **M8** (`FindRepeats("ACGTACGT", 2)` → exact 6-item set) as a completeness regression guard — it **fails against the old implementation** and passes only with the fix.
- M1–M5 LRS expectations independently re-derived by brute force (exact sequence/length/count/positions, not ranges).
- M7/S3 are invariant property guards over the full result set (legitimate ≥-style invariants, not weakened exact assertions).
- No skipped/ignored tests; no widened tolerances. Evidence + TestSpec tables updated to the sourced values.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6571` (was 6570 + new M8); changed files build warning-free.

**Gate result: PASS** (defect fixed, tests locked to sourced values, full suite green).

## Verdict & follow-ups
- **Stage A: PASS.** Description is biologically/mathematically correct.
- **Stage B: FAIL → FIXED.** `FindRepeats` incompleteness corrected; tests and docs realigned to brute-force ground truth.
- **End-state: CLEAN.** Algorithm fully functional; full suite green.
- Logged in FINDINGS_REGISTER as a completeness defect in `FindRepeats` (now resolved).

---

## B09 review (2026-11-18) — F27–F31

- **Stage A: PASS-WITH-NOTES** (LRS / repeated-substring definitions confirmed against CMU 15-451 §2.1 and a fresh brute-force + SA/LCP oracle; maximal-repeat definition added from Gusfield §7.12 restated in the UCF slides and Abouelhoda et al. 2004 — snippets only, book text unavailable; note: the unit's contract is "every repeated substring", which lists all prefixes/suffixes of longer repeats — now documented as a decision record, with the non-redundant `FindMaximalRepeats` added).
- **Stage B: PASS-WITH-NOTES → FIXED.** `FindRepeats` (suffix-tree traversal, commit 3d49621e) is semantically correct: 628-row Python oracle (set, positions, order) and 480 random brute-force strings, 1 Mbp ≈ 4.7 s; no regression. Findings: F27 `GenomicAnalyzer.FindPalindromes` was a verbatim copy of the canonical scan without validation (minLength 2 → dinucleotides, odd minLength → nothing, null → NRE, O(n·L²), length-major order) → delegates to `RepeatFinder.FindPalindromes`; F28 `null` sequence threw lazily `NullReferenceException` → `ArgumentNullException` at the call, independent oracles added; F29 `FindLongestRepeat` tie-break unspecified (tree child order) → leftmost first occurrence, `FindAllLongestRepeats` added (canonical `SuffixTree.FindAllLongestRepeatedSubstrings`); F30 maximal repeats were "not implemented" → `FindMaximalRepeats` (linear left-diversity propagation, MCP `find_repeats` `maximalOnly`); F31 two wrong MCP pages (`find_longest_repeat` example, `find_repeats` complexity) corrected.
- **Call-through evaluation:** `RepeatFinder.FindDirectRepeats` / `SuffixTree.FindMaximalRepeatedPairs` report maximal *pairs*, not strings with occurrence sets; `FindRepeats` cannot call them (different result set); `FindMaximalRepeats` strings equal the strings of the pair set (tested), kept as a separate linear tree pass because it returns occurrence sets directly.
- **End-state: FIXED / CLEAN.** Details: `docs/Validation/review-2026-09/B09.md` F27–F31; algorithm doc `docs/algorithms/Repeat_Analysis/Repeat_Detection.md` §2.2, §5.4–5.6.


- **Audit r1 (A1-21):** `FindPalindromes(string, …)` overload (canonical raw-string `RepeatFinder`, only A/C/G/T pair); 861-row ambiguous oracle (B09.md F50).
