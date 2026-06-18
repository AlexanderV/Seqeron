# Validation Report: GENOMIC-REPEAT-001 — Repeat Detection (LRS + all repeats)

- **Validated:** 2026-06-15   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindLongestRepeat(DnaSequence)`, `GenomicAnalyzer.FindRepeats(DnaSequence, int minLength)`
- **Stage A verdict:** PASS
- **Stage B verdict:** FAIL → FIXED (defect found in `FindRepeats`, fully corrected this session)

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
