# Validation Report: PAT-EXACT-001 — Exact Pattern Search (Suffix Tree)

- **Validated:** 2026-06-24   **Area:** Pattern Matching
- **Canonical method(s):** `SuffixTree.FindAllOccurrences(string / ReadOnlySpan<char>)`, `SuffixTree.Contains(...)`, `SuffixTree.CountOccurrences(...)`; wrappers `GenomicAnalyzer.FindMotif`, `MotifFinder.FindExactMotif`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

Note: the task header loosely described this unit as "naive / KMP / Boyer-Moore". The TestSpec
is authoritative and specifies exact pattern matching via a **suffix tree** (build O(n), match
O(m), report all leaves under the match node). Validated against suffix-tree references.

### Sources opened & what they confirm
- **Wikipedia — "Suffix tree"** (https://en.wikipedia.org/wiki/Suffix_tree, fetched 2026-06-24):
  all z occurrences of a pattern of length m are found in **O(m + z)** by descending to the
  match node and enumerating the leaves below it; overlapping occurrences are naturally
  reported. The BANANA example shows "ANA" at suffix starts 1 and 3 (0-indexed).
- **Rosalind SUBS** (https://rosalind.info/problems/subs/, fetched 2026-06-24): sample
  s = "GATATATGCATATACTT", t = "ATAT" → output **2 4 10**, explicitly **1-indexed**
  ("position = number of symbols to its left, including itself") and explicitly **overlapping**
  (positions 2 and 4 overlap). 0-indexed equivalent: **{1, 3, 9}**.
- **Gusfield (1997), Algorithms on Strings, Trees and Sequences** (reference): the canonical
  suffix-tree exact-matching procedure; "mississippi"/"issi" is the standard worked example.

### Formula / semantics check
- **All occurrences incl. overlapping** reported (leaf enumeration under match node). ✔
- **Index base:** 0-based start positions; spec converts Rosalind 1-based → 0-based explicitly. ✔
- **Empty pattern:** returns [0..n-1]. Not addressed by Wikipedia; spec sources it to formal-
  language theory (ε is a substring of every string at every position) — a standard, documented,
  defensible convention, not implementation-defined. (PASS-with-note level; no source contradicts it.)
- **Pattern not found / pattern longer than text:** empty result. ✔
- **Pattern = whole text:** [0]. ✔
- **Null pattern:** ArgumentNullException (implementation contract). ✔
- **Complexity:** O(m + z) search, O(m) contains, O(m) count via precomputed LeafCount. ✔

### Independent hand cross-checks (0-indexed)
- "GATATATGCATATACTT" / "ATAT" → {1, 3, 9} (= Rosalind {2,4,10}, overlapping). ✔
- "AAAA" / "AA" → {0, 1, 2} (overlapping convention confirmed). ✔
- "banana" / "ana" → {1, 3}; "a" → {1, 3, 5}; "na" → {2, 4}. ✔
- "mississippi" / "issi" → {1, 4}; "i" → {1, 4, 7, 10}. ✔
- "abracadabra" / "a" → {0, 3, 5, 7, 10}. ✔
- empty pattern in "abc" → {0, 1, 2}; pattern-longer-than-text → {}; whole text → {0}. ✔

All match the spec tables and external sources. **Stage A PASS.**

## Stage B — Implementation

### Code path reviewed
- `src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs`
  - `MatchPatternCore` (l.39): walks edges matching P (scalar <8 chars, SIMD `SequenceEqual` ≥8);
    bounds-checked vs `_text.Length` (l.54) so pattern-longer-than-text → `matched=false`.
  - `FindAllOccurrences` (string l.81 / span l.148): empty → `BuildAllStartPositions(n)` = [0..n-1];
    not matched → `Array.Empty<int>()`; else `CollectLeaves` from match node. Null guarded.
  - `CountOccurrences` (l.96/162): empty → n; matched → `node.LeafCount`; else 0.
  - `Contains` (l.23/141): empty → true; else `MatchPatternCore(...).matched`.
  - `CollectLeaves` (l.108): DFS; each leaf → `startPosition = _text.Length + 1 - suffixLength`,
    guarded by `startPosition < _text.Length` to exclude the implicit terminator-only leaf. The
    `+1` correctly accounts for the appended `$`; a suffix of length L (incl. `$`) starts at n+1−L,
    and the terminator-only suffix (L=1) yields n, excluded. Returns ALL leaves = all occurrences
    incl. overlapping.
- Wrappers: `MotifFinder.FindExactMotif` (l.24) and `GenomicAnalyzer.FindMotif` (l.164) both
  uppercase-normalize and delegate to `FindAllOccurrences`; both return **empty for null/empty
  motif** (documented wrapper behavior W3, intentionally distinct from SuffixTree's empty-pattern
  semantics).

### Formula realised correctly?
Yes. 0-based positions; overlapping matches captured because each distinct suffix start under the
match node is a distinct occurrence. The exhaustive test `FindAll_AllSubstrings_MatchLinearSearch`
cross-checks every substring of "mississippi" against `string.IndexOf` with `pos++` (overlapping
reference), proving INV-2/INV-5 including overlaps.

### Cross-verification table recomputed vs code (tests executed)
Ran `SuffixTree.Tests` Search filter: **52 passed / 0 failed**. All spec MUST/SHOULD cases present
and green: banana {1,3}/{1,3,5}/{2,4}, mississippi/issi {1,4}, Rosalind/ATAT {1,3,9},
AAAA/AA {0,1,2}, abracadabra/a {0,3,5,7,10}, empty→[0,1,2], full→[0], not-found / longer-than-text
→ empty, null→throws.

### Variant/delegate consistency
- Span overloads verified identical to string overloads (`FindAll_SpanOverload_MatchesStringOverload`,
  slice/char-array span tests, plus parallel Contains/Count span tests).
- INV-3 (`Count == FindAll.Count`) and INV-4 (`Contains == Count>0`) tested.

### Test quality audit
Assertions check exact sourced position lists (not "no throw" / tautology), are deterministic, and
cover all Stage-A edge cases. Wrapper uppercase/empty behavior covered in Genomics tests.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN** — no defects; no code changes.
- Tests: `SuffixTree.Tests` Search filter = 52 passed / 0 failed.
</content>
</invoke>
