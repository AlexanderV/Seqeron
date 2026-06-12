# Validation Report: PAT-EXACT-001 — Exact Pattern Search (Suffix Tree)

- **Validated:** 2026-06-12   **Area:** Pattern Matching
- **Canonical method(s):** `SuffixTree.FindAllOccurrences(string / ReadOnlySpan<char>)`, `SuffixTree.Contains(...)`, `SuffixTree.CountOccurrences(...)`; wrappers `GenomicAnalyzer.FindMotif`, `MotifFinder.FindExactMotif`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Suffix tree"** (https://en.wikipedia.org/wiki/Suffix_tree): exact matching of pattern P (length m) finds **all z occurrences** by descending to the match node and **enumerating the leaves of the subtree below it**, in **O(m + z)** time ("Find all z occurrences of the patterns ... in O(m+z) time"). Confirms "ana" in "banana" → starting positions **1 and 3** (0-indexed), i.e. overlapping occurrences are included.
- **Rosalind SUBS** (https://rosalind.info/problems/subs/): for s = "GATATATGCATATACTT", t = "ATAT" the answer is **2, 4, 10 (1-indexed)**; positions overlap (the simple non-overlapping scan would miss 4). 0-indexed: **1, 3, 9**.
- **Gusfield (1997), Algorithms on Strings, Trees and Sequences**: standard suffix-tree exact-matching procedure (build O(n), match O(m), report all leaves = all occurrences incl. overlapping); "mississippi"/"issi" is the canonical worked example.

### Semantics confirmed
- **All occurrences incl. overlapping** are reported (leaf enumeration under match node). ✔
- **Index base:** 0-based start positions. Spec converts Rosalind's 1-based output to 0-based explicitly. ✔
- **Empty pattern:** returns all positions [0..n-1]. Wikipedia does not address this; spec sources it to formal-language theory (ε is a substring of every string at every position) — a standard, defensible convention, **documented**, not implementation-defined. PASS-with-note level only (not a divergence from any source that contradicts it).
- **Pattern not found / pattern longer than text:** empty result. ✔
- **Pattern = whole text:** [0]. ✔
- **Null pattern:** ArgumentNullException (implementation contract). ✔
- **Complexity:** O(m + z) search, O(m) contains, O(m) count via precomputed LeafCount — matches Wikipedia. ✔

### Independent hand-computed cross-checks (0-indexed)
- "AAAA" / "AA" → {0,1,2} (overlapping). ✔
- "banana" / "ana" → {1,3}; "a" → {1,3,5}; "na" → {2,4}. ✔
- "mississippi" / "issi" → {1,4}; "i" → {1,4,7,10}. ✔
- "GATATATGCATATACTT" / "ATAT" → {1,3,9} (= Rosalind {2,4,10}). ✔
- "abracadabra" / "a" → {0,3,5,7,10}. ✔

All match the spec's tables. **Stage A PASS.**

## Stage B — Implementation

### Code path reviewed
- `src/SuffixTree/Algorithms/SuffixTree/SuffixTree.Search.cs`
  - `MatchPatternCore` (l.39): walks edges matching P (scalar < 8 chars, `SequenceEqual`/SIMD ≥ 8); bounds-checked against `_text.Length` so pattern-longer-than-text returns `matched=false`.
  - `FindAllOccurrences` (l.81 string / l.148 span): empty → `BuildAllStartPositions(_text.Length)` = [0..n-1]; not matched → `Array.Empty<int>()`; else `CollectLeaves` from the match node. Null guarded via `EnsureNotNull`.
  - `CountOccurrences` (l.96/162): empty → `_text.Length`; matched → `node.LeafCount`; else 0.
  - `Contains` (l.23/141): empty → true; else `MatchPatternCore(...).matched`.
  - `CollectLeaves` (l.108): DFS over subtree; each leaf yields `startPosition = _text.Length + 1 - suffixLength`, guarded by `startPosition < _text.Length` to exclude the implicit terminator leaf. Returns ALL leaves = all occurrences incl. overlapping.
- Wrappers: `GenomicAnalyzer.FindMotif` (GenomicAnalyzer.cs l.113) and `MotifFinder.FindExactMotif` (MotifFinder.cs l.24) normalize to uppercase and delegate to `FindAllOccurrences`; both return **empty for null/empty motif** (documented wrapper behavior W3, intentionally different from SuffixTree's empty-pattern semantics).

### Formula realised correctly?
Yes. Index base is 0-based; overlapping matches are inherently captured because each distinct suffix start under the match node is a distinct occurrence. The exhaustive test `FindAll_AllSubstrings_MatchLinearSearch` cross-checks every substring of "mississippi" against `string.IndexOf` with `pos++` (overlapping reference), proving INV-2/INV-5 hold including overlaps.

### Cross-verification table recomputed vs code (tests executed)
All spec MUST/SHOULD cases (M1–M17, S1–S7) are present and green: banana {1,3}, mississippi/issi {1,4}, Rosalind/ATAT {1,3,9}, AAAA/AA {0,1,2}, abracadabra/a {0,3,5,7,10}, empty→[0,1,2], full string→[0], not-found/longer-than-text→empty, null→throws.

### Variant/delegate consistency
- Span overloads verified identical to string overloads (`FindAll_SpanOverload_MatchesStringOverload`, slice/char-array span tests).
- INV-3 (`Count == FindAll.Count`) and INV-4 (`Contains == Count>0`) tested.

### Test quality audit
Assertions check exact sourced position lists (not "no throw"), deterministic, and cover all Stage-A edge cases. Wrapper empty/uppercase behavior covered in Genomics.Tests.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN** — no defects; no code changes required.
- Tests: `SuffixTree.Tests` search filter = 70 passed / 0 failed; `Seqeron.Genomics.Tests` motif/pattern filter = 203 passed / 0 failed.
