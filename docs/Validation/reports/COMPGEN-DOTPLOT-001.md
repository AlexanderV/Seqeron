# Validation Report: COMPGEN-DOTPLOT-001 — Dot Plot Generation (word-match / k-tuple dot matrix)

- **Validated:** 2026-06-16   **Area:** Comparative
- **Canonical method(s):** `ComparativeGenomics.GenerateDotPlot(string sequence1, string sequence2, int wordSize = 10, int stepSize = 1)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one test-quality defect found and fixed in-session)

## Stage A — Description

### Sources opened this session (independently retrieved, not trusting repo Evidence)

| Source | URL | What it confirmed |
|--------|-----|-------------------|
| Huttley TIB — Dotplot worked example | https://gavinhuttley.github.io/tib/seqcomp/dotplot.html | Match rule `if X[i] == Y[j] then matches[i,j]=0` (k=1). Worked example X=`AGCGT`, Y=`AT`: matches at **(1,1)** and **(5,2)** in the page's 1-based display ⇒ 0-based **(0,0)** and **(4,1)**. Matrix indexed [first-seq-index, second-seq-index] ⇒ confirms x=seq1, y=seq2. |
| EMBOSS dottup manual | https://www.bioinformatics.nl/cgi-bin/emboss/help/dottup | "places where words (tuples) of a specified length have an **exact match** in both sequences." No scoring matrix. Word-size trade-off (longer ⇒ less noise, less sensitive; shorter ⇒ more sensitive, more noise) quoted verbatim. |
| EMBOSS dottup manpage (xenial) | https://manpages.ubuntu.com/manpages/xenial/man1/dottup.1e.html | "Displays a wordmatch dotplot of two sequences." **wordsize default = 10.** |
| Wikipedia — Dot plot (bioinformatics) | https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) | Dot drawn when residues match at the same location. "The main diagonal represents the sequence's alignment with itself; lines off the main diagonal represent similar or repetitive patterns." Noise reduction via tuples (a tuple of 3 = three residues in a row). **Cites Gibbs & McIntyre 1970.** |

### Formula check
The implementation realises the match relation
D = { (i, j) : A[i..i+w−1] = B[j..j+w−1] }, 0 ≤ i ≤ n−w, 0 ≤ j ≤ m−w (case-insensitive),
which is exactly the dottup exact-word-match rule (Source: dottup manual) and, for w=1, the
single-residue rule X[i]==Y[j] (Source: Huttley). Confirmed correct.

### Edge-case semantics
- Sequence shorter than wordSize ⇒ no length-w word ⇒ empty (dottup word formation). ✓
- null/empty ⇒ empty. ✓ (degenerate, consistent with "no words to compare")
- Self-comparison ⇒ full main diagonal (Wikipedia). ✓
- Disjoint alphabets ⇒ empty (dot-on-match rule). ✓
- Default wordSize = 10 (dottup manpage). ✓
- wordSize ≤ 0 / stepSize ≤ 0: not a published dottup case; a non-positive window is undefined, so
  throwing `ArgumentOutOfRangeException` is a reasonable sibling-convention contract (documented as
  an implementation decision, INV-5). PASS-WITH-NOTES-level only; not a biology claim.

### Independent cross-checks (numbers traced to sources / hand computation)
1. **AGCGT vs AT, k=1** → 0-based **{(0,0),(4,1)}**. Source: Huttley page (retrieved this session; 1-based (1,1),(5,2)). Hand re-derivation: A∈AGCGT at {0}, A∈AT at {0} ⇒ (0,0); T at x=4, T∈AT at {1} ⇒ (4,1). No other equal characters. ✓
2. **ACGTACGT self, w=4** → hand-computed from the (sourced) exact-word rule over x=0..4: x=0/4 "ACGT" occurs at {0,4}; x=1 "CGTA"→{1}; x=2 "GTAC"→{2}; x=3 "TACG"→{3} ⇒ **{(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)}**. ✓ (all overlapping occurrences reported, per dottup "all places")
3. **ACGT self, w=1** → distinct residues ⇒ only equal pairs are (i,i) ⇒ exact set **{(0,0),(1,1),(2,2),(3,3)}** (full main diagonal, nothing else). Source: Wikipedia main-diagonal statement + hand computation. ✓

### Findings / divergences (Stage A)
None affecting correctness. Two documented presentation/contract choices (axis orientation x=seq1/y=seq2; case-insensitive folding; non-positive-window throws) are reasonable and labelled as decisions, not as sourced biology. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:1169-1207`
- Eager validation: `wordSize <= 0` / `stepSize <= 0` → `ArgumentOutOfRangeException` (thrown before iterator, so it surfaces immediately, not on first MoveNext).
- Iterator: null/empty or shorter-than-word ⇒ `yield break`. Builds `SuffixTree` over `ToUpperInvariant(seq2)`, slides word `seq1.Substring(i,w).ToUpperInvariant()` for i = 0, step, … ≤ n−w, yields (i, j) for every j from `FindAllOccurrences`.
- `SuffixTree.FindAllOccurrences` (`SuffixTree.Search.cs:81`) collects **all** leaves under the matched node ⇒ all overlapping start positions. Verified it returns every occurrence, not just the first.

### Formula realised correctly?
Yes. Direct enumeration of D; case folding via `ToUpperInvariant`; stepSize sampling on seq1 only;
default wordSize constant = 10 (`DefaultDotPlotWordSize`). Matches the validated Stage-A description.

### Cross-verification table recomputed vs code (via tests, full suite green)
| Case | Sourced/hand value | Code output | Match |
|------|--------------------|-------------|-------|
| AGCGT vs AT, w=1 | {(0,0),(4,1)} | same | ✓ |
| ACGTACGT self, w=4 | {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} | same | ✓ |
| ACGT self, w=1 | {(0,0),(1,1),(2,2),(3,3)} (exact) | same | ✓ |
| AAAA vs CCCC, w=1 | ∅ | ∅ | ✓ |
| ACG vs ACGT, w=4 | ∅ | ∅ | ✓ |
| null/empty | ∅ | ∅ | ✓ |
| w=0,−1 / step=0,−1 | throws AOORE | throws AOORE | ✓ |
| ACGTACGT self, w=4 step=4 | {(0,0),(0,4),(4,0),(4,4)} | same | ✓ |
| acgtacgt vs ACGTACGT, w=4 | M2 set (case-insensitive) | same | ✓ |
| ACGTACGTAC self, default w=10 | {(0,0)} | same | ✓ (locks default=10) |

### Variant/delegate consistency
Single public method; no delegates/`*Fast` variants. Default-parameter path now covered by S3.

### Test quality audit (HARD gate)
- **Defect found & fixed:** M3 (`GenerateDotPlot_SelfComparison_ContainsMainDiagonal`) used
  `Is.SupersetOf` where the **exact** match set is known (ACGT has distinct residues ⇒ the full set
  is precisely the main diagonal). A superset assertion would still pass against an implementation
  that emitted spurious off-diagonal dots — a green-washing weakness. **Tightened to `Is.EquivalentTo`
  {(0,0),(1,1),(2,2),(3,3)}** (sourced value).
- **Coverage added:** S3 (`GenerateDotPlot_DefaultWordSize_IsTenAndMatchesFullWord`) — exercises the
  default-parameter path and locks the documented EMBOSS dottup default wordSize=10 with the
  hand-verified exact set {(0,0)} for a length-10 self-comparison. Previously no test exercised the
  default value.
- Remaining assertions use exact `Is.EquivalentTo` / `Is.Empty` / exact exception type — sourced
  values, not code echoes. C1 is an explicit invariant/property test (INV-4 bound + multiples-of-step
  + (0,0) member); acceptable as a property check, no known single exact set is being weakened.
- Every Stage-A branch covered: match rule (M1/M2/S2/S3), self-diagonal (M3), disjoint (M4),
  word>seq (M5), null/empty (M6), invalid wordSize (M7), invalid stepSize (M8), stepSize sampling
  (S1), case-insensitivity (S2), default wordSize (S3), O(n×m)/multiples property (C1).
- **Honest green:** FULL unfiltered suite `dotnet test` = **6606 passed, 0 failed**, 0 skipped (the
  one Skipped is the pre-existing benchmark guard, not a dot-plot test). `dotnet build` 0 errors;
  the changed test file builds warning-free (the 4 NUnit2007 warnings are pre-existing in an
  unrelated file `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects (Stage B)
One test-quality defect (M3 weak assertion) — **fixed in-session**. One coverage gap (default
wordSize untested) — **closed in-session** (S3). No implementation defect.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches authoritative sources (dottup, Wikipedia, Huttley, Gibbs &
  McIntyre 1970), every non-trivial expected value traces to a source retrieved this session or to a
  hand computation from a sourced rule.
- **Stage B: PASS-WITH-NOTES.** Implementation is correct; one weak test assertion was strengthened
  to the sourced exact set and one default-path test added. No code change needed.
- **End-state: CLEAN.** All defects fixed; build + full suite green (6606/0).
