# Dot Plot Generation

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics / Sequence Comparison |
| Test Unit ID | COMPGEN-DOTPLOT-001 |
| Related Projects | Seqeron.Genomics.Analysis, SuffixTree |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

A dot plot is a graphical method for comparing two biological sequences: one sequence is laid along the x-axis and the other along the y-axis, and a dot is drawn at every position pair where the sequences match. Regions of similarity appear as diagonal runs of dots, while insertions/deletions break the diagonal and repeats add extra diagonals [3]. This implementation produces a *word-match (k-tuple) dot plot* in the style of EMBOSS `dottup`: a dot at (x, y) marks an **exact** match of a length-`wordSize` word starting at `sequence1[x]` and `sequence2[y]` [2]. It is an exact, deterministic, specification-driven enumeration (no scoring matrix, no approximate matching).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The dot-matrix ("diagram") method was introduced by Gibbs & McIntyre (1970) for comparing amino-acid and nucleotide sequences [1]. It remains a standard first-look tool for detecting similarity, repeats, and rearrangements between two sequences without committing to a single global alignment [3].

### 2.2 Core Model

Let A = sequence1 (length n) and B = sequence2 (length m), and let w = `wordSize`. Define the match relation

> D = { (i, j) : A[i..i+w−1] = B[j..j+w−1] }, for 0 ≤ i ≤ n−w and 0 ≤ j ≤ m−w,

where equality is character-by-character (case-insensitive). The dot plot is the set D. For w = 1 this reduces to the classic single-residue rule "place a dot at (i, j) iff X[i] == Y[j]" [5]. Using w > 1 (a *tuple*) suppresses random single-residue matches because the chance of w consecutive matches is much lower [3][2]. The implementation samples i in steps of `stepSize` along A.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | (x, y) ∈ output ⇔ A[x..x+w−1] = B[y..y+w−1] (case-insensitive) | Direct enumeration of the match relation D [2][5] |
| INV-02 | Self-comparison (A = B) contains every (i, i), 0 ≤ i ≤ n−w (full main diagonal) | A word always matches itself; "the main diagonal represents the sequence's alignment with itself" [3] |
| INV-03 | Empty output when either sequence is null/empty or shorter than w | No length-w word can be formed [2][3-manpage] |
| INV-04 | All overlapping occurrences are reported (every matching y per word), each x is a multiple of `stepSize` | `SuffixTree.FindAllOccurrences` returns all start positions; i advances by `stepSize` |

### 2.5 Comparison with Related Methods

| Aspect | dottup-style (this algorithm) | dotmatcher |
|--------|-------------------------------|------------|
| Match criterion | exact word (tuple) match | windowed score vs substitution matrix ≥ threshold [2] |
| Speed / sensitivity | fast, less sensitive [2] | slower, tunable sensitivity |
| Parameters | wordSize, stepSize | window size, threshold, scoring matrix |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence1 | string | required | First sequence; positions → x (0-based) | null/empty ⇒ empty result |
| sequence2 | string | required | Second sequence; positions → y (0-based) | null/empty ⇒ empty result |
| wordSize | int | 10 | Exact-match word (tuple) length [3-manpage] | must be > 0 |
| stepSize | int | 1 | Sampling step along sequence1 | must be > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | IEnumerable<(int x, int y)> | Lazily yielded match coordinates; x is a 0-based start position in sequence1, y a 0-based start position in sequence2 |

### 3.3 Preconditions and Validation

0-based indexing on both axes. Null/empty either input, or a sequence shorter than `wordSize`, yields an empty enumeration (no exception). Comparison is case-insensitive (both sequences are upper-cased before matching; T/U are not normalized). `wordSize ≤ 0` or `stepSize ≤ 0` throws `ArgumentOutOfRangeException` (a non-positive window is undefined for a word-match dot plot).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `wordSize > 0` and `stepSize > 0`.
2. Return empty if either sequence is null/empty or shorter than `wordSize`.
3. Build a suffix tree over upper-cased sequence2.
4. For each start i = 0, stepSize, 2·stepSize, … ≤ n−w, extract the word A[i..i+w−1] and query the suffix tree for all occurrences j in B; yield (i, j) for each.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Default wordSize = 10** mirrors EMBOSS `dottup` [3-manpage]; **stepSize = 1** is the exhaustive sampling of every start position [1].
- **Data structure — repository SuffixTree:** sequence2 is indexed once (O(m)); each word is then located in O(w + occurrences). This is the suffix tree's strong case (many exact-match queries against one fixed text) and is used here in place of a naive O(n·m) scan; the returned match set is identical to the naive definition.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| GenerateDotPlot | O(m + (n/stepSize)·w + K) | O(m) | m = build cost; K = total dots reported; worst case K = O(n·m) (e.g., self-comparison of a homopolymer), so the bound is O(n·m) overall |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.GenerateDotPlot(sequence1, sequence2, wordSize = 10, stepSize = 1)`: enumerates exact word-match coordinates as `(x, y)` pairs.

### 5.2 Current Behavior

The public method validates parameters eagerly, then delegates to a private iterator (`GenerateDotPlotIterator`) so that `ArgumentOutOfRangeException` is thrown immediately rather than being deferred until first enumeration. Matching is case-insensitive via `ToUpperInvariant`. **Search-reuse decision:** the repository `SuffixTree` is used for word location — exact-match occurrence enumeration with many queries against one text fits the suffix tree's O(m) build / O(w+k) query profile exactly; a naive scan was rejected as algorithmically inferior while producing the same output.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact word-match dot placement: (x, y) reported iff the length-`wordSize` words match [2][5].
- Word-size noise/sensitivity trade-off via the `wordSize` parameter [2][3].
- Default `wordSize` = 10 [3-manpage].
- Full main diagonal under self-comparison [3] (INV-02).

**Intentionally simplified:**

- Case folding: both inputs are upper-cased so matching is case-insensitive; consequence: lower-case soft-masked residues still participate in matches (the EMBOSS tools likewise compare uppercased sequence).

**Not implemented:**

- Substitution-matrix scored windows with a threshold (the EMBOSS `dotmatcher` model); users needing scored similarity should rely on a scoring-based aligner — no current alternative in this class.
- Rendering/plotting to an image; this method returns coordinates only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | x = sequence1, y = sequence2 | Assumption | Axis orientation is a presentation choice; transposing inputs transposes the plot | accepted | ASM-style note in Evidence |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty result | no words to compare (INV-03) |
| sequence shorter than wordSize | empty result | no length-w word can start [2] |
| identical sequences | full main diagonal + any internal-repeat diagonals | INV-02 [3] |
| disjoint alphabets | empty result | dot only on exact match [5] |
| wordSize ≤ 0 or stepSize ≤ 0 | ArgumentOutOfRangeException | window undefined |

### 6.2 Limitations

Exact matching only — no mismatches, gaps, or scoring. Small `wordSize` produces off-diagonal chance dots (documented noise, not a defect [2]). Output size is O(n·m) in the worst case (highly repetitive input). T/U are not normalized, so an RNA word will not match its DNA counterpart.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// k=1 dot plot of AGCGT (x) vs AT (y): dots where characters match.
var dots = ComparativeGenomics.GenerateDotPlot("AGCGT", "AT", wordSize: 1).ToList();
// dots == { (0,0)  // A==A , (4,1)  // T==T }
```

**Numerical walk-through:** For `AGCGT` vs `AT`, wordSize 1, the only equal characters are the A's at (x=0, y=0) and the T's at (x=4, y=1) [5].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_GenerateDotPlot_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_GenerateDotPlot_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [COMPGEN-DOTPLOT-001-Evidence.md](../../../docs/Evidence/COMPGEN-DOTPLOT-001-Evidence.md)
- Related algorithms: [Synteny_Block_Detection](./Synteny_Block_Detection.md)

## 8. References

1. Gibbs AJ, McIntyre GA. 1970. The Diagram, a Method for Comparing Sequences. Its Use with Amino Acid and Nucleotide Sequences. *Eur. J. Biochem.* 16(1):1–11. https://doi.org/10.1111/j.1432-1033.1970.tb01046.x
2. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite — `dottup` word-match dot plot. *Trends Genet.* 16(6):276–277. Manual: https://www.bioinformatics.nl/cgi-bin/emboss/help/dottup ; manpage (default wordsize 10): https://manpages.ubuntu.com/manpages/xenial/man1/dottup.1e.html
3. Wikipedia contributors. Dot plot (bioinformatics). https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) (accessed 2026-06-14).
4. Huttley G. Topics in Bioinformatics — Dotplot. https://gavinhuttley.github.io/tib/seqcomp/dotplot.html (accessed 2026-06-14).
