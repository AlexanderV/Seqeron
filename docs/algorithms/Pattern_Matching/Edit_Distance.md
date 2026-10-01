# Edit Distance (Levenshtein Distance)

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-APPROX-002 |
| Related Projects | N/A |
| Implementation Status | Complete |
| Last Reviewed | 2026-10-01 |

## 1. Overview

Edit distance measures the minimum number of insertions, deletions, and substitutions required to transform one string into another. In this repository the Levenshtein distance and the Sellers end-position search run on the Myers (1999) bit-parallel engine (global form of Hyyrö 2003, multi-word blocks as in edlib), with the Wagner–Fischer column kernel kept as the reference and as the engine of the window search `FindWithEdits` and of the traceback (`GetEditAlignment`, per-hit alignments in edlib CIGAR convention). The Damerau variants — optimal string alignment (restricted) and true Damerau–Levenshtein (Lowrance–Wagner 1975; linear-space engine of Zhao & Sahni 2020) — are separate methods. Weighted (non-unit, additive) insertion/deletion/substitution costs are available through `EditCosts` overloads of `EditDistance`, `GetEditAlignment`, `GetEditAlignmentLinearSpace` and `FindEditEndPositions` (rapidfuzz `Levenshtein` weights semantics).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Levenshtein distance is the standard edit-distance metric for strings when insertions, deletions, and substitutions all cost one unit. In bioinformatics it models sequence differences produced by insertion, deletion, and substitution events and provides a more appropriate notion of distance than Hamming distance when indels are possible. Sources: Levenshtein (1966), Wagner & Fischer (1974), Wikipedia (Levenshtein distance, Edit distance), Berger et al. (2021), Navarro (2001).

### 2.2 Core Model

The recursive definition preserved from the original document is:

$$
lev(a,b) = \begin{cases}
|a| & \text{if } |b| = 0 \\
|b| & \text{if } |a| = 0 \\
lev(tail(a), tail(b)) & \text{if } head(a) = head(b) \\
1 + \min\begin{cases}
lev(tail(a), b) \\
lev(a, tail(b)) \\
lev(tail(a), tail(b))
\end{cases} & \text{otherwise}
\end{cases}
$$

Approximate search (the *k differences* problem, Sellers 1980; Navarro 2001 §5.1) uses the same recurrence with a free start in the text: `C[0, j] = 0`, `C[i, 0] = i`, and a match ends at text position `j` whenever `C[m, j] ≤ k`, where `C[m, j] = min_i ed(P, T[i..j])`. `FindEditEndPositions(...)` returns exactly these `(j, C[m, j])` pairs.

**Bit-parallel engine (Myers 1999; Hyyrö 2003).** Each DP column is encoded by its vertical deltas `Δv ∈ {−1, 0, +1}` as two bit-vectors `Pv`/`Mv` (64 rows per word; ⌈m/64⌉ words). One column step is edlib's `calculateBlock` (Myers' Advance_Block with a horizontal input delta `hin`): `Xv = Eq | Mv; Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq; Ph = Mv | ~(Xh | Pv); Mh = Pv & Xh`, then shift in `hin` and `Pv' = Mh | ~(Xv | Ph); Mv' = Ph & Xv`. The delta entering row 0 is `+1` for global distance (`C[0, j] = j`, Hyyrö's NW form) and `0` for Sellers search (`C[0, j] = 0`); each word's outgoing delta feeds the next word. The score `C[m, j]` is tracked from the horizontal delta at row `m` (bit `(m − 1) mod 64` of the last word). The results are exactly the DP values (tests: exhaustive over `{A,C}^≤6` pairs, random incl. m > 64 and non-ASCII symbols).

**Traceback / CIGAR.** The alignment is recovered from the stored DP columns (query = rows, target = columns). Operations follow edlib (`edlib.h`): `=` match, `X` mismatch, `I` insertion to the query (query character without target counterpart), `D` deletion from the query (target character without query counterpart); EXTENDED CIGAR uses `=`/`X`, STANDARD CIGAR `M`. Tie-break (deterministic): walking back from `(m, n)`, take the diagonal whenever `C[i−1, j−1] + [q_i ≠ t_j] = C[i, j]`, else `I` when `C[i−1, j] + 1 = C[i, j]`, else `D`. edlib's traceback tries `I`, then `D`, then the diagonal; a Python traceback with edlib's order reproduced edlib's CIGAR on 2000/2000 random pairs, so the two differ only on co-optimal ties. Diagonal-first guarantees that when an equal-length window has `ed = Hamming` the returned path is the ungapped Hamming path (on the main diagonal `C[i, i] = HD(prefix_i)` for every `i`, so the diagonal test always succeeds).

**Damerau variants.** Optimal string alignment (OSA, restricted edit distance; Damerau 1964, Boytsov 2011) adds `d[i−2, j−2] + 1` when `a_i = b_{j−1}` and `a_{i−1} = b_j` (no substring edited twice; not a metric: OSA(CA, ABC) = 3 > OSA(CA, AC) + OSA(AC, ABC) = 2). True Damerau–Levenshtein (Lowrance & Wagner 1975) allows insertions/deletions between transposed characters via the last-occurrence table `da`: `d[k−1, l−1] + (i−k−1) + 1 + (j−l−1)`; DL(CA, ABC) = 2. `DL ≤ OSA ≤ Levenshtein`. `FindWithEdits(...)` reports every window `T[i..i+len)` (len ∈ `[max(1, m − k), m + k]`) with `ed(P, window) ≤ k`; the set of its window end positions (with the minimum distance per end) equals the Sellers set.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `d(a, b) = 0` iff the strings are identical under the exact comparison used by `EditDistance(...)` | Zero cost is incurred only when every aligned character matches and lengths agree |
| INV-02 | `d(a, b) >= |len(a) - len(b)|` | At least the length difference must be repaired by insertions or deletions |
| INV-03 | `d(a, b) <= max(len(a), len(b))` | One string can be transformed into the other by deletions plus substitutions |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `s1`, `s2` | `string` | required | Strings compared by `EditDistance(...)` | Null input throws `ArgumentNullException` |
| `[FindWithEdits(string)] sequence` | `string` | required | Sequence searched by `FindWithEdits(...)` | Null or empty input yields no matches |
| `[FindWithEdits(DnaSequence)] sequence` | `DnaSequence` | required | Sequence searched through the typed wrapper | Null input throws because the wrapper dereferences `sequence.Sequence` |
| `pattern` | `string` | required | Pattern compared against variable-length windows | Null or empty input yields no matches |
| `maxEdits` | `int` | required | Maximum allowed edit distance | Negative values throw `ArgumentOutOfRangeException`; values up to `int.MaxValue` are valid (window bound computed in `long`) |
| `[GetEditAlignment] query`, `target` | `string` | required | Strings aligned globally (query = pattern/rows) | Null input throws `ArgumentNullException`; case-sensitive |
| `[GetEditAlignmentLinearSpace] query`, `target` | `string` | required | Same, linear space (Hirschberg) | Null input throws `ArgumentNullException`; case-sensitive |
| `[OptimalStringAlignmentDistance / DamerauLevenshteinDistance] s1`, `s2` | `string` | required | Strings compared | Null input throws `ArgumentNullException`; case-sensitive |
| `costs` (`EditCosts(Insertion, Deletion, Substitution)`) or `insertionCost`, `deletionCost`, `substitutionCost` | `int` × 3 | unit overloads = (1, 1, 1) | Additive costs, rapidfuzz `weights=(insertion, deletion, substitution)`: transforming s1 (query / pattern) into s2 (target / text window), an insertion adds a character of s2 ('D' column), a deletion removes a character of s1 ('I' column) | Each ≥ 0 (else `ArgumentOutOfRangeException`); a result above `int.MaxValue` throws `OverflowException` |
| `[FindEditEndPositions(…, EditCosts)] maxCost` | `int` | required | Maximum weighted cost | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `distance` | `int` | Levenshtein distance between two strings |
| `Position` | `int` | Start of a matching window in `FindWithEdits(...)` |
| `MatchedSequence` | `string` | Window whose edit distance is within threshold |
| `Distance` | `int` | Observed edit distance |
| `MismatchType` | `MismatchType` | `Substitution` when the hit's alignment has no indel (⇔ equal-length window with edit distance = Hamming distance); otherwise `Edit` |
| `MismatchPositions` | `IReadOnlyList<int>` | Pattern-relative indices of the substituted (`X`) pattern characters of the hit's alignment (= the Hamming mismatch indices for `Substitution` hits) |
| `Alignment` | `EditAlignment?` | Hit alignment (pattern = query, window = target); null for Hamming results |
| `EditAlignment` | record | `Distance` (weighted cost for the `EditCosts` overloads), `Operations` (`=`/`X`/`I`/`D` per column), `Cigar` (EXTENDED), `StandardCigar`, `AlignedQuery`, `AlignedTarget` (`-` for gaps), `SubstitutionPositions`, `HasIndels` |

### 3.3 Preconditions and Validation

`EditDistance(...)` throws `ArgumentNullException` when either string is null. `FindWithEdits(string, ...)` returns no matches when the sequence or pattern is null or empty and throws `ArgumentOutOfRangeException` when `maxEdits < 0`. The `DnaSequence` wrapper overload dereferences `sequence.Sequence` directly and therefore throws when `sequence` is null. The approximate-search method uppercases both sequence and pattern before scanning, but the core `EditDistance(...)` routine compares characters as-is.

## 4. Algorithm

### 4.1 High-Level Steps

The DP entry points share one column kernel (`AdvanceColumn`) of the unit-cost Wagner–Fischer DP; they differ only in the top cell of each column. The two distance-only entry points run on the Myers bit-vector engine (`MyersBitVector`), which produces the same column scores.

1. `EditDistance`: Myers/Hyyrö global engine (row-0 delta +1) with the shorter string as bit-vector pattern; return `C[m, n]`. Reference: `EditDistanceDp` (two DP columns, top cell = `j`).
2. `FindEditEndPositions` (Sellers): uppercase inputs; Myers engine with row-0 delta 0 (free text start); yield `(j, C[m, j])` when `C[m, j] ≤ k`. Reference: `FindEditEndPositionsDp` (top cell = 0).
3. `FindWithEdits`: uppercase inputs; for each start `i`, run the start-anchored DP (top cell = window length) over `T[i..i+m+k)`, keeping the columns — one pass gives `ed(P, T[i..i+len))` for every length; for windows with `len ≥ max(1, m − k)` and distance `≤ k` trace the alignment back through the kept columns and yield the hit; stop early once the column minimum exceeds `k` (Ukkonen 1985 cut-off; column minima never decrease).
4. `GetEditAlignment`: full DP matrix (top cell = `j`), then the diagonal-first traceback.
5. `OptimalStringAlignmentDistance`: three rolling rows with the adjacent-transposition case; `DamerauLevenshteinDistance`: Zhao & Sahni (2020) linear-space algorithm as transcribed from rapidfuzz-cpp `damerau_levenshtein_distance_zhao` — rows `R` (current), `R1` (previous), `FR[j]` = `H[k−1, j−2]` saved at the last row k with `a_k = b_j`, `T` = `H[i−2, l−1]` saved at the last column l with `a_i = b_l`, and a last-row-per-symbol table; when `a_i ≠ b_j` only the transpositions with `j − l = 1` (`FR[j] + (i − k)`) or `i − k = 1` (`T + (j − l)`) can beat the Levenshtein moves. Reference: Lowrance–Wagner full matrix with sentinel row/column `m + n` and the `da` table (internal `DamerauLevenshteinDistanceFullMatrix`, test oracle).
6. Weighted costs (`EditCosts`): the weighted Wagner–Fischer recurrence `D[i, j] = min(D[i−1, j] + del, D[i, j−1] + ins, D[i−1, j−1] + [a_i ≠ b_j]·sub)`, `D[i, 0] = i·del`, `D[0, j] = j·ins` (rapidfuzz-cpp `generalized_levenshtein_wagner_fischer`), one 64-bit column kernel (`AdvanceWeightedColumn`) shared by `EditDistance` (two columns over the shorter string; uniform costs c → c × the Myers distance, as rapidfuzz), `GetEditAlignment` (full matrix + the same diagonal → `I` → `D` traceback, so unit costs reproduce the unit path), `GetEditAlignmentLinearSpace` (Hirschberg; the unit overload is this with (1, 1, 1)) and weighted Sellers `FindEditEndPositions` (top cell 0).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `EditDistance` | `O(⌈min(m,n)/64⌉ × max(m,n))` | `O(⌈min(m,n)/64⌉ + σ·⌈min(m,n)/64⌉)` | Myers bit-parallel (σ = distinct pattern symbols) |
| `FindEditEndPositions` | `O(⌈p/64⌉ × s)` | `O(σ·⌈p/64⌉)` | Sellers (1980) semi-global search, Myers engine |
| `FindWithEdits` | `O(s × p × (p + e))` worst case, plus `O(p + len)` traceback per hit | `O(p × (p + e))` | `s` = sequence length, `p` = pattern length, `e` = `maxEdits`; start-anchored DP columns kept for the traceback |
| `GetEditAlignment` | `O(m × n)` | `O(m × n)` | Full Wagner–Fischer matrix + traceback |
| `GetEditAlignmentLinearSpace` | `O(m × n)` (≈ 2× the full DP) | `O(m + n)` | Hirschberg (1975) divide and conquer over the shared column kernel |
| `OptimalStringAlignmentDistance` | `O(m × n)` | `O(n)` | Three rolling rows |
| `DamerauLevenshteinDistance` | `O(m × n)` | `O(n + σ)` | Zhao & Sahni linear space (σ = distinct symbols of s1) |
| `EditDistance(…, EditCosts)` | `O(m × n)` (uniform costs: Myers) | `O(min(m, n))` | Weighted Wagner–Fischer |
| `GetEditAlignment(…, EditCosts)` | `O(m × n)` | `O(m × n)` | Weighted full matrix + traceback |
| `GetEditAlignmentLinearSpace(…, EditCosts)` | `O(m × n)` | `O(m + n)` | Weighted Hirschberg |
| `FindEditEndPositions(…, EditCosts)` | `O(p × s)` | `O(p)` | Weighted Sellers DP |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ApproximateMatcher.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs)

- `ApproximateMatcher.EditDistance(string, string)`: Levenshtein distance (Myers/Hyyrö bit-parallel engine).
- `ApproximateMatcher.GetEditAlignment(string, string)`: Optimal global alignment (`EditAlignment`, edlib CIGAR).
- `ApproximateMatcher.GetEditAlignmentLinearSpace(string, string)`: Optimal global alignment in O(m + n) space (Hirschberg 1975; Myers & Miller 1988).
- `ApproximateMatcher.OptimalStringAlignmentDistance(string, string)`: Restricted Damerau (OSA) distance.
- `ApproximateMatcher.DamerauLevenshteinDistance(string, string)`: True Damerau–Levenshtein distance (Lowrance–Wagner 1975), O(n + σ) space (Zhao & Sahni 2020); internal `DamerauLevenshteinDistanceFullMatrix` = Lowrance–Wagner oracle (tests only).
- `ApproximateMatcher.EditDistance(string, string, int insertionCost, int deletionCost, int substitutionCost)` / `EditDistance(string, string, EditCosts)`: weighted Levenshtein distance (rapidfuzz `Levenshtein.distance(s1, s2, weights=…)`).
- `ApproximateMatcher.GetEditAlignment(string, string, EditCosts)` / `GetEditAlignmentLinearSpace(string, string, EditCosts)`: weighted optimal alignment (full matrix / Hirschberg).
- `ApproximateMatcher.FindEditEndPositions(string, string, int maxCost, EditCosts)`: weighted Sellers end positions.
- internal `EditDistanceDp` / `FindEditEndPositionsDp`: Wagner–Fischer references for the Myers engine (tests only).
- `ApproximateMatcher.FindWithEdits(string, string, int)`: All windows within `maxEdits` (start-anchored DP per start).
- `ApproximateMatcher.FindEditEndPositions(string, string, int)`: Sellers (1980) end positions with minimum distance.
- `ApproximateMatcher.FindWithEdits(DnaSequence, string, int)`: Typed wrapper over the string implementation.

### 5.2 Current Behavior

The core `EditDistance(...)` method is case-sensitive because it compares characters directly. `FindWithEdits(...)` uppercases both the sequence and pattern before scanning and distinguishes substitution-only matches from general edits by comparing the edit distance to the canonical `SequenceExtensions.HammingDistance` on equal-length windows. Every hit carries its alignment (`Alignment`); `MismatchPositions` are the pattern-relative indices of its substituted characters (empty when the only edits are indels). The `DnaSequence` overload is a thin wrapper over the string implementation and does not add its own null guard.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Levenshtein distance with insertion, deletion, and substitution costs of one.
- Space-optimized Wagner-Fischer dynamic programming.
- Approximate search by accepting windows with edit distance at most `maxEdits`.
- Sellers (1980) k-differences search (`FindEditEndPositions`), cross-checked against edlib infix (HW) mode and Navarro's `survey`/`surgery` example.

- Myers (1999) bit-parallel edit distance in the global form of Hyyrö (2003), multi-word blocks, transcribed from edlib's `calculateBlock`; identical to the DP reference (exhaustive + random tests) and to rapidfuzz/edlib (3160 random pairs, lengths 0–300, incl. non-ASCII; 600 Sellers cases vs Python DP and edlib HW best-distance/end locations).
- Traceback with edlib's operation/CIGAR convention; deterministic diagonal-first tie-break (documented above). Cross-check vs edlib `align(q, t, mode='NW', task='path')`: 2000 random pairs — distance and path validity 2000/2000, CIGAR identical in 116, and identical in 115/115 pairs whose optimal path is unique; FindWithEdits: 13019 hits in 500 random cases — windows = brute force, CIGAR replays with cost = distance, 7461 identical to edlib NW on (pattern, window).
- Optimal string alignment and true Damerau–Levenshtein distances; rapidfuzz `OSA`/`DamerauLevenshtein` and jellyfish agree on 4507 pairs (classic CA/ABC: OSA 3, DL 2). Note: `pyxDamerauLevenshtein` returns 3 for CA/ABC — it implements OSA, not unrestricted DL.

- Linear-space traceback (B05 follow-up, 2026-09-30): `GetEditAlignmentLinearSpace` — Hirschberg (1975) [10] for the unit-cost edit distance (Myers & Miller 1988 [11]). Split the query at h = ⌊m/2⌋; forward scores F[j] = ed(q[0..h), t[0..j)) and reverse scores R[j] = ed(q[h..m), t[j..n)) from the shared weighted column kernel `AdvanceWeightedColumn` with costs (1, 1, 1) (no second DP; values identical to `AdvanceColumn`); split the target at the smallest j minimising F[j] + R[j]; recurse. Base cases: empty side → all `D` / all `I`; one query character → `=` at its last occurrence in the target, else `X` against the last target character, all other target characters `D`. **Why a separate method:** the diagonal-first traceback of `GetEditAlignment` is defined on the full forward matrix, which a divide-and-conquer split never materialises, so the same co-optimal path cannot be guaranteed in linear space — measured: identical paths in 1757/3501 random pairs. Guarantees instead: optimal distance and a valid script. Cross-check: 3501 pairs (3000 random, lengths 0–120 over {AC, ACGT, 10 letters, non-ASCII}, half of them mutated copies; 500 exhaustive-small ≤ 4×4 over {A,C}; one 3000 × 3300 pair): distance = edlib 1.3 NW `editDistance` = `GetEditAlignment` = `EditDistance` 3501/3501, CIGAR replays with cost = distance 3501/3501.

- Linear-space unrestricted Damerau–Levenshtein (B05 audit group B, 2026-10-01): Zhao, C.; Sahni, S. (2020) [12] (paper hosts blocked; algorithm read from rapidfuzz-cpp `rapidfuzz/distance/DamerauLevenshtein_impl.hpp`, `damerau_levenshtein_distance_zhao`, which cites it). Cross-check: rapidfuzz 3.14.6 `DamerauLevenshtein.distance` 7620/7620 pairs (4000 random lengths 0–120 over {AC, ACGT, 10 letters, non-ASCII Greek/Cyrillic/€, mixed with 中/ü}, half mutated copies with transpositions; 3600 small pairs over {A,B,C} ≤ 4; 20 pairs up to 2500 × 2500) and jellyfish `damerau_levenshtein_distance` 7600/7600; = the Lowrance–Wagner full matrix exhaustively for all 364² strings ≤ 5 over {A,B,C}.
- Weighted additive costs (B05 audit group B, 2026-10-01): semantics of rapidfuzz `Levenshtein.distance(s1, s2, weights=(insertion, deletion, substitution))` (rapidfuzz-cpp `Levenshtein_impl.hpp`: `generalized_levenshtein_wagner_fischer`, `levenshtein_distance` uniform/InDel shortcuts; cache[i] = i·delete_cost, so deletions remove s1 characters and insertions add s2 characters). Cross-check: `EditDistance` (both overloads) = rapidfuzz on 4010 pairs (weights 0–6, uniform, up to 1000; lengths 0–80 plus 10 pairs ≤ 1500); `GetEditAlignment` and `GetEditAlignmentLinearSpace` on 3000 pairs — distance = rapidfuzz and the operations replay both strings with weighted cost = distance 3000/3000 (rapidfuzz `editops` is unit-cost only, so weighted path validity is checked by replay); weighted Sellers on 800 cases (8330 hits) = Python brute force `min_i Levenshtein.distance(p, t[i..j], weights)`. Unit costs return exactly the unit methods' results and paths (tests).

**Not implemented:**

- Affine gap costs and substitution matrices — use the pairwise aligners in `SequenceAligner`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `""` vs `"abc"` | Distance `3` | Three insertions or deletions are required |
| Empty string sequence or pattern in `FindWithEdits(...)` | Returns no matches | Explicit source guard |
| Null `DnaSequence` input in `FindWithEdits(...)` | Throws | The typed wrapper dereferences `sequence.Sequence` |
| `maxEdits < 0` | Throws `ArgumentOutOfRangeException` | Invalid threshold |
| Equal-length strings | Levenshtein distance is at most the Hamming distance | Substitutions alone can realize the Hamming path |

### 6.2 Limitations

`GetEditAlignment` keeps an `O(m·n)` matrix (for long inputs use `GetEditAlignmentLinearSpace`, O(m + n) space, same distance, possibly a different co-optimal path); `DamerauLevenshteinDistance` runs in O(n + σ) space (Zhao & Sahni) and returns the distance only (no transposition-aware traceback); `FindWithEdits` keeps `O(p·(p+e))` DP columns per start. Among co-optimal alignments exactly one (diagonal-first) is returned; it may differ from edlib's (I-first) path. The core distance methods are also case-sensitive, so callers who need normalized comparisons must uppercase or otherwise normalize inputs before calling it directly.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

`kitten -> sitting` has edit distance `3`:

1. `kitten -> sitten`
2. `sitten -> sittin`
3. `sittin -> sitting`

## 8. References

1. Levenshtein, V.I. (1966). "Binary codes capable of correcting deletions, insertions, and reversals." Soviet Physics Doklady, 10(8): 707–710.
2. Wagner, R.A.; Fischer, M.J. (1974). "The String-to-String Correction Problem." Journal of the ACM, 21(1): 168–173.
3. Sellers, P.H. (1980). "The theory and computation of evolutionary distances: Pattern recognition." Journal of Algorithms, 1(4): 359–373.
4. Ukkonen, E. (1985). "Finding approximate patterns in strings." Journal of Algorithms, 6(1): 132–137.
4a. Myers, G. (1999). "A fast bit-vector algorithm for approximate string matching based on dynamic programming." Journal of the ACM, 46(3): 395–415.
4b. Hyyrö, H. (2003). "A bit-vector algorithm for computing Levenshtein and Damerau edit distances." Nordic Journal of Computing, 10(1): 29–39.
4c. Šošić, M.; Šikić, M. (2017). "Edlib: a C/C++ library for fast, exact sequence alignment using edit distance." Bioinformatics, 33(9): 1394–1395; source opened: raw.githubusercontent.com/Martinsos/edlib/master/edlib/src/edlib.cpp (`calculateBlock`, `obtainAlignmentTraceback`, `edlibAlignmentToCigar`) and `edlib/include/edlib.h` (EDLIB_EDOP_*, EDLIB_CIGAR_*).
4d. Damerau, F.J. (1964). "A technique for computer detection and correction of spelling errors." Communications of the ACM, 7(3): 171–176.
4e. Lowrance, R.; Wagner, R.A. (1975). "An extension of the string-to-string correction problem." Journal of the ACM, 22(2): 177–183.
4f. Boytsov, L. (2011). "Indexing methods for approximate dictionary searching: comparative analyses." ACM Journal of Experimental Algorithmics, 16: 1.1.
5. Navarro, G. (2001). "A guided tour to approximate string matching." ACM Computing Surveys, 33(1): 31–88.
6. Berger, B.; Waterman, M.S.; Yu, Y.W. (2021). "Levenshtein Distance, Sequence Comparison and Biological Database Search." IEEE Transactions on Information Theory, 67(6): 3287–3294.
7. Rosetta Code - Levenshtein Distance: https://rosettacode.org/wiki/Levenshtein_distance
8. Wikipedia - Levenshtein Distance: https://en.wikipedia.org/wiki/Levenshtein_distance
9. Wikipedia - Edit Distance: https://en.wikipedia.org/wiki/Edit_distance
10. Hirschberg, D.S. (1975). "A linear space algorithm for computing maximal common subsequences." Communications of the ACM, 18(6): 341–343 (bibliographic record via WebSearch/Semantic Scholar, 2026-09-30).
11. Myers, E.W.; Miller, W. (1988). "Optimal alignments in linear space." CABIOS, 4(1): 11–17 (bibliographic record via WebSearch, 2026-09-30).
12. Zhao, C.; Sahni, S. (2020). "Linear space string correction algorithm using the Damerau-Levenshtein distance." BMC Bioinformatics, 21(Suppl 1), doi:10.1186/s12859-019-3184-8; predecessor Zhao & Sahni (2019) "String correction using the Damerau-Levenshtein distance", BMC Bioinformatics 20(Suppl 11):277 (bibliographic records via WebSearch 2026-10-01; biomedcentral/springer/scispace/PMC blocked). Implementation transcribed from rapidfuzz-cpp `rapidfuzz/distance/DamerauLevenshtein_impl.hpp` (raw.githubusercontent.com, main, 2026-10-01).
13. rapidfuzz-cpp `rapidfuzz/distance/Levenshtein_impl.hpp` (weighted `generalized_levenshtein_wagner_fischer`, `LevenshteinWeightTable`), raw.githubusercontent.com main, 2026-10-01; rapidfuzz 3.14.6 (PyPI).
