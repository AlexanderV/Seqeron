# Direct Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-DIRECT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (exact maximal pairs; reverse-complement pairs; k-mismatch repeats; supermaximal repeats) |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Direct repeat detection identifies nucleotide sequences that recur in the same 5'→3' orientation at multiple genomic positions [1][2]. Unlike inverted repeats, the downstream copy preserves the original sequence rather than its reverse complement. The repository implements exact direct-repeat discovery in `RepeatFinder.FindDirectRepeats`, reporting **maximal repeated pairs** — the convention of the reference tools MUMmer `repeat-match` [5] and REPuter [6] and of Gusfield's maximal-pair definition [7]. Spacing between copies is configurable, so adjacent tandem-like repeats, separated direct repeats and (with negative spacing) overlapping copies can all be reported. Three companion enumerations share the same maximal-pair engine (§4.4): reverse-complement maximal pairs (`FindReverseComplementRepeats`, repeat-match without `-f` / Vmatch `-p`), maximal k-mismatch degenerate repeats (`FindApproximateDirectRepeats`, REPuter / Vmatch `-h`) and supermaximal repeats (`FindSupermaximalRepeats`, Gusfield §7.12.1 / Vmatch `-supermax`).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A direct repeat consists of two or more identical sequence copies oriented in the same 5'→3' direction [1]. A canonical representation is:

```text
5' TTACG------TTACG 3'
3' AATGC------AATGC 5'
```

where `------` is an intervening spacer that may be zero bases long [1]. Direct repeats are biologically associated with transposable-element boundaries, homologous recombination, replication-slippage-mediated deletions, and some regulatory architectures [1][2][3]. The legacy reference set also notes that tandem trinucleotide expansions are a special case of direct-repeat biology and underlie disorders such as Huntington's disease, Fragile X syndrome, spinocerebellar ataxias, Friedreich's ataxia, and myotonic dystrophy [4].

### 2.2 Core Model

For a sequence $S$ of length $n$ (0-based) and positions $i < j$, let $L(i,j)$ be the length of the longest common
prefix of $S[i..)$ and $S[j..)$, where only A/C/G/T match (any other symbol never matches — MUMmer `mummer -n` [5]).
The pair $(i, j, L)$ with $L = L(i,j) \ge 1$ is a **maximal repeated pair** [7] when it is

- **right-maximal** — $L$ is the full common-prefix length (the copies cannot be extended to the right), and
- **left-maximal** — $i = 0$ or $S[i-1] \ne S[j-1]$ (or the left neighbour is not A/C/G/T).

This is exactly the test in MUMmer `repeat-match` (`Data[i-1] == Data[j-1] || Data[i+n] == Data[j+n]` → skip) [5].
The method reports every maximal pair with

$$
\text{minLength} \le L \le \text{maxLength}, \qquad \text{Spacing} = j - i - L \ge \text{minSpacing}.
$$

A maximal repeat longer than `maxLength` is **not** reported (it is not truncated into sub-windows).
Each exact direct repeat is thus reported once, at its full extent, instead of once per nested window.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `RepeatSequence` is identical at `FirstPosition` and `SecondPosition`. | $L$ is a common-prefix length. |
| INV-02 | `Spacing = SecondPosition - FirstPosition - Length`. | Constructed from the coordinates. |
| INV-03 | When `minSpacing ≥ 0`, reported copies do not overlap. | Spacing filter. |
| INV-04 | Each `(FirstPosition, SecondPosition)` pair occurs at most once; `FirstPosition < SecondPosition`. | $L(i,j)$ is unique per pair. |
| INV-05 | Every result is left- and right-maximal. | Definition [5][7]. |
| INV-06 | Widening `[minLength, maxLength]` or lowering `minSpacing` yields a superset. | Pure filters on a fixed pair set. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search (case-insensitive). | `DnaSequence` overload throws on `null`; `string` overload returns no results for `null` or empty. |
| `minLength` | `int` | `5` | Minimum maximal-repeat length. | ≥ 2 (both overloads, eager `ArgumentOutOfRangeException`). |
| `maxLength` | `int` | `50` | Maximum maximal-repeat length; longer repeats are not reported. | ≥ `minLength`. |
| `minSpacing` | `int` | `1` | Minimum `Spacing = j − i − L`. | Any value: `0` admits abutting copies, negative admits overlapping copies, `int.MinValue` = all maximal pairs (repeat-match `-f`). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `FirstPosition` | `int` | 0-based position of the first repeat copy. |
| `SecondPosition` | `int` | 0-based position of the downstream repeat copy. Results are sorted by (`FirstPosition`, `SecondPosition`). |
| `RepeatSequence` | `string` | Exact repeated sequence shared by both copies. |
| `Length` | `int` | Length of the repeat sequence in bases. |
| `Spacing` | `int` | Number of nucleotides between the two copies. |

### 3.3 Preconditions and Validation

Both overloads throw `ArgumentOutOfRangeException` eagerly (at the call) when `minLength < 2` or `maxLength < minLength`; the `DnaSequence` overload throws `ArgumentNullException` for `null`. Positions are 0-based.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the input (string overload); encode A/C/G/T as 0–3 and every other symbol at position p as the unique code 4 + p (never matches).
2. Build the suffix array (prefix doubling) and the Kasai LCP array — the shared helpers `SequenceComplexity.BuildSuffixArray` / `BuildLcpArray` (also used by the LZ76 factorization).
3. Traverse the lcp-interval tree bottom-up with a stack (Abouelhoda et al. 2004 [8]; Gusfield 1997 §7.12.3 [7]). Each interval keeps its suffix start positions in one linked list per left-character class (A, C, G, T, or "unique" for position 0 / non-ACGT left neighbour).
4. When a child interval is merged into its parent of string depth ℓ and `minLength ≤ ℓ ≤ maxLength`, emit every pair (p from the child, q from the parent's lists) whose left classes differ (or are "unique"): these pairs have LCP exactly ℓ (right-maximal) and are left-maximal. Keep those with `Spacing ≥ minSpacing`.
5. Concatenate the child's lists into the parent; finally sort results by (`FirstPosition`, `SecondPosition`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Suffix array + LCP | `O(n log² n)` | `O(n)` | Prefix doubling with comparison sort; Kasai LCP `O(n)`. |
| Maximal-pair enumeration | `O(n + z)` | `O(n)` | `z` = maximal pairs with length in `[minLength, maxLength]` (pairs failing the spacing filter are still enumerated). |

Measured: 1 Mb random DNA with planted repeats, `minLength` 16 → 1.1 s (MUMmer repeat-match 0.9 s), identical output.

### 4.4 Enumeration variants (B04 audit WP2)

All three reuse the same suffix-array/LCP maximal-pair engine (`EnumerateMaximalPairs`, shared helpers `SequenceComplexity.BuildSuffixArray` / `BuildLcpArray`); A/C/G/T only, case-insensitive, eager validation, `null`/empty → empty.

**Reverse-complement maximal pairs** — `FindReverseComplementRepeats(seq, minLength = 5, maxLength = 50, minSpacing = 1)`. Definition (Vmatch manual, Appendix A [10]): a palindromic match of `S` with itself, `S[i..i+L) = revcomp(S[k..k+L))`, `i ≤ k`, maximal = extendable neither outward (`S[i−1]` vs `S[k+L]`) nor inward (`S[i+L]` vs `S[k−1]`). `i = k` is allowed (reverse palindromes such as `GAATTC`). Engine run on `S · # · revcomp(S)` with per-(strand, left character) lists, only cross-strand pairs emitted, each pair kept in the orientation `i ≤ k` (as `repeat-match.cc` `List_Matches` does). Coordinates: `FirstPosition = i`, `SecondPosition = k`, both 0-based forward-strand starts (Vmatch `-p` prints these); MUMmer `repeat-match` prints `i + 1` and `k + L` (1-based **last** base of copy 2, where it starts on the reverse strand) followed by `r` — convert with `k = Start2 − L`. `Spacing = k − i − L` (the hairpin loop; negative for overlapping copies). One `(i, k)` may carry several lengths on different anti-diagonals (T₇ facing A₇), so results are ordered by (i, k, L). **repeat-match quirk** (documented, not reproduced): it terminates `S` and `revcomp(S)` with the same `$`, so a suffix of `revcomp(S)` equal to a suffix of `S` shares one leaf and that leaf's pairs are lost (e.g. `ATAAT` → (0,0,2) `AT` missing; 31/3 000 random cases); it aborts on a sequence equal to its own reverse complement. Appending one `N` to the FASTA removes both effects (0 mismatches).

**Maximal k-mismatch (degenerate) repeats** — `FindApproximateDirectRepeats(seq, minLength = 10, maxMismatches = 1, maxLength = ∞, minSpacing = 1, excludeContained = false)`. Definition (REPuter [6][9]; Vmatch Appendix A [10]): equal-length copies `S[i..i+L)`, `S[j..j+L)`, `i < j`, Hamming distance ≤ k (non-ACGT = mismatch, the Vmatch wildcard rule); maximal = cannot be extended left or right on its diagonal without exceeding k (a window may end on a mismatch). `excludeContained = true` applies Vmatch's literal Appendix A reading — "not contained in another k-mismatch match" across diagonals — which is exactly what `vmatch -h k -allmax` prints (and what Vmatch does **not** apply to exact repeats: `vmatch -l` keeps (0,2,6) inside (0,1,7) in A₈; hence the per-diagonal default, which for k = 0 equals `FindDirectRepeats`). Algorithm = REPuter seed-and-extend made complete: every such repeat of length ≥ m contains an exact maximal pair of length ≥ ⌊m/(k+1)⌋ (pigeonhole; Vmatch `-seedlength` rule); for each seed the first k + 1 mismatches left and right give all maximal windows containing it (splits a + b = k, or both sides at the boundary); duplicates from several seeds removed. Cross-diagonal containment: windows on one diagonal are never nested, so a containing window on diagonal d′ (|d − d′| ≤ Lmax − L) is found by binary search. Requires `maxMismatches < minLength`. Cost O(n log² n + s·k + z) for s seeds; small ⌊m/(k+1)⌋ makes s grow like n²·4^−⌊m/(k+1)⌋ (1 Mb, m = 30, k = 2: 0.8 s).

**Supermaximal repeats** — `FindSupermaximalRepeats(seq, minLength = 5)`. Definition (Gusfield 1997 §7.12.1 [7]; Vmatch `-supermax`): a maximal repeat that never occurs as a substring of any other maximal repeat. Gusfield Theorem 7.12.4 on the suffix array: an lcp-interval whose children are all singletons (a local maximum of the LCP array) whose suffixes have pairwise distinct left characters (position 0 / non-ACGT left neighbour = distinct). One record per repeated string with every occurrence (ascending; Vmatch prints each position pair instead). O(n log² n).

### 4.5 Degenerate repeats: k-differences and palindromic (B04 audit WP8)

`FindDegenerateRepeats(seq, minLength = 10, maxDifferences = 1, distance = Edit, reverseComplement = false, maxLength = ∞, minSpacing = 1)` → `DegenerateRepeatResult(FirstPosition, FirstLength, SecondPosition, SecondLength, Distance, Spacing, FirstCopy, SecondCopy, IsReverseComplement)` covers the four Vmatch degenerate modes `vmatch [-p] -l m (-h|-e) k -allmax` (REPuter [9]).

- **Definitions** (Vmatch manual App. A [10]): a match `(l, i, r, j)` pairs `u = S[i..i+l)` with `w = S[j..j+r)`; direct `u ≈ w` (`i < j`), palindromic `u ≈ revcomp(w)` (`i ≤ j`). Hamming: `l = r`, `d_H ≤ k` ("k-mismatch match"); edit: unit-cost `d_E ≤ k` (mismatch, insertion, deletion; "k-differences match"). Contained: `i′ ≤ i ≤ i+l ≤ i′+l′` and `j′ ≤ j ≤ j+r ≤ j′+r′`; maximal = not contained in another match of the same kind; both `l, r ≥ minLength` (`-l`). Wildcards always mismatch.
- **Vmatch conventions** (source `kurtz/extendED.c`, `mcontain.c`, Vmatch 2.3.1 [11]): direct edit matches pass `acceptmatch` — the right instance must not be embedded in the left one and overlapping instances need a non-overlapping part `(j − i) + (j + r) − (i + l)` larger than the distance (removes trivial self-alignments like `S[0..n)` vs `S[1..n)`); palindromic matches are compared in both orientations (the mirror `(r, j, l, i)` is the same pair of strings) and reported with `i ≤ j` (for `i = j`, `l ≠ r` both orientations).
- **Algorithm**: pigeonhole seeds (exact maximal pairs ≥ ⌊m/(k+1)⌋; direct engine of §4.1, palindromic engine of §4.4 in both orientations). Hamming: the (k+1)-th-mismatch windows of §4.4 on `S` vs `S` or `S` vs `revcomp(S)`, then cross-diagonal containment. Edit: greedy furthest-reaching fronts (Ukkonen 1985 / Myers 1986 [12], as Vmatch `frontSEP.c` [11]) left and right of each seed, every combination of front points with a + b ≤ k is a candidate with distance min(a + b); per-seed then global containment (max segment tree over the u-intervals, output-sensitive). For a maximal match the rerouting argument (an optimal alignment through any error-free block can follow the seed diagonal) makes min(a + b) = `d_E`.
- **Vmatch shortcut not reproduced by default** (opt-in since WP15, §4.6): Vmatch stops a left extension that crosses another exact match ≥ the seed length (`evalentrybackward`, "seed … detected while scanning"), assuming the match is found from that seed; for edit matches the other seed's alignment differs, so stock Vmatch misses some maximal matches (and then prints contained ones). The same Vmatch release built from source with that one test disabled (`VM_NOPRUNE`) and a brute force of the definition agree with this method on every case (Evidence §WP8).
- Cost: O(n log² n + s·k⁴ + c log c) for s seeds and c candidates; 1 Mb random + planted repeats (m 30, k 2): direct 0.8 s, palindromic 1.8 s (Hamming 0.9 / 1.7 s).

### 4.6 Vmatch output modes and stock compatibility (B04 audit WP15)

`FindDegenerateRepeats(seq, minLength, maxDifferences, distance, reverseComplement, maxLength, minSpacing, DegenerateRepeatReporting reporting, bool vmatchCompatible = false)` (and `FindApproximateDirectRepeats(…, excludeContained, reporting, vmatchCompatible)` for Hamming/direct). The MCP tools `find_degenerate_repeats` and `find_approximate_direct_repeats` expose these as `reporting` = `allMaximal` | `bestPerSeed` and `vmatchCompatible`.

- **`BestPerSeed`** is Vmatch's default output when `-allmax`, `-best` and `-complete` are all absent. The manual (`virtman.tex`) says: "for each seed a best match, i.e. one with a minimum E-value is output … no limit on the number of matches".
  - **Extension.** Each exact maximal seed is extended exactly as Vmatch does it. Hamming follows `hammingextend` / `extendmismatchesleft/right` (`kurtz/extendHD.c`). Edit follows `editextend` (`kurtz/extendED.c`) over the greedy fronts of `frontSEP.c`/`front.gen`, including their conventions (minus-infinity = −max(ulen, vlen), the band once p > min(ulen, vlen), the same-position rule).
  - **Candidate order.** Candidates are visited in Vmatch's order: distance, left errors, left diagonal, right diagonal.
  - **Choice of the best candidate.** It follows `cmpmatches` (`include/extcmp.c`): smaller E-value first, then higher identity `100(1 − d/len)`, then longer `len`. A full tie keeps the later candidate. `len` is the common length for Hamming and the longer instance for edit.
  - **E-values.** They follow `kurtz/evalues.c` (Kurtz et al. ISMB 2000). The Hamming table P(l, k) is built incrementally with match probability 1/4, stops at 1e-300 and uses multiplier 1. Edit distances are scaled by `averagequot[k]`.
  - **Output.** There is one row per seed, so the same repeat can repeat, as in Vmatch's output. Palindromic rows with i > j are dropped, as `fetchpositions` does.
- **`vmatchCompatible = true`** reproduces stock Vmatch 2.3.1.
  - **Left-extension shortcut.** A left extension stops where it scans another exact match ≥ the seed length. For edit this is `evalentrybackward`; for Hamming it is `extendmismatchesleft`, which returns h − 1.
  - **`-allmax` edit output.** Seeds are taken in Vmatch's order: the `processleafedge`/`processbranch` emission order of `Vmengine/vmatfind.c`, and for `-p` the query start. The match container keeps the distance of the **first** seed that reaches a match, because `matchcontainer` never replaces a stored match with an identical one. So stock prints `(0,0,8,8,3)` for `CAAATTTT -p -l 8 -e 3`, although d_E = 2.
  - **Default (`false`).** The default keeps the complete extension and reports the edit distance.
  - **Hamming `-allmax`.** The shortcut never changes the maximal set.
- **Cross-check against real Vmatch** (Evidence §WP15): 0 mismatching rows, compared as multisets.
  - The reference is stock `vmatch` for `vmatchCompatible = true`. For `false` it is the same release built from source with both shortcuts switched off (`VM_NOPRUNE`, `VM_NOPRUNE_H`).
  - The cases cover all 16 combinations of {h, e} × {direct, -p} × {allmax, best} × {stock, complete}.
  - Inputs: 6 000 + 2 000 (k ≤ 4) inputs of 8–50 bp, 300 of 100–1 500 bp, and 1 Mb and 200 kb inputs.
  - **One documented exception, in the default `allMaximal` edit mode.** That mode reports the true edit distance, whereas Vmatch prints its first-seed label: 1 case in 6 000, `CAAATTTT -p`.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindDirectRepeats(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindDirectRepeats(string, int, int, int)`: Uppercases raw string input and yields results for non-empty strings.
- `RepeatFinder.FindReverseComplementRepeats(DnaSequence|string, int, int, int)` → `ReverseComplementRepeatResult` (§4.4).
- `RepeatFinder.FindApproximateDirectRepeats(DnaSequence|string, int, int, int, int, bool)` → `ApproximateDirectRepeatResult` (§4.4).
- `RepeatFinder.FindSupermaximalRepeats(DnaSequence|string, int)` → `SupermaximalRepeatResult` (§4.4).
- `RepeatFinder.FindDegenerateRepeats(DnaSequence|string, int, int, ApproximateRepeatDistance, bool, int, int)` → `DegenerateRepeatResult` (§4.5).
- `RepeatFinder.FindDegenerateRepeats(DnaSequence|string, …, DegenerateRepeatReporting, bool vmatchCompatible = false)` and `FindApproximateDirectRepeats(string, …, bool excludeContained, DegenerateRepeatReporting, bool vmatchCompatible = false)`: Vmatch best-per-seed output and stock-Vmatch compatibility (§4.6).
- MCP: `find_direct_repeats` wraps `FindDirectRepeats`; the variants are exposed as `find_reverse_complement_repeats`, `find_approximate_direct_repeats`, `find_degenerate_repeats` and `find_supermaximal_repeats` (Analysis server, review-2026-09 B04 F49).

### 5.2 Current Behavior

See §4. The previous implementation (until 2026-09) enumerated every `(i, j, len)` window for every `len ∈ [minLength, maxLength]` via `Substring` + suffix-tree `FindAllOccurrences` + an `(i, j, len)` hash set: one repeat of length L produced O(L²) nested hits (`AAAAAATTTTAAAAAA`, 4–6 → 14 hits instead of 5 maximal pairs), negative `minSpacing` produced self-pairs `(i, i)`, N runs were reported as repeats, and the cost was `O(r · n · (m + k))`. Replaced by the maximal-pair enumeration (B04 F11).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact same-orientation repeats [1][2] reported as maximal repeated pairs [7], identical to MUMmer `repeat-match -f` [5] (0 mismatches on 2 200 random cases incl. 12.8 M pairs, and 50 kb–1 Mb genomes; see Evidence).
- Only A/C/G/T match (MUMmer `-n`) [5].

- Reverse-complement maximal pairs identical to `repeat-match` (without `-f`) and Vmatch `-p` (§4.4; Evidence).
- Maximal k-mismatch repeats identical to an independent brute force of the definition and, with `excludeContained`, to `vmatch -h k -allmax` (§4.4; Evidence).
- Supermaximal repeats identical to Vmatch `-supermax` and a brute force of Gusfield's definition (§4.4; Evidence).
- Degenerate repeats (`-e k`, `-p -e k`, `-p -h k`, `-h k -allmax`) identical to a brute force of the Vmatch App. A definitions (6 000 cases) and to Vmatch 2.3.1 with its left-extension shortcut disabled (6 800 cases + 1 Mb); stock Vmatch differs only in edit mode, only by that shortcut (§4.5; Evidence).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The raw-string overload did not mirror the range validation of the `DnaSequence` overload. | Deviation | `minLength = 0` produced an O(n²) blow-up. | resolved (REP-DIRECT-001 fuzzing) | Both overloads validate eagerly. |
| 2 | Every nested `(i, j, len)` window reported; self-pairs for negative `minSpacing`; N matched N. | Defect | O(L²) redundant hits per repeat; invalid (i, i) pairs; assembly gaps reported as repeats. | resolved (B04 F11, 2026-09-30) | Maximal pairs per MUMmer repeat-match; `i < j` by construction; ACGT-only. |
| 3 | `maxLength` is a filter on maximal repeats (repeat-match has no maximum). | Assumption | A repeat longer than `maxLength` is absent from the output. | documented | Pass `int.MaxValue` to disable. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / `null` string | Empty result. | Contract. |
| No repeat ≥ `minLength` | Empty result. | No qualifying maximal pair. |
| Adjacent copies, `minSpacing = 0` | `ACGTAACGTA` → (0,5,5). | Spacing 0 allowed. |
| Three copies with distinct flanks | All three pairs. | Each pair is maximal. |
| Periodic copies (`ACGTATTACGTATTACGTA`) | Spacing ≥ 1 → (0,14,5) only; overlapping (0,7,12) with negative `minSpacing`. | (0,7) and (7,14) windows are nested in (0,7,12). |
| Homopolymer A^n | Pairs (0, j, n − j) only. | Only position 0 is left-maximal. |
| N / IUPAC symbols | Never match; split copies. | MUMmer `-n`. |
| Lower-case input | Same as upper-case. | Case-insensitive (repeat-match lower-cases). |
| Maximal repeat longer than `maxLength` | Not reported. | Filter, no truncation. |

### 6.2 Limitations

`FindDirectRepeats` / `FindReverseComplementRepeats`: one record per position pair (a repeat present in c copies yields c(c−1)/2 pairs, as in repeat-match); `FindSupermaximalRepeats` groups occurrences per string. `FindApproximateDirectRepeats` uses the Hamming distance (mismatches only); k-differences (edit-distance, Vmatch `-e`) and approximate palindromic (`-p -h` / `-p -e`) repeats are `FindDegenerateRepeats` (§4.5). Small ⌊m/(k+1)⌋ makes the number of seeds grow like n²·4^−⌊m/(k+1)⌋ (as in Vmatch). No biological annotation (LTR, recombination substrate, etc.).

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Transposable elements: direct repeats and long terminal repeats mark some mobile-element architectures [1][2].
- Recombination and genome instability: same-orientation repeats are hotspots for repeat-mediated rearrangements and deletions [2][3].
- Gene regulation: some regulatory elements contain repeated same-orientation sequence motifs [1].
- Trinucleotide-repeat disease context: tandem direct-repeat expansions include the pathogenic motifs reported for Huntington's disease (`CAG`, `HTT`), Fragile X syndrome (`CGG`, `FMR1`), spinocerebellar ataxias (`CAG`, various genes), Friedreich's ataxia (`GAA`, `FXN`), and myotonic dystrophy (`CTG/CCTG`, `DMPK/ZNF9`) [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_DirectRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_DirectRepeat_Tests.cs), [RepeatFinder_RepeatVariants_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_RepeatVariants_Tests.cs); heavy tier `Properties/RepDirectVariantsProperties.cs`, `Metamorphic/RepDirectVariantsMetamorphicTests.cs`, `Fuzzing/RepDirectVariantsFuzzTests.cs`
- Test spec: [REP-DIRECT-001.md](../../../tests/TestSpecs/REP-DIRECT-001.md)
- Related snapshot tests: [RepeatSnapshotTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/RepeatSnapshotTests.cs)
- Evidence: [REP-DIRECT-001-Evidence.md](../../Evidence/REP-DIRECT-001-Evidence.md)

## 8. References

1. Wikipedia. 2026. Direct repeat. Wikipedia. https://en.wikipedia.org/wiki/Direct_repeat
2. Wikipedia. 2026. Repeated sequence (DNA). Wikipedia. https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)
3. Ussery DW, Wassenaar TM, Borini S. 2009. Computing for Comparative Microbial Genomics. Chapter 8.
4. Richard GF. 2021. Trinucleotide repeat expansions and human disease. PMC8145212. https://pmc.ncbi.nlm.nih.gov/articles/PMC8145212/
5. Kurtz S, Phillippy A, Delcher AL, Smoot M, Shumway M, Antonescu C, Salzberg SL. 2004. Versatile and open software for comparing large genomes. Genome Biol 5:R12. Source: https://github.com/mummer4/mummer (`src/tigr/repeat-match.cc`, `src/essaMEM/mummer.cpp`).
6. Kurtz S, Schleiermacher C. 1999. REPuter: fast computation of maximal repeats in complete genomes. Bioinformatics 15(5):426–427.
7. Gusfield D. 1997. Algorithms on Strings, Trees, and Sequences. Cambridge University Press. §7.12 (maximal pairs, maximal repeats).
8. Abouelhoda MI, Kurtz S, Ohlebusch E. 2004. Replacing suffix trees with enhanced suffix arrays. J Discrete Algorithms 2:53–86.
9. Kurtz S, Choudhuri JV, Ohlebusch E, Schleiermacher C, Stoye J, Giegerich R. 2001. REPuter: the manifold applications of repeat analysis on a genomic scale. Nucleic Acids Res 29(22):4633–4642.
10. Kurtz S. The Vmatch large scale sequence analysis software — a manual (Vmatch 2.3.1, ISC licence; Debian/Ubuntu `vmatch` source package `vstree-2.3.1/src/doc/virtman.tex`): options `-p`, `-h`, `-allmax`, `-seedlength`, `-supermax`; Appendix A "Basic Notions" (palindromic match, k-mismatch match, maximality by containment, supermaximal repeat).
11. Vmatch 2.3.1 source (`vstree-2.3.1`, Ubuntu `vmatch_2.3.1+dfsg.orig.tar.xz`): `kurtz/extendED.c` (`editextend`, `acceptmatch`), `kurtz/frontSEP.c` (greedy fronts, `evalentrybackward` seed shortcut), `kurtz/mcontain.c` (`matchcontainer`), `Vmengine/fself.c` (per-seed then global containment for `-allmax`); options `-e`, `-allmax`, `-seedlength` (`virtman.tex`).
12. Ukkonen E. 1985. Algorithms for approximate string matching. Information and Control 64:100–118 (furthest-reaching diagonal fronts). Myers EW. 1986. An O(ND) difference algorithm and its variations. Algorithmica 1:251–266.
