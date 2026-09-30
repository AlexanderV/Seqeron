# Direct Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-DIRECT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Complete (exact repeats) |
| Last Reviewed | 2026-09-30 |

## 1. Overview

Direct repeat detection identifies nucleotide sequences that recur in the same 5'→3' orientation at multiple genomic positions [1][2]. Unlike inverted repeats, the downstream copy preserves the original sequence rather than its reverse complement. The repository implements exact direct-repeat discovery in `RepeatFinder.FindDirectRepeats`, reporting **maximal repeated pairs** — the convention of the reference tools MUMmer `repeat-match` [5] and REPuter [6] and of Gusfield's maximal-pair definition [7]. Spacing between copies is configurable, so adjacent tandem-like repeats, separated direct repeats and (with negative spacing) overlapping copies can all be reported.

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

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindDirectRepeats(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindDirectRepeats(string, int, int, int)`: Uppercases raw string input and yields results for non-empty strings.

### 5.2 Current Behavior

See §4. The previous implementation (until 2026-09) enumerated every `(i, j, len)` window for every `len ∈ [minLength, maxLength]` via `Substring` + suffix-tree `FindAllOccurrences` + an `(i, j, len)` hash set: one repeat of length L produced O(L²) nested hits (`AAAAAATTTTAAAAAA`, 4–6 → 14 hits instead of 5 maximal pairs), negative `minSpacing` produced self-pairs `(i, i)`, N runs were reported as repeats, and the cost was `O(r · n · (m + k))`. Replaced by the maximal-pair enumeration (B04 F11).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact same-orientation repeats [1][2] reported as maximal repeated pairs [7], identical to MUMmer `repeat-match -f` [5] (0 mismatches on 2 200 random cases incl. 12.8 M pairs, and 50 kb–1 Mb genomes; see Evidence).
- Only A/C/G/T match (MUMmer `-n`) [5].

**Not implemented (different algorithms, not simplifications):**

- Reverse-complement (palindromic) repeats — `repeat-match` without `-f`; use `FindInvertedRepeats` for inverted repeats.
- Approximate (mismatch/indel) repeats (REPuter degenerate repeats [6]); use alignment-based tools.
- Supermaximal repeats / repeat families (one record per repeated string rather than per pair).

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

Exact repeats only (no mismatches/indels); forward strand only; one record per position pair (a repeat present in c copies yields c(c−1)/2 pairs, as in repeat-match). No biological annotation (LTR, recombination substrate, etc.).

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Transposable elements: direct repeats and long terminal repeats mark some mobile-element architectures [1][2].
- Recombination and genome instability: same-orientation repeats are hotspots for repeat-mediated rearrangements and deletions [2][3].
- Gene regulation: some regulatory elements contain repeated same-orientation sequence motifs [1].
- Trinucleotide-repeat disease context: tandem direct-repeat expansions include the pathogenic motifs reported for Huntington's disease (`CAG`, `HTT`), Fragile X syndrome (`CGG`, `FMR1`), spinocerebellar ataxias (`CAG`, various genes), Friedreich's ataxia (`GAA`, `FXN`), and myotonic dystrophy (`CTG/CCTG`, `DMPK/ZNF9`) [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_DirectRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_DirectRepeat_Tests.cs)
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
