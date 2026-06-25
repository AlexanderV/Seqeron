# Direct Repeat Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-DIRECT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Direct repeat detection identifies nucleotide sequences that recur in the same 5'→3' orientation at multiple genomic positions [1][2]. Unlike inverted repeats, the downstream copy preserves the original sequence rather than its reverse complement. The repository implements exact direct-repeat discovery in `RepeatFinder.FindDirectRepeats` by combining length and position enumeration with suffix-tree occurrence lookup. Spacing between copies is configurable, so adjacent tandem-like repeats and separated direct repeats can both be reported.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A direct repeat consists of two or more identical sequence copies oriented in the same 5'→3' direction [1]. A canonical representation is:

```text
5' TTACG------TTACG 3'
3' AATGC------AATGC 5'
```

where `------` is an intervening spacer that may be zero bases long [1]. Direct repeats are biologically associated with transposable-element boundaries, homologous recombination, replication-slippage-mediated deletions, and some regulatory architectures [1][2][3]. The legacy reference set also notes that tandem trinucleotide expansions are a special case of direct-repeat biology and underlie disorders such as Huntington's disease, Fragile X syndrome, spinocerebellar ataxias, Friedreich's ataxia, and myotonic dystrophy [4].

### 2.2 Core Model

For a sequence $S$, repeat length $L$, first position $i$, and second position $j$, a direct-repeat pair is present when:

$$
S[i..i+L) = S[j..j+L) \quad \text{and} \quad j > i + L - 1 + \text{minSpacing}
$$

The implementation reports spacing as:

$$
	ext{Spacing} = j - i - L
$$

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `RepeatSequence` is identical at `FirstPosition` and `SecondPosition`. | Results are emitted only from suffix-tree occurrence matches of the same extracted pattern. |
| INV-02 | `Spacing = SecondPosition - FirstPosition - Length`. | Spacing is constructed directly from the stored coordinates. |
| INV-03 | When `minSpacing > 0`, reported copies do not overlap. | The filter requires `j > i + len - 1 + minSpacing`. |
| INV-04 | Each `(FirstPosition, SecondPosition, Length)` tuple is unique. | A hash set suppresses duplicate result keys. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | The `DnaSequence` overload throws on `null`; the `string` overload returns no results for `null` or empty input. |
| `minLength` | `int` | `5` | Minimum repeat length to test. | The `DnaSequence` overload rejects values below `2`. |
| `maxLength` | `int` | `50` | Maximum repeat length to test. | The `DnaSequence` overload rejects values below `minLength`. |
| `minSpacing` | `int` | `1` | Minimum number of bases between copies. | `0` is allowed and enables adjacent repeat copies. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `FirstPosition` | `int` | 0-based position of the first repeat copy. |
| `SecondPosition` | `int` | 0-based position of the downstream repeat copy. |
| `RepeatSequence` | `string` | Exact repeated sequence shared by both copies. |
| `Length` | `int` | Length of the repeat sequence in bases. |
| `Spacing` | `int` | Number of nucleotides between the two copies. |

### 3.3 Preconditions and Validation

`FindDirectRepeats(DnaSequence, ...)` throws `ArgumentNullException` when `sequence` is `null`, throws `ArgumentOutOfRangeException` when `minLength < 2`, and throws `ArgumentOutOfRangeException` when `maxLength < minLength`. The raw-string overload uppercases non-empty input and yields no results for `null` or empty strings; it now ALSO mirrors the numeric range checks of the `DnaSequence` overload (`minLength < 2` and `maxLength < minLength` both throw `ArgumentOutOfRangeException`), validating eagerly at the call rather than only on enumeration. Positions are 0-based in all returned results.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase when the raw-string overload is used.
2. Build a suffix tree for the full sequence.
3. For each candidate repeat length from `minLength` to `maxLength`, enumerate start positions where two copies plus the required spacing could still fit.
4. Extract the candidate repeat at the current start position and query the suffix tree for all its occurrences.
5. Keep downstream occurrences that satisfy the spacing constraint, suppress duplicate `(start, occurrence, length)` tuples, and emit a `DirectRepeatResult` for each qualifying pair.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Suffix-tree construction | `O(n)` | `O(n)` | The full input sequence is indexed once before candidate enumeration. |
| Repeat detection | `O(r × n × (m + k))` | `O(d)` plus suffix tree | `r` is the number of tested lengths, `m` is the candidate repeat length, `k` is the number of suffix-tree occurrences returned for that candidate, and `d` is the number of dedup keys retained. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs)

- `RepeatFinder.FindDirectRepeats(DnaSequence, int, int, int)`: Validating overload for `DnaSequence` input.
- `RepeatFinder.FindDirectRepeats(string, int, int, int)`: Uppercases raw string input and yields results for non-empty strings.

### 5.2 Current Behavior

The core implementation builds a suffix tree with `global::SuffixTree.SuffixTree.Build(seq)` and uses `FindAllOccurrences(repeat)` to locate repeated candidates. For a single starting position and repeat length, it can emit multiple result pairs if the same sequence occurs at multiple later positions that satisfy the spacing filter. Duplicate `(FirstPosition, SecondPosition, Length)` tuples are removed with a hash set. The raw-string overload normalizes case with `ToUpperInvariant()`, while the `DnaSequence` overload preserves the already-normalized `DnaSequence.Sequence` content.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact direct-repeat detection with same-orientation copies rather than reverse complements [1][2].
- Configurable spacing so both adjacent and separated direct repeats can be represented [1].
- Pairwise reporting of repeat coordinates, repeat sequence, length, and spacing for each qualifying repeat pair.

**Intentionally simplified:**

- Exact matching only, with no mismatch or gap tolerance; **consequence:** approximate, degenerate, or interrupted direct repeats are not reported.
- Raw repeat-pair output only, without higher-level annotation of transposable elements or regulatory context; **consequence:** downstream interpretation remains the caller's responsibility.

**Not implemented:**

- Approximate repeat scoring or specialized LTR/transposon annotation; **users should rely on:** `RepeatFinder.FindDirectRepeats` only for exact direct-repeat pairs.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The raw-string overload did not mirror the range validation performed by the `DnaSequence` overload. | Deviation | A degenerate `minLength = 0` on the raw-string surface produced a zero-length candidate whose suffix-tree lookup matched every position, blowing the result set up with `O(n²)` spurious empty-/single-base "repeats". | resolved (REP-DIRECT-001 fuzzing) | Fixed by mirroring the `DnaSequence` overload's `minLength < 2` and `maxLength < minLength` guards onto the raw-string overload (hoisted into an eager wrapper so the exception surfaces at the call). Both overloads now validate identically; the MCP `find_direct_repeats` tool that forwards raw user input is no longer exposed to the blow-up. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns empty enumerable. | The raw-string overload yields nothing for empty input, and a too-short indexed sequence produces no qualifying windows. |
| Sequence too short for two copies | Returns empty enumerable. | The inner position loop requires room for two copies plus `minSpacing`. |
| No repeats found | Returns empty enumerable. | No candidate pattern has a qualifying downstream occurrence. |
| Adjacent repeats with `minSpacing = 0` | Reported when identical copies abut exactly. | The filter permits `j = i + len` when `minSpacing` is zero. |
| Multiple later occurrences of the same repeat | One result per qualifying pair. | The occurrence loop emits every downstream match that satisfies the spacing filter. |
| Case variation in raw strings | Normalized to uppercase before detection. | The string overload applies `ToUpperInvariant()`. |

### 6.2 Limitations

The algorithm reports only exact direct repeats and does not score partial similarity, mismatches, or indels. It also does not perform biological annotation of why a repeat is present, so users must interpret whether a reported direct repeat corresponds to an LTR, a tandem expansion, a recombination substrate, or another structure. Numeric validation is now symmetric across both overloads (`minLength < 2` and `maxLength < minLength` throw `ArgumentOutOfRangeException` on the raw-string surface as well as the `DnaSequence` surface).

## 7. Examples and Related Material

### 7.2 Related Use Cases

- Transposable elements: direct repeats and long terminal repeats mark some mobile-element architectures [1][2].
- Recombination and genome instability: same-orientation repeats are hotspots for repeat-mediated rearrangements and deletions [2][3].
- Gene regulation: some regulatory elements contain repeated same-orientation sequence motifs [1].
- Trinucleotide-repeat disease context: tandem direct-repeat expansions include the pathogenic motifs reported for Huntington's disease (`CAG`, `HTT`), Fragile X syndrome (`CGG`, `FMR1`), spinocerebellar ataxias (`CAG`, various genes), Friedreich's ataxia (`GAA`, `FXN`), and myotonic dystrophy (`CTG/CCTG`, `DMPK/ZNF9`) [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_DirectRepeat_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/RepeatFinder_DirectRepeat_Tests.cs)
- Test spec: [REP-DIRECT-001.md](../../../tests/TestSpecs/REP-DIRECT-001.md)
- Related snapshot tests: [RepeatSnapshotTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/RepeatSnapshotTests.cs)

## 8. References

1. Wikipedia. 2026. Direct repeat. Wikipedia. https://en.wikipedia.org/wiki/Direct_repeat
2. Wikipedia. 2026. Repeated sequence (DNA). Wikipedia. https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)
3. Ussery DW, Wassenaar TM, Borini S. 2009. Computing for Comparative Microbial Genomics. Chapter 8.
4. Richard GF. 2021. Trinucleotide repeat expansions and human disease. PMC8145212. https://pmc.ncbi.nlm.nih.gov/articles/PMC8145212/
