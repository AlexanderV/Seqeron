# Palindrome Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Repeat Analysis |
| Test Unit ID | REP-PALIN-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

DNA palindrome detection identifies sequences that are equal to their own reverse complement [1]. This is a biological notion of palindrome, distinct from textual symmetry on a single strand, and it is central to restriction-enzyme recognition-site discovery [1][2]. The repository exposes a validating implementation in `RepeatFinder.FindPalindromes` and a second, lighter-weight implementation in `GenomicAnalyzer.FindPalindromes` with a different return type. Both implementations perform exact reverse-complement comparison over even-length windows.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A DNA palindrome reads the same in the 5'→3' direction on the forward strand as it does in the 5'→3' direction on the reverse strand [1]. For example, `GAATTC` is the EcoRI recognition sequence:

```text
5'-GAATTC-3'
3'-CTTAAG-5'
     ↓ read reverse strand 5'→3'
5'-GAATTC-3'
```

DNA palindromes are biologically important because they are common restriction-enzyme recognition sites, can form cruciform and hairpin structures, and are associated with chromosome fragility when long self-complementary tracts occur in genomes [1][2].

### 2.2 Core Model

For a candidate sequence $S$, palindrome detection uses the criterion:

$$
\operatorname{Palindrome}(S) \iff S = \operatorname{ReverseComplement}(S)
$$

For an even-length sequence $S = s_1 s_2 \dots s_{2n}$, the equivalent basewise condition is:

$$
s_i = \operatorname{complement}(s_{2n + 1 - i}) \quad \text{for all } i \in [1, 2n]
$$

with complement mapping `A↔T` and `G↔C` [1]. Biological DNA palindromes therefore require even length so that every base pairs with a complementary partner across the symmetry axis.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported palindrome satisfies `Sequence = ReverseComplement(Sequence)`. | Both implementations emit results only after exact reverse-complement comparison. |
| INV-02 | All reported lengths are even. | The scan steps through candidate lengths in increments of `2`. |
| INV-03 | Reported positions are within sequence bounds. | Candidate windows are enumerated only when `start <= seq.Length - len`. |
| INV-04 | `Length` equals the actual sequence length of the reported palindrome. | `PalindromeResult` stores the scanned candidate and its length directly. |

### 2.5 Comparison with Related Concepts

| Concept | Relation |
|---------|----------|
| Inverted repeat | A palindrome is the zero-loop special case of an inverted repeat. |
| Restriction site analysis | Palindrome detection is a structural prerequisite for recognizing many Type II restriction-enzyme sites [2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to search. | `RepeatFinder` throws on `null` `DnaSequence` input and yields no results for `null` or empty raw strings. |
| `minLength` | `int` | `4` | Minimum palindrome length. | Must be even and at least `4` in `RepeatFinder`. |
| `maxLength` | `int` | `12` | Maximum palindrome length. | Must be at least `minLength` in `RepeatFinder`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | 0-based start position of the palindrome. |
| `Sequence` | `string` | Palindromic DNA sequence. |
| `Length` | `int` | Length of the palindrome in bases. |

`GenomicAnalyzer.FindPalindromes` returns `PalindromeInfo` instead of `PalindromeResult`, but it still reports the exact palindrome sequence and its start position.

### 3.3 Preconditions and Validation

`RepeatFinder.FindPalindromes(DnaSequence, ...)` throws `ArgumentNullException` on `null` sequence input, throws `ArgumentOutOfRangeException` when `minLength < 4` or `minLength` is odd, and throws `ArgumentOutOfRangeException` when `maxLength < minLength`. The raw-string overload applies the same length validation, yields no results for `null` or empty strings, and uppercases non-empty input before scanning. `GenomicAnalyzer.FindPalindromes` performs the same scan logic but does not apply those explicit argument checks.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `minLength` and `maxLength` in `RepeatFinder`.
2. Normalize raw-string input to uppercase when needed.
3. For each even candidate length from `minLength` through `maxLength`, enumerate every start position where a full window fits.
4. Extract the candidate subsequence, compute its reverse complement, and compare for exact equality.
5. Emit a palindrome result for every matching window.

### 4.2 Reference Table

Common restriction-enzyme palindrome examples cited in the legacy documentation are:

| Enzyme | Sequence | Length | Source Organism |
|--------|----------|--------|-----------------|
| EcoRI | `GAATTC` | 6 | *Escherichia coli* |
| BamHI | `GGATCC` | 6 | *Bacillus amyloliquefaciens* |
| HindIII | `AAGCTT` | 6 | *Haemophilus influenzae* |
| TaqI | `TCGA` | 4 | *Thermus aquaticus* |
| AluI | `AGCT` | 4 | *Arthrobacter luteus* |
| NotI | `GCGGCCGC` | 8 | *Nocardia otitidis* |
| SmaI | `CCCGGG` | 6 | *Serratia marcescens* |
| EcoRV | `GATATC` | 6 | *Escherichia coli* |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Palindrome detection | `O(n × r × m)` | `O(m)` | `r` is the number of even candidate lengths scanned and `m` is the cost of reverse-complement comparison for a candidate window. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation locations:** [RepeatFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs), [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `RepeatFinder.FindPalindromes(DnaSequence, int, int)`: Validating canonical API returning `PalindromeResult`.
- `RepeatFinder.FindPalindromes(string, int, int)`: Raw-string overload with uppercase normalization.
- `GenomicAnalyzer.FindPalindromes(DnaSequence, int, int)`: Alternate API returning `PalindromeInfo`.

### 5.2 Current Behavior

`RepeatFinder` enforces `minLength >= 4`, requires even `minLength`, requires `maxLength >= minLength`, and uppercases raw-string input before scanning. `GenomicAnalyzer.FindPalindromes` performs the same even-length stepwise scan and reverse-complement equality test but does not add explicit null or range validation before it dereferences `sequence.Sequence`. Both implementations report overlapping palindromes when different even lengths match at the same or neighboring positions.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact biological palindrome detection using `Sequence == ReverseComplement(Sequence)` [1].
- Even-length scanning appropriate for self-complementary DNA palindromes [1].
- Detection of restriction-enzyme-style palindromic sites such as EcoRI, BamHI, and HindIII when those windows are present [2][4].

**Intentionally simplified:**

- Exact DNA matching only, without IUPAC degeneracy or ambiguity-code handling; **consequence:** degenerate recognition motifs are not reported unless the sequence resolves to an exact palindrome.

**Not implemented:**

- Restriction-enzyme catalog lookup and cleavage semantics; **users should rely on:** [RestrictionAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/RestrictionAnalyzer.cs) for restriction-site analysis beyond structural palindrome discovery.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `GenomicAnalyzer.FindPalindromes` does not perform the explicit validation that `RepeatFinder.FindPalindromes` applies. | Deviation | Invalid lengths or `null` input can fail differently depending on the API used. | accepted | The algorithmic core is equivalent, but the validation surfaces differ. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns empty enumerable from string overloads. | No candidate window fits. |
| Sequence shorter than `minLength` | Returns empty enumerable. | No valid palindrome window can be formed. |
| No palindromes present | Returns empty enumerable. | No candidate equals its reverse complement. |
| Entire sequence is palindromic | Reported at the corresponding start position if length constraints permit it. | The window equality test succeeds for the full sequence. |
| Overlapping palindromes of different lengths | Both can be reported. | The scan checks every even window independently. |
| Odd `minLength` | `RepeatFinder` throws `ArgumentOutOfRangeException`. | Biological palindromes require even length, and the validating API enforces it. |
| `minLength < 4` | `RepeatFinder` throws `ArgumentOutOfRangeException`. | The validating API excludes trivial 2-bp matches. |
| `maxLength < minLength` | `RepeatFinder` throws `ArgumentOutOfRangeException`. | Ordered bounds are required. |

### 6.2 Limitations

The algorithm detects only exact DNA palindromes and does not interpret ambiguous IUPAC symbols, cleavage positions, or enzyme-specific recognition context. Long palindromic tracts may have secondary-structure consequences, but the implementation reports only the local exact windows and does not model cruciform energetics or enzyme digestion behavior.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RepeatFinder_Palindrome_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RepeatFinder_Palindrome_Tests.cs)
- Test spec: [REP-PALIN-001.md](../../../tests/TestSpecs/REP-PALIN-001.md)
- Related smoke tests: [GenomicAnalyzerTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GenomicAnalyzerTests.cs)

## 8. References

1. Wikipedia. 2026. Palindromic sequence. Wikipedia. https://en.wikipedia.org/wiki/Palindromic_sequence
2. Wikipedia. 2026. Restriction enzyme. Wikipedia. https://en.wikipedia.org/wiki/Restriction_enzyme
3. Rosalind. 2026. REVP: Locating Restriction Sites. Rosalind. https://rosalind.info/problems/revp/
4. REBASE. 2026. Restriction Enzyme Database. https://rebase.neb.com/
