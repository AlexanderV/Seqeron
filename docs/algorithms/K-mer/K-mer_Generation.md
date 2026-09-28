# K-mer Generation

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-GENERATE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

K-mer Generation enumerates *every* possible k-mer of length `k` over a given alphabet — the complete k-mer universe, independent of any particular sequence. For an alphabet of size `n` there are exactly `n^k` such k-mers, formed as the k-fold Cartesian product of the alphabet [1][2]. This is an exact, deterministic, combinatorial enumeration used to build frequency arrays / fixed-index k-mer tables, to initialise background models, and to iterate the space of possible motifs. It is distinct from k-mer *counting*, which extracts the k-mers actually present in a sequence.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a length-`k` string over a sequence alphabet ("substrings of length k contained within a biological sequence") [1]. For DNA the alphabet is {A,C,G,T}; for protein it is the 20 standard amino acids. Many k-mer methods need the full set of possible k-mers (e.g. to address a frequency array of size `n^k`), not only those observed in a sample.

### 2.2 Core Model

The set of all possible k-mers over alphabet `Σ` (|Σ| = n) is the k-fold Cartesian product `Σ^k`. Its cardinality is `n^k`: "there exist n^k total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)" [1]; equivalently "the possible combinations of k positions are computed as 4^k" for DNA [2]. Enumeration is the Cartesian product `product(Σ, repeat=k)` [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output count = `n^k` (4^k for default DNA alphabet) | Cardinality of the k-fold Cartesian product `Σ^k` [1][2] |
| INV-02 | All emitted k-mers are distinct (output is a set of size `n^k`) | Each k-mer is a unique length-k tuple over Σ [3] |
| INV-03 | Every k-mer has length exactly `k` and contains only alphabet characters | Direct from the definition of a k-mer over Σ [1] |
| INV-04 | Emission order is lexicographic with respect to the alphabet's own symbol order (rightmost position advances fastest); for a sorted alphabet this is ordinary sorted order | Rosalind LEXF: the alphabet is an ordered permutation a1<a2<…, s <Lex t iff the first mismatching symbol of s precedes that of t [4]; odometer ordering of the Cartesian product: "if the input's iterables are sorted, the product tuples are emitted in sorted order" [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `k` | `int` | required | k-mer length | must be > 0 |
| `alphabet` | `string` | `"ACGT"` | symbols to combine | non-empty; characters used verbatim |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IEnumerable<string>` | All `alphabet.Length^k` k-mers (distinct when the alphabet symbols are distinct); lexicographic in the alphabet's order (Rosalind LEXF) |

### 3.3 Preconditions and Validation

`k <= 0` throws `ArgumentOutOfRangeException` (a k-mer length must be positive [1]). A null or empty `alphabet` throws `ArgumentException`. The alphabet is used exactly as supplied (case-sensitive, no normalisation); the default `"ACGT"` is upper-case and already in sorted order, so default DNA output is lexicographic (INV-04). Output is deterministic.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `k > 0` and a non-empty `alphabet`.
2. Enumerate the k-fold Cartesian product with an index odometer (the CPython `itertools.product` algorithm): keep one alphabet index per position and a k-character buffer; after each emission advance the rightmost index, carrying leftwards on wrap-around.
3. Emit each completed length-`k` string.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GenerateAllKmers(k, Σ)` | O(k · n^k) | O(k) working (lazy enumeration) | n = `Σ.Length`; output size is n^k strings of length k |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.GenerateAllKmers(int k, string alphabet = "ACGT")`: returns all `n^k` k-mers.
- `KmerAnalyzer.EnumerateCartesianProduct(...)` (private): index odometer realising the Cartesian product in O(k) working space (replaced the recursive prefix extension in review 2026-09, which kept all k prefixes alive — O(k²) chars — and passed each k-mer through k nested iterators).

### 5.2 Current Behavior

The odometer advances the rightmost position fastest, i.e. odometer ordering; output is byte-identical to `itertools.product(alphabet, repeat=k)` (cross-checked 2026-09-28 for ACGT k=1,5; TAGC k=4; 20-letter protein k=3; AB k=10; AAC k=3; acgT k=3; A k=7). Enumeration is lazy (`yield return`), so callers can stream the universe without materialising all `n^k` strings. The alphabet is taken verbatim — no sorting or de-duplication is applied, so a caller passing an unsorted alphabet receives all k-mers in that alphabet's positional order (INV-04 applies only to sorted alphabets).

**Search-reuse decision:** N/A — this unit *generates* the k-mer universe; it performs no substring search against a text, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Full k-fold Cartesian product `Σ^k` with exactly `n^k` distinct k-mers [1][2][3] (INV-01, INV-02, INV-03).
- Lexicographic / odometer emission order for a sorted alphabet [3] (INV-04).

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `k = 1` | one k-mer per alphabet symbol (DNA → A,C,G,T) | n^1 = n [1] |
| single-letter alphabet, any k | exactly one k-mer (homopolymer) | 1^k = 1 [3] |
| `k <= 0` | `ArgumentOutOfRangeException` | k-mer length must be positive [1]; Rosalind LEXF takes a *positive* integer n [4]. Note: `itertools.product(Σ, repeat=0)` yields one empty tuple (|Σ|^0 = 1); the library deliberately rejects k = 0 as it is not a k-mer length |
| null/empty alphabet | `ArgumentException` | no symbols ⇒ no k-mers |
| unsorted alphabet | all n^k k-mers, lexicographic in the alphabet's order (e.g. "TAGC", k=2 → TT, TA, TG, TC, AT, …, CC) | Rosalind LEXF sample [4]; itertools.product [3] |
| repeated symbol in alphabet | n^k strings including repeats ("AAC", k=2 → AA, AA, AC, AA, AA, AC, CA, CA, CC) | itertools.product semantics [3] |

### 6.2 Limitations

Output size grows as `n^k` and becomes impractical to materialise for large `k` (e.g. DNA k≥14 exceeds 2.6×10^8 k-mers); callers should stream the lazy enumerable. The alphabet is not de-duplicated, so a repeated character yields repeated k-mers — supply a clean, sorted alphabet for lexicographic, duplicate-free output.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// All 16 (4^2) DNA 2-mers in lexicographic order:
// AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
IEnumerable<string> twoMers = KmerAnalyzer.GenerateAllKmers(2);
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_GenerateAllKmers_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_GenerateAllKmers_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [KMER-GENERATE-001-Evidence.md](../../../docs/Evidence/KMER-GENERATE-001-Evidence.md)

## 8. References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo BJ. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Python Software Foundation. 2026. *itertools — Functions creating iterators for efficient looping* (itertools.product). Python 3 Standard Library documentation. https://docs.python.org/3/library/itertools.html
4. Rosalind. *LEXF — Enumerating k-mers Lexicographically*. https://rosalind.info/problems/lexf/ (sample dataset "T A G C", n=2 and its output, read via the archived statement in github.com/mtarbit/Rosalind-Problems `e015-lexf.py`).
