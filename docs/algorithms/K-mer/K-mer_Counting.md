# K-mer Counting

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-COUNT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 |

## 1. Overview

K-mer counting extracts and counts all overlapping substrings of length `k` from a biological sequence. It is a foundational operation for genome assembly, metagenomics binning, sequence comparison, and repeat analysis. In this repository, the implementation provides string, `DnaSequence`, span-based, cancellation-aware, async, and both-strand counting surfaces. The counting logic is DNA-oriented in its framing, but the raw-string and span-based surfaces treat any uppercase symbols as literal k-mer characters rather than filtering to `A/C/G/T`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a substring of fixed length `k` contained within a sequence. For a sequence of length `L`, there are `L - k + 1` overlapping k-mers when `k <= L`. For a DNA alphabet of size 4, the number of possible distinct k-mers is `4^k`. Sources: Wikipedia (K-mer), Rosalind (K-mer Composition), Compeau et al. (2011), Marçais & Kingsford (2011).

### 2.2 Core Model

The repository uses the standard sliding-window count model:

$$
Count(kmer) = \sum_{i=0}^{L-k} \mathbf{1}(sequence[i..i+k-1] = kmer)
$$

and stores counts in a dictionary keyed by the uppercased k-mer string.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The sum of all counts equals `L - k + 1` whenever `k <= L` | The sliding window emits exactly one k-mer per valid start position |
| INV-02 | For DNA-alphabet inputs, the number of unique k-mers is at most `min(4^k, L - k + 1)` | There are at most `4^k` DNA words and at most one emitted k-mer per window |
| INV-03 | The implementation is case-insensitive | Input is normalized to uppercase before counting |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string`, `DnaSequence`, or `ReadOnlySpan<char>` | required | Sequence to analyze | Null or empty string returns an empty dictionary |
| `k` | `int` | required | K-mer length | For non-empty string input and span-based counting, `k <= 0` throws `ArgumentOutOfRangeException`; `k > sequence.Length` returns an empty dictionary |
| `cancellationToken` | `CancellationToken` | optional | Cancellation support for long-running counting | Used in cancellation-aware overloads |
| `progress` | `IProgress<double>?` | optional | Progress reporter for long-running counting | Reports values between `0.0` and `1.0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `counts` | `Dictionary<string, int>` | Maps each observed k-mer to its count |

### 3.3 Preconditions and Validation

String-based counting returns an empty dictionary for null or empty input before validating `k`. The `DnaSequence` overloads delegate to that string-based path, so an empty `DnaSequence.Sequence` also returns an empty dictionary before `k` is checked. For non-empty string input, the cancellation-aware and synchronous string overloads throw `ArgumentOutOfRangeException` when `k <= 0`, and the span-based path validates `k` before checking sequence length. When `k` is greater than the sequence length, counting returns an empty dictionary. The string and cancellation-aware overloads uppercase the input before emitting dictionary keys.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input sequence and `k`.
2. Normalize the sequence to uppercase when using string-based overloads.
3. Slide a window of length `k` across every valid start position.
4. Increment the count for the observed k-mer in the output dictionary.
5. Return the completed count map.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Entry points documented in the original file and confirmed in source:

| Method | Class | Description |
|--------|-------|-------------|
| `CountKmers(string, k)` | `KmerAnalyzer` | Canonical string-based counting |
| `CountKmers(DnaSequence, k)` | `KmerAnalyzer` | Wrapper for `DnaSequence` input |
| `CountKmersSpan(ReadOnlySpan<char>, k)` | `KmerAnalyzer` | Span-based counting variant (the class's single loop; `SequenceExtensions.CountKmersSpan` in Core is a separate B01 helper) |
| `CountKmersParallel(string, k, options, dop, ct, progress)` | `KmerAnalyzer` | Opt-in multi-core count, result identical to `CountKmers` (see Asynchronous_K-mer_Counting.md §5.3) |
| `CountKmersBothStrands(DnaSequence, k)` | `KmerAnalyzer` | Forward plus reverse-complement counting |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CountKmers` | `O(n × k)` effective | `O(u)` | `n` is sequence length, `k` is k-mer length, and each window creates and hashes a length-`k` string key |
| `CountKmersBothStrands` | `O(n × k)` effective | `O(u)` | Counts forward and reverse-complement sequences separately, with the same per-window string allocation and hashing cost |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs), [SequenceExtensions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs)

- `KmerAnalyzer.CountKmers(...)`: Canonical counting overloads for strings and `DnaSequence` values.
- `KmerAnalyzer.CountKmers(string, int, KmerCountingOptions, CancellationToken = default, IProgress<double>? = null)`: option-aware counting — literal (default), ACGT-only (Jellyfish window rule) or canonical (Jellyfish `count -C`).
- `KmerAnalyzer.DistinctKmers(string, int[, KmerCountingOptions])`: the distinct k-mer set (key set of the option-aware counter; caller-owned ordinal `HashSet<string>`).
- `KmerAnalyzer.CountKmersAsync(...)`: Async wrapper over cancellation-aware counting.
- `KmerAnalyzer.CountKmersSpan(ReadOnlySpan<char>, int)`: Span-based counting entry point.
- `KmerAnalyzer.CountKmersBothStrands(DnaSequence, int)`: Counts forward and reverse-complement k-mers and combines the totals.

### 5.2 Current Behavior

All string-based entry points uppercase the input before counting. The synchronous `CountKmers(string, int)` delegates to the cancellation-aware overload with `CancellationToken.None`, so there is a single counting loop in `KmerAnalyzer`. The synchronous `CountKmers(...)` methods are stateless, and the cancellation-aware overload checks cancellation periodically while optionally reporting progress. `CountKmersSpan(...)` copies the span once and runs the same loop, with the same contract as `CountKmers` (empty input → empty for any k; k ≤ 0 throws only for non-empty input); it no longer delegates to Core's `SequenceExtensions.CountKmersSpan`, which throws for k ≤ 0 even on an empty span and allocates one (upper-case input) or two (lower-case input) strings per window. `CountKmersBothStrands(...)` counts forward and reverse-complement sequences independently before summing counts. The single loop (`CountWindowRange`, audit round 1 WP4) looks each window up by `ReadOnlySpan<char>` through the dictionary's alternate lookup, so a string key is allocated only for the first occurrence of each distinct k-mer; runtime is the O(n·k) hashing of the windows plus O(u·k) key allocation. Measured on random 10 Mbp (Release): k = 8 — `KmerAnalyzer.CountKmers` 0.43–0.45 s / 7 MB allocated, Core `CountKmersSpan` 1.17–1.24 s / 404 MB (lower-case input: 0.45 s / 27 MB vs 1.46 s / 804 MB); k = 12 — 3.3 s vs 4.1–4.5 s; k = 21 — 3.4–3.9 s vs 4.6–4.9 s. The `DnaSequence` overloads throw `ArgumentNullException` (ParamName `dna`) for a null argument. The raw-string and span-based paths do not restrict the alphabet, so ambiguous or non-ACGT symbols are preserved as literal k-mer keys.

The option-aware overload uses the same single loop (`CountKmersCore`). With `KmerCountingOptions.AcgtOnly` it skips every window that contains a symbol other than A/C/G/T after upper-casing; the last non-ACGT index is tracked so the check is O(1) per window. This is Jellyfish's rule [5]: `mer_iterator.hpp` resets `filled_` to 0 whenever `mer_dna::code(c)` is negative, and `mer_dna.hpp` `codes[256]` maps only A/a, C/c, G/g, T/t to 0..3 (IUPAC codes, U, gaps and every other byte are negative). With `KmerCountingOptions.Canonical`, the forward ACGT-only table is folded once per distinct k-mer onto min(w, RC(w)) using the canonical `DnaSequence.GetReverseComplementString`. This is Jellyfish `count -C` [4][5]: `mer_iterator` returns `m_ < rcm_ ? m_ : rcm_`, and `mer_dna::get_canonical` does the same. Jellyfish compares 2-bit codes with A<C<G<T, the same as ordinal comparison of upper-case ACGT strings. Canonical counting is defined only over ACGT, because Jellyfish cannot encode any other base, so `Canonical = true` always applies the ACGT-only rule. The default options reproduce the literal overloads exactly.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Overlapping sliding-window k-mer extraction.
- Exact count accumulation in a dictionary keyed by the observed k-mer.
- Case-insensitive normalization for string-based input.
- ACGT-only window rule (option) — Jellyfish `mer_iterator` [5].
- Canonical k-mer representation min(w, RC(w)) (option) — Jellyfish `count -C` [4][5]; numerically identical to Jellyfish 2.3.1 (§7.3).

**Intentionally simplified:**

- Both-strand counting sums forward and reverse-complement counts rather than canonicalizing each k-mer to a single representative key; **consequence:** forward and reverse-complement words remain separate dictionary entries unless they are identical strings. Canonical keys (one entry min(w, RC(w)) per k-mer pair, Jellyfish `count -C`) are available through `KmerCountingOptions.Canonical` in `CountKmers(sequence, k, options)` (B06 audit round 1, F10); see also the kPAL ACGT-only both-strand mode (F25).
- The option-less overloads and `CountKmersSpan` do not enforce a DNA alphabet (generic definition, also used for non-DNA text); **consequence:** the usual `4^k` bound applies only to DNA input. The Jellyfish behaviour (non-ACGT windows skipped) is available as `KmerCountingOptions.AcgtOnly`, and canonical keys as `KmerCountingOptions.Canonical`. The literal mode stays the default for backward compatibility.

**Not implemented:** none (canonical collapsing and ACGT-only counting are available as options since B06 audit round 1).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null string input | Returns an empty dictionary | Explicit source guard |
| `k <= 0` with non-empty string input or span input | Throws `ArgumentOutOfRangeException` | Invalid k-mer length after the relevant input guard |
| Empty or null string input with `k <= 0` | Returns an empty dictionary | The string-based overload short-circuits before validating `k` |
| `k > sequence.Length` | Returns an empty dictionary | No valid windows exist |
| Homopolymer such as `AAAA`, `k = 2` | Returns one key with count `L - k + 1` | Every window is identical |
| String input with ambiguous or non-ACGT symbols | Counts those symbols literally after uppercasing | The string/span counting logic does not filter the alphabet |
| Same, with `AcgtOnly` or `Canonical` | Windows containing such a symbol are skipped (`ACGTNACGT`, k=4 → ACGT:2) | Jellyfish `mer_iterator` window reset [5] |
| `U` with `AcgtOnly`/`Canonical` | Treated as non-ACGT (`AAUUAA`, k=2 → AA:2) | Jellyfish `codes['U']` is negative |
| Reverse-complement palindrome with `Canonical` | Counted once per occurrence under its own key | min(w, RC(w)) = w |

### 6.2 Limitations

The current implementation uses string keys for observed k-mers and does not canonicalize reverse complements automatically. It is exact and general-purpose, but large genomes or large `k` values can still create substantial runtime and dictionary pressure because each window materializes and hashes a length-`k` key. In the default (literal) mode ambiguous symbols are retained; use `KmerCountingOptions` for the Jellyfish conventions.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

- Genome assembly via de Bruijn graph construction.
- Metagenomics binning based on k-mer signatures.
- Alignment-free sequence comparison.
- Error detection through k-mer spectra.
- Repeat analysis in repetitive regions.

### 7.3 Reference cross-check — Jellyfish 2.3.1 (B06 audit round 1)

Reference: the real Jellyfish 2.3.1 binary (Ubuntu package `jellyfish 2.3.1-3build1`), run as `jellyfish count -m k -s 10000 -t 1 [-C]` followed by `dump -c`, `stats` and `histo`. A Python replica of `mer_iterator` gave the same `dump` on all 20 rows (10 inputs × with/without `-C`). C# `CountKmers(seq, k, new KmerCountingOptions(Canonical: C, AcgtOnly: true))` equals every row (tests `KmerAnalyzer_CountingOptions_Tests`).

| Input | k | -C | Unique | Distinct | Total | Max | Notes |
|---|---|---|---|---|---|---|---|
| GAATTCACGTTGCAGGATCCATGC | 3 | yes | 7 | 14 | 22 | 3 | GCA:3 |
| GAATTCACGTTGCAGGATCCATGC | 4 | yes | 17 | 19 | 21 | 2 | ATCC:2, ATTC:2 |
| ACGTTGCATGTCGCATGATGCATGAGAGCT (BA1B) | 4 | no / yes | 17 / 16 | 21 / 20 | 27 | 3 / 4 | -C: ATGC:4 (= ATGC + GCAT 3) |
| acgtNNacgtacgRtTTGCAnA | 3 | no / yes | 6 / 2 | 8 / 5 | 11 | 3 / 5 | -C: AAA:1 ACG:5 CAA:1 GCA:2 GTA:2 |
| ACGTNACGT | 4 | either | 0 | 1 | 2 | 2 | ACGT:2 |
| AAUUAA | 2 | either | 0 | 1 | 2 | 2 | U resets |
| ATGATG | 3 | yes | 2 | 3 | 4 | 2 | ATC:1 ATG:2 TCA:1 |
| Rosalind KMER sample | 4 | no / yes | 91 / 23 | 209 / 130 | 412 | 8 / 10 | -C histo 1:23 2:34 3:27 4:25 5:7 6:5 7:4 9:3 10:2 |
| Rosalind KMER sample | 5 | no / yes | 292 / 181 | 348 / 279 | 411 | 4 / 6 | -C histo 1:181 2:74 3:18 4:3 5:2 6:1 |

Means and entropies on the `-C` tables (scipy 1.17.1 `entropy(counts, base=2)`): Rosalind k=4 mean 3.169230769230769, H 6.779144227048732; BA1B k=4 mean 1.35, H 4.1343361131944505; mixed k=3 mean 2.2, H 2.0403733936884962.

## 8. References

1. Wikipedia. "K-mer." https://en.wikipedia.org/wiki/K-mer
2. Rosalind. "K-mer Composition." https://rosalind.info/problems/kmer/
3. Compeau, P.E.C., Pevzner, P.A., Tesler, G. (2011). "How to apply de Bruijn graphs to genome assembly." Nature Biotechnology, 29(11), 987–991.
4. Marçais, G., Kingsford, C. (2011). "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers." Bioinformatics, 27(6), 764–770.
5. Jellyfish source (gmarcais/Jellyfish, master): `include/jellyfish/mer_iterator.hpp`, `include/jellyfish/mer_dna.hpp` (`codes[256]`, `get_canonical`), `sub_commands/count_main_cmdline.yaggo` (`-C, --canonical` "Count both strand, canonical representation"). https://github.com/gmarcais/Jellyfish
