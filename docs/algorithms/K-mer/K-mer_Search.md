# K-mer Search Algorithms

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-FIND-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-09-28 |

## 1. Overview

K-mer search algorithms identify k-mers of special interest within a sequence rather than returning the full count map. In this repository, the documented search surface includes the most frequent k-mers, unique k-mers, and `(L, t)` clumps. All three operations are built on exact k-mer counts, with clump finding using a sliding-window update strategy and a deduplicating result set.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A most frequent k-mer is any k-mer attaining the maximum count in a sequence, and multiple k-mers may tie for that maximum. A unique k-mer occurs exactly once. A pattern forms an `(L, t)` clump if some window of length `L` contains at least `t` occurrences of that pattern. The original document also notes that clumps can indicate biologically interesting regions such as origins of replication or other motif-rich segments. Sources: Rosalind BA1B, Rosalind BA1E, Wikipedia (K-mer).

### 2.2 Core Model

`FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` are filters over exact k-mer counts:

$$
MostFrequent = \{kmer : Count(kmer) = \max_j Count(j)\}
$$

$$
Unique = \{kmer : Count(kmer) = 1\}
$$

For clumps, the windows are the substrings `Genome[i..i+L-1]` for `i ∈ [0, |Genome| − L]`; an occurrence counts only when it lies entirely inside the window (start `p` with `i ≤ p ≤ i + L − k`), so each window holds `L − k + 1` k-mer starts and overlapping occurrences count (Rosalind BA1E; Compeau & Pevzner ch. 1). The implementation is the textbook's `BetterClumpFinding`: count the first window with `CountKmers(...)`, then per slide decrement the leaving k-mer and increment the entering one; since only the entering k-mer's count can grow, only it is tested against `t`. `FindClumpWindows` additionally reports *where* each k-mer forms a clump: the set of window starts `{ i : Genome[i..i+L−1] contains ≥ t occurrences }` as maximal runs `[FirstWindowStart, LastWindowStart]` (a run covers `Genome[FirstWindowStart .. LastWindowStart + L − 1]`). In the same pass a run opens when the entering k-mer reaches `t` and closes when the leaving k-mer falls to `t − 1` (a slide whose leaving and entering k-mers are equal changes nothing), so the first run start is the leftmost qualifying window — the window in which `BetterClumpFinding` first detects the k-mer.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every k-mer returned by `FindUniqueKmers(...)` has count exactly `1` | The method filters the count dictionary by `kvp.Value == 1` |
| INV-02 | Every k-mer returned by `FindMostFrequentKmers(...)` has the maximum observed count | The method filters by `counts.Values.Max()` |
| INV-03 | `FindClumps(...)` returns each qualifying k-mer at most once | Results are stored in a `HashSet<string>` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Sequence to search | Empty or null input yields an empty result |
| `k` | `int` | required | K-mer length | `FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` inherit `CountKmers(...)` validation; `FindClumps(...)` returns empty when `k <= 0` |
| `windowSize` | `int` | required for clumps | Window length `L` for clump detection | `FindClumps(...)` returns empty when `windowSize < k` or `windowSize > sequence.Length` |
| `minOccurrences` | `int` | required for clumps | Minimum number of occurrences `t` in a window | `FindClumps(...)` returns empty when `minOccurrences <= 0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `mostFrequent` | `IEnumerable<string>` | All k-mers tied for the maximum observed count |
| `unique` | `IEnumerable<string>` | All k-mers observed exactly once |
| `clumps` | `IEnumerable<string>` | Unique set of k-mers that satisfy the `(L, t)` clump condition |

### 3.3 Preconditions and Validation

`FindMostFrequentKmers(...)` and `FindUniqueKmers(...)` use `CountKmers(...)`, so `k <= 0` raises `ArgumentOutOfRangeException` there. `FindClumps(...)` instead treats invalid parameters as empty-result conditions: null or empty sequence, `k <= 0`, `windowSize < k`, `windowSize > sequence.Length`, or `minOccurrences <= 0` all yield no clumps. All string-based searches uppercase the sequence internally.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input sequence to uppercase.
2. For most-frequent and unique searches, compute the full k-mer count map.
3. Return either the keys tied at the maximum count or the keys with count `1`.
4. For clumps, count the first window with the canonical `CountKmers(...)`, record k-mers meeting the threshold, then slide the window by one position, updating the leaving/entering k-mer counts and testing only the entering k-mer against `t`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindMostFrequentKmers` | `O(n·k)` | `O(u)` | Builds exact counts and filters maxima (k for substring hashing) |
| `FindUniqueKmers` | `O(n)` | `O(u)` | Builds exact counts and filters singletons |
| `FindClumps` | `O(n·k)` | `O(min(L, u))` | BetterClumpFinding: one decrement + one increment + one threshold test per slide (before 2026-09 every slide rescanned the whole window map, `O(n·(L−k+1))`; E. coli k=9,L=500,t=3 went from ~5.7 s to ~0.75 s) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindMostFrequentKmers(string, int)`: Returns all maxima from the count map.
- `KmerAnalyzer.FindMostFrequentKmers(string, int, KmerCountingOptions)`: the maxima of the literal / ACGT-only / canonical (`jellyfish count -C`) table (B06 audit round 2, WP8; §7.4).
- `KmerAnalyzer.FindUniqueKmers(string, int[, KmerCountingOptions])`: Returns all singleton k-mers (option-aware overload: `jellyfish count -C` + `dump -L 1 -U 1`, see Unique_And_MinCount_Kmers.md).
- `KmerAnalyzer.FindClumps(string, int, int, int)`: Returns deduplicated clump-forming k-mers (streamed).
- `KmerAnalyzer.FindClumpWindows(string, int, int, int)`: Returns `IReadOnlyList<KmerClump>` — each clump k-mer with its maximal runs of qualifying window starts (`ClumpWindowRun(FirstWindowStart, LastWindowStart)`, inclusive, 0-based; `KmerClump.FirstWindowStart` = leftmost qualifying window), ordered by first window then ordinal k-mer (B06 audit round 1 WP3, F12). Both methods consume one private sliding pass (`ScanClumpTransitions`), so the sliding logic exists once.
- MCP (Seqeron.Mcp.Analysis): `most_frequent_kmers`, `unique_kmers`, `find_clumps` (k-mer set) and `find_clump_windows` (k-mers with their window runs; B06 audit round 1 WP5), each delegating to the method above.

### 5.2 Current Behavior

All three methods uppercase the input sequence. All three reuse the canonical `CountKmers(...)` (`FindClumps(...)` for its first window) and then `FindClumps(...)` maintains the per-window count dictionary incrementally and a `HashSet<string>` of discovered clumps, streaming each clump k-mer once in order of first detection (first-window k-mers in ordinal order). `FindClumpWindows(...)` uses the same pass and also closes runs. Result order of every method is not part of the contract (most-frequent: first-occurrence order). `FindClumps(...)` returns empty rather than throwing on invalid window or threshold parameters.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Identification of all k-mers tied at the maximum count; with `Canonical` the arg-max of `jellyfish count -C` + `dump -c` (a k-mer and its reverse complement pooled — the exact-match reverse-complement-aware frequent-words question; mismatches, Rosalind BA1J, are `ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements` in the Alignment module).
- Identification of singleton k-mers.
- Sliding-window `(L, t)` clump detection.

**Intentionally simplified:**

- (none) — the qualifying windows are reported by `FindClumpWindows(...)` (audit round 1 WP3); `FindClumps(...)` keeps the set-only textbook output.

**Not implemented:**

- Per-window multiplicity traces (count of a k-mer in every window); they are not part of the clump definition, and the runs of `FindClumpWindows(...)` already give every window in which the count is ≥ t. Callers needing raw traces can combine `FindKmerPositions(...)` with the window rule above.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns an empty result | No valid windows exist |
| `k <= 0` | Throws for count-based searches; returns empty for `FindClumps(...)` | Different validation paths in source |
| `k > sequence.Length` | Returns an empty result | No valid k-mers exist |
| `windowSize > sequence.Length` | `FindClumps(...)` returns empty | No full window exists |
| `windowSize < k` | `FindClumps(...)` returns empty | A window cannot contain a full k-mer |
| All k-mers equally frequent | `FindMostFrequentKmers(...)` returns all observed k-mers | Every k-mer is tied at the maximum count |

### 6.2 Limitations

`FindClumps(...)` reports pattern-level results; `FindClumpWindows(...)` adds the qualifying window runs. As with the underlying count-based helpers, memory usage still depends on the number of unique k-mers maintained in dictionaries and sets.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

The original document cites `ACGTTGCATGTCGCATGATGCATGAGAGCT` as a standard frequent-word example in which the most frequent 4-mers are `CATG` and `GCAT`, each appearing 3 times. It also gives `TGCA` as a `(25, 3)` clump in:

```text
gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac
```

### 7.2 Applications and Use Cases (Optional)

- Unique k-mers for marker discovery and genomic fingerprinting.
- Clump finding for motif-rich regions such as origins of replication.

### 7.3 Reference cross-check (2026-09 review)

| Dataset | Parameters | Expected (source) | Seqeron |
|---|---|---|---|
| BA1B sample `ACGTTGCATGTCGCATGATGCATGAGAGCT` | k=4 | `CATG GCAT` (Rosalind BA1B) | identical |
| BA1E sample | k=5, L=75, t=4 | `CGACA GAAGA AATGT` (Rosalind BA1E) | identical |
| BA1E statement example | k=4, t=3, L=25/22/21 | `TGCA` / `TGCA` / none (3 occurrences span 22 bp) | identical |
| E. coli genome (textbook dataset, 4,639,675 bp) | k=9, L=500, t=3 | 1904 distinct 9-mers (textbook exercise answer) | 1904, set-equal to an independent Python implementation |
| 3000 random strings (alphabets 1–4, n ≤ 40) | random k, L, t | Python brute force over all windows | 0 mismatches |

**Clump windows (`FindClumpWindows`, audit round 1 WP3).** Reference: a Python brute force that counts every window `Genome[i..i+L−1]` and compresses the qualifying `i` of each k-mer into maximal runs; it agrees with an independent occurrence-interval method (window `i` qualifies iff some t consecutive occurrences `p_j … p_{j+t−1}` satisfy `p_{j+t−1} + k − L ≤ i ≤ p_j`, unioned) on all cases below.

| Dataset | k, L, t | Reference runs (k-mer:first-last window start) | Seqeron |
|---|---|---|---|
| BA1E sample | 5, 75, 4 | `CGACA:0-6 GAAGA:0-16 AATGT:16-21` (AATGT at 21, 73, 81, 86 → [86+5−75, 21]) | identical |
| BA1B sample | 4, 30, 3 | `CATG:0-0 GCAT:0-0` | identical |
| BA1B sample | 4, 12, 2 | `GCAT:4-5,11-12 CATG:5-6,12-13 ATGA:13-14` (split runs) | identical |
| `ACACGTTTTTTTTTTACACGTGGGGGGGGGGACACGT` | 4, 15, 2 | `TTTT:0-10 GGGG:11-22` | identical |
| `AAAAAAAAAA` | 2, 4, 3 | `AA:0-6` (leaving = entering on every slide) | identical |
| 3000 random strings (alphabets AC/ACGT/ACGTN/acgt, n ≤ 60, k ≤ 5, L ≤ 30, t ≤ 5; 743 with clumps) | random | brute force | 0 mismatches (runs and order); FindClumps set equal in all 3000 |
| 2 Mbp random ACGT (seed 7) | 6, 300, 4 | occurrence-interval method | 113 k-mers, 142 runs, identical (C# 0.38 s) |

### 7.4 Canonical / ACGT-only most frequent k-mers (B06 audit round 2, WP8)

Reference: **Jellyfish 2.3.1** `jellyfish count -m k -s 10000 [-C]` + `jellyfish dump -c`, arg-max taken over the dump.

| Input | k | `-C` (canonical) | no `-C` (ACGT-only) |
|---|---|---|---|
| BA1B sample `ACGTTGCATGTCGCATGATGCATGAGAGCT` | 4 | `ATGC` (4: ATGC + GCAT; CATG palindrome 3) | `CATG GCAT` (3) |
| `GAATTCNNACGTTGCAGGATCCATGCRYacgtgcaNTTGCA` | 2 | `CA` (8) | `CA GC TG` (4) |
| same | 3 | `GCA` (7) | `TGC` (4) |
| same | 4 | `TGCA` (3) | `TGCA` (3) |
| `AAAANTTTTGGGGuCCCC` | 2 | `AA CC` (6) | `AA CC GG TT` (3) |
| `ACGTNACGTAAcgtRTT` | 3 | `ACG` (6) | — |

`FindMostFrequentKmers(sequence, k, options)` is identical on every row (`KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests`;
MCP `most_frequent_kmers` optional `canonical` / `acgtOnly`).

The E. coli genome was not re-run for the window output (the textbook file is not in the repository and no GitHub mirror path could be resolved); its 1904-set is unchanged because `FindClumps` emits exactly the run openings of the shared pass (set equality with `FindClumpWindows` is tested).

## 8. References

1. Rosalind BA1B - Find the Most Frequent Words in a String. https://rosalind.info/problems/ba1b/
2. Rosalind BA1E - Find Patterns Forming Clumps in a String. https://rosalind.info/problems/ba1e/
3. Wikipedia (K-mer). https://en.wikipedia.org/wiki/K-mer
4. Compeau P., Pevzner P. *Bioinformatics Algorithms: An Active Learning Approach*, ch. 1 (FrequentWords / BetterFrequentWords, ClumpFinding / BetterClumpFinding).
