# Unique K-mers and K-mers with Minimum Count

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-UNIQUE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 |

## 1. Overview

These operations filter the k-mer frequency spectrum of a DNA/RNA sequence. `FindUniqueKmers` returns the **unique** k-mers — those appearing exactly once (occurrence count = 1) [2][4]. `FindKmersWithMinCount` returns **recurrent** k-mers — those whose overlapping occurrence count is at least a threshold `minCount` (Count(Text, Pattern) ≥ t) [3] — paired with their counts and ordered by count descending (ties by ordinal k-mer). An overload adds an inclusive upper bound `maxCount`; together the bounds are the `-L/--lower-count` and `-U/--upper-count` filters of `jellyfish dump` [4] (KMC `-ci`/`-cx` [5]). Both are exact, deterministic, combinatorial operations derived directly from the k-mer definition; neither is heuristic or probabilistic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a substring of length k contained within a biological sequence [1]. K-mers are extracted with a sliding window of step 1, so adjacent k-mers overlap by k−1 characters and a sequence of length L yields **L − k + 1** total k-mers [1]. The frequency spectrum distinguishes *total* k-mers (with duplicates), *distinct* k-mers (each different string once), and *unique* k-mers (frequency exactly 1) [2].

### 2.2 Core Model

For a sequence `Text` and k-mer `Pattern`, `Count(Text, Pattern)` is the number of overlapping occurrences of `Pattern` in `Text` [3]. Define the multiset of all length-k substrings and group by string to obtain counts c(P) = Count(Text, P).

- **Unique k-mers:** `{ P : c(P) = 1 }` — "Unique k-mers are those that appear only once" [2].
- **K-mers with minimum count:** `{ (P, c(P)) : c(P) ≥ minCount }`, the recurrent k-mers (Count ≥ t) [3].
- **K-mers in a count range:** `{ (P, c(P)) : minCount ≤ c(P) ≤ maxCount }` — `jellyfish dump` skips a k-mer when `val < lower_count || val > upper_count` (`sub_commands/dump_main.cc`; defaults lower 0, upper 2^64) [4]; KMC excludes k-mers occurring fewer than `-ci` or more than `-cx` times [5]. Unique = range [1, 1] (`jellyfish dump -L 1 -U 1`).

**Terminology (shared with KMER-STATS-001):** *unique* = count exactly 1 = Jellyfish `stats` "Unique" = `KmerStatistics.SingletonKmers`; *distinct* = every different k-mer once = Jellyfish "Distinct" = `KmerStatistics.DistinctKmers` (legacy `UniqueKmers` field). The range filter is one private helper (`SelectByCountRange`) shared by `FindKmersWithMinCount`, `FindUniqueKmers` and `AnalyzeKmers(seq, k, lowerCount, upperCount)`, so `FindKmersWithMinCount(seq,k,L,U)` lists exactly the k-mers that `AnalyzeKmers(seq,k,L,U)` summarises.

Worked example (ATCGATCAC, k=3) [2]: 7 total, 6 distinct, 5 unique = {TCG, CGA, GAT, TCA, CAC}; ATC occurs twice (c=2) so it is distinct but not unique.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each k-mer in `FindUniqueKmers` has c(P) = 1 | Filter predicate is `c(P) = 1` [2] |
| INV-02 | unique ⊆ distinct, so &#124;unique&#124; ≤ &#124;distinct&#124; ≤ L − k + 1 | Total count bound [1]; unique is a frequency-1 subset of distinct [2] |
| INV-03 | Each pair in `FindKmersWithMinCount` has c(P) ≥ minCount and the reported count equals the overlapping occurrence count | Filter predicate `c(P) ≥ t` over exact counts [3] |
| INV-04 | `FindKmersWithMinCount` output is ordered by count non-increasing, ties by ascending ordinal k-mer; `FindUniqueKmers` is in ascending ordinal order | Documented deterministic order (Jellyfish dump emits hash order [4]; ours = `dump -c \| sort -k2,2nr -k1,1`) |
| INV-06 | Range [L, U] output = k-mers summarised by `AnalyzeKmers(seq,k,L,U)` (count = Distinct, Σ count = Total, #count-1 = SingletonKmers) | Same predicate as Jellyfish `dump`/`stats` [4], one shared helper |
| INV-05 | With minCount ≤ 1 the returned keys equal the distinct k-mer set | Every observed k-mer has c(P) ≥ 1 [2][3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Sequence to analyze | null/empty → empty result; upper-cased internally (case-insensitive) |
| k | int | required | K-mer length | Must be > 0; k > L → empty result |
| minCount | int | required (`FindKmersWithMinCount`) | Inclusive minimum occurrence threshold | ≤ 1 (incl. ≤ 0) selects all distinct k-mers |
| maxCount | int | `int.MaxValue` (4-arg overload) | Inclusive maximum occurrence threshold (Jellyfish `-U`) | negative → `ArgumentOutOfRangeException`; `maxCount < minCount` → empty |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (FindUniqueKmers) | IEnumerable\<string\> | K-mers with occurrence count = 1, ascending ordinal order |
| (FindKmersWithMinCount) | IEnumerable\<(string Kmer, int Count)\> | K-mers in the count range with their counts, ordered by Count descending, then ordinal k-mer |

### 3.3 Preconditions and Validation

Input is upper-cased (T/U not normalised; standard string k-mers). Indexing is 0-based over the character string. Null or empty sequence returns an empty result; k > sequence length returns empty (L − k + 1 ≤ 0). k ≤ 0 throws `ArgumentOutOfRangeException` for a non-empty sequence (null/empty returns empty for any k, as in `CountKmers`). All exceptions are thrown at the call, not on enumeration.

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the k-mer count map via `CountKmers` (single O(n) pass, overlapping window, step 1) [1].
2. `FindKmersWithMinCount(seq, k, min, max)`: keep entries with min ≤ count ≤ max via the shared `SelectByCountRange` [4]; order by count descending, then ordinal k-mer. The 3-argument overload uses max = `int.MaxValue`.
3. `FindUniqueKmers`: `FindKmersWithMinCount(seq, k, 1, 1)` keys [2][4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindUniqueKmers | O(n·k + u log u) | O(d·k) | n = length, d = distinct, u = unique k-mers; substring extraction is O(k) per position; ordinal sort of the output |
| FindKmersWithMinCount | O(n·k + d log d) | O(d·k) | dominated by counting; sort over d distinct k-mers |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.FindUniqueKmers(string, int)`: returns k-mers with count = 1 (ordinal order).
- `KmerAnalyzer.FindKmersWithMinCount(string, int, int)`: returns (k-mer, count) with count ≥ minCount.
- `KmerAnalyzer.FindKmersWithMinCount(string, int, int, int)`: returns (k-mer, count) with minCount ≤ count ≤ maxCount (Jellyfish `dump -L/-U`).
- MCP: `unique_kmers` and `kmers_with_min_count` (Seqeron.Mcp.Analysis; optional `maxCount`) delegate.

### 5.2 Current Behavior

Both methods delegate counting to `KmerAnalyzer.CountKmers`, inheriting its null/empty/k-bounds handling and upper-casing. Before the 2026-10 review, both methods returned dictionary enumeration order (ties in `FindKmersWithMinCount`, all of `FindUniqueKmers`), which .NET does not specify; the order is now fixed as documented above. The repository **suffix tree** was evaluated and not used: both methods need the full frequency spectrum (count of every distinct k-mer), which a single linear hash-map pass over L − k + 1 windows yields in O(n·k); a suffix tree adds construction overhead without benefit for this exhaustive-count workload (it shines for many exact-match queries against one text, not for enumerating all k-mer frequencies).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Unique k-mers = occurrence count exactly 1 [2].
- Min-count filter = overlapping occurrence count ≥ t [3].
- Overlapping, step-1 counting (L − k + 1 total) [1].
- Count-descending ordering of recurrent k-mers (most-frequent-first) [3]; ties ordinal.
- Inclusive lower/upper count filters = `jellyfish dump -L/-U` [4], KMC `-ci/-cx` [5].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Canonical (reverse-complement-merged) k-mer counting; users should rely on `KmerAnalyzer.CountKmersBothStrands` for strand-aware counts.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty sequence | empty result | No k-mers exist |
| k > L | empty result | L − k + 1 ≤ 0 [1] |
| k ≤ 0 | ArgumentOutOfRangeException | k-mer length must be positive [1] |
| Homopolymer (AAAAA, k=3) | no unique k-mers | AAA has count 3 > 1 [2] |
| minCount ≤ 1 | all distinct k-mers | every k-mer has count ≥ 1 [2][3] |
| minCount > max count | empty | no k-mer meets the threshold [3] |
| maxCount < minCount | empty | Jellyfish `dump -L 3 -U 2` prints nothing [4] |
| maxCount < 0 | ArgumentOutOfRangeException | Jellyfish bounds are unsigned [4]; same rule as `AnalyzeKmers` upperCount |

### 6.2 Limitations

Counts forward-strand string k-mers only (no reverse-complement merging); does not normalise RNA/DNA (T vs U treated as different characters). Non-ACGT symbols are counted as literal k-mer characters (Jellyfish skips them; see K-mer_Counting §5.3).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var unique = KmerAnalyzer.FindUniqueKmers("ATCGATCAC", 3).ToHashSet();
// { "TCG", "CGA", "GAT", "TCA", "CAC" }  — ATC (count 2) excluded

var recurrent = KmerAnalyzer.FindKmersWithMinCount("ACGTACGT", 4, 2).ToList();
// [ ("ACGT", 2) ]

var band = KmerAnalyzer.FindKmersWithMinCount("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 2, 2).ToList();
// [ ("ATGA", 2), ("TGCA", 2) ]   — jellyfish dump -L 2 -U 2 (CATG, GCAT have count 3)
```

### 7.2 Reference cross-check (Jellyfish `dump` replica)

Python replica of `dump_main.cc` over `collections.Counter`, sorted by (−count, k-mer); C# equals every row (tests `FindKmersWithMinCount_Range_MatchesJellyfishDumpReplica`, `FindUniqueKmers_MatchesJellyfishDumpU1_InOrdinalOrder`).

| Sequence | k | -L | -U | Output |
|----------|---|----|----|--------|
| ACGTTGCATGTCGCATGATGCATGAGAGCT (BA1B sample) | 4 | 1 | 1 | 17 k-mers: ACGT AGAG AGCT ATGC ATGT CGCA CGTT GAGA GAGC GATG GTCG GTTG TCGC TGAG TGAT TGTC TTGC (= Jellyfish stats Unique 17, KMER-STATS-001) |
| same | 4 | 2 | ∞ | CATG:3 GCAT:3 ATGA:2 TGCA:2 |
| same | 4 | 2 | 2 | ATGA:2 TGCA:2 |
| same | 4 | 3 | 3 | CATG:3 GCAT:3 |
| GTAGAGCTGT | 2 | 1 | 1 | CT GA GC TA TG |
| GTAGAGCTGT | 2 | 2 | ∞ | AG:2 GT:2 |
| ATCGATCAC | 3 | 1 | 1 | CAC CGA GAT TCA TCG |
| AAAACGTAAA | 2 | 2 | ∞ / 2 | AA:5 / (none) |

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_FindUniqueAndMinCount_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_FindUniqueAndMinCount_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [KMER-UNIQUE-001-Evidence.md](../../../docs/Evidence/KMER-UNIQUE-001-Evidence.md)

## 8. References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo B, et al. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach* (2nd ed.). Active Learning Publishers. https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637
4. Marçais G, Kingsford C. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics* 27(6):764–770. Jellyfish source `sub_commands/dump_main.cc`, `dump_main_cmdline.yaggo`, `stats_main.cc` (https://github.com/gmarcais/Jellyfish, opened via raw.githubusercontent.com).
5. Kokot M, Długosz M, Deorowicz S. 2017. KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* 33(17):2759–2761. CLI usage `kmc_CLI/kmc.cpp` (https://github.com/refresh-bio/KMC): `-ci<value>` exclude k-mers occurring less than value times, `-cx<value>` exclude k-mers occurring more than value times.
