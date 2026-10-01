# Both-Strand K-mer Counting

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer |
| Test Unit ID | KMER-BOTH-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Both-strand k-mer counting reports, for every k-mer observed in a DNA sequence, the number of times it occurs across **both** strands of the double-stranded molecule: its overlapping occurrences on the forward strand plus its overlapping occurrences on the reverse-complement strand. Because double-stranded DNA carries the same information on the complementary strand read in the opposite direction, a strand-aware profile sums each k-mer's count with that of its reverse complement [1][2]. This yields a strand-symmetric frequency table and is the additive ("balance") variant of strand handling [1], as opposed to the canonical collapsing used by tools such as Jellyfish and Mash [4][5]. The computation is exact and deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

DNA is double-stranded; a read or k-mer can originate from either strand, and the reverse-complement strand read 5'→3' presents the reverse complement of the forward sequence [6]. Strand-aware analyses therefore relate a k-mer `w` to its reverse complement `RC(w)`. Two conventions exist: (a) **canonical collapsing**, keeping only the lexicographically smaller of `w` and `RC(w)` as a single key [4][5]; (b) **additive balancing**, keeping every k-mer but summing complementary contributions [1][2]. This unit implements convention (b).

### 2.2 Core Model

For a sequence S of length L and k-mer length k, single-strand counting yields `forward[w]` = number of overlapping positions i (0 ≤ i ≤ L−k) where S[i..i+k) = w; there are L − k + 1 such positions [4]. The reverse-complement strand read 5'→3' is `RC(S)`, and by inversion symmetry the count of `w` on that strand equals the count of `RC(w)` on the forward strand [3]. Hence the both-strand count is:

> count[w] = forward[w] + forward[RC(w)]   [1][3]

equivalently obtained by counting k-mers in S and in RC(S) and summing per key. kPAL describes this exactly as balancing a profile "by adding the values of each k-mer to its reverse complement" [1][2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | count[w] = forward[w] + forward[RC(w)] | inversion symmetry: RC-strand count of w = forward count of RC(w) [3] |
| INV-02 | Σ_w count[w] = 2·(L − k + 1) for L ≥ k, else 0 | each strand contributes L − k + 1 windows [4] |
| INV-03 | count[w] = count[RC(w)] (strand-symmetric profile) | the sum is symmetric under w ↔ RC(w) [1][3] |
| INV-04 | palindromic w (RC(w) = w) ⇒ count[w] = 2·forward[w] | both contributions land on one key [3] |
| INV-05 | every count is a positive integer | counts are occurrence tallies [4] |

### 2.5 Comparison with Related Methods

| Aspect | Both-strand additive (this) | Canonical collapsing (Jellyfish -C / Mash) |
|--------|------------------------------|---------------------------------------------|
| Keys | every observed k-mer | only lexicographically smaller of {w, RC(w)} [4][5] |
| Count of w | forward[w] + forward[RC(w)] | occurrences of canonical(w) over both strands |
| w vs RC(w) | both present, equal counts | merged into one key |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string / DnaSequence | required | DNA sequence | case-insensitive (upper-cased internally); IUPAC bases complemented per repository `GetComplementBase` |
| k | int | required | k-mer length | k > 0 |
| options | KmerCountingOptions | literal | `AcgtOnly`: skip windows containing a non-ACGT symbol (kPAL / Jellyfish rule); `Canonical` is rejected | `Canonical = false` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | Dictionary<string,int> | each observed k-mer → summed forward + reverse-complement occurrence count |

### 3.3 Preconditions and Validation

- Null or empty sequence ⇒ empty dictionary.
- k > L (so L − k + 1 ≤ 0) ⇒ empty dictionary.
- k ≤ 0 ⇒ `ArgumentOutOfRangeException` (inherited from `CountKmers`).
- `options.Canonical = true` ⇒ `ArgumentException` (ParamName `options`): see §5.3.
- Input is upper-cased internally (case-insensitive); the reverse complement of recognized bases is uppercase.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count overlapping k-mers of the forward sequence (`CountKmers(sequence, k, options)`; literal or ACGT-only windows).
2. Copy the forward table; for every forward entry (w, c) add c to the entry of RC(w) (kPAL `balance()`).
   Windows of RC(S) are exactly the reverse complements of the windows of S, so this equals counting `RC(S)` and
   summing per key; a palindrome receives its own count twice.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CountKmersBothStrands | O(n·k + d·k) | O(d·k) | n = L − k + 1 windows (one counting pass); d = number of distinct forward k-mers, each reverse-complemented once |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.CountKmersBothStrands(string, int)`: additive both-strand counting; = the options overload with literal options.
- `KmerAnalyzer.CountKmersBothStrands(string, int, KmerCountingOptions)`: the implementation (literal or ACGT-only; B06 audit round 2, WP8).
- `KmerAnalyzer.CountKmersBothStrands(DnaSequence, int)`: delegates to the string overload.

### 5.2 Current Behavior

Counts forward k-mers via `CountKmers(sequence, k, options)` and adds each entry's count to the entry of its reverse complement (`DnaSequence.GetReverseComplementString` of the k-mer), which equals counting k-mers of the reverse-complement string and summing per key. The reverse-complement helper handles IUPAC ambiguity codes and is case-insensitive. A **suffix tree was not used**: this is a single linear two-pass tally over the sequence and its reverse complement (count-all-windows, not occurrence-enumeration of a query pattern), so the suffix-tree's O(m) post-construction query advantage does not apply; a direct O(n) scan is optimal here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- count[w] = forward[w] + forward[RC(w)] — kPAL "balance" (sum of each k-mer and its reverse complement) [1][2], grounded by inversion symmetry [3].
- Grand total 2·(L − k + 1) [4]; strand-symmetric profile [1]; palindrome doubling [3].
- Reference-confirmed (review 2026-09): kPAL's own `Profile.balance()` (`kpal/klib.py`, LUMC/kPAL master: `if i < i_rc: counts[i] += counts[i_rc]; counts[i_rc] += temp` / `elif i == i_rc: counts[i] += counts[i]`) was run on `ATGGC`/k=2, `ACGT`/k=2, `AAA`/k=2, `ATGC`/k=4 and `GAATTCACGTTGCAGGATCCATGC`/k=3,4,6; every balanced profile equals this method's output (and `Counter(S) + Counter(Bio.Seq.reverse_complement(S))`). Palindromes are doubled by kPAL exactly as here.
- Relation to canonical counting (Jellyfish `-C`) [4]: for a non-palindromic canonical k-mer c, count[c] = count[RC(c)] = canonical_count[c]; for a reverse-complement palindrome, count[w] = 2·canonical_count[w] (e.g. k=4 on the reference sequence above: AATT, ACGT, CATG, GATC, TGCA → 2 here, 1 under `-C`). Palindromes exist only for even k.

- ACGT-only counting (`KmerCountingOptions.AcgtOnly`, B06 audit round 2, WP8) is kPAL's counting rule: `Profile.from_sequences` (`kpal/klib.py`) splits every sequence on the regular expression `[^AaCcGgTt]` and counts only the k-mers inside each part, then `balance()` is applied. kPAL (LUMC/kPAL master, `klib.py` + `metrics.py`, run from source) gives exactly this method's ACGT-only output on `GAATTCNNACGTTGCAGGATCCATGCRYacgtgcaNTTGCA` k = 3 (26 keys, Σ 56, GCA = TGC = 7) and k = 4 (31 keys, Σ 48, TGCA = 6, AATT = 2), `ACGTNACGTAAcgtRTT` k = 3 (ACG = CGT = 6, Σ 18) and `AAAANTTTTGGGGuCCCC` k = 2 (AA = CC = GG = TT = 6, CA = TG = 1). The kPAL forward profile also equals Jellyfish 2.3.1 `count` (no `-C`) + `dump -c` on every one of these inputs (`u` is a non-ACGT symbol in both tools).
- Canonical option (decision, WP8): rejected with `ArgumentException`. The canonical count already is the both-strand count — Jellyfish `count -C` keys forward[w] + forward[RC(w)] by min(w, RC(w)) and counts a palindrome once (`KmerAnalyzer.CountKmers(sequence, k, Canonical)`). "Both strands + canonical" would be this balanced table restricted to canonical keys, which differs from `-C` exactly on palindromes (k = 4 on the sequence above: AATT 2 vs 1, ACGT 4 vs 2) — not 2× the canonical counts for non-palindromes, and a combination that neither kPAL (no canonical profiles) nor Jellyfish defines.

**Intentionally simplified:**

- (none)

**Not implemented:** none in this unit. Canonical collapsing (the lexicographically smaller of w and RC(w) as a single key, Jellyfish `-C` [4]) is provided by `KmerAnalyzer.CountKmers(sequence, k, new KmerCountingOptions(Canonical: true))` (KMER-COUNT-001, B06 audit round 1; equal to Jellyfish 2.3.1 `count -C`, K-mer_Counting §7.3). For ACGT input, canonical_count[c] = count[c] here for a non-palindromic c, and count[w]/2 for a palindrome w; canonical mode also skips non-ACGT windows, whereas this method keeps them (below).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / null sequence | empty dictionary | no windows (L − k + 1 ≤ 0) [4] |
| k > L | empty dictionary | L − k + 1 ≤ 0 [4] |
| k = L | one window per strand | boundary of the window formula [4] |
| Palindromic k-mer (RC(w)=w) | count doubled on one key | INV-04 [3] |
| k ≤ 0 | ArgumentOutOfRangeException | API contract (sibling `CountKmers`) |
| `AcgtOnly`, sequence with N/IUPAC/U | windows containing a non-ACGT symbol skipped on both strands; Σ = 2 × all-ACGT windows | kPAL `from_sequences` split rule |
| `Canonical = true` | ArgumentException | §5.3 decision (use `CountKmers` with `Canonical`) |

### 6.2 Limitations

Additive counting double-counts every k-mer's information relative to a single canonical key; it is not interchangeable with canonical k-mer sets used by sketching tools. Ambiguity codes (N, R, Y, …) are kept as literal k-mer keys (inherited from `CountKmers`, KMER-COUNT-001) and complemented by the canonical IUPAC `SequenceExtensions.GetComplementBase` (N↔N, R↔Y, …), so e.g. `AAN`/k=2 → {AA, AN, NT, TT}; kPAL and Jellyfish instead skip every k-mer containing a non-ACGT base, which is the `AcgtOnly` option (`AAN`/k=2 → {AA: 1, TT: 1}). The MCP tool accepts non-ACGT input only with `acgtOnly: true`. Non-IUPAC characters (gaps etc.) pass through the reverse-complement helper unchanged. `U` is not DNA: the forward strand keeps `U` while its complement is `A`, so strand symmetry (INV-03) does not hold for U-containing input (the MCP tool rejects non-DNA input via `RequireDna`).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var counts = KmerAnalyzer.CountKmersBothStrands("ATGGC", 2);
// counts = { AT:2, TG:1, GG:1, GC:2, CC:1, CA:1 }
```

**Numerical walk-through:** S = ATGGC (L=5), k=2. Forward 2-mers {AT,TG,GG,GC}. RC(ATGGC)=GCCAT, RC-strand 2-mers {GC,CC,CA,AT}. Summed: AT 1+1=2, TG 1, GG 1, GC 1+1=2, CC 1, CA 1. Grand total 8 = 2·(5−2+1). Check INV-01 for TG: forward[TG]+forward[RC(TG)=CA] = 1+0 = 1. ✓

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_CountKmersBothStrands_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmersBothStrands_Tests.cs) — covers INV-01..INV-05; [KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests.cs) — kPAL ACGT-only values, canonical rejection
- Evidence: [KMER-BOTH-001-Evidence.md](../../../docs/Evidence/KMER-BOTH-001-Evidence.md)
- Related algorithms: [K-mer_Generation](./K-mer_Generation.md), [K-mer_Counting](../K-mer/K-mer_Counting.md)

## 8. References

1. Anvar SY, et al. 2014. Determining the quality and complexity of next-generation sequencing data without a reference genome. Genome Biology, 15:555. https://doi.org/10.1186/s13059-014-0555-3
2. kPAL documentation — Methodology. https://kpal.readthedocs.io/en/latest/method.html
3. Shporer S, Chor B, Rosset S, Horn D. 2016. Inversion symmetry of DNA k-mer counts: validity and deviations. BMC Genomics. https://pmc.ncbi.nlm.nih.gov/articles/PMC5006273/
4. Marçais G, Kingsford C. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27(6):764–770. https://doi.org/10.1093/bioinformatics/btr011
5. Ondov BD, et al. Mash — canonical k-mer explanation, GitHub issue #45. https://github.com/marbl/Mash/issues/45
6. Clavijo BJ. 2018. BioInfoLogics — k-mer counting, part I: Introduction. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
