---
type: concept
title: "Both-strand (additive / kPAL-balance) k-mer counting"
tags: [analysis, algorithm]
mcp_tools:
  - count_kmers_both_strands
sources:
  - docs/Evidence/KMER-BOTH-001-Evidence.md
  - docs/algorithms/K-mer/Both_Strand_Kmer_Counting.md
source_commit: c85acb02f865639aee5a5c5cc2c2257a1295fe7c
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-both-001-evidence
      evidence: "Test Unit ID: KMER-BOTH-001; Algorithm: K-mer counting over both strands (docs/algorithms/K-mer/Both_Strand_Kmer_Counting.md)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:asynchronous-kmer-counting
      source: kmer-both-001-evidence
      evidence: "Both operations of the K-mer family are built on the same synchronous CountKmers sliding-window primitive (KmerAnalyzer)"
      confidence: high
      status: current
---

# Both-strand (additive / kPAL-balance) k-mer counting

The **K-mer** family's strand-aware count: `KmerAnalyzer.CountKmersBothStrands` returns a
`Dictionary<string,int>` giving, for every k-mer observed in a DNA sequence, its occurrence
count summed across **both** strands of the double-stranded molecule — its forward-strand
occurrences **plus** its reverse-complement-strand occurrences. This is the **additive
("balance")** strand convention (kPAL / Anvar et al. 2014), *not* the canonical-collapsing
convention of Jellyfish `-C` or Mash. Validated under test unit **KMER-BOTH-001** (record
[[kmer-both-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Core model

For a sequence `S` of length `L` and k-mer length `k` (`k > 0`), single-strand counting yields
`forward[w]` = the number of overlapping start positions `i` (`0 ≤ i ≤ L−k`) where `S[i..i+k) = w`
— there are **`L − k + 1`** such windows. The reverse-complement strand read 5'→3' is `RC(S)`, and
by **inversion symmetry** (Shporer et al. 2016) the count of `w` on that strand equals the count of
`RC(w)` on the forward strand. Hence:

> **count[w] = forward[w] + forward[RC(w)]**

equivalently: count k-mers in `S` and in `RC(S)`, then sum per key. kPAL states this exactly as
balancing a profile "by adding the values of each k-mer to its reverse complement" (Anvar et al.
2014) — it enforces strand symmetry, the both-strand view of double-stranded DNA.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `count[w] = forward[w] + forward[RC(w)]` | inversion symmetry: RC-strand count of `w` = forward count of `RC(w)` |
| INV-02 | `Σ count[w] = 2·(L − k + 1)` for `L ≥ k`, else 0 | each strand contributes `L − k + 1` windows |
| INV-03 | `count[w] = count[RC(w)]` (strand-symmetric profile) | the sum is symmetric under `w ↔ RC(w)` |
| INV-04 | palindromic `w` (`RC(w) = w`, e.g. `AT`, `GC`, `ACGT`) ⇒ `count[w] = 2·forward[w]` | both contributions land on the **same** key |
| INV-05 | every count is a positive integer | counts are occurrence tallies |

## Additive vs canonical (the two strand-aware conventions)

Both relate a k-mer `w` to `RC(w)`; they differ in what they emit. This unit implements the
**additive** convention only — canonical collapsing is explicitly **not implemented** (a future
canonical unit or external Jellyfish `-C` / Mash is the route for that view).

| Aspect | Additive both-strand (this) | Canonical collapsing (Jellyfish `-C` / Mash) |
|--------|------------------------------|-----------------------------------------------|
| Keys | every observed k-mer | only the lexicographically smaller of `{w, RC(w)}` |
| Count of `w` | `forward[w] + forward[RC(w)]` | occurrences of `canonical(w)` over both strands |
| `w` vs `RC(w)` | both present, equal counts | merged into one key |

## Contract and edge cases

- **Inputs:** `sequence` (`string` or `DnaSequence`; null/empty → empty result; upper-cased
  internally so counting is **case-insensitive**; IUPAC bases complemented via the repository
  `GetReverseComplementString` / `GetComplementBase` helper), `k` (`int`, `> 0`).
- **Output:** `Dictionary<string,int>` mapping each observed k-mer to its summed forward +
  reverse-complement count. The `DnaSequence` overload delegates to the `string` overload.
- **Edge cases:** empty/null sequence → empty dictionary; `k > L` (`L − k + 1 ≤ 0`) → empty
  dictionary; `k = L` → one window per strand; `k ≤ 0` → `ArgumentOutOfRangeException` (inherited
  from `CountKmers`). Non-IUPAC characters pass through the reverse-complement helper unchanged.

## Worked oracles

- `CountKmersBothStrands("ATGGC", 2)` → `{ AT:2, TG:1, GG:1, GC:2, CC:1, CA:1 }`. Forward
  `{AT,TG,GG,GC}`; `RC(ATGGC)=GCCAT` → `{GC,CC,CA,AT}`; summed. Grand total `8 = 2·(5−2+1)`.
  INV-01 check for `TG`: `forward[TG] + forward[RC(TG)=CA] = 1 + 0 = 1` ✓.
- `CountKmersBothStrands("ACGT", 2)` → `{ AC:2, CG:2, GT:2 }` — every 2-mer here is a
  reverse-complement palindrome (`RC(ACGT)=ACGT`), so each key is doubled (INV-04).
- `CountKmersBothStrands("AAA", 2)` → `{ AA:2, TT:2 }` — non-palindromic (`RC(AAA)=TTT`); the two
  keys carry equal counts (INV-03).

## Complexity and implementation

`O(n·k)` time (`n = L − k + 1` windows per strand, **two** linear passes — over `S` and `RC(S)`),
`O(d·k)` space (`d` = distinct k-mers). A **suffix tree was evaluated and not used**: full-spectrum
both-strand counting is a linear count-all-windows tally, not a fixed-pattern occurrence query, so a
suffix tree adds `O(n)` construction with no benefit.

## Deviations and assumptions

- **No source contradictions; deviations = None.** kPAL "balance" (Anvar et al. 2014) and inversion
  symmetry (Shporer et al. 2016) give the identical `forward[w] + forward[RC(w)]` semantics, and the
  Jellyfish/Mash canonical wording is cited only to *contrast* the not-implemented collapsing mode.
- **ASSUMPTION (API-shape, non-correctness-affecting):** empty/short input (`k > L`) → empty
  dictionary, and `k ≤ 0` → `ArgumentOutOfRangeException` — resolved to match the sibling
  `CountKmers` / `GenerateAllKmers` contract in this repository (sources define `k` as a positive
  substring length but prescribe no exception type or empty-boundary result).

## Relation to other units

- Built on the same synchronous sliding-window primitive `KmerAnalyzer.CountKmers` (KMER-COUNT-001,
  not yet ingested) that also underlies [[asynchronous-kmer-counting]] — this operation calls it
  once on `S` and once on `RC(S)`. The `L − k + 1` count model is shared; a future sync-count concept
  will own that definition and this page will link to it.
- The same k-mer machinery feeds the K-mer family's other operations (frequency analysis,
  [[k-mer-generation|generation]], statistics) and downstream [[de-bruijn-graph-assembly]],
  [[kmer-spectrum-error-correction]], and [[kmer-jaccard-similarity]]. Note [[k-mer-generation]]
  enumerates the full `n^k` universe (sequence-independent) rather than counting observed windows,
  and [[k-mer-positions]] indexes *where* one k-mer occurs (an ordered position list) rather than
  tallying *how many* — the inverse index to this counting operation, sharing the same
  case-insensitive upper-casing convention.
