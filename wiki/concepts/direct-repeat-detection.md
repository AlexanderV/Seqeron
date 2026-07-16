---
type: concept
title: "Direct repeat detection (same-orientation dispersed pairs)"
tags: [annotation, algorithm]
mcp_tools:
  - find_direct_repeats
sources:
  - docs/algorithms/Repeat_Analysis/Direct_Repeat_Detection.md
source_commit: 6b489eae876ae047fe066f85843b56e5359e969f
created: 2026-07-16
updated: 2026-07-16
---

# Direct repeat detection (same-orientation dispersed pairs)

Finding two or more identical sequence copies that recur in the **same 5'→3'
orientation** at different genomic positions — the *same-orientation* member of the
Seqeron `RepeatFinder` repeat family (test unit **REP-DIRECT-001**,
`RepeatFinder.FindDirectRepeats` in
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs`). Unlike an
inverted repeat, the downstream copy preserves the original sequence rather than its
reverse complement; an intervening **spacer** may be zero or many bases long
(`5' TTACG——TTACG 3'`).

This is a *distinct* operation from the two neighbouring repeat concepts and should be
linked, not re-derived:

- [[repetitive-element-detection]] — the repeats/tandem family anchor covers *tandem*
  (head-to-tail, adjacent) copies, *inverted* repeats (reverse-complement arms), and
  RepeatMasker-class assignment. Direct repeats are same-orientation **dispersed**
  pairs with configurable spacing, not any of those.
- [[longest-repeated-substring]] — enumerates **any** recurring substring within one
  sequence (`GenomicAnalyzer.FindLongestRepeat`/`FindRepeats`) and reports a position
  *list* per repeat. Direct-repeat detection instead reports **pairs**
  `(FirstPosition, SecondPosition)` filtered by an explicit spacing constraint.

See [[test-unit-registry]] for how the unit is tracked and
[[algorithm-validation-evidence]] for the artifact pattern.

## Core model

For a sequence `S`, repeat length `L`, first position `i`, and second position `j`, a
direct-repeat pair is present when the two windows are identical **and** the second
copy clears the required gap:

```text
S[i..i+L) = S[j..j+L)   and   j > i + L − 1 + minSpacing
```

Reported spacing is `Spacing = j − i − L` (nucleotides strictly between the two copies).
With `minSpacing = 0` the filter permits `j = i + L`, i.e. two copies that abut exactly
(the adjacent, tandem-like case); `minSpacing > 0` guarantees the copies do not overlap.

## API contract

`FindDirectRepeats` has two overloads over the same exact-copy model:

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `sequence` | `DnaSequence` or `string` | required | `DnaSequence` overload throws on `null`; `string` overload uppercases (`ToUpperInvariant`) and yields nothing for `null`/empty. |
| `minLength` | `int` | `5` | Rejected below `2`. |
| `maxLength` | `int` | `50` | Rejected below `minLength`. |
| `minSpacing` | `int` | `1` | `0` allowed → enables abutting copies. |

Each result (`DirectRepeatResult`) carries `FirstPosition`, `SecondPosition` (both
**0-based**), `RepeatSequence`, `Length`, and `Spacing`.

**Validation is now symmetric across both overloads.** The `DnaSequence` overload
throws `ArgumentNullException` for a null sequence and `ArgumentOutOfRangeException`
when `minLength < 2` or `maxLength < minLength`. The raw-string overload was a
[[research-grade-limitations|documented deviation]] (REP-DIRECT-001 fuzzing): it did
**not** mirror those range checks, so a degenerate `minLength = 0` produced a
zero-length candidate whose suffix-tree lookup matched every position and blew the
result set up to `O(n²)` spurious empty/single-base "repeats" — reachable through the
MCP `find_direct_repeats` tool that forwards raw user input. The fix hoists the
`minLength < 2` / `maxLength < minLength` guards into an eager wrapper on the raw-string
overload so the exception surfaces at the call site; both overloads now validate
identically.

## Algorithm and complexity

1. Normalize case (raw-string overload only).
2. Build a [[suffix-tree]] over the full sequence once (`SuffixTree.Build`, `O(n)`).
3. For each candidate length `minLength..maxLength`, enumerate start positions with room
   for two copies plus the spacing.
4. Extract the candidate at the start and query `FindAllOccurrences` for every
   occurrence.
5. Keep each downstream occurrence satisfying the spacing filter, suppress duplicate
   `(FirstPosition, SecondPosition, Length)` tuples with a hash set, and emit a result
   per qualifying pair.

Suffix-tree construction is `O(n)` time/space; detection is
`O(r × n × (m + k))` where `r` = number of tested lengths, `m` = candidate length, and
`k` = occurrences returned per candidate.

## Invariants (test oracles)

- `RepeatSequence` is identical at `FirstPosition` and `SecondPosition` (results come
  only from suffix-tree occurrence matches of the same extracted pattern).
- `Spacing == SecondPosition − FirstPosition − Length` (constructed from the stored
  coordinates).
- With `minSpacing > 0`, reported copies do not overlap (`j > i + len − 1 + minSpacing`).
- Each `(FirstPosition, SecondPosition, Length)` tuple is unique (hash-set dedup).
- A single start position yields **multiple** result pairs when the same repeat occurs
  at several later positions that all clear the spacing filter.

## Edge cases and scope

Empty sequence, a sequence too short for two copies plus `minSpacing`, or no qualifying
downstream occurrence → empty enumerable. Adjacent identical copies are reported only
when `minSpacing = 0`.

**Exact matching only** — no mismatch, gap, or indel tolerance, so approximate,
degenerate, or interrupted direct repeats are not reported. Output is **raw repeat-pair
coordinates**, with no higher-level biological annotation (LTR/transposon boundaries,
replication-slippage deletion substrates, trinucleotide-expansion loci such as
`CAG`/`CGG`/`GAA`/`CTG`, or regulatory context); interpreting *why* a reported repeat is
present is the caller's responsibility. A documented Framework/Simplified
[[research-grade-limitations|limitation]], not an invented constraint.
