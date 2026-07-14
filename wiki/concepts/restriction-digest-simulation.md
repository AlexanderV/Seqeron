---
type: concept
title: "Restriction digest simulation (sequence partitioning into fragments + restriction map + end compatibility)"
tags: [moldesign, restriction]
mcp_tools:
  - find_all_restriction_sites
  - find_restriction_sites
  - get_enzyme
sources:
  - docs/algorithms/MolTools/Restriction_Digest_Simulation.md
source_commit: 0f959309abe7855c9429aeb0059af26eb694597e
created: 2026-07-14
updated: 2026-07-14
---

# Restriction digest simulation

The **RESTR** family member that models DNA fragmentation *after* cleavage — the
downstream complement of [[restriction-enzyme-filtering]] (which selects enzymes by
cut properties) and of restriction-site *detection* (which locates cut sites in a
sequence). Where filtering operates on enzyme metadata, digest simulation operates
on a **target sequence**: it partitions the sequence into fragments, builds a
restriction map, and tests whether two enzymes leave ligatable ends. Implemented in
`RestrictionAnalyzer` (`Digest`, `GetDigestSummary`, `CreateMap`, `AreCompatible`,
`FindCompatibleEnzymes`); status **Simplified**. This is a primary per-algorithm
spec (Test Unit ID N/A), synthesized against the
[[algorithm-validation-evidence|validation-evidence]] pattern used across MolTools.

## Core model — cut positions to fragments

Given a sequence of length `L` and forward-strand cut positions `c1 < c2 < … < ck`,
the digest fragments are the **half-open intervals** partitioning `[0, L)`:

```
[0, c1), [c1, c2), …, [ck, L)     →  k+1 fragments for k cuts
```

The algorithm: validate input → find restriction sites per requested enzyme → keep
**forward-strand cut positions only** (so a palindromic site — which reads the same on
both strands — is not double-counted) → sort + deduplicate cuts → emit fragments
between adjacent boundaries. `Digest(...)` inserts the boundary list
`[0, cuts…, L]`, yielding one `DigestFragment` per interval.

Each `DigestFragment` carries its subsequence, start, length, fragment number, and
**one representative** `LeftEnzyme`/`RightEnzyme` per boundary. When several enzymes
cut at the *same* coordinate the position is kept but boundary provenance collapses to
a single enzyme name (a documented simplification).

## Summary and map surfaces

- `GetDigestSummary(...)` → `DigestSummary`: total fragment count, fragment sizes
  **sorted descending**, largest/smallest, average size, enzyme list. `O(n + k log k)`.
- `CreateMap(...)` → `RestrictionMap`: sequence length, per-site records, positions
  grouped by enzyme, forward-strand `TotalSites`, `UniqueCutters` (enzymes with exactly
  one distinct forward-strand position), and `NonCutters` (requested enzymes absent from
  the grouped map). Passing **zero** enzyme names scans the full built-in catalog.

## End compatibility

`AreCompatible(enzyme1, enzyme2)` decides whether two enzymes leave ligatable ends
— the downstream reason a user filters by cutter class before a cloning/ligation step
(the compatibility fact [[restriction-enzyme-filtering]] mentions but does not compute):

| Condition | Result |
|-----------|--------|
| Both produce **blunt** ends | Compatible |
| Overhang **types** differ | Not compatible |
| Overhang types match **and** overhang sequences match | Compatible |
| Overhang sequences differ | Not compatible |

Worked examples from the spec: **BamHI + BglII** compatible (both leave a `GATC`
5' overhang), **EcoRV + SmaI** compatible (both blunt), **EcoRI + PstI** not
compatible. The relation is symmetric — `AreCompatible(A, B) == AreCompatible(B, A)`.
`FindCompatibleEnzymes()` enumerates all compatible pairs in the built-in catalog.
An unknown enzyme name yields `false` rather than throwing.

## Key invariants and contract

- **Fragment count:** `k` forward-strand cuts → `k + 1` fragments.
- **Fragment-sum principle:** the sum of fragment lengths **equals** the original
  sequence length (adjacent half-open boundaries partition the sequence). This is the
  spec's experimentally-checkable invariant (gel band sizes must sum to the input).
- **Descending sizes:** `DigestSummary.FragmentSizes` are sorted largest-first.
- **Compatibility symmetry:** `AreCompatible` is order-independent.
- **No-cut case:** when no site is found, `Digest` returns a **single fragment equal to
  the whole sequence** (explicit special case).
- **Boundary enzymes:** first fragment `LeftEnzyme = null`, last fragment
  `RightEnzyme = null` (no cut before position 0 / after the final boundary).
- **Zero-length fragments are never emitted** (yielded only when `length > 0`).
- **Validation:** null `sequence` → `ArgumentNullException`; `Digest` with no enzyme
  names → `ArgumentException`; `CreateMap` accepts zero names (scans all).

## Scope and limitations

A **sequence-partitioning simulation only** — it reports virtual fragments,
statistics, and maps. It does **not** simulate gel electrophoresis / migration,
incomplete (partial) digestion, methylation blocking, or circular-DNA topology, and it
depends on the built-in restriction-enzyme catalog. Digestion uses forward-strand cut
positions only, so strand-paired palindromic site reports are collapsed to one cut for
fragmentation purposes.

## Sources

`docs/algorithms/MolTools/Restriction_Digest_Simulation.md` (spec), which cites
Wikipedia *Restriction digest / Restriction enzyme / Restriction map*, the Addgene
*Restriction Digest Protocol*, Roberts (1976) *Restriction endonucleases*, and REBASE.
See [[restriction-enzyme-filtering]] for the enzyme-selection sibling and the
blunt/sticky end-type partition that underlies the compatibility rules.
