---
type: concept
title: "Contig merging (suffix–prefix overlap collapse / superstring merge)"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
  - docs/algorithms/Extended_Assembly/Contig_Merging.md
  - docs/Validation/reports/ASSEMBLY-MERGE-001.md
source_commit: 6abf4edca8f18ac8c0d17c25f3949d7c1dea135d
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-merge-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-MERGE-001 ... Contig Merging (suffix–prefix overlap merge / superstring collapse) ... MergeContigs(contig1, contig2, overlapLength)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:de-bruijn-graph-assembly
      source: assembly-merge-001-evidence
      evidence: "OLC stages: Overlap – Build overlap graph; Layout – Bundle stretches of the overlap graph into contigs; Consensus ... an overlap is a suffix/prefix match — the same overlap primitive the de Bruijn / Eulerian formulation avoids by using fixed k-mers"
      confidence: medium
      status: current
---

# Contig merging (suffix–prefix overlap collapse)

**Contig merging** is the low-level primitive that collapses two overlapping strings into one by
removing a single copy of their shared region. It is the **merge** step behind greedy
shortest-common-superstring (SCS) assembly and the **layout/consensus** stitching of
Overlap-Layout-Consensus (OLC). This is the anchor for the assembly **MERGE** family, validated
under test unit **ASSEMBLY-MERGE-001** (`MergeContigs`). The literature-traced validation record is
[[assembly-merge-001-evidence]], and the independent two-stage re-validation verdict (Stage A ✅ PASS /
Stage B ✅ PASS / State CLEAN, suite 6529/0, zero code change) is [[assembly-merge-001-report]];
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the artifact
pattern.

## The overlap definition (source-traced)

Traced verbatim to Langmead's JHU SCS/OLC notes and MIT 7.91J Lecture 6:

- **Overlap** — a length-`l` **suffix of X** exactly matches a length-`l` **prefix of Y**
  (`l` given). An overlap "exists when a suffix of X of length ≥ l exactly matches a prefix of Y".
- **`suffixPrefixMatch(x, y, k)`** — the discovery primitive returns the length of the longest
  suffix of `x` (length ≥ k) that matches a prefix of `y`, else `0`. It guards
  `if len(x) < k or len(y) < k: return 0`, so an overlap can never exceed `min(|X|, |Y|)`.
- **Longest match reported** — when several overlaps exist, only the longest suffix/prefix match is
  used; merging with that length removes exactly one copy of the overlapping region.

## The merge / collapse operation

- **Superstring by collapse** — "without requirement of 'shortest', it's easy: just concatenate
  them." Collapsing keeps a **single copy** of the overlapping region, so the merge is shorter than
  plain concatenation by the overlap length.
- **Length invariant** — `|merge(X, Y, l)| = |X| + |Y| − l`. The implementation
  `MergeContigs(contig1, contig2, overlapLength)` removes `overlapLength` characters from the front
  of `contig2` and appends the remainder to `contig1`.
- **Overlap 0 → plain concatenation** — when no length-≥-`l` match exists, `suffixPrefixMatch`
  returns 0 and the merge is exactly `X + Y`.

## Worked oracles (published)

| Inputs | Overlap `l` | Merged result | Length |
|--------|-------------|---------------|--------|
| `BAA` + `AAB` | 2 (suffix `AA` = prefix `AA`) | `BAAB` | 4 = 3 + 3 − 2 |
| `AAA`+`AAB`+`ABB`+`BBB`+`BBA`, each overlap 2 | chained | `AAABBBA` | 7 (full greedy-SCS) |
| `BAA` + `AAB` | 0 (no overlap) | `BAAAAB` | 6 (= plain concatenation) |

## Assumptions (from the artifact)

- **Caller-supplied overlap length is trusted, not re-verified.** `MergeContigs` is a collapse
  primitive; verifying that the collapsed region truly matched is the separate responsibility of
  `FindOverlap` / `FindAllOverlaps`. An unverified length changes only whether the region matched,
  not the merge arithmetic — an **API-contract** boundary, not a scoring/correctness parameter.
- **Out-of-range overlap (≤ 0 or > min length) → plain concatenation.** Follows directly from the
  two source facts (overlap 0 = "no overlap → concatenate"; a valid overlap is bounded by
  `min(|X|, |Y|)`); it is the only behavior consistent with the suffix/prefix definition.

## Relation to the other assembly formulations

Contig merging is the **overlap-based** primitive: it collapses variable-length suffix/prefix
matches, which is what [[overlap-layout-consensus-assembly|OLC]] and greedy SCS do. The
[[de-bruijn-graph-assembly]] approach sidesteps
this pairwise-overlap search by using fixed-length `k`-mers and an Eulerian walk instead (the
overlap search that makes the Hamiltonian/SCS formulation NP-hard is exactly what DBG avoids). The
Consensus stage that finishes an OLC layout is [[consensus-sequence]]. [[scaffolding]] hands off to
this primitive on its **negative-gap** case: a negative inter-contig distance estimate means the two
contigs should overlap, so if the overlap is found they are merged here rather than joined by an
`N`-gap.

No contradictions among the sources — Langmead's SCS notes, Langmead's OLC notes, and MIT 7.91J all
state the identical suffix-of-X / prefix-of-Y overlap definition; each corroborates the others.
