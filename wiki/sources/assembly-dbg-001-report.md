---
type: source
title: "Validation report: ASSEMBLY-DBG-001 (de Bruijn graph assembly — k-mer graph + Eulerian-walk reconstruction)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-DBG-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-DBG-001.md
source_commit: 131c8e266fdd08713526890d833f52901b803517
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-DBG-001

The two-stage **validation write-up** for test unit **ASSEMBLY-DBG-001** (de Bruijn graph
assembly — build a de Bruijn graph from k-mers and traverse it into contigs), validated
2026-06-15 in a fresh context. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm description
and the shipped code. The algorithm itself is summarized in the concept
[[de-bruijn-graph-assembly]] (anchor of the assembly DBG family), and the wider campaign is
[[validation-and-testing]]. Distinct from [[assembly-dbg-001-evidence]] (the pre-implementation
evidence artifact, sourced from `docs/Evidence/`) — this is the independent re-validation verdict.

Canonical methods under test:
`SequenceAssembler.BuildDeBruijnGraph(reads, k)` (public (k-1)-mer out-adjacency multigraph) and
`SequenceAssembler.AssembleDeBruijn(reads, parameters)` (Eulerian-walk reconstruction → contigs +
stats), with private helpers `ReconstructContigs` / `ChooseEulerStart` / `SpellEulerianWalk`
(iterative Hierholzer).

## Verdict

**Stage A: ✅ PASS · Stage B: 🟡 PASS-WITH-NOTES · State: ✅ CLEAN.** No implementation defect;
the code realises the validated formula exactly for every cross-check. The only issue was a
**test-coverage gap** (three documented but untested branches), found and closed in-session with
sourced-value tests. Full unfiltered suite **6497 passed, Failed: 0** (1 pre-existing skip,
`MFE_Benchmark_AllScenarios`, unrelated); DBG fixture 17 → 20 tests; `dotnet build` 0 errors
(4 pre-existing warnings in an untouched file).

## Stage A — description (algorithm faithfulness)

Theory checked against sources opened this session, independent of repo artifacts:

- **Langmead (JHU) — "De Bruijn Graph assembly" lecture notes** (PDF fetched and extracted
  page-by-page with `pypdf` this session). Every Stage-A claim was read directly off the slides:
  node = (k-1)-mer, edge = k-mer directed prefix→suffix, directed multigraph (repeated k-mer ⇒
  multiedge, p.8); balanced (indeg=outdeg) / semi-balanced (differ by 1) and "Eulerian iff ≤2
  semi-balanced, rest balanced" (p.10, citing Jones & Pevzner §8.8); perfect-sequencing-is-
  Eulerian (p.15); the `chop` bound `i ∈ [0, |r|−k]` and prefix/suffix split (p.16); the spelling
  rule `path[0] + concat(last-char of path[i>0])` (p.18-19); worked walks and the k=4-correct /
  k=3-mis-ordered `to_every…` repeat case (p.21-22); coverage-gap → disconnected-but-per-
  component-Eulerian (p.24).
- **Jones & Pevzner (2004), Theorems 8.1/8.2** — Eulerian path iff ≤2 semibalanced vertices, all
  others balanced; exactly the condition `ChooseEulerStart` realises.
- **Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987-991** — publication identity + thesis
  (k-mers as edges between (k-1)-mer nodes; assembly = polynomial Eulerian-path vs NP-complete
  Hamiltonian/OLC), used for framing only.

**Formula check.** Node = (k-1)-mer, edge = k-mer (prefix→suffix), multigraph; chop bound
`i ∈ [0,|r|−k]` so #k-mers per read = |r|−k+1; balanced/semi-balanced + Euler condition; spelling
`p₀ + last-char(pᵢ, i>0)` — each matches a source verbatim. Edge semantics sourced: read < k ⇒ 0
k-mers (p.16); k < 2 rejected (nodes would be empty, follows from the (k-1)-mer definition);
repeated k-mer ⇒ multiedge (p.8); disconnected graph ⇒ one Eulerian walk per component (p.24).
Two assumptions are flagged as not source-prescribed: walk selection among multiple Eulerian
walks (ASSUMPTION-1, consistent with p.21 — only invariants asserted, never a specific wrong
walk) and empty/null reads → empty result (ASSUMPTION-2, trivial identity).

**Independent cross-check (numbers).** Re-derived an independent Hierholzer DBG assembler in
Python (separate from repo code) and reproduced every published value: `AAABBBA` k=3 → `AAABBBA`
(5 = 7−3+1 edges); `AAABBBBA` k=3 → 6 edges, BBB twice ⇒ two BB→BB (multiedge); `a_long_long_
long_time` k=5 → verbatim (Langmead p.18 walk); `to_every…season` k=4 → input exactly (p.22);
`ATGGCGTGCA` k=4 → itself (all 4-mers distinct ⇒ unique walk). **Note:** Langmead's printed
k=3 mis-ordered `to_every…` comes from his particular arbitrary walk; under this library's
deterministic insertion-order Hierholzer walk the k=3 input is recovered — the spec correctly
anticipates this (ASSUMPTION-1) and the k=3 test asserts only branching structure, not a string.
Stage A = **PASS**, no divergences.

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs`:
`AssembleDeBruijn` L104-125 (null/empty guard → build → reconstruct → `MinContigLength` filter →
stats); `BuildDeBruijnGraph` L366-402 (k<2 throw, null guard, chop bound `i ≤ read.Length−k`,
prefix/suffix split, multimap append); `ReconstructContigs` L413-455 (degree tallies, weakly-
connected components via `BuildUndirectedAdjacency` L457-474 / `CollectComponent` L476-498, per-
component walk); `ChooseEulerStart` L506-524 (semi-balanced out−in=1 source, else lexicographically
smallest node); `SpellEulerianWalk` L532-569 (iterative Hierholzer, edges consumed in insertion
order, reverse, spell).

**Formula realised correctly (evidence).** Chop bound `for (int i = 0; i <= read.Length - k; i++)`
= Langmead's `range(0, len−(k−1))`; `prefix = read.Substring(i, k-1)`, `suffix = read.Substring
(i+1, k-1)` = `kmer[:-1], kmer[1:]`; multimap append preserves edge multiplicity; `ChooseEulerStart`
returns the unique out−in=1 node (in a valid Eulerian-path graph exactly one exists, so DFS-order
iteration is safe); `SpellEulerianWalk` builds `walk[0] + concat(last char of walk[i>0])`;
`adjacency` is a per-call mutable copy consumed by a per-node cursor, and components are edge-
disjoint, so sharing one `adjacency` across the per-component loop is correct.

**Cross-verification (real `AssembleDeBruijn` via a temporary console harness + NUnit fixture) —
all match the independent Python reference:**

| Input | k | Reference / code output | Source |
|-------|---|-------------------------|--------|
| `AAABBBA` | 3 | `AAABBBA` | Langmead p.11 |
| `a_long_long_long_time` | 5 | `a_long_long_long_time` | Langmead p.18 |
| `to_every…season` | 4 | input verbatim | Langmead p.22 |
| `ATGGCGTGCA` | 4 | `ATGGCGTGCA` | construction+spelling |
| `{ATGGCGTGCA, GATTACAGGTC}` | 4 | two contigs = each read | Langmead p.24 (disconnected) |
| same, `MinContigLength=11` | 4 | `{GATTACAGGTC}` only | contract §3.1 (filter) |
| `AAABBBA` edges | 3 | 5 edges | Langmead p.7 |
| `AAABBBBA` BB→BB | 3 | 2 (multiedge), 6 edges | Langmead p.8 |

**Variant/delegate consistency.** `AssembleDeBruijn` consumes `BuildDeBruijnGraph`; no
`*Fast`/instance variants. A determinism test confirms identical output across runs (lexicographically-
smallest start + insertion-order edges).

**Test-quality audit (HARD gate: PASS).** Every exact expected value (M1-M7 strings/edges,
multiedge count, edge counts, disconnected contigs, filter) traces to Langmead (verified from the
fetched PDF) or the independent Hierholzer reference — not to code output; a deliberately-wrong
walk/spelling would fail M4-M7 and the new disconnected test. Assertions use exact `Is.EqualTo`
on strings/edge lists/node sets/counts; S3 (INV-05) uses `Contains` only because INV-05 is
*defined* as a substring property. No skips, no widened tolerances, no expected values bent to
match output.

## Findings

- **No algorithm defect.** State ✅ CLEAN; the code already produced the correct sourced values
  for every case, including the three previously untested branches.
- **Test-coverage gap fixed in-session (Stage B defect).** The original 17-test fixture covered
  all M/S/C spec rows but not three documented Stage-A behaviours/branches: (1) **disconnected
  graph → one contig per weakly-connected component** (Langmead p.24; the `BuildUndirectedAdjacency`/
  `CollectComponent`/per-component loop); (2) the **`MinContigLength` filter** (contract §3.1; the
  `Where(c => c.Length >= MinContigLength)` branch — all prior tests used `MinContigLength:1`);
  (3) **`BuildDeBruijnGraph(null, k)`** null guard (contract §3.3). Added three exact-value tests
  with externally-derived oracles — `AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent`
  (`{ATGGCGTGCA, GATTACAGGTC}` k=4 → exactly those two contigs), `AssembleDeBruijn_MinContigLength_
  FiltersShortContigs` (same reads, MinContigLength=11 → only the length-11 `GATTACAGGTC` survives,
  the length-10 contig dropped), and `BuildDeBruijnGraph_NullReads_ProducesEmptyGraph`.
- **No follow-ups.**
