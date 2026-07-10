---
type: concept
title: "De Bruijn graph assembly (k-mer graph + Eulerian-walk reconstruction)"
tags: [assembly, algorithm]
mcp_tools:
  - assemble_de_bruijn
sources:
  - docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
  - docs/algorithms/Assembly/De_Bruijn_Graph_Assembly.md
  - docs/Validation/reports/ASSEMBLY-DBG-001.md
source_commit: 131c8e266fdd08713526890d833f52901b803517
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-dbg-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-DBG-001 ... De Bruijn graph genome assembly — graph construction (BuildDeBruijnGraph) and Eulerian-walk reconstruction (AssembleDeBruijn)"
      confidence: high
      status: current
---

# De Bruijn graph assembly

**De Bruijn graph (DBG) assembly** reconstructs a genome from short reads by turning it into a
graph problem: chop reads into `k`-mers, make each distinct **(k-1)-mer** a node, and add a
directed edge for every `k`-mer from its left (k-1)-mer prefix to its right (k-1)-mer suffix.
Reconstructing the genome is then an **Eulerian walk** — a walk that uses each edge exactly once
— and spelling that walk recovers the superstring. This is the anchor for the assembly **DBG**
family, validated under test unit **ASSEMBLY-DBG-001** (`BuildDeBruijnGraph` +
`AssembleDeBruijn`). The literature-traced validation record is [[assembly-dbg-001-evidence]]; the independent
two-stage re-validation verdict (Stage A PASS / Stage B PASS-WITH-NOTES / CLEAN — a closed
test-coverage gap, no code defect) is [[assembly-dbg-001-report]];
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
artifact pattern.

## Graph construction (source definitions)

Traced to Langmead's DBG lecture notes and the Jones & Pevzner textbook:

- **k-mer / (k-1)-mer** — a `k`-mer is a length-`k` substring; its left and right (k-1)-mers are
  its length-(k-1) prefix and suffix (`AAB` → `AA`, `AB`).
- **Nodes and edges** — nodes are the distinct (k-1)-mers; each `k`-mer contributes one directed
  edge (left prefix → right suffix). Each edge corresponds to exactly one length-`k` input
  substring, so an edge encodes a length-(k-2) overlap between two (k-1)-mers.
- **Directed multigraph** — repeated `k`-mers produce parallel edges (multiedges), so `G(V, E)`
  is a directed multigraph. `indegree` / `outdegree` count incoming / outgoing edges.
- **`chop` bound** — `k`-mers are `st[i:i+k]` for `i` in `range(0, len(st)-(k-1))`; a read shorter
  than `k` yields **zero** `k`-mers and is silently skipped.

## Eulerian-walk reconstruction

- **Balanced / semi-balanced** — a node is *balanced* when indegree = outdegree, *semi-balanced*
  when they differ by 1.
- **Euler's theorems (Jones & Pevzner 8.1/8.2)** — a connected graph has an **Eulerian cycle**
  iff every vertex is balanced, and an **Eulerian path** iff it has at most two semi-balanced
  vertices and all others are balanced. Fragment assembly reduces to finding an Eulerian path
  (contrast the NP-complete Hamiltonian-path / overlap formulation — see the OLC sibling).
- **Perfect sequencing ⇒ Eulerian** — with perfect coverage the graph is Eulerian: the left-end
  (k-1)-mer is the semi-balanced source (one extra out-edge), the right-end is the semi-balanced
  sink (one extra in-edge), all others balanced (unless the genome is circular).
- **Spelling the walk** — the walk is a node sequence; the genome is
  `path[0] + ''.join(node[-1] for node in path[1:])` — emit the first (k-1)-mer in full, then
  append the **last character** of every subsequent node.
- **Complexity** — an Eulerian walk is found in **O(|E|)** time (Hierholzer-style).

## Worked oracles (unique-walk inputs)

The unit asserts exact reconstruction only on inputs whose Eulerian walk is unique (no repeated
(k-1)-mer):

| Input | k | Reconstruction |
|-------|---|----------------|
| `AAABBBA` | 3 | `AAABBBA` (nodes `{AA,AB,BB,BA}`; edges `AA→AA, AA→AB, AB→BB, BB→BB, BB→BA`; walk `AA→AA→AB→BB→BB→BA`) |
| `a_long_long_long_time` | 5 | `a_long_long_long_time` (18-node walk) |
| `to_every_thing_turn_turn_turn_there_is_a_season` | 4 | exact (but **mis-ordered at k=3** — the `turn` repeat is unresolvable) |
| `ATGGCGTGCA` | 4 | `ATGGCGTGCA` (repeat-free DNA smoke case, unique path) |

Adding a `B` (`AAABBBBA`, k=3) is the canonical **multiedge** case: node `BB` gains an extra
self-loop.

## Failure modes (when the graph is not Eulerian)

- **Repeat ≥ k-1 → multiple Eulerian walks** — a (k-1)-mer occurring more than once becomes a
  shared node joining edge-disjoint cycles; several walks exist, only one is the true genome.
  Increasing `k` can resolve a repeat shorter than `k-1` (the k=4-correct / k=3-wrong `to_every…`
  case).
- **Coverage gaps → disconnected graph** — omitting a `k`-mer splits the graph; each component is
  individually Eulerian but the whole is not, so assembly yields multiple contigs.
- **Uneven coverage / sequencing errors → non-Eulerian** — an extra `k`-mer copy or an error
  creates additional semi-balanced nodes (>2), violating Theorem 8.2.
- **Not practical at scale** — real graphs are non-Eulerian from errors, repeats and uneven
  coverage; the De Bruijn Superwalk Problem is NP-hard. The Eulerian framing is the clean
  theoretical model, not the production algorithm.

## Assumptions (from the artifact)

- **Walk selection among multiple Eulerian walks is unspecified.** Sources only state the true
  genome is *one* of the valid walks, not which a deterministic implementation must emit. Exact
  reconstruction is therefore asserted only on unique-walk inputs; non-unique inputs are checked
  against source-guaranteed invariants (each input `k`-mer is spelled by some contig; total
  reconstructed length equals the sum over edges) and structural facts (a repeated (k-1)-mer
  produces a branch node with out-degree ≥ 2) — never a specific wrong string.
- **Empty / null read set → empty `AssemblyResult`** (trivial identity, mirrors `AssembleOLC`).
- **Reads shorter than `k` contribute no `k`-mers** — direct consequence of the `chop` bound, not
  a separate modelling choice.

No contradictions among the sources — Langmead's notes derive from and cite the same Jones &
Pevzner Euler theorems that Compeau, Pevzner & Tesler (2011) build the assembly application on.
The [[overlap-layout-consensus-assembly|OLC (Overlap-Layout-Consensus)]] approach is the
alternative fragment-assembly formulation (Hamiltonian-path layout over an overlap graph vs DBG's
Eulerian walk over a k-mer graph); its suffix–prefix overlap/merge primitive is covered by
[[contig-merge-overlap-collapse]] and its consensus step by [[consensus-sequence]].
