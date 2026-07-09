---
type: concept
title: "Overlap-Layout-Consensus (OLC) assembly"
tags: [assembly, algorithm]
sources:
  - docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
source_commit: ad9a76ef0b2ef475bdc2d4d9e866bf22676f0f84
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: assembly-olc-001-evidence
      evidence: "Test Unit ID: ASSEMBLY-OLC-001 ... Overlap-Layout-Consensus (OLC) genome assembly — overlap detection (FindAllOverlaps) and OLC assembly (AssembleOLC)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:de-bruijn-graph-assembly
      source: assembly-olc-001-evidence
      evidence: "finding a Hamiltonian path is NP-complete ... motivates heuristic layout and the shift to de Bruijn / Eulerian formulations — OLC and DBG are the two fragment-assembly formulations (Compeau, Pevzner & Tesler 2011)"
      confidence: high
      status: current
---

# Overlap-Layout-Consensus (OLC) assembly

**Overlap-Layout-Consensus (OLC)** is one of the two canonical fragment-assembly paradigms
(the other being [[de-bruijn-graph-assembly]]). It reconstructs a genome from reads in three
named stages — **Overlap**, **Layout**, **Consensus** — over an *overlap graph* whose nodes are
whole reads. This is the anchor for the assembly **OLC** family, validated under test unit
**ASSEMBLY-OLC-001** (`FindAllOverlaps` for overlap detection + `AssembleOLC` for end-to-end
assembly). The literature-traced validation record is [[assembly-olc-001-evidence]];
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
artifact pattern.

## The three stages (source-traced)

Traced verbatim to Langmead's JHU OLC notes (p.4) and Compeau, Pevzner & Tesler (2011):

1. **Overlap** — build the **overlap graph**: each read is a node, and a directed edge A → B is
   drawn when a **suffix of A** matches a **prefix of B** with overlap length ≥ a minimum
   threshold `l`. Only the **longest** suffix/prefix match per ordered pair is reported
   (OLC p.10); edges are labeled with the overlap length. `FindAllOverlaps` is the discovery
   primitive; the per-pair overlap definition is the same one anchored in
   [[contig-merge-overlap-collapse]].
2. **Layout** — bundle non-branching stretches of the overlap graph into **contigs**. Exact
   layout is finding a **Hamiltonian path** (visit every node exactly once) through the overlap
   graph — which is **NP-complete** (Compeau 2011), so practical assemblers use heuristics:
   **transitive reduction** (remove edges inferrible from other edges, starting with those that
   skip one node, then two) and emitting **non-branching stretches** as contigs. Branching at an
   unresolvable repeat splits the layout into multiple contigs.
3. **Consensus** — for each contig, pick the most likely nucleotide at each column by
   **majority vote** (OLC p.28). This is the same column-wise operation anchored in
   [[consensus-sequence]].

## Complexity

- **Overlap via suffix tree** — `O(N + a)` where `N` is total input length and `a` = number of
  overlapping pairs.
- **All-pairs dynamic-programming overlap** — `O(d²n²) = O(N²)` where `d` = number of reads of
  length `n`, `N = dn`; worst case `a = O(d²)` (OLC p.10, p.16).

## Worked oracles (published)

| Dataset | Parameters | Expected |
|---------|-----------|----------|
| `GTACGTACGAT` 6-mers overlap graph | 6 distinct 6-mers, minOverlap 4 | exactly **12 directed edges** with overlap lengths 4/5 (`FindAllOverlaps`) — Langmead SCS p.24–25, re-derived |
| Unambiguous 5-overlap tiling | `AAAAACCCCC`, `CCCCCGGGGG`, `GGGGGTTTTT` | single contig `AAAAACCCCCGGGGGTTTTT` (length 20) |
| Non-overlapping reads | 3 reads, no edge ≥ threshold | 3 singleton contigs |
| Single suffix-prefix overlap | X=`CTCTAGGCC`, Y=`TAGGCCCTC`, l=3 | longest overlap 6 (`TAGGCC`) — OLC p.5 |

`FindAllOverlaps` never emits self-overlaps (`ReadIndex1 != ReadIndex2`) and reports the longest
suffix-prefix match ≥ minOverlap.

## Failure modes (documented)

- **NP-completeness of exact layout** — a Hamiltonian path through the overlap graph is
  NP-complete; no exact polynomial OLC layout exists in general. This is the historical
  motivation for the shift to the Eulerian [[de-bruijn-graph-assembly]] formulation.
- **Unresolvable repeats split contigs** — a repeat longer than the read length creates branching
  in the overlap graph; the layout breaks into multiple contigs rather than one (OLC p.25).
- **Spurious subgraphs from sequencing error** — error-induced mismatches create dead-end branches
  that must be pruned (OLC p.26).
- **Greedy layout is suboptimal** — greedy maximal-overlap merging (greedy-SCS) can yield a
  longer-than-minimal superstring (SCS p.57: same input gives `AAABBBA` len 7 or `AAABBABBB`
  len 9 depending on merge order); it is a heuristic, not an exact shortest-common-superstring
  solver.
- **Repeats below resolution length collapse** — reads/k-mers shorter than the repeat period
  cannot distinguish repeat copies (SCS p.58–62).

## Assumptions (from the artifact)

- **Exact-match (identity 1.0) overlap for canonical numeric cases.** The published numeric
  oracles (the `GTACGTACGAT` graph, greedy traces) are stated for error-free reads. The sources
  also discuss approximate overlaps (mismatch/gap via DP, OLC p.11–15); the repository
  `minIdentity` parameter generalizes this, and a separate test exercises the identity threshold
  (7/8 = 0.875 accepted at 0.85, rejected at 0.95).
- **Empty read set → empty `AssemblyResult`.** No source specifies empty-input behavior; the
  trivial identity case, treated as a trivially-correct edge (mirrors `AssembleDeBruijn`).

## Relation to the other assembly formulations

OLC and [[de-bruijn-graph-assembly]] are the **two fragment-assembly formulations**: OLC works on
an **overlap graph of whole reads** and its exact layout is a **Hamiltonian path** (NP-complete),
whereas DBG works on a **k-mer graph** and its reconstruction is an **Eulerian walk** (tractable,
O(|E|)). Within OLC, the **O**verlap/merge primitive is anchored in
[[contig-merge-overlap-collapse]] and the **C**onsensus stage in [[consensus-sequence]]; OLC is
the end-to-end paradigm that composes them.

No contradictions among the sources — Compeau, Pevzner & Tesler (2011) and both Langmead lecture
notes give the identical overlap-graph / Hamiltonian-path / three-stage account; the teaching
notes reflect the same primaries. The `GTACGTACGAT` edge weights and greedy traces were
re-derived from the suffix-prefix definition and match the source slides.
