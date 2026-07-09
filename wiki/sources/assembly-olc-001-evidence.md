---
type: source
title: "Evidence: ASSEMBLY-OLC-001 (Overlap-Layout-Consensus)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-OLC-001-Evidence.md
source_commit: ad9a76ef0b2ef475bdc2d4d9e866bf22676f0f84
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-OLC-001

The validation-evidence artifact for test unit **ASSEMBLY-OLC-001** — Overlap-Layout-Consensus
genome assembly: overlap detection (`FindAllOverlaps`) and end-to-end OLC assembly (`AssembleOLC`).
One instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the three-stage paradigm, overlap graph, Hamiltonian-path layout, oracles and assumptions
are summarized in [[overlap-layout-consensus-assembly]], the anchor for the assembly OLC family.
See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13; PDFs fetched with WebFetch, binary PDFs extracted
  locally with `pypdf`):
  - **Compeau, Pevzner & Tesler (2011), "How to apply de Bruijn graphs to genome assembly",
    *Nature Biotechnology* 29:987–991 (rank 1, peer-reviewed)** — the overlap graph (read = node,
    edge A→B when suffix of A matches prefix of B above a threshold), layout = **Hamiltonian path**
    (visit each node once) reconstructs the genome, and the **NP-completeness** of finding a
    Hamiltonian path (the motivation for the shift to de Bruijn / Eulerian formulations).
  - **Langmead (JHU) — "Overlap Layout Consensus assembly" notes (rank 3)** — the three stages
    (Overlap = build overlap graph / Layout = bundle stretches into contigs / Consensus = majority
    vote per column), the overlap definition (longest suffix of X = prefix of Y, length ≥ `l`,
    p.5), report-only-longest (p.10), edge weights labeled with overlap length (p.20–25), layout
    **transitive reduction** (remove edges that skip one node, then two, p.21–24), contig emission
    from non-branching stretches (p.25), and complexity (suffix-tree O(N+a); all-pairs DP
    O(d²n²)=O(N²), p.10/p.16).
  - **Langmead (JHU) — "Assembly & Shortest Common Superstring" notes (rank 3)** — the first law of
    assembly (suffix of A similar to prefix of B ⇒ possible overlap), the overlap graph worked
    example, **greedy-SCS** (each round merge the max-overlap pair; not optimal — same input yields
    `AAABBBA` len 7 or `AAABBABBB` len 9), and repeats-foil-assembly (repeats longer than read
    length are unresolvable).
- **Datasets (published oracles):**
  - `GTACGTACGAT` 6 distinct 6-mers, minOverlap 4 → exactly **12 directed edges** (overlap lengths
    4 and 5), tabulated A→B:length; edge weights match the source slide (Langmead SCS p.24–25,
    re-derived from the suffix-prefix definition).
  - Unambiguous 5-overlap tiling `AAAAACCCCC`,`CCCCCGGGGG`,`GGGGGTTTTT` → single contig
    `AAAAACCCCCGGGGGTTTTT` (length 20).
  - Single suffix-prefix overlap X=`CTCTAGGCC`, Y=`TAGGCCCTC`, l=3 → longest overlap 6 (`TAGGCC`),
    OLC p.5.
- **Corner cases / failure modes** — NP-complete exact layout (no polynomial OLC layout in
  general); repeats longer than read length branch the graph and split contigs (p.25); sequencing
  errors create spurious dead-end subgraphs to prune (p.26); report only the longest overlap per
  ordered pair (p.10); greedy layout is suboptimal (p.57); repeats below resolution length collapse
  (p.58–60).
- **Recommended coverage** — MUST: `FindAllOverlaps` on the `GTACGTACGAT` 6-mers returns exactly
  the 12 edges with correct lengths; never emits self-overlaps and reports the longest match ≥
  minOverlap; `AssembleOLC` on the 5-overlap tiling → single `AAAAACCCCCGGGGGTTTTT`; on three
  non-overlapping reads → 3 singleton contigs. SHOULD: identity threshold gates acceptance
  (7/8=0.875 accepted at 0.85, rejected at 0.95); `MinOverlap` boundary (exactly = accepted, one
  below rejected); contig-length invariant (≤ sum of read lengths, ≥ longest read). COULD: repeat
  limitation does not falsely collapse into one too-short contig.

## Assumptions (from the artifact)

Two assumption records: (1) **exact-match (identity 1.0) overlap for canonical numeric cases** —
the published numeric oracles are stated for error-free reads; the sources also cover approximate
overlaps (mismatch/gap via DP, OLC p.11–15), which the repository `minIdentity` parameter
generalizes (a separate test exercises the threshold). (2) **empty read set → empty
`AssemblyResult`** — not source-specified; the natural identity edge case.

No contradictions among the sources — Compeau, Pevzner & Tesler (2011) and both Langmead notes
give the identical overlap-graph / Hamiltonian-path / three-stage account; the teaching notes
reflect the same primaries and the re-derived numeric oracles match the source slides.
