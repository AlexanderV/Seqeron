---
type: source
title: "Evidence: ASSEMBLY-DBG-001 (De Bruijn Graph Assembly)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-DBG-001-Evidence.md
source_commit: ccacb64ccd743b9ec2d044971ec779c003f3b73c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-DBG-001

The validation-evidence artifact for test unit **ASSEMBLY-DBG-001** — de Bruijn graph genome
assembly: graph construction (`BuildDeBruijnGraph`) and Eulerian-walk reconstruction
(`AssembleDeBruijn`). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm's construction rule,
Euler theorems, spelling rule and failure modes are summarized in [[de-bruijn-graph-assembly]],
the anchor for the assembly DBG family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13):
  - **Langmead (JHU) — "De Bruijn Graph assembly" lecture notes** (rank 3 teaching notes, cite
    Jones & Pevzner for the Euler theorems; PDF extracted page-by-page with `pypdf`) — the
    k-mer / (k-1)-mer node/edge construction, directed-multigraph definition, balanced /
    semi-balanced nodes, the "Eulerian iff ≤2 semi-balanced" theorem, perfect-sequencing-is-
    Eulerian, the `path[0] + last-char-of-each-subsequent-node` spelling rule, the O(|E|)
    Hierholzer walk, and every worked example / failure mode below.
  - **Jones & Pevzner (2004), *An Introduction to Bioinformatics Algorithms*, MIT Press**
    (rank 1 textbook) — **Theorem 8.1** (Eulerian cycle iff all vertices balanced),
    **Theorem 8.2** (Eulerian path iff ≤2 semibalanced vertices, rest balanced), and §8.9
    Eulerian-path→cycle reduction; fragment assembly (SBH) reduces to finding an Eulerian path,
    linear in the number of edges.
  - **Compeau, Pevzner & Tesler (2011), *Nature Biotechnology* 29(11):987-991, DOI 10.1038/nbt.2023**
    (rank 1 paper) — the assembly application: reads as a de Bruijn graph (k-mers as edges
    between (k-1)-mer nodes) makes assembly a tractable Eulerian-path problem vs the NP-complete
    Hamiltonian-path overlap formulation. **Note:** both PDF mirrors are scanned/image-only
    (no extractable text), so only the publication metadata + abstract are cited from it; the
    behavioral theory is taken from the two text-extractable sources above.
- **Datasets** (all verified computationally with a Hierholzer Eulerian-walk reimplementation):
  - `AAABBBA`, k=3 → nodes `{AA,AB,BB,BA}`, edges `AA→AA, AA→AB, AB→BB, BB→BB, BB→BA` (degrees:
    AA in1/out2 source, BA in1/out0 sink, AB balanced, BB in2/out2 balanced ⇒ exactly two
    semi-balanced ⇒ Eulerian); walk `AA→AA→AB→BB→BB→BA` spells `AAABBBA`.
  - `a_long_long_long_time`, k=5 → 18-node walk reproduces the 21-char input verbatim.
  - `to_every_thing_turn_turn_turn_there_is_a_season` → correct at k=4; mis-ordered at k=3 (the
    `turn` repeat unresolvable at k=3).
  - `ATGGCGTGCA`, k=4 → exact (repeat-free DNA; all 4-mers distinct, unique Eulerian path).
- **Corner cases / failure modes** — repeat ≥ k-1 ⇒ multiple Eulerian walks (only one is the true
  genome), resolvable by larger k; coverage gap (omitted k-mer) ⇒ disconnected graph ⇒ multiple
  contigs; extra k-mer copy ⇒ 4 semi-balanced nodes ⇒ non-Eulerian; sequencing error ⇒ largest
  component non-Eulerian; the Eulerian framing is not practical at scale (De Bruijn Superwalk
  Problem is NP-hard — Medvedev 2007, Kingsford 2010). The multiedge case: `AAABBBBA`, k=3 adds an
  extra `BB→BB` self-loop.
- **Recommended coverage** — MUST: `BuildDeBruijnGraph(["AAABBBA"],3)` yields the exact node/edge
  set; `AssembleDeBruijn` reconstructs `AAABBBA` (k=3), `a_long_long_long_time` (k=5) and
  `to_every…` (k=4) exactly; `AAABBBBA` (k=3) produces the multiedge. SHOULD: read < k ⇒ no k-mers;
  empty/null reads ⇒ empty `AssemblyResult`; each distinct input k-mer ⇒ exactly one edge. COULD:
  every input read is a substring of the single unique-walk contig; a repeated (k-1)-mer produces a
  branch node (out-degree ≥ 2).

## Assumptions (from the artifact)

Three assumption records: (1) **walk selection among multiple Eulerian walks is unspecified** by
the sources — exact reconstruction is asserted ONLY on unique-walk inputs (no repeated (k-1)-mer);
non-unique inputs are checked against source-guaranteed invariants and the documented branch
structure, not a specific walk/string. (2) **empty / null read set → empty `AssemblyResult`**
(trivial identity, mirrors `AssembleOLC`). (3) **reads shorter than k contribute no k-mers** —
direct consequence of the `chop` bound `range(0, len(st)-(k-1))`, not a separate modelling choice.

No contradictions among the sources — Langmead's notes derive from and cite the same Jones &
Pevzner Euler theorems (8.1/8.2) that Compeau, Pevzner & Tesler (2011) build the assembly
application on. The only image-only source (Compeau 2011) is cited for metadata only, with the
behavioral theory drawn from the two text-extractable sources.
