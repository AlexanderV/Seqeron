---
type: concept
title: "Phylogenetic marker selection (single-copy core genes ranked by parsimony-informative sites)"
tags: [comparative-genomics, pan-genome, phylogenetics, algorithm]
mcp_tools:
  - select_phylogenetic_markers
sources:
  - docs/Evidence/PANGEN-MARKER-001-Evidence.md
source_commit: 955bde0590e52fa1c979f009a5965d15a3a44722
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pangen-marker-001-evidence
      evidence: "Test Unit ID: PANGEN-MARKER-001 ... Algorithm: Phylogenetic Marker Selection (single-copy core genes ranked by parsimony-informative sites)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:pan-genome-gene-clustering
      source: pangen-marker-001-evidence
      evidence: "SelectPhylogeneticMarkers operates on gene clusters; a marker is a single-copy core cluster (panX 'those gene clusters in which all strains are represented exactly once'), i.e. the clusters ClusterGenes (PANGEN-CLUSTER-001) produces"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:pan-genome-core-accessory-partition
      source: pangen-marker-001-evidence
      evidence: "Marker candidates are the single-copy CORE clusters (Roary 'at least 99% of samples'); markers are the single-copy subset of the core partition PANGEN-CORE-001 defines"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:ortholog-detection-reciprocal-best-hits
      source: pangen-marker-001-evidence
      evidence: "Roary: 'homologous groups containing paralogs are split into groups of true orthologs' using conserved gene neighborhood; phylogenetic markers must be single-copy orthologs, not paralog mixtures"
      confidence: medium
      status: current
---

# Phylogenetic marker selection (single-copy core genes ranked by parsimony-informative sites)

`SelectPhylogeneticMarkers` is the **pan-genome family**'s (`PANGEN-*`) marker-selection unit: from a
set of gene **clusters** across N genomes it picks the **single-copy core genes** that make good
**phylogenetic markers** — clusters present **once in every genome** and carrying **variable
positions** — then **ranks them by descending parsimony-informative-site (PIS) count** and caps the
result at `maxMarkers`. Its companion primitive `CountParsimonyInformativeSites` scores the
column-wise phylogenetic signal of an aligned cluster. This is the input-selection step that
precedes core-genome tree building (panX → FastTree): it does **not** build a tree, it selects and
ranks the marker alignment columns a tree would be built from. It sits **downstream** of the gene
clustering step [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001, which produces the clusters) and
selects the **single-copy subset** of the core partition [[pan-genome-core-accessory-partition]]
(PANGEN-CORE-001) defines; its single-copy-ortholog requirement is the same true-ortholog constraint
[[ortholog-detection-reciprocal-best-hits]] enforces pairwise. Validated under test unit
**PANGEN-MARKER-001**; the validation record is [[pangen-marker-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

## The single-copy core marker rule (panX / Roary)

A cluster qualifies as a phylogenetic marker only if it is a **single-copy core gene**:

- **Single-copy** — every genome contributes **exactly one** gene to the cluster. panX defines the
  marker set verbatim as *"those gene clusters in which all strains are represented exactly once"*;
  Roary's documentation states *"any cluster containing paralogous genes gets filtered out of the
  final core gene alignment."* A genome contributing **0** genes (missing) or **≥ 2** genes
  (paralog) disqualifies the cluster.
- **Core** — present in (nearly) all genomes. Roary defines core as being in *"at least 99% of
  samples"* (default `coreFraction = 0.99`), which tolerates assembly error in large datasets. For
  small N this is strict (with N = 3 and 0.99, only 3/3 is core).
- **Variable** — the cluster must carry **≥ 1 parsimony-informative site**. panX *"extracts all
  variable positions from the nucleotide alignments of all single-copy core genes"* precisely
  because fully conserved columns carry **no phylogenetic signal**; a single-copy core cluster with
  **0 PIS** (fully conserved) is excluded.

Retained markers are **ranked by descending PIS** (most informative first) and the list is capped at
`maxMarkers`, mirroring panX's use of the most informative variable positions.

## Parsimony-informative sites (Zvelebil & Baum 2008)

A **parsimony-informative site** is *"a position in the relevant set of aligned sequences at which
there are at least two different character states and each of those states occurs in at least two of
the sequences."* Equivalently, a column is informative iff it has **≥ 2 distinct states AND at least
two of those states each appear in ≥ 2 sequences**. Non-informative columns:

- **Monomorphic** — all sequences share one state → 0 signal.
- **Singleton** — exactly one sequence differs (a variant in only one row) → not informative;
  it implies the same number of evolutionary changes regardless of tree topology.

Worked column oracle — four aligned sequences `s1=AAAAA`, `s2=AAACA`, `s3=AACCG`, `s4=ACCTG`:

| Column | s1 | s2 | s3 | s4 | States | PIS? | Reason |
|--------|----|----|----|----|--------|------|--------|
| 1 | A | A | A | A | A:4 | No | monomorphic |
| 2 | A | A | A | C | A:3,C:1 | No | singleton (C once) |
| 3 | A | A | C | C | A:2,C:2 | **Yes** | 2 states, each ≥ 2 |
| 4 | A | C | G | T | 4×1 | No | 4 singletons (no state ≥ 2) |
| 5 | A | A | G | G | A:2,G:2 | **Yes** | 2 states, each ≥ 2 |

**Total PIS = 2** (columns 3 and 5). The count is a pure column-content property: symmetric to row
order and invariant under relabeling states (A↔C).

## Selection oracle (single-copy core + PIS)

Three genomes g1,g2,g3, `coreFraction = 0.99` → a marker must be in all 3 genomes, one gene each,
with PIS ≥ 1:

| Cluster | Genomes | Genes/genome | Single-copy core? | Outcome |
|---------|---------|--------------|-------------------|---------|
| paralog | g1,g2,g3 | g1 has 2 | No (paralog in g1) | excluded |
| not-core | g1,g2 | 1 each | No (absent in g3) | excluded |
| conserved | g1,g2,g3 | 1 each | Yes, but 0 PIS | excluded (no signal) |

A positive marker across **four** single-copy genomes with members `ACGT`,`ACGT`,`TCGG`,`TCGG`:
col1 A,A,T,T (A:2,T:2 → PI) and col4 T,T,G,G (PI) → **PIS = 2**, retained.

## Edge cases

- **Null / empty inputs** (no clusters, no genomes) → empty marker set, no exception.
- **Non-single-copy** (a genome with 0 or ≥ 2 genes) → excluded.
- **Below core threshold** (absent from ≥ 1 genome) → excluded.
- **Zero PIS** (fully conserved single-copy core) → excluded.

## Assumption (source-backed)

**Per-cluster member alignment = equal-length ungapped columns (assumption).** panX/Roary align each
single-copy core cluster (MAFFT/PRANK) before reading columns. This unit has **no in-repo multiple
aligner**, so PIS is counted directly over the cluster's member sequences **when they share a common
length** (treated as already aligned, position-by-position) — the same ungapped, position-wise
convention the CD-HIT global identity of [[pan-genome-gene-clustering]] uses elsewhere in
`PanGenomeAnalyzer`. When member sequences differ in length, no common alignment exists and **PIS is
0** (the cluster carries no usable column-wise signal). This affects only how the alignment is
obtained, not the parsimony-informative criterion (copied verbatim from Zvelebil 2008) or the
single-copy-core selection rule (both fully source-backed).

## Caveat

The selected markers feed a core-genome tree, but selection **does not guarantee a correct tree**:
panX notes the core-genome tree *"may not reflect true evolutionary history due to homologous
recombination affecting even core genes."* This unit picks and ranks markers; it makes no
phylogenetic-accuracy claim.

## Reference tools

Definitions trace to **Ding, Baumdicker & Neher 2018** (panX, *Nucleic Acids Research* 46(1):e5 —
single-copy core "exactly once", variable-position extraction, FastTree, recombination caveat),
**Page et al. 2015** (Roary, *Bioinformatics* 31(22):3691–3693 — 99% core rule, paralog splitting
into true orthologs), the **Roary documentation** (Sanger Pathogens — paralog clusters filtered out
of the core gene alignment), and **Zvelebil & Baum 2008** (*Understanding Bioinformatics*, Garland
Science — the parsimony-informative-site definition, via the Wikipedia "Informative site" article).
No source contradictions — the single-copy-core selection rule and the parsimony-informative
criterion are mutually consistent across panX, Roary, and Zvelebil.
