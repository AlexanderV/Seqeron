---
type: source
title: "Evidence: PANGEN-MARKER-001 (Phylogenetic marker selection — single-copy core genes ranked by parsimony-informative sites)"
tags: [validation, comparative-genomics, pan-genome, phylogenetics]
doc_path: docs/Evidence/PANGEN-MARKER-001-Evidence.md
sources:
  - docs/Evidence/PANGEN-MARKER-001-Evidence.md
source_commit: 955bde0590e52fa1c979f009a5965d15a3a44722
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PANGEN-MARKER-001

The validation-evidence artifact for test unit **PANGEN-MARKER-001** — **phylogenetic marker
selection**: from gene clusters across N genomes, select the **single-copy core genes** (present
exactly once in every genome) that carry **variable positions**, score each by
**parsimony-informative sites (PIS)**, and return them **ranked by descending PIS** capped at
`maxMarkers`. This is a **pan-genome family** (`PANGEN-*`) Evidence file and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm,
its selection rule, the parsimony-informative criterion, worked oracles, and the alignment
assumption are summarized in the dedicated concept [[phylogenetic-marker-selection]]. It selects the
single-copy subset of the core partition [[pan-genome-core-accessory-partition]] (PANGEN-CORE-001)
and consumes the clusters of [[pan-genome-gene-clustering]] (PANGEN-CLUSTER-001); the single-copy
constraint is the true-ortholog requirement [[ortholog-detection-reciprocal-best-hits]] enforces
pairwise. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Ding, Baumdicker & Neher 2018** (panX, *Nucleic Acids Research* 46(1):e5, authority 1) —
    single-copy core gene = *"those gene clusters in which all strains are represented exactly
    once"*; panX *"extracts all variable positions from the nucleotide alignments of all single-copy
    core genes"* to build a core-genome tree with FastTree; invariant columns carry no signal;
    caveat that the core tree *"may not reflect true evolutionary history due to homologous
    recombination affecting even core genes."*
  - **Page et al. 2015** (Roary, *Bioinformatics* 31(22):3691–3693, authority 1) — core = *"at least
    99% of samples"* (default 0.99, tolerates assembly error); *"homologous groups containing
    paralogs are split into groups of true orthologs"* via conserved gene neighborhood.
  - **Roary documentation** (Sanger Pathogens, authority 3) — *"any cluster containing paralogous
    genes gets filtered out of the final core gene alignment"*; the core gene alignment is built from
    per-cluster members (MAFFT/PRANK), signal read off the aligned columns.
  - **Zvelebil & Baum 2008** (*Understanding Bioinformatics*, Garland Science, authority 1 via the
    Wikipedia "Informative site" article, authority 4) — parsimony-informative site = *"a position …
    at which there are at least two different character states and each of those states occurs in at
    least two of the sequences"*; fully-conserved and singleton columns are **not** informative.
- **Corner cases / failure modes:** non-single-copy (a genome with 0 or ≥ 2 genes) excluded;
  below-core (absent from ≥ 1 genome under threshold) excluded; monomorphic column → 0 PIS; singleton
  column → not informative; two-states-each-≥2 → 1 PIS; PIS defined over aligned columns → sequences
  must share a common length.
- **Datasets (documented oracles):**
  - *PIS worked columns (Zvelebil 2008)* — s1=`AAAAA`, s2=`AAACA`, s3=`AACCG`, s4=`ACCTG`: cols 3 and
    5 informative, cols 1 (mono) / 2 (singleton) / 4 (four singletons) not → **PIS = 2**.
  - *Single-copy core selection (panX/Roary)* — 3 genomes, coreFraction 0.99: `paralog` (g1 has 2
    genes) excluded, `not-core` (absent in g3) excluded, `conserved` (single-copy core but 0 PIS)
    excluded; a positive marker across four single-copy genomes `ACGT`/`ACGT`/`TCGG`/`TCGG` → PIS = 2
    (cols 1 and 4).
- **Coverage recommendations:** MUST-test `CountParsimonyInformativeSites` = 2 on the 5-column
  alignment, the monomorphic/singleton/two-state-each-≥2 column classifications, `SelectPhylogenetic-`
  `Markers` exclusion of non-single-copy (paralog) and non-core clusters and of 0-PIS conserved
  clusters, descending-PIS ranking capped at `maxMarkers`, and null/empty → empty set no exception;
  SHOULD-test unequal-length members → PIS 0 → not selected; COULD-test PIS symmetry to row order and
  to state relabeling.

## Deviations and assumptions

One assumption, source-backed, no deviations from the literature. **Per-cluster member alignment =
equal-length ungapped columns (assumption).** panX/Roary align each single-copy core cluster
(MAFFT/PRANK) before reading columns; this unit has no in-repo multiple aligner, so PIS is counted
directly over the cluster's member sequences **when they share a common length** (treated as already
aligned, position-by-position — the same ungapped position-wise convention the CD-HIT global identity
of [[pan-genome-gene-clustering]] uses elsewhere in `PanGenomeAnalyzer`). Unequal-length members →
no common alignment → **PIS = 0**. This affects only how the alignment is obtained, not the
parsimony-informative criterion (verbatim Zvelebil 2008) or the single-copy-core selection rule
(both fully source-backed). No source contradictions.
