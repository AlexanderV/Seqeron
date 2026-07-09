---
type: concept
title: "Pan-genome gene clustering (greedy incremental, CD-HIT model)"
tags: [comparative-genomics, pan-genome, algorithm]
sources:
  - docs/Evidence/PANGEN-CLUSTER-001-Evidence.md
source_commit: cce387f34d67ce7348bec08c7cea53dfd1d5cb64
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: pangen-cluster-001-evidence
      evidence: "Test Unit ID: PANGEN-CLUSTER-001 ... Algorithm: Gene Clustering (homolog / ortholog grouping by sequence identity)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:ortholog-detection-reciprocal-best-hits
      source: pangen-cluster-001-evidence
      evidence: "Roary converts coding regions to proteins, 'iteratively pre-clustered with CD-HIT,' then all-against-all BLASTP + MCL; this clusterer produces homolog groups, conserved gene neighbourhood later splits paralogs from true orthologs (out of scope for a greedy identity clusterer)"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:genome-comparison-core-dispensable
      source: pangen-cluster-001-evidence
      evidence: "Page et al. 2015 (Roary) uses CD-HIT gene clustering as the pre-clustering step of prokaryote pan-genome analysis; gene clusters are the homolog gene families the pan-genome core/dispensable partition is built from"
      confidence: medium
      status: current
---

# Pan-genome gene clustering (greedy incremental, CD-HIT model)

`ClusterGenes` is the **pan-genome family**'s (`PANGEN-*`) gene-clustering unit: it groups a set of
genes (across one or many genomes) into **homolog / ortholog families** by **sequence identity**,
following the **greedy incremental** procedure of **CD-HIT** (Li & Godzik 2006). Each output cluster
is a family of similar sequences with one **representative** and reports a `GenomeCount` (how many
distinct genomes contributed) and an `AverageIdentity` (members vs the representative). It is the
sequence-identity, single-collection clustering step that pan-genome pipelines such as **Roary**
(Page et al. 2015) run to build gene families — distinct from the pairwise best-hit
[[ortholog-detection-reciprocal-best-hits]] and the positional
[[conserved-gene-clusters-common-intervals]], and upstream of the pan-genome core/dispensable
partition [[genome-comparison-core-dispensable]]. Validated under test unit **PANGEN-CLUSTER-001**;
the validation record is [[pangen-cluster-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

## The greedy incremental procedure (CD-HIT)

The CD-HIT Algorithm wiki describes the method this unit follows:

1. **Sort the input sequences long → short** and process them from longest to shortest.
2. The **first (longest) unassigned sequence becomes a cluster representative**.
3. Each remaining query is **compared to the representatives found before it** and is classified as
   *redundant* (joins an existing cluster) or a *new representative*, based on whether it meets the
   identity threshold against one of them.
4. **First-match (fast mode):** the query joins the **first** representative that meets the threshold
   — a *greedy* assignment, **not** a best-hit search over all representatives.

Consequences documented in the source:

- **The representative is the longest member** (input is sorted long→short and the first unassigned
  sequence becomes the representative).
- **Deterministic** for a fixed input, because the long→short processing order is fixed by a stable
  sort.
- Members are compared to the **cluster representative**, not all-pairs (the `.clstr` output reports
  `%` = identity between each member and the representative).

## Sequence-identity metric (shorter-length denominator)

Grouping is by **sequence identity**, not k-mer similarity. CD-HIT's **default global identity**
(`-c`, default 0.9; `-G 1`) is:

```
identity = (identical residues in the alignment) / (full length of the SHORTER sequence)
```

The `-G 0` **local identity** (identical residues / alignment length) is the non-default contrast.
CD-HIT's **shorter-sequence denominator** is the convention tied to this clustering model. EMBOSS
needle corroborates the **numerator** (count of identical aligned positions); its denominator
(alignment length) differs by convention.

The threshold is **inclusive** (`≥`): a sequence joins when identity meets or exceeds the cutoff, so
`idThreshold = 1.0` clusters only sequences that are exactly identical over the shorter length.

Worked identity values (CD-HIT shorter-length denominator):

| seq1 | seq2 | identical / shorter len | identity |
|------|------|-------------------------|----------|
| `ATGCATGC` | `ATGCATGC` | 8 / 8 | 1.0 |
| `ATGCATGC` | `ATGCATGG` | 7 / 8 | 0.875 |
| `ATGCATGC` | `ATGCATGCAAAA` | 8 / 8 | 1.0 |
| `ATGC` | `CGTA` | 0 / 4 | 0.0 |

## Outputs and invariants

- **Cluster** = a homolog/ortholog gene family: representative + members, each member ≥ threshold to
  the representative.
- **`GenomeCount`** = number of distinct genomes contributing members (identical genes across
  genomes give `GenomeCount` = number of distinct genomes).
- **`AverageIdentity`** = mean member-vs-representative identity; a **singleton cluster's
  `AverageIdentity` = 1.0** (a sequence is 100% identical to itself).
- Empty/null genomes → **no clusters**; a null inner gene list is skipped (input-validation
  contract).

## Documented oracle (greedy clustering)

Sequences length-sorted long→short, threshold 0.8, ungapped global identity = identical positions /
shorter length:

- `Q2 = AAAAAAAAAAAA` (len 12, longest → becomes representative)
- `R = AAAAAAAAAA` (len 10; vs Q2 10/10 over shorter 10 = 1.0 → joins Q2)
- `Q1 = AAAAAAAAAT` (len 10; vs Q2 rep 9/10 = 0.9 ≥ 0.8 → joins Q2)
- `Q3 = CCCCCCCCCC` (len 10; vs Q2 0/10 = 0.0 → starts a new cluster)

Result: **2 clusters** — `{Q2, R, Q1}` and `{Q3}`. Lowering/raising the threshold merges or splits
near-identical members.

## Deviations and assumptions (source-backed)

1. **Ungapped alignment (deviation).** CD-HIT computes identity over a *banded* alignment; this unit
   computes identity over an **ungapped positional comparison of the shared prefix**, divided by the
   shorter length. This **matches CD-HIT exactly** for sequences differing only by substitutions or
   by a length difference (no internal indels); for sequences that need internal gaps to align,
   identity may be **underestimated**. Recorded as a deviation in the algorithm doc §5.3 / §5.4 —
   implementing CD-HIT's full banded aligner is out of scope for one unit; the numerator (identical
   positions) and denominator (shorter length) are taken verbatim from the source.
2. **No paralog splitting (assumption).** Roary splits homolog groups into orthologs using **gene
   neighbourhood** (synteny); this clusterer produces **homolog groups only** (no synteny step),
   matching CD-HIT pre-clustering rather than the full Roary pipeline.

## Reference tools

Definitions trace to **Li & Godzik 2006** (CD-HIT, *Bioinformatics* 22(13):1658–1659) plus the
**CD-HIT User's Guide** (`-c`/`-G 1` global-identity formula, cluster/representative structure) and
the **CD-HIT Algorithm wiki** (greedy incremental long→short procedure, first-match fast mode);
**Page et al. 2015** (Roary, *Bioinformatics* 31(22):3691–3693) supplies the pan-genome context and
the default 95% identity for ortholog grouping; **EMBOSS needle** corroborates the
identical-aligned-positions numerator. No source contradictions — the identity conventions,
greedy procedure, and pan-genome framing are mutually consistent.
