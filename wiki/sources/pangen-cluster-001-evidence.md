---
type: source
title: "Evidence: PANGEN-CLUSTER-001 (Gene clustering — greedy incremental homolog grouping by identity)"
tags: [validation, comparative-genomics, pan-genome]
doc_path: docs/Evidence/PANGEN-CLUSTER-001-Evidence.md
sources:
  - docs/Evidence/PANGEN-CLUSTER-001-Evidence.md
source_commit: cce387f34d67ce7348bec08c7cea53dfd1d5cb64
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PANGEN-CLUSTER-001

The validation-evidence artifact for test unit **PANGEN-CLUSTER-001** — **Gene Clustering**,
the greedy-incremental grouping of genes into homolog/ortholog families by **sequence identity**
following the **CD-HIT** model (`ClusterGenes`). This is a **pan-genome family** Evidence file and
one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the algorithm, its metric, invariants, worked oracles, and deviations are summarized in
[[pan-genome-gene-clustering]]. Its comparative-genomics relatives are the pairwise
[[ortholog-detection-reciprocal-best-hits]] and the pan-genome pipeline
[[genome-comparison-core-dispensable]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **CD-HIT User's Guide** (Li & Godzik, authority 3) — `-c` sequence-identity threshold
    (default 0.9); **default global identity** (`-G 1`) = identical residues in the alignment /
    **full length of the shorter sequence**; **local identity** (`-G 0`, non-default) = identical
    residues / alignment length; each cluster has **one representative**, and the `.clstr` `%` is
    each member's identity **to the representative** (not all-pairs).
  - **CD-HIT Algorithm wiki** (authority 3) — the **greedy incremental** procedure: sort input
    long→short, the first sequence becomes a representative, each remaining query is compared to
    prior representatives and classified redundant-or-representative; **fast/first-match** mode
    groups a query into the *first* representative meeting the threshold (greedy, not best-hit).
  - **Roary** (Page et al. 2015, *Bioinformatics*, authority 1) — pan-genome workflow: proteins
    "iteratively pre-clustered with CD-HIT," then all-against-all BLASTP + MCL; **default 95%**
    identity for ortholog grouping; conserved gene **neighbourhood** later splits paralogs from
    true orthologs (out of scope for a greedy identity clusterer).
  - **EMBOSS needle** (authority 3) — percent-identity numerator = count of identical aligned
    positions; corroborates CD-HIT's numerator (its denominator differs by convention).
- **Corner cases / failure modes:** representative = the longest member (long→short sort);
  greedy first-match assignment is deterministic for fixed input; threshold **inclusive** (`≥`)
  so `idThreshold = 1.0` clusters only exact-identity (over the shorter length) sequences;
  empty-vs-non-empty identity behaviour (two empty → identical, one empty + one non-empty → 0).
- **Datasets (documented oracles):**
  - *Greedy clustering* — threshold 0.8, sequences `Q2 = A×12`, `R = A×10`, `Q1 = A×9+T`,
    `Q3 = C×10`: Q2 (longest) → representative, R (1.0) and Q1 (0.9) join it, Q3 (0.0) starts a new
    cluster → **2 clusters** `{Q2,R,Q1}` + `{Q3}`.
  - *Global identity worked values* (shorter-length denominator): `ATGCATGC`/`ATGCATGC` → 1.0,
    `ATGCATGC`/`ATGCATGG` → 0.875, `ATGCATGC`/`ATGCATGCAAAA` → 1.0, `ATGC`/`CGTA` → 0.0.
- **Coverage recommendations:** MUST-test the shorter-length global identity, greedy grouping
  above/below threshold, threshold-1.0 exactness, and longest-member representative + `GenomeCount`
  = distinct genomes; SHOULD-test empty/null inputs and singleton `AverageIdentity = 1.0`;
  COULD-test determinism across repeated calls.

## Deviations and assumptions

Two records, both source-backed. (1) **Ungapped alignment (deviation, algorithm doc §5.3/§5.4)** —
identity is computed over an ungapped positional comparison of the shared prefix / shorter length;
this matches CD-HIT exactly for substitution- or length-only differences (no internal indels) but
may **underestimate** identity when internal gaps are needed. (2) **Paralog splitting not performed
(assumption)** — `ClusterGenes` produces homolog groups only (no gene-neighbourhood/synteny step),
matching CD-HIT pre-clustering rather than the full Roary ortholog pipeline. No source
contradictions.
