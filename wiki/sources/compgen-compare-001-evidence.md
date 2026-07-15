---
type: source
title: "Evidence: COMPGEN-COMPARE-001 (End-to-end genome comparison — core/dispensable partition + syntenic fraction)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
source_commit: 8af87d1e0136942ca38a8ae9156e8ee5fee35080
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-COMPARE-001

The validation-evidence artifact for test unit **COMPGEN-COMPARE-001** — the end-to-end
two-genome comparison pipeline `CompareGenomes`, which partitions each genome's genes into a
**core (conserved)** set and a **dispensable (genome-specific)** set and reports an overall
**syntenic-gene fraction** `OverallSynteny`. This is a **Comparative-genomics** family Evidence
file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its outputs,
invariants, worked oracles, and corner cases are summarized in the dedicated concept
[[genome-comparison-core-dispensable]]. This is an **orchestrating pipeline** built on the
comparative-genomics sub-units it composes: reciprocal-best-hit ortholog detection
([[ortholog-detection-reciprocal-best-hits]], COMPGEN-ORTHO-001) and MCScanX synteny
(COMPGEN-SYNTENY-001, summarized in
[[synteny-and-rearrangement-detection]]). Its sibling COMPGEN units are
[[average-nucleotide-identity]] and [[conserved-gene-clusters-common-intervals]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Tettelin et al. 2005** (PNAS 102(39):13950–13955, PMC1216834; authority rank 1 — the paper
    that coined "pan-genome") — the **core/dispensable** model: the pan-genome is a *"core genome
    containing genes present in all strains"* plus *"a dispensable genome composed of genes absent
    from one or more strains and genes that are unique to each strain."* A gene present in every
    compared genome is **core/conserved**; a gene missing from one or more (in the pairwise case:
    present in only one) is **dispensable/genome-specific**. Presence is decided by an alignment
    gate: *"conserved if ... a minimum of 50% sequence conservation over 50% of the protein/gene
    length."*
  - **Moreno-Hagelsieb & Latimer 2008** (Bioinformatics 24(3):319–324) + **Tatusov et al. 1997**
    (rank 1) — the operational **shared-gene = reciprocal best hit** criterion: two genes are
    orthologs iff *"their protein products find each other as the best hit in the opposite
    genome."* In `CompareGenomes` the conserved/core genes are exactly the RBH ortholog pairs
    (the validated sub-unit COMPGEN-RBH-001); everything else in each genome is genome-specific.
    Coverage gate: *"coverage of at least 50% of any of the protein sequences."*
  - **Synteny — fraction of syntenic genes** (ScienceDirect Synteny overview + Wikipedia, rank 4,
    citing primaries) — a syntenic block is a run of collinear orthologs; *"The fraction of
    syntenic genes is a metric used to measure synteny conservation."* `CompareGenomes` reports
    `OverallSynteny = (genes in syntenic blocks) / min(|genome1|, |genome2|)`, clamped to ≤ 1.0.
  - **MCScanX** (Wang et al. 2012, NAR 40(7):e49) — the collinearity model behind the blocks
    (validated sub-unit COMPGEN-SYNTENY-001): non-overlapping chains scoring ≥ 250 (≥ 5 collinear
    anchors at MatchScore 50).
- **Algorithm behaviour (from the artifact):** partition = RBH pairs are core (conserved); the
  rest of each genome is that genome's dispensable/specific set. Outputs: `Conserved` count,
  `Specific1`, `Specific2`, and `OverallSynteny`. The partition is symmetric under swapping the
  two genomes (Specific1 ↔ Specific2, same Conserved).
- **Datasets (documented oracles):**
  - *One shared, one unique each* — g1 {a1=S, b1=U1}, g2 {c2=S, d2=U2}, pair a1↔c2 →
    Conserved 1, Specific1 1, Specific2 1.
  - *Disjoint content* — no ortholog pairs → Conserved 0, Specific1 2, Specific2 2 (all
    dispensable).
  - *Identical content (5 collinear + 1 unique each)* — 5 shared S₀…S₄ + one unique per genome →
    Conserved 5, Specific 1/1, `OverallSynteny = 5 / min(6,6) = 0.8333…`, 0 rearrangements
    (identity permutation).

## Deviations and assumptions

Two **ASSUMPTIONs**, both inherited/source-backed, neither a partition-logic gap:

1. **Alignment-free similarity replaces the Tettelin 50%/50% alignment gate.** Conserved genes are
   found by RBH where similarity is 5-mer-content Jaccard (identity ≥ 0.3) with k-mer coverage
   ≥ 0.5, not Needleman–Wunsch/BLAST. This maps Tettelin's "≥ 50% conservation over ≥ 50% length"
   and Moreno-Hagelsieb's "≥ 50% coverage, E ≤ 1e-6" onto alignment-free space (inherited verbatim
   from COMPGEN-RBH-001, Assumption 1). It does not change the partition logic tested here (core =
   reciprocal pairs, specific = the rest) — identical sequences pass the gate, disjoint sequences
   fail it, which is all the partition tests rely on.
2. **Minimum syntenic block size = 5 collinear anchors for `OverallSynteny`.** `OverallSynteny` is
   the fraction of genes inside MCScanX syntenic blocks, and MCScanX reports only chains scoring
   ≥ 250 (≥ 5 collinear anchors, the MCScanX default; Wang et al. 2012, validated in
   COMPGEN-SYNTENY-001). Hence `OverallSynteny` can be 0 even with a few conserved orthologs — a
   documented boundary, not a bug. Corner case: empty genomes → Conserved 0, Specific1 0,
   Specific2 0, OverallSynteny 0.

No contradictions among sources — Tettelin (core/dispensable definitions), Moreno-Hagelsieb/Tatusov
(RBH operationalisation), and the synteny sources (fraction-of-syntenic-genes metric + MCScanX block
threshold) are mutually consistent and each governs a distinct output of the pipeline.
