---
type: concept
title: "Genome comparison pipeline (core/dispensable partition + syntenic fraction)"
tags: [comparative-genomics, algorithm]
sources:
  - docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
  - docs/algorithms/Comparative_Genomics/Genome_Comparison.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-compare-001-evidence
      evidence: "Test Unit ID: COMPGEN-COMPARE-001 ... Algorithm: Comprehensive two-genome comparison — core/dispensable gene partition and the syntenic-gene fraction"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:synteny-and-rearrangement-detection
      source: compgen-compare-001-evidence
      evidence: "CompareGenomes reports OverallSynteny = (genes in syntenic blocks)/min(|g1|,|g2|); the blocks come from the MCScanX collinearity model, the validated sub-unit COMPGEN-SYNTENY-001"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:average-nucleotide-identity
      source: compgen-compare-001-evidence
      evidence: "COMPGEN-COMPARE-001 is a sibling end-to-end Comparative-genomics pipeline of COMPGEN-ANI-001; both quantify two-genome similarity — COMPARE at the gene core/dispensable level, ANI at nucleotide identity"
      confidence: medium
      status: current
---

# Genome comparison pipeline (core/dispensable partition + syntenic fraction)

`CompareGenomes` is the **end-to-end two-genome comparison pipeline** of the
**Comparative-genomics** family (`COMPGEN-*`). It partitions each genome's genes into a **core
(conserved)** set and a **dispensable (genome-specific)** set — the pairwise case of the
**pan-genome** model of Tettelin et al. 2005 — and reports an overall **syntenic-gene fraction**
`OverallSynteny`. Unlike the family's single-metric units it is an **orchestrating pipeline**: it
composes reciprocal-best-hit ortholog detection
([[ortholog-detection-reciprocal-best-hits]], COMPGEN-ORTHO-001) and MCScanX synteny
([[synteny-and-rearrangement-detection]], COMPGEN-SYNTENY-001). Its single-metric siblings are
[[average-nucleotide-identity]] (nucleotide identity of the two genomes) and
[[conserved-gene-clusters-common-intervals]] (order-free contiguous gene sets). Validated under
test unit **COMPGEN-COMPARE-001**; the validation record is [[compgen-compare-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
artifact pattern.

## The core/dispensable (pan-genome) partition

Tettelin et al. 2005 (PNAS, the paper that coined "pan-genome") split a set of genomes into a
**core genome** ("genes present in all strains") and a **dispensable genome** ("genes absent from
one or more strains and genes that are unique to each strain"). In the two-genome case this pipeline
implements:

- **Core / conserved** = genes present in *both* genomes = the **reciprocal-best-hit (RBH)
  ortholog pairs** ([[ortholog-detection-reciprocal-best-hits]]; Moreno-Hagelsieb & Latimer 2008;
  Tatusov et al. 1997): *"two genes ... are
  deemed orthologs if their protein products find each other as the best hit in the opposite
  genome."*
- **Dispensable / genome-specific** = everything else — a gene with no reciprocated best hit in
  the other genome. `Specific1` and `Specific2` are the per-genome dispensable counts.

Presence/conservation in the source is decided by an alignment gate (Tettelin: "≥ 50% conservation
over ≥ 50% length"; Moreno-Hagelsieb: "≥ 50% coverage, E ≤ 1e-6"), operationalised here in
alignment-free space (see Assumptions).

## OverallSynteny — the syntenic-gene fraction

Beyond content, the pipeline scores **gene-order conservation** as the fraction of genes falling in
syntenic (collinear) blocks (ScienceDirect Synteny overview: "The fraction of syntenic genes is a
metric used to measure synteny conservation"):

```
OverallSynteny = (genes in syntenic blocks) / min(|genome1|, |genome2|),  clamped to <= 1.0
```

The syntenic blocks come from the **MCScanX** collinearity model (Wang et al. 2012, the validated
COMPGEN-SYNTENY-001 sub-unit): non-overlapping chains scoring ≥ 250 (≥ 5 collinear anchors at
MatchScore 50). Because of that threshold `OverallSynteny` can be **0 even when conserved orthologs
exist** — fewer than 5 collinear anchors form no reported block.

## Outputs and invariants

| Output | Meaning |
|--------|---------|
| `Conserved` | count of core genes = number of RBH ortholog pairs |
| `Specific1` | genome-1 dispensable (genome-specific) gene count |
| `Specific2` | genome-2 dispensable gene count |
| `OverallSynteny` | fraction of genes in syntenic blocks, in [0, 1] |

- **Symmetry** — swapping genome1/genome2 swaps `Specific1`/`Specific2` and leaves `Conserved`
  unchanged (the RBH matching is symmetric).
- **`OverallSynteny` clamped to ≤ 1.0** (documented clamp for degenerate inputs).
- Every gene of each genome is *either* core *or* that genome's specific set (a partition).

## Documented oracles

- **One shared, one unique each** — g1 {a1=S, b1=U1}, g2 {c2=S, d2=U2}, ortholog a1↔c2 →
  `Conserved 1`, `Specific1 1`, `Specific2 1`.
- **Disjoint content** — no ortholog pairs → `Conserved 0`, `Specific1 2`, `Specific2 2` (all
  dispensable).
- **Identical content (5 collinear + 1 unique each)** — 5 shared S₀…S₄ plus one unique gene per
  genome → `Conserved 5`, `Specific 1/1`, `OverallSynteny = 5/min(6,6) = 0.8333…`, 0 rearrangements
  (identity permutation).

## Edge cases

- **All genes shared** → dispensable count 0 for both genomes (everything is core).
- **No genes shared** → `Conserved 0`, every gene is that genome's specific gene.
- **Empty genomes** → `Conserved 0`, `Specific1 0`, `Specific2 0`, `OverallSynteny 0` (no ortholog
  pairs possible).
- **One-directional best hit** → not reciprocated → not core, stays genome-specific.
- **Orthologs detected but < 5 collinear** → `OverallSynteny 0` while `Conserved > 0` (MCScanX
  block threshold, Assumption 2).

## Assumptions (source-backed, not correctness gaps)

1. **Alignment-free RBH similarity** — the conserved set uses 5-mer-content Jaccard (identity ≥ 0.3,
   k-mer coverage ≥ 0.5) in place of the Tettelin 50%/50% alignment gate (inherited from
   COMPGEN-RBH-001). The partition logic (core = reciprocal pairs, specific = the rest) is
   unchanged: identical sequences pass, disjoint sequences fail.
2. **Minimum syntenic block = 5 collinear anchors** for `OverallSynteny` (the MCScanX default;
   Wang et al. 2012). Hence `OverallSynteny` can be 0 with a few conserved orthologs.

## Reference tools

The definitions trace to **Tettelin et al. 2005** (PNAS, pan-genome core/dispensable model),
**Moreno-Hagelsieb & Latimer 2008** + **Tatusov et al. 1997** (RBH ortholog criterion), the
**fraction-of-syntenic-genes** metric (ScienceDirect/Wikipedia synteny overviews), and **MCScanX**
(Wang et al. 2012, the collinearity block model). No deviations from the sources are recorded; the
two assumptions are the inherited alignment-free RBH gate and the MCScanX block threshold, both
source-backed.
