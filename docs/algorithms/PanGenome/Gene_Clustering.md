# Gene Clustering (Homolog Grouping by Sequence Identity)

| Field | Value |
|-------|-------|
| Algorithm Group | PanGenome |
| Test Unit ID | PANGEN-CLUSTER-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Gene clustering groups genes from multiple genomes into homolog families ("ortholog groups") so that downstream pan-genome analysis can decide which gene clusters are core, accessory, or unique. This implementation follows the CD-HIT greedy incremental clustering model [1][2]: sequences are processed from longest to shortest, the longest sequence seeds a cluster as its representative, and each subsequent sequence joins the first existing representative whose global sequence identity meets the threshold, otherwise it starts a new cluster. Membership is decided by **global sequence identity** — identical residues divided by the length of the shorter sequence [2] — not by k-mer or alignment-score heuristics. The result is a deterministic, threshold-controlled partition of all input genes into clusters; it is a heuristic clusterer, not an exact optimum.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Pan-genome pipelines (e.g., Roary [3]) first reduce coding sequences into homologous gene families before classifying them as core/accessory. Roary pre-clusters with CD-HIT and refines with BLASTP + MCL at a default 95% identity [3]. The CD-HIT step is itself a complete identity-based clusterer [1][2] and is the model implemented here.

### 2.2 Core Model

**Global sequence identity** (CD-HIT default, `-G 1`) [2]:

> "global sequence identity calculated as: number of identical amino acids in alignment divided by the full length of the shorter sequence."

For two sequences `s1`, `s2` with `m = min(|s1|, |s2|)` and `k` identical residues in the (ungapped) positional alignment:

```
identity(s1, s2) = k / m            (1.0 if both empty; 0.0 if exactly one empty)
```

**Greedy incremental clustering** [1]:

> "CD-HIT ... sorts the input sequences from long to short ... The first sequence is automatically classified as the first cluster representative sequence. Then each query sequence ... is compared to the representative sequences found before it, and is classified as redundant or representative based on whether it is similar to one of the existing representative sequences."

> "In default manner (fast mode), a query is grouped into the first representative without comparing to other representatives." [1]

A query joins the first representative `r` with `identity(query, r) >= idThreshold`; otherwise it becomes a new representative [1][2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Sequences differ only by substitutions or length (no internal indels) | Identity over the ungapped shared prefix underestimates CD-HIT's gapped identity; near-homologs requiring internal gaps may be split into separate clusters |
| ASM-02 | Homolog grouping is sufficient (no paralog separation) | Paralogous copies may share a cluster; Roary's gene-neighbourhood ortholog split [3] is not performed |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every input gene is in exactly one cluster (partition) | Each sequence is either redundant (joins one cluster) or a new representative [1] |
| INV-02 | Σ cluster sizes = number of input genes | Direct consequence of INV-01 |
| INV-03 | 0 ≤ identity ≤ 1; identical = 1.0, disjoint = 0.0 | Ratio of identical positions to shorter length [2] |
| INV-04 | Members of a cluster have identity ≥ threshold to its representative | Greedy join condition [1] |
| INV-05 | The cluster representative is the longest member | Long→short processing order [1] |
| INV-06 | A singleton cluster has AverageIdentity = 1.0 | A sequence is 100% identical to itself |

### 2.5 Comparison with Related Methods

| Aspect | This implementation (CD-HIT greedy) | Roary (full pipeline) |
|--------|-------------------------------------|------------------------|
| Similarity measure | Global identity, shorter-length denominator [2] | BLASTP %identity (default 95%) [3] |
| Clustering | Greedy first-match, representative-based [1] | CD-HIT pre-cluster + MCL [3] |
| Paralog split | No | Yes (gene neighbourhood) [3] |
| Determinism | Yes (stable long→short order) | Pipeline-dependent |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genomes` | `IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>` | required | Genome id → genes | null → empty result; null inner lists skipped |
| `identityThreshold` | `double` | `0.9` (CD-HIT `-c` default [2]) | Global-identity cutoff for joining a cluster | typically in [0,1]; inclusive (`>=`) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ClusterId` | `string` | `cluster_1`, `cluster_2`, … in emission order |
| `GeneIds` | `IReadOnlyList<string>` | Member gene ids (representative first) |
| `GenomeIds` | `IReadOnlyList<string>` | Distinct contributing genome ids |
| `GenomeCount` | `int` | Number of distinct genomes (cluster occupancy) |
| `AverageIdentity` | `double` | Mean pairwise global identity within the cluster (1.0 for singletons) |
| `ConsensusSequence` | `string` | The representative (longest) member sequence |

### 3.3 Preconditions and Validation

`null` `genomes` → empty enumeration (no throw), matching the sibling `ConstructPanGenome` contract. Null inner gene lists are skipped. Empty input → empty enumeration. Identity is case-sensitive positional comparison (no T↔U normalization). Two empty sequences are defined identical (1.0); one empty + one non-empty → 0.0.

## 4. Algorithm

### 4.1 High-Level Steps

1. Flatten all genomes into `(genomeId, geneId, sequence)` records (skip null lists).
2. Order record indices by **decreasing sequence length** (stable; ties keep input order) [1].
3. For each record in that order: compare to existing representatives; join the first whose global identity ≥ `identityThreshold`; otherwise create a new cluster with this record as representative [1].
4. For each cluster, compute distinct genome ids, mean pairwise identity (1.0 for singletons), and set the representative (longest member) as the consensus.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Identity cutoff:** inclusive `>=` `identityThreshold`; default 0.9 [2].
- **Global identity:** identical positions over `min(|s1|,|s2|)` [2].
- **Representative selection:** longest member by construction (long→short order) [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ClusterGenes` | O(g·c·s) ≤ O(g²·s) | O(g) | g = total genes, c = clusters, s = sequence length; each gene compared to current representatives only (not all-pairs), so cost is bounded by clusters, worst case all-singletons = O(g²·s) |
| `CalculateSequenceIdentity` | O(s) | O(1) | positional scan of shorter length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PanGenomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs)

- `PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold)`: greedy CD-HIT clustering into `GeneCluster` records.
- `PanGenomeAnalyzer.CalculateSequenceIdentity(seq1, seq2)` (private): global sequence identity.
- `PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters)`: gene presence/absence rows derived from clusters (delegate).

### 5.2 Current Behavior

Each query is compared against cluster **representatives** (the first/longest member), not against all members, matching CD-HIT fast mode [1] and bounding cost below the naive all-pairs O(g²). Clustering is deterministic: the long→short ordering is a stable sort, so equal-length genes retain input order. The **suffix tree was not used**: the operation is identity-scored similarity over whole gene sequences (an O(s) positional comparison per pair), not exact-substring occurrence search, so the repository suffix tree (`Contains`/`FindAllOccurrences`) does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Global sequence identity = identical residues / length of shorter sequence (CD-HIT `-G 1` default) [2].
- Greedy incremental clustering: long→short order, longest-as-representative, first-match join, new-representative fallback [1].
- Inclusive identity cutoff with default 0.9 [2].

**Intentionally simplified:**

- Identity is computed over an **ungapped** positional alignment of the shared prefix; **consequence:** for homologs separated by internal insertions/deletions the identity is underestimated relative to CD-HIT's banded gapped alignment, so such pairs may be split (ASM-01).
- No CD-HIT short-word index/filter (performance only, not output): the implementation does the exact positional comparison directly.

**Not implemented:**

- Paralog-vs-ortholog separation by gene neighbourhood; **users should rely on:** the full Roary pipeline [3] for true ortholog groups (ASM-02).
- BLASTP/MCL refinement step of Roary; **users should rely on:** external Roary for large divergent pan-genomes.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Ungapped identity | Assumption | Internal indels underestimate identity | accepted | ASM-01; exact for substitution/length-only differences |
| 2 | Homolog (not ortholog) clustering | Assumption | Paralogs may share a cluster | accepted | ASM-02 |
| 3 | API name `CreatePresenceAbsenceMatrix` vs checklist `GeneratePresenceAbsenceMatrix` | Deviation | Naming only | accepted | non-correctness-affecting; sibling-consistent name kept |
| 4 | Prior k-mer Jaccard "identity" replaced by CD-HIT global identity | Deviation (fix) | Corrected misclassification of near-identical genes | fixed | conformance to [2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null `genomes` | empty enumeration | matches `ConstructPanGenome` contract |
| empty `genomes` | empty enumeration | no genes to cluster |
| genome with empty gene list | contributes nothing | nothing to add |
| both sequences empty | identity 1.0 | 0 differences over 0 positions |
| one sequence empty | identity 0.0 | no shared residues |
| `idThreshold = 1.0` | only exact-identity-over-shorter-length sequences cluster | inclusive cutoff [2] |
| singleton cluster | AverageIdentity 1.0 | self-identity (INV-06) |

### 6.2 Limitations

Heuristic (greedy, order-dependent on length) — not an optimal clustering. No internal-gap alignment, no paralog separation, no BLASTP/MCL refinement. Identity is case-sensitive and assumes a consistent alphabet across input genes.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
{
    ["g1"] = new[] { ("a", "ATGCATGCAT"), ("b", "TTTTTTTTTT") },
    ["g2"] = new[] { ("c", "ATGCATGCAT") }, // identical to a -> same cluster
};

foreach (var cluster in PanGenomeAnalyzer.ClusterGenes(genomes, identityThreshold: 0.9))
    Console.WriteLine($"{cluster.ClusterId}: occupancy={cluster.GenomeCount}");
// cluster_1: occupancy=2  (a + c)
// cluster_2: occupancy=1  (b)
```

**Numerical walk-through:** sequences `AAAAAAAAAA` (R), `AAAAAAAAAT` (Q1, 9/10=0.9), `CCCCCCCCCC` (Q3, 0/10=0.0). At threshold 0.9: R seeds cluster 1; Q1 joins (0.9 ≥ 0.9); Q3 starts cluster 2. Two clusters.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PanGenomeAnalyzer_ClusterGenes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Metagenomics/PanGenomeAnalyzer_ClusterGenes_Tests.cs) — covers INV-01..INV-06
- Evidence: [PANGEN-CLUSTER-001-Evidence.md](../../../docs/Evidence/PANGEN-CLUSTER-001-Evidence.md)
- Related algorithms: [PanGenome_Core_Accessory](../Metagenomics/PanGenome_Core_Accessory.md)

## 8. References

1. CD-HIT Algorithm wiki. weizhongli/cdhit. https://github.com/weizhongli/cdhit/wiki/1.-Algorithm
2. CD-HIT User's Guide (Li lab). `-c` / `-G` global sequence identity definition. https://vcru.wisc.edu/simonlab/bioinformatics/programs/cd-hit/cdhit-user-guide.pdf
3. Page AJ, Cummins CA, Hunt M, et al. 2015. Roary: rapid large-scale prokaryote pan genome analysis. Bioinformatics 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421
4. Li W, Godzik A. 2006. Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics 22(13):1658–1659. https://doi.org/10.1093/bioinformatics/btl158
5. EMBOSS needle manual (percent-identity convention). https://galaxy-iuc.github.io/emboss-5.0-docs/needle.html
