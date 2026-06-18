# Comprehensive Genome Comparison (Core / Dispensable Partition)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-COMPARE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`CompareGenomes` performs a comprehensive pairwise comparison of two gene-annotated genomes. It partitions genes into the **core (conserved)** set — genes that have an ortholog in the other genome — and each genome's **dispensable (genome-specific)** set — genes with no ortholog, following the pan-genome model of Tettelin et al. (2005) [1]. Shared genes are found as reciprocal best hits [2][3]; gene-order conservation is summarised by `OverallSynteny`, the fraction of genes inside syntenic blocks [4][5]. The method is an aggregator: orthologs, syntenic blocks, and rearrangements come from the validated sub-methods, and this unit defines the partition and the synteny-fraction semantics. It is specification-driven (the partition is exact given the ortholog set), heuristic only in the alignment-free similarity used to detect orthologs.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Comparing two genomes asks which genes are shared and which are unique. Tettelin et al. (2005) [1] framed this through the **pan-genome**: a **core genome** of "genes present in all strains" and a **dispensable genome** of "genes absent from one or more strains and genes that are unique to each strain." For two genomes, "present in all" means present in both; a gene present in only one is that genome's dispensable/specific gene. Operationally, two genes in different genomes are taken to be the "same" gene (orthologs) when they are reciprocal best hits [2][3].

### 2.2 Core Model

Let `O` be the set of ortholog (reciprocal-best-hit) pairs between genome 1 and genome 2. Define:

- **Core (conserved) genes:** `ConservedGenes = |O|` — each pair is one conserved gene shared by both genomes [1][2].
- **Genome-specific (dispensable) genes:** `GenomeSpecificGenes1 = |{g ∈ genome1 : g has no pair in O}|`, and symmetrically for genome 2 [1].
- **Overall synteny:** `OverallSynteny = (Σ block.GeneCount) / min(|genome1|, |genome2|)`, the **fraction of syntenic genes** [4][5], clamped to `≤ 1`. Syntenic blocks are MCScanX collinear chains scoring ≥ 250 (≥ 5 collinear anchors) [4].

Because the RBH matching maps each gene to at most one partner, `ConservedGenes + GenomeSpecificGenes_i = |genome_i|` for each genome (INV-02).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `ConservedGenes == Orthologs.Count` | Conserved set is exactly the RBH ortholog pairs [1][2]. |
| INV-02 | `ConservedGenes + GenomeSpecificGenes_i == |genome_i|` for i ∈ {1,2} | RBH is a matching; every gene is core or that genome's dispensable gene [1]. |
| INV-03 | `0 ≤ OverallSynteny ≤ 1` | It is a fraction of genes, explicitly clamped to 1 [4][5]. |
| INV-04 | Swapping genome1/genome2 keeps `ConservedGenes`, swaps `GenomeSpecificGenes1 ↔ GenomeSpecificGenes2` | The RBH matching is reciprocal/symmetric [2][3]. |

### 2.5 Comparison with Related Methods

| Aspect | CompareGenomes (this) | Single sub-method |
|--------|------------------------|-------------------|
| Output | Partition + synteny + orthologs + blocks + rearrangements | One facet only |
| Conserved-gene basis | Reciprocal best hits | RBH alone gives pairs, not the partition |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| genome1Genes | IReadOnlyList&lt;Gene&gt; | required | Genome-1 genes in chromosomal order; each carries a sequence for ortholog detection | non-null |
| genome2Genes | IReadOnlyList&lt;Gene&gt; | required | Genome-2 genes in chromosomal order | non-null |
| minOrthologIdentity | double | 0.3 | Minimum RBH similarity for a conserved (shared) gene | 0–1 |
| minSyntenicBlockSize | int | 3 | Minimum collinear anchors per syntenic block (MCScanX score ≥ 250 still applies) | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| SyntenicBlocks | IReadOnlyList&lt;SyntenicBlock&gt; | Collinear blocks from MCScanX chaining |
| Orthologs | IReadOnlyList&lt;OrthologPair&gt; | Reciprocal-best-hit pairs (the conserved genes) |
| Rearrangements | IReadOnlyList&lt;RearrangementEvent&gt; | Breakpoint-derived rearrangements |
| OverallSynteny | double | Fraction of syntenic genes, 0–1 |
| ConservedGenes | int | Core gene count (= Orthologs.Count) |
| GenomeSpecificGenes1 | int | Genome-1 dispensable gene count |
| GenomeSpecificGenes2 | int | Genome-2 dispensable gene count |

### 3.3 Preconditions and Validation

Both gene lists must be non-null (`ArgumentNullException` otherwise). Empty lists are valid and yield an all-zero result with empty collections. Ortholog detection ignores genes with null/empty sequence. Sequence comparison is case-insensitive (uppercased internally). Coordinates are the genes' own 0-based positions.

## 4. Algorithm

### 4.1 High-Level Steps

1. Find reciprocal-best-hit orthologs between the two genomes (the conserved/core genes).
2. Build the ortholog map and find syntenic blocks (MCScanX collinear chains).
3. Detect rearrangements (breakpoints of the signed gene-order permutation).
4. Partition: conserved = ortholog count; genome-specific = genes of each genome not in any pair.
5. Compute `OverallSynteny` = genes-in-blocks ÷ smaller genome size, clamped to 1.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Conserved-gene gate: RBH similarity ≥ `minOrthologIdentity` and k-mer coverage ≥ 0.5 (maps Tettelin's 50%/50% gate [1] and Moreno-Hagelsieb's ≥50% coverage gate [2]).
- Syntenic-block report rule: MCScanX score ≥ 250 (≥ 5 collinear anchors at MatchScore 50) [4] — so `OverallSynteny` is non-zero only with ≥ 5 collinear orthologs.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CompareGenomes | O(n·m·L) | O(n+m) | n,m = gene counts; L = per-pair k-mer similarity cost; dominated by all-vs-all RBH similarity. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.CompareGenomes(...)`: aggregator producing the core/dispensable partition, OverallSynteny, orthologs, blocks, and rearrangements.
- Delegates to `FindReciprocalBestHits`, `FindSyntenicBlocks`, `DetectRearrangements`.

### 5.2 Current Behavior

The conserved set is the RBH matching; `ConservedGenes == Orthologs.Count` by construction. `GenomeSpecificGenes_i` count genes of genome i absent from the matching. `OverallSynteny` sums `GeneCount` over reported syntenic blocks and divides by the smaller genome size, clamped to 1. A search/matching reuse note: ortholog detection uses alignment-free k-mer Jaccard (not the suffix tree) because it is a scoring-based all-vs-all similarity, not exact substring location; the syntenic-block and ANI/dot-plot siblings already reuse the repository suffix tree where exact matching applies.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Core genome = genes present in both genomes; dispensable = genes unique to one genome [1].
- Conserved genes = reciprocal best hits [2][3].
- OverallSynteny = fraction of syntenic genes [4][5].

**Intentionally simplified:**

- Ortholog similarity: alignment-free 5-mer Jaccard (id ≥ 0.3, coverage ≥ 0.5) instead of a BLAST/Needleman–Wunsch alignment; **consequence:** the 50%/50% conservation gate [1][2] is approximated, but identical sequences pass and disjoint sequences fail, so the partition is unaffected for clear cases.

**Not implemented:**

- Multi-genome (≥3) pan-genome accumulation curves; **users should rely on:** dedicated pan-genome tools (Roary, PEPPAN); this method is strictly pairwise.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Alignment-free ortholog similarity | Assumption | Approximates the 50%/50% gate | accepted | Inherited from COMPGEN-RBH-001 (Assumption 1). |
| 2 | OverallSynteny needs ≥5 collinear anchors | Assumption | Synteny can read 0 with few conserved orthologs | accepted | MCScanX default score ≥ 250 [4]. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Both genomes empty | All counts 0, empty collections, OverallSynteny 0 | No ortholog pairs possible. |
| No shared genes | ConservedGenes 0; all genes genome-specific | "genes unique to each strain" [1]. |
| All genes shared | GenomeSpecific 0 for both | All core [1]. |
| < 5 collinear orthologs | Conserved > 0 but OverallSynteny 0 | MCScanX block threshold [4]. |
| null genome list | ArgumentNullException | Contract. |

### 6.2 Limitations

Pairwise only (no pan-genome growth curves). Ortholog detection is alignment-free, so highly diverged but truly orthologous genes below the k-mer gate may be miscounted as genome-specific. `OverallSynteny` reflects the MCScanX reporting threshold rather than every collinear pair.

## 7. Examples and Related Material

### 7.1 Worked Example

```csharp
var g1 = new List<ComparativeGenomics.Gene> {
    new("a1", "G1", 0, 60, '+', sharedSeq),
    new("b1", "G1", 100, 160, '+', uniqueSeq1),
};
var g2 = new List<ComparativeGenomics.Gene> {
    new("c2", "G2", 0, 60, '+', sharedSeq),
    new("d2", "G2", 100, 160, '+', uniqueSeq2),
};
var r = ComparativeGenomics.CompareGenomes(g1, g2);
// r.ConservedGenes == 1, r.GenomeSpecificGenes1 == 1, r.GenomeSpecificGenes2 == 1
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_CompareGenomes_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CompareGenomes_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [COMPGEN-COMPARE-001-Evidence.md](../../../docs/Evidence/COMPGEN-COMPARE-001-Evidence.md)
- Related algorithms: [Reciprocal_Best_Hits](./Reciprocal_Best_Hits.md), [Synteny_Block_Detection](./Synteny_Block_Detection.md), [Genome_Rearrangement_Detection](./Genome_Rearrangement_Detection.md)

## 8. References

1. Tettelin H, Masignani V, Cieslewicz MJ, et al. 2005. Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: Implications for the microbial "pan-genome". *PNAS* 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102
2. Moreno-Hagelsieb G, Latimer K. 2008. Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585
3. Tatusov RL, Koonin EV, Lipman DJ. 1997. A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631
4. Wang Y, Tang H, DeBarry JD, et al. 2012. MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research* 40(7):e49. https://doi.org/10.1093/nar/gkr1293
5. Synteny — an overview. ScienceDirect Topics. https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/synteny (accessed 2026-06-14); Wikipedia, "Synteny", https://en.wikipedia.org/wiki/Synteny (accessed 2026-06-14).
