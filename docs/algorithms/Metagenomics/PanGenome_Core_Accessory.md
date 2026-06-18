# Pan-Genome Construction (Core / Accessory / Unique)

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics / PanGenome |
| Test Unit ID | PANGEN-CORE-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Pan-genome construction partitions the gene-cluster (ortholog group) repertoire of a set of genomes into the **core** genome (gene families present in essentially all genomes), the **accessory / dispensable** genome (present in some but not all), and the **unique / strain-specific** genome (present in exactly one genome) [1][2]. It additionally summarises gene-level diversity with **genome fluidity** [3] and classifies the pan-genome as **open** or **closed** using the Heaps'-law decay exponent of newly observed gene clusters per added genome [2][6]. The partitioning is exact given a cluster occupancy table; the upstream clustering and the openness fit are heuristic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The pan-genome of a clade is the union of all gene families across its genomes; the core genome is their intersection [1][2]. Because rare gene families and rare genomes are hard to estimate, gene-level similarity is also summarised by genome fluidity, a metric robust to sampling [3]. Whether continued sequencing keeps revealing new genes (open) or saturates (closed) is assessed with Heaps' law [2][6].

### 2.2 Core Model

Let there be `N` genomes and a set of gene-family clusters, each cluster having an *occupancy* = number of distinct genomes containing it.

- **Core:** present in at least `coreFraction` of the genomes, i.e. `occupancy / N ≥ coreFraction` (Roary `-cd 99`: "a gene being in at least 99% of samples" [4]). This is a fractional (percentage) test — **not** `floor(coreFraction · N)`: e.g. a 2-of-3 (66.7%) cluster is *not* core under a 0.99 threshold.
- **Unique (strain-specific / cloud):** occupancy = 1 [1][2].
- **Accessory (dispensable / shell):** all remaining clusters [2].

**Genome fluidity** [3]:

```
φ = ( 2 / (N(N−1)) ) · Σ_{k<l} (U_k + U_l) / (M_k + M_l)
```

where `U_k`, `U_l` are the numbers of gene families found only in genome `k` and only in genome `l`, and `M_k`, `M_l` are the total numbers of gene families in `k` and `l` [3]. The per-pair term `(U_k+U_l)/(M_k+M_l)` is the symmetric difference over the union (by size) of the two genomes' family sets.

**Open vs closed (Heaps' law)** [2][6]: the number of *new* gene clusters contributed by the `k`-th genome follows `n_new(k) = K · k^(−α)`. The pan-genome is **open** when `α < 1` and **closed** when `α > 1` [2][6].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Gene families (clusters) are well defined and stable across genomes | Mis-clustering shifts occupancy and therefore the core/accessory/unique partition |
| ASM-02 | The new-gene curve follows a power law in genome index | The α estimate (and open/closed call) is unreliable for non-power-law accumulation |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | core, accessory, unique are pairwise disjoint and partition all clusters | each cluster is classified by exactly one branch on its occupancy |
| INV-02 | CoreGeneCount + AccessoryGeneCount + UniqueGeneCount = TotalGenes | direct consequence of INV-01 |
| INV-03 | core ⟺ occupancy / N ≥ coreFraction (present in ≥ coreFraction of genomes); unique ⟺ occupancy = 1 | Tettelin/Roary definition [1][2][4] |
| INV-04 | 0 ≤ GenomeFluidity ≤ 1 | each pair term ∈ [0,1]; φ is their convex average [3] |
| INV-05 | identical gene content ⇒ φ = 0; pairwise-disjoint ⇒ φ = 1 | U_k = 0 for all pairs (resp. U_k + U_l = M_k + M_l) [3] |
| INV-06 | CoreFraction = CoreGeneCount / TotalGenes (0 if TotalGenes = 0) | definitional |
| INV-07 | open ⟺ Heaps decay exponent α < 1 (N ≥ 3, else Closed) | Tettelin/micropan [2][6] |

### 2.5 Comparison with Related Methods

| Aspect | Genome fluidity [3] | Pan/core size curves [2] |
|--------|---------------------|--------------------------|
| Sensitivity to rare genes/genomes | robust | sensitive |
| Output | single 0..1 dissimilarity | growth curves + open/closed call |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| genomes | `IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>` | required | genome id → its genes (id + sequence) | null or empty → empty result |
| identityThreshold | double | 0.9 | k-mer Jaccard identity for clustering (delegated to `ClusterGenes`) | 0..1 |
| coreFraction | double | 0.99 | fraction of genomes a cluster must occupy to be core | 0..1; Roary default 0.99 [4] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CoreGenes / AccessoryGenes / UniqueGenes | `IReadOnlyList<string>` | cluster IDs per partition |
| Statistics.CoreFraction | double | CoreGeneCount / TotalGenes |
| Statistics.GenomeFluidity | double | Kislyuk fluidity φ ∈ [0,1] |
| Statistics.Type | `PanGenomeType` | Open or Closed |

### 3.3 Preconditions and Validation

`null` or empty `genomes` returns an all-empty `PanGenomeResult` (no exception). Sequences are compared by the k-mer Jaccard heuristic in `ClusterGenes`; case and alphabet are taken as-is. With `N < 3` the open/closed exponent is not estimable and the result is `Closed`. With `N < 2` there are no genome pairs and fluidity is 0.

## 4. Algorithm

### 4.1 High-Level Steps

1. Cluster all genes into ortholog groups (`ClusterGenes`), giving each cluster an occupancy (distinct genome count).
2. For each cluster: core if `occupancy / N ≥ coreFraction` (present in ≥ coreFraction of genomes), else unique if occupancy = 1, else accessory.
3. Compute genome fluidity over all genome pairs from each genome's cluster-ID set [3].
4. Estimate the Heaps' law decay exponent α of new clusters per added genome and classify Open (α<1) / Closed (α≥1) [2][6].

### 4.2 Decision Rules, Scoring, Reference Tables

- coreFraction default 0.99; clustering identity default 0.9; Roary uses BLASTP identity default 95% and core 99% [4].
- α threshold = 1.0 (open below, closed at/above); minimum genomes for the fit = 3 [2][6].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ConstructPanGenome | O(g²·s) | O(g) | g = total genes, s = sequence length (all-vs-all clustering dominates) |
| Genome fluidity | O(N²·C) | O(N·C) | N genomes, C clusters; pairwise set differences |
| Heaps α fit | O(N·C) | O(C) | one accumulation pass + log-log regression |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PanGenomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs)

- `PanGenomeAnalyzer.ConstructPanGenome(genomes, identityThreshold, coreFraction)`: clusters genes, partitions into core/accessory/unique, computes fluidity and open/closed type.
- `PanGenomeAnalyzer.GetCoreGeneClusters(clusters, totalGenomes, threshold)`: core-gene identification (the Registry `IdentifyCoreGenes` referent) — filters clusters with occupancy ≥ floor(threshold·totalGenomes).
- `PanGenomeAnalyzer.CalculateGenomeFluidity` (private): Kislyuk φ.
- `PanGenomeAnalyzer.DeterminePanGenomeType` / `EstimateHeapsDecayExponent` (private): Heaps' law openness.

### 5.2 Current Behavior

Clustering uses an in-repo k-mer (k=7) Jaccard similarity, not BLAST. The repository **suffix tree was evaluated and not used**: this unit performs set-occupancy counting and arithmetic over already-formed clusters, not substring/occurrence search, so the suffix tree does not apply (exact-match occurrence enumeration is irrelevant here). The open/closed exponent is fit on a single genome ordering (dictionary order) rather than averaging over random permutations; zero-novelty steps are floored to 1 new cluster to keep the log defined.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Core = occupancy / N ≥ coreFraction (present in ≥ coreFraction of genomes); unique = occupancy 1; accessory = remainder [1][2][4].
- Genome fluidity φ = (2/(N(N−1)))·Σ_{k<l}(U_k+U_l)/(M_k+M_l) [3].
- Open ⟺ Heaps decay exponent α < 1 [2][6].

**Intentionally simplified:**

- Gene clustering: k-mer Jaccard heuristic instead of BLAST all-vs-all; **consequence:** cluster boundaries (and thus occupancy) may differ from a BLAST-based pipeline for divergent homologs.
- Heaps α fit: single dictionary-order accumulation, not averaged over random permutations; **consequence:** the α estimate (and rare borderline open/closed calls) depends on genome order.

**Not implemented:**

- Fluidity jackknife variance σ² [3]; **users should rely on:** no current alternative in-repo (variance is not returned).
- Heaps κ/γ size-curve fitting for size extrapolation; **users should rely on:** `PanGenomeAnalyzer.FitHeapsLaw` (separate method, separate unit).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Open/closed previously used unsourced `uniqueFraction > 0.1` heuristic | Deviation | wrong classification vs literature | fixed | replaced with Heaps decay-exponent criterion (α<1 ⇒ open) [2][6] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty genomes | empty result, all stats 0, Type Closed | no data |
| single genome (N=1) | all clusters core (occupancy = N), fluidity 0 | intersection over 1 genome; no pairs [1][3] |
| N < 3 | Type Closed | decay exponent not estimable [2][6] |
| identical gene content | fluidity 0 | no unique families [3] |
| pairwise-disjoint content | fluidity 1 | every family unique per pair [3] |

### 6.2 Limitations

Cluster quality is bounded by the k-mer Jaccard clusterer (no protein-level homology). The open/closed call is order-dependent and meant for small comparative sets, not large-scale population pan-genomics. Fluidity variance is not reported.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (fluidity):** A={c1,c2,c3}, B={c1,c2,c4}, C={c1,c5,c6}. Pairs: (A,B)=2/6, (A,C)=4/6, (B,C)=4/6. φ = (1/3)(2/6+4/6+4/6) = 10/18 ≈ 0.5555556 [3].

**API usage example:**

```csharp
var genomes = new Dictionary<string, IReadOnlyList<(string, string)>> { /* genomeId -> genes */ };
var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);
// result.Statistics.CoreGeneCount, .GenomeFluidity, .Type
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PanGenomeAnalyzer_ConstructPanGenome_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ConstructPanGenome_Tests.cs) — covers INV-01..INV-07
- Evidence: [PANGEN-CORE-001-Evidence.md](../../../docs/Evidence/PANGEN-CORE-001-Evidence.md)

## 8. References

1. Tettelin H, et al. 2005. Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *PNAS* 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102
2. Tettelin H, Riley D, Cattuto C, Medini D. 2008. Comparative genomics: the bacterial pan-genome. *Curr Opin Microbiol* 11(5):472–477. https://doi.org/10.1016/j.mib.2008.09.006
3. Kislyuk AO, Haegeman B, Bergman NH, Weitz JS. 2011. Genomic fluidity: an integrative view of gene diversity within microbial populations. *BMC Genomics* 12:32. https://doi.org/10.1186/1471-2164-12-32
4. Page AJ, et al. 2015. Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics* 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421
5. Lagesen K, et al. micropan: Microbial Pan-Genome Analysis — `heaps()` / `fluidity()`. CRAN. https://rdrr.io/cran/micropan/man/heaps.html
6. Wikipedia contributors. Pan-genome. https://en.wikipedia.org/wiki/Pan-genome (accessed 2026-06-13)
