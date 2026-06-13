# Phylogenetic Marker Selection

| Field | Value |
|-------|-------|
| Algorithm Group | PanGenome |
| Test Unit ID | PANGEN-MARKER-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Reference |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Phylogenetic marker selection picks, from a pan-genome's gene clusters, the genes most suitable for reconstructing strain relationships. The established practice (panX, Roary) is to keep **single-copy core genes** — clusters present in every genome with exactly one gene per genome — and to use the **variable (parsimony-informative) positions** of their alignments, because invariant columns carry no topological signal [1][2]. This implementation selects single-copy core clusters, scores each by its number of parsimony-informative sites (PIS) using the column criterion of Zvelebil & Baum (2008) [4], keeps those with at least one PIS, and returns them ranked by descending PIS. It is a deterministic, specification-driven filter; it does not build a tree.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A pan-genome partitions genes into clusters (ortholog groups). The **core genome** is the set of clusters shared by (nearly) all genomes; Roary defines core as present in **≥ 99%** of samples [2]. For phylogeny, only **single-copy core** clusters are used — panX: "gene clusters in which all strains are represented exactly once" [1]; Roary filters paralog-containing clusters out of the core gene alignment [2][3]. Within those genes, the **parsimony-informative sites** are the columns that can discriminate tree topologies [4].

### 2.2 Core Model

A multiple-sequence-alignment column is **parsimony-informative** iff it has at least two different character states and each of at least two of those states occurs in at least two sequences (Zvelebil & Baum 2008) [4]. Monomorphic columns (one state) and singleton columns (a variant in only one sequence) are excluded because they imply the same number of changes on every topology [4].

A cluster is selected as a marker iff: (a) it is single-copy core — present in all `totalGenomes` with exactly one gene per genome [1][2][3]; and (b) its member alignment has ≥ 1 parsimony-informative site (panX uses the variable positions; a fully conserved gene has none) [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ PIS ≤ alignment length | each column contributes 0 or 1 [4] |
| INV-02 | Monomorphic or singleton column contributes 0 | excluded by the definition [4] |
| INV-03 | Column with ≥ 2 states each in ≥ 2 sequences contributes exactly 1 | minimal informative pattern [4] |
| INV-04 | Every selected marker is single-copy core (all genomes, one gene each) | panX "exactly once"; Roary paralog filtering [1][2][3] |
| INV-05 | Every selected marker has ≥ 1 PIS | panX extracts variable positions [1] |
| INV-06 | Markers ordered by descending PIS; result size ≤ `maxMarkers` | most informative first [1][4] |
| INV-07 | PIS is invariant to row order and to bijective state relabeling | column-content property [4] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | This (single-copy core + PIS) | Identity-band heuristic (removed) |
|--------|-------------------------------|-----------------------------------|
| Marker definition | single-copy core gene [1][2] | clusters in an average-identity band (unsourced) |
| Informativeness measure | parsimony-informative sites [4] | consensus-sequence length (unsourced) |
| Source backing | panX, Roary, Zvelebil | none |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `alignedSequences` | `IReadOnlyList<string>` | required | rows of an alignment for PIS counting | equal length; ≥ 2 rows for any PIS |
| `genomes` | `IReadOnlyDictionary<string, IReadOnlyList<(string,string)>>` | required | genome → (gene id, sequence) | used to recover cluster member sequences |
| `coreClusters` | `IEnumerable<GeneCluster>` | required | candidate clusters | typically the core set |
| `totalGenomes` | `int` | required | total genomes in the analysis | > 0 |
| `maxMarkers` | `int` | 100 | cap on returned markers | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (PIS) | `int` | number of parsimony-informative columns |
| (markers) | `IEnumerable<GeneCluster>` | single-copy core markers, descending PIS, ≤ `maxMarkers` |

### 3.3 Preconditions and Validation

`CountParsimonyInformativeSites`: null, < 2 rows, empty rows, or rows of unequal length → 0 (no countable common alignment). Characters compared as-is (case-sensitive; no T↔U normalization). `SelectPhylogeneticMarkers`: null `genomes`/`coreClusters`, `totalGenomes` ≤ 0, or `maxMarkers` ≤ 0 → empty sequence. Clusters that are not single-copy core, or have 0 PIS, are dropped.

## 4. Algorithm

### 4.1 High-Level Steps

1. Index every gene id → sequence from `genomes`.
2. For each candidate cluster, keep only single-copy core: `GenomeCount == totalGenomes` and `GeneIds.Count == totalGenomes`.
3. Recover the cluster's member sequences and count parsimony-informative sites over them.
4. Keep clusters with PIS ≥ 1; order by descending PIS (ties by ordinal cluster id), take `maxMarkers`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Parsimony-informative column test (Zvelebil & Baum 2008 [4]): tally states in the column; informative iff ≥ 2 distinct states each occur in ≥ 2 rows.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CountParsimonyInformativeSites` | O(r · L) | O(s) | r rows, L columns, s distinct states per column |
| `SelectPhylogeneticMarkers` | O(n · g) | O(g) | n clusters, g = genomes (member alignment per cluster), matches checklist O(n × g) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PanGenomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs)

- `PanGenomeAnalyzer.CountParsimonyInformativeSites(...)`: counts parsimony-informative columns of an alignment.
- `PanGenomeAnalyzer.SelectPhylogeneticMarkers(...)`: selects single-copy core markers ranked by PIS.

### 5.2 Current Behavior

PIS is counted directly over the cluster's member sequences treated as an alignment (equal-length, ungapped, position-wise — the same convention CD-HIT global identity uses elsewhere in this class). No suffix tree is used: this is a column-wise statistic over a fixed set of equal-length rows and a per-cluster filter, not a substring/occurrence search, so the repository suffix tree is not applicable. Tie-breaking on equal PIS is by ordinal cluster id for determinism.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Single-copy core marker rule: present in all genomes with exactly one gene per genome (panX "all strains represented exactly once"; Roary paralog filtering) [1][2][3].
- Parsimony-informative column criterion: ≥ 2 states each in ≥ 2 sequences (Zvelebil & Baum 2008) [4].
- Selection of variable (≥ 1 PIS) markers only (panX extracts variable positions) [1].

**Intentionally simplified:**

- Alignment source: members are treated as already aligned (equal length); **consequence:** clusters whose members differ in length score 0 PIS and are not selected, whereas panX/Roary would first run MAFFT/PRANK. The PIS criterion itself is unchanged.

**Not implemented:**

- Tree construction (FastTree/RaxML) and codon-aware multiple alignment; **users should rely on:** external tools (panX, Roary + RAxML) for the downstream tree.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Equal-length members as alignment | Assumption | unequal-length clusters score 0 PIS | accepted | Evidence Assumption 1; no in-repo aligner |
| 2 | Old identity-band/length heuristic removed | Deviation (fix) | API signature changed; conforms to sources | fixed | unsourced thresholds were defects |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| < 2 aligned rows | PIS 0 | no state can occur in ≥ 2 rows [4] |
| Unequal-length members | PIS 0 → not selected | no common alignment (Assumption 1) |
| Fully conserved single-copy core | not selected | 0 variable positions [1] |
| Paralog (genome with ≥ 2 genes) | not selected | not single-copy [1][2][3] |
| Below-core cluster | not selected | absent from a genome [2] |
| null/empty inputs | empty result, no exception | degenerate input |

### 6.2 Limitations

Markers are read off ungapped equal-length members; indels are not modeled. Selection picks informative genes but does not guarantee a correct tree — core genes may be affected by homologous recombination [1].

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**Numerical walk-through:** alignment s1=`AAAAA`, s2=`AAACA`, s3=`AACCG`, s4=`ACCTG`. Column 3 (A,A,C,C) and column 5 (A,A,G,G) each have two states each appearing twice → parsimony-informative; columns 1 (mono), 2 (singleton C), 4 (four singletons) are not → **PIS = 2** [4].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests.cs) — covers `INV-01`–`INV-07`
- Evidence: [PANGEN-MARKER-001-Evidence.md](../../../docs/Evidence/PANGEN-MARKER-001-Evidence.md)
- Related algorithms: [Pan_Genome_Growth_Model](../PanGenome/Pan_Genome_Growth_Model.md)

## 8. References

1. Ding W, Baumdicker F, Neher RA. 2018. panX: pan-genome analysis and exploration. *Nucleic Acids Research* 46(1):e5. https://doi.org/10.1093/nar/gkx977
2. Page AJ, Cummins CA, Hunt M, et al. 2015. Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics* 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421
3. Sanger Pathogens. Roary: the pan genome pipeline (documentation). https://sanger-pathogens.github.io/Roary/
4. Zvelebil M, Baum JO. 2008. *Understanding Bioinformatics*. Garland Science (definition via https://en.wikipedia.org/wiki/Informative_site).
