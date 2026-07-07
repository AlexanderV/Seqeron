# Conserved Gene Clusters (Common Intervals)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-CLUSTER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Conserved gene clusters are groups of genes that occur together (as a contiguous block) in several genomes, a strong signal of shared regulation or functional association. This algorithm formalises a cluster as a **common interval** of the genomes' ortholog-group orderings: a set of ortholog-group labels that forms a contiguous window in *every* genome [1][2]. The result is exact (not heuristic): a set is reported iff it is an interval of all input genomes.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A genome is read as the sequence of its genes in chromosomal order; each gene is mapped to an ortholog-group label so that homologous genes across genomes share a label. Genomes then become permutations (or, with paralogs, sequences) over the alphabet of group labels. Genes occurring in neighbouring locations in several genomes tend to be functionally related, and such clusters have highly conserved *content* even when internal *order* differs [1].

### 2.2 Core Model

An **interval** `[i,j]` of a permutation `Pk` (with `1 ≤ i < j ≤ n`) is the set of all elements located between positions `i` and `j` inclusive [1, §2]. A **common interval** of a family `P = {P1,…,PK}` is "a set of integers that is an interval of each `Pk`, `k ∈ [K]`" [1, Definition 1, citing Uno & Yagiura [3]]. Equivalently, a set of group labels is a conserved cluster iff in every genome some contiguous window contains exactly that set of labels (and nothing else) [1][2].

For sequences with repeated labels (paralogs/duplications), an interval is `I = Set(T[i..j])` for some window, and a common interval is a set that is an interval of every sequence — any matching window ("location") suffices [2, §2, Example 1].

Worked Example 1 [1]: for `P1 = Id7 = (1 2 3 4 5 6 7)` and `P2 = (7 2 1 3 6 4 5)`, the common intervals are `{1,2}, {1,2,3}, {3,4,5,6}, {4,5}, {4,5,6}, {1..6}, {1..7}`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported cluster is a label set contiguous in **every** genome | Common-interval definition [1, Def. 1] |
| INV-02 | Every cluster has size ≥ `minClusterSize` and ≥ 2 | Interval defined only for `i < j` [1, §2] |
| INV-03 | A set contiguous in some but not all genomes is excluded | "interval of *each* Pk" [1, Def. 1] |
| INV-04 | Fewer than 2 genomes ⇒ empty result | Common interval is a family (K ≥ 2) notion [1, §2] |
| INV-05 | Output is deterministic and order-independent | Enumeration over well-defined sets [3] |

### 2.5 Comparison with Related Methods

| Aspect | Common intervals (this) | Gene teams [Bergeron et al. 2002] |
|--------|-------------------------|-----------------------------------|
| Gaps inside a cluster | none (strict contiguity) | up to δ intervening genes allowed |
| Output | exact set of conserved clusters | maximal δ-teams |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genomes` | `IReadOnlyList<IReadOnlyList<Gene>>` | required | Genomes, each a gene list in chromosomal order | non-null |
| `orthologGroups` | `IReadOnlyDictionary<string,string>` | required | gene id → ortholog-group id | non-null |
| `minClusterSize` | `int` | 3 | Min distinct groups per cluster | raised to ≥ 2 internally |
| `maxGap` | `int` | 2 | Retained for API/MCP compatibility; not used in the strict model | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IEnumerable<IReadOnlyList<string>>` | Conserved clusters; each is the sorted list of its ortholog-group labels; clusters ordered by size then lexicographically |

### 3.3 Preconditions and Validation

`genomes` and `orthologGroups` must be non-null (else `ArgumentNullException`). Genes whose id is absent from `orthologGroups` are treated as window-breaking non-members. With fewer than two genomes the result is empty. Group labels are compared ordinally (case-sensitive). `maxGap` is ignored (Evidence Assumption 1).

## 4. Algorithm

### 4.1 High-Level Steps

1. Map each genome to its ordered sequence of ortholog-group labels (missing genes → non-member sentinel).
2. Enumerate every contiguous window of genome 0 with no sentinel; collect each window's distinct label set of size ≥ `minClusterSize` as a candidate.
3. For each candidate set, test whether it is an interval of *every* genome (some window's set equals it exactly).
4. Report the sets common to all genomes, sorted deterministically.

### 4.2 Decision Rules

A window may not contain a foreign group: an interval is the set of *all* elements of the window [1, §2], so any non-member label encountered while extending a window terminates that location.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindConservedClusters | O(n² · K · m) | O(n²) | n = genes per genome, K = genomes, m = avg cluster size; candidates from genome-0 windows verified against all genomes. The simple quadratic check (cf. Uno & Yagiura's O(n²) LHP [3]) is adequate for the small gene-cluster inputs in scope; the O(n+K) RC algorithm [3] is not needed. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.FindConservedClusters(genomes, orthologGroups, minClusterSize, maxGap)`: returns the common-interval clusters.
- `ComparativeGenomics.IsIntervalOf(...)` (private): tests whether a label set is an interval of one genome.

### 5.2 Current Behavior

Candidates are generated only from genome 0 (every common interval is necessarily an interval of genome 0, so this is complete). Cluster identity is the *set* of labels (deduplicated, ordinal-sorted); duplicate windows collapse to one cluster. Output is sorted by size then joined-label order for determinism. The repository suffix tree was evaluated and **not** used: the operation is exact set-membership over contiguous windows of small label sequences, not substring occurrence search, so a suffix tree does not fit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Common interval = a set that is an interval of every genome [1, Def. 1; 3].
- Interval = set of all elements of a contiguous window; foreign elements excluded [1, §2; 2, §2].
- Size ≥ 2 (interval requires `i < j`); whole set is a trivial common interval [1, §2].
- Sequence/duplicate handling: any matching window location suffices [2, Example 1].

**Intentionally simplified:**

- Enumeration uses the simple O(n²) quadratic scheme rather than Uno & Yagiura's O(n+K) RC algorithm; **consequence:** identical output, higher asymptotic cost (acceptable for small gene-cluster inputs) [3].

**Not implemented:**

- Gene-teams δ-gap clusters (Bergeron, Corteel & Raffinot 2002); **users should rely on:** the strict common-interval model here; the gapped extension is out of scope (source not retrievable this session — Evidence Assumption 1).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `maxGap` parameter | Assumption | Parameter present but unused; clusters are strictly contiguous | accepted | API/MCP compatibility; Evidence Assumption 1 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| < 2 genomes | empty | Family notion, K ≥ 2 [1, §2] |
| gene without ortholog group | breaks windows; never a member | interval = set of all window elements [1, §2] |
| repeated group label | any matching window counts | sequence common interval [2] |
| no conserved set ≥ minClusterSize | empty | no common interval |
| null genomes / orthologGroups | `ArgumentNullException` | defensive contract |

### 6.2 Limitations

Strict-contiguity only (no gaps); does not implement gene teams. Ortholog grouping quality is the user's responsibility. Quadratic in genome length.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Two genomes, ortholog groups labelled "1".."7"; expects the common intervals of Example 1.
var clusters = ComparativeGenomics.FindConservedClusters(genomes, orthologGroups, minClusterSize: 2);
// => {1,2},{4,5},{1,2,3},{4,5,6},{3,4,5,6},{1..6},{1..7}
```

**Numerical walk-through:** P1=(1 2 3 4 5 6 7), P2=(7 2 1 3 6 4 5). The set {3,4,5,6} is positions 3–6 in P1 and positions 4,6,5,? — in P2 the labels 3,6,4,5 occupy positions 4–7, a contiguous window with exactly {3,4,5,6}, so it is a common interval. The set {2,3} is contiguous in P1 but 2 (pos 2) and 3 (pos 4) are not adjacent in P2, so it is excluded.

### 7.2 Performance Baseline

The algorithm is O(n²·K·m). Measured baseline (informal, debug build, .NET on Apple Silicon, June 2026): the full 12-test fixture — including a 3-genome × 12-gene property-based randomised case (`FindConservedClusters_RandomGenomes_*`) — runs in ~7 ms total. Inputs in scope (operon-scale clusters, tens to low hundreds of genes) are well within interactive latency; the O(n+K) RC algorithm [3] would only matter for whole-genome ortholog orderings.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_FindConservedClusters_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ComparativeGenomics_FindConservedClusters_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [COMPGEN-CLUSTER-001-Evidence.md](../../../docs/Evidence/COMPGEN-CLUSTER-001-Evidence.md)
- Related algorithms: [Reciprocal_Best_Hits](../Comparative_Genomics/Reciprocal_Best_Hits.md)

## 8. References

1. Bui-Xuan B-M, Habib M, Paul C. 2013. MinMax-Profiles: A Unifying View of Common Intervals, Nested Common Intervals and Conserved Intervals of K Permutations. arXiv:1304.5140. https://arxiv.org/abs/1304.5140
2. Didier G, Schmidt T, Stoye J, Tsur D. 2013. Extending Common Intervals Searching from Permutations to Sequences. arXiv:1310.4290. https://arxiv.org/abs/1310.4290
3. Uno T, Yagiura M. 2000. Fast Algorithms to Enumerate All Common Intervals of Two Permutations. Algorithmica 26(2):290–309. https://doi.org/10.1007/s004539910014
4. Heber S, Stoye J. 2001. Finding All Common Intervals of k Permutations. CPM 2001, LNCS 2089:207–218. https://doi.org/10.1007/3-540-48194-X_19
