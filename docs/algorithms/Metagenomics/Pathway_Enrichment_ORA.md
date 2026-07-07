# Metabolic Pathway Enrichment (Over-Representation Analysis)

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-PATHWAY-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Pathway over-representation analysis (ORA) tests whether a query gene set (e.g. predicted or
differential genes from a metagenome) contains more members of a given metabolic pathway than would be
expected by chance, given a background gene universe. For each caller-supplied pathway it returns the
right-tail hypergeometric p-value P(X ≥ overlap). This is an exact (non-heuristic, non-randomized)
probability computed from the hypergeometric distribution [1][2]. Pathway-to-gene definitions are NOT
hard-coded — KEGG/MetaCyc are large curated databases — so the caller supplies the pathway membership
and the implementation supplies the statistics; this is a *Framework* algorithm.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Functional annotation assigns genes to pathways/terms (GO, KEGG, MetaCyc). ORA then asks, per pathway,
whether the query is enriched for that pathway relative to a background distribution, ranking pathways
by significance [1].

### 2.2 Core Model

The number of query genes falling in a pathway follows a hypergeometric distribution (sampling without
replacement). With N = background size, M = pathway size, n = query size, x = overlap, the
over-representation (upper-tail) p-value is [1][2]:

```
P(X ≥ x) = 1 − Σ_{i=0}^{x−1}  C(M,i) · C(N−M, n−i) / C(N, n)
```

The hypergeometric is preferred over the binomial because it models sampling without replacement, which
is more accurate for finite gene universes [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Genes are drawn uniformly at random from the background without replacement | p-values mis-calibrated if genes are not exchangeable (e.g. length/expression bias) |
| ASM-02 | Background universe is the correct reference set | A wrong/over-broad background inflates or deflates significance |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ p-value ≤ 1 | It is a probability; result is clamped to [0, 1] [1] |
| INV-02 | x = 0 ⇒ p-value = 1 | Empty upper sum: 1 − 0 = 1 [2] |
| INV-03 | N, M, or n ≤ 0 ⇒ p-value = 1 | No success can be drawn from a degenerate population |
| INV-04 | P(X ≥ x) is invariant under swapping M ↔ n | Symmetry of the hypergeometric distribution |
| INV-05 | Pathways returned ascending by p-value | Significance ranking [1] |

### 2.5 Comparison with Related Methods

| Aspect | Hypergeometric ORA | GSEA (rank-based) |
|--------|--------------------|-------------------|
| Input | Discrete query set + background | Full ranked gene list |
| Threshold dependence | Requires a query cut-off | Threshold-free |
| Test | Exact hypergeometric tail | Enrichment-score permutation |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| queryGenes | IEnumerable&lt;string&gt; | required | Genes of interest | de-duplicated by ordinal equality |
| pathwayDatabase | IReadOnlyDictionary&lt;string, IReadOnlyCollection&lt;string&gt;&gt; | required | Pathway id → member genes (caller-supplied) | each value de-duplicated |
| backgroundGenes | IEnumerable&lt;string&gt;? | null | Background universe | if null/empty → union of pathway members |
| x, bigN, bigM, n | int | required (HypergeometricUpperTail) | overlap, N, M, n | integers |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Pathway | string | Pathway id |
| Overlap | int | Query genes in the pathway (x) |
| PathwaySize | int | Pathway members within the background (M) |
| QuerySize | int | Distinct query genes (n) |
| BackgroundSize | int | Background universe size (N) |
| PValue | double | Right-tail hypergeometric P(X ≥ Overlap) |

### 3.3 Preconditions and Validation

`queryGenes` and `pathwayDatabase` null → `ArgumentNullException`. Empty pathway database → empty list.
Gene matching is ordinal, case-sensitive. The query is unioned into the background so the query is part
of the sampled universe. `HypergeometricUpperTail` returns 1.0 for degenerate inputs (x ≤ 0 or any of
N, M, n ≤ 0).

## 4. Algorithm

### 4.1 High-Level Steps

1. Build the distinct query set and the background set (defaulting to the union of pathway members when
   none supplied); union the query into the background.
2. For each pathway: intersect its members with the background (M), count overlap with the query (x).
3. Compute P(X ≥ x) via `HypergeometricUpperTail` using log-Gamma binomial coefficients.
4. Sort pathways ascending by p-value.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Binomial coefficients are evaluated in log-space via the Lanczos approximation of ln Γ to avoid
overflow for large N. `LogChoose(n,k) = ln Γ(n+1) − ln Γ(k+1) − ln Γ(n−k+1)`, returning −∞ (term = 0)
when k ∉ [0, n], which correctly zeroes infeasible partial tables (sampling-without-replacement
constraint).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindPathwayEnrichment | O(P · (M̄ + x̄)) + O(P log P) sort | O(N) | P pathways; M̄ avg members; x̄ avg query size |
| HypergeometricUpperTail | O(x) | O(1) | x summation terms, each O(1) log-Gamma |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.FindPathwayEnrichment(queryGenes, pathwayDatabase, backgroundGenes?)`: per-pathway ORA, sorted by p-value.
- `MetagenomicsAnalyzer.HypergeometricUpperTail(x, bigN, bigM, n)`: core right-tail probability P(X ≥ x).

### 5.2 Current Behavior

The query set is always unioned into the background, and pathway members are intersected with the
background before counting, so M and overlap are measured relative to the actual sampled universe.
Computation is in log-space; the final p is clamped to [0, 1]. This is not a substring/pattern-search
task, so the repository suffix tree is **not applicable** (operations are set membership and arithmetic).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The right-tail hypergeometric p-value `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)` [1][2].
- Sampling without replacement (hypergeometric, not binomial) [1].
- Per-pathway p-values ranked ascending [1].

**Intentionally simplified:**

- Default background = union of pathway members when none is supplied; **consequence:** the user must
  pass an explicit background for a genome-wide N (the cited sources require a caller-defined background).

**Not implemented:**

- Multiple-testing correction (FDR/Benjamini-Hochberg); **users should rely on:** applying their own
  correction over the returned p-values.
- Built-in KEGG/MetaCyc pathway definitions; **users should rely on:** supplying pathway membership.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default background union | Assumption | Changes N when no background passed | accepted | ASM-02; tested in S1 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Overlap x = 0 | p-value = 1.0 | Empty upper sum (INV-02) [2] |
| Degenerate N/M/n ≤ 0 | p-value = 1.0 | No success drawable (INV-03) |
| Empty pathway database | empty result list | Nothing to test |
| Null query or database | ArgumentNullException | Contract |
| Duplicate query genes | counted once | Set semantics |

### 6.2 Limitations

No multiple-testing correction; no built-in pathway databases; assumes gene exchangeability (no
length/abundance weighting). Not a search algorithm — suffix tree not applicable.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through (PNNL §8.2 [2]):** N = 8000 background genes, a gene set of M = 400, a query of
n = 100, overlap x = 20. `HypergeometricUpperTail(20, 8000, 400, 100)` = 7.884747×10⁻⁸ ≈ 7.88×10⁻⁸,
matching the published value.

**API usage example:**

```csharp
var query = new[] { "g1", "g2", "g3" };
var pathways = new Dictionary<string, IReadOnlyCollection<string>>
{
    ["P1"] = new[] { "g1", "g2", "g3", "g4", "g5" }
};
var background = new[] { "g1","g2","g3","g4","g5","g6","g7","g8","g9","g10" };
var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background);
// results[0].PValue == HypergeometricUpperTail(3, 10, 5, 3) == 1.0/12
```

### 7.2 Applications and Use Cases

- **Metabolic capability profiling:** identify which pathways are over-represented among the genes a
  metagenome encodes, relative to a reference gene universe.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Metagenomics/MetagenomicsAnalyzer_FindPathwayEnrichment_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [META-PATHWAY-001-Evidence.md](../../../docs/Evidence/META-PATHWAY-001-Evidence.md)
- Related algorithms: [Functional_Prediction](./Functional_Prediction.md)

## 8. References

1. Boyle EI, Weng S, Gollub J, Jin H, Botstein D, Cherry JM, Sherlock G. 2004. GO::TermFinder — open
   source software for accessing Gene Ontology information and finding significantly enriched Gene
   Ontology terms associated with a list of genes. *Bioinformatics* 20(18):3710–3715.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3037731/ (DOI: 10.1093/bioinformatics/bth456)
2. PNNL Computational Mass Spectrometry. Proteomics Data Analysis in R/Bioconductor, §8.2
   Over-Representation Analysis. https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
