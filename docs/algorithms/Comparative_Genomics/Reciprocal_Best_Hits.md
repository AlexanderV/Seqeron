# Reciprocal Best Hits (RBH)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-RBH-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Reciprocal Best Hits (RBH), also called bidirectional or symmetrical best hits, is the standard operational heuristic for identifying orthologous genes between two genomes. Two genes — one in each genome — are deemed orthologs when each is the other's best hit in the opposite genome [1]. It is the pairwise special case of the mutually-consistent best-hit triangles used to build COGs [2]. RBH is a heuristic (not a phylogenetic reconstruction): it trades some sensitivity for speed and high precision, and is the seed step of many ortholog pipelines.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Orthologs are homologous genes separated by a speciation event; they tend to retain the same function, so functional annotation can be transferred along orthologous links [2]. Distinguishing orthologs from paralogs (genes separated by duplication) is the core problem of comparative genomics. RBH approximates orthology by similarity: the closest reciprocal match across genomes is taken as the orthologous counterpart [1].

### 2.2 Core Model

Let genome G1 have genes {a} and G2 have genes {b}, each with a similarity score s(x, y) to candidates in the other genome. Define the best hit of x into a genome as the candidate maximizing the score among those passing a significance gate. Per Moreno-Hagelsieb & Latimer (2008): "two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome" [1]. Formally, (a, b) is an RBH pair iff:

- `bestHit(a → G2) = b` and `bestHit(b → G1) = a`.

The best hit is chosen by sorting candidates "from highest to lowest bit-score, then, if the bit-scores were identical, from smallest to highest E-values" [1] — i.e. maximum score with a deterministic tie-break. A candidate qualifies only if it passes the significance gate: "maximum E-value threshold of 1×10⁻⁶" and "coverage of at least 50% of any of the protein sequences in the alignments" [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every returned pair is reciprocal (a's best hit is b and b's best hit is a). | RBH definition [1] |
| INV-02 | The result is a matching: no gene of either genome appears in two pairs. | each gene has at most one best hit; a pair requires both directions to agree [1] |
| INV-03 | Each pair carries the actual hit identity, coverage and alignment length of the winning hit. | the hit metrics are those of the best hit, not placeholders [1] |
| INV-04 | Deterministic and order-independent for fixed input. | the best hit is unique via the score + deterministic tie-break [1] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Reciprocal Best Hits | One-directional best hit |
|--------|----------------------|--------------------------|
| Orthology criterion | mutual (both directions agree) | single direction |
| False positives | low (reciprocity filters them) | high (asymmetric hits) |
| Output | partial matching (≤ min(\|G1\|,\|G2\|) pairs) | up to \|G1\| hits |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genome1Genes` | `IReadOnlyList<Gene>` | required | genes of genome 1 | non-null; genes need a non-empty `Sequence` to participate |
| `genome2Genes` | `IReadOnlyList<Gene>` | required | genes of genome 2 | non-null |
| `minIdentity` | `double` | 0.3 | minimum similarity score for a qualifying hit (significance gate, mapped from the 1×10⁻⁶ E-value gate) [1] | 0–1 |
| `minCoverage` | `double` | 0.5 | minimum coverage fraction for a qualifying hit (≥ 50% coverage gate) [1] | 0–1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Gene1Id` | `string` | id of the genome-1 gene |
| `Gene2Id` | `string` | id of its reciprocal best hit in genome 2 |
| `Identity` | `double` | similarity score of the hit (5-mer Jaccard; see ASM-01) |
| `Coverage` | `double` | shared-k-mer coverage of the shorter sequence |
| `AlignmentLength` | `int` | min length of the two sequences |

### 3.3 Preconditions and Validation

Null genome lists throw `ArgumentNullException`. Genes whose `Sequence` is null or empty are skipped (similarity is undefined without a sequence). An empty genome (no sequence-bearing genes) yields no pairs. Sequences shorter than k = 5 produce a similarity of 0 and so never qualify. Comparison is case-insensitive (k-mers upper-cased).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; drop genes without a sequence.
2. For each genome-1 gene, compute the best qualifying hit into genome 2 (max score, deterministic tie-break, gated by `minIdentity`/`minCoverage`).
3. For each genome-2 gene, compute the id of its best qualifying hit back into genome 1.
4. Emit (a, b) when a's best hit is b and b's best hit is a (reciprocity).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Best-hit selection: maximum similarity score; ties broken by larger coverage, then by ordinal gene id, making the winner unique (mapping the bit-score → E-value tie-break of [1] onto the alignment-free score). Significance gate: `identity ≥ minIdentity (0.3)` and `coverage ≥ minCoverage (0.5)` [1]. Similarity is a 5-mer Jaccard index (alignment-free; ASM-01).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindReciprocalBestHits` | O(n × m × L) | O((n + m) × L) | n, m = gene counts; L = k-mer set construction per pair; all-vs-all best-hit scan in both directions |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.FindReciprocalBestHits(genome1Genes, genome2Genes, minIdentity, minCoverage)`: the canonical RBH entry point.
- `ComparativeGenomics.FindOrthologs(...)`: delegates to `FindReciprocalBestHits` (orthology by RBH is the same criterion).
- `ComparativeGenomics.FindBestHit(...)` (private): selects the single best qualifying hit with deterministic tie-break.

### 5.2 Current Behavior

The similarity used to rank candidates is an alignment-free 5-mer Jaccard index (`CalculateSequenceSimilarity`), because the Analysis project does not reference an alignment engine (ASM-01). `Identity` is the Jaccard score, `Coverage` is the fraction of the shorter sequence's k-mers that are shared, and `AlignmentLength` is the minimum of the two sequence lengths. The result is returned as a fully materialized list, so enumeration is side-effect-free and repeatable.

**Search-infrastructure decision (suffix tree):** not used. RBH ranks candidates by a scoring function over k-mer sets (set intersection/union), not by exact-substring occurrence; the repository suffix tree answers exact-match occurrence queries and does not fit a scoring-based all-vs-all best-hit search.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Reciprocity criterion: a pair is returned iff each gene is the other's best hit [1].
- Best hit = maximum score with a deterministic tie-break [1].
- Significance/coverage gate: `minIdentity` and `minCoverage` (≥ 50% coverage) filter qualifying hits [1].

**Intentionally simplified:**

- Similarity metric: 5-mer Jaccard instead of BLAST bit-score; **consequence:** the ranking is alignment-free, so for near-identical candidates the tie-winner can differ from a bit-score ranking (ASM-01); identical sequences still score 1.0.
- E-value gate mapped to a minimum-similarity threshold; **consequence:** no statistical E-value is computed; the user supplies `minIdentity` instead.

**Not implemented:**

- Soft-masking and Smith-Waterman alignment recommended by [1] for improved RBH recovery; **users should rely on:** the documented `minIdentity`/`minCoverage` gates, or an external alignment-based pipeline for production ortholog calling.
- In-paralog clustering (InParanoid/OrthoMCL extensions); **users should rely on:** `ComparativeGenomics.FindParalogs` (COMPGEN-ORTHO-001) for within-genome best hits.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | 5-mer Jaccard ranking instead of bit-score | Assumption | tie ordering for near-identical candidates | accepted | ASM-01; reciprocity/gate/tie-break are source-backed |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null genome list | `ArgumentNullException` | repository contract |
| Empty genome | empty result | a pair needs one gene from each genome [1] |
| Gene without sequence | skipped | similarity undefined |
| One-directional best hit | excluded | reciprocity required [1][2] |
| Sub-threshold pair | excluded | significance/coverage gate [1] |
| Sequence shorter than k = 5 | similarity 0 → never qualifies | k-mer set empty |

### 6.2 Limitations

RBH is a heuristic and can miss orthologs lost by gene-family expansion/contraction, and can pair the wrong member of a recently expanded family. It is alignment-free here, so it does not detect orthology across distant homologs with low k-mer overlap. It identifies one-to-one orthologs only; many-to-many orthologous relationships and in-paralog co-clustering are out of scope.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var g1 = new[] { new ComparativeGenomics.Gene("a1", "G1", 0, 14, '+', "ACGTACGTACGTAC") };
var g2 = new[] { new ComparativeGenomics.Gene("b1", "G2", 0, 14, '+', "ACGTACGTACGTAC") };
var rbh = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();
// rbh[0] => (Gene1Id "a1", Gene2Id "b1", Identity 1.0, Coverage 1.0, AlignmentLength 14)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_FindReciprocalBestHits_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindReciprocalBestHits_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [COMPGEN-RBH-001-Evidence.md](../../../docs/Evidence/COMPGEN-RBH-001-Evidence.md)
- Related algorithms: [Ortholog_Identification](./Ortholog_Identification.md)

## 8. References

1. Moreno-Hagelsieb G, Latimer K. 2008. Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585
2. Tatusov RL, Koonin EV, Lipman DJ. 1997. A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631 (method described in NCBI Handbook, "The Clusters of Orthologous Groups (COGs) Database", https://www.ncbi.nlm.nih.gov/books/NBK21090/)
3. Ondov BD, et al. 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 17:132. https://doi.org/10.1186/s13059-016-0997-x
