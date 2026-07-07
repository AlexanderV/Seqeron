# Ortholog Identification (Reciprocal Best Hits) and Paralog Identification

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-ORTHO-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Identifies orthologous gene pairs between two genomes and paralogous gene pairs within one genome. Orthologs are detected by the **Reciprocal Best Hit (RBH/BBH)** criterion: two genes are orthologs iff each is the other's best hit in the opposite genome [3]. Paralogs (recent in-paralogs) are detected as within-genome mutual best hits [1][4]. The method is a heuristic for the evolutionary relationships defined by Fitch (1970) [1]; it is alignment-free in this implementation (best hits are ranked by a k-mer similarity score, see §5.2 and §5.4) but the reciprocity logic is exact.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Fitch (1970) defined two subclasses of homology: genes whose divergence is due to **speciation** are **orthologous**; genes whose divergence is due to **gene duplication** within an organism's history are **paralogous** [1]. Orthologs tend to retain function across species, which is why ortholog detection underpins cross-species functional annotation transfer [2]. Operationally, orthologs between two genomes are approximated by reciprocal best hits [2][3], and recent paralogs (in-paralogs) by within-genome best hits [4].

### 2.2 Core Model

Let G1, G2 be gene sets with sequences and let `sim(x, y)` be a symmetric similarity score.

- **Best hit:** `BH(x, T) = argmax_{y ∈ T} sim(x, y)` among targets `y` whose hit qualifies (similarity ≥ minIdentity and coverage ≥ minCoverage); ties broken deterministically [3].
- **Ortholog (RBH):** a pair `(g1 ∈ G1, g2 ∈ G2)` is an ortholog pair iff `BH(g1, G2) = g2` **and** `BH(g2, G1) = g1` [3]. A one-directional best hit (`BH(g1,G2)=g2` but `BH(g2,G1)≠g1`) is **not** an ortholog [2].
- **Paralog (in-paralog):** within a single genome G, a pair `(p, q ∈ G, p ≠ q)` is a paralog pair iff `BH(p, G\{p}) = q` **and** `BH(q, G\{q}) = p` [1][4].

Coverage gate: a hit qualifies only if it covers ≥ 50% of (the shorter) sequence [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every ortholog pair is reciprocal | RBH definition is symmetric by construction [3] |
| INV-02 | RBH output is a matching: no gene appears in two pairs | best hit is unique per gene (deterministic argmax) [3] |
| INV-03 | Every paralog pair is an unordered pair of distinct within-genome genes that are mutual best hits | in-paralog definition [1][4] |
| INV-04 | Output is deterministic and order-independent | deterministic best-hit ranking with ordinal tie-break [3] |
| INV-05 | Pairs below minIdentity or minCoverage are excluded | qualification gate [3] |

### 2.5 Comparison with Related Methods

| Aspect | RBH (this method) | COG / OrthoMCL clustering |
|--------|-------------------|---------------------------|
| Genomes per run | 2 | ≥ 3 (graph clustering) [2][5] |
| Output | pairwise matching | ortholog groups (with in-paralogs) [2][5] |
| Paralogs | reported separately by `FindParalogs` | merged into groups [2][5] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| genome1Genes / genome2Genes / genes | `IReadOnlyList<Gene>` | required | genes, each with a non-empty `Sequence` | non-null; genes without a sequence are skipped |
| minIdentity | `double` | 0.3 | minimum similarity score to qualify | [0, 1] |
| minCoverage | `double` | 0.5 | minimum coverage fraction to qualify | [0, 1]; 0.5 per [3] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IEnumerable<OrthologPair>` | reciprocal best-hit pairs; for `FindParalogs`, each unordered within-genome pair reported once |
| OrthologPair.Identity | `double` | k-mer Jaccard similarity of the pair (ranking score) |
| OrthologPair.Coverage | `double` | shared-k-mer fraction of the shorter sequence |

### 3.3 Preconditions and Validation

Null gene lists throw `ArgumentNullException`. Genes whose `Sequence` is null/empty are skipped. Sequences are upper-cased for k-mer extraction (case-insensitive). Empty genome (no sequenced genes) → no orthologs; genome with < 2 sequenced genes → no paralogs. Sequences shorter than the k-mer length (5) contribute no k-mers and never qualify.

## 4. Algorithm

### 4.1 High-Level Steps

1. Drop genes without a sequence.
2. For each gene, compute its best qualifying hit in the opposite set (orthologs) or in the rest of the same set (paralogs).
3. Keep a pair only when the best-hit relation is reciprocal.
4. For paralogs, report each unordered pair once.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Similarity score `sim` = Jaccard of 5-mer sets (alignment-free, cf. Mash [6]).
- Coverage = shared 5-mers / 5-mers of the shorter sequence; gate ≥ 0.5 [3].
- Best-hit tie-break: higher score, then higher coverage, then smaller ordinal gene id — guarantees a unique, deterministic best hit [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| FindOrthologs | O(n·m·L) | O((n+m)·L) | n,m = gene counts; L = sequence length (k-mer build + intersect) |
| FindParalogs | O(n²·L) | O(n·L) | all within-genome pairs |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.FindOrthologs(genome1Genes, genome2Genes, minIdentity, minCoverage)`: RBH ortholog pairs.
- `ComparativeGenomics.FindParalogs(genes, minIdentity, minCoverage)`: within-genome mutual-best-hit paralog pairs.
- `ComparativeGenomics.FindBestHit(...)` (private): deterministic best-qualifying-hit selection.

### 5.2 Current Behavior

`FindOrthologs` was corrected from a prior one-directional best-hit implementation (which violated INV-01) to the reciprocal criterion. `FindParalogs` is new. Best hits are ranked by an alignment-free 5-mer Jaccard similarity because the `Seqeron.Genomics.Analysis` project does not reference an alignment/bit-score implementation; the reciprocity, coverage gate, and thresholds follow the cited sources.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Reciprocal best-hit ortholog criterion [3]; symmetrical best hits [2].
- Within-genome mutual-best-hit (in-paralog) detection [1][4].
- ≥ 50% coverage qualification gate [3].
- Deterministic best-hit selection with tie-break [3].

**Intentionally simplified:**

- Best-hit ranking score: alignment-free 5-mer Jaccard similarity instead of a BLAST bit-score [3]; **consequence:** the *score values* differ from bit-scores, and for sequences with identical 5-mer multisets the tie-break (not alignment) decides the winner. Reciprocity, ordering, and thresholds are preserved (see §5.4 and Evidence Assumption 1).

**Not implemented:**

- Multi-genome ortholog groups with merged in-paralogs (COG/OrthoMCL) [2][5]; **users should rely on:** running pairwise RBH per genome pair, or an external COG/OrthoMCL tool.
- Outparalog vs in-paralog discrimination relative to a specific speciation node [4]; **users should rely on:** `FindParalogs` reports closest within-genome relatives only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | k-mer Jaccard as best-hit score | Assumption | score values are not bit-scores; tie-break decides identical-k-mer cases | accepted | Evidence Assumption 1; order-preserving for test datasets |
| 2 | RBH does not reuse the suffix tree | Deviation | n/a | accepted | best hit is a scoring search (argmax similarity), not exact-match enumeration; the suffix tree fits exact occurrence search, not similarity ranking |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null gene list | `ArgumentNullException` | contract §3.3 |
| empty genome | no orthologs | a pair needs one gene from each genome |
| single-gene genome | no paralogs | a pair needs two within-genome genes |
| gene without sequence | skipped | similarity undefined |
| one-directional best hit | not returned | RBH requires reciprocity [2][3] |
| sub-threshold pair | excluded | qualification gate [3] |

### 6.2 Limitations

Pairwise only (two genomes). Alignment-free score is a heuristic proxy for sequence identity; does not detect outparalogs vs in-paralogs, fusion/fission genes, or many-to-many orthology. Very short sequences (< 5 nt) never qualify.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var g1 = new[] { new ComparativeGenomics.Gene("a1", "G1", 0, 14, '+', "ACGTACGTACGTAC") };
var g2 = new[] { new ComparativeGenomics.Gene("b1", "G2", 0, 14, '+', "ACGTACGTACGTAC") };
var orthologs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();
// orthologs => one pair a1 <-> b1 (mutual best hits, Jaccard 1.0)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_FindOrthologs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ComparativeGenomics_FindOrthologs_Tests.cs) — covers INV-01..INV-05
- Evidence: [COMPGEN-ORTHO-001-Evidence.md](../../../docs/Evidence/COMPGEN-ORTHO-001-Evidence.md)
- Related algorithms: [Synteny_Block_Detection](../Comparative_Genomics/Synteny_Block_Detection.md)

## 8. References

1. Fitch WM. 1970. Distinguishing homologous from analogous proteins. *Systematic Zoology* 19(2):99–113. https://pmc.ncbi.nlm.nih.gov/articles/PMC3178060/
2. Tatusov RL, Koonin EV, Lipman DJ. 1997. A genomic perspective on protein families. *Science* 278(5338):631–637. https://doi.org/10.1126/science.278.5338.631
3. Moreno-Hagelsieb G, Latimer K. 2008. Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324. https://doi.org/10.1093/bioinformatics/btm585
4. Remm M, Storm CEV, Sonnhammer ELL. 2001. Automatic clustering of orthologs and in-paralogs from pairwise species comparisons. *J. Mol. Biol.* 314(5):1041–1052. https://doi.org/10.1006/jmbi.2000.5197
5. Li L, Stoeckert CJ, Roos DS. 2003. OrthoMCL: identification of ortholog groups for eukaryotic genomes. *Genome Res.* 13(9):2178–2189. https://doi.org/10.1101/gr.1224503
6. Ondov BD, et al. 2016. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 17:132. https://doi.org/10.1186/s13059-016-0997-x
