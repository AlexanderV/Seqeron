# Taxonomic Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-CLASS-001 |
| Related Projects | N/A |
| Implementation Status | Implemented (Kraken k-mer / LCA / RTL) |
| Last Reviewed | 2026-06-17 |

## 1. Overview

K-mer based taxonomic classification assigns taxonomic labels to query sequences by comparing their
constituent k-mers against a reference database. The repository implements the **Kraken** algorithm
(Wood and Salzberg, 2014): a taxonomy tree with a lowest-common-ancestor (LCA) operation, a
canonical-k-mer → taxon database in which a k-mer shared by several reference taxa is stored as the
**LCA of those taxa**, and a per-read **root-to-leaf (RTL) maximum-weight-path** classifier over the
taxonomy. Reads whose k-mers hit no database entry are left unclassified (reported as the taxonomy
root).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Exact k-mer classification is the standard fast metagenomic strategy for read assignment against
reference genomes: each read is decomposed into fixed-length substrings compared to a precomputed
lookup structure, trading alignment sensitivity for speed and indexing efficiency (Wood and Salzberg,
2014).

### 2.2 Core Model — Kraken

**Database build.** A database record is "a k-mer and the LCA of all organisms whose genomes contain
that k-mer." As reference sequences are processed, "if a k-mer from a sequence has had its LCA value
previously set, then the LCA of the stored value and the current sequence's taxon is calculated" and
stored. K-mers are stored in **canonical** form:

$$
canonical(kmer) = \min(kmer, reverse\_complement(kmer))
$$

reflecting the double-stranded nature of DNA.

**Per-read classification.** The set of taxa hit by a read's k-mers, **together with their ancestors
in the taxonomy tree, form the *classification tree*** — a pruned subtree. "Each node in the
classification tree is weighted with the number of k-mers … that mapped to the taxon associated with
that node." "Each root-to-leaf (RTL) path … is scored by calculating the sum of all node weights
along the path. The maximum scoring RTL path … is the classification path," and its leaf is the
read's assigned label. Tie-break: "if there are multiple maximally scoring paths, the LCA of all
those paths' leaves is selected." Reads with no k-mer hits "are left unclassified."

**Confidence (Kraken 2 C/Q).** The score reported per read is `C/Q`, where `C` is the number of
k-mers mapped to a taxon **in the clade rooted at the assigned label**, and `Q` is the number of
non-ambiguous k-mers queried.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output count equals input read count | `ClassifyReads` yields exactly one result per input tuple |
| INV-02 | `0 <= Confidence <= 1` | Confidence is `C/Q` with `0 <= C <= Q`, else `0` |
| INV-03 | `MatchedKmers (C) <= TotalKmers (Q)` | The assigned clade's k-mers are a subset of the queried k-mers |
| INV-04 | Reads with no k-mer hits → root / unclassified | No hits ⇒ empty classification tree ⇒ `TaxonId = TaxonomyTree.RootId` |
| INV-05 | The assigned taxon lies in the taxonomy tree | The RTL leaf / its LCA are nodes of the supplied tree |
| INV-06 | A k-mer shared by several taxa maps to their LCA | DB build folds the stored taxon with each new owner via `Lca` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[BuildKmerDatabase] referenceSequences` | `IEnumerable<(int TaxonId, string Sequence)>` | required | Labeled references | TaxonId must exist in the tree; sequences shorter than `k` are skipped |
| `[BuildKmerDatabase] taxonomy` | `TaxonomyTree` | required | Tree used for the LCA of shared k-mers | — |
| `[BuildKmerDatabase] k` | `int` | `31` | K-mer length | Must be positive |
| `[ClassifyReads] reads` | `IEnumerable<(string Id, string Sequence)>` | required | Reads to classify | Empty / shorter-than-`k` reads → unclassified |
| `[ClassifyReads] kmerDatabase` | `IReadOnlyDictionary<string, int>` | required | Canonical-k-mer → taxon-id map | Values should be taxa present in the tree |
| `[ClassifyReads] taxonomy` | `TaxonomyTree` | required | Tree for parent chains + LCA | — |
| `[ClassifyReads] k` | `int` | `31` | K-mer length | Must be positive |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[BuildKmerDatabase] return value` | `Dictionary<string, int>` | Canonical-k-mer → taxon-id (LCA of owning taxa) |
| `ReadId` | `string` | Input read identifier |
| `TaxonId` | `int` | Assigned taxon = leaf of the max-scoring RTL path (root id when unclassified) |
| `TaxonName` | `string` | Assigned taxon's name |
| `Rank` | `string` | Assigned taxon's rank label |
| `RtlScore` | `int` | Sum of node weights on the winning root-to-leaf path |
| `Confidence` | `double` | `C/Q` |
| `MatchedKmers` | `int` | `C` — k-mers in the clade rooted at the assigned label |
| `TotalKmers` | `int` | `Q` — non-ambiguous k-mers queried |
| `Kingdom`…`Species` | `string` | The standard ranks read off the assigned taxon's lineage (empty where absent) |

### 3.3 Preconditions and Validation

Both methods uppercase input internally and skip any k-mer containing a non-`A/C/G/T` character
(such k-mers are not counted in `Q`). Empty or shorter-than-`k` reads return an unclassified result
with zero counts. Null arguments throw `ArgumentNullException`; a non-positive `k` throws
`ArgumentOutOfRangeException`; a reference taxon absent from the tree throws `KeyNotFoundException`.
The `TaxonomyTree` constructor validates that there is exactly one self-parented root, no duplicate
ids, and that every non-root parent is present.

## 4. Algorithm

### 4.1 High-Level Steps

**Database build:** for each reference, enumerate canonical length-`k` k-mers (skipping ambiguous
ones); insert with the reference's taxon, or, if the k-mer is already present, replace the stored
taxon with its LCA against the new reference's taxon.

**Per-read classification:** enumerate the read's non-ambiguous canonical k-mers (counting `Q`),
look each up in the database, and tally per-taxon hit counts. If there are no hits, assign the root.
Otherwise, determine the classification-tree leaves (hit taxa that are not a proper ancestor of
another hit taxon), score each leaf by the sum of hit weights on its root path, take the
maximum-scoring leaf (or the LCA of the tied maximal leaves), and report it with `C` = the assigned
clade's k-mer count and `C/Q`.

### 4.2 Data Structures

`TaxonomyTree` stores `TaxonNode(Id, Name, Rank, ParentId)` records keyed by id; `Lca(a,b)` collects
`a`'s ancestor set and walks `b` upward to the first shared node (root otherwise); `Lca(IEnumerable)`
folds the pairwise LCA. The database is a `Dictionary<string,int>`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `BuildKmerDatabase` | `O(Σ m · h)` | `O(u)` | `m` = ref length, `h` = tree height (LCA), `u` = unique k-mers |
| `ClassifyReads` | `O(r · (l + leaves · h))` | `O(distinct hit taxa)` | `r` = reads, `l` = read length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs),
[TaxonomyTree.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/TaxonomyTree.cs)

- `TaxonomyTree` / `TaxonNode`: taxonomy data model with parent chains and `Lca`.
- `MetagenomicsAnalyzer.BuildKmerDatabase(IEnumerable<(int,string)>, TaxonomyTree, int)`: LCA k-mer database.
- `MetagenomicsAnalyzer.ClassifyReads(IEnumerable<(string,string)>, IReadOnlyDictionary<string,int>, TaxonomyTree, int)`: RTL classification.

### 5.2 Conformance to Theory / Spec

**Implemented (faithful to Kraken):** canonical k-mer indexing; database k-mer → **LCA of owning
taxa**; per-read classification tree weighted by k-mer count; **maximum-scoring root-to-leaf path**
with **LCA-of-leaves** tie-break; no-hit → unclassified (root); `C/Q` confidence; ambiguous-k-mer
exclusion from `Q`.

**Not implemented:** the minimizer / spaced-seed and compact-hash index of Kraken 2 (this is the
exact-k-mer Kraken 1 model); paired-end concatenation; the `--confidence` re-walk threshold.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| Empty / shorter-than-`k` read | Unclassified (root), `Q = 0` |
| No database matches | Unclassified (root), `Q` still recorded |
| K-mer with non-`ACGT` | Skipped (not counted in `Q`) |
| K-mer shared by several taxa (DB build) | Stored as the LCA of those taxa |
| Read split equally between sibling species | Assigned their genus (LCA of tied leaves) |
| Read split equally across genera | Assigned their family (LCA of tied leaves) |

### 6.2 Limitations

This is the exact-k-mer Kraken 1 model; it does not use minimizers/spaced seeds and depends on the
caller-supplied reference database and taxonomy. Choice of `k` trades specificity against
sensitivity.

## 7. Examples and Related Material

- [Taxonomic_Profile.md](./Taxonomic_Profile.md) — downstream aggregation of per-read classifications.
- [Alpha_Diversity.md](./Alpha_Diversity.md) — diversity metrics after profiling.

## 8. References

1. Wood, D. E., and S. L. Salzberg. 2014. Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 15:R46. doi:10.1186/gb-2014-15-3-r46.
2. Kraken Documentation. https://ccb.jhu.edu/software/kraken/
3. Wikipedia contributors. Metagenomics. Wikipedia. https://en.wikipedia.org/wiki/Metagenomics
