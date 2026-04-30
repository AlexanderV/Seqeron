# Taxonomic Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-CLASS-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

K-mer based taxonomic classification assigns taxonomic labels to query sequences by comparing their constituent k-mers against a reference database. This style of classifier is widely used in metagenomics because exact k-mer matching is fast and can scale to many short reads. The repository implements a canonical-k-mer workflow with separate database-building and read-classification entry points. The current implementation is a simplified winner-take-all classifier: it selects the taxon with the largest number of matching k-mers and reports confidence as the fraction of non-ambiguous query k-mers assigned to that winning taxon.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Exact k-mer classification is a standard metagenomic strategy for rapid read assignment against reference genomes. By decomposing each read into fixed-length substrings and comparing them to a precomputed lookup structure, the method trades alignment sensitivity for speed and database indexing efficiency (Wood and Salzberg, 2014).

### 2.2 Core Model

The core classification model has two stages:

1. Build a reference map from canonical k-mers to taxonomic labels.
2. For each query read, count how many canonical query k-mers map to each taxon and choose the taxon with the highest count.

The canonical-k-mer transform is:

$$
canonical(kmer) = \min(kmer, reverse\_complement(kmer))
$$

which identifies a k-mer with its reverse complement, reflecting the double-stranded nature of DNA. The repository reports confidence as:

$$
Confidence = \frac{MatchedKmers}{TotalKmers}
$$

where `TotalKmers` counts non-ambiguous query k-mers and `MatchedKmers` counts only those supporting the winning taxon.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output count equals input read count | `ClassifyReads` yields exactly one `TaxonomicClassification` result per input tuple |
| INV-02 | `0 <= Confidence <= 1` | Confidence is the ratio `MatchedKmers / TotalKmers` when `TotalKmers > 0`, else `0` |
| INV-03 | `MatchedKmers <= TotalKmers` | The winning taxon count is one category within the total counted query k-mers |
| INV-04 | Reads with no qualifying or matching k-mers are reported as `Unclassified` | The implementation returns `Kingdom = "Unclassified"` when the taxon count map is empty |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[BuildKmerDatabase] referenceGenomes` | `IEnumerable<(string TaxonId, string Sequence)>` | required | Reference genomes used to create the lookup database | Sequences shorter than `k` are skipped |
| `[BuildKmerDatabase] k` | `int` | `31` | K-mer length for database construction | Longer `k` increases specificity; shorter `k` increases sensitivity |
| `[ClassifyReads] reads` | `IEnumerable<(string Id, string Sequence)>` | required | Reads to classify | Sequences shorter than `k` return `Unclassified` immediately |
| `[ClassifyReads] kmerDatabase` | `IReadOnlyDictionary<string, string>` | required | Canonical k-mer to taxonomy-string mapping | Query k-mers are looked up exactly after canonicalization |
| `[ClassifyReads] k` | `int` | `31` | K-mer length for classification | Must match the database's intended k-mer length to be meaningful |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[BuildKmerDatabase] return value` | `Dictionary<string, string>` | Canonical k-mer database keyed by canonical DNA k-mer string |
| `ReadId` | `string` | Input read identifier |
| `Kingdom` | `string` | Kingdom assignment or `Unclassified` when no classification is made |
| `Phylum` | `string` | Parsed phylum assignment |
| `Class` | `string` | Parsed class assignment |
| `Order` | `string` | Parsed order assignment |
| `Family` | `string` | Parsed family assignment |
| `Genus` | `string` | Parsed genus assignment |
| `Species` | `string` | Parsed species assignment |
| `Confidence` | `double` | Winning-taxon support divided by total non-ambiguous query k-mers |
| `MatchedKmers` | `int` | Number of query k-mers assigned to the winning taxon |
| `TotalKmers` | `int` | Number of query k-mers considered after skipping ambiguous k-mers |

### 3.3 Preconditions and Validation

Both database construction and classification uppercase input sequences internally. K-mers containing any character outside `A/C/G/T` are skipped. Query reads that are empty or shorter than `k` return an `Unclassified` result with zero counts. When no database entry matches any qualifying query k-mer, the method still records `TotalKmers` but returns `Unclassified`. Taxonomy strings are parsed into rank fields by splitting on `|` or `;` in the fixed order `kingdom, phylum, class, order, family, genus, species`.

## 4. Algorithm

### 4.1 High-Level Steps

1. For database construction, enumerate every length-`k` substring of each reference genome.
2. Skip any reference k-mer containing characters outside `A/C/G/T`.
3. Canonicalize each reference k-mer by comparing it with its reverse complement.
4. Insert the canonical k-mer into the lookup database if it is not already present.
5. For read classification, enumerate canonical query k-mers after uppercasing and skipping ambiguous substrings.
6. Count taxon hits per canonical query k-mer.
7. Select the taxon with the highest hit count, compute confidence as `MatchedKmers / TotalKmers`, and parse the taxonomy string into rank fields.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The reference database is a dictionary from canonical k-mer string to taxonomy string. Classification is winner-take-all: the method orders taxon hit counts descending and selects the first entry. Canonicalization uses the lexicographically smaller of the forward k-mer and its reverse complement. The default `k` value is `31`, matching the older document's stated default and the source signature.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `BuildKmerDatabase` | `O(n * m)` | `O(k * u)` | `n` = reference genomes, `m` = average genome length, `u` = unique canonical k-mers |
| `ClassifyReads` | `O(r * l)` | `O(1)` per read plus taxon hit counts | `r` = reads, `l` = average read length |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.BuildKmerDatabase(IEnumerable<(string TaxonId, string Sequence)>, int)`: Builds the canonical k-mer lookup table from reference genomes.
- `MetagenomicsAnalyzer.ClassifyReads(IEnumerable<(string Id, string Sequence)>, IReadOnlyDictionary<string, string>, int)`: Classifies reads by canonical k-mer hit counts.

### 5.2 Current Behavior

The current database builder keeps only the first taxon seen for a canonical k-mer; later references with the same canonical k-mer do not overwrite the existing mapping. Classification skips ambiguous query k-mers instead of attempting approximate matching. The chosen taxon is the one with the largest hit count rather than a lowest-common-ancestor reconciliation over multiple taxa. Taxonomy strings are parsed only by position into the seven fixed ranks stored by `TaxonomicClassification`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact k-mer matching against a prebuilt reference database.
- Reverse-complement-aware canonical k-mer representation.
- Confidence reporting as the fraction of counted query k-mers supporting the chosen classification.

**Intentionally simplified:**

- Shared canonical k-mers are stored with the first observed taxon instead of an ambiguity-aware or lowest-common-ancestor representation; **consequence:** database build order can affect the label assigned to shared k-mers.
- Read classification selects the taxon with the highest hit count directly; **consequence:** conflicting k-mer evidence is not reconciled with a Kraken-style lowest-common-ancestor rule.

**Not implemented:**

- Lowest-common-ancestor resolution over a taxonomic tree; **users should rely on:** no current alternative in this class.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | First-taxonomy retention for duplicate canonical k-mers | Deviation | Shared reference k-mers do not preserve all candidate taxa | accepted | The dictionary insertion check keeps the earliest mapping only |
| 2 | Winner-take-all classification | Assumption | Mixed evidence across taxa is compressed to one label without explicit uncertainty decomposition | accepted | Confidence records winning support but not a full evidence distribution |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | Returns `Unclassified` with `Confidence = 0` and zero counts | Guard clause for empty input |
| Sequence length shorter than `k` | Returns `Unclassified` | No valid length-`k` substring can be extracted |
| Empty database or no database matches | Returns `Unclassified` with `TotalKmers` recorded | The taxon count map remains empty |
| Query k-mer containing non-`ACGT` characters | K-mer is skipped | The implementation filters ambiguous k-mers before lookup |

### 6.2 Limitations

This implementation is a simplified exact-match classifier. It does not model phylogenetic relationships, reconcile shared k-mers through a taxonomic tree, or provide approximate matching for divergent reads. Choice of `k` still trades specificity against sensitivity: smaller `k` increases false positives, while larger `k` can miss divergent or short reads.

## 7. Examples and Related Material

- [Taxonomic_Profile.md](./Taxonomic_Profile.md) documents the downstream aggregation of per-read classifications into a sample profile.
- [Alpha_Diversity.md](./Alpha_Diversity.md) documents diversity metrics that can be computed after taxonomic profiling.

## 8. References

1. Wood, D. E., and S. L. Salzberg. 2014. Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 15:R46. doi:10.1186/gb-2014-15-3-r46.
2. Kraken Documentation. https://ccb.jhu.edu/software/kraken/
3. Wikipedia contributors. Metagenomics. Wikipedia. https://en.wikipedia.org/wiki/Metagenomics

