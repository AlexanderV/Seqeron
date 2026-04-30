# Taxonomic Profile Generation

| Field | Value |
|-------|-------|
| Algorithm Group | Metagenomics |
| Test Unit ID | META-PROF-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Taxonomic profiling aggregates read-level classifications into a community-level summary of taxon abundance. In metagenomics, this converts per-read labels into a sample profile that can be interpreted biologically and summarized with diversity indices. The repository implementation generates abundance maps for selected taxonomic ranks and computes Shannon and Simpson diversity from species-level abundances. It is simplified in that the output record stores only kingdom, phylum, genus, and species abundances even though the upstream classification record contains additional intermediate ranks.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Taxonomic profiling summarizes the composition of a microbial community by aggregating classifications over taxonomic ranks. This kind of profile is a common downstream representation in metagenomics because it translates large numbers of read-level assignments into interpretable abundance distributions and diversity summaries (Segata et al., 2012).

### 2.2 Core Model

For a taxon $i$ with count $c_i$, relative abundance is computed as:

$$
abundance_i = \frac{c_i}{\sum_{j=1}^{n} c_j}
$$

where the denominator excludes reads treated as unclassified by the profiling workflow.

The repository then computes species-level Shannon diversity (Shannon, 1948):

$$
H = -\sum_{i=1}^{S} p_i \ln(p_i)
$$

and Simpson concentration (Simpson, 1949):

$$
D = \sum_{i=1}^{S} p_i^2
$$

from the species-abundance map.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TotalReads` equals the number of input classification records | The method materializes the full input and counts it before filtering |
| INV-02 | `ClassifiedReads` excludes records whose `Kingdom` is `Unclassified` or empty | Those records are filtered before abundance maps are built |
| INV-03 | `KingdomAbundance` sums to `1` when `ClassifiedReads > 0`; lower-rank maps do so only when every retained read has that rank populated | Each stored map divides retained rank counts by `classifiedReads`, but empty phylum/genus/species values are skipped before counting |
| INV-04 | Shannon and Simpson diversity are computed only from species-level abundances | The implementation passes `SpeciesAbundance.Values` into the diversity helpers |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `classifications` | `IEnumerable<TaxonomicClassification>` | required | Per-read taxonomic assignments to aggregate into a profile | Reads with `Kingdom == "Unclassified"` or empty kingdom are excluded from abundance denominators |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `KingdomAbundance` | `IReadOnlyDictionary<string, double>` | Relative abundances by kingdom |
| `PhylumAbundance` | `IReadOnlyDictionary<string, double>` | Relative abundances by phylum |
| `GenusAbundance` | `IReadOnlyDictionary<string, double>` | Relative abundances by genus |
| `SpeciesAbundance` | `IReadOnlyDictionary<string, double>` | Relative abundances by species |
| `ShannonDiversity` | `double` | Shannon diversity computed from species abundances |
| `SimpsonDiversity` | `double` | Simpson concentration computed from species abundances |
| `TotalReads` | `int` | Total number of input classification records |
| `ClassifiedReads` | `int` | Number of records retained after the unclassified filter |

### 3.3 Preconditions and Validation

Empty input returns a profile with `TotalReads = 0`, `ClassifiedReads = 0`, empty abundance maps, and zero diversity metrics. Records whose `Kingdom` is `Unclassified` or empty are excluded entirely from abundance denominators. Empty phylum, genus, and species values are excluded from their rank-specific maps even for otherwise classified reads.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize the input classifications and count them as `TotalReads`.
2. Filter out records whose kingdom is `Unclassified` or empty.
3. Count classified reads by kingdom, phylum, genus, and species.
4. Convert each rank count dictionary into relative abundances using `ClassifiedReads` as the denominator.
5. Compute Shannon and Simpson diversity from the species-abundance values.
6. Return the aggregated `TaxonomicProfile` record.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The implementation stores four abundance maps: kingdom, phylum, genus, and species. Empty rank strings are skipped when building phylum, genus, and species maps. To avoid division by zero, the code uses `1` as an internal denominator when `ClassifiedReads = 0`, but the corresponding count dictionaries are empty in that case, so the emitted abundance maps remain empty.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GenerateTaxonomicProfile` | `O(n)` | `O(t)` | `n` = number of classifications, `t` = number of unique taxa retained across stored ranks |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MetagenomicsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs)

- `MetagenomicsAnalyzer.GenerateTaxonomicProfile(IEnumerable<TaxonomicClassification>)`: Aggregates read-level classifications into a `TaxonomicProfile`.

### 5.2 Current Behavior

The current implementation filters reads on the `Kingdom` field only and then aggregates four stored ranks: kingdom, phylum, genus, and species. Although `TaxonomicClassification` includes class, order, and family fields, those ranks are not exposed in the `TaxonomicProfile` output record. Because empty phylum, genus, and species values are skipped before counting but the denominator remains `ClassifiedReads`, lower-rank abundance maps can sum to less than `1` when some retained reads are missing that rank. Diversity metrics are computed exclusively from species-level abundances, and Shannon uses the natural logarithm through the shared helper method.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Relative abundance as count divided by the total number of retained classified reads.
- Shannon diversity from species-level proportions.
- Simpson concentration from species-level proportions.

**Intentionally simplified:**

- The output profile stores kingdom, phylum, genus, and species abundances only; **consequence:** class, order, and family summaries must be computed elsewhere even though the upstream classification record contains those ranks.
- Diversity is computed only at the species level; **consequence:** the API does not expose analogous diversity metrics for higher taxonomic ranks.

**Not implemented:**

- (none)

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty input | `TotalReads = 0`, `ClassifiedReads = 0`, empty abundance maps, zero diversity | Early return through empty count dictionaries and zero diversity helpers |
| All reads unclassified | `ClassifiedReads = 0`, empty abundance maps, zero diversity | All records are removed by the kingdom filter |
| Single species after filtering | `ShannonDiversity = 0`, `SimpsonDiversity = 1` | One-category species distribution |
| Missing phylum, genus, or species values | Missing rank is excluded from that rank's abundance map, so that map can sum to less than `1` | Empty keys are filtered before each map is materialized while the denominator remains `ClassifiedReads` |

### 6.2 Limitations

The profile is only as informative as the upstream classifications and the stored rank set. Because the current output record omits class, order, and family abundance maps, this API is not a full rank-by-rank profile container. It also computes diversity only from species-level abundances rather than exposing diversity summaries for every rank.

## 7. Examples and Related Material

- [META-PROF-001](../../../tests/TestSpecs/META-PROF-001.md) documents the repository's taxonomic profile test specification.
- [Taxonomic_Classification.md](./Taxonomic_Classification.md) documents the upstream read-classification algorithm.
- [Alpha_Diversity.md](./Alpha_Diversity.md) documents the diversity metrics reused here.
- [Beta_Diversity.md](./Beta_Diversity.md) documents pairwise between-sample diversity metrics.

## 8. References

1. Shannon, C. E. 1948. A Mathematical Theory of Communication. Bell System Technical Journal.
2. Simpson, E. H. 1949. Measurement of Diversity. Nature.
3. Segata, N., et al. 2012. Metagenomic microbial community profiling using unique clade-specific marker genes. Nature Methods. doi:10.1038/nmeth.2066.
4. Wikipedia contributors. Metagenomics. Wikipedia. https://en.wikipedia.org/wiki/Metagenomics
