# Taxonomic Classification

## Overview

K-mer based taxonomic classification is a computational method for assigning taxonomic labels to DNA sequences, typically short reads from metagenomic studies. The approach uses exact k-mer matching against a pre-built database of reference genomes to rapidly classify sequences.

## Algorithm Description

### Core Principle

The algorithm extracts all k-mers (substrings of length k) from a query sequence and queries each against a database mapping k-mers to taxonomic identifiers. Classification is determined by the taxon with the most matching k-mers.

### Process

1. **Database Construction** (`BuildKmerDatabase`)
   - For each reference genome with known taxonomy:
     - Extract all k-mers of length k
     - Compute canonical form (lexicographically smaller of k-mer or reverse complement)
     - Map canonical k-mer → taxon ID

2. **Read Classification** (`ClassifyReads`)
   - For each query read:
     - Extract all k-mers of length k
     - Convert each to uppercase for case-insensitive matching
     - Query database for each k-mer
     - Count k-mer hits per taxon
     - Classify to taxon with highest hit count
     - Calculate confidence as MatchedKmers / TotalKmers

### Canonical K-mers

To reduce database size by ~50%, the algorithm uses canonical k-mers:

```
canonical(kmer) = min(kmer, reverse_complement(kmer))
```

This exploits the double-stranded nature of DNA where a sequence and its reverse complement represent the same biological entity.

## Complexity

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| BuildKmerDatabase | O(n × m) | O(k × u) |
| ClassifyReads | O(r × l) | O(1) per read |

Where:
- n = number of reference genomes
- m = average reference genome length
- k = k-mer length
- u = unique k-mers in database
- r = number of reads
- l = average read length

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| k | 31 | K-mer length. Longer k increases specificity, shorter k increases sensitivity |

### K-mer Length Selection

- **k=31**: Standard for metagenomic classification (Kraken default)
- **k<20**: Higher sensitivity, more false positives
- **k>35**: Higher specificity, may miss divergent sequences

## Output Structure

### TaxonomicClassification Record

| Field | Type | Description |
|-------|------|-------------|
| ReadId | string | Input read identifier |
| Kingdom | string | Taxonomic kingdom (or "Unclassified") |
| Phylum | string | Taxonomic phylum |
| Class | string | Taxonomic class |
| Order | string | Taxonomic order |
| Family | string | Taxonomic family |
| Genus | string | Taxonomic genus |
| Species | string | Taxonomic species |
| Confidence | double | MatchedKmers / TotalKmers (0–1) |
| MatchedKmers | int | Count of k-mers found in database |
| TotalKmers | int | Total k-mers extracted from read |

## Edge Cases

| Condition | Behavior |
|-----------|----------|
| Empty sequence | Returns Unclassified with Confidence=0 |
| Sequence length < k | Returns Unclassified (cannot extract k-mers) |
| No database matches | Returns Unclassified with TotalKmers recorded |
| Empty database | All reads return Unclassified |
| Non-ACGT characters | Skipped during k-mer extraction |

## Invariants

1. Output count equals input read count
2. 0 ≤ Confidence ≤ 1
3. MatchedKmers ≤ TotalKmers
4. TotalKmers = max(0, sequence_length - k + 1)
5. If MatchedKmers = 0, Kingdom = "Unclassified"

## Implementation Notes

### Taxonomy String Format

The implementation expects taxonomy strings in hierarchical format:
```
Kingdom|Phylum|Class|Order|Family|Genus|Species
```

Alternative delimiter: semicolon (`;`)

### Case Handling

All sequences are converted to uppercase internally for matching.

## References

1. Wood DE, Salzberg SL. Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 2014, 15:R46. DOI: 10.1186/gb-2014-15-3-r46

2. Wikipedia: Metagenomics — https://en.wikipedia.org/wiki/Metagenomics

3. Kraken Documentation — https://ccb.jhu.edu/software/kraken/

## See Also

- [META-PROF-001](./Taxonomic_Profile.md) — Taxonomic Profile Generation
- [META-ALPHA-001](./Alpha_Diversity.md) — Alpha Diversity Calculation
