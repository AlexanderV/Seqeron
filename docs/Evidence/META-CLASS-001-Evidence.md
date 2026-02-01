# Evidence: META-CLASS-001 — Taxonomic Classification

## Overview

| Field | Value |
|-------|-------|
| Test Unit ID | META-CLASS-001 |
| Algorithm | K-mer Based Taxonomic Classification |
| Area | Metagenomics |
| Methods | `ClassifyReads`, `BuildKmerDatabase` |
| Date | 2026-02-01 |

---

## Authoritative Sources

### Primary Sources

1. **Wikipedia — Metagenomics**
   - URL: https://en.wikipedia.org/wiki/Metagenomics
   - Relevance: Defines metagenomics, k-mer based classification, and taxonomic profiling
   - Key excerpts:
     - "Binning is the process of associating a particular sequence with an organism"
     - "Similarity-based binning methods such as BLAST are used to rapidly search for phylogenetic markers"
     - "Composition based binning methods use intrinsic features of the sequence, such as oligonucleotide frequencies"

2. **Kraken Documentation (CCB Johns Hopkins University)**
   - URL: https://ccb.jhu.edu/software/kraken/
   - Relevance: Canonical implementation of k-mer based metagenomic classification
   - Key excerpts:
     - "Kraken is a system for assigning taxonomic labels to short DNA sequences"
     - "Kraken aims to achieve high sensitivity and high speed by utilizing exact alignments of k-mers"
     - Uses k=31 as default k-mer size
     - Canonical k-mers used to reduce database size (forward or reverse complement, whichever is lexicographically smaller)

3. **Wood & Salzberg (2014)** — Primary Reference
   - Citation: Wood DE, Salzberg SL. Kraken: ultrafast metagenomic sequence classification using exact alignments. Genome Biology 2014, 15:R46.
   - DOI: 10.1186/gb-2014-15-3-r46
   - Relevance: Foundational paper describing k-mer classification algorithm

---

## Algorithm Characteristics

### From Kraken Documentation

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| Default k | 31 | Standard for genomic classification |
| Classification method | Exact k-mer matching | Maps k-mers to taxonomy database |
| Canonical k-mers | Yes | Min(kmer, reverse_complement(kmer)) lexicographically |
| Confidence | MatchedKmers / TotalKmers | Fraction of k-mers with database hits |

### Classification Algorithm

Per Kraken/Wikipedia:
1. Extract all k-mers from a read
2. Query each k-mer against the database
3. Count k-mer hits per taxon
4. Classify to taxon with most k-mer hits
5. Reads with no hits → "Unclassified"

---

## Test Datasets from Sources

### Kraken Accuracy Dataset
- Source: https://ccb.jhu.edu/software/kraken/dl/accuracy.tgz
- Contains: Three FASTA files with 10,000 simulated reads each
- Read length: 100 bp
- Error rate: 2.1% SNP, 1.1% indel

### Documented Test Cases

| Case | Source | Expected Behavior |
|------|--------|-------------------|
| Empty sequence | Implementation standard | Return "Unclassified", no error |
| Sequence shorter than k | Kraken manual | Cannot extract k-mers, "Unclassified" |
| No matching k-mers | Kraken | Return "Unclassified" |
| Multiple taxon matches | Kraken (LCA) | Classify to highest-count taxon |
| All-N sequence | DUST filtering docs | Handle gracefully |

---

## Edge Cases and Corner Cases

### From Literature

1. **Low-complexity sequences** (Morgulis et al. 2006)
   - Runs of single nucleotides (e.g., AAAA...)
   - Can cause false positive hits
   - Should use DUST filtering in production

2. **Reads shorter than k**
   - Cannot produce any k-mers
   - Must return Unclassified

3. **Canonical k-mer normalization**
   - Forward: ATGCGATCGATCGA
   - Reverse complement: TCGATCGATCGCAT
   - Use lexicographically smaller

4. **Empty input**
   - Empty read list → empty output
   - Empty database → all reads Unclassified

---

## Implementation-Specific Notes

### Current Seqeron.Genomics Implementation

From `MetagenomicsAnalyzer.cs`:

1. **Taxonomy string format**: Semicolon or pipe-delimited hierarchy
   - Example: "Bacteria|Proteobacteria|Gamma|Escherichia|coli"
   - Parsed into kingdom, phylum, class, order, family, genus, species

2. **Canonical k-mer**: Uses `DnaSequence.GetReverseComplementString()` for reverse complement

3. **Confidence calculation**: `MatchedKmers / TotalKmers`

4. **Case handling**: Converts sequences to uppercase internally

---

## Testing Methodology

### From Literature (Kraken accuracy evaluation)

1. **Precision/Sensitivity at genus level**
   - Measure correct classifications vs total classifications
   - Measure correct classifications vs total reads

2. **Simulated metagenomes**
   - Known composition for ground truth
   - Variable error rates

### Required Test Categories

| Category | Justification | Source |
|----------|---------------|--------|
| Empty/null inputs | Robustness | Standard |
| Short sequences | Edge case | Kraken manual |
| No matches | Expected behavior | Kraken |
| Exact matches | Core functionality | Kraken |
| Multiple reads | Batch processing | Kraken |
| Canonical k-mers | Correctness | Kraken paper |
| Mixed case input | Robustness | Standard |

---

## Invariants

1. **Output count invariant**: |output| = |input reads|
2. **Confidence range**: 0 ≤ Confidence ≤ 1
3. **K-mer count invariant**: TotalKmers = max(0, len(sequence) - k + 1)
4. **Matched k-mers bound**: MatchedKmers ≤ TotalKmers
5. **Unclassified criteria**: If MatchedKmers = 0, Kingdom = "Unclassified"
6. **Canonical k-mer uniqueness**: For any k-mer, exactly one canonical form exists

---

## Open Questions

None identified. Algorithm behavior is well-documented by Kraken.

---

## Summary

Evidence is sufficient to define comprehensive tests for META-CLASS-001. The k-mer based classification algorithm is well-established with canonical reference (Kraken). Key test areas:
- Input validation (empty, short, null)
- Core classification logic
- Canonical k-mer handling
- Confidence calculation
- Batch processing invariants
