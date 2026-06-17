# Evidence: META-CLASS-001 — Taxonomic Classification

## Overview

| Field | Value |
|-------|-------|
| Test Unit ID | META-CLASS-001 |
| Algorithm | K-mer Based Taxonomic Classification |
| Area | Metagenomics |
| Methods | `ClassifyReads`, `BuildKmerDatabase`, `TaxonomyTree` (+ `Lca`) |
| Date | 2026-06-17 |

> **Classifier scope:** This unit implements the **faithful Kraken algorithm** (Wood & Salzberg 2014). `BuildKmerDatabase` maps each canonical k-mer to the **lowest common ancestor (LCA) of all reference taxa that contain it** ("a k-mer and the LCA of all organisms whose genomes contain that k-mer"; subsequent owners fold via LCA). `ClassifyReads` builds the per-read **classification tree** (hit taxa + ancestors, weighted by k-mer count), assigns the leaf of the **maximum-scoring root-to-leaf (RTL) path**, and breaks ties by the **LCA of the maximally-scoring leaves**; reads with no hits are **unclassified** (taxonomy root). Confidence is Kraken 2's **C/Q** (C = k-mers in the clade rooted at the assigned label, Q = non-ambiguous k-mers queried). `TaxonomyTree.Lca` provides the LCA used in both DB build and tie-break. (Prior to enhancement C1 this was a flat best-hit classifier with no LCA; that wording is now obsolete.)

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

2. **Kraken 1 Manual (CCB Johns Hopkins University)**
   - URL: https://ccb.jhu.edu/software/kraken/MANUAL.html
   - Relevance: Canonical implementation of k-mer based metagenomic classification
   - Key excerpts:
     - "Kraken is a system for assigning taxonomic labels to short DNA sequences"
     - "Kraken does not query k-mers containing ambiguous nucleotides (non-ACGT)"
     - Uses k=31 as default k-mer size
     - Canonical k-mers used in both database building AND querying
     - Confidence score: C/Q where C = clade k-mers, Q = non-ambiguous k-mers queried

3. **Kraken 2 Manual (GitHub Wiki)**
   - URL: https://github.com/DerrickWood/kraken2/wiki/Manual
   - Relevance: Updated Kraken implementation confirming scoring formula
   - Key excerpt:
     - "A sequence label's score is a fraction C/Q, where C is the number of k-mers mapped to LCA values in the clade rooted at the label, and Q is the number of k-mers in the sequence that lack an ambiguous nucleotide (i.e., they were queried against the database)."

4. **Wood & Salzberg (2014)** — Primary Reference
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
| Confidence | C / Q | C = k-mers in the clade rooted at the assigned label; Q = non-ambiguous k-mers queried |
| Ambiguous filtering | Yes | K-mers with non-ACGT characters skipped during classification |

### Classification Algorithm (Kraken — Wood & Salzberg 2014)

**Database build.** For each labeled reference, set each contained canonical k-mer's stored taxon to
the **LCA of the previously-stored taxon (if any) and this reference's taxon** — "if a k-mer … has
had its LCA value previously set, then the LCA of the stored value and the current sequence's taxon
is calculated."

**Per-read classification:**
1. Extract all k-mers from a read; skip those with ambiguous (non-ACGT) nucleotides (not counted in Q).
2. Canonicalize each k-mer: `min(kmer, reverse_complement(kmer))`; query against the database.
3. Tally per-taxon k-mer hit counts; if none, the read is **unclassified** (root).
4. The hit taxa + their ancestors form the **classification tree**; each node is weighted by its
   k-mer count.
5. Score every **root-to-leaf (RTL)** path as the sum of node weights along it; the **maximum-scoring
   path** is the classification path and its **leaf** is the assigned taxon.
6. Tie-break: if several paths share the maximum score, assign the **LCA of their leaves**.
7. Confidence = **C/Q**, C = k-mers mapped to a taxon in the clade rooted at the assigned label, Q =
   non-ambiguous k-mers queried.

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
| Multiple taxon matches | Kraken RTL / LCA | Assign the leaf of the max-scoring root-to-leaf path; equal-score paths → LCA of their leaves |
| K-mer shared by several taxa (DB) | Kraken LCA | Stored as the LCA of the owning taxa |
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

3. **Confidence calculation**: `C / Q` per Kraken — C = k-mers in the clade rooted at the assigned label, Q = non-ambiguous k-mers

4. **Ambiguous k-mer filtering**: K-mers with non-ACGT characters skipped (not counted in TotalKmers)

5. **Canonical lookup in ClassifyReads**: Read k-mers canonicalized before database lookup

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
3. **K-mer count invariant**: TotalKmers (Q) = count of non-ambiguous k-mers (equals len - k + 1 for all-ACGT sequences)
4. **Matched k-mers bound**: MatchedKmers (C) ≤ TotalKmers (Q); C = k-mers in the clade rooted at the assigned label
5. **Unclassified criteria**: no k-mer hits ⇒ TaxonId = root (`TaxonomyTree.RootId`), C = 0
6. **Canonical k-mer uniqueness**: For any k-mer, exactly one canonical form exists
7. **Confidence formula**: Confidence = C / Q (Kraken 2)
8. **DB-build LCA**: a canonical k-mer owned by several taxa is stored as their LCA
9. **RTL assignment**: assigned taxon = leaf of the maximum-scoring root-to-leaf path; ties → LCA of tied leaves
10. **LCA correctness**: siblings → parent; ancestor/descendant → ancestor; self → self; disjoint → root

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
