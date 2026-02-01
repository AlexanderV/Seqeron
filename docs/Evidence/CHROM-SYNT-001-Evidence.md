# Evidence: CHROM-SYNT-001 - Synteny Analysis

**Test Unit ID:** CHROM-SYNT-001  
**Area:** Chromosome Analysis  
**Date:** 2026-02-01  
**Status:** Complete

---

## 1. Authoritative Sources

### 1.1 Primary Sources

| Source | Type | URL / Reference | Key Information |
|--------|------|-----------------|-----------------|
| Wikipedia - Synteny | Encyclopedia | https://en.wikipedia.org/wiki/Synteny | Definition, concepts, computational detection |
| Wikipedia - Comparative genomics | Encyclopedia | https://en.wikipedia.org/wiki/Comparative_genomics | Synteny blocks, rearrangement detection |
| Wikipedia - Chromosomal rearrangement | Encyclopedia | https://en.wikipedia.org/wiki/Chromosomal_rearrangement | Inversion, translocation, deletion, duplication |
| Wang et al. (2012) | Primary Literature | MCScanX: Nucleic Acids Res. 40(7):e49 | Algorithm for synteny detection |
| Goel et al. (2019) | Primary Literature | SyRI: Genome Biology 20:277 | Synteny and rearrangement identification |
| Liu et al. (2018) | Primary Literature | BMC Bioinformatics 19(1):26 | Systematic evaluation of synteny inference |

### 1.2 Tool References

| Tool | Description | Reference |
|------|-------------|-----------|
| MCScanX | Detection and evolutionary analysis of gene synteny and collinearity | Wang et al. (2012) |
| SyRI | Synteny and Rearrangement Identifier for whole-genome assemblies | Goel et al. (2019) |
| MUMmer | Whole-genome alignment suite | Marçais et al. (2018) |

---

## 2. Core Concepts

### 2.1 Synteny Definition

**From Wikipedia (Synteny):**
> In genomics, synteny more commonly refers to colinearity, i.e. conservation of blocks of order within two sets of chromosomes that are being compared with each other. These blocks are referred to as syntenic blocks.

**From Wikipedia (Comparative genomics):**
> Synteny blocks are more formally defined as regions of chromosomes between genomes that share a common order of homologous genes derived from a common ancestor.

### 2.2 Synteny Block Properties

Based on literature:
- **Collinearity:** Genes maintain relative order between species
- **Strand orientation:** '+' (forward/same) or '-' (inverted/opposite)
- **Gene count:** Minimum number of genes required to define a block
- **Gap tolerance:** Maximum allowed gap between consecutive ortholog pairs
- **Sequence identity:** Conservation level within the block

### 2.3 Chromosomal Rearrangement Types

**From Wikipedia (Chromosomal rearrangement):**
- **Inversion:** Segment is reversed in orientation (detected by strand change within same chromosome)
- **Translocation:** Segment moves to a different chromosome
- **Deletion:** Segment is removed
- **Duplication:** Segment is copied

---

## 3. Algorithm Characteristics

### 3.1 FindSyntenyBlocks Algorithm

**Input:**
- List of ortholog pairs with positions in both genomes
- `minGenes`: Minimum genes to form a block (typical default: 3-5)
- `maxGap`: Maximum gap between consecutive genes (in megabases)

**Process:**
1. Group ortholog pairs by chromosome pairs
2. Sort by position in reference genome
3. Identify collinear runs (same relative order)
4. Merge consecutive collinear segments respecting gap constraints
5. Output blocks meeting minimum gene threshold

**Output:**
- Synteny blocks with coordinates, strand, gene count, identity

### 3.2 DetectRearrangements Algorithm

**Input:**
- List of synteny blocks

**Process:**
1. Sort blocks by reference chromosome and position
2. Compare adjacent blocks:
   - Different target chromosome → Translocation
   - Same target chromosome, different strand → Inversion
3. Output detected rearrangement events

**Output:**
- Rearrangement events with type, position, size

---

## 4. Edge Cases and Corner Cases

### 4.1 FindSyntenyBlocks

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty input | Return empty | Implementation |
| Fewer genes than minGenes | Return empty | Definition |
| All genes collinear (forward) | Single block with '+' strand | Definition |
| All genes collinear (reverse) | Single block with '-' strand | Definition |
| Gap exceeds maxGap | Break into separate blocks | Definition |
| Multiple chromosome pairs | Separate blocks per pair | Definition |
| Single gene | Return empty (below minGenes default) | Definition |
| Two genes (minGenes=3) | Return empty | Definition |

### 4.2 DetectRearrangements

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty input | Return empty | Implementation |
| Single block | Return empty (no adjacent pairs) | Definition |
| All same chromosome, same strand | Return empty (no rearrangements) | Definition |
| Strand change (same chromosome) | Inversion detected | Wikipedia |
| Chromosome change | Translocation detected | Wikipedia |

---

## 5. Test Data Sets

### 5.1 Collinear Forward Block

```
Ortholog pairs (all on chr1 → chrA, forward order):
Gene1: (chr1, 1000-2000) → (chrA, 1000-2000)
Gene2: (chr1, 3000-4000) → (chrA, 3000-4000)
Gene3: (chr1, 5000-6000) → (chrA, 5000-6000)
Gene4: (chr1, 7000-8000) → (chrA, 7000-8000)

Expected: 1 block, strand '+', GeneCount ≥ 3
```

### 5.2 Inverted Block

```
Ortholog pairs (reverse order in target):
Gene1: (chr1, 1000-2000) → (chrA, 8000-9000)
Gene2: (chr1, 3000-4000) → (chrA, 6000-7000)
Gene3: (chr1, 5000-6000) → (chrA, 4000-5000)
Gene4: (chr1, 7000-8000) → (chrA, 2000-3000)

Expected: 1 block, strand '-'
```

### 5.3 Translocation Detection

```
Synteny blocks:
Block1: chr1:1000-50000 → chrA:1000-50000, strand '+'
Block2: chr1:60000-100000 → chrB:1000-40000, strand '+'

Expected: Translocation detected (chrA → chrB)
```

### 5.4 Inversion Detection

```
Synteny blocks:
Block1: chr1:1000-50000 → chrA:1000-50000, strand '+'
Block2: chr1:60000-100000 → chrA:60000-100000, strand '-'

Expected: Inversion detected (strand change on same chromosome)
```

---

## 6. Invariants

### 6.1 FindSyntenyBlocks Invariants

- **I1:** All returned blocks have GeneCount ≥ minGenes
- **I2:** Block coordinates are valid (Start ≤ End) for both species
- **I3:** Strand is either '+' or '-'
- **I4:** SequenceIdentity is in range [0, 1]
- **I5:** All genes in a block belong to the same chromosome pair

### 6.2 DetectRearrangements Invariants

- **I1:** Rearrangement Type is a recognized value ("Inversion", "Translocation")
- **I2:** Position1 is always set (non-null)
- **I3:** For translocations, Chromosome2 differs from source block's target chromosome

---

## 7. Known Limitations

1. **Gap parameter scale:** Implementation uses `maxGap * 1000000` suggesting megabase units
2. **Placeholder identity:** Current implementation returns fixed 0.9 identity (not calculated)
3. **Deletion/Duplication not detected:** Only inversions and translocations are identified
4. **Single chromosome pair per block:** Breaks at chromosome boundaries

---

## 8. References

1. Wang Y, Tang H, et al. (2012). MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. Nucleic Acids Research. 40(7):e49. doi:10.1093/nar/gkr1293

2. Goel M, Sun H, et al. (2019). SyRI: Finding genomic rearrangements and local sequence differences from whole-genome assemblies. Genome Biology. 20:277. doi:10.1186/s13059-019-1911-0

3. Liu D, Hunt M, Tsai IJ (2018). Inferring synteny between genome assemblies: a systematic evaluation. BMC Bioinformatics. 19(1):26. doi:10.1186/s12859-018-2026-4

4. Wikipedia contributors. Synteny. Wikipedia, The Free Encyclopedia. Accessed 2026-02-01.

5. Wikipedia contributors. Comparative genomics. Wikipedia, The Free Encyclopedia. Accessed 2026-02-01.

6. Wikipedia contributors. Chromosomal rearrangement. Wikipedia, The Free Encyclopedia. Accessed 2026-02-01.
