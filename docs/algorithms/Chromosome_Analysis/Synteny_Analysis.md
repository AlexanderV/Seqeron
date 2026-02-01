# Synteny Analysis

**Algorithm Group:** Chromosome Analysis  
**Implementation:** `ChromosomeAnalyzer.FindSyntenyBlocks`, `ChromosomeAnalyzer.DetectRearrangements`  
**Date:** 2026-02-01

---

## 1. Overview

Synteny analysis identifies conserved blocks of gene order between two genomes and detects chromosomal rearrangements that have occurred during evolution. This is a fundamental technique in comparative genomics for understanding genome evolution and establishing orthology relationships.

---

## 2. Biological Context

### 2.1 Synteny Definition

**Synteny** (Greek: σύν "together" + ταινία "band") refers to the physical co-localization of genetic loci on the same chromosome. In modern genomics, the term more commonly refers to **collinearity**—conservation of gene order between chromosomes of different species.

**Synteny blocks** are regions where multiple genes maintain their relative order and orientation, indicating descent from a common ancestor without major intervening rearrangements.

### 2.2 Biological Significance

- **Evolutionary studies:** Reveals genome organization changes during speciation
- **Orthology assignment:** Conserved synteny supports gene homology predictions
- **Functional genomics:** Syntenic regions often share regulatory mechanisms
- **Comparative mapping:** Enables cross-species gene prediction

### 2.3 Chromosomal Rearrangements

Four major types of chromosomal rearrangements break synteny:

| Type | Description | Detection Signature |
|------|-------------|---------------------|
| **Inversion** | Segment reversed in orientation | Strand change within same chromosome |
| **Translocation** | Segment moved to different chromosome | Chromosome change between adjacent blocks |
| **Deletion** | Segment removed | Gap in synteny coverage |
| **Duplication** | Segment copied | Multiple synteny hits |

---

## 3. Algorithm Description

### 3.1 FindSyntenyBlocks

#### Input
- **orthologPairs:** Collection of ortholog gene pairs with coordinates in both genomes
  - Format: `(Chr1, Start1, End1, Gene1, Chr2, Start2, End2, Gene2)`
- **minGenes:** Minimum number of genes to form a valid block (default: 3)
- **maxGap:** Maximum gap between consecutive genes (in megabases, default: 10)

#### Algorithm

```
1. Group ortholog pairs by chromosome pair (Chr1, Chr2)
2. For each chromosome pair group:
   a. Sort pairs by position in reference genome (Start1)
   b. Initialize block tracking variables
   c. For each consecutive pair:
      - Determine direction (forward if pos2 increasing, reverse otherwise)
      - Check if collinear (same direction as current block)
      - Check if gap within maxGap threshold
      - If collinear and within gap: extend current block
      - Otherwise: emit current block (if ≥ minGenes), start new block
   d. Emit final block if ≥ minGenes
3. Yield all valid synteny blocks
```

#### Output
- **SyntenyBlock** records containing:
  - Species1Chromosome, Species1Start, Species1End
  - Species2Chromosome, Species2Start, Species2End
  - Strand ('+' or '-')
  - GeneCount
  - SequenceIdentity

#### Complexity
- Time: O(n log n) for sorting, O(n) for scanning
- Space: O(n) for storing pairs

### 3.2 DetectRearrangements

#### Input
- **syntenyBlocks:** Collection of synteny blocks from FindSyntenyBlocks

#### Algorithm

```
1. Sort blocks by Species1Chromosome, then Species1Start
2. For each consecutive pair of blocks (current, next):
   a. If same Species1 chromosome:
      - If same Species2 chromosome but different Strand:
        → Emit Inversion event
      - If different Species2 chromosome:
        → Emit Translocation event
3. Yield all detected rearrangement events
```

#### Output
- **ChromosomalRearrangement** records containing:
  - Type ("Inversion" or "Translocation")
  - Chromosome1, Position1 (location in reference)
  - Chromosome2, Position2 (for translocations)
  - Size (optional)
  - Description

#### Complexity
- Time: O(n log n) for sorting, O(n) for detection
- Space: O(n) for storing blocks

---

## 4. Implementation Notes

### 4.1 Current Implementation

The `Seqeron.Genomics.ChromosomeAnalyzer` implementation:

1. **Gap scaling:** The maxGap parameter is multiplied by 1,000,000 (megabase units)
2. **Identity placeholder:** Returns fixed 0.9 identity (actual calculation not implemented)
3. **Rearrangement types:** Only Inversion and Translocation are detected
4. **Block boundaries:** Coordinates span from first to last gene in block

### 4.2 Data Structures

```csharp
public readonly record struct SyntenyBlock(
    string Species1Chromosome,
    int Species1Start,
    int Species1End,
    string Species2Chromosome,
    int Species2Start,
    int Species2End,
    char Strand,
    int GeneCount,
    double SequenceIdentity);

public readonly record struct ChromosomalRearrangement(
    string Type,
    string Chromosome1,
    int Position1,
    string? Chromosome2,
    int? Position2,
    int? Size,
    string? Description);
```

### 4.3 Limitations

1. **Deletion/duplication detection:** Not implemented in current version
2. **Sequence identity:** Not calculated from actual sequences
3. **Complex rearrangements:** Only pairwise block comparisons performed

---

## 5. Usage Examples

### 5.1 Finding Synteny Blocks

```csharp
var orthologPairs = new List<(string, int, int, string, string, int, int, string)>
{
    ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
    ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
    ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
};

var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
    orthologPairs, 
    minGenes: 3, 
    maxGap: 10);

foreach (var block in blocks)
{
    Console.WriteLine($"{block.Species1Chromosome}:{block.Species1Start}-{block.Species1End} " +
                      $"-> {block.Species2Chromosome}:{block.Species2Start}-{block.Species2End} " +
                      $"Strand: {block.Strand}, Genes: {block.GeneCount}");
}
```

### 5.2 Detecting Rearrangements

```csharp
var syntenyBlocks = ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs);
var rearrangements = ChromosomeAnalyzer.DetectRearrangements(syntenyBlocks);

foreach (var r in rearrangements)
{
    Console.WriteLine($"{r.Type} at {r.Chromosome1}:{r.Position1}");
}
```

---

## 6. References

1. Wang Y, Tang H, et al. (2012). MCScanX: a toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research*. 40(7):e49.

2. Goel M, Sun H, et al. (2019). SyRI: Finding genomic rearrangements and local sequence differences from whole-genome assemblies. *Genome Biology*. 20:277.

3. Liu D, Hunt M, Tsai IJ (2018). Inferring synteny between genome assemblies: a systematic evaluation. *BMC Bioinformatics*. 19(1):26.

4. Wikipedia contributors. "Synteny." *Wikipedia, The Free Encyclopedia*.

5. Wikipedia contributors. "Comparative genomics." *Wikipedia, The Free Encyclopedia*.

6. Wikipedia contributors. "Chromosomal rearrangement." *Wikipedia, The Free Encyclopedia*.
