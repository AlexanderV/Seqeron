# Centromere Analysis

## Overview

Centromere analysis involves identifying and characterizing the centromeric region of chromosomes. The centromere is a specialized DNA sequence that links sister chromatids during cell division and serves as the attachment point for spindle fibers via the kinetochore.

## Biological Background

### Centromere Structure

The centromere creates two chromosome arms:
- **p arm** (short arm): Named from French "petit" (small)
- **q arm** (long arm): Named as q follows p in the alphabet

### Classification by Position

Based on Levan et al. (1964) nomenclature:

| Classification | Arm Ratio (p/q) | Description |
|----------------|-----------------|-------------|
| Metacentric    | 1.0 - 1.7       | Arms approximately equal length |
| Submetacentric | 1.7 - 3.0       | Arms close but unequal |
| Acrocentric    | 3.0 - 7.0       | One arm much shorter |
| Telocentric    | > 7.0           | Centromere at or very near end |

### Molecular Characteristics

Human centromeres are characterized by:
- **Alpha-satellite DNA**: ~171 bp tandem repeats (alphoid sequences)
- **High repeat content**: Repetitive DNA organized in higher-order repeat units
- **Heterochromatin**: Constitutive heterochromatin packaging
- **Low GC variability**: More uniform GC content compared to gene-rich regions

## Algorithm Description

### Implementation: `ChromosomeAnalyzer.AnalyzeCentromere`

The algorithm uses a heuristic sliding-window approach:

1. **Window Scanning**: Scan sequence with overlapping windows (default 100kb)
2. **Repeat Content Estimation**: Calculate k-mer (k=15) frequency to estimate repetitiveness
3. **GC Variability**: Measure variance in GC content across sub-windows
4. **Scoring**: Score = RepeatContent × (1 - GCVariability)
5. **Boundary Extension**: Extend detected region while repeat content remains high
6. **Classification**: Determine type based on position relative to chromosome length

### Position-Based Classification Logic

```
Position ratio = centromere_midpoint / chromosome_length

< 0.15 → Acrocentric
0.15 - 0.35 → Submetacentric
0.35 - 0.65 → Metacentric
0.65 - 0.85 → Submetacentric
> 0.85 → Acrocentric
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chromosomeName` | string | - | Identifier for the chromosome |
| `sequence` | string | - | DNA sequence to analyze |
| `windowSize` | int | 100,000 | Sliding window size in bp |
| `minAlphaSatelliteContent` | double | 0.3 | Minimum repeat content threshold |

## Output

```csharp
public readonly record struct CentromereResult(
    string Chromosome,
    int? Start,
    int? End,
    int Length,
    string CentromereType,
    double AlphaSatelliteContent,
    bool IsAcrocentric);
```

## Complexity

- **Time**: O(n) where n = sequence length
- **Space**: O(k) for k-mer dictionary within each window

## Limitations

1. **Heuristic approach**: Does not use reference alpha-satellite databases
2. **Single candidate**: Returns only the best-scoring region
3. **Synthetic sequences**: May produce false positives on artificial repetitive sequences

## References

1. Wikipedia - Centromere: https://en.wikipedia.org/wiki/Centromere
2. Wikipedia - Karyotype: https://en.wikipedia.org/wiki/Karyotype
3. Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes". Hereditas. 52(2): 201-220.
4. Mehta GD, Agarwal MP, Ghosh SK (2010). "Centromere identity: a challenge to be faced". Molecular Genetics and Genomics. 284(2): 75-94.

## Implementation Notes

- The implementation uses `EstimateRepeatContent` with 15-mer analysis
- GC variability is calculated over 1kb sub-windows
- Boundary extension uses 70% of the threshold to allow gradual transition zones
- The algorithm is designed for computational analysis, not clinical diagnostics
