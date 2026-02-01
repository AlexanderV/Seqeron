# Telomere Analysis

## Overview

Telomere analysis detects and measures telomeric repeat sequences at chromosome ends. Telomeres are protective caps consisting of tandem repeats that prevent chromosome degradation and fusion.

## Biological Background

### Definition
A **telomere** is a region of repetitive nucleotide sequences at chromosome ends that protects coding DNA from progressive degradation during replication.

### Structure
- **Vertebrate canonical repeat:** TTAGGG (5'→3')
- **Reverse complement:** CCCTAA (found at 5' chromosome end)
- **Repeat unit:** 6 base pairs
- **Human telomere length:** 5,000–15,000 bp at birth

### Chromosome End Orientation
| End | Repeat Sequence | Direction |
|-----|-----------------|-----------|
| 5' (start) | CCCTAA | Toward chromosome interior |
| 3' (terminus) | TTAGGG | Toward chromosome end |

**Source:** Wikipedia (Telomere), Meyne et al. (1989)

---

## Algorithm: AnalyzeTelomeres

### Purpose
Detect telomeric repeats at both chromosome ends and measure their length and purity.

### Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| chromosomeName | string | required | Chromosome identifier |
| sequence | string | required | DNA sequence to analyze |
| telomereRepeat | string | "TTAGGG" | Canonical repeat unit |
| searchLength | int | 10000 | Maximum distance from ends to search |
| minTelomereLength | int | 500 | Minimum length to report as present |
| criticalLength | int | 3000 | Threshold for critically short flag |

### Return Value: TelomereResult
| Field | Type | Description |
|-------|------|-------------|
| Chromosome | string | Chromosome name |
| Has5PrimeTelomere | bool | Whether 5' telomere meets minimum threshold |
| TelomereLength5Prime | int | Measured length at 5' end |
| Has3PrimeTelomere | bool | Whether 3' telomere meets minimum threshold |
| TelomereLength3Prime | int | Measured length at 3' end |
| RepeatPurity5Prime | double | Fraction of matching bases at 5' |
| RepeatPurity3Prime | double | Fraction of matching bases at 3' |
| IsCriticallyShort | bool | True if any telomere below critical threshold |

### Algorithm Steps
1. Convert sequence and repeat to uppercase
2. Compute reverse complement of telomere repeat (CCCTAA for TTAGGG)
3. Search 5' end region for reverse complement repeats
4. Search 3' end region for forward repeats
5. For each end:
   - Iterate through repeat-sized windows
   - Calculate similarity to expected repeat
   - Continue while similarity ≥ 70%
   - Accumulate length and matching bases
6. Calculate purity = matching bases / total bases
7. Determine flags based on thresholds

### Complexity
**Time:** O(n) where n = min(searchLength, sequence.Length)  
**Space:** O(1) auxiliary

---

## Algorithm: EstimateTelomereLengthFromTSRatio

### Purpose
Convert qPCR T/S ratio to estimated telomere length in base pairs.

### Background
The T/S ratio method (Cawthon, 2002) measures the ratio of telomere repeat copy number to a single-copy gene. This ratio is proportional to average telomere length.

**Source:** Cawthon (2002), DOI: 10.1093/nar/30.10.e47

### Formula
```
estimatedLength = referenceLength × (tsRatio / referenceRatio)
```

### Parameters
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| tsRatio | double | required | Measured T/S ratio |
| referenceRatio | double | 1.0 | Reference sample T/S ratio |
| referenceLength | double | 7000 | Reference sample telomere length (bp) |

### Return Value
Estimated telomere length in base pairs (double).

### Complexity
**Time:** O(1)  
**Space:** O(1)

---

## Implementation Notes

### Current Implementation
Located in `ChromosomeAnalyzer.cs` (lines 230–340).

### Repeat Matching Tolerance
The implementation uses 70% similarity threshold, allowing:
- ~2 mismatches per 6-bp repeat
- Detection of aged/divergent telomeres
- Robustness to sequencing errors

### Edge Cases Handled
1. **Empty sequence:** Returns no telomere, critically short = true
2. **Sequence shorter than repeat:** Returns no telomere
3. **No telomeric repeats:** Returns lengths = 0, Has*Telomere = false
4. **Partial repeats:** Only complete repeat units counted

---

## Species-Specific Repeats

| Organism | Repeat | Notes |
|----------|--------|-------|
| Vertebrates | TTAGGG | Conserved across all vertebrates |
| Arabidopsis | TTTAGGG | 7-bp repeat |
| Tetrahymena | TTGGGG | Original discovery organism |
| S. cerevisiae | Variable | Irregular repeats |

**Source:** Wikipedia (Telomere), TeloBase database

---

## Clinical Significance

- **Normal human range:** 5,000–15,000 bp
- **Critical threshold:** ~3,000 bp (triggers senescence)
- **Association:** Shorter telomeres correlate with aging and disease risk

**Source:** Rossiello et al. (2022), Nature Cell Biology

---

## References

1. Wikipedia - Telomere. https://en.wikipedia.org/wiki/Telomere
2. Meyne J, Ratliff RL, Moyzis RK (1989). Conservation of the human telomere sequence (TTAGGG)n among vertebrates. PNAS 86(18):7049-53.
3. Cawthon RM (2002). Telomere measurement by quantitative PCR. Nucleic Acids Res 30(10):e47.
4. Blackburn EH, Gall JG (1978). A tandemly repeated sequence at the termini of the extrachromosomal ribosomal RNA genes in Tetrahymena. J Mol Biol 120(1):33-53.
