# Probe Validation

## Overview

Probe validation is a computational process for assessing the quality and specificity of hybridization probes. It ensures that probes will bind specifically to their intended target sequences with minimal off-target hybridization. This is critical for applications such as FISH, DNA microarrays, qPCR, and other molecular biology techniques.

## Algorithm Description

### Core Algorithm: `ValidateProbe`

The `ValidateProbe` method assesses a probe sequence against reference sequences:

1. **Normalize Input**
   - Convert probe sequence to uppercase for consistent comparison

2. **Off-Target Analysis**
   - Search each reference sequence for approximate matches within tolerance
   - Count total number of matching sites across all references
   - Multiple hits indicate potential cross-hybridization issues

3. **Self-Complementarity Check**
   - Calculate the degree of self-complementarity
   - High self-complementarity (>30%) may form hairpins or dimers

4. **Secondary Structure Detection**
   - Check for potential hairpin formation
   - Secondary structures can reduce hybridization efficiency

5. **Calculate Specificity Score**
   - Unique match (1 hit): specificity = 1.0
   - Multiple hits: specificity = 1.0 / hitCount
   - No match: specificity = 0.0 (probe doesn't match target)

6. **Return Validation Result**
   - IsValid: true if issues count is 0 or (offTargetHits ≤ 1 and selfComp ≤ 0.4)
   - SpecificityScore: calculated as above
   - OffTargetHits: total count of matching sites
   - SelfComplementarity: fraction of self-complementary bases
   - HasSecondaryStructure: boolean indicating hairpin potential
   - Issues: list of detected problems

### Complexity

- **Time**: O(n × g × m) where n = probe length, g = reference count, m = reference lengths
- **Space**: O(1) for validation result (excluding input sequences)

### CheckSpecificity Algorithm (Suffix Tree Optimization)

The `CheckSpecificity` method provides fast specificity checking using a pre-built suffix tree:

```csharp
// O(m) specificity check per probe using suffix tree
var positions = genomeIndex.FindAllOccurrences(probeSequence);
return positions.Count == 1 ? 1.0 : 1.0 / positions.Count;
```

This provides O(m) lookup time instead of O(n × m) linear search, where m = probe length.

## Scientific Background

### Cross-Hybridization

Cross-hybridization occurs when a probe binds to sequences other than the intended target. This is a major concern in probe design because:

- **Mismatch Tolerance**: Hybridization can occur with 1-5 base pair mismatches depending on conditions
- **Stringency**: Higher stringency (temperature, salt) reduces non-specific binding but may also reduce specific binding
- **Probe Length**: Longer probes have higher specificity but may have more off-target sites

Source: Wikipedia (Hybridization probe), Wikipedia (DNA microarray)

### Specificity Factors

| Factor | Effect on Specificity | Source |
|--------|----------------------|--------|
| Unique sequence | High specificity (score = 1.0) | Implementation |
| Multiple genome hits | Reduced specificity (score = 1/hits) | Implementation |
| Self-complementarity >30% | May form hairpins | General practice |
| Low-complexity regions | Higher cross-hybridization risk | Wikipedia (DNA microarray) |

### Off-Target Detection Methods

Modern off-target detection approaches include:

1. **Approximate Matching**: Allow mismatches within a threshold
2. **BLAST-like algorithms**: For database searching (O(n) heuristic)
3. **Suffix tree indexing**: For O(m) exact or approximate matching

Source: Wikipedia (BLAST), Altschul et al. (1990)

### Validation Parameters

| Parameter | Default | Description | Source |
|-----------|---------|-------------|--------|
| maxMismatches | 3 | Maximum allowed mismatches for off-target detection | Implementation (aligns with CRISPR off-target tolerance) |
| Self-complementarity threshold | 0.3 | Maximum acceptable self-complementarity | Implementation |
| Secondary structure check | true | Whether to check for hairpin potential | Implementation |

## Implementation Notes

### ValidateProbe Method

```csharp
public static ProbeValidation ValidateProbe(
    string probeSequence,
    IEnumerable<string> referenceSequences,
    int maxMismatches = 3)
```

**Parameters:**
- `probeSequence`: The probe sequence to validate
- `referenceSequences`: Reference sequences to search for off-target sites
- `maxMismatches`: Maximum mismatches for approximate matching (default: 3)

**Returns:** `ProbeValidation` record with validation results

### CheckSpecificity Method

```csharp
public static double CheckSpecificity(
    string probeSequence,
    ISuffixTree genomeIndex)
```

**Parameters:**
- `probeSequence`: The probe sequence to check
- `genomeIndex`: Pre-built suffix tree index of the genome

**Returns:** Specificity score (0.0 to 1.0)

## Invariants

1. **Specificity Range**: 0.0 ≤ specificityScore ≤ 1.0 (Source: Implementation)
2. **Self-Complementarity Range**: 0.0 ≤ selfComplementarity ≤ 1.0 (Source: Mathematical definition)
3. **Off-Target Count Non-Negative**: offTargetHits ≥ 0 (Source: Implementation)
4. **Unique Probe Maximum Specificity**: offTargetHits == 1 → specificityScore == 1.0 (Source: Implementation)
5. **No Match Zero Specificity**: offTargetHits == 0 → specificityScore == 0.0 (Source: Implementation)

## Related Methods

- `DesignProbes`: Designs probes with optional genome-wide specificity check
- `FindApproximateMatches`: Internal method for approximate string matching

## References

1. Wikipedia: Hybridization probe - https://en.wikipedia.org/wiki/Hybridization_probe
2. Wikipedia: DNA microarray - https://en.wikipedia.org/wiki/DNA_microarray
3. Wikipedia: Off-target genome editing - https://en.wikipedia.org/wiki/Off-target_genome_editing
4. Wikipedia: BLAST (biotechnology) - https://en.wikipedia.org/wiki/BLAST_(biotechnology)
5. Altschul et al. (1990) - Basic local alignment search tool, J. Mol. Biol.
6. Amann R, Ludwig W (2000) - Ribosomal RNA-targeted nucleic acid probes for studies in microbial ecology, FEMS Microbiology Reviews
