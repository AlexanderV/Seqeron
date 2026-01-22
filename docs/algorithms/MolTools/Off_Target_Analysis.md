# Off-Target Analysis

## Overview

Off-target analysis predicts potential unintended cleavage sites for CRISPR guide RNAs (gRNAs). When a guide RNA targets a genomic sequence, similar sequences elsewhere in the genome may also be cleaved, leading to off-target mutations. This analysis identifies these potential sites and scores their likelihood of causing off-target activity.

## Algorithm

### Core Concept

The algorithm searches for genomic sites that:
1. Have a valid PAM sequence for the CRISPR system
2. Have sequence similarity to the guide RNA (within allowed mismatches)
3. Are not exact matches (which would be on-targets)

### Method: `FindOffTargets`

```
Input:
  - guideSequence: The 20bp (typically) guide RNA sequence
  - genome: The genomic sequence to search
  - maxMismatches: Maximum allowed mismatches (0-5)
  - systemType: CRISPR system (SpCas9, SaCas9, Cas12a, etc.)

Output:
  - Collection of OffTargetSite records with position, mismatches, and score

Algorithm:
1. Find all PAM sites in the genome (both strands)
2. For each PAM site:
   a. Extract the target sequence (guide length upstream/downstream of PAM)
   b. Count mismatches between guide and target
   c. If 0 < mismatches ≤ maxMismatches:
      - Calculate off-target score based on mismatch positions
      - Record mismatch positions
      - Yield the off-target site
```

### Position-Dependent Scoring

The off-target score reflects the likelihood of cleavage at the off-target site. Mismatches in different regions have different effects:

**Seed Region** (PAM-proximal, typically last 10-12bp for Cas9):
- Mismatches in the seed region reduce off-target activity but are still concerning
- Implementation weights seed mismatches higher (5 points per mismatch)

**PAM-Distal Region** (first 8-10bp):
- Mismatches here are more tolerated by the Cas9 protein
- Implementation weights these lower (2 points per mismatch)

### Method: `CalculateSpecificityScore`

```
Input:
  - guideSequence: The guide RNA sequence
  - genome: The genomic sequence to analyze
  - systemType: CRISPR system type

Output:
  - Specificity score (0-100, higher = more specific)

Algorithm:
1. Find all off-targets with up to 4 mismatches
2. If no off-targets: return 100 (maximum specificity)
3. Sum the off-target scores for all found sites
4. Return max(0, 100 - totalPenalty)
```

## Complexity

- **Time**: O(n × m) where n = genome length, m = guide length
- **Space**: O(k) where k = number of off-targets found

In practice, complexity may be higher due to PAM site enumeration on both strands.

## Evidence Base

### Key Scientific Findings

| Finding | Source |
|---------|--------|
| Off-target mutations occur with 3-5 base mismatches | Wikipedia: Off-target genome editing |
| Seed sequence (10-12nt at PAM-proximal) is critical for specificity | Hsu et al. (2013), Fu et al. (2013) |
| PAM-proximal mismatches (8-14bp) define single-base specificity | Hsu et al. (2013) Nature Biotechnology |
| 2+ mismatches in PAM-proximal region considerably reduce activity | Hsu et al. (2013) |
| 3+ interspaced or 5+ concatenated mismatches eliminate detectable cleavage | Hsu et al. (2013) |
| NAG PAM can also cause off-target activity at ~20% efficiency | Hsu et al. (2013) |

### Validation Studies

The scoring model is based on empirical data from:
- >700 guide RNA variants tested (Hsu et al. 2013)
- >100 predicted off-target loci assessed
- Deep sequencing quantification of indel frequencies

## Implementation Notes

### CRISPR System Variations

| System | PAM | Guide Length | PAM Position | Seed Region |
|--------|-----|--------------|--------------|-------------|
| SpCas9 | NGG | 20bp | After target | Last 12bp |
| SaCas9 | NNGRRT | 21bp | After target | Last 12bp |
| Cas12a | TTTV | 23bp | Before target | First 12bp |

### Seed Region Calculation

For Cas9 (PAM after target):
- Seed is positions guide.Length - 12 to guide.Length - 1

For Cas12a (PAM before target):
- Seed is positions 0 to 11

## Limitations

1. **Simplified Scoring**: The implementation uses a basic position-weighted model. More sophisticated models (CFD score, MIT specificity score) incorporate base-pair-specific weights.

2. **PAM Flexibility**: Some off-targets occur at non-canonical PAMs (NAG) with reduced efficiency. The implementation uses strict PAM matching for the specified system.

3. **Chromatin Context**: In vivo off-target activity is affected by chromatin accessibility, which cannot be predicted from sequence alone.

4. **Bulge Mismatches**: The current implementation only considers base mismatches, not insertions/deletions (bulges) between guide and target.

## Usage Example

```csharp
// Find off-targets for a guide
string guide = "ACGTACGTACGTACGTACGT";
var genome = new DnaSequence("ACGT...long genome...ACGT");

var offTargets = CrisprDesigner.FindOffTargets(
    guide, 
    genome, 
    maxMismatches: 3,
    CrisprSystemType.SpCas9);

foreach (var ot in offTargets)
{
    Console.WriteLine($"Position: {ot.Position}, Mismatches: {ot.Mismatches}");
    Console.WriteLine($"Mismatch positions: {string.Join(",", ot.MismatchPositions)}");
    Console.WriteLine($"Off-target score: {ot.OffTargetScore}");
}

// Calculate overall specificity
double specificity = CrisprDesigner.CalculateSpecificityScore(
    guide, genome, CrisprSystemType.SpCas9);
Console.WriteLine($"Specificity: {specificity}/100");
```

## References

1. Hsu PD, Scott DA, Weinstein JA, et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. *Nature Biotechnology*, 31(9):827-832. doi:10.1038/nbt.2647

2. Fu Y, Foden JA, Khayter C, et al. (2013). High-frequency off-target mutagenesis induced by CRISPR-Cas nucleases in human cells. *Nature Biotechnology*, 31(9):822-826. doi:10.1038/nbt.2623

3. Wikipedia. Off-target genome editing. https://en.wikipedia.org/wiki/Off-target_genome_editing
