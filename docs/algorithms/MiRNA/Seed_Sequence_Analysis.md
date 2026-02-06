# Seed Sequence Analysis

## Algorithm Overview

MicroRNA (miRNA) seed sequence analysis identifies and compares the critical 5' region of mature miRNAs that determines target specificity. The seed region is the primary determinant of miRNA-mRNA interaction in animals.

## Definitions

### Seed Region
- **Position**: Nucleotides 2–8 (1-indexed) from the 5' end of the mature miRNA (7 nucleotides)
- **Core seed**: Positions 2–7 (6 nucleotides), as defined by Bartel (2009) and TargetScan
- **Extended seed**: Positions 2–8 (7 nucleotides), commonly used in computational tools

### miRNA Family
miRNAs sharing the same seed sequence belong to the same functional family. Family members are expected to target overlapping sets of mRNAs (Griffiths-Jones et al., 2006; TargetScan).

## Implementation

### `GetSeedSequence(string miRnaSequence) → string`
- **Input**: Mature miRNA sequence (RNA, uppercase or mixed case)
- **Output**: 7-nucleotide seed (positions 2-8), uppercased
- **Complexity**: O(1)
- **Edge cases**: Returns empty string for null, empty, or sequences shorter than 8 nt
- **Note**: Implementation uses the extended seed (positions 2-8, 7 nt) consistent with TargetScan's matching approach

### `CreateMiRna(string name, string sequence) → MiRna`
- **Input**: miRNA name and sequence (RNA or DNA)
- **Output**: MiRna record with extracted seed, name, and positional metadata
- **DNA handling**: Converts T→U before seed extraction
- **Case handling**: Normalizes to uppercase

### `CompareSeedRegions(MiRna mirna1, MiRna mirna2) → SeedComparison`
- **Input**: Two MiRna records
- **Output**: SeedComparison record containing:
  - Number of matches between seed sequences
  - Number of mismatches (Hamming distance)
  - Boolean indicating same family membership (identical seeds)
- **Complexity**: O(k) where k = seed length (7)
- **Semantics**: Character-by-character comparison of pre-extracted seed sequences

## Seed Site Type Hierarchy (context for understanding)
From TargetScan and Bartel (2009), canonical site types in order of effectiveness:
1. **8mer**: Seed match (pos 2-8) + A at position 1 → strongest repression
2. **7mer-m8**: Seed match (pos 2-8)
3. **7mer-A1**: Core seed (pos 2-7) + A at position 1
4. **6mer**: Core seed (pos 2-7) only

## Sources
- Bartel DP (2009). "MicroRNAs: Target recognition and regulatory functions." Cell 136(2):215-233. PMID: 19167326
- Lewis BP, Burge CB, Bartel DP (2005). "Conserved seed pairing." Cell 120(1):15-20. PMID: 15652477
- Grimson A et al. (2007). "MicroRNA targeting specificity in mammals." Molecular Cell 27(1):91-105. PMID: 17612493
- Agarwal V et al. (2015). "Predicting effective microRNA target sites." eLife 4:e05005. PMID: 26267216
- miRBase: https://mirbase.org/

## Implementation-Specific Notes
- The library uses the 7-nucleotide extended seed (positions 2-8) as its standard seed representation
- Seed comparison uses Hamming distance (no gaps, fixed-length comparison)
- ASSUMPTION: CompareSeedRegions handles empty seeds gracefully (returns 0 matches, 0 mismatches, not same family)
