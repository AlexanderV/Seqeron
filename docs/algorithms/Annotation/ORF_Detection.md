# ORF Detection

## Algorithm Overview

Open Reading Frame (ORF) detection identifies potential protein-coding regions in DNA sequences. An ORF is defined as a sequence that starts with a start codon and ends with a stop codon, without any intervening stop codons in the same reading frame.

## Biological Background

### Definition
- **Open Reading Frame (ORF)**: A continuous stretch of codons beginning with a start codon (usually ATG/AUG) and ending with a stop codon (TAA, TAG, TGA), with no stop codons in between [Wikipedia, Rosalind].
- In prokaryotes, only one of six possible reading frames is typically "open" (actively translated) [Wikipedia].
- The presence of an ORF does not guarantee the region is translated [Deonier et al., 2005].

### Six-Frame Translation
Since DNA is read in groups of three nucleotides (codons), there are three reading frames on each strand. The double helix has two antiparallel strands, yielding six total reading frames [Wikipedia, Rosalind]:
- **Forward strand**: Frames +1, +2, +3
- **Reverse complement**: Frames -1, -2, -3

### Start Codons
- **Primary**: ATG (Methionine) - most common start codon
- **Alternative start codons** (prokaryotes): GTG, TTG (less common)
- NCBI ORF Finder supports both 'ATG only' and 'ATG and alternative initiation codons' modes [NCBI ORF Finder].

### Stop Codons (Standard Genetic Code)
- TAA (ochre)
- TAG (amber)
- TGA (opal)

### Minimum Length Considerations
- Some authors require ≥100 codons for a meaningful ORF [Claverie, 1997].
- Others use ≥150 codons [Deonier et al., 2005].
- Short ORFs (sORFs) < 100 codons can still produce functional micropeptides [Wikipedia].

## Algorithm Description

### Basic ORF Detection Algorithm
```
For each of 6 reading frames:
    For each position in frame (step 3):
        If start codon found:
            Track ORF start position
        If stop codon found:
            If valid ORF start exists AND length ≥ minLength:
                Report ORF (start, end, frame, sequence, protein)
            Clear ORF tracking
```

### Key Parameters
| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| minLength | Minimum ORF length (amino acids) | 30-150 aa |
| searchBothStrands | Search reverse complement | true |
| requireStartCodon | Require ATG/GTG/TTG start | true/false |
| geneticCode | Translation table | Standard (1), Bacterial (11) |

## Invariants

1. **Start Codon**: Every ORF must begin with a valid start codon (ATG, GTG, TTG) when requireStartCodon=true
2. **Stop Codon**: Every ORF must end with a stop codon (TAA, TAG, TGA)
3. **Frame Integrity**: ORF nucleotide length must be divisible by 3
4. **No Internal Stops**: No stop codons between start and end in the same frame
5. **Length Constraint**: ORF protein length ≥ minLength
6. **Coordinate Validity**: 0 ≤ start < end ≤ sequence.Length

## Implementation Notes

### GenomeAnnotator.FindOrfs (Canonical)
- Searches all 6 reading frames when `searchBothStrands=true`
- Supports alternative start codons (GTG, TTG)
- Returns `OpenReadingFrame` records with:
  - Start/End positions
  - Frame number (±1, ±2, ±3)
  - IsReverseComplement flag
  - DNA sequence
  - Translated protein sequence

### GenomeAnnotator.FindLongestOrfsPerFrame
- Returns dictionary with frame number → longest ORF
- Keys: 1, 2, 3 (forward) and -1, -2, -3 (reverse complement)

### GenomicAnalyzer.FindOpenReadingFrames (Alternate)
- Wrapper with simpler signature
- Searches only for ATG start codons
- Returns `OrfInfo` records

## Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| Empty sequence | Return empty collection |
| No start codon | Return empty (when requireStartCodon=true) |
| No stop codon | Return empty OR truncated ORF at end (implementation-specific) |
| Very short sequence (< 3bp) | Return empty |
| Lowercase input | Should be handled (case-insensitive) |
| N characters | Skip or handle as unknown |
| Overlapping ORFs | Report all valid ORFs independently |
| Nested ORFs | Report outer and inner if both meet criteria |

## Testing Methodology

### Test Categories
1. **Canonical ORF detection**: ATG...TAA/TAG/TGA patterns
2. **Alternative start codons**: GTG, TTG recognition
3. **Six-frame search**: All reading frames including reverse complement
4. **Minimum length filtering**: Enforce minLength threshold
5. **Edge cases**: Empty, short, no-start, no-stop sequences
6. **Invariant verification**: Coordinate and sequence integrity

### Reference Test Cases

#### Rosalind ORF Problem Dataset
```
>Rosalind_99
AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG

Expected proteins:
MLLGSFRLIPKETLIQVAGSSPCNLS
M
MGMTPRLGLESLLE
MTPRLGLESLLE
```

## References

1. Wikipedia. "Open reading frame." https://en.wikipedia.org/wiki/Open_reading_frame
2. Rosalind. "Open Reading Frames." https://rosalind.info/problems/orf/
3. NCBI ORF Finder. https://www.ncbi.nlm.nih.gov/orffinder/
4. Deonier R, Tavaré S, Waterman M (2005). Computational Genome Analysis: an introduction. Springer-Verlag. p. 25.
5. Claverie JM (1997). "Computational methods for the identification of genes in vertebrate genomic sequences." Human Molecular Genetics 6(10): 1735-44.
6. Sieber P, Platzer M, Schuster S (2018). "The Definition of Open Reading Frame Revisited." Trends in Genetics 34(3): 167-170.
