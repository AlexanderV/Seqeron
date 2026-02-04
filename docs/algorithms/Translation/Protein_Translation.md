# Protein Translation

## Overview

Protein translation is the process of converting a nucleotide sequence (DNA or RNA) into a protein sequence (amino acid chain) using the genetic code. The `Translator` class provides methods for translating sequences in different reading frames, finding Open Reading Frames (ORFs), and performing six-frame translations.

---

## Theoretical Background

### The Translation Process (Source: Wikipedia)

Translation reads nucleotides in groups of three called **codons**. Each codon maps to either:
- One of 20 standard amino acids
- A stop signal (termination)

**Key properties:**
- Direction: 5' → 3' (amino-to-carboxyl)
- Triplet code: Each codon = 3 nucleotides
- Non-overlapping: Codons read sequentially without overlap
- Degenerate: Multiple codons can encode the same amino acid

### Reading Frames (Source: Wikipedia)

A **reading frame** defines which triplets are grouped as codons:
- **Frame 0**: Start reading from position 0
- **Frame 1**: Start reading from position 1 (skip 1 nucleotide)
- **Frame 2**: Start reading from position 2 (skip 2 nucleotides)

For double-stranded DNA:
- **6 total frames**: 3 forward (+1, +2, +3) and 3 reverse complement (-1, -2, -3)

### Open Reading Frames (Source: Wikipedia)

An **ORF** is a sequence that:
1. Begins with a start codon (typically AUG)
2. Continues through a region with length divisible by 3
3. Ends with a stop codon (UAA, UAG, UGA)

ORFs are used in gene prediction. Short ORFs (<100 codons) may still produce functional peptides.

---

## Algorithm

### Basic Translation

**Input:**
- `sequence`: DNA or RNA string
- `geneticCode`: Translation table (default: Standard)
- `frame`: Reading frame offset (0, 1, or 2)
- `toFirstStop`: Stop at first stop codon (boolean)

**Output:**
- `ProteinSequence`: Amino acid sequence

**Process:**
```
1. Convert DNA to RNA (T → U)
2. Skip 'frame' nucleotides from start
3. For each triplet from position 'frame':
   a. Extract codon (3 nucleotides)
   b. Translate codon to amino acid using genetic code
   c. If toFirstStop and amino acid is '*': stop
   d. Append amino acid to result
4. Return protein sequence
```

**Complexity:**
- Time: O(n/3) where n = sequence length
- Space: O(n/3) for output

### Six-Frame Translation

**Input:**
- `dna`: DNA sequence
- `geneticCode`: Translation table

**Output:**
- Dictionary<int, ProteinSequence> with keys: +1, +2, +3, -1, -2, -3

**Process:**
```
1. Forward strand translation:
   - Frame +1: Translate(dna, frame=0)
   - Frame +2: Translate(dna, frame=1)
   - Frame +3: Translate(dna, frame=2)
2. Reverse complement strand:
   - revComp = ReverseComplement(dna)
   - Frame -1: Translate(revComp, frame=0)
   - Frame -2: Translate(revComp, frame=1)
   - Frame -3: Translate(revComp, frame=2)
```

### ORF Finding

**Input:**
- `dna`: DNA sequence
- `geneticCode`: Translation table
- `minLength`: Minimum ORF length in amino acids
- `searchBothStrands`: Include reverse complement

**Output:**
- Enumerable of `OrfResult` (StartPosition, EndPosition, Frame, Protein)

**Process:**
```
1. For each strand (forward, optionally reverse complement):
   a. For each frame (0, 1, 2):
      i. Scan for start codons
      ii. Accumulate amino acids until stop codon
      iii. If length >= minLength: yield ORF result
      iv. Continue scanning for next start codon
```

**Complexity:**
- Time: O(n) per strand/frame = O(6n) for full search
- Space: O(k) where k = number of ORFs found

---

## Implementation Notes

### Current Implementation (Seqeron.Genomics.Core.Translator)

**Methods:**
| Method | Description |
|--------|-------------|
| `Translate(DnaSequence, ...)` | Translate DNA to protein |
| `Translate(RnaSequence, ...)` | Translate RNA to protein |
| `Translate(string, ...)` | Translate sequence string |
| `TranslateSixFrames(DnaSequence, ...)` | All 6 reading frames |
| `FindOrfs(DnaSequence, ...)` | Find all ORFs |

**Input Handling:**
- DNA input: T automatically converted to U
- Case-insensitive: Lowercase converted to uppercase
- Empty/null: Empty sequence returns empty protein; null throws ArgumentNullException
- Invalid frame: Throws ArgumentOutOfRangeException for frame > 2 or < 0

**Alternative Genetic Codes:**
- Standard (Table 1) - default
- VertebrateMitochondrial (Table 2) - AGA/AGG = Stop
- YeastMitochondrial (Table 3) - CUx = Threonine
- Bacterial (Table 11) - Alternative starts

---

## Test Considerations

### Invariants
1. `Translate(seq, frame=0)` equals `TranslateSixFrames(seq)[+1]`
2. Translated protein length ≤ floor((seqLength - frame) / 3)
3. ORF protein always starts with start codon amino acid (M in standard code)
4. ORF position: endPosition - startPosition + 1 = nucleotideLength

### Edge Cases
- Empty sequence → empty protein
- Sequence shorter than 3 nucleotides → empty protein
- No start codon → no ORFs found
- All stop codons → empty protein with toFirstStop
- minLength > possible ORF length → empty results

---

## References

1. Wikipedia: "Translation (biology)"
2. Wikipedia: "Reading frame"  
3. Wikipedia: "Open reading frame"
4. NCBI: "The Genetic Codes"
5. Lodish H et al. (2007). Molecular Cell Biology, 6th ed.
