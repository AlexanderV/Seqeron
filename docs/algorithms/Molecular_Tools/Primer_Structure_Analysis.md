# Primer Structure Analysis

## Overview

Primer structure analysis evaluates PCR primers for secondary structure formation and self-complementarity issues that can reduce amplification efficiency.

## Algorithm Category

- **Test Unit ID:** PRIMER-STRUCT-001
- **AlgorithmGroup:** Molecular_Tools (MolTools)
- **Area:** Primer Design / Secondary Structure Detection

## Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| `HasHairpinPotential(seq, minStemLength, minLoopLength)` | Detects self-complementary regions that can form hairpin structures | O(n²) for <100bp, O(n) for ≥100bp* |
| `HasPrimerDimer(primer1, primer2, minComp)` | Detects 3' end complementarity between primers | O(n) |
| `Calculate3PrimeStability(seq)` | Calculates ΔG of the 3' pentamer using nearest-neighbor thermodynamics | O(1) |
| `FindLongestHomopolymer(seq)` | Finds longest mononucleotide repeat (e.g., AAAA) | O(n) |
| `FindLongestDinucleotideRepeat(seq)` | Finds longest dinucleotide repeat (e.g., ATATAT) | O(n) |

*Uses suffix tree optimization for long sequences (≥100bp)

## Theory

### Hairpin Structures (Stem-Loop)

Hairpins form when a primer has self-complementary regions separated by a loop. Key properties:
- **Minimum stem length:** Typically 4 bp (stable enough to interfere with PCR)
- **Minimum loop length:** 3 nucleotides (sterically required)
- **Stability:** GC-rich stems are more stable than AT-rich

*Source: Wikipedia (Stem-loop), Primer3 Manual (PRIMER_MAX_HAIRPIN_TH)*

### Primer-Dimer Formation

Primer-dimers form when primers have complementary 3' ends, allowing self-extension:
- 3' end complementarity is critical (extension occurs here)
- High GC content at 3' increases dimer stability
- Forms 30-50 bp products visible on gel electrophoresis

*Source: Wikipedia (Primer dimer), Primer3 Manual (PRIMER_MAX_SELF_END)*

### 3' End Stability (ΔG Calculation)

The stability of the 3' terminal bases determines extension efficiency:
- More negative ΔG = more stable binding
- Calculated using nearest-neighbor thermodynamic parameters
- Last 5 bases are typically evaluated (Primer3 standard)
- SantaLucia (1998) nearest-neighbor parameters:
  - Most stable 5mer: GCGCG (ΔG = -6.86 kcal/mol)
  - Least stable 5mer: TATAT (ΔG = -0.86 kcal/mol)

*Source: SantaLucia (1998) PNAS 95:1460-65, Primer3 Manual (PRIMER_MAX_END_STABILITY)*

### Homopolymer Runs

Consecutive identical nucleotides that cause:
- Polymerase slippage
- Mispriming
- Primer3 default max: 5 bp

*Source: Primer3 Manual (PRIMER_MAX_POLY_X)*

### Dinucleotide Repeats

Alternating dinucleotide patterns (e.g., ATATAT, GCGCGC):
- Contribute to mishybridization
- Loop formation potential
- Should be avoided in primer design

*Source: Wikipedia (Primer molecular biology - PCR primer design section)*

## Implementation Notes

### Hairpin Detection Algorithm

The implementation uses two strategies based on sequence length:

#### Simple Algorithm (for sequences <100bp)
O(n²) nested loop approach optimized for short primers:
1. Extract fragment of minStemLength at position i
2. Search for complementary sequence at position j ≥ i + minStemLength + minLoopLength
3. Compare fragment against reverse of target for complementarity

#### Suffix Tree Algorithm (for sequences ≥100bp)
O(n) approach using the SuffixTree library:
1. Build suffix tree on sequence: O(n)
2. Compute reverse complement of sequence
3. For each position p in revComp, search for stems via `FindAllOccurrences`: O(m + k)
4. Check loop constraint: j ≥ i + stemLength + loopLength
5. Position mapping: position p in revComp corresponds to position (n - p - stemLength) in original

**Break-even point:** ~100bp based on suffix tree construction overhead vs O(n²) iterations.
For typical PCR primers (18-25bp), the simple algorithm is faster.

### Primer-Dimer Detection Algorithm

Focuses on 3' end complementarity:
1. Extract last 8 bases of primer1
2. Compute reverse complement of primer2
3. Compare 5' of revcomp(primer2) with 3' of primer1
4. Count complementary pairs

### 3' Stability Calculation

Uses simplified nearest-neighbor ΔG values:
- Extracts last 5 bases
- Sums ΔG for each dinucleotide pair
- Values based on SantaLucia (1998) unified model

## Nearest-Neighbor ΔG Values Used

| Dinucleotide | ΔG (kcal/mol) |
|--------------|---------------|
| AA/TT | -1.00 |
| AT | -0.88 |
| TA | -0.58 |
| CA/TG | -1.45 |
| GT/AC | -1.44 |
| CT/AG | -1.28 |
| GA/TC | -1.30 |
| CG | -2.17 |
| GC | -2.24 |
| GG/CC | -1.84 |

*Values from SantaLucia (1998)*

## Edge Cases

### Hairpin Detection
- Sequence too short: Returns false (minimum = 2×minStemLength + 3)
- No self-complementary regions: Returns false
- Perfect palindrome: May or may not form hairpin depending on loop

### Primer-Dimer
- Empty primer: Returns false
- Self-dimer: Same primer compared with itself
- Non-complementary 3' ends: Returns false

### Homopolymer
- Empty sequence: Returns 0
- All unique bases: Returns 1
- All same base: Returns sequence length

### Dinucleotide Repeat
- Sequence < 4 bp: Returns 0
- No repeats: Returns 1
- Perfect repeat: Returns (length/2)

## References

1. Wikipedia - Primer (molecular biology): PCR primer design section
2. Wikipedia - Primer dimer: Mechanism of formation
3. Wikipedia - Stem-loop (Hairpin loop): Formation and stability
4. Wikipedia - Nucleic acid thermodynamics: Nearest-neighbor method
5. SantaLucia JR (1998) "A unified view of polymer, dumbbell and oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95:1460-65
6. Primer3 Manual (primer3.org): PRIMER_MAX_HAIRPIN_TH, PRIMER_MAX_SELF_END, PRIMER_MAX_END_STABILITY, PRIMER_MAX_POLY_X
