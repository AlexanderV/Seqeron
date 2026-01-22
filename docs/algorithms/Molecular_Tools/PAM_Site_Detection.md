# PAM Site Detection Algorithm

## Overview

PAM (Protospacer Adjacent Motif) site detection is the process of identifying DNA sequences adjacent to potential CRISPR target sites that are required for Cas nuclease binding and cleavage.

## Biological Background

### Definition

A protospacer adjacent motif (PAM) is a 2–6 base pair DNA sequence immediately following the DNA sequence targeted by the Cas9 nuclease in the CRISPR bacterial adaptive immune system [Wikipedia: PAM].

PAM is an essential targeting component which distinguishes bacterial self from non-self DNA, preventing the CRISPR locus from being targeted and destroyed by the CRISPR-associated nuclease [Wikipedia: PAM].

### PAM Requirements

Cas9 will not successfully bind to or cleave the target DNA sequence if it is not followed by the PAM sequence [Wikipedia: PAM, Jinek et al. 2012].

## PAM Sequences by CRISPR System

| System | PAM Sequence | PAM Location | Guide Length | Notes |
|--------|--------------|--------------|--------------|-------|
| SpCas9 | NGG | 3' of target | 20 bp | Canonical Cas9 from *Streptococcus pyogenes* |
| SpCas9-NAG | NAG | 3' of target | 20 bp | Lower efficiency variant |
| SaCas9 | NNGRRT | 3' of target | 21 bp | *Staphylococcus aureus* Cas9 |
| Cas12a (Cpf1) | TTTV | 5' of target | 23 bp | V = A, C, or G |
| AsCas12a | TTTV | 5' of target | 23 bp | *Acidaminococcus sp.* |
| LbCas12a | TTTV | 5' of target | 24 bp | *Lachnospiraceae bacterium* |
| CasX | TTCN | 5' of target | 20 bp | Compact Cas protein |

**Sources:** Wikipedia (CRISPR, PAM), Zetsche et al. 2015 (Cas12a), Jinek et al. 2012 (SpCas9)

### IUPAC Codes in PAM Patterns

| Code | Meaning | Nucleotides |
|------|---------|-------------|
| N | Any | A, C, G, T |
| R | Purine | A, G |
| V | Not T | A, C, G |

**Source:** IUPAC-IUB nomenclature (1970)

## Algorithm

### PAM Detection Process

1. **Input validation**: Verify sequence is valid DNA
2. **Forward strand search**:
   - Scan sequence for PAM pattern matches using IUPAC matching
   - For each PAM match at position `i`:
     - Calculate target region based on PAM position (before or after)
     - Verify target region is within sequence bounds
     - Extract target sequence of appropriate guide length
3. **Reverse strand search**:
   - Generate reverse complement of input sequence
   - Repeat scanning process on reverse complement
   - Convert positions back to forward strand coordinates
4. **Return all found PAM sites with metadata**

### Position Conventions

- **PAM after target (Cas9 systems)**: `[Target 20bp][NGG]`
  - Target region: positions `PAM_pos - guideLength` to `PAM_pos - 1`
  
- **PAM before target (Cas12a systems)**: `[TTTV][Target 23bp]`
  - Target region: positions `PAM_pos + PAM_length` to `PAM_pos + PAM_length + guideLength - 1`

### Complexity

- **Time**: O(n) where n = sequence length
- **Space**: O(k) where k = number of PAM sites found

## Implementation Notes

### Current Implementation

The `CrisprDesigner.FindPamSites()` method in this library:

1. Searches both forward and reverse strands
2. Uses IUPAC pattern matching via `IupacHelper.MatchesIupac()`
3. Returns `PamSite` records containing:
   - Position (on forward strand)
   - PAM sequence (actual nucleotides matched)
   - Target sequence
   - Target start position
   - Strand orientation
   - CRISPR system information

### GetSystem Method

Returns `CrisprSystem` record with:
- Name, PAM sequence pattern, guide length
- PAM position relative to target (before/after)
- Description of the system

## Edge Cases

| Case | Expected Behavior |
|------|-------------------|
| Empty sequence | Return empty collection |
| Null sequence | Throw ArgumentNullException |
| Sequence shorter than PAM + guide | Return empty (no valid sites) |
| No PAM matches | Return empty collection |
| PAM at sequence boundary | Include if target fits within bounds |
| Multiple overlapping PAMs | Return all distinct PAM sites |
| Both strands have matches | Return sites from both strands |
| Lowercase input | Case-insensitive matching |

## References

1. Wikipedia: Protospacer adjacent motif. https://en.wikipedia.org/wiki/Protospacer_adjacent_motif
2. Wikipedia: CRISPR. https://en.wikipedia.org/wiki/CRISPR
3. Jinek M, et al. (2012). "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity". Science 337(6096):816-821.
4. Zetsche B, et al. (2015). "Cpf1 is a single RNA-guided endonuclease of a class 2 CRISPR-Cas system". Cell 163(3):759-771.
5. Anders C, et al. (2014). "Structural basis of PAM-dependent target DNA recognition by the Cas9 endonuclease". Nature 513(7519):569-573.
