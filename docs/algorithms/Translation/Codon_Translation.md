# Codon Translation

**Algorithm Group:** Translation  
**Test Unit ID:** TRANS-CODON-001  
**Complexity:** O(1) per codon lookup

---

## Overview

Codon translation is the fundamental process of mapping three-nucleotide sequences (codons) to amino acids according to the genetic code. This is a core operation in molecular biology that underlies protein synthesis.

---

## Theoretical Background

### The Genetic Code (Source: Wikipedia, NCBI)

The genetic code is the set of rules by which information encoded in genetic material (DNA or RNA) is translated into proteins. Key properties:

- **Triplet code**: Each codon consists of exactly 3 nucleotides
- **Non-overlapping**: Codons are read sequentially without overlap
- **Degenerate**: Multiple codons can encode the same amino acid (64 codons → 20 amino acids + stop)
- **Nearly universal**: Most organisms use the same standard code with minor variations

### Codon Composition

- **64 possible codons**: 4³ combinations of A, U/T, G, C
- **61 sense codons**: Encode 20 amino acids
- **3 stop codons**: UAA (ochre), UAG (amber), UGA (opal)

### Alternative Genetic Codes (Source: NCBI Translation Tables)

| Table | Name | Key Differences |
|-------|------|-----------------|
| 1 | Standard | Universal default |
| 2 | Vertebrate Mitochondrial | AGA/AGG=Stop, AUA=Met, UGA=Trp |
| 3 | Yeast Mitochondrial | CUU/CUC/CUA/CUG=Thr, AUA=Met, UGA=Trp |
| 11 | Bacterial/Plastid | Same as standard, different start codons |

---

## Implementation Notes

### Current Implementation (Seqeron.Genomics.Core.GeneticCode)

The implementation provides:

1. **Translate(codon)**: Maps a 3-character codon to single-letter amino acid
2. **IsStartCodon(codon)**: Checks if codon can initiate translation
3. **IsStopCodon(codon)**: Checks if codon terminates translation
4. **GetCodonsForAminoAcid(aa)**: Reverse lookup for degeneracy
5. **GetByTableNumber(n)**: Factory for supported genetic codes

### Input Normalization

- DNA codons automatically converted to RNA (T→U)
- Case-insensitive matching (AUG = aug = Aug)

### Supported Genetic Codes

| Property | Table Number | Name |
|----------|--------------|------|
| `GeneticCode.Standard` | 1 | Standard |
| `GeneticCode.VertebrateMitochondrial` | 2 | Vertebrate Mitochondrial |
| `GeneticCode.YeastMitochondrial` | 3 | Yeast Mitochondrial |
| `GeneticCode.BacterialPlastid` | 11 | Bacterial, Archaeal and Plant Plastid |

---

## Invariants

1. **Bijective mapping**: Each codon maps to exactly one amino acid
2. **Complete coverage**: All 64 codons must be defined
3. **Stop codon representation**: Stop codons return '*' character
4. **DNA/RNA equivalence**: ATG and AUG produce identical results
5. **Case insensitivity**: Case variations produce identical results

---

## Edge Cases

| Case | Input | Expected |
|------|-------|----------|
| Standard start | AUG | M |
| Stop codons | UAA, UAG, UGA | * |
| DNA input | ATG | M |
| Lowercase | aug | M |
| Mixed case | AuG | M |
| Invalid length | AU | ArgumentException |
| Invalid codon | XYZ | ArgumentException |
| Null input | null | ArgumentException |

---

## Codon Degeneracy

Amino acids have variable numbers of synonymous codons:

| # Codons | Amino Acids |
|----------|-------------|
| 1 | Met (M), Trp (W) |
| 2 | Phe (F), Tyr (Y), His (H), Gln (Q), Asn (N), Lys (K), Asp (D), Glu (E), Cys (C) |
| 3 | Ile (I) |
| 4 | Val (V), Pro (P), Thr (T), Ala (A), Gly (G) |
| 6 | Leu (L), Ser (S), Arg (R) |

---

## References

1. Wikipedia. "Genetic code." https://en.wikipedia.org/wiki/Genetic_code
2. Wikipedia. "Start codon." https://en.wikipedia.org/wiki/Start_codon
3. Wikipedia. "Stop codon." https://en.wikipedia.org/wiki/Stop_codon
4. NCBI. "The Genetic Codes." https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Crick FH (1968). The origin of the genetic code. J Mol Biol. 38(3):367-79.
