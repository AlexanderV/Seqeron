# Evidence: TRANS-CODON-001 - Codon Translation

**Test Unit ID:** TRANS-CODON-001  
**Area:** Translation  
**Date:** 2026-02-04

---

## Authoritative Sources

### Primary Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia: Genetic code](https://en.wikipedia.org/wiki/Genetic_code) | Encyclopedia | Standard codon table, degeneracy, 64 codons → 20 amino acids + 3 stop codons |
| [Wikipedia: Start codon](https://en.wikipedia.org/wiki/Start_codon) | Encyclopedia | AUG as universal start codon, alternative starts (GUG, UUG) in bacteria |
| [Wikipedia: Stop codon](https://en.wikipedia.org/wiki/Stop_codon) | Encyclopedia | Three stop codons (UAA/ochre, UAG/amber, UGA/opal), nomenclature |
| [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) | Official Reference | Complete translation tables 1-33, organism-specific variations |

### Reference Datasets from NCBI Translation Tables

#### Table 1: Standard Genetic Code
```
AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Stop codons: UAA, UAG, UGA
Start codons: AUG (primary), alternative: GUG, UUG (rare in eukaryotes)
```

#### Table 2: Vertebrate Mitochondrial Code
```
AAs  = FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
Starts = ----------**--------------------MMMM----------**---M------------
Differences from Standard:
- AGA, AGG → Stop (not Arg)
- AUA → Met (not Ile)
- UGA → Trp (not Stop)
Start codons: AUG, AUA, AUU, AUC
```

#### Table 3: Yeast Mitochondrial Code
```
AAs  = FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ----------**----------------------MM---------------M------------
Differences from Standard:
- CUU, CUC, CUA, CUG → Thr (not Leu)
- AUA → Met (not Ile)
- UGA → Trp (not Stop)
Start codons: AUG, AUA
```

#### Table 11: Bacterial, Archaeal and Plant Plastid Code
```
AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M------------MMMM---------------M------------
Same codon table as Standard, but alternative start codons are common:
Start codons: AUG, GUG, UUG
```

---

## Documented Corner Cases

### From Wikipedia/NCBI

| Corner Case | Expected Behavior | Source |
|-------------|-------------------|--------|
| Codon length ≠ 3 | Error/Exception | Definition of codon (NCBI) |
| Unknown codon (e.g., NNN) | Error/Exception | Standard genetic code definition |
| DNA vs RNA input | Both should work (T↔U conversion) | Implementation decision |
| Case sensitivity | Case-insensitive (AUG = aug = AuG) | Common convention |
| Stop codon → '*' | Standard representation | NCBI format |
| Start codon as non-initiator | Translates to M (not fMet) | Wikipedia: Start codon |

### Degeneracy (Codon Redundancy)

| Amino Acid | # Codons | Example Codons |
|------------|----------|----------------|
| Met (M) | 1 | AUG |
| Trp (W) | 1 | UGG |
| Phe (F) | 2 | UUU, UUC |
| Leu (L) | 6 | UUA, UUG, CUU, CUC, CUA, CUG |
| Ser (S) | 6 | UCU, UCC, UCA, UCG, AGU, AGC |
| Arg (R) | 6 | CGU, CGC, CGA, CGG, AGA, AGG |
| Stop (*) | 3 | UAA, UAG, UGA |

---

## Test Datasets

### Complete Standard Codon Table (64 codons)

| Codon | AA | Codon | AA | Codon | AA | Codon | AA |
|-------|----|----|----|----|----|----|----| 
| UUU | F | UCU | S | UAU | Y | UGU | C |
| UUC | F | UCC | S | UAC | Y | UGC | C |
| UUA | L | UCA | S | UAA | * | UGA | * |
| UUG | L | UCG | S | UAG | * | UGG | W |
| CUU | L | CCU | P | CAU | H | CGU | R |
| CUC | L | CCC | P | CAC | H | CGC | R |
| CUA | L | CCA | P | CAA | Q | CGA | R |
| CUG | L | CCG | P | CAG | Q | CGG | R |
| AUU | I | ACU | T | AAU | N | AGU | S |
| AUC | I | ACC | T | AAC | N | AGC | S |
| AUA | I | ACA | T | AAA | K | AGA | R |
| AUG | M | ACG | T | AAG | K | AGG | R |
| GUU | V | GCU | A | GAU | D | GGU | G |
| GUC | V | GCC | A | GAC | D | GGC | G |
| GUA | V | GCA | A | GAA | E | GGA | G |
| GUG | V | GCG | A | GAG | E | GGG | G |

### Verification Codons for Alternative Genetic Codes

#### Vertebrate Mitochondrial (Table 2) Differences
| Codon | Standard | Vert Mito |
|-------|----------|-----------|
| AGA | R | * |
| AGG | R | * |
| AUA | I | M |
| UGA | * | W |

#### Yeast Mitochondrial (Table 3) Differences
| Codon | Standard | Yeast Mito |
|-------|----------|------------|
| CUU | L | T |
| CUC | L | T |
| CUA | L | T |
| CUG | L | T |
| AUA | I | M |
| UGA | * | W |

---

## Known Failure Modes

| Failure Mode | Cause | Expected Behavior |
|--------------|-------|-------------------|
| Empty codon | Invalid input | ArgumentException |
| Codon too short (2 chars) | Invalid format | ArgumentException |
| Codon too long (4+ chars) | Invalid format | ArgumentException |
| Null codon | Invalid input | ArgumentNullException |
| Invalid nucleotide (X, Z) | Unknown codon | ArgumentException |
| Invalid table number | Unsupported code | ArgumentException |

---

## Testing Methodologies

### Recommended Test Categories (NCBI/Wikipedia-based)

1. **Complete coverage of standard table** - All 64 codons
2. **Start codon identification** - AUG + alternatives by table
3. **Stop codon identification** - UAA, UAG, UGA + table-specific
4. **Alternative genetic code verification** - Key differences from standard
5. **DNA/RNA normalization** - T↔U equivalence
6. **Case insensitivity** - Mixed case handling
7. **Reverse lookup (amino acid → codons)** - Degeneracy verification
8. **Input validation** - Invalid codons, lengths, nulls

---

## Implementation Notes

### Current Implementation (GeneticCode.cs)
- Supports 4 genetic codes: Standard (1), Vertebrate Mitochondrial (2), Yeast Mitochondrial (3), Bacterial/Plastid (11)
- Normalizes DNA to RNA (T→U) internally
- Returns '*' for stop codons
- Throws ArgumentException for invalid codons

---

## References

1. Wikipedia contributors. "Genetic code." Wikipedia, The Free Encyclopedia.
2. Wikipedia contributors. "Start codon." Wikipedia, The Free Encyclopedia.
3. Wikipedia contributors. "Stop codon." Wikipedia, The Free Encyclopedia.
4. Elzanowski A, Ostell J. "The Genetic Codes." NCBI Taxonomy.
5. Nirenberg M, Matthaei JH (1961). The dependence of cell-free protein synthesis in E. coli upon naturally occurring or synthetic polyribonucleotides.
6. Crick FH (1968). The origin of the genetic code. J Mol Biol.
