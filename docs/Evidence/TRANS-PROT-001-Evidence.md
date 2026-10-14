# Evidence: TRANS-PROT-001 - Protein Translation

**Test Unit ID:** TRANS-PROT-001  
**Area:** Translation  
**Date:** 2026-02-04

---

## Primary Sources

### 1. Wikipedia: Translation (biology)
**URL:** https://en.wikipedia.org/wiki/Translation_(biology)

**Key Information:**
- Translation is the process where proteins are produced using RNA molecules as templates
- Nucleotides are read in triplets (codons), each resulting in addition of one specific amino acid
- The matching from nucleotide triplet to amino acid is called the genetic code
- Codons are read in the 5'→3' direction
- The start codon (typically AUG) encodes methionine in eukaryotes and archaea
- Stop codons (UAA, UAG, UGA) terminate translation
- Termination depends on release factors (eRF1) that recognize all three stop codons

**Relevance:** Defines the fundamental process of translating nucleotide sequences to protein sequences.

---

### 2. Wikipedia: Reading Frame
**URL:** https://en.wikipedia.org/wiki/Reading_frame

**Key Information:**
- A reading frame is a specific choice out of the possible ways to read the sequence of nucleotides
- A single strand has three possible reading frames (starting at position 0, 1, or 2)
- Double-stranded DNA has six possible reading frames (three forward, three reverse complement)
- Reading frames are oriented 5' to 3'
- mRNA is single-stranded and contains three possible reading frames, of which only one is translated

**Relevance:** Defines the concept of reading frames which is essential for understanding how frame parameter affects translation.

---

### 3. Wikipedia: Open Reading Frame
**URL:** https://en.wikipedia.org/wiki/Open_reading_frame

**Key Information:**
- An ORF is a reading frame that has the potential to be transcribed and translated
- Requires a continuous sequence starting with a start codon and ending with a stop codon
- ORFs are used as evidence for gene prediction
- Short ORFs (sORFs) < 100 codons can still produce functional peptides
- Six-frame translation: DNA has six possible frame translations (3 forward + 3 reverse complement)
- A stop codon is expected approximately once every 21 codons in random sequence

**Relevance:** Defines ORF finding algorithm requirements and six-frame translation concept.

---

### 4. NCBI Translation Tables
**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

**Key Information:**
- Multiple translation tables exist for different organisms
- Standard code (Table 1): Universal genetic code; start codons AUG, UUG, CUG; stop codons UAA, UAG, UGA
- Vertebrate Mitochondrial (Table 2): AUA→M, UGA→W, AGA→*, AGG→*; start codons AUG, AUA, AUU, AUC, GUG; stop codons UAA, UAG, AGA, AGG
- Yeast Mitochondrial (Table 3): CUU/CUC/CUA/CUG→T, AUA→M, UGA→W; start codons AUG, AUA, GUG; stop codons UAA, UAG
- Bacterial/Archaeal/Plant Plastid (Table 11): Same amino acid translation as standard; start codons AUG, GUG, UUG, CUG, AUU, AUC, AUA; stop codons UAA, UAG, UGA

**Relevance:** Documents alternative genetic codes that should be supported in translation.

---

### 5. UniProt P01308 — Human Insulin
**URL:** https://www.uniprot.org/uniprot/P01308

**Key Information:**
- Human preproinsulin is 110 amino acids: `MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN`
- B chain (positions 25–54): `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` (30 amino acids)
- Gene: INS, Chromosome 11
- Organism: Homo sapiens

**Relevance:** Provides verified biological protein sequence for translation correctness test.

---

## Extracted Test Requirements

### From Wikipedia (Translation biology)
1. **Codon-based reading**: Sequence must be read in triplets
2. **Start codon recognition**: AUG initiates translation in standard code
3. **Stop codon termination**: UAA, UAG, UGA terminate translation
4. **5'→3' direction**: Translation proceeds from 5' to 3' end

### From Wikipedia (Reading Frame)
1. **Three forward frames**: Frame 0, 1, 2 for single strand
2. **Six frames total**: Forward and reverse complement strands
3. **Frame offset**: Each frame shifts reading by one nucleotide

### From Wikipedia (Open Reading Frame)
1. **ORF definition**: Starts with start codon, ends with stop codon
2. **Minimum length filtering**: ORFs can be filtered by minimum amino acid length
3. **Both strands**: ORFs should be searchable on both strands

---

## Edge Cases from Sources

1. **Empty/null sequences** - Edge case for input validation
2. **Sequence length not divisible by 3** - Incomplete final codon
3. **No start codon found** - ORF finder returns empty
4. **Stop codon before minimum length** - ORF filtered out
5. **Multiple ORFs in same sequence** - Multiple results returned
6. **Alternative genetic codes** - Different amino acid mappings
7. **Case insensitivity** - Input should handle mixed case
8. **DNA vs RNA input** - T→U conversion handled automatically
9. **Selenocysteine (U) / pyrrolysine (O)** (review 2026-09, B02 F32) - Biopython 1.88 `IUPACData.extended_protein_letters` = "ACDEFGHIKLMNPQRSTVWYBXZJUO" includes U and O (three-letter Sec/Pyl, `protein_letters_1to3_extended`); `ProteinSequence` accepts them (`ValidCharacters` = these 26 + `*` = 27). Reference values: Biopython `molecular_weight(s, "protein")` MKUAO 735.7975, MSKU 514.4988, MSKO 601.759, ACDUEFGOHK 1293.3102 (EMBOSS 6.6 pepstats 735.80 / 514.50 / 601.76 / 1293.32); EMBOSS 6.6 `iep` pI MSKU = MSKO = MSK = 9.1500, ACDUEFGOHK = ACDEFGHK = 5.4305 (U/O neutral); Biopython `IsoelectricPoint("MKUAO")` = `IsoelectricPoint("MKA")` (pi 8.500108…, charge_at_pH(7) 0.49935568654234397); Biopython `ProteinAnalysis("MKUAO").gravy()` raises `KeyError: 'U'` — Seqeron skips U/O like X/B/Z/J (= gravy("MKA") = −0.0667). Standard-table translation never produces U/O (UGA/UAG → `*`).
10. **Nucleotide `X` (opt-in)** (review 2026-09, B02 F33) - Biopython 1.88 treats nucleotide `X` as any base (`IUPACData.ambiguous_dna_values['X'] == 'GATC'`): `Seq("GCX").translate()` = `A`, `Seq("ATGNNN").translate()` = `MX` but `Seq("ATGXXX").translate()` and `Seq("TAX").translate()` raise `TranslationError` (possible-stop codon; X not in `GATCRYWSMKHBVDN`). `Translator.Translate(string, code, frame, toFirstStop, allowNucleotideX: true)` reproduces this; default (4-argument) overload still rejects X. Differential: 600 random DNA/RNA/mixed-case strings with X × 27 tables × to_stop (oracle `TestData/GeneticCode/biopython_nucleotide_x.tsv`, 206 X-containing strings translate, rest raise) → 600/600 identical. `DnaSequence`/`RnaSequence` cannot contain X, so `FindOrfs`/`TranslateSixFrames` are unaffected.
11. **Residue-class Mole% (`ProteinSequence.TypePercentage`)** (review 2026-09, B02 F47) - EMBOSS 6.6.0 `pepstats` "Property / Residues / Number / Mole%" table: classes from `Eamino.dat` columns (Tiny A C G S T U; Small A B C D G N P S T U V; Aliphatic A I J L V; Aromatic F H W Y; Non-polar A C F G I J L M P U V W Y; Polar D E H K N O Q R S T Z; Charged B D E H K O R Z; Basic H K O R; Acidic B D E Z — the printed labels omit U/O/J); counts by `embIepCompC` (B → D `(int)(0.5+nB·5.5/9.8)`, rest N; Z → E `(int)(0.5+nZ·6.0/9.9)`, rest Q); Mole% = 100·Number/Residues, Residues = length with one terminal `*` dropped by the reader (`MKRHDE*` → 6, `AK**` → 3, `AAXXHK` → 6). pepstats values: insulin A (21) Tiny 8/38.095, Small 11/52.381, Aliphatic 5/23.810, Aromatic 2/9.524, Non-polar 12/57.143, Polar 9/42.857, Charged 2/9.524, Basic 0, Acidic 2/9.524; LACI_ECOLI (360) 120/33.333, 197/54.722, 137/38.056, 21/5.833, 200/55.556, 160/44.444, 69/19.167, 37/10.278, 32/8.889; `MBBBZZZZXJUO*KH*` (15) 1/6.667, 4/26.667, 1/6.667, 1/6.667, 3/20.000, 10/66.667, 7/46.667, 3/20.000, 4/26.667.

---

## Test Datasets

### From Wikipedia Genetic Code Table
```
Standard Translation Examples:
- ATG (DNA) / AUG (RNA) → M (Methionine)
- GCT/GCC/GCA/GCG → A (Alanine)  
- TAA/TAG/TGA (DNA) → * (Stop)

Six-Frame Translation:
- Forward strand: frames +1, +2, +3
- Reverse complement: frames -1, -2, -3
```

### Known Biological Sequences
```
Human Insulin B Chain (UniProt P01308, positions 25–54):
- DNA: TTCGTGAACCAGCACCTGTGCGGCTCCCACCTGGTGGAAGCTCTGTACCTGGTGTGTGGGGAGCGTGGCTTCTTCTACACACCCAAGACC
- Protein: FVNQHLCGSHLVEALYLVCGERGFFYTPKT
- Length: 30 amino acids
```

---

## Implementation-Specific Notes

### Current Implementation (Translator class)
Based on code review:
1. Supports DNA and RNA input sequences
2. Automatic T→U conversion for DNA input
3. Frame parameter: 0, 1, or 2 (throws for invalid values)
4. Optional `toFirstStop` parameter to terminate at stop codon
5. Six-frame translation returns dictionary with keys -3 to +3 (excluding 0)
6. ORF finding with configurable minimum length and strand search
7. Supports alternative genetic codes via `GeneticCode` parameter

---

## Verification Methodology

All four implemented genetic code tables were verified codon-by-codon against the NCBI Genetic Codes page (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, last update Sep. 23, 2024). Verification method:

1. Each NCBI table encodes 64 codons in the AAs line using the positional encoding:
   - Base1: T×16, C×16, A×16, G×16
   - Base2: (TTTTCCCCAAAAGGGG)×4
   - Base3: (TCAG)×16
2. All 64 codon→amino acid mappings in `CreateStandardCodonTable()` match NCBI Table 1.
3. Differences applied in Tables 2, 3, 11 match the NCBI "Differences from the Standard Code" sections.
4. Start and stop codon sets for each table match the NCBI Starts/AAs lines.

Insulin B chain test data verified against UniProt P01308 (Homo sapiens preproinsulin, Swiss-Prot reviewed entry).

---

## References

1. Wikipedia contributors. "Translation (biology)." Wikipedia, The Free Encyclopedia.
2. Wikipedia contributors. "Reading frame." Wikipedia, The Free Encyclopedia.
3. Wikipedia contributors. "Open reading frame." Wikipedia, The Free Encyclopedia.
4. Elzanowski, A. & Ostell, J. "The Genetic Codes." NCBI Taxonomy Group. Last updated Sep. 23, 2024.
5. UniProt Consortium. "P01308 · INS_HUMAN." UniProt Knowledgebase (Swiss-Prot).
6. Lodish H, et al. (2007). Molecular Cell Biology, 6th ed.
7. Pierce BC (2012). Genetics: A Conceptual Approach.
