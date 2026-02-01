# Gene Prediction Algorithm

## Overview

Gene prediction (gene finding) is the computational process of identifying regions of genomic DNA that encode genes. In prokaryotes, this primarily involves identifying Open Reading Frames (ORFs) combined with analysis of upstream regulatory signals such as the Shine-Dalgarno (SD) sequence.

## Biological Background

### Prokaryotic Gene Structure

In prokaryotes, genes have:
- **Promoter sequences**: -35 box (TTGACA) and -10 box (TATAAT, Pribnow box)
- **Ribosome Binding Site (RBS)**: Contains the Shine-Dalgarno sequence upstream of start codon
- **Start codon**: Usually ATG, but also GTG and TTG
- **Continuous ORF**: No introns in prokaryotes
- **Stop codon**: TAA, TAG, or TGA

### Shine-Dalgarno Sequence

The Shine-Dalgarno (SD) sequence is a ribosome-binding site in bacterial and archaeal mRNA, generally located ~8 bases upstream of the start codon AUG.

**Key characteristics** (Wikipedia, Shine & Dalgarno 1975):
- Consensus sequence: **AGGAGG** (or AGGAGGU in E. coli)
- Shorter variants: GGAGG, AGGAG, GAGG, AGGA
- Located 4-15 nucleotides upstream of start codon
- Optimal spacing: 5-13 bp from SD to start codon (Chen et al. 1994)
- Functions by base-pairing with 3' end of 16S rRNA

**RBS distance constraints** (Chen et al. 1994, Laursen et al. 2005):
- Minimum distance: 4 bp
- Maximum distance: ~15 bp
- Optimal distance: 5-9 bp (species-dependent)

## Algorithm Description

### Simple ORF-Based Gene Prediction

The implementation uses an ORF-based approach:

1. **Find all ORFs** in all six reading frames (3 forward + 3 reverse complement)
2. **Filter by minimum length** (e.g., 100 amino acids)
3. **Require start codon** (ATG, GTG, or TTG)
4. **Terminate at stop codon** (TAA, TAG, TGA)
5. **Generate gene annotations** with positional information

### RBS Detection (FindRibosomeBindingSites)

1. For each ORF found in the sequence:
2. Search upstream region (default: 20 bp window)
3. Look for SD motif matches: AGGAGG, GGAGG, AGGAG, GAGG, AGGA
4. Validate distance constraints (4-15 bp to start codon)
5. Score based on motif length (longer = higher score)

## Implementation Details

### Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| `PredictGenes(dna, minOrfLength, prefix)` | ORF-based gene prediction | O(n) |
| `FindRibosomeBindingSites(dna, window, minDist, maxDist)` | SD sequence detection | O(n × m) where m = ORF count |

### Parameters

**PredictGenes:**
- `dnaSequence`: Input DNA sequence
- `minOrfLength`: Minimum ORF length in amino acids (default: 100)
- `prefix`: Prefix for gene IDs (default: "gene")

**FindRibosomeBindingSites:**
- `upstreamWindow`: Search window size (default: 20 bp)
- `minDistance`: Minimum SD-to-start distance (default: 4 bp)
- `maxDistance`: Maximum SD-to-start distance (default: 15 bp)

### Output

**GeneAnnotation record:**
- GeneId: Unique identifier (e.g., "gene_0001")
- Start/End: 0-based positions
- Strand: '+' or '-'
- Type: "CDS"
- Product: Description (default: "hypothetical protein")
- Attributes: Frame, protein length, translation

## Edge Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty sequence | Returns empty | Implementation |
| No start codon | No ORF found | Definition |
| No stop codon | ORF extends to sequence end (if requireStartCodon=false) | Implementation |
| Short ORFs | Filtered by minOrfLength | Parameter |
| Overlapping ORFs | All reported independently | Implementation |
| Reverse strand | Coordinates adjusted to forward strand | Implementation |

## Limitations

1. **Simple heuristic**: Not comparable to GLIMMER, GeneMark, or other sophisticated gene finders
2. **No intron handling**: Prokaryotic model only
3. **No codon bias analysis**: Does not use codon adaptation index
4. **No promoter scoring**: RBS detection only, no -10/-35 box scoring integration
5. **No training**: Static motif patterns, not trained on organism-specific data

## References

1. Wikipedia: Gene prediction - https://en.wikipedia.org/wiki/Gene_prediction
2. Wikipedia: Shine-Dalgarno sequence - https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence
3. Wikipedia: Ribosome-binding site - https://en.wikipedia.org/wiki/Ribosome-binding_site
4. Shine J, Dalgarno L (1975). "Determinant of cistron specificity in bacterial ribosomes". Nature 254:34-38.
5. Chen H et al. (1994). "Determination of the optimal aligned spacing between the Shine-Dalgarno sequence and the translation initiation codon". Nucleic Acids Research 22(23):4953-4957.
6. Laursen BS et al. (2005). "Initiation of Protein Synthesis in Bacteria". Microbiology and Molecular Biology Reviews 69(1):101-123.
7. Stormo GD et al. (1982). "Characterization of translational initiation sites in E. coli". Nucleic Acids Research 10(9):2971-2996.
