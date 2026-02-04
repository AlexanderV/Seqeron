# Evidence Document: CODON-RARE-001

## Test Unit
**ID:** CODON-RARE-001  
**Title:** Rare Codon Detection  
**Method Under Test:** `CodonOptimizer.FindRareCodons(string, CodonUsageTable, double)`

---

## Sources Consulted

### Primary Sources

1. **Wikipedia - Codon usage bias**
   - URL: https://en.wikipedia.org/wiki/Codon_usage_bias
   - Key findings:
     - Rare codons correspond to low-abundance tRNAs
     - Translation elongation rates are slower at rare codons
     - Consecutive rare codons can inhibit translation
     - Rare codons have biological significance in protein folding and translational regulation

2. **GenScript - GenRCA Rare Codon Analysis Tool**
   - URL: https://www.genscript.com/tools/rare-codon-analysis
   - Key findings:
     - Comprehensive evaluation of codon usage preferences
     - Supports 31 codon preference indices
     - Reference: Fan et al. (2024) BMC Bioinformatics
     - Uses Kazusa Codon Usage Database as reference

3. **Kazusa Codon Usage Database**
   - URL: https://www.kazusa.or.jp/codon/
   - Key findings:
     - 35,799 organisms with codon usage data
     - Data sourced from NCBI GenBank
     - Standard reference for codon frequencies

4. **Shu et al. (2006) - PMC6032470**
   - Title: "Inhibition of Translation by Consecutive Rare Leucine Codons in E. coli"
   - Key findings:
     - CUA (Leu) is a rare codon in E. coli with frequency 0.04
     - AGA, AGG (Arg) are rare codons in E. coli (0.07, 0.04)
     - CGA (Arg) is rare in E. coli (0.07)
     - Five consecutive rare CUA codons cause ~3-fold inhibition of translation
     - Rare codons at 5' end have stronger effect than internal positions

5. **Sharp & Li (1987) - PMC340524**
   - Title: "The codon adaptation index - a measure of directional synonymous codon usage bias"
   - Key findings:
     - CAI measures directional codon usage bias
     - Rare codons have low relative adaptiveness values
     - Frequency threshold approach is valid for identifying rare codons

---

## Algorithm Behavior

### Definition of Rare Codons
A codon is considered "rare" when its usage frequency falls below a specified threshold relative to the codon usage table of the target organism.

### Standard E. coli K12 Rare Codons (frequency < 0.10)
Based on the implementation's EColiK12 table:
- **AGA** (Arg): 0.07
- **AGG** (Arg): 0.04
- **CGA** (Arg): 0.07
- **CUA** (Leu): 0.04
- **AUA** (Ile): 0.11 (borderline)
- **UAG** (Stop): 0.09

### Default Threshold
The implementation uses a default threshold of 0.15, which identifies codons occurring at less than 15% relative frequency within their synonymous codon family.

---

## Edge Cases Documented

### From Sources

1. **Empty sequence** → returns empty enumerable (no codons to analyze)
2. **Sequence not divisible by 3** → incomplete codons are ignored
3. **No rare codons present** → returns empty enumerable
4. **All codons are rare** → returns all positions
5. **Mixed DNA/RNA input** → T is converted to U internally
6. **Threshold = 0** → all codons reported as rare (except those with freq 0)
7. **Threshold = 1** → all codons reported as rare
8. **Unknown codons** → codons not in table get frequency 0, always reported

### Implementation-Specific

1. **Position reporting** → nucleotide position (0-indexed * 3), not codon index
2. **Amino acid reporting** → translated amino acid included in result tuple
3. **Frequency reporting** → actual frequency from table included in result

---

## Test Dataset Recommendations

### From Literature

1. **E. coli rare arginine cluster**: "AGA AGG" - documented inhibition effect
2. **E. coli rare leucine**: "CUA CUA CUA CUA CUA" - 3-fold inhibition
3. **Known rare codon genes**: Human genes expressed in E. coli often contain problematic rare codons

### Constructed Test Cases

| Test Case | Sequence | Expected Rare Codons | Rationale |
|-----------|----------|---------------------|-----------|
| Single rare | AUGAGA | AGA at pos 3 | Basic detection |
| Multiple rare | AUGAGAAGG | AGA at 3, AGG at 6 | Multiple hits |
| No rare | AUGCUG | None | Common codons only |
| All rare | AGAAGGCGA | All positions | Stress test |
| Mixed | CUGCUA | CUA at pos 3 | CUG common, CUA rare |

---

## Invariants

1. **Position validity**: All reported positions are multiples of 3
2. **Position bounds**: 0 ≤ position ≤ sequence.Length - 3
3. **Codon validity**: Reported codon length is always 3
4. **Frequency range**: 0 ≤ frequency ≤ 1.0
5. **Threshold behavior**: All reported codons have frequency < threshold
6. **Amino acid consistency**: Reported AA matches standard genetic code
7. **Determinism**: Same input produces same output

---

## Known Limitations

1. Implementation does not consider codon context (consecutive rare codons)
2. Does not weight by position (5' vs internal)
3. Single threshold for all amino acids (no per-AA analysis)
4. Does not report tRNA availability, only frequency

---

## References

- Sharp PM, Li WH. The codon adaptation index. Nucleic Acids Res. 1987;15(3):1281-1295.
- Shu P, et al. Inhibition of translation by consecutive rare leucine codons in E. coli. Gene Expr. 2006;13(2):97-106.
- Fan K, et al. GenRCA: a user-friendly rare codon analysis tool. BMC Bioinformatics. 2024;25:309.
- Plotkin JB, Kudla G. Synonymous but not the same: causes and consequences of codon bias. Nat Rev Genet. 2011;12(1):32-42.
