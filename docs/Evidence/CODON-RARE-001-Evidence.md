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
     - AGA, AGG (Arg) are rare codons in E. coli (0.04, 0.02) - Kazusa MG1655
     - CGA (Arg) is rare in E. coli (0.06)
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
Based on Kazusa MG1655 (species=316407):
- **AGA** (Arg): 0.04
- **AGG** (Arg): 0.02
- **CGA** (Arg): 0.06
- **CUA** (Leu): 0.04
- **AUA** (Ile): 0.07 (borderline)
- **UAG** (Stop): 0.07

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

---

# Addendum (2026-06-24): Rare-Codon Cluster / Run Detection

Closes the documented limitation that the per-codon `FindRareCodons` does not detect
consecutive rare-codon **clusters / runs / ramps**. Two published, citable cluster-detection
methods were retrieved this session and implemented as opt-in additions
(`CalculateMinMaxProfile`, `FindRareCodonClusters`). The per-codon default is unchanged.

## Online Sources (retrieved 2026-06-24)

### %MinMax — Clarke & Clark (2008), "Rare Codons Cluster"

**Authority rank:** 1 (peer-reviewed, PLoS ONE).
**Retrieved via:** WebSearch "Clarke Clark 2008 %MinMax codon usage sliding window …" →
WebFetch of the PMC fulltext https://pmc.ncbi.nlm.nih.gov/articles/PMC2565806/ and the PLoS
ONE article page https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003412 ;
summed-window form cross-checked against the %MinMax 2018 paper
https://pmc.ncbi.nlm.nih.gov/articles/PMC5734269/ .

**Key Extracted Points (verbatim / close paraphrase):**

1. **Per-amino-acid quantities:** For the *j*th codon of the *i*th amino acid with *n*
   synonymous codons, the algorithm uses Xij (the actual codon usage frequency), Xmax,i
   ("the usage frequency for the most common codon encoding this amino acid"), Xmin,i ("the
   usage frequency for the least common codon encoding this amino acid"), and Xavg,i ("an
   average usage frequency … calculated for each residue by summing the individual codon
   frequencies and dividing by the number of codons"). These are **per-amino-acid
   comparisons among synonymous codons only**, matching the per-family relative fractions
   already stored in `CodonUsageTable.CodonFrequencies`.
2. **Window formula (2018 summed form):**
   If Σ Xij > Σ Xavg,i then `%Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i)` ;
   if Σ Xij < Σ Xavg,i then `%Min = Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i)` (×100).
3. **Window size:** "%MinMax results are typically averaged over an 18-codon sliding window"
   (the 2008 paper tested 5–30 codons; the 2018 paper states z = 17 typical). Default chosen
   here: **18 codons** (Clarke & Clark 2008).
4. **Cluster interpretation:** "clusters of predominantly rare codons appear as negative
   (%Min) peaks"; −100% = a window encoded using only the rarest synonymous codons, +100% =
   only the most common, 0% = equal to the mean. Each window yields either a %Min or a %Max
   value, not both.

### Sherlocc — Chartier, Gaudreault, Najmanovich (2012)

**Authority rank:** 1 (peer-reviewed, Bioinformatics) + 3 (reference implementation).
**Retrieved via:** WebSearch "Chartier Najmanovich Sherlocc rare codon cluster …" → WebFetch
of the Bioinformatics article page
https://academic.oup.com/bioinformatics/article/28/11/1438/265476 (DOI and window/cluster
rule) and the tool README https://raw.githubusercontent.com/mtthchrtr/sherlocc/master/README.md .

**Key Extracted Points (verbatim):**

1. **Detection window:** "a seven codon-wide window" moves across the sequence/alignment,
   computing average codon usage frequencies at each position.
2. **Cluster rule:** a "seven position-wide window … containing at least four pause positions
   out of seven" is reported as a rare-codon cluster. (README: "Any window with 4 'slow'
   positions or more are considered as rare codon clusters … a centered 7 positions wide
   window.")
3. **'Slow'/pause position:** "Positions with a codon usage frequency average equal or lower
   than the specified threshold … considered as 'slow' positions." Implemented here as the
   same per-codon rare criterion as `FindRareCodons` (frequency strictly below the supplied
   `rareThreshold`, default 0.15).
4. **DOI:** 10.1093/bioinformatics/bts149.

## Documented Corner Cases (cluster methods)

1. **Single-codon amino acids (Met AUG, Trp UGG):** no synonymous spread, so Xmax = Xmin =
   Xavg = Xij and the codon contributes 0 to both %MinMax numerator and denominator — the
   window value stays defined (no divide-by-zero / NaN).
2. **Window longer than the sequence:** no windows / no clusters (Clarke & Clark windowing).
3. **Isolated rare codons:** flagged by per-codon `FindRareCodons` but NOT reported as a
   Sherlocc cluster unless ≥ 4 fall inside a 7-codon window — this is precisely the
   capability the unit adds.
4. **Overlapping qualifying windows:** merged into one maximal cluster so a long rare run is
   reported once (implementation choice; Sherlocc reports cluster regions, not raw windows).

## Worked Examples (E. coli K12, Arg family: CGU 0.38, CGC 0.40, CGA 0.06, CGG 0.10, AGA 0.04, AGG 0.02; Xavg = 1.00/6 = 0.16667, Xmax = 0.40, Xmin = 0.02)

| Input (window) | windowSize | Method | Expected | Derivation |
|----------------|-----------|--------|----------|------------|
| AGA·AGA·AGA | 3 | %MinMax | −86.363636…% | −(0.16667−0.04)/(0.16667−0.02)×100 |
| CGC·CGC·CGC | 3 | %MinMax | +100% | (0.40−0.16667)/(0.40−0.16667)×100 |
| CUG·AGA | 2 | %MinMax | +36.470588…% | (0.54−0.33333)/((0.50−0.16667)+(0.40−0.16667))×100 |
| 7× AGA | 7/≥4 | Sherlocc RCC | 1 cluster, codons 0–6, 7 rare | 7 rare ≥ 4 |
| 3 AGA + 4 CGC | 7/≥4 | Sherlocc RCC | no cluster | 3 rare < 4 |
| 4 AGA + 3 CGC | 7/≥4 | Sherlocc RCC | 1 cluster, 4 rare | exactly 4 ≥ 4 |

## References (cluster methods)

- Clarke TF, Clark PL. 2008. Rare Codons Cluster. PLoS ONE 3(10):e3412.
  https://doi.org/10.1371/journal.pone.0003412
- Rodriguez A, Wright G, Emrich S, Clark PL. 2018. %MinMax: A versatile tool for calculating
  and comparing synonymous codon usage and its impact on protein folding. Protein Science.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5734269/
- Chartier M, Gaudreault F, Najmanovich R. 2012. Large-scale analysis of conserved rare codon
  clusters suggests an involvement in co-translational molecular recognition events.
  Bioinformatics 28(11):1438–1445. https://doi.org/10.1093/bioinformatics/bts149
- Sherlocc reference implementation (Chartier & Najmanovich).
  https://github.com/mtthchrtr/sherlocc (README accessed 2026-06-24)
