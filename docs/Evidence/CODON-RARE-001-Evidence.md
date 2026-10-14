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
     - AGA, AGG (Arg) are rare codons in E. coli (0.04, 0.02) - Kazusa K-12 W3110 (species=316407)
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
Based on Kazusa K-12 W3110 (species=316407):
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

---

# Addendum (2026-09-28, review campaign B02)

## Sources opened this session
- **python-codon-tables 0.1.18** (PyPI wheel, Edinburgh Genome Foundry; tables downloaded from
  Kazusa `showcodon.cgi?aa=1&style=N&species=<taxid>`): `e_coli_316407.csv`,
  `s_cerevisiae_4932.csv`, `h_sapiens_9606.csv`. All 64 relative frequencies of each built-in
  preset (`EColiK12`, `Yeast`, `Human`) are identical. Kazusa 316407 is E. coli K-12 **W3110**
  (earlier text said MG1655; corrected).
- **Kane JF (1995)** Curr Opin Biotechnol 6(5):494–500 — WebSearch snippet: "The rarest codons
  AGG, AGA, CUA, AUA, CGA, and CCC of E. coli". With the W3110 fractions all six are < 0.15.
- **CHARMING (Wright … Clark 2022, Protein Sci, doi:10.1002/pro.4223)** reference code from the
  %MinMax authors' lab: `raw.githubusercontent.com/wrightgs/CHARMING/main/CHARMING.py`
  (`calculateMinMax`), `README.md`, `ScerCUB.txt`, `EcolCUB.txt`. Findings:
  1. The window formula is Σ-based exactly as implemented (averages over the window cancel).
  2. The synonymous families include the stop family `'*': ['TGA','TAA','TAG']`.
  3. The codon usage input is **frequency per thousand codons** ("TTT 26.1 / TCT 23.5 / TAT
     18.8", i.e. Kazusa 4932 per-thousand), not per-amino-acid fractions. The 2026-06 claim that
     the published quantities "match the per-family relative fractions already stored in
     `CodonUsageTable.CodonFrequencies`" is therefore corrected: relative fractions reproduce
     the reference only for single-residue windows (yeast CUG·AGA, w=2: 33.99 reference vs
     58.62 with fractions).
- Not reachable (proxy): journals.plos.org, academic.oup.com, codons.org, onlinelibrary.wiley.com,
  arxiv.org, bmcgenomics.biomedcentral.com, rdrr.io. The Sherlocc GitHub source fetch was not
  permitted in this session, so the Sherlocc "equal or lower than the threshold" wording quoted
  above (2026-06) could not be re-checked against code (see open item below).

## Ambiguous triplets
EMBOSS `ajCodSetTripletsS` ("Skips triplets with ambiguity codes and any incomplete triplet at the
end") is the codon-set contract already adopted by `CodonUsageAnalyzer.CountCodons`. A triplet
such as NNN has no row in a Kazusa/cusp table, so it has no usage frequency; it is skipped by all
three rare-codon screens (frame preserved) instead of being reported as a frequency-0 "rare codon".

## Open item (2026-09) — resolved 2026-10 (review FIN-B02 WP7, A1-1)
Sherlocc README (quoted 2026-06): slow positions have frequency "equal or lower than the specified
threshold"; `FindRareCodonClusters` used strict `<` (shared with `FindRareCodons`). With the
E. coli preset this matters at the 0.15 default for UCU/UCC/UCG/AGU/GUA/GGG (exactly 0.15).

# Review 2026-10 (FIN-B02 WP7) — Sherlocc source; per-thousand presets

## Sherlocc 1.0 reference implementation (opened)
`git clone https://github.com/mtthchrtr/sherlocc` (commit b3b4d15987c5f8205d1957db3a757dbd390d5049,
2020-06-25; files `sherlocc.pl` 1223 lines, `README.md`, `files/cut`); README also via
`raw.githubusercontent.com/mtthchrtr/sherlocc/master/README.md`. Verbatim:
- Threshold prompt (l. 26–28): `print "\nEnter Sherlocc threshold:\n"; chomp($avrg_lim = <STDIN>);`
  — README: "A rare codon frequency threshold (i.e.: 14,15,16,etc)".
- Frequencies (stage 2, l. 449–522): per-organism rows of `files/cut`
  (`org_id=…|transl_tbl=…` then `TTT5.5|TCT19.2|…` — CUTG frequency per thousand codons, keyed by
  Kazusa/NCBI taxid); each aligned residue gets its codon's per-thousand value, gaps/unknown → `x`.
- Pause (stage 3, l. 775–825): for every alignment column `i`, sum the values of columns
  `i−3 … i+3` over all rows (sequences); out-of-range columns and `x` cells are not counted;
  `$tresh=($row*(7))*(0.75)`; `if($freq_count >= $tresh)` → `$avrg=$freq_sum/$freq_count`;
  **`if($avrg <= $avrg_lim){ #We have a pause`**; otherwise (or column not representative) not a pause.
  README: "a centered 7 positions wide window … windows occupied by alignment gaps … by 25% or more
  are excluded … positions with a codon usage frequency average equal or lower than the specified
  threshold … considered as 'slow' positions".
- Clusters (l. 834–960): `for($i=0; $i<=($#col_type-6); $i++)` count pauses in columns `i … i+6`;
  `if($win_pause_count >= 4)` → cluster `first pause – last pause` of that window; consecutive
  clusters merged only when `this_start >= start && this_start <= stop` (overlap); each reported
  with the mean column average over its span (`$cluster_avrg`), 1-based.
- No P-value / randomisation / shuffling code in `sherlocc.pl` (grep `rand|shuffl|pval|simul`:
  none); the conservation aspect is the averaging over all family members of the Pfam alignment
  plus the 75 % representativeness filter. P-values are only in the paper (not opened here).

Decision: `FindRareCodonClusters` (TestSpec cites Chartier 2012 as its canonical source) now uses the
source comparison `≤` by default; the previous strict `<` stays available (`inclusiveThreshold:
false`). `FindRareCodons` is not a Sherlocc unit (sources: Kazusa / GenRCA / Shu 2006, TestSpec M06
strict `<`), so its default stays strict, with an opt-in `inclusiveThreshold: true`. The window
averaging / alignment / per-thousand scale differences are recorded as a kept simplification
(algorithm doc §5.3) and are the subject of audit item A1-2.

## Per-thousand presets (A1-3)
- %MinMax input scale: CHARMING `calculateMinMax` sums the table values directly; its tables are
  per-thousand (`ScerCUB.txt` "TTT 26.1"), so per-thousand is the sourced default when available.
- Yeast (Kazusa 4932): CHARMING `ScerCUB.txt` (commit 2c75d1b, one decimal). Per-amino-acid ratios
  vs python-codon-tables 0.1.18 `s_cerevisiae_4932`: 62/64 equal at 2 dp; stop family UAA 0.478 vs
  0.47, UAG 0.217 vs 0.23 (one-decimal rounding of 1.1/0.5/0.7 — e.g. 1.05/0.52/0.68 give
  0.47/0.23/0.30).
- Human (Kazusa 9606): PyPI `codon-optimizer` 0.2.0 `codon_optimizer/reference_seq/Homo_sapiens_codon_frequency.csv`
  (AmAcid, Codon, Number, /1000, Fraction; tool help cites kazusa.or.jp/codon); Σ Number =
  40 662 582 = the Kazusa 9606 codon total; ratios equal python-codon-tables `h_sapiens_9606` at all
  64 codons (2 dp); stored as 1000·n/40 662 582 (its /1000 column agrees within 0.005).
- E. coli (Kazusa 316407): per-thousand column not obtainable — kazusa.or.jp, genscript.com,
  detaibio.com, dnahive.fda.gov blocked (proxy / WebFetch ENOTFOUND); GitHub code search not
  available; python-codon-tables (0.1.3.3, 0.1.9, 0.1.12, 0.1.18, GitHub source),
  codon-harmony, CodonU, pycodon, freqgen, codonopt, codon-tools, dnachisel, kuma
  (`kuma_core/kuro/resources/codon_tables/ecoli.json`, fractions only) carry fractions only.
  EMBOSS 6.6 `Eecoli.cut` (CUTG146, 5045 CDS) differs from the 316407 fractions at 9 codons,
  CHARMING `EcolCUB.txt` at 43 — different entries, not used. `EColiK12.PerThousand` = null.
- Independent check: Python port of the formula and CHARMING's own `calculateMinMax` (isolated)
  over the same tables — yeast `CUGAGA` w2 33.99209486166008 (CHARMING 33.99); 20-codon gene w18
  12.68791340950089 / 17.028571428571425 / 15.991062879029688 (12.69/17.03/15.99); human gene w18
  28.796581754809075 / 23.264123403923747 / 16.195449171123403 (28.8/23.26/16.2); human `UUACGU` w2
  −96.40404591634571 (−96.4); human fractions (CodonFrequencies basis) 19.650655021834133 /
  13.878080415045439 / 13.803548247511891, `UUACGU` −94.54545454545456.

# Review 2026-10 (FIN-B02 WP8) — Sherlocc alignment analysis (A1-2)

## Source re-read in full
`sherlocc.pl` (commit b3b4d15, 1223 lines) read end to end. Stage 1 (l. 95–445) retrieves, per Pfam
row, UniProt/ENA records over HTTP and writes `u_id/org_id/pfam_ori/pfam_aa/transl_tbl/coding_seq`
(l. 369–374). Stage 2 (l. 447–537) gives every aligned residue the per-thousand usage of its codon in
its OWN organism's CUTG row (`$cu_hash{$org_id}`, l. 497), non-residue/unknown → `x`. Stage 3:
- Averaging + filter (l. 772–819): `$tresh=($row*(7))*(0.75)`; for column `i`, cells of columns
  `i−3…i+3` of all rows are summed (`$freq_count`), out-of-range / `x` cells skipped;
  `if($freq_count >= $tresh)` → `$avrg=$freq_sum/$freq_count`, pause `if($avrg <= $avrg_lim)`; else
  the column is "NA" (not representative). This is the whole "conservation" filter — no separate
  per-sequence or identity filter exists.
- Clusters (l. 838–854): windows `i…i+6`, `$win_pause_count >= 4` → `first pause – last pause`.
- Merge (l. 894–961): next cluster extends the current one when `this_start` ∈ [start, stop]; when
  `this_start > stop` a new cluster starts (abutting clusters stay separate).
- Cluster mean (l. 875–879): `$freq_sum+=$col_freq[$k]` over the span — an NA column inside a span
  is an undefined value (0). Middle point (l. 1137–1140) `stop − (len/2)`, decimals stripped;
  fraction (l. 1154) `total_cluster_length / alignment_length`; `clusters_freq_avrg` = mean of cluster
  means.
- FoldIndex block (l. 683–744, 997–1103): Uversky/Prilusky unstructured-region prediction — a separate
  output, not ported.
- P-values: none (no `rand`, shuffle, P-value or significance code; `grep -n -i
  "p-val|pval|random|shuffl|signific|rand("` → only the comment "significant windows (clusters)" l. 837/848).

## Oracle run
perl 5.38 + libwww-perl (apt) in an isolated copy; stages 2+3 on hand-written stage-1 records (the
only change to the script: l. 15 `cleanup("stage1")` and l. 17 `mkdir "stage1"` removed so the records
survive; stage 1 needs live UniProt/ENA). `files/cut`: `org_id=4932` / `org_id=9606` rows with the
Yeast / Human preset per-thousand values. Results (`results/out/results`, HTML "Average" row):
- A (yeast, 3×20, threshold 15): pauses 9–12 (1-based); `clusters=9-12|9.83928571428571`,
  `total_nb_of_pauses=4`.
- B (yeast, 2×18, threshold 15, two fully gapped columns): column 10 (1-based) "NA";
  `clusters=6-13|5.265625` (the NA column counted as 0), `total_nb_of_pauses=7`.
- C (yeast/human/yeast, 3×34, threshold 16): `clusters=4-10|9.20962052101093,21-27|9.30926559372343`,
  `clusters_freq_avrg=9.25944305736718`, `total_nb_of_pauses=14`.
An independent Python port of l. 772–961 reproduces every 2-dp column average and every cluster
value (unrounded values used in the tests at 1e-12).

## Cluster P-values — BLOCKED
The paper (Chartier, Gaudreault & Najmanovich 2012, doi:10.1093/bioinformatics/bts149, PMID 22467916)
reports P-values (WebSearch snippet: P22 tailspike cluster "P-value of 1.72×10^−8 at a threshold of
15"), but the model is not in `sherlocc.pl` and the full text could not be opened: academic.oup.com
(…/article/28/11/1438/265476 and /266932), pmc.ncbi.nlm.nih.gov, pubmed.ncbi.nlm.nih.gov,
eutils.ncbi.nlm.nih.gov, www.ebi.ac.uk/europepmc, api.semanticscholar.org, web.archive.org,
ui.adsabs.harvard.edu, journals.plos.org, arxiv.org — proxy `CONNECT tunnel failed, response 403` /
WebFetch `ENOTFOUND`; WebSearch snippets give no null model. No P-value is reported.

## F48 — genetic-code awareness of the rare-codon APIs (review 2026-09, audit A2-1)

Sources opened: dnachisel 3.2.16 `biotools/biotables.get_backtranslation_table` and
`builtin_specifications/EnforceTranslation.py` (`genetic_table`), Biopython 1.88
`Bio.Data.CodonTable` tables 2 / 4 / 11 / 27 / 28 / 31.

- `FindRareCodons`: rarity is a per-codon usage-table lookup, so the *selection* is code-independent;
  the genetic code only names the amino acid. Added optional `geneticCode` (default Standard):
  on `AUGCGUAGACGCAGGUGAAUAUAA` at threshold 0.30 the reported positions are identical under tables 1
  and 2, while `UGA` is reported as `*` under table 1 and `W` under table 2, and `AGA` as `*` under
  table 2.
- `CalculateMinMaxProfile`: the %MinMax sums are per-family (Xavg, Xmax, Xmin), so the profile is
  code-dependent. Added optional `geneticCode` (default Standard, the Clark-lab reference
  behaviour); under table 2 the Arg family is CGU/CGC/CGA/CGG (AGA/AGG are stops) and Met/Trp become
  two-codon families, which changes the profile while the window count is unchanged.
- `FindRareCodonClusters` / `FindConservedRareCodonClusters` (Sherlocc): no translation and no
  synonymous family is involved — pause marking compares a per-codon (or per-column average) usage
  value with the threshold — so these methods are genetic-code independent and no parameter was
  added. Documented in the XML docs and `Rare_Codon_Detection.md`.
