# CODON-CAI-001: Codon Adaptation Index (CAI) Calculation Evidence

## Test Unit
- **ID:** CODON-CAI-001
- **Area:** Codon
- **Algorithm:** CAI Calculation

## Sources

### Primary Sources

1. **Wikipedia - Codon Adaptation Index**
   - URL: https://en.wikipedia.org/wiki/Codon_Adaptation_Index
   - Key concepts:
     - CAI is "the most widespread technique for analyzing codon usage bias"
     - Measures "deviation of a given protein coding gene sequence with respect to a reference set of genes"
     - Used as "a quantitative method of predicting the level of expression of a gene based on its codon sequence"
     - Reference set ideally composed of highly expressed genes

2. **Sharp, P.M. & Li, W.H. (1987)** - Original CAI Paper
   - Title: "The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications"
   - Journal: Nucleic Acids Research, 15(3): 1281-1295
   - PMC: 340524, PMID: 3547335
   - DOI: 10.1093/nar/15.3.1281

3. **Xia, X. (2007)** — "An Improved Implementation of Codon Adaptation Index", Evolutionary Bioinformatics 3:53-58
   (PMC2684136; this entry was previously mis-attributed to "Jansen, Bauer & Stadler 2003", which is a different NAR paper)
   - Retrieved 2026-06-24 from PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC2684136/; statement re-confirmed 2026-09-28 via search snippet of journals.sagepub.com/doi/full/10.1177/117693430700300028
   - Quotes the original Sharp & Li (1987) rule **verbatim** and gives the reason:
     > "The original paper proposing CAI (Sharp and Li, 1987) specifically stated that codon
     > families containing a single codon (e.g. AUG and UGG in the standard genetic code)
     > should be excluded in computing CAI."
     > "Note that, for such codons (e.g. AUG and UGG in the standard genetic code), their
     > corresponding w value will always be 1 regardless of codon usage bias of the gene. If a
     > gene happens to use a high proportion of methionine and tryptophan, then it will have a
     > high CAI value even if its codon usage is not at all biased."
   - This is the authoritative basis for the **single-codon amino-acid exclusion** implemented as
     the `excludeSingleCodonAminoAcids` mode of `CalculateCAI` (default `true` since review 2026-09).

### Mathematical Definition (from Wikipedia)

**Relative Adaptiveness (w_i):**
```
w_i = f_i / max(f_j)  where i,j ∈ [synonymous codons for amino acid]
```

Where:
- f_i = observed frequency of codon i
- max(f_j) = frequency of the most frequent synonymous codon for that amino acid

**CAI Calculation:**
```
CAI = (∏_{i=1}^{L} w_i)^{1/L}
```

Equivalent to:
```
CAI = exp((1/L) × Σ ln(w_i))
```

Where L = number of codons (excluding stop codons per implementation)

### Key Properties (from Sharp & Li 1987 / Wikipedia)

1. **Range:** 0 ≤ CAI ≤ 1
   - CAI = 1 when all codons are the most frequent for their amino acids
   - CAI = 0 when any codon has zero frequency in the reference set, or when input is empty
   - 0 < CAI < 1 for typical genes

2. **Geometric Mean:** CAI uses geometric mean, which is sensitive to low values
   - A single rare codon significantly lowers overall CAI

3. **Single-Codon Amino Acids:** Methionine (AUG) and Tryptophan (UGG) have w=1.0 always
   - Only one codon exists, so it's always the "most frequent"
   - **Canonical rule (Sharp & Li 1987; Xia 2007):** such single-codon families
     **should be EXCLUDED** from the CAI geometric mean, because w≡1 regardless of bias and
     including them inflates CAI for Met/Trp-rich genes (see Source 3, retrieved verbatim).
   - Default `CalculateCAI(seq, table)` *excludes* Met/Trp (review 2026-09, F12);
     `excludeSingleCodonAminoAcids: false` scores them with w=1.0 (EMBOSS-style opt-in).

4. **Stop Codons:** Excluded from CAI calculation
   - Source: Sharp & Li (1987) — stop codons do not encode amino acids

### Documented Edge Cases

1. **Empty Sequence:**
   - No codons to evaluate → CAI = 0 (by convention)
   - Source: Implementation convention

2. **Sequence with only Met/Trp:**
   - No scored codon → CAI = 0 (CodonW 1.4.4 binary: `ATGTGG` → 0.000); opt-in inclusion → 1.0

3. **All Optimal Codons:**
   - Every codon is the most frequent for its amino acid → CAI = 1.0
   - Source: CAI definition

4. **All Rare Codons:**
   - Low w values for all codons → CAI approaches 0
   - Source: Geometric mean properties

5. **Different Organisms:**
   - Same sequence has different CAI values for different organisms
   - Source: Codon usage varies by organism

### Implementation Notes (review 2026-09)

`CodonOptimizer.CalculateCAI` delegates to the canonical `CodonUsageAnalyzer.CalculateCai` core:
stop codons and single-codon families (genetic-code dependent) not scored; `w = value / family max`;
`w < 0.0001 → 0.01` (CodonW `cai_out`); family without data not scored; non-nucleotide triplets skipped
frame-preservingly; DNA/RNA any case. The former `1e-6` clamp and the former silent drop of w=0 codons
(CodonUsageAnalyzer) were unsourced and are removed (F13).

### Reference implementations opened (2026-09-28)

- **CodonW 1.4.4** (original tarball, compiled): `codon_us.c` `cai_out` — "Non-synonymous codons and
  termination codons (genetic code dependent) are excluded … these codons have fitness of zero (<.0001)
  are adjusted to 0.01"; w-building: "if a codon is absent then adjust its frequecy to 0.5";
  `codonW.h` `cai[]` E. coli w = Biopython `SharpEcoliIndex` (identical 61 values).
- **seqinr** `R/cai.R`, `man/cai.Rd` (raw.githubusercontent.com/cran/seqinr): excludes stops and
  singulets, `zero.threshold = 0.0001, zero.to = 0.01` ("default is from Bulmer (1988)"), "intended to
  work exactly as in the program codonW".
- **Biopython 1.88** `Bio/SeqUtils/__init__.py` `CodonAdaptationIndex`: 0.5 for codons absent from the
  reference ("Following the description in the original paper"); `calculate` skips ATG/TGG; scores stop
  codons present in the index (quirk, not followed). **Biopython 1.79** `CodonUsage.py` `cai_for_gene`
  divides by `cai_length - 1.0` (bug, not followed); `CodonUsageIndices.SharpEcoliIndex`.
- **EMBOSS** `ajcod.c` `ajCodCalcCaiSeq`: L = all codons (Met/Trp and stops included), w = 0 codons
  contribute nothing to the sum but count in L — basis of the `false` opt-in only.
- Sharp & Li 1987 full text: not reachable (academic.oup.com / PMC blocked); rules taken from the above
  and the Xia 2007 snippet.

### Numerical cross-check (2026-09-28)

- Python port of CodonW `cai_out` vs CodonW binary: 208 genes × 8 codes (NCBI 1,2,3,4,5,6,9,10) × 4 w
  tables (E. coli built-in + 3 random with ~15 % zeros) = 6 656 cases, 0 mismatches (3-dp rounding).
- C# vs port: 832 inputs (random DNA/RNA/lower/IUPAC) × 27 NCBI tables, max |Δ| 5.6e-16; CodonOptimizer
  (Standard): 832/832. Before the fix: CodonUsageAnalyzer 556/832 and CodonOptimizer 556/832
  (zero-w tables) mismatching; default CodonOptimizer mode 698/832.
- C# vs Biopython 1.88 `CodonAdaptationIndex` (40 indices from small reference sets, 400 genes): max |Δ| 0.

## Test Datasets

### Reference Codon Tables (Kazusa MG1655, species=316407)

**E. coli K12 Leucine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CUG | 0.50 | 1.00 (optimal) |
| UUA | 0.13 | 0.26 |
| UUG | 0.13 | 0.26 |
| CUU | 0.10 | 0.20 |
| CUC | 0.10 | 0.20 |
| CUA | 0.04 | 0.08 (rare) |

**E. coli K12 Arginine Codons:**
| Codon | Frequency | Relative Adaptiveness |
|-------|-----------|----------------------|
| CGC | 0.40 | 1.00 (optimal) |
| CGU | 0.38 | 0.95 |
| CGG | 0.10 | 0.25 |
| CGA | 0.06 | 0.15 |
| AGA | 0.04 | 0.10 (rare) |
| AGG | 0.02 | 0.05 (rare) |

### Hand-Calculated Test Cases

**Test Case 1: Single Met (AUG)**
- Not scored (single-codon family) → CAI = 0; opt-in inclusion: w = 1 → 1.0

**Test Case 2: CUG-CCG-ACC (E. coli)**
- CUG: w = 0.50/0.50 = 1.0 (Leu optimal)
- CCG: w = 0.53/0.53 = 1.0 (Pro optimal)
- ACC: w = 0.44/0.44 = 1.0 (Thr optimal)
- CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0

**Test Case 3: CUA-CCA-ACA (E. coli rare)**
- CUA: w = 0.04/0.50 = 0.08 (Leu rare)
- CCA: w = 0.19/0.53 = 0.3585 (Pro suboptimal)
- ACA: w = 0.13/0.44 = 0.2955 (Thr suboptimal)
- CAI = (0.08 × 0.3585 × 0.2955)^(1/3) = 0.1980

### Hand-Calculated Test Cases — Exclusion Mode (`excludeSingleCodonAminoAcids: true`, the default; "Inclusive" = opt-in `false`)

All values E. coli K12 (Kazusa species=316407); CUA(Leu) w = 0.04/0.50 = 0.08; AUG(Met) and UGG(Trp) excluded.

**Test E1: `AUGUGG` (Met + Trp only)**
- Both codons excluded → no scored codons (L=0) → CAI = 0 (contrast: inclusive default = 1.0).

**Test E2: `AUGCUACUA` (Met + 2×CUA)**
- Inclusive default: exp((ln1 + ln0.08 + ln0.08)/3) = 0.18566355334451112
- Exclusive: Met dropped → exp((ln0.08 + ln0.08)/2) = **0.08** exactly.

**Test E3: `AUGUGGCUA` (Met + Trp + CUA)**
- Inclusive default: exp((ln1 + ln1 + ln0.08)/3) = 0.43088693800637673
- Exclusive: Met+Trp dropped → exp(ln0.08 / 1) = **0.08** exactly.

**Test E4: `CUGCUA` (no single-codon AA)**
- Inclusive == Exclusive = (1.0 × 0.08)^(1/2) = 0.28284271247461906 (flag has no effect).

## Assumptions

- Empty sequence / no scored codon returns 0 (CodonW prints 0.000).
- A family with no reference data is not scored.
- All codon frequency tables verified against Kazusa database (March 2026).

## References

- Sharp, P.M. & Li, W.H. (1987). Nucleic Acids Res. 15(3):1281-1295
- Xia, X. (2007). Evol. Bioinform. 3:53-58 (PMC2684136)
- Peden, J.F. (1999). CodonW 1.4.4; Bulmer, M. (1988). J. Evol. Biol. 1:15-26
- Wikipedia: Codon Adaptation Index
- Kazusa Codon Usage Database: https://www.kazusa.or.jp/codon/
