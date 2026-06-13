# Evidence Artifact: CODON-STATS-001

**Test Unit ID:** CODON-STATS-001
**Algorithm:** Codon Usage Statistics (GetStatistics: codon counts, RSCU, ENC, GC1/GC2/GC3, GC3s, total codons; CalculateCai; E. coli / human reference tables)
**Date Collected:** 2026-06-13

---

## Online Sources

### Sharp & Li (1987) — Codon Adaptation Index (primary paper)

**URL:** https://doi.org/10.1093/nar/15.3.1281 (PubMed: https://pubmed.ncbi.nlm.nih.gov/3547335/)
**Accessed:** 2026-06-13 (retrieved via WebSearch "Sharp Li 1987 codon adaptation index CAI formula"; the full-text PDF at academic.oup.com/nar/article-pdf/15/3/1281 returns a tokened redirect and could not be opened directly, so the w-value table was taken from the Biopython reproduction below and the formula from Wikipedia/seqinr/CodonW, all retrieved this session)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **CAI definition:** CAI is the geometric mean of the relative-adaptiveness (w) values of the codons of a gene. Range 0–1 (1 = uses only the most-frequent synonym of every amino acid).
2. **Relative adaptiveness:** `w_i = f_i / f_max`, the ratio of the usage of codon i to the usage of the most-abundant synonymous codon for the same amino acid.
3. **Reference set:** w values for E. coli were derived from 27 very highly expressed genes (ribosomal proteins, outer-membrane proteins, elongation factors); reported in Table 1 ("Values of RSCU and w for codons in very highly expressed genes from E.coli and yeast").

### Wikipedia — Codon Adaptation Index (formula transcription, cites Sharp & Li 1987)

**URL:** https://en.wikipedia.org/wiki/Codon_Adaptation_Index
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 4 (Wikipedia citing the primary Sharp & Li 1987)

**Key Extracted Points:**

1. **Weight formula (verbatim):** `w_i = f_i / max(f_j)` where `max(f_j)` is "the frequency of the most frequent synonymous codon for that amino acid."
2. **CAI formula (verbatim):** `CAI = (∏ w_i)^(1/L)`, equivalently `CAI = exp[(1/L) Σ ln(w_i)]`, with L the number of codons.

### seqinr `cai` reference documentation (R reference implementation)

**URL:** https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, seqinr R package)

**Key Extracted Points:**

1. **w:** "The relative adaptiveness (w) of each codon is the ratio of the usage of each codon, to that of the most abundant codon for the same amino acid."
2. **CAI:** "The CAI index is defined as the geometric mean of these relative adaptiveness values"; computed via natural-log summation.
3. **Excluded codons:** "Non-synonymous codons and termination codons are excluded from the calculation (genetic code dependent)" — i.e. single-codon amino acids Met (ATG) and Trp (TGG) plus the stop codons.

### CodonW codon-usage indices documentation (Peden) — CAI and GC3s

**URL:** https://codonw.sourceforge.net/Indices.html
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, CodonW)

**Key Extracted Points:**

1. **CAI exclusions (verbatim):** "Non-synonymous codons and termination codons (dependent on genetic code) are excluded."

### Peden (1999) — *Analysis of Codon Usage* thesis (definitive CodonW reference)

**URL:** https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf
**Accessed:** 2026-06-13 (WebFetch returned the binary PDF; extracted with `pdftotext`. Definition located at §1.8.2.1.3, line 2063, and the user-manual option 6, line 3257)
**Authority rank:** 1–3 (M.Sc./Ph.D. thesis defining the CodonW indices)

**Key Extracted Points:**

1. **GC3s definition (verbatim, §1.8.2.1.3):** "The index GC3s, is the frequency of G or C nucleotides present at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)."
2. **GC3s computation (verbatim, manual option 6):** "This option calculates the fraction of codons synonymous at the third codon position, and having guanine or cytosine at that third position."
3. **ENc range / clamp (verbatim):** "20 (when effectively only a single codon is used for each amino-acid) and 61 (when codons are used randomly). If ENc is greater than 61 ... it is corrected to 61."

### EMBOSS `cusp` documentation — GC position statistics

**URL:** https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation, EMBOSS)

**Key Extracted Points:**

1. **GC fields (verbatim labels):** "#Coding GC" (overall), "#1st letter GC" (codon position 1), "#2nd letter GC" (position 2), "#3rd letter GC" (position 3).
2. **Fraction field:** "the proportion of usage of the codon among its redundant set" (used to confirm RSCU-style per-amino-acid normalization).

### Biopython `SharpEcoliIndex` (reference implementation of Sharp & Li 1987 w values)

**URL:** https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py
**Accessed:** 2026-06-13 (fetched via `curl` from the `biopython-179` tag; the file was removed from `master`)
**Authority rank:** 3 (reference implementation, Biopython v1.79)

**Key Extracted Points:**

1. **Provenance (file header, verbatim):** "Codon adaption indxes, including Sharp and Li (1987) E. coli index ... from Sharp & Li, Nucleic Acids Res. 1987."
2. **Sample w values:** `GCT=1, GCC=0.122, GCA=0.586, GCG=0.424` (Ala); `CTG=1, CTC=0.037, CTA=0.007, CTT=0.042, TTA=0.02, TTG=0.02` (Leu); `CGT=1, CGC=0.356, CGA=0.004, CGG=0.004, AGA=0.004, AGG=0.002` (Arg); `TTT=0.296, TTC=1` (Phe); `ATG=1`, `TGG=1`. (Full 64-codon table transcribed into `EColiOptimalCodons`; stop codons added as 0.0.)

### Kazusa Codon Usage Database — *Homo sapiens* [gbpri]

**URL:** https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 5 (curated database; underlying method Nakamura et al. 2000)

**Key Extracted Points:**

1. **Dataset scale:** 93,487 coding sequences, 40,662,582 codons.
2. **Sample per-thousand frequencies:** CTG 39.6, GTG 28.1, GCC 27.7, AAG 31.9, CAG 34.2, GAG 39.6, ATG 22.0, TGG 13.2. (Full 64-codon per-thousand table used to derive `HumanOptimalCodons` RSCU; see Dataset below.)

### CodonW GC3s + ENc plot reference (Wright 1990 relationship, viral codon study)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7596632/
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 4 (peer-reviewed paper applying CodonW; confirms exclusion set)

**Key Extracted Points:**

1. **59-codon synonymous set (verbatim):** RSCU/GC3s are computed over "the 59 synonymous codons (excluding the codons of AUG, UGG and the three stop codons)" — confirms the Met/Trp/stop exclusion in Peden's GC3s definition.

---

## Documented Corner Cases and Failure Modes

### From seqinr / CodonW

1. **Single-codon amino acids and stop codons:** excluded from CAI (and from GC3s). A sequence containing only Met/Trp/stop codons has no scorable codons → CAI is undefined (this implementation returns 0).
2. **Zero-frequency codons:** seqinr/EMBOSS substitute a small value (0.01, Bulmer 1988) to avoid `ln(0)`. This implementation instead skips codons whose relative adaptiveness is 0 (so an entirely-zero gene yields CAI 0). Documented as a deviation in the algorithm doc.

### From Peden thesis

1. **GC3s denominator empty:** a sequence with no synonymous codons (only Met/Trp/stop) has an empty GC3s denominator → GC3s reported as 0.

---

## Test Datasets

### Dataset: E. coli relative adaptiveness (Sharp & Li 1987 via Biopython SharpEcoliIndex)

**Source:** Sharp & Li (1987) NAR 15:1281–1295, Table 1; reproduced in Biopython v1.79 `SharpEcoliIndex`.

| Codon (AA) | w |
|-----------|----|
| CTG (Leu) | 1.000 |
| CTC (Leu) | 0.037 |
| CTA (Leu) | 0.007 |
| GCT (Ala) | 1.000 |
| GCC (Ala) | 0.122 |
| GCA (Ala) | 0.586 |
| CGT (Arg) | 1.000 |
| AGG (Arg) | 0.002 |
| TTC (Phe) | 1.000 |
| TTT (Phe) | 0.296 |

### Dataset: Human RSCU (derived from Kazusa H. sapiens per-thousand frequencies)

**Source:** Kazusa Codon Usage Database, Homo sapiens [gbpri] (accessed 2026-06-13). RSCU_j = n·x_j / Σ_k x_k.

| Codon (AA) | per-1000 | RSCU |
|-----------|----------|------|
| CTG (Leu, n=6) | 39.6 | 2.3713 |
| GCC (Ala, n=4) | 27.7 | 1.5988 |
| ATG (Met, n=1) | 22.0 | 1.0000 |
| TGG (Trp, n=1) | 13.2 | 1.0000 |

### Dataset: Worked CAI / GC examples (derived from the formulas above)

| Input (DNA) | Quantity | Expected | Derivation |
|-------------|----------|----------|------------|
| `CTGATCGTTGCTCGTAAA` | CAI vs E. coli w | 1.0 | all six codons have w=1 → geo mean 1 |
| `GCTGCC` | CAI vs E. coli w | 0.34928498393146 | √(1×0.122) (Ala GCT w=1, GCC w=0.122) |
| `CTAATAGTC` | CAI vs E. coli w | 0.01114947479545 | ∛(0.007×0.003×0.066) |
| `ATGTGGTAA` | CAI vs E. coli w | 0 | Met+Trp+stop only → no scorable codon |
| `GCCGCA` | GC3s | 50.0 | 2 Ala codons; 3rd = C(GC), A(not) → 1/2 |
| `ATGGCA` | GC3s | 0.0 | Met excluded; Ala GCA 3rd = A → 0/1 |
| `ATGGCA` | GC3 (all) | 50.0 | positions 3 = G(ATG), A → 1/2 (shows GC3 ≠ GC3s) |
| `CTGGTTAAA` | GC1 / GC2 / GC3 | 66.666… / 0.0 / 33.333… | per-position counts over 3 codons |

---

## Assumptions

1. **ASSUMPTION: GC3s reported as a percentage.** CodonW reports GC3s as a fraction in [0,1]; this implementation reports it as a percentage (×100) for consistency with the existing GC1/GC2/GC3 fields (which follow EMBOSS cusp percentage style). Non-correctness-affecting unit/labeling choice; documented in the algorithm doc. The synonymous-codon *subset* used in the numerator/denominator is exactly per Peden.
2. **ASSUMPTION: zero-w codons are skipped rather than floored to 0.01.** Sharp & Li / Bulmer floor missing codons to 0.01; this implementation skips codons whose relative adaptiveness is 0. For the supplied reference tables no synonymous codon has w=0, so CAI on real CDS is unaffected; only a gene using a codon entirely absent from the reference differs. Documented as a deviation.

---

## Recommendations for Test Coverage

1. **MUST Test:** CAI of an all-optimal E. coli sequence equals 1.0 — Evidence: Sharp & Li 1987 (w=1 for the preferred synonym).
2. **MUST Test:** CAI equals the geometric mean of w for a small hand-derived sequence (`GCTGCC` → √0.122) — Evidence: Sharp & Li 1987 formula.
3. **MUST Test:** CAI excludes Met, Trp and stop codons (`ATGTGGTAA` → 0) — Evidence: seqinr / CodonW exclusion rule.
4. **MUST Test:** GC3s counts only synonymous third positions, excluding Met/Trp/stop (`ATGGCA` GC3s=0 while GC3=50) — Evidence: Peden §1.8.2.1.3.
5. **MUST Test:** GC1/GC2/GC3 per-position percentages (`CTGGTTAAA`) — Evidence: EMBOSS cusp position GC.
6. **MUST Test:** `EColiOptimalCodons` reproduces Sharp & Li 1987 w values (CTG=1, GCC=0.122, …) — Evidence: Biopython SharpEcoliIndex.
7. **MUST Test:** `HumanOptimalCodons` reproduces Kazusa-derived RSCU (CTG≈2.3713, GCC≈1.5988, ATG=1) — Evidence: Kazusa H. sapiens.
8. **SHOULD Test:** null `DnaSequence` / null reference throw `ArgumentNullException`; empty string returns zeroed statistics / CAI 0 — Rationale: documented input contract.
9. **COULD Test:** `OverallGc` equals (GC1+GC2+GC3)/3 — Rationale: derived convenience property.

---

## References

1. Sharp PM, Li W-H. (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15(3):1281–1295. https://doi.org/10.1093/nar/15.3.1281
2. Peden JF. (1999). Analysis of Codon Usage (Ph.D. thesis, University of Nottingham; CodonW reference). https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf
3. CodonW — Codon usage indices. https://codonw.sourceforge.net/Indices.html
4. Charif D, Lobry JR. seqinr `cai` function documentation. https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html
5. EMBOSS `cusp` application documentation. https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp
6. Biopython v1.79, Bio.SeqUtils.CodonUsageIndices (`SharpEcoliIndex`). https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py
7. Nakamura Y, Gojobori T, Ikemura T. (2000). Codon usage tabulated from international DNA sequence databases. Nucleic Acids Research 28(1):292. Kazusa Codon Usage Database, Homo sapiens [gbpri]. https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
8. Sharp PM, Tuohy TMF, Mosurski KR. (1986). Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Research 14(13):5125–5143. https://doi.org/10.1093/nar/14.13.5125
9. Wright F. (1990). The 'effective number of codons' used in a gene. Gene 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9

---

## Change History

- **2026-06-13**: Initial documentation (CODON-STATS-001).
