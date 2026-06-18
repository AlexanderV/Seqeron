# Evidence Artifact: SEQ-DINUC-001

**Test Unit ID:** SEQ-DINUC-001
**Algorithm:** Dinucleotide Analysis (frequencies, observed/expected relative abundance, codon frequencies)
**Date Collected:** 2026-06-13

---

## Online Sources

### Karlin S. — "Pervasive properties of the genomic signature" (PMC126251)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC126251/
**Accessed:** 2026-06-13 (fetched via WebFetch of the PMC article page)
**Authority rank:** 1 (peer-reviewed; restates the Karlin & Burge 1995 genomic-signature model)

**Key Extracted Points:**

1. **Odds-ratio (relative abundance) definition:** Verbatim from the *Hypothesis formulation* section — "The assessment of bias in dinucleotide relative abundance begins with the 'odds ratios' r_xy = f_xy/f_x f_y where f_x denotes the (normalized) frequency of nucleotide (base) x and f_xy is the frequency of dinucleotide (base step) xy in the leading strand." Hence ρ_XY = f_XY / (f_X · f_Y).
2. **No-bias value:** r_xy = 1.0 corresponds to the dinucleotide frequency expected under statistical independence of neighbouring bases (no bias); r_xy > 1 over-represented, r_xy < 1 under-represented.
3. **Frequencies are normalized:** f_x and f_xy are normalized (relative) frequencies, i.e. counts divided by the respective totals, not raw counts.

### Tsirigos & Rigoutsos / Karlin & Burge thresholds (MBE 19(6):964, "Relative Abundance of Dinucleotides in Transposable Elements")

**URL:** https://academic.oup.com/mbe/article/19/6/964/1095097
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed; applies and cites the Karlin & Burge 1995 criterion)

**Key Extracted Points:**

1. **Formula restated:** ρ*_XY = f*_XY / (f*_X f*_Y), with f*_XY the dinucleotide frequency and f*_X, f*_Y the individual base frequencies.
2. **Classification thresholds (verbatim):** "the XY dinucleotide was considered to be underrepresented if ρ*XY ≤ 0.78 and overrepresented if ρ*XY ≥ 1.23", attributed to Karlin and Burge (1995). (Used only for interpretation/documentation; the methods under test return the raw ratio, not a classification.)

### Gardiner-Garden M, Frommer M (1987) CpG O/E formula (search-derived; confirmed by retrieved restatements)

**URL:** https://link.springer.com/article/10.1007/BF00162972 (Gardiner-Garden & Frommer cited context) + WebSearch result text
**Accessed:** 2026-06-13 (WebSearch query "CpG observed expected ratio formula Gardiner-Garden Frommer 1987")
**Authority rank:** 1 / 4 (primary paper; formula restated by retrieved sources)

**Key Extracted Points:**

1. **CpG O/E (alternative normalization):** O/E = (#CpG/N) / ((#C/N)·(#G/N)), N = number of bases in the segment. This is the same odds-ratio shape as Karlin's ρ but normalizes the dinucleotide count by N rather than by (N−1). The repository follows the Karlin normalization (dinucleotide frequency over the (N−1) dinucleotide positions); the difference is the N/(N−1) factor and is recorded as a documented modeling choice, not an error.

### Kazusa Codon Usage Database — CUTG readme

**URL:** https://www.kazusa.or.jp/codon/readme_codon.html
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 5 (curated database) / 2 (de-facto standard tabulation)

**Key Extracted Points:**

1. **Codon frequency definition (verbatim):** "The table shows frequency (per thousand) and count for each codon as a sum of all CDS's of the organism." Each coding sequence is read as consecutive **non-overlapping** triplets within a reading frame; a codon's frequency is its count divided by the total number of counted codons.
2. **Ambiguous codons excluded:** codons containing ambiguous bases are excluded from the count. (The repository excludes any triplet not over {A,T,G,C}.)

---

## Documented Corner Cases and Failure Modes

### From Karlin (PMC126251)

1. **Statistical-independence baseline:** the expected frequency f_X·f_Y assumes positional independence; the ratio is undefined/uninformative when a constituent base is absent (f_X = 0 ⇒ expected = 0 ⇒ division by zero). The repository returns ratio 0 for such dinucleotides (expected = 0 guard).

### From Kazusa CUTG

1. **Reading frame:** counting starts at the frame offset and advances in steps of 3; trailing 1–2 leftover bases are not a codon and are ignored.
2. **Non-ACGT triplets:** excluded from counts.

---

## Test Datasets

### Dataset: Hand-derived odds-ratio example (sequence `ATGCGCGT`)

**Source:** Direct application of ρ_XY = f_XY/(f_X f_Y) [Karlin, PMC126251]. All values exact rationals.

| Parameter | Value |
|-----------|-------|
| Sequence | `ATGCGCGT` (length 8) |
| Mononucleotide counts | A=1, T=2, G=3, C=2 (total 8) |
| Dinucleotide counts (7 positions) | AT=1, TG=1, GC=2, CG=2, GT=1 |
| f_GC | 2/7 ; f_G=3/8, f_C=2/8 |
| ρ_GC = (2/7)/((3/8)(2/8)) | 64/21 = 3.047619047619048 |
| ρ_CG = (2/7)/((2/8)(3/8)) | 64/21 = 3.047619047619048 |
| ρ_AT = (1/7)/((1/8)(2/8)) | 32/7 = 4.571428571428571 |
| ρ_TG = (1/7)/((2/8)(3/8)) | 32/21 = 1.523809523809524 |
| ρ_GT = (1/7)/((3/8)(2/8)) | 32/21 = 1.523809523809524 |

### Dataset: Dinucleotide frequency example (same sequence)

**Source:** normalized frequency = count/(N−1) [Karlin: f_xy is a normalized frequency].

| Dinucleotide | Frequency |
|--------------|-----------|
| GC | 2/7 = 0.2857142857142857 |
| CG | 2/7 = 0.2857142857142857 |
| AT | 1/7 = 0.14285714285714285 |
| TG | 1/7 = 0.14285714285714285 |
| GT | 1/7 = 0.14285714285714285 |

### Dataset: Codon frequency example

**Source:** non-overlapping triplets per reading frame [Kazusa CUTG readme].

| Input | Frame | Codons counted | Expected frequencies |
|-------|-------|----------------|----------------------|
| `ATGATGAAA` | 0 | ATG, ATG, AAA (3) | ATG=2/3, AAA=1/3 |
| `ATGATGAAA` | 1 | TGA, TGA (2) | TGA=1.0 |

---

## Assumptions

1. **ASSUMPTION: Karlin (N−1) normalization vs Gardiner-Garden (N) normalization** — The repository computes dinucleotide frequency as count/(N−1) (Karlin odds-ratio convention, PMC126251) rather than count/N (Gardiner-Garden CpG convention). Both are published; the ratio is the same odds-ratio shape and differs only by the factor N/(N−1) → 1 for long sequences. This is a documented modeling choice, not an unresolved correctness gap, because the chosen convention is itself authoritative (Karlin).
2. **ASSUMPTION: U treated as a fifth base in `CalculateDinucleotideRatios`** — The single-base frequency denominator includes U for RNA inputs (matches `CalculateNucleotideComposition`). Not contradicted by any source; RNA dinucleotide signatures are defined the same way with U replacing T.

---

## Recommendations for Test Coverage

1. **MUST Test:** ρ_XY = f_XY/(f_X f_Y) for a multi-base DNA sequence with exact rational values (GC, CG, AT, TG, GT). — Evidence: Karlin PMC126251; MBE 19(6):964.
2. **MUST Test:** dinucleotide frequencies sum semantics and exact values for the same sequence. — Evidence: Karlin (normalized frequency).
3. **MUST Test:** codon frequencies for frame 0 and frame 1 with exact fractions; trailing bases ignored; non-ACGT excluded. — Evidence: Kazusa CUTG readme.
4. **MUST Test:** ρ for a perfectly uniform/independent sequence approximates 1.0 (no-bias baseline). — Evidence: Karlin (r=1 ⇒ no bias).
5. **SHOULD Test:** edge cases null/empty/length-<2 (ratios, freqs) and length-<3 (codons) return empty. — Rationale: documented input guards.
6. **SHOULD Test:** division-by-zero guard — dinucleotide whose constituent base is absent yields ratio 0. — Rationale: expected=0 guard.
7. **COULD Test:** lowercase and RNA (U) inputs handled (case-insensitive, U counted). — Rationale: implementation normalizes case and supports U.

---

## References

1. Karlin S. (1998). Global dinucleotide signatures and analysis of genomic heterogeneity / "Pervasive properties of the genomic signature". PMC. https://pmc.ncbi.nlm.nih.gov/articles/PMC126251/
2. Karlin S., Burge C. (1995). Dinucleotide relative abundance extremes: a genomic signature. Trends in Genetics 11(7):283-290. https://doi.org/10.1016/S0168-9525(00)89076-9 (criterion ρ≤0.78 / ρ≥1.23 retrieved from MBE 19(6):964: https://academic.oup.com/mbe/article/19/6/964/1095097)
3. Gardiner-Garden M., Frommer M. (1987). CpG islands in vertebrate genomes. J Mol Biol 196(2):261-282. https://doi.org/10.1016/0022-2836(87)90689-9 (CpG O/E formula restated via https://link.springer.com/article/10.1007/BF00162972)
4. Nakamura Y., Gojobori T., Ikemura T. — Codon Usage Database (CUTG), Kazusa. https://www.kazusa.or.jp/codon/readme_codon.html

---

## Change History

- **2026-06-13**: Initial documentation.
