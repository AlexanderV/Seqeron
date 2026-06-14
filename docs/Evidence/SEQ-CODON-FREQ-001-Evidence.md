# Evidence Artifact: SEQ-CODON-FREQ-001

**Test Unit ID:** SEQ-CODON-FREQ-001
**Algorithm:** Codon Frequencies (non-overlapping in-frame triplet usage)
**Date Collected:** 2026-06-14

---

## Online Sources

### Kazusa Codon Usage Database (CUTG) — Announcement / README

**URL:** https://www.kazusa.or.jp/codon/readme_codon.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5 (well-maintained bioinformatics database; the canonical codon-usage convention)

**Key Extracted Points:**

1. **Metric definition:** Retrieved text states "The table shows frequency (**per thousand**) and count for each codon as a sum of all CDS's of the organism." Codon usage is therefore a count of each codon aggregated over coding sequences, reported both as a raw count and as a per-thousand frequency.
2. **Aggregation method:** "the frequency (per thousand) of codon use in each organism was calculated by summing up the numbers of codons used" — i.e. frequency = (codon count / total codons) scaled to per-thousand.
3. **Ambiguous-base handling:** "Codons containing ambiguous base were excluded from count." Triplets with any non-ACGT base do not contribute to either the per-codon count or the total.

### EMBOSS `cusp` application documentation

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation, EMBOSS suite)

**Key Extracted Points:**

1. **Five output columns:** Codon, AA, Fraction, Frequency, Number.
2. **Fraction** = "the proportion of usage of the codon among its redundant set, i.e. the set of codons which code for this codon's amino acid" — a *per-amino-acid* proportion. This is a **different** metric from the one under test.
3. **Frequency** = "the expected number of codons, given the input sequence(s), per 1000 bases" (per-thousand, extrapolated for short inputs).
4. **Number** = the actual count of that codon observed in the input sequences.
5. **Worked dataset (verbatim sample output, all 64 codons):** the sum of the Number column is 386 codons. Selected rows: CGC Number=22, Frequency=56.995; GGC Number=23, Frequency=59.585; GCC Number=18, Frequency=46.632; TGA(*) Number=1, Frequency=2.591. The relationship Number/386 × 1000 = Frequency-per-thousand holds exactly (22/386×1000 = 56.995; 23/386×1000 = 59.585), confirming that the count/total-codons fraction under test equals Kazusa/cusp per-thousand frequency ÷ 1000.

### Wikipedia — Codon usage bias (for cited primaries only)

**URL:** https://en.wikipedia.org/wiki/Codon_usage_bias
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4

**Key Extracted Points:**

1. **Definition:** "Codon usage bias refers to differences in the frequency of occurrence of synonymous codons in coding DNA." Establishes that codon frequency is the unit of measurement.
2. **Primaries cited:** Sharp & Li (1987, CAI), Ikemura (1981, frequency of optimal codons). These are downstream indices built on codon counts; not needed for the raw frequency definition, which is covered by Kazusa.

### Nakamura, Gojobori, Ikemura (2000) — citation verification

**URL:** https://pubmed.ncbi.nlm.nih.gov/10592250/ (fetched via WebFetch)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, the paper underlying the Kazusa CUTG database)

**Key Extracted Points:**

1. **Full citation confirmed:** Nakamura Y, Gojobori T, Ikemura T (2000). "Codon usage tabulated from international DNA sequence databases: status for the year 2000." Nucleic Acids Research 28(1):292. DOI: 10.1093/nar/28.1.292.
2. The PMC full-text abstract (https://pmc.ncbi.nlm.nih.gov/articles/PMC102460/, fetched) confirms the database tabulates codon usage over complete CDSs; the per-codon counting/per-thousand convention is given in the Kazusa README above.

---

## Documented Corner Cases and Failure Modes

### From Kazusa CUTG README

1. **Ambiguous / non-ACGT codon:** excluded from the count entirely (does not affect counts or total).

### From the count/total definition (CUTG aggregation)

1. **Trailing partial codon (length not a multiple of 3 from the frame start):** only complete non-overlapping triplets are counted; the 1–2 leftover bases are not a codon and are ignored.
2. **Reading frame offset:** counting starts at the chosen frame index; the same sequence yields a different codon multiset per frame.
3. **No valid codon (all triplets ambiguous, or sequence shorter than 3):** total = 0; an empty frequency table is the only well-defined result (no division by zero).

---

## Test Datasets

### Dataset: EMBOSS cusp sample output (count/total cross-check)

**Source:** EMBOSS `cusp` documentation, https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html (accessed 2026-06-14)

| Parameter | Value |
|-----------|-------|
| Total codons (Σ Number column) | 386 |
| CGC count | 22 |
| CGC per-thousand frequency | 56.995 |
| CGC count/total fraction (metric under test) | 22/386 = 0.0569948… = 56.995 / 1000 |
| GGC count | 23 |
| GGC count/total fraction | 23/386 = 0.0595854… = 59.585 / 1000 |

### Dataset: Hand-derived small sequences (exact rationals)

**Source:** Direct application of the CUTG count/total definition.

| Input | Frame | Codons read | Expected frequencies |
|-------|-------|-------------|----------------------|
| `ATGATGAAA` | 0 | ATG, ATG, AAA | ATG = 2/3, AAA = 1/3 |
| `ATGATGAAA` | 1 | TGA, TGA (AA leftover) | TGA = 1.0 |
| `ATGNNNAAA` | 0 | ATG, NNN(excluded), AAA | ATG = 1/2, AAA = 1/2 |
| `ATGAA` | 0 | ATG (AA leftover) | ATG = 1.0 |
| `atgaaa` | 0 | ATG, AAA (lowercase normalized) | ATG = 1/2, AAA = 1/2 |

---

## Assumptions

1. **ASSUMPTION: empty result for zero valid codons.** Kazusa specifies ambiguous codons are excluded and that frequency is count/total, but does not explicitly define the result when *no* valid codon exists (total = 0). Returning an empty frequency table is the only value consistent with "count / total" (a frequency over zero codons is undefined; reporting none is the conservative, non-fabricating choice). This is non-correctness-affecting for any input containing at least one valid codon and matches the implementation's guard. Recorded as the single open assumption.

---

## Recommendations for Test Coverage

1. **MUST Test:** Frame-0 exact frequencies on a known multiset (`ATGATGAAA` → ATG=2/3, AAA=1/3) — Evidence: Kazusa CUTG count/total.
2. **MUST Test:** Reading-frame offset changes the codon multiset (`ATGATGAAA` frame 1 → TGA=1.0) — Evidence: CUTG non-overlapping triplets from frame start.
3. **MUST Test:** Frequencies over counted codons sum to 1.0 (INV-02) — Evidence: count/total normalization.
4. **MUST Test:** Non-ACGT triplet excluded (`ATGNNNAAA` → ATG=1/2, AAA=1/2) — Evidence: CUTG "ambiguous codons excluded".
5. **MUST Test:** Cross-check against the EMBOSS cusp dataset relationship (count/total = per-thousand ÷ 1000) on a constructed multiset reproducing a cusp ratio — Evidence: cusp sample output.
6. **SHOULD Test:** Trailing 1–2 bases ignored (`ATGAA` → ATG=1.0) — Rationale: documented non-overlapping-triplet remainder rule.
7. **SHOULD Test:** Case-insensitive normalization (`atgaaa` equals `ATGAAA`) — Rationale: implementation upper-cases input; codons are case-independent.
8. **SHOULD Test:** All-ambiguous sequence yields empty table (total = 0, no division by zero) — Rationale: documented zero-codon corner case.
9. **COULD Test:** null / empty / length < 3 return empty — Rationale: standard guard contract shared with sibling methods.

---

## References

1. Nakamura Y, Gojobori T, Ikemura T (2000). Codon usage tabulated from international DNA sequence databases: status for the year 2000. Nucleic Acids Research 28(1):292. https://doi.org/10.1093/nar/28.1.292
2. Kazusa DNA Research Institute. Codon Usage Database (CUTG) — Announcement / README. https://www.kazusa.or.jp/codon/readme_codon.html (accessed 2026-06-14)
3. Rice P, Longden I, Bleasby A (2000). EMBOSS: The European Molecular Biology Open Software Suite — `cusp` application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html (accessed 2026-06-14)
4. Wikipedia. Codon usage bias. https://en.wikipedia.org/wiki/Codon_usage_bias (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
