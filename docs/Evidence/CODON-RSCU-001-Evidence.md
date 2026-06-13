# Evidence Artifact: CODON-RSCU-001

**Test Unit ID:** CODON-RSCU-001
**Algorithm:** Relative Synonymous Codon Usage (RSCU) and codon counting
**Date Collected:** 2026-06-13

---

## Online Sources

### Sharp, Tuohy & Mosurski (1986) — "Codon usage in yeast: cluster analysis…" (NAR 14(13):5125-5143)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC311530/ (PMID 3526280; DOI 10.1093/nar/14.13.5125)
**Accessed:** 2026-06-13 (WebSearch "Sharp Tuohy Mosurski 1986 Codon usage in yeast cluster analysis RSCU formula"; WebFetch of the PMC article page)
**Authority rank:** 1 (peer-reviewed primary paper that introduced RSCU)

**Key Extracted Points:**

1. **Origin of RSCU:** This paper introduces the Relative Synonymous Codon Usage (RSCU) measure used to quantify synonymous codon bias across yeast genes (abstract/title retrieved from the PMC page and Oxford Academic record). Full citation confirmed: *Nucleic Acids Res.* 14(13):5125-5143.
2. **Retrieval note:** the PMC page renders the methods pages as images, so the verbatim equation could not be transcribed from this page; the explicit formula was obtained from the reference implementations below, which cite this same paper.

### Suzuki et al. / LIRMM "RSCU RS" page — explicit formula

**URL:** https://www.lirmm.fr/~rivals/rscu/
**Accessed:** 2026-06-13 (WebSearch + WebFetch)
**Authority rank:** 3 (reference implementation / academic tool documentation)

**Key Extracted Points:**

1. **Explicit formula (verbatim):** "RSCU_i = X_i / ((1/N_i) ∑ X_j), where X_i is the number of occurrences … of codon i." The sum runs over the N_i synonymous codons of the family codon i belongs to (its degeneracy).
2. **Interpretation:** RSCU values indicate whether a codon is used more (>1) or less (<1) than expected under uniform usage within its synonymous group; the no-bias value is 1.

### GenomicSig (CRAN) `RSCU` reference function — explicit indexed formula

**URL:** https://rdrr.io/cran/GenomicSig/man/RSCU.html
**Accessed:** 2026-06-13 (WebSearch + WebFetch)
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Explicit indexed formula (verbatim):** "For the j-th codon of amino acid i, let x_{i,j} denote the number of occurrences … RSCU_{i,j} = (n_i × x_{i,j}) / (Σ_{j=1}^{n_i} x_{i,j})", where n_i is the number of codons that code for amino acid i and x_{i,j} is the count of codon j. This is algebraically identical to the LIRMM form.
2. **Bounds:** RSCU is "comprised between 0 and the number of synonymous codons for that amino acid" → range [0, n_i]; max n_i reached when only one synonymous codon is used.

### seqinr `uco` reference function (CRAN) — definition and references

**URL:** https://search.r-project.org/CRAN/refmans/seqinr/html/uco.html and https://www.rdocumentation.org/packages/seqinr/versions/4.2-8/topics/uco
**Accessed:** 2026-06-13 (WebSearch + WebFetch)
**Authority rank:** 3 (Biopython-equivalent reference implementation in R)

**Key Extracted Points:**

1. **Definition (verbatim):** "RSCU values are the number of times a particular codon is observed, relative to the number of times that the codon would be observed for a uniform synonymous codon usage (i.e. all the codons for a given amino-acid have the same probability)."
2. **No-bias / interpretation:** "In the absence of any codon usage bias, the RSCU values would be 1.00"; values <1 underused, >1 overused.
3. **Primary reference cited:** Sharp, P.M., Tuohy, T.M.F., Mosurski, K.R. (1986) — confirming the formula traces to the 1986 paper.

### cubar `est_rscu` reference function (CRAN) — zero-count handling

**URL:** https://rdrr.io/cran/cubar/man/est_rscu.html
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 3 (reference implementation documentation)

**Key Extracted Points:**

1. **Pseudocount (verbatim):** `pseudo_cnt` is a "Numeric pseudo count added to avoid division by zero when few sequences are available for RSCU calculation (default: 1)." This shows the 0/0 case (an entirely absent synonymous family) is implementation-defined; classic Sharp et al. RSCU has no pseudocount, so the repository's choice (return 0 for absent families) is a documented convention, not a divergence from the canonical formula for present families.
2. **Stop codons:** `incl_stop` default FALSE (stop codons excluded by default); the repository instead treats the three standard stop codons as a synonymous family of size 3, which does not affect RSCU of any amino-acid codon.

### Begomovirus codon-usage paper (PMC2528880) — independent restatement of definition

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/
**Accessed:** 2026-06-13 (WebFetch)
**Authority rank:** 1 (peer-reviewed) / restatement

**Key Extracted Points:**

1. **Definition (verbatim):** "Relative synonymous codon usage (RSCU) is defined as the ratio of the observed frequency of codons to the expected frequency given that all the synonymous codons for the same amino acids are used equally." Attributed to Sharp & Li 1986.

---

## Documented Corner Cases and Failure Modes

### From GenomicSig / LIRMM

1. **Bounds:** RSCU ∈ [0, n_i]. A codon never used in its family → 0; a codon that is the only one used → n_i.
2. **Family-sum invariant:** because RSCU_{i,j} = n_i·x_{i,j}/Σx, the RSCU values within one present family sum to n_i (Σ_j n_i·x_{i,j}/Σx = n_i).

### From cubar / seqinr

1. **Absent family (0/0):** when no codon of a family appears, the denominator is 0; canonical RSCU is undefined. cubar avoids this with a pseudocount; the repository returns 0 for every codon of an absent family.
2. **Single-codon families (Met=ATG, Trp=TGG):** n_i = 1, so RSCU = 1 whenever present (no bias possible).

### From repository contract (CountCodons)

1. **Reading frame / triplets:** codons are non-overlapping triplets from offset 0; trailing 1–2 bases are ignored.
2. **Non-ACGT triplets excluded:** any triplet containing a character outside {A,C,G,T} is not counted (string overload uppercases first).

---

## Test Datasets

### Dataset: Leu (6-fold) worked example — sequence `CTGCTGCTGCTA`

**Source:** direct application of RSCU_{i,j} = (n_i × x_{i,j}) / Σ x [GenomicSig; LIRMM; Sharp et al. 1986]. Leu has n_i = 6 codons {TTA,TTG,CTT,CTC,CTA,CTG}.

| Parameter | Value |
|-----------|-------|
| Sequence | `CTGCTGCTGCTA` (4 Leu codons) |
| Counts | CTG=3, CTA=1, others=0; family total = 4 |
| RSCU(CTG) = 6·3/4 | 4.5 |
| RSCU(CTA) = 6·1/4 | 1.5 |
| RSCU(TTA)=RSCU(TTG)=RSCU(CTT)=RSCU(CTC) | 0.0 |
| Σ over the 6 Leu codons | 6.0 (= n_i) |

### Dataset: Phe (2-fold) worked example — sequence `TTTTTTTTC`

**Source:** RSCU_{i,j} = (n_i × x_{i,j}) / Σ x. Phe has n_i = 2 {TTT,TTC}.

| Parameter | Value |
|-----------|-------|
| Sequence | `TTTTTTTTC` (TTT×2, TTC×1) |
| RSCU(TTT) = 2·2/3 | 1.3333333333333333 (4/3) |
| RSCU(TTC) = 2·1/3 | 0.6666666666666666 (2/3) |
| Σ over Phe | 2.0 (= n_i) |

### Dataset: Unbiased 2-fold — sequence `TTTTTC`

**Source:** equal usage ⇒ no bias ⇒ RSCU = 1 [seqinr; LIRMM "no-bias value is 1"].

| Parameter | Value |
|-----------|-------|
| RSCU(TTT) = 2·1/2 | 1.0 |
| RSCU(TTC) = 2·1/2 | 1.0 |

### Dataset: Single-codon family — sequence `ATGATG`

**Source:** n_i = 1 for Met ⇒ RSCU = 1 [GenomicSig bounds; cubar single-codon note].

| Parameter | Value |
|-----------|-------|
| RSCU(ATG) | 1.0 |

### Dataset: CountCodons frame/exclusion — `ATGAAATGA`, `ATGAA`, `ATGNNNAAA`

**Source:** non-overlapping triplets; trailing bases ignored; non-ACGT excluded (repository contract; consistent with Kazusa CUTG codon-counting convention).

| Input | Codons counted | Notes |
|-------|----------------|-------|
| `ATGAAATGA` | ATG=1, AAA=1, TGA=1 | 3 full triplets |
| `ATGAA` | ATG=1 | trailing `AA` ignored |
| `ATGNNNAAA` (string overload) | ATG=1, AAA=1 | `NNN` triplet excluded |

---

## Assumptions

1. **ASSUMPTION: Absent-family 0/0 handling returns 0** — The canonical formula is undefined when a synonymous family has zero observed codons (0/0). cubar resolves this with a pseudocount (default 1); the repository instead returns 0 for every codon of an absent family. This only affects families that do not occur in the input (where no canonical value exists), never the RSCU of an observed codon. Documented as a convention, not a correctness gap for present families.
2. **ASSUMPTION: Stop codons treated as a 3-fold family** — The repository groups TAA/TAG/TGA as one synonymous family (degeneracy 3) and computes RSCU for them like any family. Reference tools commonly exclude stop codons; this choice does not change RSCU for any amino-acid codon.

---

## Recommendations for Test Coverage

1. **MUST Test:** RSCU for a 6-fold Leu family (`CTGCTGCTGCTA`) with exact values 4.5, 1.5, 0.0 and family sum = 6. — Evidence: GenomicSig formula; LIRMM; Sharp et al. 1986.
2. **MUST Test:** RSCU for a 2-fold Phe family (`TTTTTTTTC`) = 4/3 and 2/3; sum = 2. — Evidence: GenomicSig formula.
3. **MUST Test:** unbiased family (`TTTTTC`) ⇒ RSCU = 1.0 for both codons. — Evidence: seqinr/LIRMM no-bias value.
4. **MUST Test:** single-codon family (`ATGATG`) ⇒ RSCU(ATG) = 1.0. — Evidence: GenomicSig bounds; cubar.
5. **MUST Test (CountCodons):** non-overlapping triplets, repeated codons, trailing-base truncation, non-ACGT exclusion (string overload). — Evidence: repository contract / Kazusa convention.
6. **SHOULD Test:** null DnaSequence throws ArgumentNullException; empty string returns empty dictionary (both methods). — Rationale: documented input guards.
7. **SHOULD Test (invariant):** within any present family, RSCU values sum to n_i. — Evidence: derivation from formula (GenomicSig).
8. **COULD Test:** lowercase input handled by string overload (case-insensitive). — Rationale: implementation uppercases.

---

## References

1. Sharp P.M., Tuohy T.M.F., Mosurski K.R. (1986). Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. *Nucleic Acids Research* 14(13):5125-5143. https://doi.org/10.1093/nar/14.13.5125 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC311530/)
2. Suzuki H. et al. / LIRMM. RSCU RS — Measuring the bias in codon usage. https://www.lirmm.fr/~rivals/rscu/
3. GenomicSig (CRAN). RSCU: Relative Synonymous Codon Usage. https://rdrr.io/cran/GenomicSig/man/RSCU.html
4. Charif D., Lobry J.R. seqinr — `uco`: Codon usage indices. https://search.r-project.org/CRAN/refmans/seqinr/html/uco.html
5. cubar (CRAN). `est_rscu`: Estimate Relative Synonymous Codon Usage. https://rdrr.io/cran/cubar/man/est_rscu.html
6. (Restatement) Analysis of synonymous codon usage and evolution of begomoviruses. PMC2528880. https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/

---

## Change History

- **2026-06-13**: Initial documentation.
