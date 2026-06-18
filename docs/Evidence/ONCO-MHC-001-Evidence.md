# Evidence Artifact: ONCO-MHC-001

**Test Unit ID:** ONCO-MHC-001
**Algorithm:** MHC-Peptide Binding Classification (length filtering + affinity/%rank thresholds)
**Date Collected:** 2026-06-14

---

## Online Sources

### NetMHCpan-4.1 / NetMHCIIpan-4.0 (Reynisson et al. 2020)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7319546/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Nucleic Acids Research*)
**How retrieved:** WebSearch query `NetMHCpan-4.1 Reynisson 2020 strong binder weak binder %rank threshold 0.5 2 IC50` → opened the PMC full text with WebFetch.

**Key Extracted Points:**

1. **Class I %Rank thresholds (verbatim):** "by default, %Rank < 0.5% and %Rank < 2% thresholds are considered for detecting SBs and WBs for class I" (SB = strong binder, WB = weak binder).
2. **Class II %Rank thresholds (verbatim):** "and %Rank < 2% and %Rank < 10%, for SBs and WBs for class II".
3. **Class I peptide length range (verbatim):** "for class I, the length range goes from 8 to 14 amino acids, default is 8–11".
4. **Citation extracted from the page:** Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. *Nucleic Acids Res*. 2020;48(W1):W449–W454. doi:10.1093/nar/gkaa379.

---

### Sette et al. (1994) — class I affinity vs immunogenicity

**URL:** https://pubmed.ncbi.nlm.nih.gov/7527444/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Journal of Immunology*)
**How retrieved:** WebSearch query `Sette 1994 Journal of Immunology relationship affinity MHC class I immunogenicity 500 nM CTL response threshold` → opened the PubMed abstract page with WebFetch.

**Key Extracted Points:**

1. **Affinity threshold (verbatim from abstract):** "an affinity threshold of approximately 500 nM (preferably 50 nM or less) apparently determines the capacity" (of a peptide epitope to elicit a CTL response). This is the primary biological basis for the 50 nM (strong) and 500 nM (binder) IC50 cutoffs.
2. **Citation:** Sette A, Vitiello A, Reherman B, Fowler P, Nayersina R, Kast WM, Melief CJ, Oseroff C, Yuan L, Ruppert J, Sidney J, del Guercio MF, Southwood S, Kubo RT, Chesnut RW, Grey HM, Chisari FV. The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *J Immunol*. 1994 Dec 15;153(12):5586–92. PMID: 7527444.

---

### IEDB — selecting thresholds for MHC binding predictions

**URL:** https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated bioinformatics resource; IEDB)
**How retrieved:** WebSearch query `IEDB MHC class I binding prediction strong binder IC50 < 50 nM intermediate < 500 nM threshold` returned this IEDB help article; the article body is served behind a 403 to WebFetch, so the verbatim statement below is the text returned in the IEDB search-result snippet (the direct page fetch returned HTTP 403 — recorded here for auditability).

**Key Extracted Points:**

1. **IC50 affinity tiers (verbatim from snippet):** "Peptides with IC50 values <50 nM are considered high affinity, <500 nM intermediate affinity and <5000 nM low affinity."
2. **Binder cutoff (verbatim from snippet):** "an absolute binding affinity (IC50) threshold of 500 nM identifies strong binders".

---

### Roomp, Antes & Lengauer (2010) — corroborating the 500 nM binder cutoff

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2836306/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *BMC Bioinformatics*)
**How retrieved:** WebSearch (IEDB threshold query, above) listed this PMC article; opened it with WebFetch (after following the ncbi→pmc redirect) to corroborate the IC50 = 500 nM binder demarcation independently of IEDB.

**Key Extracted Points:**

1. **500 nM binder demarcation (verbatim):** "Any peptides annotated in IEDB as binders with IC50 values greater than 500 nM, and peptides annotated as non-binders with IC50 values less than 500 nM were discarded." — confirms 500 nM as the binder/non-binder demarcation in a peer-reviewed source.
2. **Citation:** Roomp K, Antes I, Lengauer T. Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics*. 2010 Feb 17;11:90. doi:10.1186/1471-2105-11-90.

---

### IEDB — MHC class II binding tool description (peptide length range)

**URL:** https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated bioinformatics resource; IEDB)
**How retrieved:** WebSearch query `IEDB MHC class II peptide length 13-25 binding prediction core 9-mer length range` returned this IEDB article; page body is 403 to WebFetch, so the verbatim statement below is from the search-result snippet (direct fetch returned HTTP 403 — recorded for auditability).

**Key Extracted Points:**

1. **Class II length range (verbatim from snippet):** "Peptides binding to MHC class II molecules can vary considerably, and typically range between 13 and 25 amino acids long."
2. **Binding core (verbatim from snippet):** "The peptide binding core ... is usually nine amino acids long (9-mer)."

---

## Documented Corner Cases and Failure Modes

### From Reynisson et al. (2020)

1. **Class I length out of range:** lengths below 8 or above 14 are not valid class I peptides; the canonical neoantigen search default is 8–11. A peptide whose length is outside the class's accepted range is not a valid binder candidate.

### From Sette et al. (1994) / IEDB

1. **Boundary semantics:** the cutoffs are stated as strict "<" inequalities (`<50 nM`, `<500 nM`). A peptide at exactly 50 nM is therefore NOT strong by the strict tier text; a peptide at exactly 500 nM is NOT a binder by the strict tier text. The tests encode the strict-inequality reading from the source.

---

## Test Datasets

### Dataset: NetMHCpan-4.1 class I %Rank decision points

**Source:** Reynisson et al. (2020), PMC7319546 — "%Rank < 0.5% and %Rank < 2% ... for SBs and WBs for class I".

| Parameter | Value |
|-----------|-------|
| Strong binder %Rank cutoff (class I) | < 0.5 |
| Weak binder %Rank cutoff (class I) | < 2.0 |
| %Rank = 0.5 (boundary) | NOT strong (strict `<`) → Weak (0.5 < 2.0) |
| %Rank = 2.0 (boundary) | NOT weak (strict `<`) → NonBinder |
| %Rank = 0.4 | Strong |
| %Rank = 1.0 | Weak |
| %Rank = 5.0 | NonBinder |

### Dataset: IEDB / Sette class I IC50 (nM) decision points

**Source:** Sette et al. (1994), PMID 7527444; IEDB threshold tiers (high <50, intermediate <500).

| Parameter | Value |
|-----------|-------|
| Strong binder IC50 cutoff | < 50 nM |
| Weak (intermediate) binder IC50 cutoff | < 500 nM |
| IC50 = 10 nM | Strong |
| IC50 = 50 nM (boundary) | NOT strong → Weak (50 < 500) |
| IC50 = 200 nM | Weak |
| IC50 = 500 nM (boundary) | NOT weak → NonBinder |
| IC50 = 1000 nM | NonBinder |

### Dataset: peptide length validity by MHC class

**Source:** Reynisson et al. (2020) (class I 8–14, default 8–11); IEDB class II tool description (13–25).

| Parameter | Value |
|-----------|-------|
| Class I accepted length range | 8–11 (canonical neoantigen search default) |
| Class II accepted length range | 13–25 |
| Class I length 7 | invalid (too short) |
| Class I length 9 | valid |
| Class I length 12 | invalid (above canonical default range) |
| Class II length 15 | valid |
| Class II length 12 | invalid |

---

## Assumptions

1. **ASSUMPTION: Class I canonical length range = 8–11.** The source gives 8–14 as the full class I range with 8–11 as the default. This unit adopts 8–11 as the accepted class I range to match the existing `OncologyAnalyzer.MhcClassIMin/MaxPeptideLength` constants (ONCO-NEO-001) and the pVACtools canonical neoantigen search. This is a documented default, not an invented value; callers can pass an explicit range. It affects `IsValidPeptideLength` output for lengths 12–14.

---

## Recommendations for Test Coverage

1. **MUST Test:** IC50 classification at and around both cutoffs (10/50/200/500/1000 nM) yields Strong/Weak/NonBinder per the strict-`<` tiers — Evidence: Sette 1994; IEDB tiers.
2. **MUST Test:** %Rank classification at and around both cutoffs (0.4/0.5/1.0/2.0/5.0) yields Strong/Weak/NonBinder for class I — Evidence: Reynisson 2020.
3. **MUST Test:** class II %Rank cutoffs (2.0 / 10.0) classify SB/WB — Evidence: Reynisson 2020.
4. **MUST Test:** peptide length validity for class I (8–11) and class II (13–25) at boundaries — Evidence: Reynisson 2020; IEDB class II tool description.
5. **SHOULD Test:** invalid inputs rejected — negative/NaN/infinite IC50, negative/NaN %Rank, %Rank > 100 — Rationale: %Rank is a percentile in [0,100]; IC50 is a positive concentration (invariant IC50 > 0).
6. **COULD Test:** a combined "classify candidate" helper that gates on length validity then affinity — Rationale: mirrors how a caller wires ONCO-NEO-001 windows to a supplied affinity.

---

## References

1. Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. *Nucleic Acids Research* 48(W1):W449–W454. https://doi.org/10.1093/nar/gkaa379 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC7319546/)
2. Sette A, Vitiello A, Reherman B, et al. (1994). The relationship between class I binding affinity and immunogenicity of potential cytotoxic T cell epitopes. *Journal of Immunology* 153(12):5586–5592. https://pubmed.ncbi.nlm.nih.gov/7527444/
3. Roomp K, Antes I, Lengauer T (2010). Predicting MHC class I epitopes in large datasets. *BMC Bioinformatics* 11:90. https://doi.org/10.1186/1471-2105-11-90 (PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC2836306/)
4. IEDB. What thresholds (cut-offs) should I use for MHC class I and II binding predictions. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094152371-What-thresholds-cut-offs-should-I-use-for-MHC-class-I-and-II-binding-predictions
5. IEDB. T Cell Epitopes - MHC Class II Binding Prediction Tools Description. Accessed 2026-06-14. https://help.iedb.org/hc/en-us/articles/114094151731-T-Cell-Epitopes-MHC-Class-II-Binding-Prediction-Tools-Description

---

## Change History

- **2026-06-14**: Initial documentation.
