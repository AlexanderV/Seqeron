# Evidence Artifact: ANNOT-CODONUSAGE-001

**Test Unit ID:** ANNOT-CODONUSAGE-001
**Algorithm:** Relative Synonymous Codon Usage (RSCU)
**Date Collected:** 2026-06-13

---

## Online Sources

### LIRMM — "RSCU RS: Measuring the bias in codon usage" (Rivals et al., Université de Montpellier)

**URL:** https://www.lirmm.fr/~rivals/rscu/
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference-implementation / methods page from an academic computational-biology group, restating the Sharp & Li definition)
**Retrieved by:** WebSearch query `Relative Synonymous Codon Usage RSCU formula Sharp Li 1986 definition`, then WebFetch of the URL above.

**Key Extracted Points:**

1. **Verbatim formula:** "For an amino acid i, let n_i denote the number of codons that code amino acid i. For the j-th codon of amino acid i, let x_{i,j} denote the number of occurrences of codon j. Then the RSCU for codon j of amino acid i is determined using the following formula: RSCU_{i,j} = (n_i * x_{i,j}) / Σ(x_{i,j})" — the denominator sums occurrences over all synonymous codons of amino acid i.
2. **Symbol n_i:** "the number of codons that code amino acid i" (the size of the synonymous family).
3. **Symbol x_{i,j}:** "the number of occurrences of codon j."
4. **Range:** RSCU values are "comprised between 0 and the number of synonymous codons for that amino acid."
5. **No-bias value:** when synonymous codons are used uniformly each codon yields RSCU = 1.0 (n_i identical codons → n_i/n_i = 1).

### PMC2528880 — "Analysis of synonymous codon usage and evolution of begomoviruses" (peer-reviewed, PMC)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper applying and citing the Sharp & Li RSCU definition)
**Retrieved by:** Same WebSearch as above, then WebFetch of the URL.

**Key Extracted Points:**

1. **Definition:** "Relative synonymous codon usage (RSCU) is defined as the ratio of the observed frequency of codons to the expected frequency given that all the synonymous codons for the same amino acids are used equally."
2. **No-bias value:** RSCU = 1.0 indicates no codon usage bias; values above 1.0 indicate a preferred codon, values below 1.0 indicate an under-represented codon.
3. **Attribution:** the definition is attributed to Sharp and Li (1986).

### Sharp & Li (1986) — original primary citation (Nucleic Acids Research)

**URL:** https://academic.oup.com/nar/article/14/19/7737/2385389 (DOI https://doi.org/10.1093/nar/14.19.7737)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary)
**Retrieved by:** WebSearch query `Sharp Li 1986 "codon usage in regulatory genes" Nucleic Acids Research RSCU DOI`, then WebFetch of the article landing page.

**Key Extracted Points:**

1. **Citation confirmed:** Sharp P.M., Li W.-H. (1986). "Codon usage in regulatory genes in Escherichia coli does not reflect selection for 'rare' codons." Nucleic Acids Research 14(19):7737–7749, DOI 10.1093/nar/14.19.7737 — the work that introduces RSCU as a codon-bias measure (full mathematical body behind the OUP paywall; the formula itself is taken verbatim from the LIRMM page and corroborated by PMC2528880 above).

### CodonU — reference implementation (`rscu_comp.py` / `internal_comp.py`, SouradiptoC/CodonU)

**URL:** https://github.com/SouradiptoC/CodonU/blob/master/CodonU/analyzer/internal_comp.py (function `rscu`, lines ~190–219)
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)
**Retrieved by:** `gh api repos/SouradiptoC/CodonU/contents/CodonU/analyzer/internal_comp.py?ref=master` (base64-decoded), grep for `def rscu`.

**Key Extracted Points:**

1. **Verbatim formula (code):** `rscu_dict[codon] = counts[codon] / ((len(_syn_codons[codon]) ** -1) * (sum(counts[_codon] for _codon in _syn_codons[codon])))`. This is algebraically `n_i * counts[codon] / Σ(synonymous counts)`, identical to the LIRMM formula.
2. **Docstring:** "Calculates relative synonymous codon usage (RSCU) value for a given nucleotide sequence according to Sharp and Li (1987)."
3. **Sense codons only:** iterates `unambiguous_dna_by_id[genetic_code].forward_table`, i.e. the Biopython forward table, which contains the 61 sense codons (stop codons excluded).
4. **Multiple references:** the function accepts a list of reference sequences and pools all codon occurrences (`Counter`) before computing RSCU — RSCU is computed on the aggregate codon counts across the whole reference set, not per sequence.
5. **Pseudocount (CAI-specific, NOT base RSCU):** it sets any zero-count codon to 0.5, quoting Sharp & Li (1987) p.1285 — this is the CAI adjustment to avoid log(0), distinct from the plain RSCU definition (LIRMM and PMC2528880 define RSCU with raw counts and allow RSCU = 0). The base-RSCU implementation in this unit does NOT apply the 0.5 pseudocount.

### NCBI Genetic Codes — Standard code (translation table 1)

**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
**Accessed:** 2026-06-13
**Authority rank:** 2 (official specification)
**Retrieved by:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Synonymous families (table 1):** the full 64-codon mapping was extracted. Stop codons are TAA, TAG, TGA. Single-codon amino acids are Met (ATG) and Trp (TGG) → their RSCU is always 1.0 (n_i = 1). Leucine (TTA, TTG, CTT, CTC, CTA, CTG), Arginine (CGT, CGC, CGA, CGG, AGA, AGG) and Serine (TCT, TCC, TCA, TCG, AGT, AGC) are the three six-codon families (n_i = 6).

---

## Documented Corner Cases and Failure Modes

### From CodonU (`rscu` function)

1. **Aggregation over multiple sequences:** counts are pooled across all reference sequences before computing RSCU.
2. **Sense codons only:** stop codons are excluded from the codon families (`forward_table`).
3. **Zero family count:** if a whole synonymous family is unobserved, the denominator is 0; the base RSCU is undefined for that family (the CAI pseudocount path is a separate, CAI-only convention).

### From LIRMM

1. **Single-codon amino acids:** Met and Trp have n_i = 1, so RSCU is always exactly 1.0.
2. **Range:** RSCU is bounded in [0, n_i]; RSCU = 1 means no bias for that codon.

---

## Test Datasets

### Dataset: Leucine-only worked example (derived from the LIRMM/CodonU formula)

**Source:** Formula from LIRMM (https://www.lirmm.fr/~rivals/rscu/); NCBI table 1 Leu family (n_i = 6).

CDS = `CTTCTTCTGTTA` → codons CTT, CTT, CTG, TTA (all Leucine). Family counts: CTT=2, CTG=1, TTA=1, TTG=0, CTC=0, CTA=0; Σ = 4; n_i = 6.

| Codon | x_{i,j} | RSCU = 6·x/4 |
|-------|---------|--------------|
| CTT   | 2       | 3.0          |
| CTG   | 1       | 1.5          |
| TTA   | 1       | 1.5          |
| TTG   | 0       | 0.0          |
| CTC   | 0       | 0.0          |
| CTA   | 0       | 0.0          |

Invariant check: Σ RSCU over the family = 3.0 + 1.5 + 1.5 + 0 + 0 + 0 = 6.0 = n_i.

### Dataset: Uniform-usage example (no bias → RSCU = 1.0)

**Source:** LIRMM no-bias property; PMC2528880 ("RSCU = 1.0 indicates no codon usage bias").

CDS = `TTTTTC` → codons TTT, TTC (Phenylalanine, n_i = 2). Counts TTT=1, TTC=1; Σ = 2.
RSCU(TTT) = 2·1/2 = 1.0, RSCU(TTC) = 2·1/2 = 1.0.

### Dataset: Single-codon amino acid (Met)

**Source:** NCBI table 1 (ATG = Met, only codon); LIRMM range property.

CDS = `ATGATG` → RSCU(ATG) = 1·2 / 2 = 1.0 (always 1.0 for n_i = 1).

---

## Assumptions

1. **ASSUMPTION: Genetic code defaults to Standard (NCBI table 1).** When no table is specified, the Standard code is used. Justification: table 1 is the default in Biopython/CodonU and the universal default for nuclear genes; the method exposes an overload to pass a `GeneticCode` for non-standard tables, so this is an API default, not a correctness gap.
2. **ASSUMPTION: zero-count family → RSCU 0 for every codon in that family.** The base RSCU definition (LIRMM, raw counts) leaves RSCU undefined when Σ = 0. To avoid division by zero and keep the output total over codons, an unobserved family is reported as 0.0 for each member (no codon was used, so none is "preferred"). The CAI 0.5 pseudocount is intentionally NOT applied because that is a CAI-specific convention (Sharp & Li 1987), not part of plain RSCU. This affects only families with no observations.

---

## Recommendations for Test Coverage

1. **MUST Test:** Leucine worked example returns the exact RSCU values 3.0/1.5/1.5/0/0/0 — Evidence: LIRMM formula + NCBI Leu family.
2. **MUST Test:** Uniform usage (TTT=TTC=1) returns RSCU 1.0 for both — Evidence: LIRMM/PMC2528880 no-bias = 1.0.
3. **MUST Test:** Single-codon amino acid (ATG/Met) returns RSCU 1.0 regardless of count — Evidence: NCBI table 1, n_i = 1.
4. **MUST Test:** Counts pooled across multiple input sequences — Evidence: CodonU aggregates over the reference list.
5. **MUST Test:** Stop codons excluded from output — Evidence: CodonU uses forward_table (sense codons only).
6. **MUST Test:** Σ RSCU over a synonymous family equals n_i when the family is observed — Evidence: algebraic identity of the formula (LIRMM).
7. **SHOULD Test:** lower-case input handled case-insensitively — Rationale: consistency with sibling annotator methods.
8. **SHOULD Test:** partial trailing codon (length not multiple of 3) ignored — Rationale: matches existing GetCodonUsage and reference impl reading-frame stepping.
9. **COULD Test:** null / empty input handled by documented failure mode — Rationale: edge-case completeness.

---

## References

1. Sharp P.M., Li W.-H. (1986). Codon usage in regulatory genes in Escherichia coli does not reflect selection for 'rare' codons. Nucleic Acids Research 14(19):7737–7749. https://doi.org/10.1093/nar/14.19.7737
2. Rivals E. et al. RSCU RS: Measuring the bias in codon usage (LIRMM, Université de Montpellier). https://www.lirmm.fr/~rivals/rscu/ (accessed 2026-06-13).
3. Analysis of synonymous codon usage and evolution of begomoviruses. PMC2528880. https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/ (accessed 2026-06-13).
4. SouradiptoC. CodonU, `CodonU/analyzer/internal_comp.py`, function `rscu`. https://github.com/SouradiptoC/CodonU/blob/master/CodonU/analyzer/internal_comp.py (accessed 2026-06-13).
5. NCBI. The Genetic Codes — Standard Code (transl_table=1). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
