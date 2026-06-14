# Evidence Artifact: ONCO-HLA-001

**Test Unit ID:** ONCO-HLA-001
**Algorithm:** HLA allele nomenclature parsing/validation + allele-specific HLA LOH (LOHHLA) classification
**Date Collected:** 2026-06-15

---

## Online Sources

### WHO HLA Nomenclature — "Naming Alleles" (IPD-IMGT/HLA, hla.alleles.org)

**URL:** https://hla.alleles.org/pages/nomenclature/naming_alleles/
**Accessed:** 2026-06-15
**Authority rank:** 2 (official nomenclature standard of the WHO Nomenclature Committee for Factors of the HLA System)
**How retrieved:** WebSearch query `hla.alleles.org nomenclature naming HLA allele four fields colon separated allele group specific protein` → WebFetch of the returned canonical page URL above.

**Key Extracted Points:**

1. **Name structure:** An allele name has the form `HLA-[Gene]*[Field1]:[Field2][:Field3][:Field4]`. The gene name is separated from the numeric fields by an asterisk `*`; the numeric fields are separated by colons.
2. **Field 1 (type / allele group):** "The digits before the first colon describe the type, which often corresponds to the serological antigen carried by an allotype."
3. **Field 2 (specific HLA protein / subtype):** "The next set of digits are used to list the subtypes, numbers being assigned in the order in which DNA sequences have been determined." Field 2 denotes the specific HLA protein.
4. **Field 3 (synonymous coding substitutions):** "Alleles that differ only by synonymous nucleotide substitutions … within the coding sequence are distinguished by the use of the third set of digits."
5. **Field 4 (non-coding differences):** "Alleles that only differ by sequence polymorphisms in the introns, or in the 5′ or 3′ untranslated regions … are distinguished by the use of the fourth set of digits."
6. **Minimum designation:** "All alleles receive at least a four digit name, which corresponds to the first two sets of digits, longer names are only assigned when necessary." → a valid allele name has **at least two fields** (Field1:Field2).
7. **Expression suffixes:** an optional trailing letter encodes expression status — `N` Null (not expressed), `L` Low cell-surface expression, `S` Secreted (soluble) only, `C` Cytoplasm only, `A` Aberrant expression, `Q` Questionable expression.

### Marsh et al. (2010) — "Nomenclature for factors of the HLA system, 2010", Tissue Antigens

**URL:** https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1399-0039.2010.01466.x
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed publication of the WHO Nomenclature Committee, the originating reference for the colon-delimited HLA naming convention)
**How retrieved:** WebSearch query `Marsh 2010 Nomenclature factors HLA system 2010 Tissue Antigens allele field naming HLA-A*02:01` → the Wiley record was returned and identifies the article (Tissue Antigens, 2010, 75(4):291–455).

**Key Extracted Points:**

1. **Colon-delimited fields:** Each HLA allele name has a unique number corresponding to up to four sets of digits separated by colons; the length of the designation depends on the allele sequence relative to its nearest relative.
2. **Minimum two fields:** all alleles receive at least a four-digit name corresponding to the first two sets of digits (consistent with the hla.alleles.org statement above).

### McGranahan et al. (2017) — "Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution", Cell (LOHHLA)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper, *Cell* 171(6):1259–1271; introduces the LOHHLA method)
**How retrieved:** WebSearch query `McGranahan 2017 Cell LOHHLA allele-specific HLA loss of heterozygosity copy number allelic imbalance` → WebFetch of the PMC full-text URL above (after following the ncbi.nlm.nih.gov → pmc.ncbi.nlm.nih.gov 301 redirect).

**Key Extracted Points:**

1. **Per-allele copy number:** LOHHLA infers, for each HLA gene, the allele-specific copy number of both homologous alleles from tumor coverage relative to germline (logR) and B-allele frequencies (BAF) at polymorphic sites that distinguish the two alleles.
2. **Loss threshold (verbatim):** "A copy number < 0.5, is classified as subject to loss, and thereby indicative of LOH." → an HLA allele whose estimated copy number is **< 0.5** is the *lost* allele.
3. **Allelic imbalance guard (verbatim):** "To avoid over-calling LOH, we also calculate a p value relating to allelic imbalance for each HLA gene. Allelic imbalance is determined if **p < 0.01** using the paired Student's t-Test." → HLA LOH is called only when there is significant allelic imbalance (paired t-test **p < 0.01**) AND one allele has copy number < 0.5.

### LOHHLA reference implementation (mskcc/lohhla, `LOHHLAscript.R`)

**URL:** https://raw.githubusercontent.com/mskcc/lohhla/master/LOHHLAscript.R
**Accessed:** 2026-06-15
**Authority rank:** 3 (reference implementation accompanying the McGranahan 2017 paper)
**How retrieved:** WebFetch of the raw GitHub source file URL above.

**Key Extracted Points:**

1. **Paired t-test for imbalance (verbatim):** `PairedTtest <- t.test(tmpOut[,2],tmpOut[,4],paired=TRUE)` and `PVal <- PairedTtest$p.value` — confirms allelic imbalance is assessed by a paired Student's t-test and the p-value drives the call.
2. **Per-allele copy-number variables:** `HLA_type1copyNum_withoutBAF`, `HLA_type1copyNum_withBAF`, `HLA_type2copyNum_withoutBAF`, `HLA_type2copyNum_withBAF` — confirms LOHHLA reports an estimated copy number for each of the two homologous HLA alleles, which is the input our classifier consumes.

---

## Documented Corner Cases and Failure Modes

### From hla.alleles.org / Marsh 2010

1. **Two-field minimum:** a designation with only a first field (e.g. `HLA-A*02`) is not a complete allele name; a valid allele has at least two fields.
2. **Trailing expression suffix:** a single letter (N/L/S/C/A/Q) may follow the last numeric field and must be accepted; any other trailing letter is invalid.
3. **Variable field count:** 2, 3, or 4 fields are all valid; more than 4 is not.

### From McGranahan 2017 (LOHHLA)

1. **No-LOH despite low copy:** if allelic imbalance is not significant (paired t-test p ≥ 0.01), LOH is *not* called even if a raw allele copy number dips below 0.5 — this is the explicit "avoid over-calling LOH" guard.
2. **Both alleles retained:** when both homologous alleles have copy number ≥ 0.5, the locus is heterozygous-retained (no LOH).
3. **Homozygous locus:** the method estimates copy number per *distinct* allele; a locus with two identical alleles cannot be assessed for allele-specific loss (no polymorphic sites distinguish the homologs).

---

## Test Datasets

### Dataset: HLA nomenclature examples (WHO / hla.alleles.org)

**Source:** hla.alleles.org "Naming Alleles"; Marsh et al. (2010), Tissue Antigens 75(4):291–455.

| Input | Valid? | Gene | Field1 | Field2 | Field3 | Field4 | Suffix |
|-------|--------|------|--------|--------|--------|--------|--------|
| `HLA-A*02:01` | yes | A | 02 | 01 | — | — | — |
| `HLA-B*07:02:01` | yes | B | 07 | 02 | 01 | — | — |
| `HLA-C*07:02:01:03` | yes | C | 07 | 02 | 01 | 03 | — |
| `HLA-A*24:02:01:02L` | yes | A | 24 | 02 | 01 | 02 | L |
| `HLA-A*02` | no (only one field; min is two) | — | — | — | — | — | — |
| `A*02:01` | no (missing `HLA-` prefix) | — | — | — | — | — | — |
| `HLA-A*02:01:01:01:01` | no (five fields; max is four) | — | — | — | — | — | — |
| `HLA-A*02:01X` | no (X is not a valid expression suffix) | — | — | — | — | — | — |

### Dataset: HLA allele-specific LOH (LOHHLA thresholds, McGranahan 2017)

**Source:** McGranahan et al. (2017), Cell 171(6):1259–1271 (PMC5720478): allele copy number < 0.5 → loss; allelic-imbalance paired t-test p < 0.01.

| Case | Allele1 CN | Allele2 CN | Imbalance p | HLA LOH? | Lost allele |
|------|-----------|-----------|-------------|----------|-------------|
| Clear loss of allele 2 | 1.8 | 0.30 | 0.001 | yes | allele 2 |
| Clear loss of allele 1 | 0.10 | 1.50 | 0.0005 | yes | allele 1 |
| Both retained | 1.10 | 0.90 | 0.30 | no | — |
| Low CN but no significant imbalance | 1.60 | 0.40 | 0.05 | no (p ≥ 0.01 guard) | — |
| Boundary CN exactly 0.5 | 1.50 | 0.50 | 0.001 | no (0.5 is not < 0.5) | — |
| Boundary p exactly 0.01 | 1.70 | 0.40 | 0.01 | no (0.01 is not < 0.01) | — |

---

## Assumptions

1. **ASSUMPTION: lost-allele tie-break** — McGranahan 2017 defines a *lost* allele by copy number < 0.5 but does not specify behaviour if both alleles were < 0.5 with significant imbalance (biologically a homozygous deletion, not allele-specific LOH). We classify LOH only when exactly one allele is < 0.5; if both are < 0.5 we report it as `HomozygousLoss` (not allele-specific LOH). Only the label is affected; the two source-exact thresholds (0.5, 0.01) are unchanged.

---

## Recommendations for Test Coverage

1. **MUST Test:** Parse each valid nomenclature example into the correct gene + field tuple + suffix. — Evidence: hla.alleles.org "Naming Alleles"; Marsh 2010.
2. **MUST Test:** Reject malformed names (missing `HLA-` prefix, single field, five fields, invalid suffix letter, non-numeric field). — Evidence: hla.alleles.org two-field minimum / four-field maximum / suffix set.
3. **MUST Test:** Call HLA LOH iff one allele CN < 0.5 AND imbalance p < 0.01; verify the lost allele and the two boundary cases (CN = 0.5, p = 0.01). — Evidence: McGranahan 2017 verbatim thresholds.
4. **SHOULD Test:** Both-retained and no-significant-imbalance cases return no LOH. — Rationale: the explicit over-calling guard.
5. **COULD Test:** Both alleles < 0.5 → HomozygousLoss label. — Rationale: documented assumption boundary.

---

## References

1. WHO Nomenclature Committee for Factors of the HLA System. Naming Alleles. IPD-IMGT/HLA. https://hla.alleles.org/pages/nomenclature/naming_alleles/ (accessed 2026-06-15).
2. Marsh SGE, Albert ED, Bodmer WF, et al. (2010). Nomenclature for factors of the HLA system, 2010. Tissue Antigens 75(4):291–455. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1399-0039.2010.01466.x
3. McGranahan N, Rosenthal R, Hiley CT, et al. (2017). Allele-Specific HLA Loss and Immune Escape in Lung Cancer Evolution. Cell 171(6):1259–1271. https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/ (DOI: 10.1016/j.cell.2017.10.001)
4. mskcc/lohhla. LOHHLAscript.R (reference implementation). https://raw.githubusercontent.com/mskcc/lohhla/master/LOHHLAscript.R (accessed 2026-06-15).

---

## Change History

- **2026-06-15**: Initial documentation.
