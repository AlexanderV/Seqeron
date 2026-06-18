# Evidence Artifact: ONCO-SIG-004

**Test Unit ID:** ONCO-SIG-004
**Algorithm:** Mutational Process Classification (SBS exposure → active mutational processes)
**Date Collected:** 2026-06-14

---

## Online Sources

### COSMIC Mutational Signatures — SBS (Single Base Substitution) catalogue

**URL:** https://cancer.sanger.ac.uk/signatures/sbs/
**Accessed:** 2026-06-14
**Authority rank:** 5 (curated reference database; the canonical registry of SBS signatures and their proposed aetiologies; underpinned by Alexandrov et al. 2020)

**Key Extracted Points (proposed aetiology quoted verbatim from the per-signature pages):**

(Retrieved by fetching the SBS index page and asking for the proposed aetiology text of each named signature; the page returned the quoted aetiology strings below.)

1. **SBS1 aetiology:** "Spontaneous deamination of 5-methylcytosine (clock-like signature)" — maps to the **Aging / clock-like** process.
2. **SBS5 aetiology:** "Unknown (clock-like signature)" — also a **clock-like** process; commonly grouped with SBS1 as age-correlated clock-like activity.
3. **SBS2 aetiology:** "Activity of APOBEC family of cytidine deaminases" — **APOBEC** process.
4. **SBS13 aetiology:** "Activity of APOBEC family of cytidine deaminases" — **APOBEC** process.
5. **SBS4 aetiology:** "Tobacco smoking" — **Tobacco smoking** process.
6. **SBS7a aetiology:** "Ultraviolet light exposure"; **SBS7b/7c/7d** likewise "Ultraviolet light exposure" — **UV** process.
7. **SBS6 aetiology:** "Defective DNA mismatch repair" — **MMR deficiency** process.
8. **SBS15 aetiology:** "Defective DNA mismatch repair" — **MMR deficiency** process.
9. **SBS20 aetiology:** "Concurrent POLD1 mutations and defective DNA mismatch repair" — **MMR deficiency** process (involves defective MMR).
10. **SBS26 aetiology:** "Defective DNA mismatch repair" — **MMR deficiency** process.

### deconstructSigs (Rosenthal et al. 2016) — relative-contribution cutoff for signature presence

**URL (paper):** https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
**URL (reference source):** https://github.com/raerose01/deconstructSigs/blob/master/R/whichSignatures.R
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed: *Genome Biology* 2016, 17:31) / 3 (reference implementation source)

**Key Extracted Points:**

1. **Weights are normalized relative contributions:** "Finally, the weights W are normalized between 0 and 1 …" — each signature's reported activity is its fraction of the total reconstructed mutation burden (a relative contribution / proportion). (Retrieved by fetching the PMC article and asking how contributions are reported and normalized.)
2. **Presence cutoff = 6%:** "… and any signature with Wᵢ < 6% is excluded." A signature is reported as **present/active** only when its normalized relative contribution Wᵢ ≥ 0.06; below that it is set to zero. (Same fetch.)
3. **Reference-implementation cutoff (verbatim source):** the `whichSignatures` function declares the default parameter `signature.cutoff = 0.06`, and applies it with the line `weights[weights < signature.cutoff ] <- 0`. (Retrieved by fetching the raw `whichSignatures.R`; both the default value and the cutoff line were returned verbatim.)
4. **False-negative rate of the cutoff:** the 6% exclusion threshold "only resulted in 38 instances where a signature was incorrectly excluded for a false negative rate of 1.4%" across 2,646 simulated signatures — i.e. the 6% cutoff is the source-recommended, empirically calibrated threshold for declaring a signature absent. (Retrieved from the PMC article methods/results.)

### Alexandrov et al. (2020) — The repertoire of mutational signatures in human cancer (primary aetiology reference)

**URL:** https://pubmed.ncbi.nlm.nih.gov/32025018/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed: *Nature* 578(7793):94–101, 2020; DOI 10.1038/s41586-020-1943-3)

**Key Extracted Points:**

1. **Primary signature catalogue:** the paper characterises 4,645 whole-genome and 19,184 exome sequences and identifies 81 single-base-substitution, doublet-base-substitution and indel mutational signatures — the published basis behind the COSMIC SBS aetiology assignments cited above. (Retrieved via search of the PubMed record; citation and scope confirmed.)

---

## Documented Corner Cases and Failure Modes

### From deconstructSigs (Rosenthal et al. 2016)

1. **Sub-cutoff signatures dropped:** any signature whose normalized contribution is below 6% is forced to zero, so the surviving contributions can sum to **less than 1** (the missing mass is attributed to "unknown"). A process must therefore be declared active only from contributions that survive the 6% cutoff.
2. **Tie / multiple active processes:** a tumour can show several active processes simultaneously (e.g. APOBEC SBS2+SBS13 plus aging SBS1); classification reports the full active set, not a single label.

### From COSMIC / Alexandrov 2020

1. **Unknown-aetiology signatures:** several common signatures (e.g. SBS5 "Unknown (clock-like)") have only a partially known aetiology; SBS5 is nonetheless a clock-like process. Signatures absent from the caller-supplied label→process map contribute to none of the named processes.

---

## Test Datasets

### Dataset: Hand-derived normalized-contribution classification

**Source:** Derived directly from the deconstructSigs 6% cutoff rule (Rosenthal et al. 2016) and the COSMIC SBS→aetiology map. Reference signature **profiles** are caller-supplied and NOT fabricated here; only the per-signature exposures (activities) and their COSMIC labels are used.

| Signature | Raw exposure | Normalized contribution | ≥ 0.06 cutoff? | Process |
|-----------|--------------|--------------------------|----------------|---------|
| SBS2  | 50 | 50/100 = 0.50 | yes | APOBEC |
| SBS13 | 30 | 30/100 = 0.30 | yes | APOBEC |
| SBS1  | 15 | 15/100 = 0.15 | yes | Aging |
| SBS4  | 5  | 5/100 = 0.05  | **no** (0.05 < 0.06) | Tobacco smoking (excluded) |

Total raw = 100. Per-process active contribution (sum of surviving signatures): APOBEC = 0.50 + 0.30 = 0.80; Aging = 0.15; Tobacco = 0 (SBS4 below cutoff). Active processes = {APOBEC, Aging}; **dominant process = APOBEC (0.80)**.

---

## Assumptions

1. **ASSUMPTION: Per-process aggregation by summation.** When more than one signature maps to the same process (e.g. SBS2 and SBS13 → APOBEC; SBS6/15/20/26 → MMR deficiency), the process's relative contribution is the **sum** of its member signatures' surviving normalized contributions. COSMIC defines the per-signature aetiologies but does not prescribe an aggregation rule; summation is the natural and standard interpretation of additive relative contributions (deconstructSigs weights are additive fractions of one reconstruction). This affects only the per-process total, not the per-signature cutoff decision.
2. **ASSUMPTION: Cutoff applied to per-signature normalized contribution (not per-process total).** Following deconstructSigs exactly, the 6% cutoff is applied to each individual signature's normalized contribution; surviving contributions are then aggregated by process. (deconstructSigs operates per signature; processes are a downstream grouping.)

---

## Recommendations for Test Coverage

1. **MUST Test:** Normalized contributions computed as exposureᵢ / Σ exposure. — Evidence: deconstructSigs "weights W are normalized between 0 and 1" (Rosenthal 2016).
2. **MUST Test:** A signature with normalized contribution < 0.06 is excluded (set to 0); ≥ 0.06 is retained. — Evidence: `weights[weights < signature.cutoff] <- 0`, `signature.cutoff = 0.06` (deconstructSigs source).
3. **MUST Test:** Boundary at exactly 0.06 is retained (cutoff uses strict `<`). — Evidence: `weights < signature.cutoff` (strict less-than, deconstructSigs source).
4. **MUST Test:** SBS labels map to the correct process per COSMIC aetiology (SBS1/5→Aging, SBS2/13→APOBEC, SBS4→Tobacco, SBS7a–d→UV, SBS6/15/20/26→MMR). — Evidence: COSMIC SBS aetiology strings.
5. **MUST Test:** Per-process contribution = sum of surviving member-signature contributions; dominant process = max per-process contribution. — Evidence: additive weights (deconstructSigs) + hand-derived dataset above.
6. **MUST Test:** Active-process set excludes processes all of whose signatures fell below the cutoff. — Evidence: hand-derived dataset (Tobacco excluded).
7. **SHOULD Test:** All-zero exposures ⇒ no active processes, no dominant process. — Rationale: degenerate input; Σ = 0 so no normalization possible.
8. **SHOULD Test:** Signature label not present in the process map contributes to no process. — Rationale: COSMIC unknown/unmapped signatures.
9. **COULD Test:** Custom cutoff parameter overrides the 0.06 default. — Rationale: deconstructSigs exposes `signature.cutoff` as a parameter.

---

## References

1. Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C. (2016). deconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. *Genome Biology* 17:31. https://doi.org/10.1186/s13059-016-0893-4
2. deconstructSigs reference implementation, `whichSignatures.R` (default `signature.cutoff = 0.06`; `weights[weights < signature.cutoff] <- 0`). https://github.com/raerose01/deconstructSigs/blob/master/R/whichSignatures.R
3. COSMIC Mutational Signatures — SBS (proposed aetiologies). Wellcome Sanger Institute. https://cancer.sanger.ac.uk/signatures/sbs/
4. Alexandrov LB, Kim J, Haradhvala NJ, et al. (2020). The repertoire of mutational signatures in human cancer. *Nature* 578(7793):94–101. https://doi.org/10.1038/s41586-020-1943-3

---

## Change History

- **2026-06-14**: Initial documentation.
