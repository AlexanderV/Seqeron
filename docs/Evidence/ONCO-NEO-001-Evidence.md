# Evidence Artifact: ONCO-NEO-001

**Test Unit ID:** ONCO-NEO-001
**Algorithm:** Neoantigen Candidate Peptide Window Generation (somatic missense mutation)
**Date Collected:** 2026-06-14

---

## Online Sources

### pVACtools (Hundal et al. 2020), Cancer Immunology Research 8(3):409–420

**URL:** https://aacrjournals.org/cancerimmunolres/article/8/3/409/469797 (search result summary, AACR);
discovery via WebSearch "pVACtools neoantigen peptide window generation 8-11 mer mutant wildtype agretope Hundal 2020"
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper) / 3 (the tool is also a reference implementation)

**How retrieved:** WebSearch query above returned the AACR article and the biorxiv preprint; the AACR
full-text and biorxiv PDF returned HTTP 403, so the extracted facts below are taken from the indexed
full-text summary returned by WebSearch (verbatim phrases quoted) and corroborated by the readthedocs and
ProGeo-neo sources that WERE fetched in full.

**Key Extracted Points:**

1. **Class I / Class II peptide lengths:** "The tools predict the strongest MHC-binding peptides (8–11-mer
   for Class I MHC and 13–25-mer for Class II) using one or more prediction algorithms." (verbatim)
2. **Supported variant mechanisms:** pVACtools "supports identification of altered peptides from different
   mechanisms, including point mutations, inframe and frameshift insertions and deletions, and gene fusions."
3. **Agretopicity (mutant vs wild-type):** predictions are ranked by, among others, "differential binding
   compared to the wild type peptide (agretopicity index score)"; the tools "can filter predictions where the
   predicted binding affinities are lower than each corresponding wild-type peptide affinity, indicating that
   the mutant version of the peptide is a stronger binder."

### pVACtools readthedocs — Optional Downstream Analysis Tools / Features

**URL:** https://pvactools.readthedocs.io/en/latest/pvacseq/optional_downstream_analysis_tools.html and
.../pvacseq/features.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (official reference-implementation documentation)

**How retrieved:** WebFetch of the readthedocs RST source
(https://raw.githubusercontent.com/griffithlab/pVACtools/master/docs/pvacseq/optional_downstream_analysis_tools.rst)
and the features.html page.

**Key Extracted Points:**

1. **Flanking construction:** "The `flanking_sequence_length` positional parameter controls how many amino
   acids will be included on either side of the mutation." (verbatim) and "The alteration in the VCF (e.g. a
   somatic missense SNV) will be centered in the protein sequence returned (if possible)." (verbatim)
2. **Class I lengths used in benchmarking:** "we have run predictions for all class I algorithms supported by
   pVACtools on 100,000 reference peptides each in lengths 8-11" (verbatim, features.html).

### ProGeo-neo (Li et al. 2020), BMC Medical Genomics 13:52

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7118832/ ; DOI https://doi.org/10.1186/s12920-020-0683-4
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper)

**How retrieved:** WebSearch "neoantigen prediction missense mutation peptide windows spanning mutated residue
MHC class I 8-11mers" → PMC link; WebFetch of the PMC article and the BMC DOI page.

**Key Extracted Points:**

1. **21-mer ±10-flank construction (verbatim):** "Sequences corresponding to each of the coding missense
   mutations that would cause amino-acid substitutions were translated into a 21-mer amino acid fasta
   sequence, with 10 amino acids flanking the substituted amino acid on each side."
2. **Window extraction:** from these 21-mers the binding predictor extracts and evaluates the shorter 8–11-mer
   epitopes; a 21-mer with the substitution centred (10 flanks each side) contains exactly every 8–11-mer
   window that overlaps the mutated residue.

### NetMHCpan-4.1 service / NetMHCpan-4.0 (Jurtz et al. 2017)

**URL:** https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/ ;
Jurtz et al. (2017) J. Immunol. 199(9):3360–3368, https://doi.org/10.4049/jimmunol.1700893
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation) / 1 (paper)

**How retrieved:** WebSearch "netMHCpan 4.1 peptide length 8 9 10 11 mers" → DTU service page; WebFetch of the
service page.

**Key Extracted Points:**

1. **Length options (verbatim list from the service):** "8mer peptides 9mer peptides 10mer peptides 11mer
   peptides 12mer peptides 13mer peptides 14mer peptides"; "Predictions can be made for peptides of any
   length." Class I canonical neoantigen search uses 8–11 (per pVACtools).
2. **Dominant length:** "most presented MHC class I ligands are of length 9 amino acids" (NetMHCpan-4.1
   description), with length preference varying by allele — motivating enumeration of 8–11 rather than 9 only.

### Wells et al. 2020 (TESLA consortium), Cell 183(3):818–834

**URL:** discovery via WebSearch; DOI https://doi.org/10.1016/j.cell.2020.09.015
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper)

**How retrieved:** WebSearch "...TESLA Wells 2020" returned the citation and the differential-agretopic-index
description.

**Key Extracted Points:**

1. **Differential agretopicity:** "The differential agretopic index (DAI) — the ratio of peptide binding
   affinity of germline peptide to matched mutated peptide — has been shown to be a strong signal of point
   mutation immunogenicity." Confirms the mutant peptide must be paired with its matched wild-type (germline)
   counterpart at the same coordinates — the agretope — which this unit produces.

---

## Documented Corner Cases and Failure Modes

### From ProGeo-neo (Li et al. 2020) / pVACtools (Hundal et al. 2020)

1. **Mutation near a protein terminus:** when the substitution is within fewer than k−1 residues of either end,
   fewer than k windows of length k exist (the 21-mer is built "if possible"; pVACtools centres "if possible").
   The set is the windows that still fit within the protein while spanning the mutation.
2. **Non-coding / no amino-acid change:** variants that do not produce an amino-acid substitution yield no
   candidate peptide (pVACtools generates peptides only from protein-altering variants).

### From NetMHCpan (Jurtz et al. 2017)

1. **Peptide shorter than 8 / longer than 11:** outside the class I canonical search band; not enumerated here.

---

## Test Datasets

### Dataset: Worked window enumeration (derived from the windowing definition)

**Source:** Windowing definition — ProGeo-neo (Li et al. 2020) + pVACtools class I lengths 8–11 (Hundal 2020).
Wild-type protein and mutation chosen for this unit; expected windows derived by hand from the definition
"every length-k window of the mutant protein that spans the mutated residue", NOT from code output.

Wild-type protein `MKTAYIAKQRSTVWLNDEFGH` (length 21), missense `Y5C` (position 5, Y→C).
Mutant protein `MKTACIAKQRSTVWLNDEFGH`.

| Parameter | Value |
|-----------|-------|
| Protein length L | 21 |
| Mutation position p (1-based) | 5 |
| Wild-type residue | Y |
| Mutant residue | C |

For length k, valid 0-based start s ∈ [max(0, p−1−k+1), min(p−1, L−k)]; count = last − first + 1.

| k | first start (1-based) | last start (1-based) | window count | mutant windows |
|---|----------------------|----------------------|--------------|----------------|
| 8 | 1 (s0=0) | 5 (s0=4) | 5 | MKTACIAK, KTACIAKQ, TACIAKQR, ACIAKQRS, CIAKQRST |
| 9 | 1 | 5 | 5 | MKTACIAKQ … CIAKQRSTV |
| 10 | 1 | 5 | 5 | MKTACIAKQR … CIAKQRSTVW |
| 11 | 1 | 5 | 5 | MKTACIAKQRS … CIAKQRSTVWL |

Total candidate peptides over k = 8..11: 5 + 5 + 5 + 5 = **20**.
Each wild-type peptide equals its mutant peptide with `C`→`Y` at the mutation offset.

### Dataset: Terminal mutation (truncated windows)

Wild-type protein `MKTAYIAKQRSTVWLNDEFGH` (length 21), missense `M1V` (position 1, M→V).
For k=9, valid start s ∈ [max(0, −8)=0, min(0, 12)=0] → 1 window only: mutant `VKTAYIAKQ`, WT `MKTAYIAKQ`,
mutation offset 0. (Only one window spans a terminal residue per length.)

---

## Assumptions

1. **ASSUMPTION: Single-residue substitution only.** This unit implements candidate-peptide windowing for a
   somatic missense SNV (one amino-acid substitution). Frameshift / indel / fusion neopeptides are out of
   scope (pVACtools supports them through a different translation step). Documented, not a correctness gap for
   the missense case.
2. **ASSUMPTION: Binding affinity out of scope.** IC50 / binding-rank scoring requires a trained MHC model
   (NetMHCpan weights) and is caller-supplied (ONCO-MHC-001). Not fabricated here. Source-backed decision.

---

## Recommendations for Test Coverage

1. **MUST Test:** interior missense — all k∈[8,11] windows spanning the mutation are produced (5 per length,
   20 total), correct mutant/WT pairing and mutation offsets. — Evidence: ProGeo-neo (Li 2020), pVACtools (Hundal 2020).
2. **MUST Test:** mutant peptide differs from WT only at the mutation offset, and the WT residue is the
   original. — Evidence: agretope pairing (Wells 2020; Hundal 2020).
3. **MUST Test:** terminal mutation produces only the windows that fit (truncated count). — Evidence: ProGeo-neo "if possible".
4. **MUST Test:** single length k=9 yields exactly k windows for a sufficiently interior mutation. — Evidence: windowing definition.
5. **SHOULD Test:** null protein → ArgumentNullException; empty protein, mutant==wildtype, invalid length range,
   out-of-range position → exceptions. — Rationale: documented failure modes / input validation.
6. **COULD Test:** protein shorter than some requested k skips that length but still returns shorter windows. — Rationale: bounds handling.

---

## References

1. Hundal J, Kiwala S, McMichael J, et al. (2020). pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens. Cancer Immunol Res. 8(3):409–420. https://doi.org/10.1158/2326-6066.CIR-19-0401
2. Li Y, Wang G, Tan X, et al. (2020). ProGeo-neo: a customized proteogenomic workflow for neoantigen prediction and selection. BMC Med Genomics. 13:52. https://doi.org/10.1186/s12920-020-0683-4
3. Jurtz V, Paul S, Andreatta M, et al. (2017). NetMHCpan-4.0: Improved Peptide-MHC Class I Interaction Predictions Integrating Eluted Ligand and Peptide Binding Affinity Data. J Immunol. 199(9):3360–3368. https://doi.org/10.4049/jimmunol.1700893
4. NetMHCpan-4.1 web service. DTU Health Tech. https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/ (accessed 2026-06-14)
5. Wells DK, van Buuren MM, Dang KK, et al. (2020). Key Parameters of Tumor Epitope Immunogenicity Revealed Through a Consortium Approach Improve Neoantigen Prediction (TESLA). Cell. 183(3):818–834. https://doi.org/10.1016/j.cell.2020.09.015

---

## Change History

- **2026-06-14**: Initial documentation.
