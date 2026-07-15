---
type: source
title: "Evidence: CODON-RSCU-001 (Relative Synonymous Codon Usage — RSCU + codon counting)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-RSCU-001-Evidence.md
sources:
  - docs/Evidence/CODON-RSCU-001-Evidence.md
source_commit: 458398bca4eee7e7fa828acbe182e07695db5e28
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-RSCU-001

The validation-evidence artifact for test unit **CODON-RSCU-001** (Relative Synonymous Codon
Usage — RSCU **and codon counting**). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is summarized in
[[relative-synonymous-codon-usage]]. This is the **second RSCU evidence unit** — it validates the
same `RSCU_{i,j} = n_i·x_{i,j}/Σ x` measure as [[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]]
but adds the supporting `CountCodons` operation and a broader reference-implementation panel. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13) — **Sharp, Tuohy & Mosurski (1986)**, *Nucleic Acids
  Res.* 14(13):5125-5143 (PMID 3526280, DOI 10.1093/nar/14.13.5125), named as the **primary paper
  that introduced RSCU** (authority rank 1; PMC renders the methods equations as images, so the
  formula was transcribed from the reference implementations below that cite it). The **LIRMM /
  Suzuki et al. RSCU RS** page for the verbatim `RSCU_i = X_i / ((1/N_i) Σ X_j)` formula and the
  no-bias = 1 property. **GenomicSig** (CRAN) `RSCU` for the equivalent indexed form
  `RSCU_{i,j} = (n_i · x_{i,j}) / Σ x_{i,j}` and the explicit **[0, n_i] bounds**. **seqinr** `uco`
  (Charif & Lobry) for the observed-over-expected-under-uniform definition, the 1.00-no-bias reading,
  and the Sharp/Tuohy/Mosurski 1986 citation. **cubar** `est_rscu` for the `pseudo_cnt` (default 1)
  zero-division guard and the `incl_stop` default-FALSE convention. **PMC2528880** (begomovirus codon
  usage) as an independent restatement of the definition (there attributed to "Sharp & Li 1986").
- **Algorithm spec** — for codon `j` of amino acid `i` with family size `n_i`,
  `RSCU = n_i · x_{i,j} / Σ x_{i,j}`; equivalently observed frequency ÷ frequency expected under
  uniform synonymous usage. Uniform ⇒ 1.0, >1 over-used, <1 under-used, bounded **[0, n_i]** (max
  `n_i` when only one synonymous codon is used); the family's RSCU values sum to `n_i`.
- **`CountCodons` contract** — codons are **non-overlapping triplets from offset 0**; trailing 1–2
  bases are ignored; any triplet containing a non-`{A,C,G,T}` character is **excluded** (string
  overload uppercases first). Kazusa CUTG counting convention.
- **Datasets** — Leu 6-fold `CTGCTGCTGCTA` (CTG×3, CTA×1) ⇒ RSCU(CTG)=**4.5**, RSCU(CTA)=**1.5**,
  the other four = 0, Σ = 6 = n_i; Phe 2-fold `TTTTTTTTC` (TTT×2, TTC×1) ⇒ **4/3** and **2/3**,
  Σ = 2; unbiased `TTTTTC` ⇒ **1.0 / 1.0**; single-codon Met `ATGATG` ⇒ **1.0**; `CountCodons`
  frame/exclusion cases `ATGAAATGA`→3 triplets, `ATGAA`→ATG only (trailing `AA` dropped),
  `ATGNNNAAA`→ATG+AAA (`NNN` excluded).
- **Recommended coverage** — MUST tests for the Leu / Phe / unbiased / single-codon RSCU values and
  the `CountCodons` frame/exclusion behaviour; SHOULD tests for null → `ArgumentNullException`,
  empty → empty dictionary, and the Σ-over-family = n_i invariant; COULD test for lowercase input.

## Assumptions (from the artifact)

1. **Absent-family 0/0 ⇒ RSCU 0.** When a synonymous family has zero observed codons the canonical
   denominator is 0 (undefined). cubar avoids this with a pseudocount (default 1); Seqeron instead
   returns 0.0 for every codon of an absent family. Only affects families that never occur — never a
   real observed RSCU value. Documented as a convention, not a divergence from the canonical formula.
2. **Stop codons treated as a synonymous family of size 3.** This artifact states the repository
   groups TAA/TAG/TGA as one degeneracy-3 family and computes their RSCU like any other family
   (reference tools such as cubar/seqinr instead exclude stops via `incl_stop = FALSE`). Either way
   this **does not change the RSCU of any amino-acid codon**.

## Contradictions and cross-page nuances

- **No source contradictions** — the LIRMM, GenomicSig, seqinr and Sharp/Tuohy/Mosurski forms are
  algebraically identical; cubar's pseudocount is an explicitly documented zero-division convention.
- **Cross-page nuance (flagged, not a source conflict):** [[relative-synonymous-codon-usage]] — built
  from [[annot-codonusage-001-evidence|ANNOT-CODONUSAGE-001]] — describes stop codons as **excluded**
  from the families (Biopython `forward_table` convention); **this** artifact says the repository
  treats the three stops as a **degeneracy-3 synonymous family**. The two descriptions agree on the
  only observable consequence (no effect on amino-acid-codon RSCU) but differ on the framing of stop
  handling; noted on the concept page.
- **Citation nuance:** this unit names the RSCU-introducing paper as **Sharp, Tuohy & Mosurski
  (1986)**, NAR 14(13):5125-5143; the concept and ANNOT-CODONUSAGE-001 wrote "Sharp & Li (1986)".
  Both 1986 attributions circulate in the literature (the begomovirus restatement here also says
  "Sharp & Li 1986"); recorded on the concept as the more precise primary attribution.
