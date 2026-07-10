---
type: concept
title: "Variant effect annotation (VEP-style consequence + Sequence-Ontology IMPACT)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/VARIANT-ANNOT-001-Evidence.md
  - docs/Evidence/VARIANT-CALL-001-Evidence.md
source_commit: 5b4dd805db54d51bae30445a884e122fc4d97bd5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: variant-annot-001-evidence
      evidence: "Test Unit ID: VARIANT-ANNOT-001 — Variant Annotation — functional impact / consequence prediction (VEP / Sequence Ontology)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:genetic-code-translation
      source: variant-annot-001-evidence
      evidence: "VariationEffect.pm peptide predicates compare ref_pep/alt_pep obtained by translating the reference and alternate codons through NCBI Standard code table 1 (gc.prt); synonymous/missense/stop_gained/stop_lost are all decided from those translations"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:somatic-variant-calling-tumor-normal
      source: variant-annot-001-evidence
      evidence: "Variant annotation is the interpretation layer downstream of variant calling — it takes an already-called variant and predicts its functional consequence"
      confidence: medium
      status: current
---

# Variant effect annotation (VEP-style consequence + Sequence-Ontology IMPACT)

**Variant effect annotation** maps an already-*called* sequence variant to **what it does** — a
standardized **functional consequence** and an **impact severity** — reproducing the **Ensembl
Variant Effect Predictor (VEP)**. It is the **interpretation** step of the variant-analysis family:
calling produces the variant, annotation says whether it is a `missense_variant`, a `stop_gained`, a
`frameshift_variant`, an `intron_variant`, and so on, each carrying its **Sequence Ontology (SO)** term,
accession, and an **IMPACT** class of HIGH / MODERATE / LOW / MODIFIER. Validated as test unit
**VARIANT-ANNOT-001** ([[variant-annot-001-evidence]]); see [[test-unit-registry]] for how units are
tracked and [[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The OverlapConsequence predicate system

VEP models each consequence as an **`OverlapConsequence`**: a boolean **predicate** plus a fixed
record of `SO_term`, `SO_accession`, `impact`, and `rank`. Crucially, **IMPACT is a stored property of
the term, not computed** — `frameshift_variant` is HIGH, `missense_variant` MODERATE, `synonymous_variant`
LOW, `intron_variant` MODIFIER, by table lookup (`Constants.pm`). The **rank** (1 = most severe) drives
**most-severe selection**: when several consequences apply across transcripts, the reported term is the
**lowest-rank** one — VEP's rule that "only the most specific term that applies under any given parent
term is assigned" (McLaren 2016).

| IMPACT | Representative SO terms (rank) |
|--------|------------------------------|
| **HIGH** | `transcript_ablation`(1), `splice_acceptor_variant`(2), `splice_donor_variant`(3), `stop_gained`(4), `frameshift_variant`(5), `stop_lost`(6), `start_lost`(7) |
| **MODERATE** | `inframe_insertion`(11), `inframe_deletion`(12), `missense_variant`(13), `protein_altering_variant`(14) |
| **LOW** | `splice_region_variant`(16), `stop_retained_variant`(21), `synonymous_variant`(22) |
| **MODIFIER** | `coding_sequence_variant`(23), `5_prime_UTR_variant`(25), `3_prime_UTR_variant`(26), `intron_variant`(28), `upstream_gene_variant`(32), `downstream_gene_variant`(33), `intergenic_variant`(40) |

## The coding-consequence engine translates codons

The coding-substitution consequences are decided by **translating the reference and alternate codons**
through the [[genetic-code-translation|standard genetic code]] (NCBI `transl_table` 1) and comparing the
resulting peptides `ref_pep` / `alt_pep`. This is why the concept **`depends_on`** the codon table —
change the table and the `*`-codon classification changes. The verbatim `VariationEffect.pm` predicates:

- **`synonymous_variant`** — `alt_pep == ref_pep`, not a retained stop, and **neither peptide contains `X`**
  (an ambiguous/untranslatable codon is never called synonymous).
- **`missense_variant`** — `ref_pep != alt_pep` **and** same length, **but false** if `start_lost`,
  `stop_lost`, `stop_gained`, or `partial_codon` (those override it).
- **`stop_gained`** — `alt_pep` gains a `*` the `ref_pep` lacks (premature stop).
- **`stop_lost`** — `ref_pep` has a `*`, `alt_pep` does not.
- **`frameshift_variant`** — purely **length-based**: `|alt_len − ref_len| mod 3 ≠ 0`, regardless of
  sequence content.
- **`inframe_insertion` / `inframe_deletion`** — indel length divisible by 3, alt codon longer / shorter
  than ref codon.

**Precedence:** the HIGH-rank terms dominate — `stop_gained` overrides `missense`; `start_lost`
(rank 7) overrides both `missense_variant` and `inframe_insertion`. A substitution at a reference stop
codon (TAA/TAG/TGA) that yields a sense codon is `stop_lost`.

**Worked oracles.** `GAA→GTA` (E→V) → missense/MODERATE · `TTA→TTG` (L→L) → synonymous/LOW ·
`CAA→TAA` (Q→*) → stop_gained/HIGH · `TAA→CAA` (*→Q) → stop_lost/HIGH · `ATG→ATC` (M→I at CDS start)
→ start_lost/HIGH · `AC→A` (Δ−1) → frameshift/HIGH · `A→ATTT` (Δ+3) → inframe_insertion/MODERATE ·
`ATTT→A` (Δ−3) → inframe_deletion/MODERATE.

## Where it sits in the variant-analysis family

Annotation is **downstream of variant calling** and orthogonal to pathogenicity/tier classification.
The upstream germline caller is [[germline-variant-calling-snp-indel]] (VARIANT-CALL-001 —
SNP/indel detection from a reference↔query alignment, with Ti/Tv); on the oncology side, calling is
[[somatic-variant-calling-tumor-normal]] (tumor-vs-normal VAF classification). This unit is the
germline-style VEP consequence layer that consumes the variants either caller produces.
It **reuses** the [[genetic-code-translation|genetic code]] table that also powers
[[open-reading-frame-detection|ORF detection]] and [[codon-optimization]] — here to translate the two
codons whose peptides the predicates compare.

## Scope and assumptions

- **Standard genetic code (table 1) only** — VEP uses transcript-specific tables; non-standard organism
  codes are out of scope (they would change `*`-codon classification). Human nuclear transcripts, VEP's
  dominant use, are table 1.
- **Single-codon SNV peptide comparison** — ref/alt peptides span only the codon(s) the variant overlaps;
  long-range effects (NMD, downstream re-initiation) are not modelled.

Research-grade, not for clinical use. Reference sources — **McLaren 2016** (VEP, *Genome Biology* 17:122),
Ensembl **`Constants.pm`** (consequence/IMPACT/rank table) + **`VariationEffect.pm`** (predicates), NCBI
**`gc.prt`** (Standard code), and the **Sequence Ontology** (Eilbeck 2005) — full record in
[[variant-annot-001-evidence]]. No source contradictions.
