---
type: source
title: "Evidence: VARIANT-ANNOT-001 (variant annotation — VEP-style functional consequence / Sequence-Ontology impact)"
tags: [validation, annotation]
doc_path: docs/Evidence/VARIANT-ANNOT-001-Evidence.md
sources:
  - docs/Evidence/VARIANT-ANNOT-001-Evidence.md
source_commit: dfa368702c4b3153cdb3b2bb48877ac663f4e019
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: VARIANT-ANNOT-001

The validation-evidence artifact for test unit **VARIANT-ANNOT-001** — **variant annotation**: mapping
a sequence variant to its **functional consequence** (missense / synonymous / stop_gained / stop_lost /
start_lost / frameshift / inframe indel / …) and the associated **Sequence-Ontology (SO) term + IMPACT
class** (HIGH / MODERATE / LOW / MODIFIER), reproducing the **Ensembl Variant Effect Predictor (VEP)**.
One instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the synthesizing concept is [[variant-effect-annotation-vep]]. [[test-unit-registry]] tracks the unit.

This is the **interpretation / annotation** layer of the variant-analysis family — it takes an already-
*called* variant and predicts *what it does*, sitting downstream of variant **calling**
([[somatic-variant-calling-tumor-normal]] on the oncology side) and orthogonal to
pathogenicity/tier classification. Its consequence engine **translates the reference and alternate
codons** through the [[genetic-code-translation|standard genetic code]] (NCBI `transl_table` 1) and
compares the resulting peptides.

## What this file records

- **Algorithm:** the VEP **OverlapConsequence predicate system** — each consequence is a boolean
  predicate over the reference/alternate peptide (or allele length), and every term carries a fixed
  `SO_term`, `SO_accession`, `impact`, and `rank`. Two source files pin it down: `Constants.pm`
  (the consequence → IMPACT → rank table) and `VariationEffect.pm` (the predicate logic).

- **Authoritative sources:**
  - **McLaren W. et al. (2016)** *Genome Biology* 17:122, PMC4893825 (rank 1, primary) — VEP; SO terms
    defined with the Sequence Ontology; **hierarchy = "only the most specific term that applies under
    any given parent term is assigned"**; IMPACT classification over the SO consequence types.
  - **Ensembl `Constants.pm`** (ensembl-variation release/110, rank 3, reference impl) — the ordered
    `OverlapConsequence` table (rank, SO_term, IMPACT, SO_accession); **IMPACT is a stored property of
    the term, not computed separately**.
  - **Ensembl `VariationEffect.pm`** (rank 3, reference impl) — the verbatim consequence predicates.
  - **NCBI Genetic Codes `gc.prt`**, Standard code table 1 (rank 2) — the codon→AA table used to
    translate ref/alt codons (drives synonymous/missense/stop_gained/stop_lost).
  - **Eilbeck K. et al. (2005)** the Sequence Ontology (SO accessions).

- **The consequence → IMPACT → rank table (verbatim from `Constants.pm`, abridged):**
  HIGH: `transcript_ablation`(1), `splice_acceptor_variant`(2), `splice_donor_variant`(3),
  `stop_gained`(4), `frameshift_variant`(5), `stop_lost`(6), `start_lost`(7); MODERATE:
  `inframe_insertion`(11), `inframe_deletion`(12), `missense_variant`(13), `protein_altering_variant`(14);
  LOW: `splice_region_variant`(16), `stop_retained_variant`(21), `synonymous_variant`(22); MODIFIER:
  `coding_sequence_variant`(23), `5_prime_UTR_variant`(25), `3_prime_UTR_variant`(26), `intron_variant`(28),
  `upstream_gene_variant`(32), `downstream_gene_variant`(33), `intergenic_variant`(40). Lower rank = more severe.

- **The peptide predicates (verbatim from `VariationEffect.pm`):**
  - `synonymous_variant`: `alt_pep == ref_pep` **and** not stop_retained **and** neither peptide contains `X`.
  - `missense_variant`: false if start_lost/stop_lost/stop_gained/partial_codon; else `ref_pep != alt_pep` **and** `len(ref_pep) == len(alt_pep)`.
  - `stop_gained`: `alt_pep` contains `*` that `ref_pep` does not.
  - `stop_lost`: `ref_pep` contains `*`, `alt_pep` does not.
  - `frameshift`: `abs(allele_len − var_len) % 3 != 0` (indel length **not** a multiple of 3).
  - `inframe_insertion` / `inframe_deletion`: `len(alt_codon) > / < len(ref_codon)` with the indel length divisible by 3.

## Oracles and datasets (all independent of the library code)

- **Standard-code codon substitutions (NCBI table 1 + VEP predicates):**
  `GAA→GTA` (E→V) → **missense_variant** MODERATE; `TTA→TTG` (L→L) → **synonymous_variant** LOW;
  `CAA→TAA` (Q→*) → **stop_gained** HIGH; `TAA→CAA` (*→Q) → **stop_lost** HIGH;
  `ATG→ATC` (M→I at CDS start) → **start_lost** HIGH.
- **Coding indels (length rule):** `AC→A` (Δ−1) → **frameshift_variant** HIGH;
  `A→ATTT` (Δ+3) → **inframe_insertion** MODERATE; `ATTT→A` (Δ−3) → **inframe_deletion** MODERATE.

## Precedence and corner cases (from `VariationEffect.pm` + `gc.prt`)

- **Most-severe selection:** across transcripts/consequences the reported term is the **lowest-rank**
  (most severe) one — VEP's "most specific term under any parent" plus the rank ordering.
- **stop_gained overrides missense** — a substitution introducing a premature stop is stop_gained, never missense.
- **start_lost overrides coding substitution** — `missense_variant`/`inframe_insertion` return 0 when start_lost (a HIGH rank-7 term wins over the MODERATE substitution terms).
- **Frameshift is purely length-based** — any coding indel with `|alt−ref| mod 3 ≠ 0` is frameshift regardless of sequence content.
- **Ambiguous codon (`X`) excludes synonymous** — an untranslatable/ambiguous codon is never reported synonymous.
- **Stop in reference codon** — a substitution at a reference stop (TAA/TAG/TGA) yielding a sense codon is stop_lost.

## Assumptions

1. **Standard genetic code (table 1) only.** VEP uses transcript-specific codon tables; this unit
   translates with the NCBI Standard code. Non-standard organism codes would change `*`-codon
   classification (out of scope; human nuclear transcripts — the dominant VEP use — use table 1).
2. **Single-codon SNV peptide comparison.** Ref/alt peptides are compared over the codon(s) directly
   overlapped by the variant; long-range effects (NMD, downstream re-initiation) are not modelled.

**No source contradictions** — McLaren 2016, `Constants.pm`, `VariationEffect.pm`, and NCBI `gc.prt`
agree on the SO term set, the IMPACT/rank mapping, the peptide predicates, and the most-specific-term
hierarchy. Research-grade, not for clinical use.
