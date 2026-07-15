---
type: source
title: "Evidence: VARIANT-CALL-001 (variant calling — SNP/indel detection from reference↔query alignment + Ti/Tv)"
tags: [validation, annotation]
doc_path: docs/Evidence/VARIANT-CALL-001-Evidence.md
sources:
  - docs/Evidence/VARIANT-CALL-001-Evidence.md
source_commit: 5b4dd805db54d51bae30445a884e122fc4d97bd5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: VARIANT-CALL-001

The validation-evidence artifact for test unit **VARIANT-CALL-001** — **variant calling**: detecting
**SNPs / insertions / deletions** by comparing a query sequence against a reference and classifying each
substitution as a **transition or transversion** (with the **Ti/Tv** ratio). One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing
concept is [[germline-variant-calling-snp-indel]]. [[test-unit-registry]] tracks the unit.

This is the **calling / detection** layer of the variant-analysis family — it *produces* the variants
that downstream units interpret. It sits **upstream of** variant **annotation**
([[variant-effect-annotation-vep]], VARIANT-ANNOT-001), and is the **germline, reference↔query**
counterpart to the oncology **tumor-vs-normal** caller [[somatic-variant-calling-tumor-normal]]
(ONCO-SOMATIC-001). Method: `CallVariantsFromAlignment` runs `SequenceAligner.GlobalAlign` on
reference vs query, then walks the aligned columns emitting one `Variant` per differing column.

## What this file records

- **Algorithm:** alignment-column variant detection. Three variant classes (Danecek 2011): **SNP**
  (single-base substitution), **Insertion**, **Deletion**. Substitutions are further classified
  **Transition** (purine↔purine A↔G or pyrimidine↔pyrimidine C↔T) vs **Transversion** (purine↔pyrimidine:
  A↔C, A↔T, G↔C, G↔T), case-insensitively, and reduced to the **Ti/Tv ratio**.

- **Authoritative sources:**
  - **VCFv4.3 specification** (samtools/hts-specs, 2024, rank 2) — the field grammar for serialized
    output (`ToVcfLines`): **POS is 1-based** ("1st base has position 1"); REF/ALT bases ∈
    A,C,G,T,N case-insensitive, `ALT` matches `^([ACGTNacgtn]+|\*|\.)$`; **indel padding base** — a pure
    insertion/deletion must carry the base *before* the event (or after, at contig position 1). Canonical
    §1.1 examples: `G→A` simple SNP; `GTC→G` 2-base deletion (TC); `GTC→GTCT` 1-base insertion (T), both
    anchored at the preceding `G`.
  - **Danecek P. et al. (2011)** *Bioinformatics* 27(15):2156, PMID 21653522 (rank 1) — the VCF reference
    paper: "VCF is a generic format for storing DNA polymorphism data such as SNPs, insertions, deletions
    and structural variants" (confirms the three classes; a variant is a *difference from reference*).
  - **Tan A., Abecasis G.R., Kang H.M. (2015)** *Bioinformatics* 31(13):2202 (rank 1) — a VCF entry is
    **normalized ⟺ left-aligned + parsimonious**; the same biological indel has multiple representations.
    Bounds the correctness claim: an alignment-derived indel is **not** guaranteed canonical.
  - **Collins D.W., Jukes T.H. (1994)** *Genomics* 20(3):386, PMID 8034311 (rank 1) — transitional bias:
    transition rate (1.71×10⁻⁹) > transversion rate (1.22×10⁻⁹), the basis for **Ti/Tv > 0.5** expectations;
    "≈2/3 of SNPs are transitions."
  - **Wikipedia Transition (genetics) / Transversion** (rank 4, citing Futuyma 2013) — the enumerated
    transition/transversion classification table.

## Oracles and datasets (independent of the library code)

- **VCFv4.3 §1.1 records:** SNP `G→A` @14370 (a transition) · deletion `GTC→G` @1234567 (2 bases TC) ·
  insertion `GTC→GTCT` @1234567 (1 base T).
- **Ti/Tv classification:** A→G, G→A, C→T, T→C → **Transition**; A→C, A→T, G→C, G→T, C→A, T→A →
  **Transversion**.
- **Ti/Tv ratio:** `{A→G, A→C}` → Ti 1 / Tv 1 → **1.0**; `{A→G, C→T}` → Ti 2 / Tv 0 → **0** by repo
  convention (undefined ratio mapped to 0).

## Corner cases and failure modes

- **Identical sequences → zero variants** (a variant is a difference from reference — Danecek 2011).
- **Empty allele forbidden in VCF** — a pure insertion/deletion carries a padding base in serialized VCF;
  the in-memory caller instead uses a `"-"` gap sentinel per column (internal representation, ASM-01).
- **Non-unique indel placement** — in low-complexity / repeated regions the alignment producing the indel
  column is not unique, so the *reported position* is implementation-dependent unless normalized (Tan 2015).
- **Case-insensitive** REF/ALT and Ti/Tv classification.
- **Mismatched aligned lengths → `ArgumentException`**; empty input → empty; non-SNP variant classifies
  as `Other` (Ti/Tv defined only for SNPs).

## Assumptions

1. **Internal gap-sentinel indel representation.** `CallVariantsFromAlignment` reports indels with the
   `"-"` gap character and a **0-based** `Position`, not the VCF padded-allele **1-based** form. The
   padding/1-based rule governs only *serialized VCF* (`ToVcfLines`, out of scope here); the in-memory
   `Variant` model is an implementation choice, internally consistent with sibling methods. Documented, not
   changed.
2. **Indels not left-aligned / parsimony-normalized.** Per Tan 2015 the canonical form needs
   left-alignment + parsimony; the caller reports the indel at the column `GlobalAlign` produces, with no
   normalization pass. Affects **position** in repeated regions only, not variant **counts/types**. Tests
   assert counts/types/alleles on unambiguous inputs and position only where the alignment is unique.
3. **Ti/Tv with zero transversions returns 0.** The mathematically-undefined `#Tv = 0` case is mapped to
   **0** (existing contract) rather than throwing or `+∞`; no source mandates a sentinel.

**No source contradictions** — VCFv4.3, Danecek 2011, Tan 2015, and Collins & Jukes 1994 agree on the
variant classes, the padding/normalization conventions, and the transition/transversion bias. The two
representation gaps (gap-sentinel vs VCF padding; un-normalized indel position) are flagged as documented
implementation choices bounded to positional (not count/type) correctness. Research-grade, not for clinical
use.
