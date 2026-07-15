---
type: source
title: "Evidence: VARIANT-INDEL-001 (indel detection — FindInsertions/FindDeletions + directional length invariant + normalization theory)"
tags: [validation, annotation]
doc_path: docs/Evidence/VARIANT-INDEL-001-Evidence.md
sources:
  - docs/Evidence/VARIANT-INDEL-001-Evidence.md
source_commit: 744b87f09d5f0949dd2f3eec7ec5a2c84e389606
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: VARIANT-INDEL-001

The validation-evidence artifact for test unit **VARIANT-INDEL-001** — **indel detection**: the
insertion/deletion facet of the alignment-column variant caller. `FindInsertions` / `FindDeletions`
are **filters over** `CallVariants` (the SNP/indel caller validated as VARIANT-CALL-001), returning
only the indel columns of a reference↔query global alignment. One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[germline-variant-calling-snp-indel]] (this unit **enriches** that page — it is the indel facet of the
same caller, not a separate algorithm). [[test-unit-registry]] tracks the unit.

## What this file records

- **Algorithm:** identify **insertions** and **deletions** between a reference and a query via the
  gap columns of `SequenceAligner.GlobalAlign`. In the per-column in-memory model: an **insertion** is a
  column with `ReferenceAllele = "-"` (gap sentinel) + a one-base `AlternateAllele`, `Type = Insertion`;
  a **deletion** is a column with a one-base `ReferenceAllele` + `AlternateAllele = "-"`, `Type =
  Deletion`. `FindInsertions` returns insertions only; `FindDeletions` returns deletions only. A
  multi-base indel is reported as consecutive per-base indel columns.

- **Directional length invariant** (the load-bearing correctness claim): **insertion ⇒ ALT longer than
  REF**; **deletion ⇒ REF longer than ALT**. In serialized VCF this is REF=`C`/ALT=`CA` (insertion) and
  REF=`TC`/ALT=`T` (deletion); in the in-memory model it manifests as which side carries the `"-"` gap
  sentinel.

- **Authoritative sources:**
  - **VCFv4.3 specification** (samtools/hts-specs, rank 2) — verbatim indel grammar: single-base
    insertion of A after pos 3 → REF=`C`, ALT=`CA`; single-base deletion of C at pos 3 → REF=`TC`,
    ALT=`T`; the **padding-base rule** (a pure insertion/deletion must carry the base *before* the event,
    or after at contig position 1); **POS is 1-based**; REF alphabet A,C,G,T,N. §1.1 microsatellite
    record `GTC → G,GTCT` — one 2-base deletion (TC), one 1-base insertion (T), both anchored at `G`.
  - **Tan A., Abecasis G.R., Kang H.M. (2015)** *Bioinformatics* 31(13):2202, DOI 10.1093/bioinformatics/btv112,
    PMID 25701572 (rank 1) — a VCF entry is **normalized ⟺ left-aligned + parsimonious**; **left-aligned**
    = smallest base position among equivalent entries; **parsimonious** = shortest allele length.
    **Algorithm 1** = right-trim shared nucleotides (re-pad left when an allele empties), then left-trim
    shared leading nucleotides while all alleles length ≥2. The same indel has multiple representations;
    an alignment-column indel is **not** guaranteed canonical → bounds the *position* claim in
    repeated/low-complexity regions only (not counts/types).
  - **minimal_representation** (E. Minikel, `normalize.py`, rank 3) — reference implementation of Tan 2015
    Algorithm 1: **suffix-then-prefix trimming** (`while ref[-1]==alt[-1]` re-padding left + decrementing
    pos, then `while len>1 and ref[0]==alt[0]` incrementing pos). Worked cases: CFTR p.F508del
    `(7,117199646,CTT,-) → (7,117199644,ATCT,A)` (3-base deletion); BRCA2 `(13,32914438,T,-) →
    (13,32914437,GT,G)` (1-base deletion) — both confirm deletion ⇒ len(REF)>len(ALT) and left-anchor
    padding of empty alleles.
  - **PharmCAT — Variant Normalization** (rank 3) — worked left-alignment in a tandem repeat: `POS
    97740414 AATGA→A` shifts to `POS 97740410 GATGA→G`; count (one deletion) and type unchanged, only
    position/spelling. Corroborates Tan 2015: indel position is alignment/normalization dependent in repeats.

## Oracles and datasets (independent of the library code)

- **VCFv4.3 indel examples:** insertion of A after pos 3 → REF=C, ALT=CA (ALT>REF) · deletion of C at
  pos 3 → REF=TC, ALT=T (REF>ALT) · microsat deletion 2 bp `GTC→G` · microsat insertion 1 bp `GTC→GTCT`.
- **Alignment-derived indel columns (unique alignment, no repeat):** `ATGCAT`→`ATGTCAT` = 1-base
  insertion of `T` (Insertion, `-`, `T`); `ATGTCAT`→`ATGCAT` = 1-base deletion of `T` (Deletion, `T`,
  `-`); `ATGC`→`ATGC` = no event.
- **minimal_representation trimming:** `(7,117199646,CTT,-)→(7,117199644,ATCT,A)`;
  `(13,32914438,T,-)→(13,32914437,GT,G)` — both verify deletion ⇒ len(REF)>len(ALT).

## Corner cases and failure modes

- **Empty allele forbidden in serialized VCF; padding base required** — a pure insertion (REF empty) or
  deletion (ALT empty) carries a left anchor base in `ToVcfLines` output; the in-memory caller uses the
  `"-"` gap sentinel per column instead (ASM-01, representation shape only).
- **Non-unique indel placement** — in a repeated/low-complexity region the alignment producing the indel
  column is not unique, so reported position is implementation-dependent unless normalized (ASM-02, Tan
  2015). Affects **position** only, not counts or types.
- **Filters return their class only** — `FindInsertions` yields no deletions/SNPs; `FindDeletions` yields
  no insertions/SNPs. Identical sequences → no indels. A substitution-only input → no indels. Null
  reference/query → `ArgumentNullException` (propagated from `CallVariants`).

## Assumptions

1. **Internal gap-sentinel representation for indels (ASM-01).** Indel columns use `"-"` for the absent
   allele and a **0-based** `Position`, not the VCF padded-allele/1-based form (which `ToVcfLines`
   produces, out of scope). Representation shape, not a detection decision — counts and types are
   unchanged. Consistent with the sibling VARIANT-CALL-001 contract.
2. **Indels not left-aligned / parsimony-normalized (ASM-02).** The caller reports the indel at the
   `GlobalAlign` column with no normalization pass. Correctness-affecting for **position** in repeats
   only; tests assert counts/types/alleles generally and exact position only where the alignment is
   provably unique.
3. **0-based in-memory `Position`.** 1-based POS is produced only by `ToVcfLines`; matches the sibling
   contract. Not a source-governed value for the in-memory model.

**No source contradictions** — VCFv4.3, Tan 2015, minimal_representation, and PharmCAT agree on the
directional length invariant, the padding-base rule, and that indel position is normalization-dependent
in repeats. The two representation gaps (gap-sentinel vs VCF padding; un-normalized position) are
documented implementation choices bounded to positional (not count/type) correctness. Research-grade,
not for clinical use.
