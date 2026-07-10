---
type: source
title: "Evidence: VARIANT-SNP-001 (SNP detection — FindSnps / FindSnpsDirect + Hamming-mismatch enumeration + Ti/Tv)"
tags: [validation, annotation]
doc_path: docs/Evidence/VARIANT-SNP-001-Evidence.md
sources:
  - docs/Evidence/VARIANT-SNP-001-Evidence.md
source_commit: b0f0fbc85f38f574f0ce5c1ffef41c0bbd45dddb
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: VARIANT-SNP-001

The validation-evidence artifact for test unit **VARIANT-SNP-001** — **SNP detection**: the
single-nucleotide-substitution facet of the reference↔query variant caller. `FindSnps` (alignment-based)
and `FindSnpsDirect` (positional / Hamming-style) return **only the substitution columns** — SNPs, no
indels. One instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the synthesizing concept is [[germline-variant-calling-snp-indel]] (this unit **enriches** that
page — it is the SNP facet of the same caller, not a separate algorithm, mirroring how
[[variant-indel-001-evidence]] folded in the indel facet). [[test-unit-registry]] tracks the unit.

## What this file records

- **Algorithm:** identify **single-nucleotide substitutions** between a reference and a query. Two entry
  points:
  - **`FindSnps`** — alignment-based; runs the same reference↔query global alignment as the caller and
    reports the SNP columns only (both bases present and differing; `Type = SNP`), never insertions or
    deletions.
  - **`FindSnpsDirect`** — positional / **Hamming-style**; compares the two sequences index-by-index and
    reports each mismatched position as one SNP. Over two **equal-length** sequences this is exactly the
    enumeration of the **Hamming-mismatch positions** (each mismatched index = one substitution). No gap
    handling — substitutions only.
  - In-memory per SNP: `Position = i` (**0-based** mismatched index), `ReferenceAllele = ref[i]`,
    `AlternateAllele = query[i]`, with `ref[i] ≠ query[i]`.

- **Transition / transversion classification** (case-insensitive): each SNP is a **transition** (A↔G,
  C↔T — within a ring class) or a **transversion** (A↔C, A↔T, G↔C, G↔T — across ring classes); **4
  possible transversions vs 1 transition per base**, yet transitions occur more often (transitional bias).
  Ti/Tv is defined only for substitutions.

- **Authoritative sources:**
  - **VCFv4.3 specification** (samtools/hts-specs, rank 2) — §1.1 canonical **simple-SNP** record
    `20 14370 rs6054257 G A …` = a single-base REF `G` substituted by single-base ALT `A`, no padding
    base; **POS is 1-based** ("1st base has position 1"); **REF alphabet A,C,G,T,N, case-insensitive** ⇒
    base comparison and Ti/Tv classification are case-insensitive; a position where **REF == ALT is not a
    variant** (a SNP is a *substitution*).
  - **Wikipedia — Transversion** (rank 4; primary Futuyma, *Evolution*, 2013) — a transversion is a point
    mutation swapping a purine (A/G) for a pyrimidine (T/C) or vice versa; the four transversions are
    A↔C, A↔T, G↔C, G↔T; four transversions but one transition possible per base, transitions still more
    frequent (Ti/Tv bias).
  - **Hamming Distance as a Concept in DNA Molecular Recognition** (Acharya 2017, *ACS Omega*, PMC5410656,
    rank 1; cites Bierbrauer, *Introduction to Coding Theory*, 2016) — "the number of positions that two
    codewords of the same length differ is the Hamming distance." ⇒ `FindSnpsDirect` over two same-length
    sequences enumerates exactly the Hamming-mismatch positions; the distance is **defined only for equal
    lengths**.
  - **Collins DW, Jukes TH (1994)** *Genomics* 20(3):386, DOI 10.1006/geno.1994.1192, PMID 8034311 (rank
    1) — transitional silent-substitution rate 1.71×10⁻⁹ yr⁻¹ > transversional 1.22×10⁻⁹ yr⁻¹; grounds
    that Ti vs Tv is a real asymmetric classification (Ti/Tv > 0.5), not any expected numeric test value.

## Oracles and datasets (independent of the library code)

- **VCFv4.3 §1.1 simple SNP:** POS 14370, REF=`G`, ALT=`A` — single-base substitution G→A. In the 0-based
  per-position model: `Position = i`, `ReferenceAllele = ref[i]`, `AlternateAllele = query[i]`.
- **Positional substitution (Hamming mismatches):** `ATGC`→`ATTC` = {2} → G→T · `AAAA`→`TGTA` = {0,1,2}
  → A→T, A→G, A→T · `ATGC`→`ATGC` = {} (Hamming distance 0, no SNP).
- **Ti/Tv classification of alleles:** A→G Transition (purine↔purine) · C→T Transition (pyrimidine↔
  pyrimidine) · A→C Transversion (purine↔pyrimidine) · G→T Transversion (purine↔pyrimidine).

## Corner cases and failure modes

- **REF == ALT is not a variant** — a column/position whose query base equals the reference base is a
  match; the detector skips equal columns (VCFv4.3: a SNP is a substitution).
- **Case insensitivity** — `a`/`g` classify identically to `A`/`G` (VCFv4.3 REF alphabet is
  case-insensitive).
- **Equal-length precondition for `FindSnpsDirect`** — the Hamming distance is defined only for equal
  lengths (PMC5410656). When inputs differ in length, positional comparison is well-defined only over the
  **common prefix** (`min(reference.Length, query.Length)`); the trailing region of the longer sequence
  is **indel territory** (VARIANT-INDEL-001), not a SNP — out of scope here.
- **Input validation** — null inputs to `FindSnps` throw `ArgumentNullException`; empty inputs to
  `FindSnpsDirect` yield an empty set. `FindSnps` on a substitution-only input reports SNPs only, no
  indels. Every emitted variant is `Type = SNP` with distinct ref/alt alleles.

## Assumptions

1. **Unequal-length `FindSnpsDirect` compares only the common prefix (ASM-01).** The Hamming distance is
   defined for equal-length strings only (PMC5410656); the contract iterates `min(reference.Length,
   query.Length)` and reports substitutions over the common prefix. The trailing region of the longer
   input is not reported (it is indel territory, VARIANT-INDEL-001). No source mandates substitution
   semantics beyond the common length — this prefix behavior is the documented contract, not a defect;
   tests assert it explicitly.
2. **0-based in-memory `Position` (ASM-02).** VCF serialization POS is 1-based (spec field 4), but the
   in-memory `Variant.Position` is **0-based**; the 1-based POS is produced only by `ToVcfLines` (out of
   scope). Matches the sibling VARIANT-CALL-001 / VARIANT-INDEL-001 contract; internally consistent, not a
   source-governed value for the in-memory model.

**No source contradictions** — VCFv4.3 (SNP = single-base substitution, case-insensitive alphabet),
Wikipedia/Futuyma (the four transversions), PMC5410656 (Hamming mismatch = substitution enumeration, equal
lengths only), and Collins & Jukes 1994 (transitional bias) agree. The two representation choices
(common-prefix behavior on unequal lengths; 0-based in-memory position) are documented contracts bounded
to scope/representation, not counts or types. Research-grade, not for clinical use.
