---
type: concept
title: "Codon Usage Statistics (aggregation view + positional GC)"
tags: [annotation, algorithm]
sources:
  - docs/algorithms/Codon/Codon_Usage_Statistics.md
  - docs/Evidence/CODON-STATS-001-Evidence.md
  - docs/Validation/reports/CODON-STATS-001.md
source_commit: 06fc6f57e4f9c4947db8a2c97df2b931367b92b6
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:relative-synonymous-codon-usage
      source: codon-usage-statistics
      evidence: "GetStatistics bundles the RSCU table (RSCU_j = n·x_j/Σx_k, Sharp/Tuohy/Mosurski 1986) as one field of the aggregated CodonUsageStatistics record"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:effective-number-of-codons
      source: codon-usage-statistics
      evidence: "GetStatistics reports Enc (Wright 1990, 20–61) alongside the codon counts and positional GC as one aggregated statistic"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:codon-adaptation-index
      source: codon-usage-statistics
      evidence: "CodonUsageAnalyzer.CalculateCai(sequence, referenceRscu) is the same geometric-mean-of-w=f/max_synonym_f CAI, re-validated as one piece of CODON-STATS-001"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-usage-statistics
      evidence: "Test Unit ID: CODON-STATS-001 — Codon Usage Statistics (CodonUsageAnalyzer.GetStatistics + CalculateCai)"
      confidence: high
      status: current
---

# Codon Usage Statistics (aggregation view + positional GC)

The **aggregation / reporting view** of the codon-usage-bias family: a single call
(`CodonUsageAnalyzer.GetStatistics`) that bundles every standard descriptor of synonymous codon
bias for one coding DNA sequence into one `CodonUsageStatistics` record, plus a separate
`CalculateCai` against a reference set. It is the umbrella unit **CODON-STATS-001**
(`Seqeron.Genomics.MolTools`); the same unit's validation artifacts are already ingested as the
source pages [[codon-stats-001-evidence|CODON-STATS-001 evidence]] and
[[codon-stats-001-report|CODON-STATS-001 report]] (two-stage: Stage A **PASS-WITH-NOTES** / Stage B
**PASS** / End state **✅ CLEAN**). See [[test-unit-registry]] for how the unit is tracked and
[[algorithm-validation-evidence]] for the shared evidence pattern.

Most of what it computes is owned by dedicated family concepts — this page synthesizes the **two
things unique to the aggregation view**: (1) the *bundle contract* (what one pass returns), and
(2) the **positional GC block** (GC1/GC2/GC3, GC3s, OverallGc) that no other codon concept covers.

## What the bundle returns

`GetStatistics(string | DnaSequence)` reads the sequence in **frame 1, non-overlapping triplets
from offset 0** (trailing < 3 nt ignored; case-insensitive; any codon containing a non-`{A,C,G,T}`
character is **skipped**, not an error), and returns:

| Field | Owned/defined by | Note |
|-------|------------------|------|
| `CodonCounts` | (this unit) | per-codon in-frame occurrence counts |
| `Rscu` | [[relative-synonymous-codon-usage\|RSCU]] | `RSCU_j = n·x_j/Σ_k x_k` (Sharp et al. 1986) |
| `Enc` | [[effective-number-of-codons\|ENC / Nc]] | Wright 1990, bounded [20, 61] |
| `TotalCodons` | (this unit) | number of valid in-frame ACGT codons (INV-05) |
| `Gc1`, `Gc2`, `Gc3` | (this unit — positional GC) | % G/C at codon position 1 / 2 / 3 |
| `Gc3s` | (this unit — positional GC) | % G/C at **synonymous** third positions |
| `OverallGc` | (this unit) | `(Gc1+Gc2+Gc3)/3` (record-derived, INV-06) |

`CalculateCai(sequence, referenceRscu)` is the same **[[codon-adaptation-index|CAI]]** geometric
mean of relative adaptiveness `w_i = f_i / max_j f_j` (Sharp & Li 1987), returned separately in
`[0, 1]`. Contract corners: a `null` `DnaSequence` **or** `null` reference table throws
`ArgumentNullException`; a `null`/empty *string* returns a **zeroed** record (CAI 0).

## Positional GC content (the distinctive block)

These four descriptors are the part of the family this unit alone provides — the GC composition
resolved *by codon position*, a standard codon-bias axis (base composition drives synonymous bias
independently of the RSCU/CAI/ENC measures).

- **GC1 / GC2 / GC3** — the fraction (×100, so 0–100) of in-frame codons carrying **G or C at codon
  position 1 / 2 / 3** respectively ("1st/2nd/3rd letter GC", EMBOSS `cusp` convention). INV-04:
  each ∈ [0, 100].
- **GC3s** — "the frequency of G or C at the third position of **synonymous** codons", i.e.
  counting **only** codons whose amino acid is degenerate (degeneracy > 1). It therefore
  **excludes Met (ATG), Trp (TGG) and the stop codons** (Peden 1999, CodonW §1.8.2.1.3; INV-03).
  This is what separates GC3s from GC3: GC3 counts all third positions, GC3s only the ones that can
  vary synonymously.
- **OverallGc** = `(GC1+GC2+GC3)/3` — a record-derived mean, *not* the whole-sequence GC%.

**Worked contrast (`ATGGCA`):** Met (ATG) is excluded from GC3s; the only synonymous codon left is
Ala (GCA), whose third base `A` is not G/C ⇒ **GC3s = 0/1 = 0%**. Over *all* third positions the
`G` of ATG plus the `A` of GCA give **GC3 = 1/2 = 50%**. The gap between them is exactly the
Met/Trp/stop exclusion of Peden §1.8.2.1.3.

**Deviation (display units).** GC3s is reported as a **percentage** (0–100), i.e. **100×** the
CodonW *fraction*; the synonymous subset counted is identical. Documented, accepted
(CODON-STATS-001 §5.4).

## CAI inside the aggregation (one re-validated deviation)

The bundled `CalculateCai` re-validates the standalone [[codon-adaptation-index|CAI]] unit, and
carries CODON-STATS-001's only standalone-CAI note: **zero-frequency codons are *skipped*** (not
scored) rather than floored — where the dedicated CAI unit clamps to `1e-6` and Bulmer 1988 uses a
`0.01` floor. Consequence: a gene using a codon entirely **absent** from the reference yields a
slightly higher CAI than EMBOSS/seqinr would. **No effect with the two bundled reference tables** —
neither `EColiOptimalCodons` (Sharp & Li 1987 *w* values, transcribed from Biopython
`SharpEcoliIndex`) nor `HumanOptimalCodons` (Kazusa *H. sapiens* [gbpri] RSCU, 93,487 CDS) contains
a synonymous *w* of 0. Because CAI rescales each reference value by its family maximum, passing the
*w* table (max 1.0) and passing the RSCU table are **equivalent** (the family-max cancels).

## Place in the codon-usage family

This is the **aggregation node** of the family: it *reports* the measures the dedicated concepts
*define*. **[[relative-synonymous-codon-usage|RSCU]]** is the per-codon normalization at the base;
**[[codon-adaptation-index|CAI]]** and **[[effective-number-of-codons|ENC / Nc]]** are the two
single-number whole-gene summaries (reference-based vs reference-free); **[[codon-optimization]]**
is the family's *rewriting* actuator; **[[rare-codon-analysis]]** localizes the low-`w` codons;
**[[codon-usage-comparison]]** is the *raw* count-table + distribution-distance end. What is unique
to CODON-STATS-001 is the **bundle-in-one-pass contract** plus the **positional GC block
(GC1/GC2/GC3/GC3s)** that none of the bias measures expose.

## Scope / limitations

DNA alphabet only (no IUPAC ambiguity, no T↔U conversion, no RNA). **Frame 1 only.** No per-amino-
acid "optimal codon" lists, Fop, CBI, or GC skew (use CodonW/EMBOSS or future units). CAI quality
depends on the supplied reference (ASM-01): the bundled human table is whole-genome Kazusa RSCU, not
a curated highly-expressed set, so its absolute CAI values are **descriptive**, not
expression-predictive.
