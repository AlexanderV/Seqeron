---
type: concept
title: "Tumor ploidy estimation (length-weighted mean segment CN) + whole-genome-doubling detection"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-PLOIDY-001-Evidence.md
source_commit: 57c2be1ccf184e08702ece85a7ad5a5d5618388c
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-ploidy-001-evidence
      evidence: "Test Unit ID: ONCO-PLOIDY-001 ... Algorithm: Tumor Ploidy Estimation (length-weighted mean segment copy number) and Whole-Genome-Doubling detection"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:allele-specific-copy-number-ascat
      source: onco-ploidy-001-evidence
      evidence: "Van Loo et al. PNAS 2010 (ASCAT) is cited as the n-scale ploidy source; ASCAT outputs a final tumour `ploidy` field, whereas this unit computes the same average ploidy as a post-hoc length-weighted mean of already-called allele-specific segment total CN."
      confidence: high
      status: current
---

# Tumor ploidy estimation (length-weighted mean segment CN) + whole-genome-doubling detection

Two scalar/binary **genome-state summaries** computed **post-hoc from already-called allele-specific
copy-number segments** (the shared `AlleleSpecificSegment` record — Major CN, Minor CN, chr, start,
end): the tumor's **average ploidy ψ** and a per-sample **whole-genome-doubling (WGD)** flag. Validated
under test unit **ONCO-PLOIDY-001**; the literature-traced record is [[onco-ploidy-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## Distinct from the ASCAT joint fit

This unit is **not** the [[allele-specific-copy-number-ascat|ASCAT joint purity/ploidy grid fit]]
(ONCO-ASCAT-001). ASCAT *fits* ρ and ψ **jointly from raw per-locus logR/BAF** by minimising a
length-weighted integer-closeness objective over a grid — an inference that produces the segments in
the first place. This unit runs **downstream** of that (or any) segmentation: given the resulting
allele-specific segments, it computes ψ as a **deterministic closed-form weighted mean** and applies a
**rule-based** WGD classifier. Same n-scale scalar (Van Loo 2010's ASCAT `ploidy` field), different
method — a post-hoc summary, not a grid search. It reads the same segment substrate as the other
genomic-scar oncology units [[loss-of-heterozygosity-detection]] and
[[homologous-recombination-deficiency-score]].

## 1. Average ploidy ψ — length-weighted mean total CN (`EstimatePloidy`)

Per **Patchwork** (Genome Biology 2010, PMC4053982), average tumour ploidy is "the average total copy
number of all genomic segments weighted by segment length":

```
ψ = Σ (CN_i · L_i) / Σ L_i          # CN_i = total copy number = Major + Minor ; L_i = End − Start
```

The averaged per-segment quantity is the **total** copy number (Major + Minor), not an allele-specific
value. ψ is reported on the **n-scale** (Van Loo 2010: **2n = diploid**; a pure copy-number-2 genome →
ψ = 2.0 exactly), and **elevated ploidy > 2.7n marks aneuploidy** (near-triploid basal-like genomes).

**Worked oracle.** Segments CN 2 / 4 / 3 at lengths 100 / 100 / 50 Mb → Σ(CN·L) = 750,000,000,
Σ(L) = 250,000,000 → **ψ = 3.0**. It is genuinely length-*weighted*: a long CN-2 segment plus a short
CN-4 segment weights toward 2, not toward the plain segment mean of 3.

## 2. Whole-genome doubling — the facets-suite / Bielski rule (`DetectWholeGenomeDoubling`)

WGD is a **per-sample binary genome-state** classification (Bielski et al., Nature Genetics 2018,
PMID 30013179 — WGD in ~30% of 9,692 advanced cancers), operationalised by the **facets-suite
`is_genome_doubled`** reference rule: WGD is called when the **autosome-restricted fraction of the
genome at major copy number ≥ 2 is strictly greater than 0.5**.

```
frac_elevated_mcn = Σ length[ mcn ≥ 2  AND  chrom ∈ 1..22 ]  /  autosomal_genome
wgd               = frac_elevated_mcn > 0.5                    # STRICT >, exactly-half is NOT doubled
mcn               = tcn − lcn                                  # major CN = total − minor (the larger allele)
autosomal_genome  = Σ chrom_info$size[ chr ∈ 1..22 ]          # a REFERENCE chromosome-size table
```

Three load-bearing details:

- **Major, not total, CN.** The elevated-fraction numerator uses **major** allele CN ≥ 2, so a balanced
  1:1 segment (total 2) has major = 1 and is **NOT** elevated, whereas a 2:0 (copy-neutral LOH) or 2:1
  segment **IS**. A genome entirely of 1:1 segments is never doubled.
- **Reference-genome denominator, autosomes only.** The denominator is a **reference chromosome-size
  table** (`ReferenceGenome { GRCh38, GRCh37 }`, default GRCh38), **not** the interrogated segments —
  so a small fully-amplified region cannot bias the call (100 Mb / 2.875 Gb ≈ 0.035 → false). Only
  autosomal (chr1–22) segments contribute; sex chromosomes / contigs are excluded, so a single-X male
  is not biased. The embedded tables equal the UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes` values
  (Ensembl GRCh38.p14 cross-verified): **Σ(chr1–22) = 2,875,001,522 bp (GRCh38)**,
  **2,881,033,286 bp (GRCh37)**.
- **Strict `> 0.5`.** Exactly half the autosomal genome at MCN ≥ 2 is **not** doubled — it must be
  strictly more than half. GRCh38 half = 1,437,500,761 bp: half+1 → true, exactly half → false,
  half−1 → false.

**Legacy overload.** `DetectWholeGenomeDoublingFromSuppliedLength` retains the pre-fix semantics —
denominator = Σ supplied segment length (60% at MCN ≥ 2 → true, exactly 50% → false, empty → throws) —
kept for back-compatibility; the reference-genome denominator is the corrected default.

## Corner cases and failure modes

- **Empty / invalid segments (ψ):** no segments → Σ(L) = 0 → ψ undefined (division by zero), rejected; a
  segment with Length ≤ 0 or negative CN corrupts the weighted mean → invalid input.
- **WGD boundary:** strict `> 0.5` on the GRCh38 autosomal genome (true at half+1, false at exactly
  half, false at half−1).
- **WGD major-CN gate:** all-1:1 autosomes → 0.0 → not doubled; chrX/chrY amplification only → 0.0 (sex
  chromosomes excluded); "chr"-prefixed autosome names recognised.
- **WGD empty input:** returns false (numerator 0 over a fixed reference denominator) — unlike the legacy
  overload, which throws on an empty supplied-length denominator.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for two post-segmentation genome-state
summaries. **Not for clinical or diagnostic use.** No source contradictions: Patchwork (the
length-weighted ploidy definition), Van Loo 2010 (the n-scale and > 2.7n aneuploidy threshold),
Bielski 2018 + facets-suite `is_genome_doubled` (the ≥ 50% / major-CN ≥ 2 WGD binary rule), and
UCSC / Ensembl (the reference autosome lengths) each cover a disjoint part and agree. The only
non-numeric decisions are input-shape reuse (`AlleleSpecificSegment` = Major + Minor CN) and the
resolved reference-genome WGD denominator — every chromosome length is the published UCSC value.
