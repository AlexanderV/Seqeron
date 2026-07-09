---
type: source
title: "Evidence: ONCO-PLOIDY-001 (tumor ploidy = length-weighted mean segment CN + whole-genome-doubling detection)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-PLOIDY-001-Evidence.md
sources:
  - docs/Evidence/ONCO-PLOIDY-001-Evidence.md
source_commit: 57c2be1ccf184e08702ece85a7ad5a5d5618388c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-PLOIDY-001

The validation-evidence artifact for test unit **ONCO-PLOIDY-001** — **Tumor Ploidy Estimation**
(`EstimatePloidy`, the length-weighted mean of per-segment total copy number) and **Whole-Genome-Doubling
(WGD) detection** (`DetectWholeGenomeDoubling` + legacy `DetectWholeGenomeDoublingFromSuppliedLength`).
The **twenty-seventh ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in its own concept, [[tumor-ploidy-estimation-and-whole-genome-doubling]];
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Patchwork (Genome Biology 2010, PMC4053982)** (rank 1, primary — the ploidy definition) — verbatim:
    "The average ploidy, PloidyTum, is the average total copy number of all genomic segments weighted by
    segment length", i.e. **ψ = Σ(CN_i · L_i) / Σ(L_i)** over segments, where the averaged per-segment
    quantity is the **total** copy number (Major + Minor), not an allele-specific value.
  - **ASCAT — Van Loo et al. (PNAS 2010, PMID 20837533)** (rank 1) — establishes the **n-scale**
    (2n = diploid): the abstract reports "aneuploidy (>2.7n) in 45% of the cases" and near-triploid basal-like
    genomes, so ψ is on the n-scale and elevated whole-genome ploidy (> 2.7n) marks aneuploidy. ASCAT outputs
    a final tumour `ploidy` field.
  - **facets-suite `R/copy-number-scores.R` — `is_genome_doubled`** (rank 3, MSKCC reference impl, verbatim R,
    encoding the rank-1 Bielski 2018 definition) — WGD is called when the **autosome-restricted fraction of the
    genome with major copy number ≥ 2 is strictly greater than 0.5**:
    `frac_elevated_mcn = sum(length[mcn>=2 & chrom %in% 1:22]) / autosomal_genome`, `wgd = frac_elevated_mcn > 0.5`,
    with `mcn = tcn − lcn` (major = total − minor) and denominator
    `autosomal_genome = sum(chrom_info$size[chr %in% 1:22])` — a **reference chromosome-size table**, NOT the
    interrogated segments.
  - **Bielski et al. (Nature Genetics 2018, PMID 30013179)** (rank 1) — WGD is a per-sample binary genome-state
    classification ("whole-genome doubling in the tumors of nearly 30% of 9,692 patients"), operationalised by
    the facets-suite ≥50%/MCN≥2 rule.
  - **UCSC `hg38.chrom.sizes` / `hg19.chrom.sizes`** (rank 5, Ensembl GRCh38.p14 cross-verified) — the canonical
    autosome lengths; **Σ(chr1–22) = 2,875,001,522 bp (GRCh38)** and **2,881,033,286 bp (GRCh37)** — the WGD
    denominators, selected by a `ReferenceGenome { GRCh38, GRCh37 }` parameter (default GRCh38).

- **Documented corner cases / failure modes:**
  - **Ploidy (`EstimatePloidy`):** empty segment set → Σ(L)=0 → undefined (division by zero), rejected; a
    segment with Length ≤ 0 or negative CN corrupts the weighted mean → invalid input; a pure copy-number-2
    (all 1:1) genome has ψ exactly 2.0.
  - **WGD (`DetectWholeGenomeDoubling`):** the fraction is over **autosomes (chr1–22) only** (sex chromosomes /
    contigs excluded, so a single-X male does not bias the call); the threshold is **strict `> 0.5`** (exactly
    half at MCN ≥ 2 is NOT doubled); doubling uses the **major** allele CN ≥ 2, not total CN ≥ 2 — a 1:1
    (total 2) segment has major = 1 and is NOT elevated, whereas 2:0 (LOH) or 2:1 IS; a small fully-amplified
    region (e.g. 100 Mb) is well under 0.5 of the 2.875 Gb reference genome, so **no supplied-segment bias**;
    "chr"-prefixed autosome names are recognised.

- **Datasets (deterministic worked oracles):**
  - **Length-weighted ploidy worked example** (Patchwork): segments CN 2 / 4 / 3 at lengths 100 / 100 / 50 Mb →
    Σ(CN·L) = 750,000,000, Σ(L) = 250,000,000 → **ψ = 3.0**.
  - **Pure-diploid identity:** all 1:1 (total 2) segments → **ψ = 2.0** exactly.
  - **WGD against the GRCh38 reference genome** (denominator 2,875,001,522 bp; half = 1,437,500,761 bp):
    1,437,500,762 bp at MCN ≥ 2 → **true**; exactly 1,437,500,761 bp (half) → **false** (strict `>`);
    1,437,500,760 bp (half − 1) → **false**; 100 Mb fully amplified → 0.035 → **false**; all-1:1 autosomes
    (major 1) → 0.0 → **false**; chrX/chrY amplified only → 0.0 → **false**. Legacy
    `DetectWholeGenomeDoublingFromSuppliedLength` (denominator = Σ supplied length): 60% at MCN ≥ 2 → true,
    exactly 50% → false.

- **Coverage recommendations:** MUST test `EstimatePloidy` on the 3-segment worked example → 3.0, on pure-diploid
  → 2.0, that it is length-weighted (long CN-2 + short CN-4 weights toward 2), and that it rejects empty /
  Length ≤ 0 / negative-CN input; MUST test `DetectWholeGenomeDoubling` flips at the 0.5 boundary against the
  GRCh38 autosomal genome (true at half+1, false at exactly half, false at half−1), uses **major** CN not total CN
  (all-1:1 → not doubled), that the embedded GRCh38/GRCh37 tables equal the UCSC values (sums 2,875,001,522 /
  2,881,033,286 bp), and that a small amplified region / sex-chromosome-only elevation is NOT WGD; SHOULD test
  shared invalid-input validation (empty → false over a fixed denominator; legacy overload 60%→true / 50%→false /
  empty→throws); COULD test the GRCh37 selector disagreeing with GRCh38 near the boundary.

## Deviations and assumptions

- **ASSUMPTION — input shape.** Per-segment total copy number is supplied as the existing `AlleleSpecificSegment`
  record (shared with ONCO-LOH-001 / ONCO-HRD-001): total CN = Major + Minor, length = End − Start. Patchwork
  defines ploidy on per-segment total CN, and Major + Minor **is** that total and is also required to evaluate
  the major-CN ≥ 2 WGD rule. An input-shape reuse decision, not an invented numeric constant.
- **RESOLVED (2026-06-22) — WGD denominator.** The WGD fraction denominator is now the **reference autosomal
  genome length** (embedded UCSC hg38 / hg19 tables, Ensembl-cross-verified, `ReferenceGenome` selector, default
  GRCh38), replacing the earlier supplied-segment-length denominator, per facets-suite `autosomal_genome`. Only
  autosomal (chr1–22) segments contribute to the numerator. The legacy supplied-segment-length behaviour remains
  available via `DetectWholeGenomeDoublingFromSuppliedLength`. No invented constants — every chromosome length is
  the published UCSC value.

No source contradictions — Patchwork (ploidy definition), Van Loo 2010 (n-scale / aneuploidy threshold),
Bielski 2018 + facets-suite (the WGD binary rule), and UCSC/Ensembl (reference chromosome lengths) each cover a
disjoint part and agree.
