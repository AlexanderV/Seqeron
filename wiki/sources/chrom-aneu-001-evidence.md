---
type: source
title: "Evidence: CHROM-ANEU-001 (Aneuploidy detection — copy number from read depth)"
tags: [validation, chromosome]
doc_path: docs/Evidence/CHROM-ANEU-001-Evidence.md
sources:
  - docs/Evidence/CHROM-ANEU-001-Evidence.md
source_commit: b7e487d18a1de294e6480341479959136380a546
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CHROM-ANEU-001

The validation-evidence artifact for test unit **CHROM-ANEU-001** — aneuploidy detection: estimating
per-bin copy number from sequencing read depth and classifying whole-chromosome ploidy against the
normal disomic (CN=2) baseline. This is the first **Chromosome-analysis** family Evidence file and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
the depth→copy-number model, the classification rule, the confidence formula, the documented oracles,
and the two limitations are summarized in [[aneuploidy-detection]], the anchor for the chromosome
copy-number/ploidy family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Wikipedia — "Aneuploidy"** (encyclopedia) — the working definition: "an abnormal number of
    chromosomes in a cell", e.g. a human somatic cell with 45 or 47 instead of 46; mosaicism
    (variable CN across cells) and partial aneuploidy (from translocations).
  - **Wikipedia — "Copy number variation"** (encyclopedia) — the copy-number terminology ladder:
    nullisomy (0), monosomy (1), disomy (2, normal), trisomy (3), tetrasomy (4), pentasomy (5).
  - **Griffiths et al. (2000)** — *An Introduction to Genetic Analysis*, 7th ed. (textbook) — genetic
    background for ploidy nomenclature.
  - **Santaguida & Amon (2015)** — *Nat Rev Mol Cell Biol* 16(8):473–485, doi:10.1038/nrm4025 —
    chromosome mis-segregation / aneuploidy review.
  - **McCarroll & Altshuler (2007)** — *Nat Genet* 39(7 Suppl):S37–42, doi:10.1038/ng2080 — CNV and
    disease-association background.
- **Algorithm behaviour (from the implementation):**
  - **Copy number from depth:** `logRatio = log2(observedDepth / medianDepth)`, then
    `copyNumber = round(2^logRatio × 2)`, clamped to `[0, 10]`. The depth→CN table anchors the model
    (ratio 0.5 → CN 1 monosomy, 1.0 → CN 2 disomy, 1.5 → CN 3 trisomy, 2.0 → CN 4 tetrasomy).
  - **Whole-chromosome classification:** the dominant CN must span ≥ a `minFraction` of bins (default
    80%); CN 0→Nullisomy, 1→Monosomy, 2→Normal, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, >5→"Copy number = N".
  - **Confidence:** `1 − min(1, |expected − observed|)` where `expected = copyNumber/2` and
    `observed = 2^logRatio`; ∈ [0, 1].
- **Datasets (documented oracles):**
  - **Clinical examples** — Down (chr21 trisomy), Edwards (chr18 trisomy), Patau (chr13 trisomy),
    Turner (chrX monosomy 45,X), Klinefelter (chrX XXY / CN 3).
  - **Confidence at exact CN boundaries (S1 test)** — depth ratios 0.0/0.5/1.0/1.5/2.0 all yield
    confidence 1.0 because expected = observed exactly at each integer-CN ratio.
- **Edge cases** — mosaicism handled via `minFraction`; empty input → empty enumerable; zero or
  negative median depth → empty (guards division by zero); multiple chromosomes grouped by name;
  depth aggregated into bins by position/binSize; output bins ordered.

## Assumptions / limitations (from the artifact)

The file's §7 records two documented limitations (its §6 fallback is "not applicable — sufficient
authoritative sources found"):

1. **Sex chromosomes are not special-cased.** X/Y are treated like autosomes against the CN=2
   baseline, so a normal male (1 copy of X) would be flagged monosomic. Documented as a limitation;
   it does not affect the autosome-focused detection the unit targets.
2. **Partial aneuploidy.** Whole-chromosome classification requires a consistent CN across ≥80% of
   bins; sub-chromosomal (regional) CN changes from translocations are detected at the per-bin level
   but do not trigger a whole-chromosome call.

No contradictions among the sources — the two Wikipedia articles supply the aneuploidy definition and
the copy-number terminology ladder; the depth→CN model and the confidence formula are implementation
definitions the sources do not contradict.
