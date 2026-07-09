---
type: source
title: "Evidence: META-RESIST-001 (antibiotic-resistance gene detection — ResFinder-style identity/coverage dual-threshold best-match screen)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-RESIST-001-Evidence.md
sources:
  - docs/Evidence/META-RESIST-001-Evidence.md
source_commit: c81ef58a4d2d8fd4a9ceb1e322d2c6a1ee237cfc
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-RESIST-001

The validation-evidence artifact for test unit **META-RESIST-001** — **antibiotic-resistance gene
detection**: screen assembled contigs against a caller-supplied reference database of resistance genes and
report the best-matching gene per contig when its BLAST percent identity and its coverage of the reference
gene both clear user-selectable thresholds (`MetagenomicsAnalyzer.FindAntibioticResistanceGenes`). This is
the **seventh ingested unit of the Metagenomics family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own concept,
[[antibiotic-resistance-gene-detection]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, no contradictions):**
  - **Zankari et al. 2012 — original ResFinder paper** (J Antimicrob Chemother 67(11):2640–2644;
    accessed 2026-06-13, authority rank 1) — the core method: "all genes from the ResFinder database were
    BLASTed against the assembled genome, and the best-matching genes were given as output"; a gene is
    reported only if it covers **≥ 2/5 of the reference gene's length**; default ID 100%, user-selectable;
    %ID = "the percentage of nucleotides that are identical between the best-matching resistance gene in
    the database and the corresponding sequence in the genome"; the study operated at **ID = 98%** because
    lower thresholds "give too much noise (e.g. fragments of genes)".
  - **genomicepidemiology/resfinder README** (rank 3, reference implementation) — default identity
    threshold `-t` **0.80** (`CGE_RESFINDER_GENE_ID`) and default breadth-of-coverage threshold `-l`
    **0.60** (`CGE_RESFINDER_GENE_COV`); coverage = proportion of the reference gene covered by the
    alignment, 0–1.
  - **Sci Rep 2023 13:15543 (carbapenem-resistant K. pneumoniae)** (rank 1) — operating thresholds
    "98% identity" and "60% coverage"; the 60% coverage floor exists so genes "lying on the edge of a
    contig or spread over two contigs are not missed, due to non-perfect assembly".
  - **JAC 2016 71(9):2484–2492 (AMR-method benchmarking)** (rank 1) — confirms ResFinder operated at
    98% identity / 60% coverage, coverage defined relative to the reference gene length.
  - **Heng Li 2018, "On the definition of sequence identity"** (rank 3) — BLAST identity = matching bases
    ÷ alignment columns; in a gapless alignment there are no gap columns, so the denominator = number of
    aligned positions (worked example 43/50 = 86%).
  - **CARD Resistance Gene Identifier (RGI)** (rank 5) — a "Perfect" RGI match is 100% identical to the
    reference over its entire length (identity = 1.0 at full coverage is the unambiguous top of scale);
    RGI ranks candidates by bit score and reports the best hit, corroborating the single best-matching-gene
    output convention.

- **Extracted formulas & rules:** identity = matches / gapless-alignment-window length `w`;
  coverage = `w / m` (m = reference gene length); reporting rule = report reference gene for a contig only
  if `identity ≥ idThreshold` **AND** `coverage ≥ covThreshold`; best-matching gene per contig =
  max identity, ties broken by greater coverage.

- **Documented corner / failure cases:**
  - Sub-threshold identity (below ~98% in the study) → gene fragments / noise → must be **rejected**,
    not reported.
  - Coverage below the floor (2/5 → 60% in later versions) → not reported.
  - Genes split across contigs / at contig edges are exactly why the coverage is measured against the
    **reference** length and the floor is < 1 — a hit need not span the full reference.

- **Datasets (documented oracles, derived arithmetically from the identity/coverage formulas** — the
  detector is generic, taking caller-supplied reference genes; CARD/ResFinder tables are not embedded):
  - **Exact full-length:** contig `AAACGTACGT` contains `CGTACGT` (m=7) → %ID = 7/7 = **1.0**,
    coverage = 7/7 = **1.0**.
  - **One mismatch, full length:** `CGTTCGT` vs `CGTACGT` → %ID = 6/7 ≈ **0.857142857**, coverage = **1.0**.
  - **Partial (contig edge), 4 of 7 bases:** contig ends with `CGTA`, ref `CGTACGT` → %ID = 4/4 = **1.0**,
    coverage = 4/7 ≈ **0.571428571**.

## Recommended test coverage (from the Evidence file)

MUST: exact full-length → ID 1.0 / cov 1.0 / reported; single-mismatch full-length → ID 6/7 / cov 1.0;
contig-edge partial hit scored by **reference** length, passes if ≥ coverage threshold; below identity
threshold → not reported; below coverage threshold → not reported; best-matching gene only (highest
identity) reported per contig; default thresholds equal ResFinder values (0.90 ID / 0.60 cov). SHOULD:
null/empty/invalid-threshold inputs raise the documented exceptions. COULD: tie-break by coverage when
identity ties.

## Deviations and assumptions

- **ASSUMPTION (ASM-01) — gapless (ungapped) alignment model.** The detector locates the best **ungapped**
  alignment (reference slid across the contig, no indels). ResFinder uses full **gapped** BLAST, but the
  BLAST identity formula (matches / alignment columns) is **identical for the gapless case** (no gap
  columns), so the formula is used exactly as cited. This affects output **only** for genes whose true
  alignment to the contig requires insertions/deletions; for substitution-only divergence and contig-edge
  truncation (the cases the coverage floor targets) it is exact. Documented as a scope simplification.

No source contradictions — Zankari's best-match/coverage rules, the ResFinder README defaults, the Sci Rep
/ JAC operating thresholds, Heng Li's identity denominator, and CARD's perfect-match / best-hit convention
are mutually consistent.

## Note — evidence 0.80/0.98 vs implementation 0.90 default

The evidence file records ResFinder's operational thresholds as **0.80** (GitHub README default `-t`) and
**0.98** (the value the Zankari study *selected* to suppress fragment noise), with 0.60 coverage
throughout. The Seqeron implementation ships a **0.90** identity default (a mid-point between the 0.80
service default and the 0.98 study threshold) with the same **0.60** coverage default; both are named,
user-selectable constants (`DefaultResistanceIdentityThreshold` / `DefaultResistanceCoverageThreshold`), so
the choice does not change the algorithm — it only sets the out-of-box operating point. Recorded here so
the 0.90 default is not mistaken for a source value; see [[antibiotic-resistance-gene-detection]].
