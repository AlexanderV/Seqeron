---
type: source
title: "Evidence: META-BIN-001 MarkerQC addendum (CheckM single-copy marker-gene completeness/contamination)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-BIN-001-MarkerQC-Evidence.md
sources:
  - docs/Evidence/META-BIN-001-MarkerQC-Evidence.md
source_commit: c21e0be5032ffdc39a0e53405a4c7fbd09482958
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-BIN-001 (MarkerQC addendum)

The **marker-gene quality-metrics addendum** to test unit **META-BIN-001** — the CheckM-style
**single-copy marker-gene completeness & contamination** now built on top of the existing
TNF/coverage binning. It **builds the honest residual** the base
[[meta-bin-001-evidence|META-BIN-001 evidence]] left unbuilt: the base binning path,
its proxy completeness (length ratio) and contamination (GC std-dev), and the default `BinContigs`
entry point are **unchanged**. The method is synthesized in the [[metagenomic-binning]] concept;
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the artifact
pattern. See `docs/Evidence/META-BIN-001-MarkerQC-Evidence.md`.

## What this file records

The two published CheckM equations (Parks et al. 2015) are implemented verbatim over **collocated
marker sets** `M`, with the exact arithmetic pinned against the CheckM reference implementation
(`MarkerSet.genomeCheck` in `Ecogenomics/CheckM`):

```
Completeness  = 100 · ( Σ_{s∈M} |s ∩ G_M| / |s| ) / |M|      (Eq. 1)
Contamination = 100 · ( Σ_{s∈M} Σ_{g∈s} C_g / |s| ) / |M|    (Eq. 2), C_g = N−1 for a gene seen N≥1×, else 0
```

- A **multi-copy marker** counts **once** toward `present` (completeness uses `|s ∩ G_M|`, not the
  copy count) and contributes `N−1` to contamination. Completeness ≈ fraction of unique single-copy
  genes present; contamination ≈ how many are present in multiple copies.
- Grouping into **marker sets** (rather than independent markers) down-weights correlated,
  consistently-collocated genes — the metric is averaged per-set then over the number of sets.

## Bundled marker sets (licence-gated)

Markers feed the CheckM formula through `EstimateBinQualityFromMarkerCounts` /
`EstimateBinQualityFromMarkers` / `DetectMarkers`; profiles load via `LoadMarkerHmms` (caller-supplied)
or the bundled loaders:

| Bundle | Markers | Basis |
|--------|---------|-------|
| 9 universal ribosomal Pfams | PF00318/00177/00410/00380/00338/00411/00203/00687/00297 (S2/S7/S8/S9/S10/S11/S19/L1/L3) | Xu et al. 2022 USCGs; each its own singleton set |
| bac120 Pfam subset (`LoadBundledBacterialMarkerHmms` / `BundledBacterialMarkerSets`) | 6 (PF00380, PF00410, PF00466, PF01025, PF02576, PF03726) | GTDB bac120 (Parks 2018) |
| ar122 Pfam subset (`LoadBundledArchaealMarkerHmms` / `BundledArchaealMarkerSets`) | 35 | GTDB ar122 |

- **Pfam profiles are CC0** (public domain) → embedded as-is (HMMER3/f, `ALPH amino`, `GA1`
  per-sequence gathering threshold). Distinct Pfam union bac120 ∪ ar122 = 39 accessions.
- **TIGRFAM-defined bac120/ar122 members are CC BY-SA 4.0 (share-alike) → NOT bundled**; the
  caller supplies them via `LoadMarkerHmms`. Only the CC0 Pfam subsets ship.

## Oracles

- **Hand-derived synthetic bin (pins Eqs. 1–2):** 3 sets — s1={A,B} A=1/B=0, s2={C,D,E} C=2/D=1/E=1,
  s3={F} F=1 → `comp=0.5+1.0+1.0=2.5` ⇒ **Completeness = 250/3 = 83.333…%**; `cont=0+1/3+0` ⇒
  **Contamination = 100/9 = 11.111…%**. MarkersPresent=5, MarkerSetCount=3, MarkerCount=6.
- **HMM detection integration:** E. coli uS8 (UniProt P0A7W7) vs bundled PF00410 → **≈176.38 bits ≥
  GA1 24.0 ⇒ hit (count 1)**; all 8 other ribosomal families **< 0 bits ⇒ no hit**. Over the 9
  singleton sets, exactly one fully present ⇒ **Completeness = 100/9 ≈ 11.111%, Contamination = 0**.
- **Domain-level detection:** E. coli GrpE (P09372) vs bac120 PF01025 (GA1 25.8) → 1 of 6 sets ⇒
  **100/6 ≈ 16.667%**; uS8 vs ar122 PF00410 → 1 of 35 ⇒ **≈2.857%**, contamination 0. Family-level
  specificity (GrpE matches only PF01025; uS8 only PF00410).

## Corner cases and assumptions

- **Empty marker set** excluded (no `÷|s|` div-by-zero); **|M|=0** → both metrics 0 (undefined-guard).
  Missing marker `C_g=0`. Triplicated marker `C_g=N−1=2`.
- **ASSUMPTION — marker present iff Viterbi bit score ≥ GA1.** The numeric gate (GA1) is sourced from
  the HMM file, but this engine's **glocal (whole-sequence) Plan7 Viterbi** log-odds score differs
  from HMMER's local + null2-corrected `hmmsearch`, so absolute bit scores diverge (a documented
  engine difference, consistent with `ProteinMotifFinder.FindDomainsByHmm`). The true-positive
  separation is nonetheless decisive (uS8 vs PF00410 = +176 bits vs all others < 0). The exact-formula
  tests feed counts directly and don't depend on this.
- **ASSUMPTION — bundled markers as singleton sets.** The bundled universal families are distinct,
  non-collocated Pfams, so each is its own set (|s|=1). CheckM's operon-based collocation grouping
  needs the full lineage DB — left as the remaining residual.

Sources: Parks et al. 2015 (CheckM, *Genome Res* 25:1043) + `MarkerSet.genomeCheck` reference impl +
Xu et al. 2022 (universal single-copy genes) + EMBL-EBI InterPro Pfam HMM API (CC0) + Parks et al.
2018 / GTDB-Tk bac120/ar122 marker lists + Pfam/TIGRFAM licence docs + UniProt true-positive proteins.
No source contradictions.
