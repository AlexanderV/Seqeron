---
type: source
title: "Evidence: GENOMIC-ORF-001 (Open Reading Frame detection, six-frame ATG→stop)"
tags: [validation, annotation]
doc_path: docs/Evidence/GENOMIC-ORF-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-ORF-001-Evidence.md
source_commit: b8df572dba1d19e639aaf1a2ac9e35e2d68aea4e
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-ORF-001

The validation-evidence artifact for test unit **GENOMIC-ORF-001** — **Open Reading
Frame (ORF) detection**, six-frame ATG→first-in-frame-stop enumeration
(`GenomicAnalyzer.FindOpenReadingFrames`). It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm, its contract, invariants, worked oracles, and corner cases are synthesized
in [[open-reading-frame-detection]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Rosalind "ORF"** (rank 4, worked-example authority) — ORF definition (start
    codon → stop codon, no internal stop), six reading frames (3 forward + 3 reverse
    complement), start AUG/ATG, stops UAA/UAG/UGA (TAA/TAG/TGA), return all **distinct**
    protein candidates translated **until** a stop; nested ORFs sharing a stop are both
    reported.
  - **Wikipedia "Open reading frame"** (rank 4) — "spans of DNA between start and stop
    codons", six frame translations on double-stranded DNA, length divisible by three,
    minimal-length convention (e.g. 100 codons).
  - **NCBI ORFfinder** (rank 5, tool) — start-codon options ("ATG only" default),
    minimal-ORF-length selectable in **nucleotides** (30/75/150/300/600), inclusive
    lower bound (length ≥ threshold).
  - **NCBI Genetic Codes transl_table=1** (rank 2, official spec) — standard start ATG
    (Met), stops TAA/TAG/TGA all "*".
- **Datasets (documented oracles):**
  - *Rosalind ORF sample* (`Rosalind_99`): the 94-nt input → the **4 distinct** proteins
    `MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`, `MTPRLGLESLLE` (last two share a
    stop); independent re-derivation (standard code, every-ATG-to-first-in-frame-stop,
    both strands, distinct) reproduces exactly these four.
  - *Single canonical forward ORF*: `ATGAAAAAATAA` (12 nt) → one ORF, `Sequence`
    `ATGAAAAAATAA`, position 0, frame 1, protein candidate `MKK` (stop excluded).

## Deviations and assumptions

**Deviations: none** (one *fixed* pre-existing deviation is noted — a greedy scan that
missed nested ORFs sharing a stop was corrected to match the Rosalind sample). Three
**assumptions** (all source-anchored, none affecting *which* ORFs are detected):

1. **ORF `Sequence`/`Length` includes the terminating stop codon** (so `Length % 3 == 0`,
   boundaries explicit; Wikipedia "length divisible by three … bounded by stop codons").
   The reported protein candidate **excludes** the stop (Rosalind translates "until a stop").
2. **`minLength` measured in nucleotides**, inclusive lower bound (NCBI ORFfinder nt
   options); default 100 nt retained from the existing public API, any value honored.
3. **Standard start codon ATG only** (NCBI ORFfinder default "ATG only"); alternative
   initiation codons (GTG/TTG) and non-standard codes are Not Implemented for this unit.

Recommended coverage (MUST): Rosalind six-frame → exactly the 4 distinct proteins;
`ATGAAAAAATAA` → one ORF / pos 0 / frame 1 / protein `MKK`; nested ORFs sharing a stop
both reported; ATG with no in-frame stop → no ORF; reverse-complement-only ORF detected
with `IsReverseComplement=true`; minLength excludes shorter / includes exactly-at-threshold;
invariants (starts ATG, ends stop, length %3==0). SHOULD: lowercase input handled;
all three stops recognized. COULD: empty/too-short → empty. No contradictions among sources.

## Scope note (distinct from the annotation ORF finder)

This unit validates the **`GenomicAnalyzer.FindOpenReadingFrames`** surface only
(ATG-only, six-frame always, `minLength` in **nucleotides**). It is **not** the
annotation-layer `GenomeAnnotator.FindOrfs` (test unit ANNOT-ORF-001), which recognizes
the prokaryotic start set ATG/GTG/TTG, measures `minLength` in **amino acids**, and
takes `searchBothStrands`/`requireStartCodon` flags — a deliberately non-contract-equivalent
sibling. Only `docs/algorithms/Analysis/Open_Reading_Frame_Detection.md` is validated here.
