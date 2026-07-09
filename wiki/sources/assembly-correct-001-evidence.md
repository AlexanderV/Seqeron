---
type: source
title: "Evidence: ASSEMBLY-CORRECT-001 (K-mer Spectrum Read Error Correction)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-CORRECT-001-Evidence.md
source_commit: 4167da78fdda4821db88c34df210958d7da43cde
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-CORRECT-001

The validation-evidence artifact for test unit **ASSEMBLY-CORRECT-001** (two-sided
k-mer-spectrum read error correction). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm's trusted/untrusted
model and correction rule are summarized in [[kmer-spectrum-error-correction]], the anchor for the
assembly CORRECT family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13, quoted verbatim):
  - **Musket** (Liu, Schmidt & Maskell 2013, *Bioinformatics* 29(3):308; rank 1 paper / 3
    reference tool) — the two-sided correction rule. A k-mer is *trusted* iff its multiplicity
    exceeds a coverage cut-off (chosen automatically from "the smallest density around the valley"
    of the coverage histogram), otherwise *untrusted*; a base is *trusted* if covered by any
    trusted k-mer. Two-sided correction "aims to find a unique alternative base that makes all
    k-mers that cover position *i* trusted," evaluating both the leftmost and rightmost covering
    k-mers, and conservatively "assumes that there is at most one substitution error in any k-mer."
    **Ambiguity rule:** if more than one alternative makes both trusted, the base is left unchanged.
  - **Quake** (Kelley, Schatz & Salzberg 2010, *Genome Biology* 11:R116; rank 1) — corroborates
    the trusted/untrusted (high- vs low-coverage) split, treats bases in the right-most trusted
    k-mer as correct (removing them from the error region), searches single-base *edits* until "a
    set of corrections C makes all k-mers in the region trusted," and localizes errors to the
    intersection then union of the read's untrusted k-mers. Substitution/nucleotide-edit model.
  - **Mining statistically-solid k-mers** (Song & Florea 2018, PMC6311904; rank 1) — the same
    idea under *solid* (frequent, error-free) vs *weak* (below threshold f0) k-mers; a correction
    is accepted only when forward and backward searches meet and the shared k-mers are all solid.
- **Datasets** — (1) single-substitution worked example: reads `ACGTACGT`×3 (true) +
  `ACGTTCGT`×1 (A→T at index 4), `k=3`, cut-off `2` ⇒ spectrum ACG=7,CGT=8,GTA=3,TAC=3 trusted vs
  GTT=1,TTC=1,TCG=1 untrusted; index 4 covered only by untrusted k-mers; only base `A` makes all
  covering k-mers trusted ⇒ output `ACGTACGT`×4. (2) Ambiguity example: reads `A,A,C,C,G`, `k=1`,
  cut-off `2` ⇒ 1-mers A=2,C=2 trusted; correcting read `T` — both `A` and `C` yield a trusted
  1-mer ⇒ ambiguous ⇒ `T` unchanged.
- **Corner cases / failure modes** — ambiguous position (>1 valid alternative) ⇒ unchanged;
  ≤1-error-per-k-mer assumption (multiple close errors may be uncorrectable); no correcting set ⇒
  read left uncorrected (Quake trim/discard out of scope); the frequency cut-off cannot perfectly
  separate erroneous from correct k-mers (solid-with-error / weak-without-error — general limit of
  k-spectrum methods).
- **Recommended coverage** — MUST: single-substitution fully corrected; trusted-covered base
  never modified; ambiguous position unchanged; no valid base ⇒ unchanged; output count and
  per-read length preserved (substitution-only); null reads ⇒ `ArgumentNullException`, `kmerSize`
  < 1 ⇒ `ArgumentOutOfRangeException`. SHOULD: correct reads pass through unchanged;
  case-insensitive input. COULD: reads shorter than k contribute no k-mers and pass unchanged.

## Assumptions (from the artifact)

One assumption record: **default parameter values are non-behavioral.** Musket/Quake pick `k` and
the coverage cut-off automatically from the coverage-histogram valley; the library instead exposes
both as parameters with fixed defaults (`kmerSize=15`, `minKmerFrequency=2`). Because every
behavioral test passes `k` and cut-off explicitly, correctness for a given `(k, cut-off)` is fully
source-defined and the defaults do not affect the tested contract.

No contradictions among the three sources — Musket, Quake, and the solid-k-mer paper describe the
same trusted/untrusted (solid/weak) two-sided correction model. The only assumption is the
non-behavioral default-parameter choice above.
