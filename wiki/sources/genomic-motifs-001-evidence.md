---
type: source
title: "Evidence: GENOMIC-MOTIFS-001 (Known motif search / multi-pattern exact substring matching)"
tags: [validation, motif]
doc_path: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
source_commit: 58f37bc5de666c59a60b8e7997c0894c16768c96
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-MOTIFS-001

The validation-evidence artifact for test unit **GENOMIC-MOTIFS-001** — **Known
Motif Search**, multi-pattern **exact** substring matching (`GenomicAnalyzer.FindMotif`
over a set of known motifs). It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its
contract, invariants, worked oracles, and corner cases are synthesized in
[[known-motif-search]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Tufts COMP 150GEN exact-matching notes** (rank 1, Gusfield-derived) — the exact
    matching problem *"find all occurrences of a pattern string P of length n in a text
    string T of length m"* returning the set of **all** start positions, and the
    canonical **overlapping** rule: `P=aaa`, `T=aaaaa` → **three** occurrences
    (positions 0,1,2). A search that skips past a match (non-overlapping/greedy) under-reports.
  - **Biopython `Bio.Seq`** (rank 3, reference implementation) — `Seq.search` accepts a
    **list of multiple query motifs** and yields `(index, substring)` per hit (the
    multi-motif "known motif search" semantics); `count_overlap` confirms overlaps are
    counted (`Seq("AAAA").count_overlap("AA")` → **3** vs plain `count` → **2**); and the
    module **upper-cases** sequence data before processing.
  - **Wikipedia "Restriction site"** (rank 4, for a REBASE/NEB-held fact) — EcoRI
    recognizes the palindrome **`GAATTC`** (length 6), supplying a real biological motif
    for the fixture.
- **Datasets (documented oracles):**
  - *Overlapping homopolymer* (Gusfield): T=`AAAAA`, P=`AAA` → starts **{0,1,2}** (count 3).
  - *EcoRI motif*: T=`GAATTCAAAGAATTC`, P=`GAATTC` → starts **{0,9}** (count 2).
  - *Multi-motif set* (Biopython): T=`ACGTACGTAA`, set `{ACGT, AA, TTT}` → `ACGT`→{0,4},
    `AA`→{8}, `TTT`→(absent, omitted from the result).

## Deviations and assumptions

**Deviations: none** — exact all-occurrences (overlapping) matching as defined by
Gusfield/Biopython. Two **API-shape assumptions** (not algorithm-correctness params,
no authoritative source defines them):

1. **Empty motif** — no source defines an empty-pattern occurrence set as meaningful
   (a suffix tree's `FindAllOccurrences("")` would return every position). Repository
   convention (`FindMotif` returns `Array.Empty<int>()` for null/empty/whitespace)
   is reused: an empty/whitespace motif contributes no entry.
2. **Result key normalization** — motifs are **upper-cased** before search (consistent
   with Biopython upper-casing DNA and `DnaSequence` upper-casing on construction); the
   result is keyed by the upper-cased motif.

Recommended coverage (MUST): overlapping fully reported (`AAA` in `AAAAA`→{0,1,2});
multi-motif set → per-motif position lists, absent motifs omitted; real motif `GAATTC`
(EcoRI) at all starts; returned positions sorted ascending + 0-based (deterministic
contract). SHOULD: empty motif set → empty result; lower/mixed-case motif normalized,
key upper-cased. COULD: single-character motif counts all (consecutive/overlapping)
positions. No contradictions among sources — Gusfield and Biopython agree that all
occurrences including overlaps are reported.
