---
type: source
title: "Evidence: GENOMIC-TANDEM-001 (Tandem repeat detection, GenomicAnalyzer.FindTandemRepeats)"
tags: [validation, genomic-analysis]
doc_path: docs/Evidence/GENOMIC-TANDEM-001-Evidence.md
sources:
  - docs/Evidence/GENOMIC-TANDEM-001-Evidence.md
source_commit: 4ee1ab19359eab0c144e9a59219013e0c0f4ec91
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: GENOMIC-TANDEM-001

The validation-evidence artifact for test unit **GENOMIC-TANDEM-001** — **Tandem Repeat
Detection** (`GenomicAnalyzer.FindTandemRepeats`): an exact detector for two-or-more
contiguous copies of a repeat unit (microsatellites / STRs and larger tandems). One
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the algorithm and its place among the repeats family are synthesized in
[[repetitive-element-detection]] (tandem sub-problem). See [[test-unit-registry]] for how
units are tracked.

**Consolidation.** GENOMIC-TANDEM-001 is a **duplicate Registry entry for the identical
method** already delivered under **REP-TANDEM-001** (same class, same O(n²) brute-force
scan, canonical fixture `GenomicAnalyzer_TandemRepeat_Tests.cs`, 27 tests). It is resolved
by consolidation, not re-implementation — no new or duplicate tests were created.

## What this file records

- **Online sources** (accessed 2026-06-14):
  - **Benson, G. (1999) — Tandem Repeats Finder**, *Nucleic Acids Research* 27(2):573–580
    (rank 1, the foundational primary paper) — the formal definition ("two or more
    contiguous copies of a pattern of nucleotides"), **period** = unit length between
    matching positions of adjacent copies, **copy number** = number of copies, and the
    **k ≥ 2** minimum. Benson's TRF detects *approximate* copies; the repository detects
    only exact copies.
  - **Wikipedia "Tandem repeat"** (rank 4) — verbatim definition and the worked example
    `ATTCG ATTCG ATTCG`; unit-length classification (di-/tri-nucleotide; microsatellite/STR
    short units, minisatellite 10–60 nt, macrosatellite ≈1,000 nt); ~8% of the human genome,
    implicated in >50 lethal diseases; detection via suffix trees/arrays.
- **Datasets (documented oracles):**
  - Wikipedia worked example: `ATTCGATTCGATTCG` → unit `ATTCG`, period 5, 3 copies, total
    length 15, start 0.
  - Canonical trinucleotide fixture (REP-TANDEM-001 M1): `ATGATGATG` → unit `ATG`, 3
    copies, start 0.
- **Recommended coverage** — two MUST tests (the two oracles above), SHOULD (`minRepetitions`
  floor k ≥ 2 and `minUnitLength` threshold honored), COULD (no-repeat / empty input →
  empty). All already implemented and verified under REP-TANDEM-001.

## Corner cases and assumption

- **Exact vs approximate:** the detector reports only **exact** contiguous copies, not
  Benson's approximate copies (documented simplification, algorithm doc §5.3). Over exact
  repeats the output matches the formal definition.
- **Multiple period interpretations:** the same region (e.g. `AAAA`) satisfies the
  definition for several periods (period 1 ×4, period 2 ×2); the detector reports **each**
  unit-length interpretation meeting the threshold — it does **not** canonicalize to the
  primitive period (algorithm doc §6). This diverges from the annotation `RepeatAnalyzer`
  path ([[annot-repeat-001-evidence|ANNOT-REPEAT-001]]), which prefers the shortest period.
- **Minimum copies:** fewer than two contiguous copies is not a tandem repeat; the
  `minRepetitions` default is 2, matching the floor.

**ASSUMPTION (exact-copy restriction):** only exact contiguous copies are reported (a
Framework/Simplified limitation) — no correctness-affecting parameter is invented, since
over exact repeats the output matches the formal definition. No source contradictions.
