---
type: source
title: "Evidence: RESTR-FILTER-001 (Restriction-enzyme filtering — blunt/sticky cutters + recognition-length range)"
tags: [validation, moldesign]
doc_path: docs/Evidence/RESTR-FILTER-001-Evidence.md
sources:
  - docs/Evidence/RESTR-FILTER-001-Evidence.md
source_commit: 94748c8ca67c5c8b70dcdcf27a2c2d69087c0eaa
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RESTR-FILTER-001

The validation-evidence artifact for test unit **RESTR-FILTER-001** — **Restriction
Enzyme Filtering**: selecting enzymes from the library by end type (`GetBluntCutters()`,
`GetStickyCutters()`) and by recognition-site length (`GetEnzymesByCutLength(min, max)`
plus the single-length overload). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the contract, invariants and
the blunt/sticky partition are synthesized in [[restriction-enzyme-filtering]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Wikipedia "Sticky and blunt ends"** (rank 4) — blunt end = *"both strands terminate
    in a base pair"* (no overhang); overhang = *"a stretch of unpaired nucleotides in the
    end of a DNA molecule"* (longer overhangs = cohesive/sticky ends), on either strand →
    5' or 3'; *"blunt ends are always compatible with each other."* Establishes that
    blunt and sticky are the **only** two categories (total, disjoint partition).
  - **Wikipedia "Restriction enzyme"** (rank 4) — Type II recognition sites are *"usually
    undivided and palindromic and 4–8 nucleotides in length"*; a **center** cut → blunt,
    a **staggered** cut → sticky overhangs. Worked examples: EcoRI `G^AATTC` (5' overhang),
    SmaI `CCC^GGG` (blunt), KpnI `GGTAC^C` (3' overhang), PstI `CTGCA^G` (3' overhang).
  - **Wikipedia "List of restriction enzyme cutting sites"** (rank 4) — corroborates the
    length categories: EcoRI/BamHI/HindIII/XhoI/SacI/PstI/NdeI = 6 bp; TaqI/AluI/HaeIII =
    4 bp; NotI = 8 bp.
  - **SfiI — PMC/REBASE** (rank 1/5) — SfiI recognizes the **interrupted palindrome**
    `5'GGCCNNNN^NGGCC3'` (13-nt string, 5-base `NNNNN` degenerate spacer) — a divided
    palindrome whose length 13 lies **outside** the inclusive `[4, 8]` filter.
  - **KpnI NEB/REBASE** (rank 3) — KpnI generates a 3', 4-base overhang (`G_GTAC^C`) →
    sticky, not blunt.
  - **EcoRI NEB/REBASE** (rank 3) — EcoRI cuts `G^AATTC` leaving a 4-base 5' `AATT`
    overhang → sticky.

- **Dataset (documented oracle — library enzymes with known end type + length):**
  SmaI `CCCGGG` 6 blunt; EcoRV `GATATC` 6 blunt; AluI `AGCT` 4 blunt; HaeIII `GGCC` 4
  blunt; EcoRI `GAATTC` 6 sticky(5'); KpnI `GGTACC` 6 sticky(3'); PstI `CTGCAG` 6
  sticky(3'); NotI `GCGGCCGC` 8 sticky; TaqI `TCGA` 4 sticky; SfiI `GGCCNNNN^NGGCC` 13
  (interrupted) sticky(3', outside 4–8).

## Deviations and assumptions

**Deviations: none** — the blunt/sticky partition and the 4–8 nt undivided-site range
follow the cited definitions exactly. One **API-shape assumption** (no authoritative
source defines a software filter's boundary convention):

1. **Range bounds are inclusive** — `GetEnzymesByCutLength(min, max)` treats both bounds as
   inclusive (`min ≤ len ≤ max`), the conventional min/max meaning. The recognition-length
   values themselves are source-backed; only the boundary convention is conventional.

Recommended coverage (MUST): `GetBluntCutters()` returns only symmetric cutters (SmaI,
EcoRV, AluI, HaeIII) and excludes EcoRI/KpnI; `GetStickyCutters()` returns only
overhang cutters (EcoRI, KpnI, PstI) and excludes SmaI; the two sets are **disjoint and
their union is the full library**; `GetEnzymesByCutLength(min, max)` returns exactly the
enzymes with length in `[min, max]` inclusive (`6→6` = all 6-cutters; `9→10` = none;
`[4, 8]` = whole undivided library **except** the 13-nt interrupted SfiI). SHOULD:
`min > max` → empty; single-length overload agrees at `min == max`. COULD: negative/zero
bounds → empty. No contradictions — all sources agree end type is a blunt-or-overhang
dichotomy and the undivided-site length range is 4–8 nt.
