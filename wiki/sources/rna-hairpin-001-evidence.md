---
type: source
title: "Evidence: RNA-HAIRPIN-001 (hairpin-loop + stem free-energy — Turner 2004)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-HAIRPIN-001-Evidence.md
sources:
  - docs/Evidence/RNA-HAIRPIN-001-Evidence.md
source_commit: 8346ce2d97d95f5b806caf203fd3d1dc19271cf5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-HAIRPIN-001

The validation-evidence artifact for test unit **RNA-HAIRPIN-001** — **Hairpin Loop and Stem
Free-Energy Calculation** (Turner 2004 nearest-neighbor model). It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept
is [[rna-free-energy-turner-model]]. [[test-unit-registry]] tracks the unit.

This unit is **not** a stem-loop *enumerator* (that is [[rna-stem-loop-enumeration]], RNA-STEMLOOP-001,
which *finds* hairpins and consumes these energy terms for its tetraloop bonus) and is **not**
miRNA-specific: it is a focused deep-dive on the two hairpin/stem **energy** methods
`CalculateHairpinLoopEnergy` and `CalculateStemEnergy` on `RnaSecondaryStructure`. It therefore **overlaps and complements**
[[rna-free-energy-turner-model|RNA-ENERGY-001]] (the aggregate energy layer, which already exposes
these same methods) rather than defining a new concept, and it is **distinct** from
[[pre-mirna-hairpin-detection]] (which *finds* pre-miRNA stem-loops and only sums a hairpin
`FreeEnergy` as a downstream step). The energies attach to adjacent [[rna-base-pairing|base-pair]]
stacks of a [[rna-dot-bracket-notation|dot-bracket]] hairpin.

## What this file records

- **Hairpin free-energy formula (loop > 3 nt):** `ΔG°37 = initiation(n) + terminal mismatch +
  (UU/GA first mismatch) + (GG first mismatch) + (special GU closure) + (all-C penalty)`. The
  hairpin decomposes into a **loop component** (initiation + terminal mismatch + first-mismatch
  bonuses) and a **helix/stem component** (per-stack stacking + AU-end penalties).
- **3-nt loop formula:** `initiation(3) + all-C penalty` only — **no** sequence-dependent
  first-mismatch term for 3-nt loops.
- **Authoritative source:** NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04/), retrieved via
  Internet Archive Wayback snapshots (the live NNDB server was down on the access date). Primary
  literature: **Mathews et al. (2004)** PNAS 101:7287-7292; Turner & Mathews (2010) NAR 38:D280.

### Parameters worth pinning

- **Hairpin initiation (kcal/mol):** 3→5.4, 4→5.6, 5→5.7, 6→5.4, 7→6.0, 8→5.5, 9→6.4, 10→6.5;
  n>9 extrapolated `initiation(9) + 1.75·RT·ln(n/9)`.
- **First-mismatch bonuses:** UU or GA −0.9; GG −0.8. **Special GU closure** −2.2 (applies only to a
  `G-U` closing pair, **not** `U-G`, preceded by two Gs — a documented asymmetry).
- **All-C loop penalty:** 3-nt flat **+1.5** (C3); >3-nt linear `A·n + B` with A=+0.3/nt, B=+1.6.
- **Stem/stacking (WC ΔG°37):** e.g. `CA/GU` −2.11, `GU/CA` −2.24, `GG/CC` −3.26, `GC/CG` −3.42;
  **+0.45 per AU end**. A helix of **P** base pairs has **P−1** stacks (0 pairs → 0 stacking).
- **Special hairpins override the model** — experimental totals keyed by closing pair + loop replace
  the additive calculation: triloop `CAACG` 6.8, `GUUAC` 6.9; tetraloop `CCUCGG` 2.5; hexaloop
  `ACAGUGUU` 1.8.

## Datasets / oracles

- **NNDB Example 1 (closing A-U, 6-nt loop A…A):** loop component +4.6 (init(6) 5.4 + terminal
  mismatch A·A −0.8); helix component −6.01 (3 stacks −2.11/−2.24/−2.11 + AU-end +0.45); total
  **−1.4 kcal/mol**.
- **NNDB Example 2 (closing A-U, 5-nt loop G…G, GG first mismatch):** loop component +4.1 (init(5)
  5.7 − 0.8 terminal mismatch − 0.8 GG bonus); helix −6.01; total **−1.9 kcal/mol**.
- **Special/penalty loops:** triloop `CAACG` → 6.8; tetraloop `CCUCGG` → 2.5; all-C 3-nt loop →
  init(3) 5.4 + C3 1.5 = **6.9**; all-C 4-nt loop → init(4) 5.6 + tm(GCCC −0.7) + (0.3·4+1.6) =
  **7.7**.

## Corner cases

- **Loops < 3 nt are prohibited** (no thermodynamic value defined — prohibitive energy).
- **3-nt loops** receive no first-mismatch term (only initiation + optional all-C penalty).
- **Special-loop table overrides** the additive model when the closing pair + loop matches.
- **Empty / single base-pair stem** contributes **0** stacking energy (P−1 stacks ⇒ 0 for P ≤ 1).

## Deviations and assumptions

All correctness-affecting parameters are source-backed; the only recorded item is a **display-precision
choice** — NNDB tabulates final ΔG°37 to one decimal, the implementation rounds intermediate sums to
two decimals; tests assert with `.Within(1e-9)` against the exact arithmetic sum of the cited
parameters. **No source contradictions.**
