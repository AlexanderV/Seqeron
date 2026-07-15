---
type: source
title: "Evidence: PRIMER-TM-001-SPECIAL-LOOP (bundled special tri/tetraloop hairpin bonus tables)"
tags: [validation, primer]
doc_path: docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md
sources:
  - docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md
source_commit: c7a5a1030898a4cc5f776e764c05a46d17c80904
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PRIMER-TM-001-SPECIAL-LOOP

The validation-evidence artifact for test unit **PRIMER-TM-001** — the **special tri/tetraloop
hairpin bonus** sub-unit: the sequence-specific stability bonuses that libprimer3 `ntthal`
adds to a length-3 or length-4 hairpin loop, now **bundled** (not caller-supplied) and applied
inside the full `calc_hairpin` DP. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. This unit **completes** the
opt-in special-loop feature that [[primer-dimer-thermodynamics-tm]] previously described as an
unbundled, opt-in caller increment; the shared SantaLucia NN thermodynamics and unimolecular hairpin
Tm are synthesized there. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **libprimer3 `primer3_config/triloop.dh` + `triloop.ds`** (rank 3, vendored in primer3-py) —
    **16** triloop entries. Each line is `<loop> <value>` where `<loop>` is the **5-character**
    loop string *including the closing base pair* (closing-5′ + 3 loop nt + closing-3′;
    `readTLoop` copies 5 chars, thal.c L1161). ΔH verbatim (cal/mol): AG**·**AT loops −1500,
    CG**·**AG loops −2000, GG**·**AC loops −2000, TG**·**AA loops −1500 (four closing-pair
    families × four middle bases). **All 16 ΔS = 0.**
  - **libprimer3 `tetraloop.dh` + `tetraloop.ds`** (rank 3) — **76** tetraloop entries. Loop is
    the **6-character** string including the closing pair (`readTLoop` copies 6, thal.c L1164).
    Transcribed verbatim into `NtthalHairpin.TetraloopTable` (e.g. CGAAAG ΔH=−1100 ΔS=0;
    GGGGAC ΔH=−1100 ΔS=0; AAAAAT ΔH=500 ΔS=−650; ACTTGT ΔH=0 ΔS=4190).
  - **libprimer3 `thal.c` `calc_hairpin`** (rank 3, L2067–2146) — **how the bonus is applied**:
    loop energy starts from size-keyed `hairpinLoop{Enthalpies,Entropies}[loopSize−1]`; for
    `loopSize > 3` a `tstack2` terminal-mismatch increment is added, for `loopSize == 3` the
    closing-A·T penalty (`atPenalty{H,S}`). Then **for loopSize == 3 the triloop bonus** and
    **for loopSize == 4 the tetraloop bonus** (bsearch, key = `numSeq1 + i` starting at the
    closing-pair 5′ base → full loop string incl. closing pair) are **added** to ΔH°/ΔS°. ΔH in
    cal/mol, ΔS in cal/(K·mol); negative ΔH stabilises. Loop-length tables `loops.dh`/`loops.ds`
    (hairpin column): ΔH° = 0 all sizes, ΔS° by size (3=4=−11.28, 5=−10.64, … 30=−20.31).
    Unimolecular init (`type==4`): `dplx_init_H=0`, `dplx_init_S=−1e-11`, `RC=0`;
    `Tm = mh/(ms + (N/2 − 1)·saltCorrection) − 273.15`, `saltCorrection = 0.368·ln(mv/1000)`,
    no strand-concentration term (intramolecular).
  - **SantaLucia & Hicks (2004)** (rank 1) — the model source (Eqs. 8–9): length-3 hairpins add
    a triloop bonus + closing-A·T penalty; length-4 add a tetraloop bonus + terminal mismatch.
    Article body confirms the model; the **numeric bonus tables** come from the libprimer3 files
    (the same values the article's supplementary material tabulates).
  - **primer3-py 2.3.0 `calc_hairpin`** (rank 3) — ground-truth numbers at mv=50, dv=0, dntp=0,
    dna_conc=50 nM (see datasets below).
- **Datasets (primer3-py `calc_hairpin` oracles):**
  | Sequence | loop | dH (cal/mol) | dS (cal/K/mol) | Tm (°C) |
  |----------|------|-------------:|---------------:|--------:|
  | GGGGCGAAAGCCCC | tetraloop CGAAAG | −40900 | −114.1872884299936 | 85.03347700825856 |
  | GGGGGGGACCCCC | tetraloop GGGGAC | −34000 | −94.1872884299836 | 87.8328944728006 |
  | GGGCGAAGCCC | triloop CGAAG | −27800 | −77.68485895331574 | 84.7060915802943 |
  | GGGGGAACCCC | triloop GGAAC | −26000 | −73.18485895331571 | 82.11474153055735 |
  | GGGCTTTTGCCC | non-special 4-nt | −32400 | −94.58485895332572 | 69.39954078842845 |
  | GGGCAAAAAGCCC | non-special 5-nt | −30100 | −87.74485895332572 | 69.89004085311882 |
  | AAAAAAAAAAAA / GCGC | (none) | structure_found = False | | |
- **Corner cases / failure modes:**
  - **Only 3-nt and 4-nt loops get a special bonus** (loopSize == 3 / == 4). A ≥ 5-nt loop never
    gets one; a non-special 3/4-nt loop falls through the bsearch with no match → no bonus.
  - **Key includes the closing pair** — the right loop nt but a different closing base pair is
    NOT in the table → no bonus.
  - **No structure** (homopolymer / too short) → ntthal `no_structure` (`structure_found = False`).

## Deviations and assumptions

- **ASSUMPTION — none correctness-affecting.** The bonus tables, their keying (full loop string
  incl. closing pair), the cal/mol-ΔH / cal-K-mol-ΔS convention, and their additive application
  are taken **verbatim** from the libprimer3 config files + `thal.c`, and cross-verified to
  machine precision against primer3-py `calc_hairpin`. No value is invented or recalled.
- **Bundled vs opt-in:** this unit *bundles* the tables inside `NtthalHairpin` so the special
  bonus is applied automatically for recognised loops; the legacy Table-4 `FindMostStableHairpin`
  and its caller-supplied `loopBonusDeltaG37` path remain unchanged (no-regression guarantee).
- No contradictions among sources.

Recommended coverage (MUST): a recognised **tetraloop** CGAAAG (`GGGGCGAAAGCCCC`) reproduces
primer3 ΔH/ΔS/ΔG/Tm with the −1100 bonus present; a recognised **triloop** CGAAG (`GGGCGAAGCCC`)
reproduces primer3 with the −2000 bonus; a **non-special** 4-nt loop TTTT is unchanged
(regression) and still matches primer3; a homopolymer / too-short input → null (ntthal
no_structure). SHOULD: a second tetraloop (GGGGAC) and triloop (GGAAC) exercise a second table
row each; the legacy `FindMostStableHairpin` + caller-supplied `loopBonusDeltaG37` path stays
unchanged.
