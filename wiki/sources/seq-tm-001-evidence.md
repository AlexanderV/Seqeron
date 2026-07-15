---
type: source
title: "Evidence: SEQ-TM-001 (CalculateMeltingTemperature — scalar Wallace-rule / Marmur-Doty GC-formula melting temperature, length-dispatched)"
tags: [validation, thermodynamics, sequence-statistics, primer]
doc_path: docs/Evidence/SEQ-TM-001-Evidence.md
sources:
  - docs/Evidence/SEQ-TM-001-Evidence.md
source_commit: 52c02ee8f4642a46e7ab17988a729a45ffbe5268
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-TM-001

The validation-evidence artifact for test unit **SEQ-TM-001** — the standalone scalar
melting temperature `CalculateMeltingTemperature`, which **length-dispatches** between two classic
**%GC rule-of-thumb** Tm formulas. Synthesized on the dedicated concept
[[melting-temperature]] (the third distinct Tm surface and the canonical home of the Wallace /
Marmur-Doty formulas); one instance of the templated [[algorithm-validation-evidence|evidence
artifact]] pattern; [[test-unit-registry]] tracks the unit.

## What this file records

- **Length-dispatch** at `ThermoConstants.WallaceMaxLength = 14`:
  - **Wallace rule** (Thein & Wallace 1986), short oligos: `Tm = 4·(G+C) + 2·(A+T)`. Biopython
    `Tm_Wallace('ACGTTGCAATGCCGTA')` → **48.0** (16 nt, 8 GC: `2·8 + 4·8`).
  - **Marmur-Doty GC formula** (Marmur & Doty 1962), len ≥ 14, repo `ThermoConstants` encoding
    `Tm = 64.9 + 41·(GC − 16.4)/N`. Oracle `GCGCGCGCGCATATATATAT` (GC 10, N 20) → **51.78**.
- **Online sources** (each rank 3 unless noted):
  - **Biopython `Bio.SeqUtils.MeltingTemp`** source + 1.76 API docs — `Tm_Wallace`
    (`4·(G+C)+2·(A+T)`, "rule of thumb for 14–20 nt primers"), `Tm_GC` general form
    `Tm = A + B·%GC − C/N + salt − D·%mismatch`; valueset 1 `69.3 + 0.41·%GC − 650/N`
    (Marmur & Doty 1962 / Chester & Marshak 1993), valueset 2 "QuikChange" `81.5 + 0.41·%GC − 675/N`
    (docstring `Tm_GC('CTGCTGATXGCACGAGGTTATGG', valueset=2)` → 69.20).
  - The same file also documents the **nearest-neighbor** `Tm_NN` (DNA_NN3, Allawi & SantaLucia
    1997), **SantaLucia (1998)** oligo Tm `Tm = ΔH°/(ΔS° + R·ln(C_T/x)) − 273.15` (rank 1) and its
    **salt-correction method 5** `0.368·(N−1)·ln[Na+]` — but those belong to the full-tuple NN
    engine (see the consolidation note below), not to this scalar unit.
- **Datasets (hand-derived, no library run):** Wallace `ACGTTGCAATGCCGTA` → 48.0; Marmur-Doty
  `GCGCGCGCGCATATATATAT` → 51.78; and (via the length-dispatch, cross-checked by SEQ-SUMMARY-001)
  16-mer `ATGCATGCATGCATGC` → 43.375, 8-mer `ATGCATGC` → 24.0.
- **Corner cases:** Wallace is a rule of thumb only for 14–20 nt; the NN model is undefined for
  < 2 nt; the NN Tm depends on C_T and the self-complementary symmetry factor x (1 vs 4).

## Consolidation note (why the NN content is not synthesized here)

The Change History records that **SEQ-TM-001 is a duplicate/consolidated Registry entry** for the
two melting-temperature methods already delivered under **SEQ-THERMO-001**:
`CalculateMeltingTemperature` (the scalar Wallace/Marmur unit → [[melting-temperature]]) **and**
`CalculateThermodynamics` (the full ΔH°/ΔS°/ΔG°/Tm engine →
[[dna-duplex-nearest-neighbor-thermodynamics]]). The evidence was independently re-retrieved and
the implementation verified conformant; the unit was **consolidated, not re-implemented** (TestSpec
§7), reusing the canonical fixture
`SequenceStatistics_CalculateThermodynamics_Tests.cs`. So this page's genuinely new wiki surface is
the **scalar %GC Tm**; its nearest-neighbor datasets cross-reference the existing NN concept.

## Deviations and assumptions

- **One documented ASSUMPTION — the 14-nt dispatch boundary.** `useWallaceRule = length < 14`
  (`ThermoConstants.WallaceMaxLength`); Biopython documents Wallace only loosely as "14–20 nt",
  fixing no exact switch point. Non-correctness-affecting when the caller passes an explicit flag.
- **One documented ASSUMPTION (NN facet) — the default C_T.** The repo `CalculateThermodynamics`
  default `primerConcentration = 2.5e-7` (250 nM) vs Biopython's 50 nM; the **formula is identical**
  and `primerConcentration: 5e-8` reproduces Biopython's 60.32 exactly. A documented default
  parameter, not an invented constant.
- **Formula-variant note:** the repo Marmur-Doty constants (64.9 / −672.4) differ from Biopython
  `Tm_GC` valueset 1 (69.3 / −650); same `A + B·%GC − C/N` family, different published constants.
- **No source contradictions** — Biopython, SantaLucia 1998, and Marmur & Doty agree within the
  documented variant/parameter choices.

Recommended coverage (from the artifact): MUST — Wallace `ACGTTGCAATGCCGTA` → 48.0; Marmur-Doty
`GCGCGCGCGCATATATATAT` → 51.78; the NN oracles (60.32; DNA_NN3 ΔH/ΔS/ΔG tuple; empty/len-1 → all
zero) — the last set are the SEQ-THERMO-001 facet. SHOULD — higher [Na+] raises NN Tm;
case-insensitivity and Watson-Crick symmetry (AA == TT).
</content>
