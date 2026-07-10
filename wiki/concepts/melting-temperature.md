---
type: concept
title: "Scalar melting temperature (Wallace rule / Marmur-Doty GC formula)"
tags: [thermodynamics, sequence-statistics, primer, validation]
sources:
  - docs/Evidence/SEQ-TM-001-Evidence.md
source_commit: 52c02ee8f4642a46e7ab17988a729a45ffbe5268
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-tm-001-evidence
      evidence: "Test Unit ID: SEQ-TM-001 ... Algorithm: Melting Temperature (Wallace rule / Marmur-Doty GC formula)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:dna-duplex-nearest-neighbor-thermodynamics
      source: seq-tm-001-evidence
      evidence: "SEQ-TM-001 records the legacy %GC scalar Tm (Wallace 4(G+C)+2(A+T) / Marmur-Doty 64.9+41·(GC−16.4)/N), a rule-of-thumb carrying no ΔH°/ΔS°/ΔG°; the same evidence file's nearest-neighbor Tm datasets belong to the full-tuple SantaLucia engine SEQ-THERMO-001. Two Tm formulas for the same duplex, thermodynamic vs. %GC rule of thumb."
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:primer-dimer-thermodynamics-tm
      source: seq-tm-001-evidence
      evidence: "The Wallace / Marmur-Doty scalar Tm is the unchanged legacy default of the PCR-primer per-oligo design Tm; the opt-in NN Tm on primer-dimer-thermodynamics-tm is the thermodynamic alternative to this %GC rule of thumb."
      confidence: high
      status: current
---

# Scalar melting temperature (Wallace rule / Marmur-Doty GC formula)

The **standalone scalar melting temperature** `CalculateMeltingTemperature` (test unit
**SEQ-TM-001**, [[seq-tm-001-evidence]]) — a single **°C** number for a DNA oligo, computed by a
**length-dispatch** between two classic **%GC rule-of-thumb** formulas. This is the **third
distinct Tm surface** in the library and the **canonical home of the Wallace and Marmur-Doty
formulas** (previously only sketched as the "legacy default" on
[[primer-dimer-thermodynamics-tm]]). [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

Unlike the two nearest-neighbor engines, this is **not a thermodynamic model** — it carries **no
ΔH°/ΔS°/ΔG°**, no strand-concentration term and no salt correction. It is the fast approximation
that [[seq-summary-001-evidence|SummarizeNucleotideSequence]]'s `MeltingTemperature` field and the
primer engine's *default* per-oligo Tm both use.

## The length-dispatch and the two formulas

`CalculateMeltingTemperature` selects a formula by oligo length against
`ThermoConstants.WallaceMaxLength = 14`:

- **Short oligos (length < 14 nt) → Wallace rule** (Thein & Wallace 1986; Biopython
  `Tm_Wallace`):

  ```
  Tm = 4·(G + C) + 2·(A + T)
  ```

  A rule of thumb documented for ~14–20 nt primers. Worked oracle
  `ACGTTGCAATGCCGTA` (16 nt, 8 GC) → `2·8 + 4·8 = 48.0 °C` (Biopython `Tm_Wallace` docstring);
  8-mer `ATGCATGC` (GC 4) → `2·4 + 4·4 = 24.0 °C`.

- **Longer oligos (length ≥ 14 nt) → Marmur-Doty GC formula** (Marmur & Doty 1962), in the repo
  `ThermoConstants` encoding:

  ```
  Tm = 64.9 + 41·(GC − 16.4) / N          (equivalently 64.9 + 0.41·%GC − 672.4/N)
  ```

  Worked oracles `GCGCGCGCGCATATATATAT` (GC 10, N 20) → `64.9 + 41·(10 − 16.4)/20 = 51.78 °C`;
  16-mer `ATGCATGCATGCATGC` (GC 8) → `64.9 + 41·(8 − 16.4)/16 = 43.375 °C` (the length-dispatch
  case verified by [[seq-summary-001-evidence]]).

> **Formula-variant note.** Biopython's `Tm_GC` valueset 1 cites the sibling Marmur-Doty form
> `Tm = 69.3 + 0.41·%GC − 650/N` (and a "QuikChange" valueset 2 `81.5 + 0.41·%GC − 675/N`); the
> repository uses the **64.9 / −672.4** constants, a different published Marmur-Doty variant. The
> family and the `A + B·%GC − C/N` shape are the same; only the constants differ. The evidence
> hand-derivations use the repo's 64.9 form.

## Failure modes and contract

- **Length-dispatch boundary is 14 nt** (`length < 14` ⇒ Wallace, else GC formula). Biopython
  documents Wallace only loosely as a "rule of thumb for 14–20 nt"; the exact 14-nt switch point
  is the SEQ-TM-001 convention (`WallaceMaxLength`), the one documented **assumption** on this
  unit — non-correctness-affecting for callers that pass an explicit `useWallaceRule` flag.
- **empty / null ⇒ 0** (Wallace of an empty count is `0`).
- **case-insensitive** — counts fold `ToUpperInvariant` before the A/T/G/C tally.

## Relationship to the other Tm surfaces

Three distinct melting-temperature surfaces coexist; this page is the **legacy %GC scalar** one.
The evidence file `SEQ-TM-001-Evidence.md` is itself a **consolidated/duplicate registry entry**:
its Change History records that SEQ-TM-001 covers the two methods already delivered under
**SEQ-THERMO-001** — the scalar `CalculateMeltingTemperature` (this page) **and** the full-tuple
`CalculateThermodynamics` (the NN engine below) — so its nearest-neighbor datasets belong to that
other concept, not here.

- [[dna-duplex-nearest-neighbor-thermodynamics]] — the SEQ-family **full-tuple NN engine**
  (SEQ-THERMO-001): returns ΔH°/ΔS°/ΔG°/Tm over the 1997 DNA_NN3 table. The **thermodynamic**
  counterpart of this %GC rule of thumb; the two are the two Tm formulas
  `SequenceStatistics` exposes for the same duplex.
- [[primer-dimer-thermodynamics-tm]] — the **PCR-primer** NN engine (PRIMER-TM-001): the Wallace /
  Marmur-Doty scalar here is that engine's **unchanged legacy default** per-oligo Tm, with its
  2004-unified salt-corrected NN Tm as the opt-in thermodynamic alternative.

The RNA analogue (single-strand folding free energy) is the Turner-parameter model on
[[rna-free-energy-turner-model]] — a different alphabet, parameters, and structural assumption.
</content>
</invoke>
