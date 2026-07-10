---
type: concept
title: "DNA duplex nearest-neighbor thermodynamics (ΔH°/ΔS°/ΔG°/Tm)"
tags: [thermodynamics, sequence-statistics, validation]
sources:
  - docs/Evidence/SEQ-THERMO-001-Evidence.md
  - docs/Evidence/SEQ-TM-001-Evidence.md
source_commit: 52c02ee8f4642a46e7ab17988a729a45ffbe5268
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-thermo-001-evidence
      evidence: "Test Unit ID: SEQ-THERMO-001 ... Algorithm: DNA Duplex Thermodynamics (Nearest-Neighbor ΔH°/ΔS°/ΔG°/Tm)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: seq-thermo-001-evidence
      evidence: "Both compute the SantaLucia bimolecular nearest-neighbor Tm 'Tm = ΔH·1000/(ΔS + R·ln(C_T/x)) − 273.15'; this SEQ engine uses the DNA_NN3 / Allawi & SantaLucia 1997 table + per-terminus init_A/T (2.3,4.1) / init_G/C (0.1,−2.8), the primer engine the 2004 unified params + fixed duplex-init + terminal-A·T penalty."
      confidence: high
      status: current
---

# DNA duplex nearest-neighbor thermodynamics (ΔH°/ΔS°/ΔG°/Tm)

Computing the full **thermodynamic profile of a Watson-Crick DNA duplex** — enthalpy **ΔH°**,
entropy **ΔS°**, Gibbs free energy **ΔG°₃₇**, and melting temperature **Tm** — from the sequence,
using the **nearest-neighbor (NN) model** of duplex stability. This is the **SEQ-family full-tuple
thermodynamics engine** (test unit **SEQ-THERMO-001**, [[seq-thermo-001-evidence]]), a verbatim port
of Biopython `Bio.SeqUtils.MeltingTemp.Tm_NN` over the **DNA_NN3** parameter set (Allawi &
SantaLucia 1997). [[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]]
describes the artifact pattern.

## The nearest-neighbor model

A duplex's stability is the **sum of the ΔH°/ΔS° of its overlapping dinucleotide (NN) steps** —
each of the 10 distinct Watson-Crick stacks (`AA/TT`, `AT/TA`, `TA/AT`, `CA/GT`, `GT/CA`, `CT/GA`,
`GA/CT`, `CG/GC`, `GC/CG`, `GG/CC`) carries a tabulated (ΔH kcal/mol, ΔS cal/(mol·K)) pair
(DNA_NN3, Allawi & SantaLucia 1997) — plus a **two-terminus initiation** term. The initiation is
counted off the **two terminal bases**: `ends = seq[0] + seq[-1]`, then `init_A/T (2.3, 4.1)` is
added once per terminal A or T and `init_G/C (0.1, −2.8)` once per terminal G or C — so a duplex
with A/T at one end and G/C at the other picks up one of each. ΔG°₃₇ is the standard combination
`ΔG° = ΔH° − 310.15·ΔS°/1000` (kcal/mol).

The **melting temperature** is the SantaLucia bimolecular equation

```
Tm = (1000·ΔH°) / (ΔS° + R·ln k) − 273.15     R = 1.987 cal/(mol·K)
```

where `k = (dnac1 − dnac2/2)·1e−9`. With the default equimolar strands `dnac1 = dnac2 = 25` nM,
`k = 12.5` nM `= C_T/4` (C_T = 50 nM) — i.e. the **F = 4 symmetry factor** for two non-self-
complementary strands (MELTING's `-F`; a self-complementary duplex would use F = 1). ΔH°×1000
converts kcal→cal so the ratio's units match ΔS°.

A **[Na⁺] salt correction** (method 5, SantaLucia 1998) is folded into the entropy before Tm:
`ΔS° ← ΔS° + 0.368·(len − 1)·ln[Na⁺]` with [Na⁺] in mol/L (`len − 1` = the number of NN stacks).
Lower salt lowers ΔS° magnitude, so Tm is **monotonic in [Na⁺]** — more salt raises Tm.

## Worked oracles

- **Biopython `Tm_NN` docstring:** `CGTTCCAAAGATGTGGGCATGAGCTTAC` at dnac1 = dnac2 = 25 nM,
  Na = 50 mM → **Tm = 60.32 °C**.
- **`GCGC`** (Na = 0.05 M, C_T = 250 nM): init both ends G/C ⇒ ΔH 0.2 / ΔS −5.6; NN GC+CG+GC ⇒
  ΔH −30.2 / ΔS −76.0 ⇒ **ΔH° = −30.0**; salt ΔS 0.368×3×ln(0.05) = −3.307 ⇒ **ΔS° = −84.91**;
  **ΔG°₃₇ = −3.67**; **Tm = −18.6 °C** (a legitimate negative Tm for a 4-mer at low concentration).
- **`ATCG`** (one A/T end + one G/C end): **ΔH° = −23.6, ΔS° = −71.81** — the per-terminus init at work.

## Failure modes and contract

- **Length < 2 ⇒ `(0, 0, 0, 0)`** — the NN model has no dinucleotide step for an empty or single-base
  input, so no duplex is defined. This is the one documented API/edge-case convention (not a
  thermodynamic value); it alters nothing for valid length ≥ 2 input.
- **Case-insensitive** — input is `ToUpperInvariant`-folded before table lookup.
- **Non-self-complementary only** — the unit fixes the default **F = 4** case.

## Relationship to the other Tm surfaces

Three distinct melting-temperature surfaces coexist in the library; this page is the **full-tuple NN
thermodynamics** one:

- [[primer-dimer-thermodynamics-tm]] — the **PCR-primer** NN engine. Same SantaLucia bimolecular
  physics, but a **different parameter vintage**: the **2004 unified** SantaLucia & Hicks table
  (`NnUnifiedParams`) with a single fixed **duplex-initiation** (+0.2 / −5.7) plus a **terminal-A·T
  penalty** (+2.2 / +6.9), versus this unit's **1997 DNA_NN3** table with **per-terminus
  `init_A/T` / `init_G/C`** counted off the terminal bases. That engine also models dimers, hairpins
  and mismatches and exposes only a Tm/structure; SEQ-THERMO returns the whole **ΔH°/ΔS°/ΔG°** tuple.
  The salt term (`0.368·(N/2)·ln[Na⁺]`, N/2 = NN stacks = `len − 1`) is the **same** Eq. 5 in both.
- The **legacy Wallace / Marmur-Doty scalar Tm** — the SEQ-family [[melting-temperature]]
  (`CalculateMeltingTemperature`, SEQ-TM-001), also bundled by
  [[seq-summary-001-evidence|SummarizeNucleotideSequence]] (see [[base-composition]]) — is a `%GC`
  rule of thumb length-dispatched between `4(G+C)+2(A+T)` (short) and `64.9 + 41·(GC−16.4)/N`
  (len ≥ 14), **not** a thermodynamic model — it carries no ΔH°/ΔS°/ΔG°.

The RNA analogue — nearest-neighbor free energy for single-strand folding — is the Turner-parameter
model on [[rna-free-energy-turner-model]]; this DNA duplex engine and that RNA folding engine share
the NN-summation idea but use different alphabets, parameter sets, and (duplex vs. single-strand)
structural assumptions.
</content>
