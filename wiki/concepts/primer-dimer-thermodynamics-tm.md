---
type: concept
title: "Primer-dimer thermodynamics and nearest-neighbour melting temperature (Tm)"
tags: [primer, algorithm, validation]
sources:
  - docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
source_commit: 6c22f2b1ee758fb6d2b4c748c139d0193a4e313a
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: primer-tm-001-dimer-evidence
      evidence: "Test Unit ID: PRIMER-TM-001 (self-/hetero-dimer Tm extension) ... Algorithm: Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment"
      confidence: high
      status: current
---

# Primer-dimer thermodynamics and nearest-neighbour melting temperature (Tm)

Computing the **melting temperature (Tm)** of a DNA duplex — and, for PCR primer QC, the Tm
of the most stable **self-dimer / hetero-dimer** two primers can form — from the
**SantaLucia unified nearest-neighbour (NN) thermodynamic model**. This is the **first
ingested unit of the PRIMER (PCR primer-design / MolTools) family** (test unit
**PRIMER-TM-001**). The literature-traced record is [[primer-tm-001-dimer-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
artifact pattern.

Unlike the **RNA** nearest-neighbour folding used in [[rna-base-pairing]] and
[[pre-mirna-hairpin-detection]] (Turner 2004 parameters, single-strand secondary structure),
this unit uses the **DNA** SantaLucia & Hicks (2004) parameters and models an
**intermolecular** duplex between two oligos at a finite strand concentration.

## The nearest-neighbour Tm model

A duplex's stability is the sum of the ΔH°/ΔS° of its **10 distinct Watson-Crick NN stacks**
(SantaLucia & Hicks 2004, Table 1, at 1 M NaCl; in the repo as `NnUnifiedParams`), plus:

- **duplex initiation** ΔH° = +0.2 kcal/mol, ΔS° = −5.7 cal/(K·mol);
- a **terminal A·T penalty** ΔH° = +2.2, ΔS° = +6.9 for **each** A·T-closed helix end;
- a **symmetry** correction ΔS° = −1.4 for a self-complementary duplex.

The **bimolecular Tm** (Eq. 3) is

```
Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15    R = 1.9872 cal/(K·mol)
```

where `C_T` is total strand concentration and **x is a symmetry factor**: **x = 4** for a
non-self-complementary duplex, **x = 1** for a self-complementary (palindromic) one. In
Primer3's `thal.c` the divisor is `RC = R·ln(dna_conc/1e9)` (x=1) or `/4e9` (x=4), with
`dna_conc` in nM; a sequence is "symmetric" iff it is an **even-length reverse-complement
palindrome** (`symmetry_thermo`). Default conditions: `dna_conc = 50` nM, `temp = 310.15` K.

Entropy carries a **[Na⁺] salt correction** (Eq. 5):
`ΔS°[Na⁺] = ΔS°[1 M] + 0.368·(N/2)·ln[Na⁺]`, where N = total phosphates so N/2 = number of NN
stacks. In `thal.c` this is `saltCorrection = 0.368·ln((mv + 120·√max(0, dv − dntp))/1000)`
per stack. This makes Tm **monotonic in salt** — lower [Na⁺] gives lower Tm.

## Dimer analysis via ntthal thermodynamic alignment

Primer-dimer QC asks: over **all** ways two oligos can anneal, which duplex is the most
stable, and what is its Tm? Primer3 answers this with the `ntthal` **thermodynamic alignment**
(Untergasser et al. 2012) — a dynamic program over the SantaLucia NN model that finds the
maximum-Tm structure, exposed as `calcHomodimer` / `calcHeterodimer`. The Seqeron port
(`PrimerDesigner.CalculateDimerThermodynamicsNtthal`, delegated to by
`CalculateDimerMeltingTemperature` / `CalculateSelfDimerMeltingTemperature`) is a verbatim
port of `thal.c` and models the **full, possibly non-contiguous** optimum:

- an internal single **mismatch** (`stackmm`),
- an internal **loop** (`interior` + `tstack` + `ILAS·|Δn|`),
- a single/multi-base **bulge** (`bulge`),
- a **terminal overhang / dangling end** (`tstack2`, `dangle3`/`dangle5`).

`EnthalpyDPT(i,j)` / `EntropyDPT(i,j)` hold the best ΔH°/ΔS° of a duplex closing at WC pair
`(i,j)`; the best terminal pair over all `(i,j)` minimises ΔG, then N (stacks) from traceback
feeds the Tm formula above. The C# port reproduces primer3-py 2.3.0's `calc_homodimer` /
`calc_heterodimer` to machine precision on both contiguous and non-contiguous cases.

## Failure modes and contract

- **No stable dimer:** a homopolymer (poly-A) self-dimer forms no Watson-Crick duplex →
  `structure_found = False` → the method returns **null / NaN**.
- **Symmetry trap:** x = 1 applies only when **both** aligned oligos are palindromes; a
  self-dimer of a *non*-palindromic oligo still uses **x = 4**.
- **Invalid input** (null, < 2 bases, non-ACGT) → null / NaN.
- The only capability not modelled is the optional caller-supplied tri/tetraloop &
  terminal-mismatch **hairpin bonus tables** (a hairpin/monomer feature, not a dimer one).

## Worked oracles

- `GCGCGCGC` self-dimer (x=1): 4·GC + 3·CG stacks + init, no A·T end → ΔH° = −70.8 kcal/mol,
  ΔS° = −192.6170, **Tm = 40.0906 °C** (hand-derived from Table 1 and primer3-parity).
- `TGCATGCATG`/`CATGCATGCA` (non-palindromic, x=4): ΔH° = −74.1, ΔS° = −211.8219,
  **Tm = 25.6596 °C**.
- `GCGCATGCGC` self (2×2 internal loop): ΔH = −84.4, **Tm = 43.1572 °C** (non-contiguous DP).
