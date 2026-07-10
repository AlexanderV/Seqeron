---
type: concept
title: "Primer-dimer thermodynamics and nearest-neighbour melting temperature (Tm)"
tags: [primer, algorithm, validation]
mcp_tools:
  - primer_dimer
  - primer_melting_temperature
  - primer_melting_temperature_salt
sources:
  - docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
  - docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md
  - docs/Evidence/PRIMER-TM-001-NN-Evidence.md
  - docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md
  - docs/Evidence/SEQ-THERMO-001-Evidence.md
  - docs/Evidence/SEQ-TM-001-Evidence.md
source_commit: 52c02ee8f4642a46e7ab17988a729a45ffbe5268
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
    - predicate: relates_to
      object: concept:test-unit-registry
      source: primer-tm-001-hairpin-evidence
      evidence: "Test Unit ID: PRIMER-TM-001 (hairpin / secondary-structure Tm extension) ... Algorithm: DNA self-folding hairpin (stem + loop) MFE detection + unimolecular hairpin Tm"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: primer-tm-001-nn-evidence
      evidence: "Test Unit ID: PRIMER-TM-001 ... Algorithm: Primer melting temperature — nearest-neighbour salt-corrected Tm (opt-in); adds the ΔH°/ΔS° → Tm path and the salt corrections"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: primer-tm-001-special-loop-evidence
      evidence: "Test Unit ID: PRIMER-TM-001 (special tri/tetraloop hairpin bonus tables — bundled & applied) ... bundled sequence-specific special triloop / tetraloop stability bonuses"
      confidence: high
      status: current
---

# Primer-dimer thermodynamics and nearest-neighbour melting temperature (Tm)

Computing the **melting temperature (Tm)** of a DNA duplex — and, for PCR primer QC, the Tm
of the most stable **self-dimer / hetero-dimer** two primers can form — from the
**SantaLucia unified nearest-neighbour (NN) thermodynamic model**. This is the **first
ingested unit of the PRIMER (PCR primer-design / MolTools) family** (test unit
**PRIMER-TM-001**). The literature-traced records are [[primer-tm-001-dimer-evidence]] (the
intermolecular dimer), [[primer-tm-001-hairpin-evidence]] (the intramolecular hairpin
extension, synthesized in the [[#Intramolecular hairpin self-folding|hairpin section]] below),
and [[primer-tm-001-nn-evidence]] (the **per-oligo design Tm** path plus the salt-correction
methods, synthesized in the [[#Per-oligo design Tm and salt corrections|per-oligo section]]
below). [[primer-tm-001-special-loop-evidence]] adds the **bundled special tri/tetraloop
hairpin bonus tables** — the increment the hairpin section flags — now applied automatically
inside the ntthal hairpin DP. [[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes
the artifact pattern. The **base** PRIMER-TM-001 unit is a separate algorithm —
[[primer3-weighted-penalty-objective]], the weighted per-primer selection penalty — which
*consumes* a per-primer Tm (and self-/dimer-alignment scores) as input terms.

Unlike the **RNA** nearest-neighbour folding used in [[rna-base-pairing]] and
[[pre-mirna-hairpin-detection]] (Turner 2004 parameters, single-strand secondary structure),
this unit uses the **DNA** SantaLucia & Hicks (2004) parameters and models an
**intermolecular** duplex between two oligos at a finite strand concentration.

A **third** DNA-duplex Tm surface — the SEQ-family
[[dna-duplex-nearest-neighbor-thermodynamics]] (SEQ-THERMO-001) — computes the same SantaLucia
bimolecular Tm but over the **1997 DNA_NN3** (Allawi & SantaLucia) table with **per-terminus
`init_A/T` (2.3, 4.1) / `init_G/C` (0.1, −2.8)** counted off the terminal bases, rather than this
engine's **2004 unified** table with a fixed duplex-initiation (+0.2 / −5.7) plus a terminal-A·T
penalty (+2.2 / +6.9). That engine returns the whole **ΔH°/ΔS°/ΔG°** tuple (not just a Tm) and does
not model dimers/hairpins; the **salt correction is identical** (`0.368·(N/2)·ln[Na⁺]`, N/2 = NN
stacks = `len − 1`).

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

## Per-oligo design Tm and salt corrections

The simplest use of this model is the **per-oligo (single-strand) design Tm** — a primer's own
Tm as one input to the [[primer3-weighted-penalty-objective|penalty objective]] — computed by
summing the NN stacks of the primer against its own complement and applying Eq. 3 directly
(opt-in NN Tm; the legacy Wallace / Marmur-Doty scalar Tm — the SEQ-family
[[melting-temperature]], SEQ-TM-001 — is the unchanged default). The published
oracle is SantaLucia & Hicks (2004) p.419: ΔH°=−43.5, ΔS°=−122.5, 0.2 mM each strand →
**Tm = 35.8 °C**. Hand-derived Table-1 sums confirm it term-by-term (`GCGCGC` self x=1 →
Tm 35.0473 °C; `ATGCATGC` x=4 → 30.4338 °C; `CGCGAATTCGCG` self → 61.1452 °C at 1 M Na, 0.5 µM).

Beyond the SantaLucia Eq. 5 entropy correction, two **Owczarzy 1/Tm corrections** are ported
verbatim from Biopython `salt_correction` (methods 5–7):

- **Monovalent, Owczarzy 2004 (method 6):**
  `1/Tm[Na] = 1/Tm[1 M] + (4.29·f(GC) − 3.95)·1e−5·ln[Na⁺] + 9.40e−6·(ln[Na⁺])²` — the Tm-vs-Na⁺
  relationship is **quadratic** (non-linear) over 69 mM–1.02 M, GC-dependent via f(GC).
- **Divalent Mg²⁺/dNTP, Owczarzy 2008 (method 7):**
  `1/Tm = 1/Tm[1 M] + a + b·ln[Mg] + f(GC)·(c + d·ln[Mg]) + (1/(2(N−1)))·(e + f·ln[Mg] + g·ln[Mg]²)`;
  in the mixed regime (0.22 ≤ R < 6, R = √[Mg]/[Mon]) a, d, g are reparameterised, and free Mg²⁺
  is reduced by **dNTP chelation** (Ka = 3×10⁴). With no Mg²⁺ it falls back to the monovalent path.

The same NN engine also scores **internal single mismatches and single dangling ends** on a
per-oligo duplex (`CalculateNearestNeighborThermodynamicsMismatch`). The parameter tables are
Biopython transcriptions of the primaries — `NnInternalMismatch` from `DNA_IMM`
(Allawi & SantaLucia 1997/1998, Peyret 1999; the WC-placeholder and inosine entries dropped)
and `NnDanglingEnd` from `DNA_DE` (Bommarito 2000) — verified against SantaLucia & Hicks (2004)
Tables 2 (mismatch worked example ΔG°37 −8.32 vs −8.349 within rounding) and 3 (all 32
dangling-end ΔH° reproduced exactly). The C# path mirrors the `Tm_NN` convention exactly: the
bottom strand is the **plain complement** (not reverse-complement), each NN/mismatch key is
looked up forward then character-reversed, dangling ends strip the outer column, and the
terminal-A·T count uses the un-dotted top strand. A **non-ACGT** base has no NN parameter →
thermodynamics not computable → **null / NaN**. Default C_T = 0.5 µM (caller-overridable) and
Eq. 5's phosphate count N = 2·(length − 1) are the two documented assumptions.

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

## Intramolecular hairpin self-folding

The same SantaLucia thermodynamics power a **unimolecular** structure: a single primer folding
back on itself into a **stem + loop hairpin** (ntthal hairpin mode). QC asks for the most stable
self-fold (MFE) and its Tm. The energy model reuses the identical **`NnUnifiedParams` stem NN
stacks** (SantaLucia & Hicks 2004 Table 1) and adds a **hairpin-loop initiation** ΔG°37 read
from Table 4 by loop size (3→3.5, 4→3.5, 5→3.3, 6→4.0 … 30→6.3 kcal/mol), with **loop ΔH° = 0**
and **loop ΔS° = −ΔG°37·1000/310.15** (footnote a). Loops beyond the largest tabulated size use
the **Jacobson-Stockmayer** extrapolation `ΔG°37(n) = ΔG°37(x) + 2.44·R·310.15·ln(n/x)`.

Two features distinguish the hairpin from the dimer:

- **No bimolecular initiation, no terminal-A·T penalty.** The hairpin is stem NN stacks + loop
  only; the +0.2/−5.7 duplex-initiation term (a two-strand nucleation cost) and the Eq. 3
  terminal-A·T penalty are **excluded** — the loop-initiation term is the unimolecular nucleation
  cost. (UNAFold-consistent; the exact, sourced core.)
- **Concentration-independent Tm (Eq. 11):** `Tm = ΔH°·1000/ΔS° − 273.15`, carrying **no**
  `R·ln(C_T/x)` term — a unimolecular transition, confirmed by Vallone & Benight (1999) over
  0.5–260 µM.

**Steric floor:** loops < 3 nt are prohibited. **Length-3 / length-4** loops additionally get a
sequence-specific triloop/tetraloop stability bonus (length-3 also a closing-A·T penalty; length-4
a terminal mismatch). Those bonus values are **not** in the SantaLucia article body — they are the
libprimer3 `triloop.*`/`tetraloop.*` config tables ([[primer-tm-001-special-loop-evidence]]). The
full `NtthalHairpin` port now **bundles** them (16 triloop rows, all ΔS = 0; 76 tetraloop rows)
and applies the bonus **automatically** in `calc_hairpin` for a recognised loop — the bsearch key
is the **full loop string including the closing pair**, ΔH in cal/mol, ΔS in cal/(K·mol), added to
the loop ΔH°/ΔS°. A ≥ 5-nt loop, or a 3/4-nt loop whose closing pair or bases are not in the table,
gets no bonus. Verified to machine precision against primer3-py 2.3.0 `calc_hairpin` (e.g. tetraloop
`GGGGCGAAAGCCCC` → Tm 85.033 °C; triloop `GGGCGAAGCCC` → Tm 84.706 °C). The **legacy** Table-4
`FindMostStableHairpin` with its caller-supplied `loopBonusDeltaG37` path is unchanged — the earlier
"opt-in hairpin bonus tables" capability is now also available bundled.

**Hairpin oracles** (hand-derived from Table 1 + Table 4): `GGGCTTTTGCCC` (4-bp GGGC stem, TTTT
loop) → ΔH° = −25.8, ΔS° = −75.48486216346927, ΔG°37 = −2.3883700000000054, **Tm = 68.6404 °C**;
`GGGCAAAAAGCCC` (5-nt loop, increment 3.3) → **Tm ≈ 71.585 °C**; poly-A → no hairpin (null) / NaN.

## Failure modes and contract

- **No stable dimer:** a homopolymer (poly-A) self-dimer forms no Watson-Crick duplex →
  `structure_found = False` → the method returns **null / NaN**.
- **Symmetry trap:** x = 1 applies only when **both** aligned oligos are palindromes; a
  self-dimer of a *non*-palindromic oligo still uses **x = 4**.
- **Invalid input** (null, < 2 bases, non-ACGT) → null / NaN.
- The tri/tetraloop **hairpin bonus tables** are a hairpin/monomer feature (see
  [[primer-tm-001-special-loop-evidence]]), not a dimer one — the dimer engine deliberately
  does not model them; the hairpin path bundles and applies them.

## Worked oracles

- `GCGCGCGC` self-dimer (x=1): 4·GC + 3·CG stacks + init, no A·T end → ΔH° = −70.8 kcal/mol,
  ΔS° = −192.6170, **Tm = 40.0906 °C** (hand-derived from Table 1 and primer3-parity).
- `TGCATGCATG`/`CATGCATGCA` (non-palindromic, x=4): ΔH° = −74.1, ΔS° = −211.8219,
  **Tm = 25.6596 °C**.
- `GCGCATGCGC` self (2×2 internal loop): ΔH = −84.4, **Tm = 43.1572 °C** (non-contiguous DP).
