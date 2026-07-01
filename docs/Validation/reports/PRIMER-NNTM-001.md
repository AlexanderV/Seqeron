# Validation Report: PRIMER-NNTM-001 — Nearest-Neighbour Salt/Mismatch/Dangling-End Tm

- **Validated:** 2026-06-25   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.CalculateMeltingTemperatureNN`,
  `PrimerDesigner.CalculateMeltingTemperatureNNMismatch`
  (helpers: `CalculateNearestNeighborThermodynamics`, `CalculateNearestNeighborThermodynamicsMismatch`,
  `ApplyOwczarzy2004`, `ApplyOwczarzy2008`)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs` (lines ~478–953, 1715–1824)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_NearestNeighborTm_Tests.cs` (25 tests, all green)

## Authoritative sources opened this session
- **SantaLucia & Hicks (2004)** *Annu Rev Biophys* 33:415–440 — Table 1 (unified NN ΔH°/ΔS°,
  init +0.2/−5.7, terminal-A·T +2.2/+6.9, symmetry 0/−1.4), Eq. 3 (bimolecular Tm), Eq. 5
  (entropy salt correction). **These are the constants the code actually uses.**
- **SantaLucia (1998)** *PNAS* 95(4):1460 — original unified set (Biopython `DNA_NN3`); differs
  from the 2004 review in the *initiation* convention (init A/T = 2.3/4.1, init G/C = 0.1/−2.8,
  AA = −7.9/−22.2). The code uses the **2004** revision, not the 1998 one (see Findings).
- **Allawi & SantaLucia (1997/1998)** + **Peyret et al. (1999)** — internal single-mismatch NN
  (Biopython `DNA_IMM1`).
- **Bommarito, Peyret & SantaLucia (2000)** *NAR* 28:1929 — single dangling-end NN (Biopython `DNA_DE1`).
- **Owczarzy et al. (2004)** *Biochemistry* 43:3537 — monovalent (Na⁺) 1/Tm correction (Biopython method 6).
- **Owczarzy et al. (2008)** *Biochemistry* 47:5336 — divalent (Mg²⁺/dNTP) correction (Biopython method 7).

## Stage A — Description

**Formula.** Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15, R = 1.9872 cal/(K·mol),
x = 4 (non-self-complementary) / 1 (self-complementary). This is SantaLucia & Hicks (2004) Eq. 3
verbatim. ΔH°/ΔS° = duplex init + Σ NN stacks + per-A·T-terminus penalty + (self-comp) symmetry —
matches Eq. 1 + Table 1.

**Salt corrections.**
- *SantaLucia entropy (Eq. 5):* ΔS°[Na] = ΔS° + 0.368·(N/2)·ln[Na⁺], N = 2·(L−1) phosphates. Correct.
- *Owczarzy 2004 monovalent:* 1/Tm[Na] = 1/Tm[1M] + (4.29e-5·fGC − 3.95e-5)·ln[Na] + 9.40e-6·(ln[Na])².
  Matches the paper's Eq. and Biopython method 6.
- *Owczarzy 2008 divalent:* full a..g model with the R = √[Mg²⁺]/[Mon] regime split (R<0.22 monovalent,
  0.22≤R<6 mixed a/d/g reparameterisation, R≥6 divalent) and the dNTP·Mg²⁺ chelation quadratic.
  Matches Biopython method 7.

**Edge-case semantics.** len<2 / non-ACGT / empty / null → not-computable (null thermo, NaN Tm);
unequal-length strands → null; a stack with two adjacent mismatches (no NN parameter) → null. All defined.

**Independent cross-check (hand computation, exact numbers).**
| Quantity | Hand-derived (SantaLucia&Hicks 2004 Table 1 + Eq. 3) | Test expects |
|---|---|---|
| ATGCATGC ΔH°/ΔS° | −57.1 / −156.5 | −57.1 / −156.5 ✓ |
| GCGCGC (self-comp) ΔH°/ΔS° | −50.4 / −134.7 | −50.4 / −134.7 ✓ |
| GCGCGC Tm (no salt, x=1) | 35.04730599 °C | 35.0473059911 ✓ |
| ATGCATGC Tm (no salt, x=4) | 30.43380607 °C | 30.4338060665 ✓ |
| Worked example ΔH=−43.5, ΔS=−122.5, 0.2 mM | 35.79 °C | 35.8 ✓ |
| MM1 (CGTGAC/GCGCTG) ΔH/ΔS, Tm | −35.5 / −101.5, −6.40608793 °C | −35.5 / −101.5, −6.4060879279 ✓ |
| DE1 (AGCGCGC/.CGCGCG) ΔH/ΔS, Tm | −51.9 / −136.4, 35.80349218 °C | −51.9 / −136.4, 35.8034921829 ✓ |

**Findings.** The file header comment loosely labels the table "SantaLucia (1998) Table 1 (unified)",
but the constants used are the **SantaLucia & Hicks (2004)** revision (Biopython `DNA_NN4`), which the
per-constant comments correctly cite. The 2004 review is the authoritative restatement of the unified
parameters; using it is correct, only the one summary label is loose. Documented, not a defect → Stage A PASS.

## Stage B — Implementation

**Code realises the formula.** `CalculateNearestNeighborThermodynamics` sums init + NN stacks +
terminal-A·T (per A/T terminus) + symmetry (self-comp). `CalculateMeltingTemperatureNN` applies Eq. 3
with x from self-comp, then the selected salt correction. The mismatch/dangling path strips dangling
columns first (DNA_DE), then sums WC or internal-mismatch stacks (DNA_IMM, forward-then-reverse key),
mirroring Biopython `Tm_NN(imm_table=DNA_IMM, de_table=DNA_DE)`.

**Parameter tables — verbatim audit vs Biopython.**
- `NnUnifiedParams` = `DNA_NN4` (16 dinucleotides) — exact match.
- `NnInternalMismatch` (50 entries) = `DNA_IMM1` ACGT subset — **all 50 match** (programmatic diff: none).
- `NnDanglingEnd` (32 entries) = `DNA_DE1` — **all 32 match** (programmatic diff: none).

**Cross-verification table recomputed vs code (run) and vs external oracles.**
| Case | C# (test, run green) | Independent re-derivation | Biopython `Tm_NN` (DNA_NN4) | primer3 `calc_tm` |
|---|---|---|---|---|
| GCGCGC no-salt | 35.0473 | 35.0473 | 35.0528 (R=1.987) | — |
| ATGCATGC no-salt | 30.4338 | 30.4338 | 30.4389 (R=1.987) | — |
| ATGCATGC Owczarzy2004 50 mM | 18.1900 | 18.1900 | 18.1947 (m6) | 18.6227 (1998 tbl) |
| GCGCGC Owczarzy2004 50 mM | 28.1593 | 28.1593 | 28.1645 (m6) | — |
| GCGCGC SantaLucia-entropy 50 mM | 24.9977 | 24.9977 | — | — |
| EcoRI CGCGAATTCGCG no-salt | 61.1452 | 61.1452 | — | 61.2532 (Na=1M) |
| CGCGAATTCGCG Owczarzy2008 Na50/Mg3 mM | (C2: >no-Mg) | 55.4498 | 55.4529 (m7, mixed regime R=1.10) | — |
| MM1 (internal G·T) Tm | −6.4061 | −6.4061 | −6.3997 (R=1.987) | — |
| DE1 (5′-dangling A) Tm | 35.8035 | 35.8035 | (Biopython mishandles this input) | — |

**Tool divergence reconciled.** Biopython uses **R = 1.987**; the C# code (and SantaLucia & Hicks 2004
Eq. 3) uses **R = 1.9872**. This fully explains the uniform ~0.005–0.006 °C C#-vs-Biopython gap; the C#
value is the one matching the paper's stated constant. primer3 uses the **SantaLucia 1998** table and its
own salt model, so primer3 differs more (≤0.5 °C at 50 mM, 0.1 °C at 1 M) — expected, within the spec's
±0.5 °C oracle band. The C# code declares it follows SantaLucia & Hicks (2004) + Biopython `DNA_NN4`/method
6–7, and is validated against THAT source's exact values (all match to 1e-8 by hand).

**Variant consistency.** EQ1/EQ2 confirm a fully-paired duplex through the mismatch path equals the
perfect-match path exactly (ΔH°, ΔS°, IsSelfComplementary, Tm). C1 confirms divalent-with-no-Mg reduces
to monovalent. M8 confirms the default `saltMode` is `Owczarzy2004Monovalent`.

**Numerical robustness.** No div-by-zero on stated ranges; invalid inputs short-circuit to null/NaN.

**Test quality audit.** Expected values are hand-derived from the 2004 table + Eq. 3 and independently
reproduced here; they are not code echoes. No skips, no widened tolerances (exact 1e-8 / 1e-9 on numeric
cases), no weakened assertions. Monotonicity/ordering tests (S5, C2, MM2) are genuine physical invariants.
Coverage spans every public method and overload, both perfect-match and mismatch/dangling paths, all four
salt modes, self-/non-self-complementary, very short (6-mer) and 12-mer, varying [Na⁺]/[Mg²⁺], and every
invalid-input class. No green-washing found.

**Findings / defects.** None.

## Verdict & follow-ups
- **Stage A: ✅ PASS** (one documented label looseness "1998" vs the 2004 constants actually used — the
  body comments cite 2004 correctly; numerically correct).
- **Stage B: ✅ PASS** — NN/IMM/DE tables verbatim vs Biopython; all worked values reproduce by hand and
  vs Biopython/primer3 within the reconciled R-constant / table-version differences.
- **State: ✅ CLEAN.** Full unfiltered `dotnet test Seqeron.sln -c Debug`: Failed 0
  (Seqeron.Genomics.Tests 18737 passed). No defect logged.
