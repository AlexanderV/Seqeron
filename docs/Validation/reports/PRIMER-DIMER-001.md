# Validation Report: PRIMER-DIMER-001 — ntthal Self/Hetero-Dimer Tm

- **Validated:** 2026-06-25   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.FindMostStableDimer`, `CalculateDimerMeltingTemperature`,
  `CalculateSelfDimerMeltingTemperature` (plus `CalculateDimerThermodynamicsNtthal` and the internal
  `NtthalDimer.Run` thermodynamic-alignment engine).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
  (public surface: lines 1411–1657), `NtthalDimer.cs` (full ntthal DP engine).
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_DimerTm_Tests.cs` (27 tests).

## Authoritative sources opened this session
1. **SantaLucia J, Hicks D (2004)** *Annu Rev Biophys* 33:415–440 — unified nearest-neighbour
   ΔH°/ΔS° (Table 1), bimolecular Tm Eq. 3 (`Tm = ΔH°/(ΔS° + R·ln(C_T/x)) − 273.15`, x = 1 for
   self-complementary, x = 4 otherwise) and the entropy salt correction (0.368 coefficient, Eq. 5).
2. **Untergasser A et al. (2012)** *Nucleic Acids Res* 40:e115 — Primer3 2.0 / ntthal thermodynamic
   alignment for bimolecular duplexes (best alignment maximising stability; NN stacks + internal
   loop + dangling/terminal-mismatch tstack2 + bimolecular initiation).
3. **primer3 `thal.c` + `primer3_config/*.dh,*.ds`** (libprimer3, the parameter set ntthal loads):
   `fillMatrix`, `LSH`, `RSH`, `maxTM`, `calc_bulge_internal`, `traceback`, `calcDimer`; constants
   `R=1.9872`, `DPLX_INIT_H=200`, `DPLX_INIT_S=-5.7`, AT-penalty `(2200, 6.9)`, `RC = R·ln(dna_conc/x)`
   with `x = 1e9` (palindrome) else `4e9`, `saltCorrectS = 0.368·ln(mv/1000)` per stack.
4. **primer3-py 2.3.0** reference oracle — installed this session
   (`pip install primer3-py==2.3.0`); `calc_homodimer`/`calc_heterodimer` at the library default
   conditions (mv=50 mM, dv=0, dntp=0, dna_conc=50 nM).

## Stage A — Description

- **Formula check.** The bimolecular Tm equation, the x=1/x=4 self-complementarity convention, the
  0.368 salt-correction coefficient, the unified NN init (+0.2 kcal ΔH / −5.7 cal ΔS) and the
  terminal A·T penalty (+2.2 kcal / +6.9 cal) all match SantaLucia & Hicks (2004) exactly. The full
  ntthal model (stacks + 1×1 internal mismatch via `stackmm`/Int2 + interior-loop length params +
  bulge-loop length params + `tstack`/`tstack2` terminal stacks + 5′/3′ dangles + the −300/310.15
  internal-loop asymmetry term) matches `thal.c` line-for-line.
- **Hand-derivation (independent of the code).** From SantaLucia & Hicks Table 1:
  - **GCGCGCGC self-dimer** (palindrome → x=1): ΔH° = 0.2 + 4·(−9.8) + 3·(−10.6) = **−70.8 kcal/mol**;
    ΔS° = −5.7 + 4·(−24.4) + 3·(−27.2) + 7·0.368·ln(0.05) = **−192.61700633667505 cal/(K·mol)**;
    ΔG°37 = **−11.059835484680235 kcal/mol**; Tm = **40.09064476882935 °C**.
  - **TGCATGCATG / CATGCATGCA hetero-dimer** (x=4, one 5′-T A·T end): ΔH° = **−74.1 kcal/mol**;
    ΔS° = **−211.8218652900108**; Tm = **25.659587124835923 °C**.
  Both reproduce the test fixture's hand-derived assertions to 1e-9 and primer3-py exactly.
- **Edge-case semantics (sourced).** No complementary base pair ⇒ ntthal `structure_found=False`
  ⇒ `null` / NaN (poly-A self, poly-G vs poly-A). Self-dimer of S ≡ hetero-dimer(S, S). x=1 only
  when *both* strands are reverse-complement palindromes (ntthal `symmetry_thermo`).
- **Verdict:** PASS — the description is biologically and mathematically correct and matches the
  cited primary sources and the ntthal reference algorithm.

## Stage B — Implementation

- **Code path.** `NtthalDimer.Run` is a faithful port of `thal.c`: `initMatrix`/`fillMatrix`, the
  `LSH`/`RSH` terminal stacks, `maxTM` stack extension, `calc_bulge_internal` (bulge / 1×1 mismatch /
  general internal loop), best-terminal-pair scan, and `traceback` to count NN stacks. ntthal native
  cal/mol are converted to the library's kcal/mol convention for ΔH/ΔG in
  `CalculateDimerThermodynamicsNtthal`. `CalculateDimerMeltingTemperature` routes through the full DP;
  `CalculateSelfDimerMeltingTemperature` delegates with the same sequence twice. `FindMostStableDimer`
  is the legacy contiguous-WC-run scorer (≥2 bp) retained for the most-stable contiguous-duplex API.
- **Differential cross-check vs primer3-py 2.3.0 (mv=50, dv=0, dntp=0, dna=50 nM).** Every C#
  value reproduces primer3-py to machine precision (ΔH/ΔS/ΔG/Tm/base-pairs all identical):

  | Case | type | primer3-py Tm (°C) | C# Tm (°C) | ΔH (cal/mol) | ΔS (cal/K/mol) | Match |
  |------|------|------|------|------|------|------|
  | GCGCGCGC | homo | 40.09064 | 40.09064 | −70800 | −192.61701 | exact |
  | ACGTACGTACGT | homo | 37.62509 | 37.62509 | −92000 | −262.62672 | exact |
  | CGATCGATCG | homo | 29.66002 | 29.66002 | −78800 | −226.82187 | exact |
  | GCATGC | homo | 0.68589 | 0.68589 | −43600 | −125.81215 | exact |
  | GGGGCCCC | homo | 29.01503 | 29.01503 | −57600 | −157.21701 | exact |
  | GCGCATGCGC (2×2 internal loop) | homo | 43.15725 | 43.15725 | −84400 | −233.42187 | exact |
  | ATCGATCGATCG / CGATCGATCGAT | hetero | 32.61073 | 32.61073 | −92000 | −264.72672 | exact |
  | TGCATGCATG / CATGCATGCA | hetero | 25.65959 | 25.65959 | −74100 | −211.82187 | exact |
  | GCGCAAAGCGC / GCGCTTTGCGC (3×3 loop) | hetero | 41.88164 | 41.88164 | −92300 | −256.82429 | exact |
  | GCGCGCGC / GCGCAGCGC (1-base bulge) | hetero | 19.81250 | 19.81250 | −70800 | −205.50701 | exact |
  | GCGCACGCGC / GCGCTAGCGC (2×2 mixed) | hetero | 18.56037 | 18.56037 | −68400 | −198.31701 | exact |
  | GCGCGCAAAA / AAAAGCGCGC (terminal overhang) | hetero | 24.65474 | 24.65474 | −60000 | −165.31215 | exact |
  | AAAAAAAA (self) / GGGGGGGG·AAAAAAAA | both | not found | `null` | — | — | exact |

  The contiguous-WC optima match to machine precision (confirming the memory note that K1 verified
  this) **and** the internal-loop / bulge / terminal-overhang extension cases match exactly.
- **Edge cases in code (traced/run).** Self == hetero(S,S) (37.62509 for ACGTACGTACGT); 1-bp inputs
  rejected by the contiguous scorer (`MinDimerBasePairs=2` ⇒ `A`/`T` → null); identical seqs handled;
  null/empty/non-ACGT ⇒ null; no-complementarity ⇒ null/NaN. The full ntthal DP returns a (very
  negative-Tm, positive-ΔG) structure for forced terminal-only pairings exactly as primer3 does
  (e.g. AT/AT Tm = −217.606 in both), confirming faithful no-special-casing.
- **Variant consistency.** `CalculateSelfDimerMeltingTemperature(S) == CalculateDimerMeltingTemperature(S,S)`
  (test S1); `CalculateDimerMeltingTemperature` == `CalculateDimerThermodynamicsNtthal().TmCelsius`
  (test N7). Salt/concentration monotonicity (S2/S3) hold per Eq. 5 / R·ln(C_T/x).
- **Test quality audit.** 27 tests. Expected values trace to (a) hand-derivation from SantaLucia &
  Hicks Table 1 (asserted to 1e-9) and (b) primer3-py 2.3.0 (asserted to 1e-3, i.e. the recorded
  4-dp reference; the full-precision probe this session showed agreement far beyond that). No
  weakened assertions, no skips, no widened tolerances, no code-echoes. Every public method/overload
  and every Stage-A path (no complementarity, self==hetero, 1-bp, identical, internal loop, bulge,
  terminal overhang, invalid input) is covered.
- **Verdict:** PASS — the code faithfully realises the validated description.

## Verdict & follow-ups
- **Final: ✅ CLEAN.** No defect found. Stage A PASS, Stage B PASS. The C# ntthal engine reproduces
  primer3-py 2.3.0 to machine precision across contiguous, internal-loop, bulge and terminal-overhang
  optima; the hand-derived SantaLucia & Hicks values match to 1e-9. Full unfiltered
  `dotnet test Seqeron.sln -c Debug` green (Genomics.Tests: 18741 passed, 0 failed; 0 warnings on
  changed files). No follow-ups.
