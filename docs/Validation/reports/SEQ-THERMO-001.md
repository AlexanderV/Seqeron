# Validation Report: SEQ-THERMO-001 — DNA Duplex Thermodynamics (Nearest-Neighbor ΔH°/ΔS°/ΔG°/Tm)

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateThermodynamics(string, double, double)`; delegate `SequenceStatistics.CalculateMeltingTemperature(string, bool)` (Wallace / Marmur-Doty via `ThermoConstants`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (with retrieved numbers)

1. **Biopython `Bio.SeqUtils.MeltingTemp`** (master) — fetched raw source.
   - `DNA_NN3` (Allawi & SantaLucia 1997) verbatim: `init (0,0)`, `init_A/T (2.3, 4.1)`, `init_G/C (0.1, -2.8)`, `sym (0,-1.4)`; `AA/TT (-7.9,-22.2)`, `AT/TA (-7.2,-20.4)`, `TA/AT (-7.2,-21.3)`, `CA/GT (-8.5,-22.7)`, `GT/CA (-8.4,-22.4)`, `CT/GA (-7.8,-21.0)`, `GA/CT (-8.2,-22.2)`, `CG/GC (-10.6,-27.2)`, `GC/CG (-9.8,-24.4)`, `GG/CC (-8.0,-19.9)`. dH kcal/mol, dS cal/(mol·K).
   - `Tm_NN`: `k = (dnac1 - dnac2/2)·1e-9`; `R = 1.987`; `Tm = (1000·delta_h)/(delta_s + R·ln(k)) − 273.15`. Terminal init applied per terminal base (AT count × init_A/T, GC count × init_G/C).
   - Salt method 5 (default): `corr = 0.368·(N−1)·ln(mon)`, `mon = Na·1e-3` (mol/L).
   - Docstring example: `Tm_NN('CGTTCCAAAGATGTGGGCATGAGCTTAC')` = **60.32**.
2. **Wikipedia — Nucleic acid thermodynamics** (cites SantaLucia 1998 as primary). NN table in **kJ**: AA/TT −33.1/−92.9; AT/TA −30.1/−85.4; TA/AT −30.1/−89.1; CA/GT −35.6/−95.0; GT/CA −35.1/−93.7; CT/GA −32.6/−87.9; GA/CT −34.3/−92.9; CG/GC −44.4/−113.8; GC/CG −41.0/−102.1; GG/CC −33.5/−83.3. Terminal A/T 9.6 kJ / 17.2 J; Terminal G/C 0.4 kJ / −11.7 J. Tm = ΔH°/(ΔS° + R·ln(C_T/x)), x=4 non-self-complementary / x=1 self-complementary.
3. **Real Biopython 1.85** (installed and executed this session): `mt.Tm_NN(Seq('CGTTCCAAAGATGTGGGCATGAGCTTAC'), nn_table=DNA_NN3, dnac1=25, dnac2=25, Na=50)` = **60.32091949919129**.

### Formula check
- **NN table identity.** Converting the Wikipedia "SantaLucia 1998 unified" kJ values by ÷4.184 reproduces the kcal `DNA_NN3` table exactly: AA/TT −33.1/4.184 = −7.91 ≈ −7.9 (ΔS −92.9/4.184 = −22.2); CG −44.4/4.184 = −10.61 ≈ −10.6 (−113.8/4.184 = −27.2); GC −41.0/4.184 = −9.80 (−102.1/4.184 = −24.4). The implementation's `NearestNeighborParams` matches this set value-for-value, including all Watson-Crick mirror keys (AA=TT, CA=TG, GT=AC, CT=AG, GA=TC, GG=CC).
- **Initiation.** Terminal init applied once at each terminus (init_A/T 2.3/4.1, init_G/C 0.1/−2.8) — matches Biopython per-terminal-base logic and Wikipedia terminal A/T (9.6 kJ = 2.29 kcal ≈ 2.3) and G/C (0.4 kJ = 0.096 ≈ 0.1).
- **Tm equation.** `Tm = (1000·ΔH)/(ΔS + R·ln(C_T/4)) − 273.15`, R = 1.987 — identical to Biopython once C_T/4 is recognised as `k = (dnac1 − dnac2/2)` for equimolar non-self-complementary strands.
- **ΔG.** `ΔG = ΔH − 310.15·ΔS/1000` (37 °C) — standard Gibbs relation (INV-02).
- **Salt.** `ΔS += 0.368·(N−1)·ln[Na+]` (method 5) — matches Biopython/SantaLucia 1998.

### Edge-case semantics
- Length < 2 (empty / single base): NN model undefined → contract returns (0,0,0,0). Null guarded by `string.IsNullOrEmpty` (no throw). Documented as an API convention (not a thermodynamic value).
- Case-insensitive via `ToUpperInvariant`.
- F = 4 (non-self-complementary equimolar) is the documented in-scope case; self-complementary (F=1, sym term) and Mg²⁺ are explicitly out of scope.

### Independent cross-check (numbers)
A faithful independent Python reimplementation of the repo algorithm reproduced every expected value; real Biopython 1.85 confirms M1.

| Case | ΔH | ΔS | ΔG | Tm | Source of expectation |
|------|----|----|----|----|-----------------------|
| M1 CGTTCC…TTAC (Na 0.05, C_T 50nM) | — | — | — | **60.3** | Biopython 1.85 = 60.3209 (executed) |
| M2 GCGC (defaults) | −30.0 | −84.91 | −3.67 | −18.6 | DNA_NN3 + SantaLucia 1998 hand-derivation |
| M3 ATCG | −23.6 | −71.81 | (−1.33) | (−47.9) | two-end init derivation |
| M4 AATT | −18.4 | −59.91 | 0.18 | −75.0 | init_A/T×2 derivation |
| S2 GCGC Na 1.0 | −30.0 | −81.60 | — | −11.3 | salt-term derivation |
| C1 ATGC Wallace | — | — | — | 12 | 2(A+T)+4(G+C) |
| C2 20-mer Marmur-Doty | — | — | — | 51.78 | 64.9+41(10−16.4)/20 |

### Findings / divergences (Stage A)
None material. Note: the table commonly labelled "SantaLucia 1998 unified" (Wikipedia) and the "Allawi & SantaLucia 1997" `DNA_NN3` are numerically the same parameter set (kJ↔kcal); the doc/spec wording referencing both is accurate. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs:435-564` (NN table, constants, `CalculateThermodynamics`, `AddTerminalInitiation`) and `:569-590` (`CalculateMeltingTemperature`).
- `src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs` (Wallace 2/4, Marmur-Doty 64.9/41/16.4).

### Formula realised correctly?
Yes. Two-end initiation (`AddTerminalInitiation(upper[0])` + `(upper[^1])`), NN sum over overlapping dinucleotides via `TryGetValue`, salt `0.368·(N−1)·ln(Na)`, `ΔG = ΔH − 310.15·ΔS/1000`, `Tm = (ΔH·1000)/(ΔS + 1.987·ln(C_T/4)) − 273.15`. Rounding: ΔH/ΔS/ΔG to 2 dp, Tm to 1 dp.

### Cross-verification recomputed vs code
All 7 thermodynamic cases above recomputed independently match the code output (12 original tests green; values traced by hand-reimplementation). M1 matches real Biopython.

### Variant/delegate consistency
Delegate `CalculateMeltingTemperature` switches Wallace (len < 14) vs Marmur-Doty correctly; C1/C2 confirmed against the published formulas.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every MUST/SHOULD/COULD expectation traces to Biopython source / real Biopython execution / Wikipedia(SantaLucia 1998) / published Wallace & Marmur-Doty formulas, recomputed independently — not to code output.
- **No green-washing:** exact `.Within(1e-10)` equality assertions on known values; no weakened assertions, widened tolerances, skips, or value-fudging. The one `Is.GreaterThan` in S2 is a *monotonicity invariant* asserted **in addition to** the two exact Tm values (−18.6, −11.3), so it is not a weak substitute.
- **Coverage of logic & edges:** both public methods exercised; both formula paths of the delegate (Wallace, Marmur-Doty); both initiation branches (A/T and G/C termini, mixed and same); salt monotonicity; case-insensitivity; WC symmetry; Gibbs relation.
- **Gaps found and fixed this session:** (1) M6 (length-1) asserted only 2 of 4 fields → now asserts all four = 0. (2) No null-input coverage → added `CalculateThermodynamics_NullInput_ReturnsAllZero` and `CalculateMeltingTemperature_NullOrEmptyInput_ReturnsZero` (both guard via `string.IsNullOrEmpty`; verified no NRE).
- **Honest green:** full unfiltered suite **6521 passed, 0 failed** (1 pre-existing benchmark skipped, unrelated); `dotnet build` 0 errors. Changed test file builds warning-free.

### Findings / defects (Stage B)
- **Minor (BY-DESIGN, documented):** invalid bases (e.g. `N`, RNA `U`) are silently skipped in the NN loop (`TryGetValue`) while still counted in `N` for the salt term and possibly as a terminal-init base. Biopython instead raises on non-ACGT. This is outside the Stage-A documented scope (valid DNA, ACGT) and does not affect any valid-input result; left as-is and noted, no code change. No correctness defect for the validated domain.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: CLEAN.**
- No correctness defect. Two test-completeness gaps (M6 partial, missing null cases) were fixed in-session by adding/strengthening assertions to sourced edge behavior. One BY-DESIGN note recorded (silent non-ACGT skip).
