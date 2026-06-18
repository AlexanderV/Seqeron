# Validation Report: SEQ-TM-001 — Melting Temperature (Wallace / Marmur-Doty / Nearest-Neighbor Tm)

- **Validated:** 2026-06-16   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateThermodynamics(dnaSequence, naConcentration, primerConcentration)` (NN ΔH°/ΔS°/ΔG°/Tm) and `SequenceStatistics.CalculateMeltingTemperature(dnaSequence, useWallaceRule)` (Wallace / Marmur-Doty via `ThermoConstants`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN — **DUPLICATE-OF SEQ-THERMO-001** (consolidated; same two methods, same canonical fixture, no new production code). One test-coverage gap fixed this session.

## Duplicate determination

SEQ-TM-001's TestSpec §7, Evidence, algorithm doc §5.4, and `ALGORITHMS_CHECKLIST_V2.md` all state SEQ-TM-001 is a **duplicate Registry entry** for the exact two melting-temperature methods already delivered under **SEQ-THERMO-001** (validated `✅ / ✅ / ✅ CLEAN`, ledger #27). Confirmed by inspection: both units list the identical `CalculateThermodynamics` + `CalculateMeltingTemperature` on `SequenceStatistics`, and both point at the single canonical fixture `SequenceStatistics_CalculateThermodynamics_Tests.cs`. Per the duplicate-elimination rule, SEQ-TM-001 reuses the existing implementation and fixture rather than re-implementing. The sibling **PRIMER-TM-001** (Phase 1) independently exact-locks the shared `ThermoConstants` Wallace/Marmur-Doty formula methods. Despite the duplicate status, a **full independent validation** was performed this session against external sources (below), not a rubber-stamp.

## Stage A — Description

### Sources opened & what they confirm
- **Biopython `Bio/SeqUtils/MeltingTemp.py`** (WebFetched raw master, this session) — verbatim:
  - Wallace: `Tm = 4·(G+C) + 2·(A+T)`; docstring `Tm_Wallace('ACGTTGCAATGCCGTA')` = **48.0**.
  - GC valueset 1: `Tm = 69.3 + 0.41(%GC) − 650/N` (Marmur & Doty 1962; Chester & Marshak 1993).
  - `Tm_NN` return: `(1000*ΔH)/(ΔS + R·ln(k)) − 273.15`, **R = 1.987**, `k = (dnac1 − dnac2/2)·1e-9`, defaults dnac1=dnac2=25 nM (⇒ k = C_T/4 with C_T = 50 nM), Na = 50 mM, nn_table = DNA_NN3.
  - DNA_NN3: all 16 dinucleotides — AA −7.9/−22.2, AT −7.2/−20.4, TA −7.2/−21.3, CA −8.5/−22.7, GT −8.4/−22.4, CT −7.8/−21.0, GA −8.2/−22.2, CG −10.6/−27.2, GC −9.8/−24.4, GG −8.0/−19.9; init_A/T (2.3, 4.1), init_G/C (0.1, −2.8).
  - Salt method 5: `0.368·(N−1)·ln[Na+]` (SantaLucia 1998, PNAS 95:1460–1465).
- **UGENE / Primer3 GC method** (WebSearch, this session) — confirms the classic GC form the repo uses: `Tm = 64.9 + 41·(GC − 16.4)/N` for length ≥ 15 (50 nM primer, 50 mM Na⁺).
- **SantaLucia (1998) PNAS** & **Allawi & SantaLucia (1997) Biochemistry** — primary refs for the unified NN parameter table, the `Tm = ΔH°/(ΔS° + R·ln(C_T/x)) − 273.15` equation (x = 4 non-self-complementary), and the salt correction.

### Formula check
All three formulas in the implementation match the cited sources exactly: Wallace `2(A+T)+4(G+C)`, Marmur-Doty `64.9 + 41(GC−16.4)/N`, NN `(1000ΔH)/(ΔS + R·ln(C_T/4)) − 273.15` with R = 1.987 and the SantaLucia salt correction. The repo's NN parameter table and initiation terms are identical to Biopython DNA_NN3.

### Edge-case semantics
NN undefined for length < 2 → all-zero tuple (sourced: NN sums over dinucleotides). Wallace/GC delegate returns 0 for null/empty. Wallace is a rule of thumb only below 14 nt; the impl auto-switches to Marmur-Doty at length ≥ 14. All defined and sourced.

### Independent cross-check (numbers, recomputed from sources — not from the repo code)
A Python reimplementation built **only from the externally-sourced constants** reproduced every test value:

| Input | Method | Params | Result (sourced) |
|-------|--------|--------|------------------|
| `ATGC` | Wallace | — | 2·2+4·2 = **12.0** |
| `GCGCGCGCGCATATATATAT` (N20, GC10) | Marmur-Doty | — | 64.9+41(10−16.4)/20 = **51.78** |
| `ACGTTGCAATGCCGTA` (N16, GC8), useWallace=true | auto-switch → Marmur-Doty | — | 64.9+41(8−16.4)/16 = **43.375** |
| `GCGC` | NN | Na 0.05, C_T 250 nM | ΔH −30.0, ΔS −84.91, ΔG −3.67, Tm **−18.6** |
| `ATCG` | NN | defaults | ΔH −23.6, ΔS −71.81 |
| `AATT` | NN | defaults | ΔH −18.4, ΔS −59.91, ΔG 0.18, Tm −75.0 |
| `CGTTCCAAAGATGTGGGCATGAGCTTAC` | NN | Na 0.05, C_T 50 nM | Tm **60.3** (= Biopython published 60.32) |
| `GCGC` | NN salt | 0.05 vs 1.0 M | −18.6 vs −11.3 (monotonic ↑) |

### Findings / divergences (Stage A)
- Repo Marmur-Doty uses the **classic** `64.9 + 41(GC−16.4)/N` form, not Biopython valueset-1's `69.3 + 0.41(%GC) − 650/N`. Both are legitimate published GC forms; the classic form is independently sourced (UGENE/Primer3). Documented, not a defect.
- Repo NN default C_T = 250 nM vs Biopython's 50 nM — explicit, documented parameter; identical formula; passing C_T = 50 nM reproduces 60.32 exactly. Documented assumption, not a defect.

Stage A verdict: **PASS**.

## Stage B — Implementation

### Code path reviewed
- `SequenceStatistics.cs:502–564` — `CalculateThermodynamics`: two-end initiation, dinucleotide sum, salt correction, ΔG°37 = ΔH − 310.15·ΔS/1000, Tm equation with R = 1.987 and C_T/4. Rounding ΔH/ΔS/ΔG to 2 dp, Tm to 1 dp.
- `SequenceStatistics.cs:435–476` — NN parameter table + init/salt/R/F constants (match Biopython DNA_NN3 cell-for-cell).
- `SequenceStatistics.cs:569–590` — `CalculateMeltingTemperature`: Wallace (length < 14) / Marmur-Doty auto-switch, delegating constants to `ThermoConstants`.
- `ThermoConstants.cs:87–99` — `CalculateWallaceTm` (`2·AT + 4·GC`), `CalculateMarmurDotyTm` (`64.9 + 41(GC−16.4)/N`, length-0 guard).

### Formula realised correctly?
Yes. The code computes the validated formulas verbatim for the right inputs; constants are source-cited inline and all match the external references retrieved this session.

### Cross-verification table recomputed vs code
The actual implementation output (run via the test suite) equals the independently-sourced values in the Stage-A table for every case (M1–M7, S1–S3, C1, C2, and the new auto-switch case).

### Variant/delegate consistency
`CalculateMeltingTemperature` delegates to the same `ThermoConstants` formulas used by `PrimerDesigner.CalculateMeltingTemperature` (PRIMER-TM-001), which exact-locks them in `PrimerDesigner_MeltingTemperature_Tests.cs`. Consistent.

### Test quality audit (HARD gate)
- **Sourced, not echoed:** M1 60.3, M2 GCGC tuple, M3/M4 ATCG/AATT, C1 Wallace 12, C2 Marmur-Doty 51.78 all asserted as exact `Within` values that match the externally-sourced computation; a deliberately-wrong implementation would fail them.
- **All branches/edges:** NN both-end init (GC/AT termini), salt, ΔG Gibbs relation, empty/length-1/null, case-insensitivity, WC symmetry, salt monotonicity all covered.
- **Gap found & fixed (test-only, 0 code change):** the `CalculateMeltingTemperature` **auto-switch** branch (`useWallaceRule=true` but length ≥ 14 → Marmur-Doty) was unexercised. Added `CalculateMeltingTemperature_WallaceRequestedButTooLong_UsesMarmurDoty` asserting `ACGTTGCAATGCCGTA` (16 nt, GC8, useWallaceRule=true) → **43.375** (the sourced Marmur-Doty value, explicitly NOT the Wallace 48.0, because length 16 ≥ 14 disables Wallace). Fixture 16 → 17.
- **Honest green:** no assertion weakened, no skip, no tolerance widened, no expected value adjusted to match output.

### Findings / defects (Stage B)
No code defect. One test-coverage gap (auto-switch branch), fixed this session.

Stage B verdict: **PASS**.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: ✅ CLEAN (DUPLICATE-OF SEQ-THERMO-001).**
- **Test-quality gate: PASS** — full unfiltered suite **6613 passed, 0 failed, 0 skipped**; `dotnet build` 0 errors; changed test file warning-free. (Suite was 6612; +1 new evidence-locked test.)
- Finding logged in `FINDINGS_REGISTER.md`; ledger row #107 added.
- No production code changed; the algorithm is fully functional.
