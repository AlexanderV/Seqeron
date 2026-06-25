# Validation Report: PROBE-LNATM-001 — LNA-Adjusted NN Tm + MGB Probe Design

- **Validated:** 2026-06-25   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.CalculateMeltingTemperatureNNLna`,
  `PrimerDesigner.CalculateNearestNeighborThermodynamicsLna`, `ProbeDesigner.EvaluateMgbProbeDesign`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
- LNA-adjusted thermodynamics + Tm: `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs`
  (`CalculateNearestNeighborThermodynamicsLna` lines 1054–1093; `CalculateMeltingTemperatureNNLna`
  lines 1113–1154; McTigue increment table lines 998–1035).
- MGB qualitative design rules: `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`
  (`EvaluateMgbProbeDesign` lines 459–480).
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs`.

## Authoritative sources opened (this session)
1. **McTigue, Peterson & Kahn (2004)** Biochemistry 43:5388–5405, DOI 10.1021/bi035976d —
   sequence-dependent ΔΔH°/ΔΔS° LNA-DNA nearest-neighbour increments (32 entries).
2. **MELTING 5** data file `McTigue2004lockedmn.xml` (Dumousseau et al. 2012, BMC Bioinformatics
   13:101; EMBL-EBI) — the canonical machine-readable transcription of the McTigue (2004) set; the
   oracle used here. Retrieved verbatim from the rmelting mirror
   (`aravind-j/rmelting/inst/extdata/Data/McTigue2004lockedmn.xml`).
3. **rmelting tutorial / MELTING `mct04`** worked example `CCATT(L)GCTACC` → Tm 63.61426 °C
   (Bioconductor rmelting vignette).
4. **SantaLucia & Hicks (2004)** Annu Rev Biophys 33:415 / SantaLucia (1998) PNAS 95:1460 — base DNA
   unified NN set (the PRIMER-NNTM-001 model the LNA increments add onto).
5. **Kutyavin et al. (2000)** Nucleic Acids Res 28(2):655–661, DOI 10.1093/nar/28.2.655 — 3'-MGB
   design rules.

## Model (validated)
LNA-adjusted duplex thermodynamics = base DNA NN stack (SantaLucia unified params + initiation +
terminal-A·T + symmetry, all computed on the underlying DNA sequence) **plus** the additive McTigue
(2004) ΔΔH°/ΔΔS° increment for each NN step that contains an LNA base. A single internal LNA base
participates in **both** of its flanking dinucleotide steps (one with the LNA as the 3' base of the
step, one with it as the 5' base), so two increments are added per internal LNA. Tm then follows the
standard bimolecular equation `Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15`, R = 1.9872,
x = 4 (non-self-comp) / 1 (self-comp). Terminal LNA positions (index 0 or length−1) are **not**
parameterised by McTigue (2004) and are correctly rejected (per MELTING `isApplicable`). Empty
LNA set reduces **exactly** to the standard NN Tm (PRIMER-NNTM-001). MGB rules are **qualitative**
only (3'-attachment + 12–20mer length window); the quantitative MGB ΔTm is empirical with no closed
form and is honestly left as a residual.

## Stage A — Description

### Source quality & formula correctness
- All 32 repo increments were cross-checked **byte-for-byte** against the MELTING
  `McTigue2004lockedmn.xml` data file: parsed the XML, mapped each `<modified sequence="…">` key to
  the repo's `(step, locked-position)` entry, and compared. **0 mismatches across all 32 keys**
  (cal/mol ÷ 1000 = kcal/mol; ΔS in cal/(mol·K) verbatim). Every XML key is covered and no repo key
  is spurious.
- The XML key naming was decoded independently: `TTL/AA` → step `TT` with the **3' base** locked;
  `TLG/AC` → step `TG` with the **5' base** locked. This exactly matches the code's
  `LnaStepPosition.ThreePrime` / `FivePrime` lookup logic.
- MELTING-notation `CCATTLGCTACC` means the base **immediately before** the `L` is locked (MELTING
  uses the `Al/Cl/Gl/Tl` "L-follows-base" convention), i.e. the **second T at 0-based index 4** —
  confirming the test's `WorkedLnaIndex = 4`.

### Edge-case semantics
- Empty LNA set → standard NN Tm (sourced: additive model, increment = 0). ✓
- Terminal LNA → not computable (sourced: McTigue parameterises internal NN only). ✓
- Non-ACGT base → not computable (base DNA NN lookup fails). ✓
- Out-of-range / duplicate / unsorted positions → out-of-range null, duplicate counts once,
  order-independent (set semantics). ✓
- MGB: 12–20mer in-range true, out-of-range flagged with the cited window; 3' attachment always
  reported. ✓

### Independent cross-check (exact numbers)
Oracle: **MELTING 5 `McTigue2004lockedmn.xml`** + hand-derivation (MELTING binary/R not installable
in this environment; documented).

| Quantity | Hand-derived (this session) | Reference | Δ |
|---|---|---|---|
| Base DNA NN `CCATTGCTACC` | ΔH° = −80.8000, ΔS° = −221.7000 | SantaLucia unified set (= repo `CalculateNearestNeighborThermodynamics`) | exact |
| All-DNA Tm `CCATTGCTACC` (C=1e-4, Na=1, no salt) | 59.692264 °C | PRIMER-NNTM-001 reduction | exact |
| LNA `CCATT(L₄)GCTACC` ΔH°/ΔS° | −80.014 / −216.6 | base + TTL/AA(+2.326,+8.1) + TLG/AC(−1.540,−3.0) (MELTING XML) | exact |
| LNA Tm `CCATT(L₄)GCTACC` | 63.527594 °C | MELTING `mct04` reports **63.61426 °C** | **+0.087 °C** |
| M6 `GGGCC` LNA@1 ΔΔ | ΔΔH = −3.787, ΔΔS = −7.6 | GGL/CC(−0.943,−0.9) + GLG/CC(−2.844,−6.7) | exact |
| C1 add LNA@6 to `CCATTGCTACC` ΔΔ | ΔΔH = −0.217, ΔΔS = +3.1 | GCL/CG(−0.925,−1.1) + CLT/GA(+0.708,+4.2) | exact |

The residual **+0.087 °C** vs MELTING `mct04` is < 0.1 °C and is fully explained: MELTING's `mct04`
runs the McTigue increments on top of MELTING's own default DNA NN set, which differs slightly from
the SantaLucia (1998) unified set used here. The **LNA increment contribution is identical**; only
the base-DNA stack differs. This is the documented, sourced explanation — not a tuning fudge.

### MGB (Kutyavin 2000) cross-check
- 3'-end attachment: confirmed ("3′-MGB-ODNs are easier to prepare … MGB-modified solid supports and
  automated DNA synthesis can be used").
- 12–20mer window: confirmed ("For MGB probes this length variation is narrowed to a range of
  12–20mers"); the 12mer-MGB ≈ 27mer-unmodified Tm equivalence (Tm 66 °C vs 65 °C) justifies leaving
  the quantitative ΔTm as a residual.

**Stage A verdict: ✅ PASS.** Description and constants are correct and primary-sourced; no divergence.

## Stage B — Implementation

### Code realises the validated model
- `CalculateNearestNeighborThermodynamicsLna` computes the base DNA NN via the unchanged
  `CalculateNearestNeighborThermodynamics`, then for each step (i, i+1) adds the 5'-locked increment
  when `i` is locked and the 3'-locked increment when `i+1` is locked — i.e. one internal LNA picks
  up **both** flanking increments. Matches the McTigue/MELTING `enthalpy += lockedAcidValue` model.
- Terminal/out-of-range positions return `null` before any increment is added (lines 1068–1071).
- `CalculateMeltingTemperatureNNLna` reuses the identical Tm equation, x-factor and salt-correction
  switch as `CalculateMeltingTemperatureNN`; with an empty LNA set the thermodynamics are byte-equal
  to the perfect-match path, so it reduces **exactly** to PRIMER-NNTM-001 (test M4 asserts `.Within(1e-9)`).
- `EvaluateMgbProbeDesign` checks only the length window and emits the 3'-attachment guidance string;
  it deliberately computes **no** quantitative ΔTm.

### Cross-verification recomputed vs code
The fixture (13 tests) was run; every value matches the Stage-A hand/oracle figures above to the
asserted tolerances (`1e-9` for ΔH/ΔS and reduction; `1e-4` for Tm).

### Test-quality audit (HARD gate)
- Expected values are **sourced** (McTigue/MELTING XML increments and hand-derived Tm), not code
  echoes: M2 asserts the literal `63.527594` AND independently `|tm − 63.61426| < 0.1` against MELTING;
  M6/C1 assert literal McTigue increment sums.
- No green-washing: no `Assert.Pass`, no `Ignore`/`Skip`, no widened tolerances, no tautologies.
- Coverage: all three public methods + every Stage-A edge case (zero LNA → equals standard NN [M4];
  1 LNA [M1/M2]; multiple LNA additive [C1]; negative increment sign [M6]; terminal → null/NaN [M5];
  non-ACGT → NaN [S3]; null/empty/short → null + null-mask throws [S1]; out-of-range/dup/unsorted [S2];
  all 32 internal contexts parameterised [S4]; MGB pass/fail with reasons + 3' end [M7]; MGB null
  throws). The C1/M6 assertions write the increment sum expression literally (e.g. `-0.925 + 0.708`),
  keeping the sourced provenance visible.

### Full suite
`dotnet test Seqeron.sln -c Debug` (unfiltered): all 6 NUnit assemblies **Failed: 0**
(`Seqeron.Genomics.Tests` = 18741 passed, 0 failed). No code/tests modified, 0 warnings.

**Stage B verdict: ✅ PASS.** Implementation faithfully realises the validated model; tests are real,
sourced, and complete.

## Verdict & follow-ups
- **State: ✅ CLEAN.** No defect found in description, implementation, or tests. The LNA increment
  table is verbatim-identical to the MELTING McTigue (2004) data file; the all-DNA reduction equals
  PRIMER-NNTM-001 exactly; the MGB qualitative rules match Kutyavin (2000).
- **Documented residual (by design, not a defect):** the quantitative MGB ΔTm is empirical with no
  published closed form (Kutyavin 2000) and is intentionally not computed — already recorded in
  `docs/Validation/LIMITATIONS.md` and reflected by `EvaluateMgbProbeDesign` returning qualitative
  guidance only. The ~0.087 °C LNA-Tm offset vs MELTING `mct04` is the base-DNA-NN-set choice, < 0.1 °C,
  and fully explained.
- Oracle note: MELTING binary / R / rmelting are not installable in this environment; the canonical
  MELTING `McTigue2004lockedmn.xml` data file was retrieved and used directly as the differential
  oracle, supplemented by hand-derivation against McTigue (2004) and the published `mct04` Tm.
