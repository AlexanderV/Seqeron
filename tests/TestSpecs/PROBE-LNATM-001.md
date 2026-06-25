# Test Specification: PROBE-LNATM-001

**Test Unit ID:** PROBE-LNATM-001
**Area:** MolTools
**Algorithm:** LNA-Adjusted NN Tm + MGB Probe Design
**Status:** ☑ Validated (Stage A ✅ / Stage B ✅ / CLEAN) — 2026-06-25
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | McTigue, Peterson & Kahn (2004) Biochemistry 43:5388–5405, DOI 10.1021/bi035976d — LNA-DNA NN ΔΔH°/ΔΔS° increments (32 entries). |
| 2 | MELTING 5 `McTigue2004lockedmn.xml` (Dumousseau et al. 2012, BMC Bioinformatics 13:101) — canonical machine-readable McTigue (2004) set; the differential oracle. |
| 3 | rmelting tutorial / MELTING `mct04` worked example `CCATT(L)GCTACC` → 63.61426 °C. |
| 4 | SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 / SantaLucia (1998) PNAS 95:1460 — base DNA unified NN set (PRIMER-NNTM-001). |
| 5 | Kutyavin et al. (2000) Nucleic Acids Res 28(2):655–661, DOI 10.1093/nar/28.2.655 — 3'-MGB design rules. |

## 2. Canonical Method(s)

`PrimerDesigner.CalculateMeltingTemperatureNNLna`,
`PrimerDesigner.CalculateNearestNeighborThermodynamicsLna`,
`ProbeDesigner.EvaluateMgbProbeDesign`

- **Source files:** `PrimerDesigner.cs` (LNA Tm + increment table), `ProbeDesigner.cs` (MGB rules)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs`

## 3. Contract / Invariants

- **Additivity:** LNA-adjusted ΔH°/ΔS° = base DNA NN stack + Σ McTigue increment per LNA-containing step.
- **Reduction:** empty LNA set ⇒ result equals PRIMER-NNTM-001 exactly (`Within 1e-9`).
- **Each internal LNA contributes two increments** (3'-locked for its left step, 5'-locked for its right step).
- **Terminal/out-of-range LNA ⇒ not computable** (`null` thermo / `NaN` Tm); McTigue parameterises internal NN only.
- **Non-ACGT ⇒ not computable.** Order/duplicates of positions tolerated (set semantics).
- **Determinism.** MGB rules: boolean length-window + 3'-attachment guidance; quantitative ΔTm intentionally not computed (empirical, no closed form).

## 4. Cross-check / Differential Oracle

- **Reference:** MELTING 5 `McTigue2004lockedmn.xml` (binary/R not installable here → data file used directly) + hand-derivation.
- **Exact numbers:** all 32 increments verbatim-match the MELTING XML (0 mismatches). Base DNA NN `CCATTGCTACC` = ΔH° −80.8 / ΔS° −221.7; all-DNA Tm 59.692264 °C; LNA(L₄) ΔH° −80.014 / ΔS° −216.6, Tm 63.527594 °C vs MELTING `mct04` 63.61426 °C (Δ +0.087 °C < 0.1, base-DNA-set choice).
- **Comparison gate:** LNA Tm within 0.1 °C of MELTING `mct04`; reduction `Within 1e-9`.

## 5. Validation Checklist (restored to ☑)

- [x] Stage A: every source retrieved this session; all 32 increments cross-checked verbatim vs MELTING XML; XML key→(step,locked) decoding verified; MGB rules confirmed against Kutyavin (2000).
- [x] Stage B: implementation realises the additive McTigue+SantaLucia model; reduces to PRIMER-NNTM-001 on empty LNA set; tests assert sourced (non-echoed) values and cover all edge cases.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics.Tests 18741 passed).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the `docs/checklists/*.md`.
