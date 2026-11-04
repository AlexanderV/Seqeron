# Test Specification: PROBE-LNATM-001

**Test Unit ID:** PROBE-LNATM-001
**Area:** MolTools
**Algorithm:** LNA-Adjusted NN Tm + MGB Probe Design
**Status:** ☑ Revalidated in review campaign 2026-09 (B07 F26–F28: base set + Owczarzy 2011 model + heteroduplex symmetry) — 2026-10-01
**Last Updated:** 2026-10-01

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | McTigue, Peterson & Kahn (2004) Biochemistry 43:5388–5405, DOI 10.1021/bi035976d — LNA-DNA NN ΔΔH°/ΔΔS° increments (32 entries). |
| 2 | MELTING 5 `McTigue2004lockedmn.xml` (Dumousseau et al. 2012, BMC Bioinformatics 13:101) — canonical machine-readable McTigue (2004) set; the differential oracle. |
| 3 | rmelting tutorial / MELTING `mct04` worked example `CCATT(L)GCTACC` → 63.61426 °C. |
| 4 | SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 / SantaLucia (1998) PNAS 95:1460 — base DNA unified NN set (PRIMER-NNTM-001). |
| 5 | Kutyavin et al. (2000) Nucleic Acids Res 28(2):655–661, DOI 10.1093/nar/28.2.655 — 3'-MGB design rules. |

## 1a. Campaign 2026-09 update (supersedes §3–§4 below where they differ)

- Base DNA set is the SantaLucia (1998) unified / Allawi & SantaLucia (1997) table (Biopython DNA_NN3, MELTING `all97`) —
  the set McTigue (2004) and Owczarzy (2011) parameterised on; the former SantaLucia & Hicks (2004) base caused the
  63.528 vs 63.614 °C gap that was mis-attributed to "base-set choice".
- New model `LnaNearestNeighborModel.Owczarzy2011` (default of the 2-argument overloads, as in MELTING): single-LNA (32),
  consecutive-LNA (16) and LNA-mismatch (96) parameters; `target` (3'→5') for mismatched duplexes; DNA internal mismatches
  via DNA_IMM1. McTigue model: isolated LNAs; runs use the Owczarzy tables (MELTING behaviour).
- Oracle: melting5.jar 5.2.0 run in full precision (`Main.getMeltingResults`): rmelting test values 63.61426 / 63.48299 /
  12.94323; Owczarzy (2011) triplet set (14 duplexes, 2 µM); McTigue (2004) duplexes (5 µM, both models); Na⁺ 50 mM;
  Mg²⁺ 3 mM; DNA mismatch. 4 200 random duplexes: ΔH°/ΔS° bit-identical, same not-computable set for internal LNAs (terminal LNAs:
  MELTING 5.2.0 computes them because its terminal-LNA `isApplicable` guard never fires — we return not computable, see below).
- Reduction (empty LNA set) = Biopython `Tm_NN(nn_table=DNA_NN3)` (R 1.987): 59.833634529845824, 75.52391452117843,
  46.231595611732246 (self-comp), 16.621392992113726 (self-comp, 50 mM, method 6).
- LNA-modified strands are never self-complementary duplexes; C_T/Na⁺/Mg²⁺/dNTP/R guards throw like `CalculateMeltingTemperatureNN`.
- MGB: qualitative rules only (quantitative ΔTm BLOCKED — vendor parameters unobtainable; see B07.md).

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
- **Terminal/out-of-range LNA ⇒ not computable** (`null` thermo / `NaN` Tm); McTigue / Owczarzy parameterise internal LNAs only.
  MELTING 5.2.0 intends the same (its `isApplicable` warns "not established for terminal locked nucleic acids" and returns
  false) but the guard compares the end pair with the literal `"L"`/`"-"` and never fires, so melting5 extrapolates the
  internal-LNA doublet parameter (`CLCATTGCTACC`, 0.1 mM, 1 M Na⁺: owc11 66.65650883512683 °C) — deliberately not reproduced
  (audit round 7, A7-5, B07 F66; test `NotComputable_Cases`).
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
