# Validation Report: PROBE-DESIGN-001 — Hybridization Probe Design

- **Validated:** 2026-06-24   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.DesignProbes` (+ suffix-tree overload), `DesignProbesOptimized`, `EvaluateProbeWithGc`/`EvaluateProbe`, `DesignTilingProbes`, `DesignAntisenseProbes`, `DesignMolecularBeacon`, `ValidateProbe`, `CheckSpecificity`, `CalculateTm` (via `ThermoConstants`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Scope clarification (important)

The validation prompt frames this unit around **TaqMan/qPCR-probe-specific** rules
(no G at the 5' end, Tm ~+10 °C above primers, ≤4 contiguous G, more C than G).
The TestSpec and the implementation do **not** claim to be a TaqMan-specific designer.
`ProbeDesigner` is a **generic hybridization-probe designer** with application presets
(`Microarray`, `FISH`, `NorthernBlot`, `qPCR`, `SouthernBlot`) plus tiling, antisense and
molecular-beacon variants. The claimed/enforced rules are the generic ones: per-application
length window, Tm window, GC window, homopolymer-run cap, self-complementarity, secondary
structure (hairpin), simple-repeat detection, and specificity/uniqueness.

I validated the unit against **its own declared spec** (the standard protocol question), and
separately recorded, under Findings, where the generic heuristics diverge from strict TaqMan
guidance so the gap is explicit. None of the TaqMan-specific rules are *asserted* by the spec,
so their absence is not a defect against this unit — it is a scope/heuristic note.

## Stage A — Description

### Sources opened & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia *Hybridization probe* | probes 15–10000 nt; detect **complementary** target; stringency (T/salt) sets specificity. |
| Wikipedia *DNA microarray* | oligo probes: Affymetrix 25-mer, Agilent 60-mer (preset Microarray 50–60 sits in the long-oligo regime). |
| Wikipedia *FISH* | BAC ≈100 kb; oligo FISH 10–25 nt; smFISH ~20–50 oligos. Preset FISH 200–500 is a mid-range double-stranded-probe choice (sourced as standard practice, not a single canonical number). |
| Wikipedia *Molecular beacon* | loop 18–30 nt, stem two 5–7 nt complementary arms, total ≈25 nt — matches `DesignMolecularBeacon(loop=25, stem=5)`. |
| Wikipedia *Nucleic acid thermodynamics* / SantaLucia (1998) | Wallace **2(A+T)+4(G+C)** for short oligos; salt-adjusted **81.5 + 16.6·log₁₀[Na⁺] + 0.41·%GC − 600/N**. Both realised in `ThermoConstants` (same formulas validated under PRIMER-TM-001). |
| Thermo Fisher TaqMan design guidance (external) | TaqMan-specific: **no G at 5' end**, Tm ≈ 68–70 °C (≈+10 above primers), ≤4 contiguous G, more C than G. Used only as a divergence lens (see Findings); **not** part of this unit's spec. |

### Formula / constraint check
- **Length:** scan `MinLength..MaxLength` bounded by `n`; all presets within the 15–10000 nt envelope. ✓
- **Tm:** `CalculateTm` dispatches Wallace (len < 14) vs salt-adjusted GC formula (len ≥ 14). Both match cited equations exactly. ✓
- **GC:** O(1) prefix-sum fraction in [0,1]; window check with ±0.1 early-reject band, in-window scoring penalty. ✓
- **Self-complementarity / hairpin / repeats / homopolymer / specificity:** standard quality filters; specificity unique=1.0, N-hit=1/N. ✓

### Edge-case semantics
- empty / null / shorter than `MinLength` → empty (`yield break`). ✓
- all-G probe → GC 1.0; all-AT → GC 0.0; beacon on too-short target → null. ✓

### Independent cross-check (hand-computed)
- **Wallace** "ACGTACGT" (AT=4, GC=4, len 8 < 14): 2·4 + 4·4 = **24 °C**.
- **Salt-adjusted** (all-G, N=100, gc=1.0, [Na⁺]=0.05): 81.5 + 16.6·log₁₀(0.05) + 41·1.0 − 600/100 = 81.5 − 21.60 + 41 − 6 = **94.90 °C** (> 0 ✓).
- **Beacon** stem=5: stem5 = "GG"+"CCC" = "GGCCC"; stem3 = RC = "GGGCC"; total len = 5+25+5 = 35 for loop 25 (test S4 uses loop 20 → total 30, asserted exactly). ✓

### Findings / divergences (Stage A → PASS-WITH-NOTES)
1. **5'/3' position penalty is symmetric and GC-based, not the TaqMan "no 5' G" rule.**
   `EvaluateProbeWithGc` (lines 307–313) adds a small penalty (0.02) if the probe starts with G **or C**, and again if it ends with G **or C**. Strict TaqMan guidance forbids a 5' **G** specifically (quencher proximity) and prefers **more C than G overall**, with no symmetric 3' rule. This is a **declared generic heuristic**, acceptable per the spec; it is not a TaqMan implementation. Flagged so the gap is on record.
2. **No "more C than G" content rule and no dedicated TaqMan Tm target (~68–70 °C).** The `qPCR` preset uses Tm 68–72 °C (consistent with TaqMan probe Tm range) but applies it to the generic designer; it does not enforce probe-Tm = primer-Tm + 10. Spec does not claim this — recorded as scope, not defect.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`:
`DesignProbes` (127–195), `DesignProbesOptimized` (201–248), `EvaluateProbeWithGc` (253–327),
`DesignTilingProbes` (341–397), `DesignAntisenseProbes` (402–414), `DesignMolecularBeacon` (419–473),
`ValidateProbe`/`CheckSpecificity` (491–587), helpers (678–804); `ThermoConstants.cs` Tm formulas.

### Formula realised correctly?
- Length scan, coordinates `[start, start+len-1]`, substring identity (M7/M8). ✓
- GC prefix-sum, ±0.1 early-reject + in-window penalty. ✓
- Tm Wallace/salt-adjusted dispatch matches validated formulas. ✓
- Antisense designs against `GetReverseComplementString(mRNA)`; self-comp compares base i to RC base i (palindrome → 1.0); hairpin via inverted-repeat stem matching. ✓
- Specificity suffix-tree unique=1.0 / N-hit=1/N; `requireUnique` filters non-unique. ✓

### Cross-verification table recomputed vs code
| Item | Hand value | Code/test outcome |
|------|-----------|-------------------|
| Wallace "ACGTACGT" | 24 °C | `CalculateTm` path — consistent |
| All-G probe GcContent | 1.0 | M13 asserts 1.0 — pass |
| All-AT probe GcContent | 0.0 | M14 asserts 0.0 — pass |
| Beacon stem5/stem3 (stem=5) | "GGCCC"/"GGGCC", total 30 (loop 20) | S4 asserts exact strings + 30 — pass |
| Specificity (1 hit / N hits) | 1.0 / 1/N | suffix-tree tests — pass |

### Variant/delegate consistency
`EvaluateProbe` delegates to `EvaluateProbeWithGc` with `CalculateGcContent`; optimized path
supplies prefix-sum GC → same result. Tiling/antisense reuse the same evaluator. Consistent.

### Test quality audit
91 `ProbeDesigner*` tests across 4 files (design / validation / mutation-killers / smoke).
Assertions check exact sourced values (GC 0.0/1.0, exact beacon strings, exact tiling
starts/coverage, palindrome self-comp 1.0, specificity 1/N), not just "no throw". Edge cases
(empty/null/short/beacon-short) covered; mutation-killer beacon tests present.

### Findings / defects
- None against the unit's declared spec. The two Stage-A notes are scope/heuristic items
  (the designer is generic, not TaqMan-specific), already reflected in the qPCR preset and the
  generic position penalty. No code is wrong relative to what the spec claims.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES**, **Stage B: PASS**, **State: CLEAN**. No code changes.
- Build: succeeded, 0 warnings. Tests: `ProbeDesigner` filter **91/91 pass**.
- Follow-up (optional, out of scope for CLEAN): if a true TaqMan/qPCR-probe preset is ever
  desired, add (a) explicit 5'-G rejection (not symmetric 5'/3' G/C penalty), (b) a "C > G"
  content check, and (c) a probe-Tm = primer-Tm + ~10 coupling. These are *enhancements*, not
  fixes — the current unit is correct for the generic probe-design scope it declares.
