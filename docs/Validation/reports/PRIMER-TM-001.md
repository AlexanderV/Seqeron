# Validation Report: PRIMER-TM-001 — Primer Melting Temperature (Tm)

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `PrimerDesigner.CalculateMeltingTemperature(string)`, `PrimerDesigner.CalculateMeltingTemperatureWithSalt(string, double)`; helpers in `ThermoConstants` (`CalculateWallaceTm`, `CalculateMarmurDotyTm`, `CalculateSaltCorrection`).
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect; methods are honestly labeled simplifications, internally correct, all tests pass)

---

## Methods present in the library

The library implements **three** Tm formulas; it does **not** implement a nearest-neighbor (SantaLucia 1998) Tm path for the canonical primer Tm. SantaLucia appears only as a cited comparison source in the spec, not as a code path here. The implemented methods are:

1. **Wallace rule** (basic) — short oligos (< 14 valid bases)
2. **Marmur-Doty / basic GC% form** — long oligos (≥ 14 valid bases)
3. **Salt correction** (Schildkraut-Lifson 16.6·log10 coefficient)

`ThermoConstants` also defines an unused `CalculateSaltAdjustedTm` (81.5 + 16.6·log10[Na] + 41·%GC − 600/N form); it is not wired into the canonical Tm path and is out of scope for the two canonical methods, but its constants were checked and are standard.

---

## Stage A — Description

### Sources opened
- Wikipedia, *Nucleic acid thermodynamics* (NN two-state model + NN ΔH°/ΔS° table; no GC%/Wallace constants — confirms NN is the rigorous model but does not contradict the simpler methods used).
- OpenWetWare, *Primer Tm estimation methods* — confirms both simple formulas verbatim.
- novoprolabs / sciencecodons / Sigma-Aldrich references — confirm Wallace and basic GC% formulas and constants.
- calculatorsconversion / oligocalc references — confirm the 16.6·log10([Na+]) salt term (Schildkraut-Lifson).

### Formula check (each method)

| Method | Code formula | Authoritative source formula | Match |
|--------|--------------|------------------------------|-------|
| Wallace | `2·AT + 4·GC`, len < 14 | Tm = 2(A+T) + 4(G+C), short oligos < 14 nt (Thein & Wallace 1986) | ✅ exact |
| Marmur-Doty (basic GC%) | `64.9 + 41·(GC−16.4)/N`, len ≥ 14 | Tm = 64.9 + 41·(yG+zC−16.4)/(wA+xT+yG+zC) (basic/"modified Marmur-Doty", ≥14 nt) | ✅ exact constants 64.9 / 41 / 16.4 |
| Salt correction | `16.6·log10(Na_mM/1000)` | ΔTm = 16.6·log10([Na+]) coefficient (Schildkraut-Lifson; cited via Owczarzy 2004) | ✅ coefficient exact; see note below |

### Naming note (Stage-A divergence, documented)
There is genuine naming inconsistency in the public literature: some sources (e.g. OpenWetWare) label `2(A+T)+4(G+C)` as "Marmur-Doty" and `64.9+41(GC−16.4)/N` as the "Wallace formula" — i.e. the reverse of this spec's labels. The **formulas, constants, and length applicability are correct regardless of which label is attached**; only the names are swapped between sources. Not a defect.

### Salt-correction reference-state note (Stage-A note)
The Schildkraut-Lifson relation is a **difference** between two salt states: ΔTm = 16.6·log10([Na+]₂/[Na+]₁). The code computes an **absolute** term `16.6·log10([Na]_M)` (i.e. referenced to a 1 M Na+ standard, since Na_mM/1000 = [Na]_M). Consequently at the standard 50 mM PCR condition the term is ≈ −21.6 °C, so `WithSalt` returns a Tm *below* the base Tm. As an absolute "Tm at 50 mM" this is physically counter-intuitive (real PCR Tm at 50 mM is close to the base value, not 21 °C lower). However, the spec explicitly defines this exact term and its expected outputs (M12 ≈ 30.2, M13 ≈ 18.58, M14 ≈ 40.18), and it is a self-consistent monotone salt response (more salt → higher Tm). It is an honestly-documented simplification, not an internal error. Flagged here so the limitation is on record.

### Edge-case semantics
Empty/null → 0; only ACGT counted (N, IUPAC, RNA U ignored); all-invalid → 0; case-insensitive; threshold strictly `< 14` valid bases for Wallace. All defined and sourced in the spec (Section 3). Consistent.

### Independent cross-check (hand computations)
- "ATATATAT": 2·8+4·0 = **16.0** (spec M3) ✓
- "GCGCGCGC": 4·8 = **32.0** (M4) ✓
- "ACGTACGT": 2·4+4·4 = **24.0** (M5) ✓
- 20 bp 50% GC: 64.9+41·(10−16.4)/20 = **51.78** (M8) ✓
- 20 bp 100% GC: 64.9+41·(3.6)/20 = **72.28** (M10) ✓
- 16 bp all-A: 64.9+41·(−16.4)/16 = **22.875** (M19) ✓
- Salt 50 mM: 16.6·log10(0.05) = **−21.597** (M12) ✓

Stage A verdict: **PASS-WITH-NOTES** — every implemented formula and constant matches an authoritative source; notes are (a) literature label swap and (b) salt-correction reference state, both documented, neither a math error.

---

## Stage B — Implementation

### Code path reviewed
- `src/.../Seqeron.Genomics.MolTools/PrimerDesigner.cs:197-235` (`CalculateMeltingTemperature`, `CalculateMeltingTemperatureWithSalt`)
- `src/.../Seqeron.Genomics.Infrastructure/ThermoConstants.cs:15-122` (constants + `CalculateWallaceTm`, `CalculateMarmurDotyTm`, `CalculateSaltCorrection`)

### Formula realised correctly?
- Constants in code: `WallaceMaxLength=14`, `WallaceAtContribution=2`, `WallaceGcContribution=4`, `MarmurDotyBase=64.9`, `MarmurDotyGcCoefficient=41.0`, `MarmurDotyGcOffset=16.4`, `SaltCoefficient=16.6`. All match the validated sources exactly.
- Wallace branch: `validLength < 14` → `2·AT + 4·GC`. ✅
- Marmur-Doty branch: `64.9 + 41·(gc−16.4)/validLength`, clamped `Math.Max(0, …)`. ✅ (clamp explains INV-1; note it can mask negative raw values for extreme low-GC very-short-≥14 inputs, but those are not realistic primers and the spec's INV-1 sources the clamp).
- Salt: `baseTm + 16.6·log10(Na_mM/1000)`, rounded to 1 dp. ✅ matches spec term.
- `validLength = at + gc` (only ACGT) drives BOTH threshold and formula → non-ACGT/U correctly ignored. ✅

### Cross-verification table recomputed vs code (tests executed)
All 41 melting tests pass; hand values above equal code outputs. Key matches: M3=16, M4=32, M5=24, M8≈51.78, M10≈72.28, M19=22.875, salt-50mM correction≈−21.6.

### Variant/delegate consistency
`CalculateMeltingTemperatureWithSalt` = `CalculateMeltingTemperature` + `CalculateSaltCorrection` (INV-3 additivity holds). Threshold boundary verified: 13 bp → Wallace (M6=38), 14 bp → Marmur-Doty (M7≈37.36). Consistent.

### Test quality audit
34+ tests assert exact sourced values (not "no-throw"): exact 16.0/32.0/24.0/12.0/2.0/4.0/38.0 for Wallace; exact formula values within 0.1 for Marmur-Doty; exact salt corrections; direct `ThermoConstants_*` constant assertions; edge cases (empty/null/N/U/all-invalid/case). Deterministic. Covers all Stage-A edge cases.

Stage B verdict: **PASS-WITH-NOTES** — code faithfully realises the validated formulas with exact constants; same two documented notes as Stage A (label swap, salt reference state) carry over but are honestly disclosed in spec.

---

## Verdict & follow-ups

- **STATE: CLEAN.** No defect. All three implemented Tm methods are simple, well-published formulas, honestly labeled as simplifications (the spec explicitly states Marmur-Doty is a simplified alternative to NN and that the −7 Wallace correction is intentionally omitted). Constants match authoritative sources exactly; all worked examples reproduce; build green; 41/41 unit tests and 4461/4461 full-suite tests pass.
- **Not a LIMITED case:** the library does not *advertise* a SantaLucia NN model in code and then mis-implement it — it offers simpler models, correctly. No wrong NN constants exist to flag.
- **Follow-up suggestions (non-blocking, for future accuracy work, not defects):**
  1. Salt correction is referenced to 1 M (absolute term), giving large negative offsets at PCR salt; if an *absolute* salt-adjusted Tm is desired, switch to a difference form relative to a 50 mM reference, or adopt the full salt-adjusted GC% form already present (`CalculateSaltAdjustedTm`). Document/decide intent.
  2. For real primer-design accuracy, a SantaLucia (1998) nearest-neighbor path with the unified ΔH°/ΔS° table, R·ln(C_T/4), and the 273.15 conversion would be the rigorous upgrade — currently absent by design.

No code changed.
