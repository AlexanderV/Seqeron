# Validation Report: PROBE-DESIGN-001 — Hybridization Probe Design

- **Validated:** 2026-06-12   **Area:** MolTools
- **Canonical method(s):** `ProbeDesigner.DesignProbes`, `DesignTilingProbes`, `DesignAntisenseProbes`, `DesignMolecularBeacon`, `ScoreProbe` (via `EvaluateProbe`/`EvaluateProbeWithGc`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia *Hybridization probe* | Hybridization probes are "usually 15–10000 nucleotides long"; probes detect **complementary** target sequences; stringency (temperature/salt) controls hybridization specificity. |
| Wikipedia *DNA microarray* | Oligonucleotide microarray probes: Affymetrix **25-mer**, Agilent **60-mer**; longer probes are more specific to individual genes. |
| Wikipedia *Fluorescence in situ hybridization* | BAC probes "on the order of 100 thousand base-pairs … basis for most FISH probes"; oligo FISH probes "often 10–25 nucleotides"; smFISH/RNA-FISH uses "~20–50 oligonucleotide pairs"; "shorter probes hybridize less specifically than longer probes". |
| Wikipedia *Molecular beacon* | Loop = **18–30 nt** region complementary to target; stem = **two 5–7 nt** complementary oligos (one each terminus); total ≈ 25 nt. |
| *Nucleic acid thermodynamics* / SantaLucia (1998), Wallace rule references | Tm of short oligos via Wallace rule **Tm = 2(A+T) + 4(G+C)** (rule-of-thumb up to ~14 nt); salt-adjusted formula **Tm = 81.5 + 16.6·log₁₀[Na⁺] + 0.41·(%GC) − 600/N** for longer oligos (length factor 600 or 675 both in use). |

### Constraint ranges confirmed (probe-design constraints + accepted ranges)

- **Probe length:** application-specific, all within the 15–10000 nt envelope. Microarray 50–60 (Agilent-class 60-mer regime), qPCR 20–30, FISH 200–500, Northern 100–300, Southern 150–500, molecular-beacon loop 25 (stem 5). All consistent with sources.
- **Tm target window:** application-specific (e.g. qPCR 68–72, microarray 75–85). Tm computed via PRIMER-TM-001 logic (`ThermoConstants`): Wallace for length < 14, salt-adjusted GC formula otherwise — both standard formulas confirmed above.
- **GC content window:** 0.35–0.65 fractional bands (microarray 0.40–0.60), within mathematical 0..1 and standard practice (avoid AT-rich/GC-rich extremes).
- **Secondary structure / repeats / runs:** homopolymer-run cap, inverted-repeat hairpin detection, di/tri-nucleotide simple-repeat detection — all standard probe quality filters, sourced as general practice.
- **Specificity / uniqueness:** `CheckSpecificity` via suffix tree (unique = 1.0, N hits = 1/N); `ValidateProbe` off-target via approximate matching. Consistent with "shorter probes hybridize less specifically" / uniqueness requirement.
- **Probe-target complementarity:** a probe hybridizes to its **reverse complement** target. `DesignAntisenseProbes` correctly designs against `GetReverseComplementString(mRNA)`; `CalculateSelfComplementarity` compares base i to RC base i (palindrome → 1.0).

### Edge-case semantics
- No valid probe / empty / null / too-short → empty (`yield break`). Defined.
- GC at window edges, all-GC → 1.0, all-AT → 0.0. Mathematically defined.
- Molecular beacon on too-short target → null. Defined.

### Independent cross-check (numbers, hand-computed)
- **Wallace example:** probe "ACGTACGT" (4 AT, 4 GC, len 8 < 14) → 2·4 + 4·4 = **24 °C**.
- **Salt-adjusted example (M13, all-G, N=100, gc=1.0, [Na⁺]=0.05):** 81.5 + 16.6·log₁₀(0.05) + 41·1.0 − 600/100 = 81.5 − 21.597 + 41 − 6 = **94.90 °C** (Tm > 0 ✓).
- **Molecular beacon (S4, stemLength=5):** stem5 = "GG"+"CCC" = "GGCCC"; stem3 = RC("GGCCC") = "GGGCC"; total = 5+20+5 = **30** ✓ (matches test assertions exactly).

### Findings / divergences
- None affecting correctness. The TestSpec table label "M11 — length 50-60 bp" vs. an inline comment in the test "50-70 bp" is a stale comment only; the assertion uses the live `param.MinLength/MaxLength` (50/60), so no behavioural divergence.

## Stage B — Implementation

### Code path reviewed
- `ProbeDesigner.cs` (`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs`): `DesignProbes` (132–195), `DesignProbesOptimized` (201–248), `EvaluateProbeWithGc` (253–327), `DesignTilingProbes` (341–397), `DesignAntisenseProbes` (402–414), `DesignMolecularBeacon` (419–473), `ValidateProbe`/`CheckSpecificity` (491–587), helpers (676–804).
- `ThermoConstants.cs`: Wallace + salt-adjusted Tm (linked to PRIMER-TM-001).

### Formula realised correctly?
- **Length:** scan `MinLength..MaxLength`, bounded by `n`; probe coordinates `[start, start+len-1]`, substring matches (M7/M8 verified). ✓
- **GC:** O(1) prefix-sum GC fraction; early reject outside `[MinGc-0.1, MaxGc+0.1]`; bounds penalty in scoring. ✓
- **Tm:** `CalculateTm` dispatches Wallace (<14) vs salt-adjusted GC formula — matches validated formulas. ✓
- **Complementarity:** reverse-complement used for antisense design, self-complementarity, and hairpin stem matching. ✓
- **Specificity:** suffix-tree unique=1.0 / N-hit=1/N. ✓

### Cross-verification table recomputed vs code
| Item | Hand value | Code/test outcome |
|------|-----------|-------------------|
| Beacon length (stem5+loop+stem3) | 30 | Test S4 asserts 30, stem "GGCCC"/"GGGCC" — pass |
| All-G probe GcContent | 1.0 | M13 asserts exactly 1.0 — pass |
| All-AT probe GcContent | 0.0 | M14 asserts exactly 0.0 — pass |
| Tiling starts (len50,ov10,208bp) | 0,40,80,120; coverage 170 | M9 asserts same — pass |
| Palindrome self-comp "AACCGGTT" | 1.0 | C2 asserts 1.0 — pass |

### Variant/delegate consistency
- `EvaluateProbe` delegates to `EvaluateProbeWithGc` with `CalculateGcContent`; optimized path supplies prefix-sum GC — same result. Tiling and antisense reuse the same evaluator. Consistent.

### Test quality audit
- 29 PROBE-DESIGN tests (full Probe filter = 72, incl. validation). Assertions check exact sourced/computed values (exact GC 0.0/1.0, exact beacon string, exact tiling starts/coverage, palindrome=1.0), not just "no throw". Edge cases (empty/null/short/beacon-short) covered. Mutation-killing beacon tests present.

### Findings / defects
- None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN**. No code changes required.
- Build: succeeded, 0 warnings. Tests: Probe filter 72/72 pass; full `Seqeron.Genomics.Tests` suite **4461 passed, 0 failed**.
- Minor non-blocking: stale "50-70 bp" inline comment in the M11 test; assertion logic uses live params so behaviour is correct (left as-is to avoid churn).
