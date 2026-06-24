# Validation Report: MIRNA-PRECURSOR-001 — Pre-miRNA Hairpin Detection

- **Validated:** 2026-06-24   **Area:** MiRNA
- **Updated:** 2026-06-24 — added opt-in MFE-structure-based detection (reuses the RNA-STRUCT-001 Zuker–Stiegler folder); Status reset to ☐ for re-validation.
- **Canonical method(s):** `MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength=55, maxHairpinLength=120, matureLength=22)` (default heuristic); internal helper `AnalyzeHairpin`; energy via `CalculateHairpinEnergy`. **New opt-in:** `FindPreMiRnaHairpinsByMfe`, `AssessHairpinByMfe`, `CalculateMfeIndex` — fold the candidate with `RnaSecondaryStructure.CalculateMfeStructure` (RNA-STRUCT-001) and derive hairpin features (single terminal loop, paired-stem count, ΔG°, AMFE, MFEI) from the real MFE structure.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (+ mutation-killer tests) — 38 heuristic tests + 10 MFE-fold tests (MF1–MF10) under the `~PreMiRna` filter
- **Stage A verdict:** PASS-WITH-NOTES (heuristic); MFE-fold opt-in pending re-validation
- **Stage B verdict:** PASS-WITH-NOTES (heuristic); MFE-fold opt-in pending re-validation

> **2026-06-24 limitation fix.** The previously documented limitation — "pre-miRNA hairpin detected by
> a consecutive-complementary-pairing heuristic, not a full MFE-structure-folding caller" — is removed
> by adding an opt-in MFE path. `AssessHairpinByMfe` folds the candidate with the already-validated
> Zuker–Stiegler engine (RNA-STRUCT-001, Turner 2004) and reads the hairpin from the **actual** MFE
> dot-bracket: a single dominant hairpin (one terminal loop, nested stem, no multibranch) with stem
> base pairs ≥ 16 (Ambros 2003), terminal loop 3–25 nt (Bartel 2004), and MFEI ≥ 0.85 (Zhang 2006,
> MFEI = AMFE/(G+C)%, AMFE = 100·|ΔG°|/length). ΔG° is the engine's `CalculateMinimumFreeEnergy`
> verbatim. The MFE path **detects the real miRBase precursors the heuristic rejected** — hsa-mir-21
> (ΔG°=−35.13, 32 bp, MFEI 1.0037) and hsa-let-7a-1 (ΔG°=−34.31, 32 bp, MFEI 1.0091) — while
> rejecting non-complementary (ΔG°=0) and multibranch (5S-like ΔG°=−47.04) structures. The default
> heuristic is unchanged.

> **2026-06-24 limitation fix (cleavage-site prediction).** Added the opt-in
> `PredictDroshaDicerCleavage(sequence, basalJunction)` implementing the PUBLISHED measuring ("ruler")
> rules: Drosha cleaves **~11 bp** from the basal ssRNA–dsRNA junction (Han et al. 2006, Cell 125:887,
> "approximately 11 bp from the stem-ssRNA junction"); Dicer cleaves **~22 nt** from the
> Drosha-generated 5' end — the 5' counting rule (Park et al. 2011, Nature 475:201, "∼22 nucleotides
> from the 5' end"), fixing the mature length at 22 nt; each RNase III cut leaves a **2-nt 3' overhang**
> (Lee 2003 / Han 2006); an optional CNNC confidence flag is reported when a C-N-N-C sits 16–18 nt 3'
> of the Drosha cut (Auyeung et al. 2013). The method reports the Drosha cut, the mature (5p) and star
> (3p) coordinates, and the 2-nt overhangs. **Cross-checked** against miRBase hsa-miR-21-5p
> (MIMAT0000076): with an 11-nt lower stem prepended so the +11 ruler lands at the annotated 5p start,
> the predicted mature equals `UAGCUUAUCAGACUGAUGUUGA` exactly (DD4). Defaults and existing methods are
> unchanged. New tests DD1–DD10 in `MiRnaAnalyzer_PreMiRna_Tests.cs`. **Honest residual (data-blocked):**
> a competitive **trained** natural-vs-background precursor classifier (miRDeep2-style fitted
> probabilistic model) — requires a trained model and labelled data. Status remains ☐ (registry
> unchanged; not previously ☑).

> Note: the checklist block (`ALGORITHMS_CHECKLIST_V2.md` §MIRNA-PRECURSOR-001) lists the canonical
> method as `FindPreMiRnas` / `ValidateHairpin`, neither of which exists. The actual public API is
> `FindPreMiRnaHairpins`, which is what the TestSpec, Evidence doc, and tests use. Checklist naming is
> stale; the TestSpec/Evidence/tests are internally consistent. (Non-blocking.)

This is a re-validation in a fresh context. Findings match the prior archived report (`git show
cb113ce:…`); all numeric cross-checks were independently re-derived rather than copied.

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia: MicroRNA** (fetched 2026-06-24). Confirms: pre-miRNA hairpins are "about 70 nucleotides
  each"; mature miRNA 21–23 nt; miRNAs derive from RNA that folds back into a **stem-loop (hairpin)**;
  **5p/3p** arm nomenclature for the two mature products from opposite arms; the passenger (star)
  strand occupies the **opposite arm** and is selected against by thermodynamic stability. All match
  TestSpec §1.2.
- **Ambros et al. (2003), RNA 9:277-279** (the "uniform system" paper, via search + PMC). Confirms the
  annotation criterion that the putative miRNA is embedded in **one arm of a fold-back hairpin** with
  **~16+ complementary bases** to the opposite arm. This directly supports the spec's "mature on one
  arm" + a meaningful stem-complementarity requirement.
- **Bartel (2004)/(2009), Krol (2004)** — cited for: pre-miRNA ~60–70 nt (DB range up to ~120);
  effective stem ~18–22 bp with allowed mismatches; loop ~3–15 nt (up to ~25); ~22 nt mature from one
  arm; G:U wobbles are valid RNA stem pairs. Consistent with the spec thresholds (stem ≥18 bp, loop
  3–25, mature ~22 nt, wobble counted).
- **Turner/Xia 1998 + NNDB (Turner 2004)** — nearest-neighbor free-energy source. The hard-coded
  `StackingEnergies` (GC/CG = −3.42, GG/CC = −3.26, CG/GC = −2.36, AU/UA = −1.10, UA/AU = −1.33,
  AA/UU = −0.93, etc.) are the standard published Xia/Turner Watson-Crick stacking parameters.

### Feature / threshold check (spec/code vs published)

| Feature | Spec / code | Published | Match |
|---|---|---|---|
| Hairpin length | [55, 120] nt window; `AnalyzeHairpin` rejects n<55 | ~60–70 typical, 55–120 in miRBase | ✓ |
| Min stem | ≥ 18 consecutive bp | ~18–22 effective (Krol); ≥16 complementary (Ambros) | ✓ |
| Loop size | 3–25 nt | ~3–15, up to ~25 (Bartel) | ✓ |
| Mature | ~22 nt, 5' arm | ~22 nt, one arm (Bartel 2009) | ✓ |
| Star | opposite arm, == mature length | passenger from opposite arm | ✓ |
| G:U wobble | counted as a pair (`CanPair`) | valid in RNA stems (Krol) | ✓ |
| Dot-bracket | `(` 5'-stem, `.` loop, `)` 3'-stem | standard notation | ✓ |
| MFE | Turner 2004 NN model | NNDB | ✓ |

### Independent cross-check (numbers)

`ValidHairpin57` = `GCAUAGCUAGCUAGCUAGCUAGCUA`+`GAAAUUU`+`UAGCUAGCUAGCUAGCUAGCUAUGC` (57 nt).
Re-implemented `AnalyzeHairpin` + `CalculateHairpinEnergy` in Python this session:
- maxStem = 57/2 − 5 = 23; consecutive end-pairs = 23 → stem = **23** (≥18 ✓); loop = 57 − 46 = **11** (3–25 ✓).
- Stacking sum over 22 pairs = **−49.10**; + loop-init(11) = +6.60 → −42.50; + terminal mismatch
  `CUAG` = −1.00 → **−43.50**; both outer (G-C) and closing (C-G) pairs are WC, so no AU/GU penalty.
  **Total = −43.50 kcal/mol** — exactly matches tests M10/M11.
- **hsa-mir-21** (MI0000077, 72 nt) re-run through the full scan: **0 detections** (16 consecutive
  end-pairs at the full window < 18; no sub-window satisfies both stem≥18 and loop 3–25). Matches M18.

### Findings / divergences (Stage A)

1. **PASS-WITH-NOTES — INV-9 / M11 framing.** INV-9 reads "longer stem → more negative FreeEnergy",
   but M11 compares stem-23/loop-11 vs stem-20/loop-15 (stem *and* loop differ). The ordering still
   holds and is asserted with exact values; only the "all else equal" implication is loose. Cosmetic.
2. **Checklist method-name staleness** (see header note). Non-blocking.
3. Scope is honestly stated: a **simplified consecutive-pairing hairpin heuristic**, not a
   folding-based pre-miRNA classifier (Zuker/Nussinov/miRDeep). M18/M19 lock that real miRBase
   pre-miRNAs (hsa-mir-21, hsa-let-7a-1) are *not* detected — an honest, documented design limitation.

Stage A: biology, thresholds, and the MFE model are correct and sourced. **PASS-WITH-NOTES.**

---

## Stage B — Implementation

### Code path reviewed
- `FindPreMiRnaHairpins` (MiRnaAnalyzer.cs:629) — null/empty/short guard, `ToUpperInvariant()` + T→U,
  O(n²) window scan over [min, min(max, remaining)].
- `AnalyzeHairpin` (:664) — n<55 guard; consecutive end-pair stem (maxStem = min(n/2−5, 35)); stem≥18;
  loop ∈ [3,25]; mature = first min(matureLength, stemLength); star = last `matureEnd`; dot-bracket
  build; energy.
- `CalculateHairpinEnergy` (:575) — stacking + loop-init (+extrapolation for >30) + terminal mismatch
  (loops ≥4) + terminal AU/GU penalty (0.45).
- `CanPair` (:323) — Watson-Crick + G:U wobble, case-insensitive, DNA T→U normalized. Shared with
  target-prediction code; behaves consistently.

### Formula realised correctly?
Yes. The Python re-implementation of the structural logic and the full NN energy model reproduced the
canonical worked example to the cent (−43.50). The stacking-key construction
`seq[i]seq[i+1]/seq[n-1-i]seq[n-2-i]` correctly indexes the antiparallel partner stack; all 22 keys for
the example are present (no silent misses). Every accept/reject decision in the test data reproduced.

### Cross-verification table recomputed vs code

| Case | Expected | Recomputed (this session) | Test |
|---|---|---|---|
| ValidHairpin57 energy | −43.50 | −43.50 | M10 ✓ |
| ValidHairpin57 detections (min=55) | 2 (i=0 len57, i=1 len55) | 2 | M5 ✓ |
| mature (22) | `GCAUAGCUAGCUAGCUAGCUAG` | matches | M6 ✓ |
| star (22) | `CUAGCUAGCUAGCUAGCUAUGC` | matches | M7 ✓ |
| structure | 23×`(`+11×`.`+23×`)` | matches | M8 ✓ |
| stem<18 / loop>25 / no-comp | reject | reject | M12/M13/S4 ✓ |
| hsa-mir-21 / hsa-let-7a-1 | not detected | 0 windows (mir-21: 16 end-pairs) | M18/M19 ✓ |

### Variant/delegate consistency
Single public entry point; no `*Fast`/delegate variants. Energy tables and `CanPair` shared with
target-prediction and consistent.

### Test quality audit
Tests assert **exact** values (energies, sequences, structure string, positions, counts) — not
"no-throw" tautologies. Edge cases covered: null, empty, short, no-complementarity, stem<18, loop>25,
T→U, wobble, case-insensitive, max/min length, custom matureLength, real-miRBase non-detection;
INV-1..10 in M17. The loop-too-**small** arm (<3) is structurally unreachable because maxStem = n/2−5
forces loop ≥ 10, so only loop-too-large is testable — correctly documented. `~PreMiRna` filter =
**38 passed / 0 failed**.

### Findings / defects (Stage B)
No code defect. The implementation faithfully realises the validated (heuristic) description; the only
items are the two documentation imprecisions in Stage A and the stale checklist method name — none
affect behaviour or any test.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — biology/thresholds/MFE model correct and sourced; INV-9 "all else
  equal" framing is loose; checklist method names stale.
- **Stage B:** PASS-WITH-NOTES — code matches the validated description; worked example and all reject
  cases reproduced independently to the cent; no code change required.
- **State: CLEAN** — no code defect. Honest, well-documented scope (simplified consecutive-pairing
  hairpin heuristic; real pre-miRNAs intentionally not detected, locked by M18/M19). The Evidence
  doc's prior "8 consecutive pairs" nit for hsa-mir-21 is already corrected to **16** in
  `docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md`.
- **Tests:** `~PreMiRna` filter = 38 passed / 0 failed.
- **Code changed:** none.
- **Optional follow-ups (non-blocking):** (a) soften INV-9 wording or make M11 vary only the stem;
  (b) fix the stale checklist method names (`FindPreMiRnas`/`ValidateHairpin` → `FindPreMiRnaHairpins`).
