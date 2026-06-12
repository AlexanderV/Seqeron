# Validation Report: MIRNA-PRECURSOR-001 — Pre-miRNA Hairpin Detection

- **Validated:** 2026-06-12   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.FindPreMiRnaHairpins(sequence, minHairpinLength=55, maxHairpinLength=120, matureLength=22)`; internal helper `AnalyzeHairpin(sequence, matureLength)`; energy via `CalculateHairpinEnergy`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (25 tests)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

> Note: the checklist block (ALGORITHMS_CHECKLIST_V2.md §MIRNA-PRECURSOR-001) lists the
> canonical method as `FindPreMiRnas` / `ValidateHairpin`, neither of which exists. The actual
> public API is `FindPreMiRnaHairpins`, which is what the TestSpec, Evidence doc, and tests use.
> Checklist naming is stale; the TestSpec/Evidence/tests are internally consistent.

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia: MicroRNA** (fetched). Confirms: pre-miRNA hairpins are "about 70 nucleotides each";
  pathway pri-miRNA → Drosha/DGCR8 → pre-miRNA → Exportin-5 → Dicer → ~22 nt miRNA:miRNA* duplex;
  Drosha leaves a **2-nt 3' overhang**; mature miRNA **21–23 nt**; **stem-loop (hairpin)** with the
  mature miRNA in the stem; **5p/3p** arm nomenclature; guide vs passenger (star) strand selected by
  thermodynamic instability at the 5' end. All match the TestSpec/Evidence §1.2.
- **Bartel (2004) / Drosha-Dicer literature** (web search). Confirms: Drosha cleaves ~11 bp (one
  helical turn) into the stem; pre-miRNA ~60–70 nt (database range up to ~100–120); stem is long
  (~33 bp in Bartel) but **imperfect** with internal mismatches/bulges; terminal loop variable
  (~4–15 nt). Confirms Evidence §1.2 points 1–3.
- **Krol et al. (2004)** — cited for: stem length critical for processing, effective stem ~18–22 bp
  (with allowed mismatches), G:U wobbles valid. Matches the spec thresholds (stem ≥ 18 bp, wobble pairs counted).
- **Turner/Xia 1998 + NNDB (Turner 2004)** — confirmed as the canonical nearest-neighbor free-energy
  source. The stacking values hard-coded in `StackingEnergies` (e.g. GC/CG = −3.42, GG/CC = −3.26,
  AU/UA = −1.10, UA/AU = −1.33, AA/UU = −0.93) are the standard published Xia/Turner Watson-Crick
  stacking parameters. (NNDB HTML pages 404'd on direct fetch; values cross-checked against the
  widely-published Xia 1998 table.)

### Feature / threshold check (vs published pre-miRNA characteristics)

| Feature | Spec / code | Published | Match |
|---|---|---|---|
| Hairpin length | [55, 120] nt window; `AnalyzeHairpin` hard-rejects n<55 | ~60–70 nt typical, 55–120 in miRBase | ✓ |
| Min stem | ≥ 18 consecutive bp | ~18–22 bp effective (Krol) | ✓ |
| Loop size | 3–25 nt | ~3–15, up to ~25 (Bartel) | ✓ |
| Mature | ~22 nt, 5' arm | ~22 nt, one arm (Bartel 2009) | ✓ |
| Star | opposite arm, == mature length | passenger from opposite arm | ✓ |
| G:U wobble | counted as a pair (`CanPair`) | valid in RNA stems (Krol) | ✓ |
| Dot-bracket | `(` 5'-stem, `.` loop, `)` 3'-stem | standard notation | ✓ |
| MFE | Turner 2004 NN model | NNDB | ✓ |

### Independent cross-check (numbers) — worked example

`ValidHairpin57` = `GCAUAGCUAGCUAGCUAGCUAGCUA`+`GAAAUUU`+`UAGCUAGCUAGCUAGCUAGCUAUGC` (57 nt).
Re-implemented `AnalyzeHairpin` + `CalculateHairpinEnergy` in Python independently:
- maxStem = 57/2 − 5 = 23; consecutive end-pairs = 23 → stem = 23 (≥18 ✓); loop = 57 − 46 = **11** (3–25 ✓).
- Stacking sum over 22 pairs = **−49.10**; + loop-init(11) = +6.60 → −42.50; + terminal mismatch
  `CUAG` (closing C-G) = −1.00 → **−43.50**; outer (G-C) and closing (C-G) pairs are WC → no AU/GU
  penalty. **Total = −43.50 kcal/mol** — exactly matches test M10/M11.
- Rejection cases reproduced: `ShortStemHairpin55` (15 bp stem) rejected (stem<18); 30-nt-loop case
  rejected (loop>25); all-purine `NoComplementarity` rejected (no stem).

### Findings / divergences (Stage A)

1. **PASS-WITH-NOTES — INV-9 / M11 framing.** INV-9 says "longer stem → more negative FreeEnergy
   *all else equal*", but M11 compares stem-23/loop-11 against stem-20/loop-15 (stem **and** loop
   differ). The ordering still holds and is correctly asserted with exact values; only the
   "all else equal" wording is imprecise. Cosmetic.
2. **Evidence doc nit.** §"Dataset 3" claims hsa-mir-21 yields "only 8 consecutive pairs from ends".
   Independent recomputation gives **16** consecutive end-pairs (still < 18, so still rejected; and
   no sub-window passes both the stem≥18 and loop 3–25 constraints → 0 detections). Test outcome
   (M18 = not detected) is correct; the "8" figure in the Evidence doc is inaccurate.
3. Honest scope is correctly stated: this is a **simplified consecutive-pairing hairpin heuristic**,
   not a real RNA-folding pre-miRNA classifier (Zuker/Nussinov/miRDeep). M18/M19 explicitly lock in
   that real miRBase pre-miRNAs (hsa-mir-21, hsa-let-7a-1) are *not* detected. This is an honest,
   well-documented design limitation, not a hidden defect.

Stage A: biology and thresholds are correct and sourced. PASS-WITH-NOTES (items 1–2 are doc-level).

---

## Stage B — Implementation

### Code path reviewed
- `FindPreMiRnaHairpins` (MiRnaAnalyzer.cs:610) — null/empty/short guard, T→U + uppercase, O(n²)
  window scan over [min, min(max, remaining)].
- `AnalyzeHairpin` (:645) — n<55 guard; consecutive end-pair stem (maxStem = min(n/2−5, 35));
  stem≥18; loop ∈ [3,25]; mature = first min(matureLength, stemLength); star = last `matureEnd`;
  dot-bracket build; energy.
- `CalculateHairpinEnergy` (:556) — stacking + loop-init + terminal-mismatch + terminal AU/GU penalty.

### Formula realised correctly?
Yes. The Python re-implementation of both the structural logic and the full NN energy model
reproduced the canonical worked example to the cent (−43.50) and reproduced every accept/reject
decision in the test data. The stacking-key construction `seq[i]seq[i+1]/seq[n-1-i]seq[n-2-i]`
correctly indexes the antiparallel partner stack, and all 22 keys for the example are present in
the table (no silent misses).

### Cross-verification table recomputed vs code

| Case | Expected (spec) | Recomputed | Code/test |
|---|---|---|---|
| ValidHairpin57 energy | −43.50 | −43.50 | M10 ✓ |
| stem-20/loop-15 energy | −36.02 | (matches; uses UAGA tm −1.10 + U-A AU penalty +0.45) | M11 ✓ |
| ValidHairpin57 detections (min=55) | 2 (i=0 len57, i=1 len55) | 2 | M5 ✓ |
| mature (22) | `GCAUAGCUAGCUAGCUAGCUAG` | matches | M6 ✓ |
| star (22) | `CUAGCUAGCUAGCUAGCUAUGC` | matches | M7 ✓ |
| structure | 23×`(`+11×`.`+23×`)` | matches | M8 ✓ |
| ShortStem55 / 30-loop / no-comp | reject | reject | M12/M13/S4 ✓ |
| hsa-mir-21 / hsa-let-7a-1 | not detected | 0 windows | M18/M19 ✓ |

### Variant/delegate consistency
Single public entry point (`FindPreMiRnaHairpins`); no `*Fast`/delegate variants. `CanPair` and the
energy tables are shared with target-prediction code and behave consistently.

### Test quality audit
Tests assert **exact** values (energies, sequences, structure string, positions, counts), not just
"no-throw" — strong. Edge cases covered: null, empty, short, no-complementarity, stem<18, loop>25,
T→U, wobble, case-insensitive, max/min length, custom matureLength, real-miRBase non-detection.
INV-1..10 checked in M17. Notes:
- M13 correctly documents that the loop-too-**small** arm (<3) is structurally unreachable because
  maxStem = n/2−5 forces loop ≥ 10; only loop-too-large is testable. Honest.
- INV-9 wording vs M11 as in Stage A note 1 (cosmetic).

### Findings / defects (Stage B)
No code defect. Implementation faithfully realises the validated (heuristic) description; all
documented accept/reject and numeric values reproduced independently. The only items are the two
documentation imprecisions noted in Stage A and the stale checklist method name — none affect
behaviour or any test.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — biology/thresholds/MFE model correct and sourced; two doc-level
  inaccuracies (INV-9 "all else equal" framing; Evidence "8 pairs" should be 16 for hsa-mir-21).
- **Stage B:** PASS-WITH-NOTES — code matches the validated description; worked example and all
  reject cases reproduced to the cent; no code change required.
- **State: CLEAN** — no code defect. Honest, well-documented scope: a simplified consecutive-pairing
  hairpin heuristic (not a folding-based classifier); real pre-miRNAs intentionally not detected
  (locked by M18/M19).
- **Tests:** `~PreMiRna` filter = 28 passed / 0 failed; full suite = **4486 passed, 0 failed**.
- **Code changed:** none.
- **Optional follow-ups (non-blocking):** (a) correct the Evidence doc "8 consecutive pairs" → 16
  for hsa-mir-21; (b) soften INV-9 wording or make M11 vary only the stem; (c) fix the stale
  checklist method names (`FindPreMiRnas`/`ValidateHairpin` → `FindPreMiRnaHairpins`).
