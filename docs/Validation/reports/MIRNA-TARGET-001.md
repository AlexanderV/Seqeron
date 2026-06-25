# Validation Report: MIRNA-TARGET-001 — microRNA Target Site Prediction

- **Validated:** 2026-06-25 (fresh re-validation; supersedes the 2026-06-24 entry)   **Area:** MiRNA
- **Scope:** the OWN canonical surface of MIRNA-TARGET-001 — seed-match detection + base/default scoring + site-type classification (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer). The opt-in TargetScan context++ scorer (`ScoreTargetSiteContextPlusPlus`) belongs to MIRNA-CONTEXT-001 / MIRNA-PCT-001 and is referenced but NOT re-litigated here.
- **Canonical method(s):** `MiRnaAnalyzer.FindTargetSites(mRna, miRna, minScore)` (canonical); supporting internal `CreateTargetSite`, `CalculateTargetScore` (the "base scorer"), `CalculateDuplexEnergy`; public helpers `GetReverseComplement`, `AlignMiRnaToTarget`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs` (M-/S-/E-/CTX- regions); alignment helpers in `MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs`.
- **Stage A verdict:** PASS (✅)
- **Stage B verdict:** PASS-WITH-NOTES (🟡)
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened THIS session
- **Bartel DP (2009) Cell 136:215–233, "MicroRNAs: Target Recognition and Regulatory Functions"** — open-access mirror PMC3794896. Confirms verbatim:
  - **Seed = miRNA nucleotides 2–7** (nt8 is supplementary).
  - **6mer** = perfect WC match to miRNA nt 2–7.
  - **7mer-m8** = seed match (2–7) augmented by a WC match to miRNA nt8.
  - **7mer-A1** = seed match (2–7) augmented by an **A at target position 1**.
  - **8mer** = seed match flanked by **both** the nt8 match **and** the A at target position 1.
  - **Efficacy hierarchy (quoted): "8mer >> 7mer-m8 > 7mer-A1 >> 6mer > no site"** — 6mer barely above background.
- **Friedman/Agarwal context (PMC4062300)** — re-confirms the seed is nt 2–7, the A1 is "an adenosine opposite miRNA nucleotide 1", and m8 is "a Watson–Crick pair with miRNA nucleotide 8".
- **Lewis BP et al. (2005) Cell 120:15–20** — target sites are the **reverse complement** of the seed, sought antiparallel along the 3′UTR 5′→3′. This is the catastrophic-error checkpoint.
- **Grimson et al. (2007)** site-type efficacy weights 8mer = 0.31, 7mer-m8 = 0.161, 7mer-A1 = 0.099 — the basis of the base scorer's normalized constants.

### A1 anchor rule (critical biological check)
The A1 is an **adenosine at target position 1** (opposite miRNA nt1), recognized as an anchor by Argonaute; in the TargetScan convention it is a **pure identity test on the target base = 'A'**, independent of miRNA nt1 and NOT a Watson–Crick pair. The code tests `hasA1 = mrna[i+6] == 'A'` (identity, no complementarity) — matches the convention exactly.

### Seed reverse-complement rule
Hand computation for **hsa-let-7a-5p** `UGAGGUAGUAGGUUGUAUAGUU`:
- Seed (nt 2–8 stored as `SeedSequence`) = `GAGGUAG` (`Substring(1,7)`).
- Reverse complement `GetReverseComplement("GAGGUAG")` = `CUACCUC` (verified char-by-char: comp `CUCCAUC`, reversed `CUACCUC`).
- 6mer core (RC of nt 2–7) = `seedRC[1..7]` = `UACCUC`.
- offset-6mer pattern (RC of nt 3–8) = `seedRC[0..6]` = `CUACCU`.

### Edge-case semantics (all defined & sourced)
Empty/null mRNA or miRNA → empty; mRNA shorter than the 6-nt core → empty; verbatim seed identity (not RC) → empty; single seed mismatch / seed-region G:U wobble → no canonical site (exact-match scan); DNA T normalized to RNA U; non-ACGU (N) neither throws nor matches; minScore filter; both m8+A1 → 8mer (dominance); scores clamped to [0,1].

### Independent cross-check (exact numbers, run against the live code this session via a temporary probe)
seed `GAGGUAG`, seedRC `CUACCUC`, 6mer core `UACCUC`, offset `CUACCU`. Constructs `GGGGG <site> GGGGG`:

| Construct (let-7a) | Type | Start | End | Len | Score |
|---|---|---|---|---|---|
| `…CUACCUC A…` | Seed8mer | 5 | 12 | 8 | **0.910** |
| `…CUACCUC G…` | Seed7merM8 | 5 | 11 | 7 | **0.430** |
| `…UACCUC A…` | Seed7merA1 | 5 | 11 | 7 | **0.240** |
| `…UACCUC G…` | Seed6mer | 5 | 10 | 6 | **0.070** |
| `…CUACCU G…` | Offset6mer | 5 | 10 | 6 | **0.040** |
| `…GAGGUAG…` (verbatim seed identity) | **none** | — | — | — | — |
| short `CUAC` (4 nt) | none | — | — | — | — |

Score chain **strictly monotone** `0.910 > 0.430 > 0.240 > 0.070 > 0.040`, matching the Bartel/Grimson hierarchy 8mer >> 7mer-m8 > 7mer-A1 >> 6mer (> offset-6mer). **Verbatim seed produced ZERO sites**, confirming targeting is by reverse complement, not identity. Start/End are 0-based inclusive; `End − Start + 1 == SeedMatchLength`.

> **Note (coordinate correction):** the prior (2026-06-24) report listed the 7mer-A1 start as 6 and the 6mer at 6..11/etc.; the fresh probe shows **7mer-A1 = 5..11** and **6mer = 5..10** (the 6mer core begins at index 5 in `GGGGG`+core). The current code is correct; the earlier report's coordinates were stale. Locked by new test E-001.

**Stage A finding:** No biological or definitional error. **PASS.**

---

## Stage B — Implementation

### Code path reviewed
`FindTargetSites` (MiRnaAnalyzer.cs:156–266), `CreateTargetSite` (268–287), `GetReverseComplement` (296–319), `AlignMiRnaToTarget` (363–419), `CalculateDuplexEnergy` (428–453), `CalculateTargetScore` (455–482).

### Realises the validated description? (evidence)
- **Reverse complement, not identity** — `seedRC = GetReverseComplement(miRna.SeedSequence)`; pass-1 scans for the exact substring `sixmerCore = seedRC[1..7]`; pass-2 for `offset6Pat = seedRC[0..6]`. Verbatim-seed probe + test E-002 return empty. ✓ (no catastrophic error)
- **Site classification** — nt8 tested upstream `mrna[i-1] == seedRC[0]`; A1 tested downstream `mrna[i+6] == 'A'` (pure identity). 8mer requires **both**; 6mer neither. Antiparallel geometry correct (miRNA nt1 maps to the 3′ end of the target site). SeedMatchLengths 8/7/7/6/6 and Start/End verified by probe + E-001. ✓
- **Exact seed match** — the scan uses `Substring(i,6) != sixmerCore` (exact equality), so a single seed mismatch or a seed-region G:U wobble does NOT yield a canonical site (tests E-003). Matches Bartel 2009's "perfect seed complementarity" requirement. ✓
- **Offset-6mer (pass 2)** suppressed when overlapping a higher-priority pass-1 site (`coveredPositions`) or when it is part of a full seedRC. ✓
- **Base scorer** (`CalculateTargetScore`) — base scores 1.0 / 0.52 / 0.32 / 0.15 / 0.10 for 8mer/7m8/7A1/6mer/offset. The 0.52 and 0.32 are the Grimson (2007) weights normalized to 8mer = 1.0 (0.161/0.310 ≈ 0.519; 0.099/0.310 ≈ 0.319). A −0.01/mismatch penalty over the full-length antiparallel alignment and a >10-match 3′-supplementary +0.05 bonus are minor adjustments; output clamped [0,1]. The monotone ordering is preserved (cross-check above). ✓

### Variant/delegate consistency
`AlignMiRnaToTarget` reads the target antiparallel (`target[Length-1-i]`); WC → `|`, G:U → `:`, mismatch → space; empty → empty duplex. `CalculateDuplexEnergy` sums Turner-2004 NN stacking terms over consecutive paired positions. These are exercised by `MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs` (well-sourced, exact).

### Test quality audit
The MIRNA-TARGET-001 suite asserts **exact** site types, exact SeedMatchLengths, exact counts (`Count==1`/`==3`), strict score inequalities (M-005), score-range bounds (M-009), and now **exact 0-based-inclusive coordinates** (E-001). Expected values trace to Bartel 2009 / Lewis 2005 / hand-construction, NOT code echoes. This session ADDED 7 tests to close gate gaps:
- **E-001** exact coordinates per site type (locks the corrected 7mer-A1/6mer starts).
- **E-002** verbatim seed identity → no site (RC-not-identity checkpoint).
- **E-003** single seed mismatch / seed-region G:U wobble → no canonical site (exact-match requirement).
- **E-004** A1 + m8 both present → classified 8mer, not split (dominance).
- **E-005** mRNA shorter than the 6mer core → empty; non-ACGU (N) neither throws nor spuriously matches, and a valid 8mer among N flanks is still found.

Full class: 63 tests pass. Full unfiltered `dotnet test Seqeron.sln -c Debug`: **Seqeron.Genomics.Tests = 18799 passed, 0 failed**; whole solution green; 0 warnings on the changed file.

### Findings / notes (honestly-scoped, NOT defects)
1. **Default `Score` is a seed-type-proportional base score** (Grimson-normalized + small mismatch/3′ adjustments), not the full TargetScan context++ model. The opt-in `ScoreTargetSiteContextPlusPlus` (MIRNA-CONTEXT-001 / MIRNA-PCT-001) supplies the source-fitted context++ score and is validated under those units. Documented simplification, not a defect of THIS unit.
2. **ΔG (`FreeEnergy`) uses a stacking-only model** (no loop/bulge/coaxial terms), so short G-flanked constructs may show ΔG = 0; classification depends on `Score`, not ΔG. Documented simplification (MIRNA-PAIR-001 energy model).
3. **Registry naming inaccuracy (documentation, not code):** the catalog row (`ALGORITHMS_CHECKLIST_V2.md` §MIRNA-TARGET-001 and the method index) lists a public method `ScoreTargetSite(site)` that **does not exist** in the source — the scoring is the internal `CalculateTargetScore`, exercised through `FindTargetSites`. The algorithm is fully implemented and tested; only the registry label is imprecise. Logged in FINDINGS_REGISTER as a minor doc note. Not blocking CLEAN.

---

## Verdict & follow-ups
- **Stage A: PASS** — seed = nt 2–7, the four site-type definitions, the A1 anchor (target-base identity) rule, the m8 WC rule, the seed reverse-complement targeting rule, and the efficacy hierarchy all independently confirmed against Bartel (2009), Lewis (2005), Friedman/Agarwal, Grimson (2007).
- **Stage B: PASS-WITH-NOTES** — `FindTargetSites` targets by seed **reverse complementarity** (not identity), with correct exact-match seed scan, correct nt8/A1 classification, correct antiparallel geometry, correct 0-based inclusive coordinates, and a monotone Grimson-proportional base scorer. The three notes are intentional/documented (two model simplifications, one registry-label inaccuracy), not biological errors.
- **End-state: ✅ CLEAN** — no code defect; the test gate gaps were closed with 7 source/hand-derived edge tests; full suite green.
- Follow-up: the context++ / PCT extensions are validated under MIRNA-CONTEXT-001 / MIRNA-PCT-001.
