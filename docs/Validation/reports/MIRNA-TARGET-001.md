# Validation Report: MIRNA-TARGET-001 — microRNA Target Site Prediction

- **Validated:** 2026-06-24   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.FindTargetSites(mRna, miRna, minScore)`; supporting `AlignMiRnaToTarget`, `CreateTargetSite` (internal), `CalculateTargetScore` (internal), `GetReverseComplement`, `GetSeedSequence`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs` (39 tests)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (efficacy scoring is honestly-scoped seed-type-proportional, not full TargetScan context++; ΔG model is simplified — both documented simplifications, not defects)

---

## Stage A — Description

### Sources opened & what they confirm
- **TargetScan 8.0 help / FAQ** (targetscan.org/vert_80/docs/help.html) — confirmed verbatim the four canonical site definitions and the position-1/position-8 logic:
  - **8mer** = exact match to miRNA positions **2–8** *plus* an **A** at target position 1 (i.e. both an adenine opposite miRNA nt1 **and** WC pairing at nt8).
  - **7mer-m8** = exact match to positions **2–8** (WC pairing at nt8, no A1).
  - **7mer-A1** = exact match to positions **2–7** plus an **A** at target position 1 (no nt8 pairing).
  - **6mer** = exact match to positions **2–7** only (neither A1 nor nt8 pairing).
  - **Offset 6mer** = match to positions **3–8** (marginal efficacy; centered sites removed in TargetScan 8.0).
- **A1 rule confirmed precisely**: the A1 is an **adenosine in the mRNA** across from miRNA nt1, **independent of the miRNA nt1 identity** (it is an anchor, not a WC pair). The code's `hasA1 = mrna[i+6] == 'A'` (pure identity test, no complementarity) matches this exactly.
- **Grimson et al. (2007)** (PMID 17612493) — site-type efficacy weights **8mer = 0.31, 7mer-m8 = 0.161, 7mer-A1 = 0.099**; efficacy hierarchy 8mer > 7mer-m8 > 7mer-A1 > 6mer.
- **Bartel (2009)** (PMID 19167326) and **Agarwal et al. (2015)** (PMID 26267216) — seed = positions 2–8, antiparallel pairing, 6mer/offset-6mer marginal-but-detectable efficacy. Evidence doc citations are accurate.

### Seed-complementarity rule (the critical biological check)
Target = **reverse complement of the seed**, sought antiparallel in the mRNA 5′→3′. This is the catastrophic-error checkpoint. Hand computation for **hsa-let-7a-5p** `UGAGGUAGUAGGUUGUAUAGUU`:
- Seed (pos 2–8) = `GAGGUAG` (= `Sequence.Substring(1,7)`). ✓
- Reverse complement = `CUACCUC`.
- 6mer core (pos 2–7 RC) = `UACCUC` = `seedRC[1..7]`. ✓
- Offset-6mer (pos 3–8 RC) = `CUACCU` = `seedRC[0..6]`. ✓
- 8mer = `CUACCUCA`; 7mer-m8 = `CUACCUC`; 7mer-A1 = `UACCUCA`; 6mer = `UACCUC` — match spec/Evidence worked example exactly.

### Edge-case semantics
Empty/null → empty; mRNA shorter than seed (<6) → empty; multiple/overlapping sites independent; DNA (T) normalized to RNA (U); minScore filter; scores clamped [0,1]; 8mer dominates 7mer-m8 at same position (both `hasPos8 && hasA1`). All defined and sourced.

### Independent cross-check (numbers)
Ran the actual code via a temporary probe (since removed). seed=`GAGGUAG`, seedRC=`CUACCUC`. Reproduced hand-computed classification + coordinates and the strict score chain:

| Construct (let-7a) | Type | Start | End | Len | Score |
|---|---|---|---|---|---|
| `…CUACCUCA…` | Seed8mer | 5 | 12 | 8 | 0.910 |
| `…CUACCUCG…` | Seed7merM8 | 5 | 11 | 7 | 0.430 |
| `…G UACCUC A…` | Seed7merA1 | 6 | 12 | 7 | 0.240 |
| `…G UACCUC G…` | Seed6mer | 6 | 11 | 6 | 0.070 |
| `…G CUACCU G…` | Offset6mer | 6 | 11 | 6 | 0.040 |
| `…GAGGUAG…` (identity) | **none** | — | — | — | — |
| short `CUAC` | none | — | — | — | — |

Score chain strictly monotone (0.910 > 0.430 > 0.240 > 0.070 > 0.040) per Grimson/Bartel hierarchy. **Verbatim seed `GAGGUAG` produced zero sites** — confirming targeting is by reverse complement, not identity. Start/End are 0-based inclusive.

**Stage A finding:** No biological or definitional error. PASS.

---

## Stage B — Implementation

### Code path reviewed
`FindTargetSites` (MiRnaAnalyzer.cs:154–264), `CreateTargetSite` (266–285), `GetReverseComplement` (294–317), `AlignMiRnaToTarget` (361–417), `CalculateTargetScore` (453–480).

### Formula realised correctly? (evidence)
- **Reverse complement, not identity** — `seedRC = GetReverseComplement(miRna.SeedSequence)`; pass-1 scans for `sixmerCore = seedRC[1..7]`, pass-2 for `offset6Pat = seedRC[0..6]`. Verbatim-seed probe returns empty. ✓ (no catastrophic error)
- **Site classification** — nt8 tested upstream at `mrna[i-1] == seedRC[0]` (RC of miRNA nt8); A1 tested downstream at `mrna[i+6] == 'A'` (pure identity). Antiparallel geometry correct (miRNA nt1 at 3′ end of target site, on the mRNA at `i+6`). 8mer requires **both** (`hasPos8 && hasA1`); 6mer requires neither. SeedMatchLengths 8/7/7/6/6 and Start/End offsets verified by probe. ✓
- **Offset-6mer (pass 2)** suppressed when overlapping a higher-priority pass-1 site (`coveredPositions`) or when it is part of a full seedRC (`mrna[i+6] == seedRC[6]`) — avoids double-counting. ✓
- **Scoring** (`CalculateTargetScore`) — base scores 1.0/0.52/0.32/0.15/0.10. The 0.52 and 0.32 are exactly Grimson (2007) weights normalized to 8mer=1.0 (0.161/0.310≈0.519, 0.099/0.310≈0.319). Mismatch penalty (−0.01/mismatch) and a >10-match 3′-supplementary bonus (+0.05) are minor adjustments; output clamped to [0,1]. ✓

### Cross-verification table recomputed vs code
See the Stage-A table above — recomputed directly from the running code; every value matches the hand computation and the spec/Evidence worked example. miR-21 8mer (`AUAAGCUA`, seedRC `AUAAGCU`) likewise classifies Seed8mer (covered by S-003 test).

### Variant/delegate consistency
`AlignMiRnaToTarget` reads the target antiparallel (`target[Length-1-i]`); perfect complement → all WC `|`; G:U → `:` wobble; same-base → mismatch; empty → empty duplex. `CanPair`/`IsWobblePair` normalize DNA T→U consistently. `CalculateDuplexEnergy` sums Turner-2004 NN stacking terms over consecutive paired positions (StackingEnergies dictionary spot-checked against the documented NNDB turner04 set).

### Test quality audit
27 tests assert **exact** site types, exact SeedMatchLengths, exact counts (`Count == 1` / `== 3`), strict score inequalities (M-005), and score-range bounds — not tautologies or "no throw". They cover all 5 site types, identity-vs-RC implicitly (M-007 no-match), DNA/U handling (M-017), minScore filtering (M-010), empty inputs (M-006/M-014), multiplicity (M-008), alignment counts/wobble/mismatch (M-011/12/13), ΔG<0 for a perfect duplex (M-015), and both reference miRNAs (S-003).

### Findings / notes (honestly-scoped, NOT defects)
1. **Default efficacy `Score` is seed-type-proportional only**, not the full TargetScan context++ model. An **opt-in** `ScoreTargetSiteContextPlusPlus` provides the source-fitted context++ score (Agarwal et al. 2015, eLife 4:e05005) using the **verbatim coefficients** from the TargetScan distribution `Agarwal_2015_parameters.txt` and the feature/scaling logic ported from `targetscan_70_context_scores.pl`. As of 2026-06-24 it realises every feature derivable from the miRNA + 3′UTR: site-type Intercept, Local_AU, sRNA position-1/8 identity, target site-position-8 identity (7mer-A1/6mer), **3′ supplementary pairing `3P_score`** (faithful DP port of `get3primePairingContribution`, raw scores cross-checked against the reference perl), **`Min_dist`** (log10 distance to the nearest 3′UTR end), **`Len_3UTR`** (log10 3′UTR length), and **`Off6m`** (offset-6mer count, used raw). The data-blocked features `SPS`, `TA_3UTR`, `Len_ORF`, `ORF8m` are computed faithfully only when supplied via the optional `ContextPlusPlusInputs`; `SA` (RNAplfold partition-function accessibility — deliberately NOT approximated by the library's MFE folder) and `PCT` (multi-species conservation) remain honest residuals. The still-residual set is reported in `OmittedFeatures`, so `ContextScorePartial` is a partial CS. Default `Score` unchanged. Honestly scoped per spec Simplification #3 and `LIMITATIONS.md`.
2. **ΔG (`FreeEnergy`) for target duplexes uses a simplified stacking-only model** (no loop/bulge/coaxial terms), so short G-flanked constructs may yield small positive ΔG. Classification depends on `Score`, not ΔG; M-015 only requires ΔG<0 for a fully complementary duplex, which holds. Honestly scoped per spec Simplification #1.

---

## Verdict & follow-ups
- **Stage A: PASS** — site definitions, the A1 anchor rule, the seed reverse-complement rule, and the efficacy hierarchy all confirmed against TargetScan 8.0, Bartel (2009), Agarwal (2015), Grimson (2007).
- **Stage B: PASS-WITH-NOTES** — implementation correctly targets by seed **reverse complementarity** (not identity), with correct nt8/A1 rules, correct site-type classification, correct antiparallel geometry, correct 0-based inclusive coordinates, and Grimson-proportional scoring. The two notes are intentional, documented simplifications, not biological errors.
- **State: RE-VALIDATION PENDING (Status remains ☐ in the root registry)** — the opt-in TargetScan context++ scorer (`ScoreTargetSiteContextPlusPlus`, Agarwal 2015) now computes the additional miRNA+3′UTR-derivable features (`3P_score`, `Min_dist`, `Len_3UTR`, `Off6m`) and accepts the data-blocked features (`SPS`, `TA_3UTR`, `Len_ORF`, `ORF8m`) as optional caller inputs. Covered by evidence-based tests CTX-001..CTX-011 whose exact expected values are hand-derived from the verbatim `Agarwal_2015_parameters.txt` coefficients and the `targetscan_70_context_scores.pl` formulas (3P_score raw values cross-checked against the reference perl). No defect found in the existing path; the default Grimson-proportional `Score` is unchanged. SA and PCT remain honest residuals.
- Follow-ups: independent re-validation of the extended opt-in method; SA requires a true McCaskill partition-function accessibility (not in scope), PCT requires a multi-species alignment.
