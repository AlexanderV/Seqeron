# Validation Report: MIRNA-TARGET-001 — microRNA Target Site Prediction

- **Validated:** 2026-06-12   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.FindTargetSites(mRna, miRna, minScore)`; supporting `AlignMiRnaToTarget`, `CreateTargetSite` (internal), `CalculateTargetScore` (internal), `GetReverseComplement`, `GetSeedSequence`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs` (27 tests)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (efficacy scoring is honestly-scoped, not full TargetScan context++; ΔG model is simplified — both documented, not defects)

---

## Stage A — Description

### Sources opened & what they confirm
- **Bartel (2009) Cell 136:215–233** (PMID 19167326) and **Agarwal et al. (2015) eLife** (PMID 26267216, PMC4532895): confirm the four/five canonical site definitions verbatim:
  - **8mer** = Watson–Crick (WC) match to miRNA positions **2–8** + adenosine opposite position 1.
  - **7mer-m8** = WC match to positions **2–8** only.
  - **7mer-A1** = WC match to positions **2–7** + A opposite position 1.
  - **6mer** = WC match to positions **2–7** only.
  - **Offset 6mer** = WC match to positions **3–8** (marginal efficacy).
- **A1 rule confirmed precisely**: the A1 is an adenosine in the *mRNA* across from miRNA nucleotide 1, **independent of the miRNA nucleotide identity** (Agarwal 2015 / TargetScan). It is an anchor, not a WC pair. The implementation's `hasA1 = mrna[i+6] == 'A'` (no complementarity test) matches this exactly.
- **Efficacy hierarchy confirmed**: 8mer > 7mer-m8 > 7mer-A1 > 6mer ≈ offset-6mer (Bartel 2009 "decreasing preferential conservation and efficacy"; Grimson 2007 fold-repression).
- **Grimson et al. (2007)** (PMID 17612493): site-type efficacy weights 8mer=0.310, 7mer-m8=0.161, 7mer-A1=0.099; context determinants (site-type, 3′-supplementary pairing, local AU content, UTR position).
- **TargetScan 8.0**: centered sites removed; 3 canonical + offset-6mer retained.

### Seed-complementarity rule (the critical biological check)
Target = **reverse complement of the seed**, sought antiparallel in the mRNA 5′→3′. This is the catastrophic-error checkpoint. Hand computation for **hsa-let-7a-5p** `UGAGGUAGUAGGUUGUAUAGUU`:
- Seed (pos 2–8) = `GAGGUAG` (= `Sequence.Substring(1,7)`). ✓
- Reverse complement = `CUACCUC`. Layout `[RCpos8=C][RCpos7=U][RCpos6=A][RCpos5=C][RCpos4=C][RCpos3=U][RCpos2=C]`.
- 6mer core (pos 2–7 RC) = `UACCUC` = `seedRC[1..7]`. ✓
- Offset-6mer (pos 3–8 RC) = `CUACCU` = `seedRC[0..6]`. ✓
- 8mer site = `CUACCUCA`; 7mer-m8 = `CUACCUC`; 7mer-A1 = `UACCUCA`; 6mer = `UACCUC` — match spec/Evidence worked example exactly.

### Edge-case semantics
Empty/null → empty; mRNA shorter than seed → empty; multiple/overlapping sites reported independently; DNA (T) normalized to RNA (U); minScore filter; scores in [0,1]; 8mer dominates 7mer-m8 at same position. All defined and sourced.

### Independent cross-check (numbers)
Direct probe of the code reproduced the hand-computed classifications and the strict score chain 0.91 > 0.43 > 0.24 > 0.07 for 8mer/7m8/7A1/6mer, and offset-6mer = 0.04. A verbatim (non-RC) seed `GAGGUAG` in the mRNA produced **zero** sites — confirming targeting is by reverse complement, not identity.

**Stage A finding:** No biological or definitional error. PASS.

---

## Stage B — Implementation

### Code path reviewed
`FindTargetSites` (MiRnaAnalyzer.cs:154–264), `CreateTargetSite` (266–285), `GetReverseComplement` (294–317), `AlignMiRnaToTarget` (350–404), `CalculateTargetScore` (434–461).

### Formula realised correctly? (evidence)
- **Reverse complement, not identity** — `seedRC = GetReverseComplement(seed)`; scanning matches `sixmerCore`/`offset6Pat` derived from `seedRC`. Verbatim-seed probe returns empty. ✓ (no catastrophic error)
- **Site classification** — pos8 match tested upstream at `mrna[i-1] == seedRC[0]` (RC of miRNA pos 8); A1 tested downstream at `mrna[i+6] == 'A'`. Antiparallel geometry correct (miRNA pos1 at 3′ end of the target site). 8mer requires **both** (`hasPos8 && hasA1`); 6mer requires neither. SeedMatchLengths 8/7/7/6/6 correct. ✓
- **Offset 6mer (pass 2)** suppressed when overlapping a higher-priority pass-1 site or when it is part of a full seedRC — avoids double-counting. ✓
- **Scoring** — base scores 1.0/0.52/0.32/0.15/0.10 are Grimson (2007) weights normalized to 8mer=1.0 (0.161/0.310≈0.52, 0.099/0.310≈0.32). Honestly proportional. Mismatch penalty and a >10-match supplementary bonus are minor adjustments; output clamped to [0,1].

### Cross-verification table recomputed vs code
| Construct (let-7a) | Expected type | Code type | Len | Score |
|---|---|---|---|---|
| `…CUACCUCA…` | 8mer | Seed8mer | 8 | 0.91 |
| `…CUACCUCG…` | 7mer-m8 | Seed7merM8 | 7 | 0.43 |
| `…UACCUCA…` | 7mer-A1 | Seed7merA1 | 7 | 0.24 |
| `…UACCUC…` | 6mer | Seed6mer | 6 | 0.07 |
| `…CUACCU…` | offset 6mer | Offset6mer | 6 | 0.04 |
| `…GAGGUAG…` (identity) | none | none | — | — |
| short `CUAC` | none | none | — | — |

Score chain strictly monotone, matching Grimson/Bartel hierarchy. miR-21 8mer (`AUAAGCUA`) likewise classified Seed8mer.

### Variant/delegate consistency
`AlignMiRnaToTarget`: perfect complement → all WC matches; G:U → wobbles; same-base → mismatches; empty → empty duplex. Antiparallel reading (`target[len-1-i]`) correct. DNA→RNA normalization consistent between `FindTargetSites` and `AlignMiRnaToTarget`.

### Test quality audit
27 tests assert exact site types, exact seed-match lengths, exact counts (`Count == 1`/`== 3`), strict score inequalities, and score-range bounds — not tautologies. Cover all 5 site types, identity-vs-RC implicitly (no false matches), DNA/U handling, minScore filtering, empty inputs, multiplicity, and both reference miRNAs (let-7a, miR-21).

### Findings / notes (honestly-scoped, NOT defects)
1. **Efficacy score is seed-type proportional only**, not the full TargetScan context++ model (no real 3′-supplementary, local-AU, UTR-position integration into the headline `Score`). Neither code nor spec advertises context++; `AnalyzeTargetContext`/`CalculateSiteAccessibility` exist as clearly-separate simplified helpers. Honestly scoped per spec Deviation #3.
2. **ΔG (`FreeEnergy`) for target duplexes uses a simplified positional model**, so short G-flanked test constructs can yield small positive ΔG. Classification depends on `Score`, not ΔG; M-015 only requires ΔG<0 for a perfectly complementary duplex, which holds. Honestly scoped per spec Deviation #1.

---

## Verdict & follow-ups
- **Stage A: PASS** — site definitions, A1 rule, seed-reverse-complement rule, and efficacy hierarchy all confirmed against Bartel (2009), Agarwal (2015)/TargetScan, Grimson (2007).
- **Stage B: PASS-WITH-NOTES** — implementation correctly uses seed **reverse complementarity** (not identity), correct A1/m8 rules and site-type classification, correct antiparallel geometry, Grimson-proportional scoring. The two notes are intentional, documented simplifications, not biological errors.
- **State: CLEAN** — no defect found; no code changes required. Full suite 4486/4486 passing; 27/27 TargetPrediction tests passing.
- Follow-ups: none required.
