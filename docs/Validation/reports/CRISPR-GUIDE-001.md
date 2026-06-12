# Validation Report: CRISPR-GUIDE-001 — Guide RNA (sgRNA) Design & On-Target Quality Scoring

- **Validated:** 2026-06-12   **Area:** MolTools (CRISPR)
- **Canonical method(s):**
  - `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)`
  - `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_GuideRNA_Tests.cs`, `Properties/CrisprProperties.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect; the heuristic nature of the score is already correctly declared)

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Guide RNA"** (https://en.wikipedia.org/wiki/Guide_RNA):
  - "The length of guide sequences is typically 20 bp, but they can also range from 17 to 24 bp." → confirms canonical **20 nt** spacer for SpCas9 (range 17–24).
  - "The optimal GC content of the guide sequence should be **over 50%**. A higher GC content enhances the stability of the RNA-DNA duplex and reduces off-target hybridization."
  - sgRNA = crRNA + tracrRNA fused by a tetraloop; tracrRNA contributes the stem-loop scaffold (spacer + scaffold structure confirmed).
- **Addgene CRISPR Guide** (https://www.addgene.org/guides/crispr/):
  - gRNA = scaffold (Cas-binding) + **~20 nt spacer** (genomic target).
  - "the **seed sequence (8–10 bases at the 3′ end** of the gRNA targeting sequence) will begin to anneal"; 3′ seed mismatches inhibit cleavage. → confirms seed = PAM-proximal 8–10 nt at 3′ for Cas9.
  - PAM table: SpCas9 = **3′ NGG**, SaCas9 = **3′ NNGRRT / NNGRR(N)**, Cas12a = **5′ TTTV**. All match the implementation's `GetSystem` table.
- **Pol III termination & GC literature** (Synthego/abm/Schindele et al. 2020; Gao et al.):
  - RNA Pol III terminates at a poly-T stretch; **T4 (TTTT) is the minimal terminator**, full efficiency at ≥6 T. → TTTT detection is the correct, conservative trigger.
  - GC content recommendation commonly cited as **40–80%** (with 40–60% "optimal"); "strongest predictor of efficiency is high GC content in the first ten bases." → the implementation's 40–70% window sits validly inside the accepted 40–80% band.

### Defined design rules (verified)
| Rule | Spec/impl value | Source | Verdict |
|------|-----------------|--------|---------|
| Spacer length (SpCas9) | 20 nt (SaCas9 21, Cas12a 23) | Wikipedia/Addgene | ✅ canonical |
| PAM adjacency (NGG SpCas9) | links to CRISPR-PAM-001 via `FindPamSites` | Addgene/Wikipedia | ✅ |
| GC window | [40, 70] default; penalty `(delta)×2` per % | 40–80% accepted band | 🟡 valid subset of accepted range (not >50% per Wikipedia, but within community 40–80%) |
| Seed region | last 10 nt (3′) for Cas9, first 10 for Cas12a | Addgene 8–10 nt | ✅ upper-bound choice, documented |
| Poly-T terminator | TTTT (≥4 T) → −20 | T4 minimal terminator | ✅ conservative |
| Scaffold | fixed 76 nt SpCas9 scaffold appended | tracrRNA-derived | ✅ |

### Is the score a published predictor or a heuristic?
The score is a **simple composition heuristic** (base 100 minus GC/poly-T/self-complementarity/seed-GC/restriction-site deductions). It does **not** claim to be Doench Rule Set 2 / Azimuth, and the algorithm doc (`docs/algorithms/MolTools/Guide_RNA_Design.md` §1, §2.2, §5.3) explicitly and repeatedly states the score is "a heuristic quality measure rather than a learned or experimentally calibrated predictor" and that it does "not reproduce models such as Doench-style activity prediction." **This is the correct, non-overclaiming declaration the protocol asks for** — so it is NOT a LIMITED case.

### Worked examples (hand-checked, GC window + poly-T)
- **M-001** `ACGTACGTACGTACGTACGT`: 10 GC/20 = 50% GC (40≤50≤70, no penalty); seed last-10 `ACGTACGTAC` = 50% (in [30,80]); no TTTT. Score = 100. ✓
- **M-004** `ACGTACGTTTTTACGTACGT`: 8 GC/20 = 40% (=Min, no penalty); contains `TTTTT` ⊇ `TTTT` → −20; seed last-10 `TTACGTACGT`... = 40% (in range). Score = 100−20 = **80**. ✓
- **S-009** `ACGTACGTACGTACGTTTAC`: only 3 consecutive T → no poly-T; 45% GC. Score = 100. ✓ (confirms TTTT, not TTT, is the trigger.)
- **C-001** `GCGCGCGC` (8 nt): 100% GC → (100−70)×2 = −60; selfComp 0.3125 → ×30 = −9.375; seed (full 8 nt) 100% → −5. Score = 100−60−9.375−5 = **25.625**. ✓

### Findings / divergences
- **N-1 (minor):** GC window is 40–70%, whereas Wikipedia cites ">50% optimal" and broader literature cites 40–80%. The chosen window is a defensible subset of the accepted range and is documented; not a defect. (PASS-WITH-NOTES.)
- **N-2 (informational, not in scope of scoring correctness):** `AvoidPolyT` and `CheckSelfComplementarity` flags in `GuideRnaParameters` are stored but not honored by the scoring path. This is **already documented** in the algorithm doc (§3.1, §4.2, §5.2, §5.3) and is not asserted otherwise by any test, so it is a declared simplification rather than a hidden defect.

---

## Stage B — Implementation

### Code path reviewed
`CrisprDesigner.cs`:
- `EvaluateGuideRna(string,...)` (lines 200–284): GC content (212–242), poly-T `HasPolyT(seq,4)` (216, 449–452), self-complementarity (219, 454–471), seed GC via PAM-orientation 10 nt window (222–226), restriction-site check (266–271), `Math.Max(0, score)` floor (281).
- `DesignGuideRnas` (169–195): null/range guards (176–180), PAM extraction via `FindPamSitesCore`, region filter `IsInRegion` (309–317), per-candidate evaluation + `MinScore` filter (192).
- `GuideRnaParameters.Default` (537–542): MinGC 40, MaxGC 70, MinScore 50.
- `GuideRnaCandidate.FullGuideRna` (560–563): spacer + fixed 76 nt scaffold.

### Formula realised correctly?
Yes. Deductions match the validated table exactly:
- GC: `(MinGc−gc)×2` / `(gc−MaxGc)×2` (lines 235, 240).
- Poly-T: flat −20 when `seq.Contains("TTTT")` (247).
- Self-comp: `selfComp×30` only when `> 0.3` (252–254).
- Seed GC: flat −5 when `<30 || >80` (259–261).
- Restriction site: flat −5 (266–270).
- Floor: `Math.Max(0, score)` (281).

### Cross-verification table recomputed vs code (via tests)
All exact-value expectations (M-001…C-010, 30 tests + 2 property tests) reproduced by the suite. Spot-checked worked examples above match the code output (M-001=100, M-004=80, S-009=100, C-001=25.625, C-002 clamped to 0). PAM-adjacency tests (S-003 pos 24, C-008 pos 4, C-005 count 3, C-004 empty) confirm guides only emerge adjacent to NGG. Edge cases covered: no valid guide → empty (C-004); GC at window edges 40%/70% not penalized (S-006/S-007); below/above 30%/80% (C-007/C-009); poly-T present (M-004/S-008) vs 3-T negative (S-009); empty/null guards (M-005/M-007/M-008/M-009/C-003).

### Variant/delegate consistency
`DesignGuideRnas` internally calls `EvaluateGuideRna(PamSite,...)` → same scoring path as the public string overload; only `Position`/`IsForwardStrand` are overwritten. Consistent.

### Test quality audit
Tests assert exact sourced/derived numeric values (scores, GC%, seed GC%, issue counts and contents), are deterministic, and cover every Stage-A edge case. Property tests (GC∈[0,100], Score∈[0,100]) are non-vacuous. No tautological "no-throw-only" assertions on the scoring path.

### Findings / defects
None. Code faithfully realises the validated (heuristic) description. The two simplifications (N-2 passive flags; SpCas9 metadata reuse for unmapped systems in designed candidates) are pre-documented in the algorithm doc and do not contradict any source or test.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — all canonical design rules (20 nt spacer, NGG/NNGRRT/TTTV PAM adjacency, seed 8–10 nt at 3′, TTTT Pol III terminator, scaffold) confirmed against Wikipedia/Addgene/literature; GC window 40–70% is a documented subset of the accepted 40–80% range.
- **Stage B: PASS** — implementation matches the formula; worked examples recompute exactly; edge cases handled and tested.
- **End-state: ✅ CLEAN.** The on-target "score" is a sequence-composition **heuristic**, and it is **correctly declared as such** (not advertised as Doench/Azimuth), so the overclaim risk the protocol flags does not apply here. Build green; 47 Guide tests pass; full suite 4461/4461 pass. No code or test changes required.
