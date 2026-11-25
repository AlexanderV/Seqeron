# Validation Report: MIRNA-CLASSIFY-001 — Pre-miRNA Structure-Feature Classifier

- **Validated:** 2026-06-25   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.ClassifyPreMiRna` (logistic regression over MFE / AMFE / MFEI / GC% / %paired);
  helpers `ExtractPreMiRnaFeatures`, `ScorePreMiRnaFeatures`, `CalculateMfeIndex`, `DinucleotideShuffle`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

> Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
> Tests: `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (region "Trained natural-vs-background classifier", CL1–CL14)

---

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

- **Bonnet, Wuyts, Rouzé & Van de Peer (2004)**, *Bioinformatics* 20(17):2911–2917, "Evidence that microRNA
  precursors, unlike other non-coding RNAs, have lower folding free energies than expected." Establishes the
  **di-nucleotide-shuffle null model**: a candidate's MFE is compared to MFEs of sequences with identical
  mono- and **di-nucleotide composition** (the shuffle preserves local composition so discrimination must come
  from **structure**, not base content). Real pre-miRNAs fold to significantly lower (more negative) MFE than
  their shuffled controls. This is exactly the negative-set convention used here.
- **Zhang, Pan, Cobb & Anderson (2006)**, *Cell Mol Life Sci* 63:246–254. Defines **AMFE = 100·|ΔG°|/length**
  (adjusted MFE, i.e. MFE per 100 nt) and **MFEI = AMFE / (G+C)%**, and reports that pre-miRNA MFEI is
  characteristically high (≥ ~0.85), separating miRNA precursors from tRNA/rRNA/mRNA. Both formulas in the code
  match the paper exactly.
- **miRBase** (public-domain hairpin precursors). The 13 positives + the two named precursors (hsa-mir-21
  MI0000077, hsa-let-7a-1 MI0000060) were re-fetched from `mirbase.org/hairpin/<MI>` this session and match the
  repo strings **byte-for-byte** (see cross-check). Public-domain license confirmed.
- **Altschul & Erickson (1985)** Eulerian-walk shuffle — the algorithm `DinucleotideShuffle` implements to
  produce di-nucleotide-preserving negatives (preserves first base, last base, length, exact doublet counts).
- **Hastie, Tibshirani & Friedman, ESL 2nd ed. §4.4** — logistic regression link `P = σ(b0 + Σ b_j z_j)` and the
  0.5 Bayes cutoff. The scorer is standard standardised-feature logistic regression.

### Formula check
- **MFE (ΔG°):** the discriminative thermodynamic signal; taken from the **validated RNA-STRUCT-001 Zuker–Stiegler
  Turner-2004 folder** (`RnaSecondaryStructure.CalculateMfeStructure`) — *not* a private heuristic. ✅
- **AMFE = 100·|ΔG°|/length** (`ExtractPreMiRnaFeatures`, `CalculateMfeIndex`). Matches Zhang 2006. ✅
- **MFEI = AMFE/(G+C)%** with G+C as a **percentage** (guarded for length≤0, GC%≤0 → 0). Matches Zhang 2006. ✅
- **GC%/GC fraction:** straight G+C count over length. ✅
- **%paired = 2·(base pairs)/length** (each pair occupies two bases). Standard base-pairing-propensity feature
  (microPred / Xue 2005 family). ✅
- **Logistic score:** `z = bias + Σ w_j·(raw_j − mean_j)/std_j` over feature order
  `[FreeEnergy, Amfe, Mfei, GcContent, PairedFraction]`, then `σ(z)`. Standard standardised logistic regression. ✅

### Edge-case semantics (sourced/defined)
- null/empty → `null` (documented failure mode).
- Non-complementary / structureless candidate → ΔG°=0 ⇒ AMFE=MFEI=%paired=0 ⇒ P≈0 ⇒ **background** (the classifier
  scores any non-empty sequence; it is not a hairpin gate, unlike `AssessHairpinByMfe`). Defined and tested (CL13).
- DNA input (T) and case are normalised (T→U, upper) by the folder/shuffle. Non-ACGU symbols do not throw.
- `DinucleotideShuffle`: length<2 returned unchanged; null RNG throws `ArgumentNullException`.

### Independent cross-check (exact numbers, hand-derived this session)

miRBase sequences re-fetched and verified identical to repo: hsa-mir-21 (72 nt) and hsa-let-7a-1 (80 nt; the
mirbase page render miscounted as "81" but the string is 80 nt and identical char-for-char).

Given the validated folder's outputs for **hsa-mir-21** (ΔG° = −35.13, base pairs = 32, length = 72, G+C = 35):

| Feature | Hand-derived | Code (runtime) | Test pin |
|---|---|---|---|
| GC fraction | 35/72 = 0.4861111111111111 | 0.48611111111111116 | ✅ (1e-9) |
| AMFE | 100·35.13/72 = 48.79166666666667 | 48.79166666666667 | ✅ |
| MFEI | 48.79166.../48.6111... = 1.0037142857142858 | 1.0037142857142858 | ✅ |
| %paired | 2·32/72 = 0.8888888888888888 | 0.8888888888888888 | ✅ |

**Logistic score, hand-computed from the bundled weights** for hsa-mir-21:
standardised contributions = [+0.9134, +4.6396, +4.7950, −0.0486, +6.2291], `z = 12.18772013945276`,
**P(natural) = 0.9999949074161931** — **identical to the code's runtime output** (verified via a temporary probe
test, since removed). Classified NATURAL. ✅

**hsa-let-7a-1** (ΔG° = −34.31, bp = 32, n = 80, G+C = 34): AMFE = 42.8875, MFEI = 1.0091176470588237,
%paired = 0.8 — match the MF3/CL pins. ✅

**Di-shuffle negative (independently regenerated):**
- seed 999 shuffle of hsa-mir-21 (CL5): same dinucleotide counts as the original (verified: identical doublet
  multiset, same first/last base, length 72). Its fold degrades to ΔG° = −15.78, MFEI = 0.4509, %paired = 0.667
  ⇒ **P(natural) = 0.000522, IsNatural = False** (below threshold). Discrimination rests on structure, not
  composition. ✅
- **Independent fresh seed 12345** shuffle (also composition-matched): **P(natural) = 0.0240, IsNatural = False** —
  the background call is not seed-specific. Locked by new test CL14. ✅

**Held-out separability (Bonnet-2004 expectation):** rebuilding the exact training dataset deterministically
(13 miRBase positives + 4 di-shuffles each, fixed seed 20060101, fixed 70/30 split → 20 held-out examples) and
scoring with the bundled model gives **held-out accuracy = 1.0 and Mann–Whitney AUC = 1.0** (CL9 asserts ≥ 0.90;
actual = 1.0), confirming the claimed AUC ≈ 1.0 on the small held-out set.

### Findings / divergences (Stage A)
None affecting correctness. The classifier deliberately uses the 5 thermodynamic/composition/pairing features
(StemBasePairs/LoopSize/Length are reported but excluded as collinear/length-confounded — documented in source).
The small (13-positive) public-domain training set and the **out-of-scope miRDeep2 read-stacking signal** are
acknowledged, documented boundaries (no GPL miRDeep2 code is used). **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed
- `ExtractPreMiRnaFeatures` (MiRnaAnalyzer.cs:2268) — folds once via the validated RNA-STRUCT-001 engine, then
  computes AMFE/MFEI/GC/%paired/loop exactly per the formulas above.
- `CalculateMfeIndex` (1832) and `GcPercent` (2007) — Zhang-2006 formulas with degenerate guards.
- `ScorePreMiRnaFeatures` (2327) + `Sigmoid` (2350) — standardised-feature logistic regression; zero-std guard.
- `ClassifyPreMiRna` (2310) — threshold-parameterised boolean call; probability is threshold-independent.
- `DinucleotideShuffle` (2385) + `LastEdgesReachRoot` / `FisherYatesShuffle` — Altschul–Erickson Eulerian-walk
  shuffle with acyclic-arborescence retry and a composition-preserving fallback.

### Formula realised correctly?
Yes — the runtime feature values and the logistic probability reproduce the hand-derivation **to full double
precision** (table above). MFE comes from the validated folder (cross-checked against `CalculateMinimumFreeEnergy`
in tests MF1–MF3). No approximation substituted.

### Cross-verification (recomputed vs code)
All values in the Stage-A table were produced by the actual code at runtime and match the external/hand-derived
references exactly. Held-out AUC/accuracy reproduced = 1.0.

### Variant/delegate consistency
`ClassifyPreMiRna` delegates feature extraction to `ExtractPreMiRnaFeatures` and scoring to
`ScorePreMiRnaFeatures`; CL10 pins that a real precursor scores strictly above its shuffle; CL11 pins
threshold independence of the probability. Consistent.

### Test quality audit
Expected values trace to **miRBase** (sequences verified this session), **Zhang 2006** (AMFE/MFEI hand-derived),
**hand-computation** (logistic z and P, %paired, GC), and the **validated RNA-STRUCT-001 folder** (ΔG°/bp,
cross-checked against the engine) — **not code echoes**. The suite covers: real precursor → positive (CL3/CL4),
di-shuffle → negative (CL5), feature math pinned exactly (CL1), held-out AUC (CL9), monotonicity (CL10), threshold
(CL11), shuffle invariants (CL6–CL8), null/empty (CL2/CL12).

**Hardening added this session (test-quality gate):**
- **CL13** — `ClassifyPreMiRna_NonComplementary_IsBackground`: structureless all-purine input (ΔG°=0, %paired=0)
  is classified background (P<0.05, IsNatural=False). Covers the Stage-A non-hairpin edge case for the *classifier*
  (previously only `AssessHairpinByMfe` had a non-complementary test).
- **CL14** — `ClassifyPreMiRna_IndependentDiShuffle_IsBackground`: an **independent fresh seed** di-shuffle is also
  background while preserving dinucleotide composition — proves the negative call is structural, not seed-cherry-picked.

### Findings / defects (Stage B)
None. **Stage B: PASS.**

---

## Verdict & follow-ups

**✅ CLEAN.** Stage A PASS, Stage B PASS. Feature math (AMFE/MFEI/GC/%paired) and the logistic score reproduce
hand-derivation to full precision; MFE traces to the validated RNA-STRUCT-001 folder; miRBase positives verified
against the live database; di-shuffle negatives (two independent seeds) score below threshold; held-out AUC = 1.0.
Two hardening tests (CL13, CL14) added and locked. Full unfiltered `dotnet test Seqeron.sln -c Debug` →
Failed: 0 (Genomics 18771 passed), 0 warnings. No defects logged.

**Documented boundaries (acceptable):** small public-domain training set (13 positives); miRDeep2 read-stacking
signal is out of scope by design (no GPL code used).

---

## B13 review (2026-11-25) — Stage A PASS-WITH-NOTES, Stage B FAIL → FIXED (F23–F26)

Sources opened: `randfold.c` and squid 1.9g `shuffle.c` (`StrDPShuffle`) from github.com/eb00/randfold (master);
uShuffle 1.1.2 source (`ushuffle.c`, PyPI sdist), compiled with a replaying `randfunc`; secondary descriptions of
miPred dP (Ng & Mishra 2007; Schultes 1999: dP = |bp|/L ∈ [0, 0.5]); ESL §4.4.1 (Newton–Raphson/IRLS).
Oracles: scikit-learn 1.9.1 / scipy 1.18.1 (fit), ViennaRNA 2.8.0 `RNA.fold` (randfold MFEs), uShuffle (shuffle).

- **F23 `DinucleotideShuffle`:** the Altschul–Erickson rejection step was capped at 1000 draws and then returned the
  INPUT unshuffled. Repro `A×500 G A×500 C` (999 valid shuffles): 69/200 seeds returned the input. Now Wilson's
  arborescence as in uShuffle: 2000/2000 scripted-draw cases identical to `ushuffle.c`; χ² uniform over the 612 / 540
  valid sequences of two test strings (p 0.66 / 0.62); repro 5/2000 (expected 2).
- **F24 model provenance:** the trainer was a never-committed scratch file; its 20 000-epoch gradient ascent was an
  unconverged approximation of the penalised MLE, and B12's folder change moved 61/65 training MFEs. The old
  constants were re-derived bit-exactly (numpy, pre-B12 tree) to prove the reconstruction; the trainer is now in the
  repo (Newton–Raphson), the constants were regenerated and equal scikit-learn to 4e-12; a test re-derives them.
  Held-out: accuracy 0.90, AUC 0.984375 (was 1.0/1.0 on the old shuffles and folder).
- **F25 feature definitions:** PairedFraction = 2·|bp|/L is twice the published dP; dP added as
  `PreMiRnaFeatures.BasePairingPropensity`; docs no longer attribute raw MFE to Bonnet 2004 nor %paired to Xue 2005;
  threshold doc corrected (training prior 1 : 4, not balanced).
- **F26 Bonnet 2004 randfold test** was missing: `CalculateRandfoldPValue` (randfold.c counting, ε 1e-7) added.
