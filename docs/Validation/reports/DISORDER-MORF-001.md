# Validation Report: DISORDER-MORF-001 — MoRF (Molecular Recognition Feature) Prediction

- **Validated:** 2026-06-16   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.PredictMoRFs(string sequence, int minLength = 10, int maxLength = 70)` (internal helper `PredictDisorder` supplies per-residue scores)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened & what they confirm (all retrieved live this session)

| Source | URL | Confirms |
|--------|-----|----------|
| Mohan et al. (2006) J Mol Biol, PubMed record | https://pubmed.ncbi.nlm.nih.gov/16935303/ | MoRFs are "relatively short (10-70 residues), loosely structured protein regions within longer, largely disordered sequences"; "upon binding to their partner(s), MoRFs undergo disorder-to-order transitions"; "at least three basic types … alpha-MoRFs, beta-MoRFs, and iota-MoRFs, which form alpha-helices, beta-strands, and irregular secondary structure when bound, respectively". |
| Cheng/Oldfield et al., Biochemistry | https://pmc.ncbi.nlm.nih.gov/articles/PMC2570644/ | Heuristic "identifies short regions of order within longer regions of disorder – or 'dips'" in disorder prediction profiles; α-MoRF is a "short (around 20 residues) structural element"; candidates "30 residues or less"; order/disorder boundary "the threshold of 0.5". |
| Wikipedia, Molecular recognition feature | https://en.wikipedia.org/wiki/Molecular_recognition_feature | "small (10-70 residues) intrinsically disordered regions … that undergo a disorder-to-order transition upon binding"; "disordered prior to binding … form a common 3D structure after interacting". (Lists subtypes as α/β/irregular/complex.) |
| Campen et al. (2008) TOP-IDP, PMC | https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ | Table 2 per-residue propensities; confirmed W −0.884, I −0.486, L −0.326, E 0.736, P 0.987; scale min/max = −0.884 / 0.987 (range 1.871). |
| Oldfield et al. (2005), PubMed | https://pubmed.ncbi.nlm.nih.gov/16156658/ | MoRE = "short region that undergoes coupled binding and folding within a longer region of disorder"; exact numeric dip flank/run parameters are in the paywalled Methods (not retrievable). |

### Formula check

The description models a MoRF as a maximal interval `[s,e]` with (1) `d(i) < 0.5` for all `i` (ordered "dip"), (2) `10 ≤ e−s+1 ≤ 70`, (3) `d(s−1) ≥ 0.5` and `d(e+1) ≥ 0.5` (flanked by disorder). Score = `(0.5 − mean d)/0.5` clamped to `[0,1]`. Every constant traces to a retrieved source: threshold 0.5 (PMC2570644), length band 10–70 (Mohan 2006 / Wikipedia), per-residue scores = normalized TOP-IDP (Campen 2008). The score normalization is an honest derivation from the 0.5 threshold (the maximum possible dip depth), not a tuned constant — confirmed reasonable.

### Edge-case semantics check

Fully ordered → no flanking disorder → no MoRF; fully disordered → no dip → no MoRF; dip outside 10–70 → excluded; dip at terminus → not flanked both sides → excluded; null/empty → empty. All four corner cases are sourced (PMC2570644, Mohan 2006, Oldfield 2005) and defined, not "implementation-defined".

### Independent cross-check (numbers)

Re-derived the smoothed disorder profile *from the source TOP-IDP raw values* (not from the repo) with a standalone Python reimplementation of the window-21 mean: normalized P = 1.000, L = 0.2983, I = 0.2128. For the 25P+30L+25P construct the smoothed profile dips below 0.5 over residues **[29,50]** (length 22), mean disorder **0.362033**, score **0.275934**. For 25P+30I+25P the dip is [28,51] (length 24), mean 0.300196, score **0.399608**. These independently reproduce the test's locked values.

### Findings / divergences (Stage A)

- **N1 (PASS-WITH-NOTES):** The exact flank/run-length dip parameters live in Oldfield 2005's paywalled Methods and could not be retrieved; the unit uses a documented qualitative approximation (ordered run < 0.5, flanked by ≥1 disordered residue, 10–70 band). This is honestly recorded in the Evidence Assumption Register and the algorithm doc §5.3/§5.4. The load-bearing constants (0.5, 10–70, TOP-IDP) are all source-traceable, so this is a bounded modeling note, not a correctness error.
- **N2 (doc nit, fixed):** TestSpec §4.1/§5.6 M1 row stated coordinates "20–34" describing a 15L/20P construct that does not match the actual test (25P+30L+25P → 29–50). Corrected the spec prose to 29–50 / length 22 / score 0.275934 to match the sourced derivation and the test. No code/test behaviour change.
- **N3 (terminology, no defect):** Wikipedia lists subtypes as α/β/irregular/complex; the spec/Evidence use Mohan 2006's "ι (iota)" naming. Mohan 2006 (retrieved) uses "iota-MoRFs" verbatim, so the spec is correct; Wikipedia is merely less precise. No change.

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs:615-671` (`PredictMoRFs`), using `PredictDisorder` (`:190`) → `CalculatePerResidueScores` (window 21, `:227`) → `CalculateDisorderScore` (normalized TOP-IDP mean, `:255`). Constants `MoRFOrderThreshold = 0.5` (`:578`), `MoRFMinLength = 10` (`:584`), `MoRFMaxLength = 70` (`:589`).

### Formula realised correctly?

Yes. The scan finds maximal runs with `DisorderScore < 0.5`, applies the `[minLength, maxLength]` filter, requires `>= 0.5` flanks on both immediate sides (rejecting terminal dips, since `start > 0` / `end < count-1` are required), and emits `(0.5 − meanDisorder)/0.5` clamped to `[0,1]`. This matches the validated description exactly.

### Cross-verification table recomputed vs code (independent Python vs test expectations)

| Case | Construct | Independent recompute | Test expectation | Match |
|------|-----------|----------------------|------------------|-------|
| M1 | 25P+30L+25P | (29,50) score 0.2759341 | (29,50) 0.275934 | ✅ |
| M7 (I) | 25P+30I+25P | (28,51) score 0.3996081 | score 0.399608 | ✅ |
| M2 | 40L | ∅ | ∅ | ✅ |
| M3 | 40P | ∅ | ∅ | ✅ |
| M4 | 25P+16L+25P | ∅ (dip len 8 < 10) | ∅ | ✅ |
| M5 | 25P+95L+25P | ∅ (dip len 87 > 70) | ∅ | ✅ |
| M6 | 15L+30P | ∅ (terminal) | ∅ | ✅ |
| S1 | two L runs | (29,50)+(89,110) | (29,50)+(89,110) | ✅ |

All independently-derived values (from source TOP-IDP, not from the C# code) match the locked test expectations.

### Variant/delegate consistency

The default-parameter overload (M1/M2/…) and the custom-bound calls (S3 min/max) are the same method with different arguments; both exercised and consistent. `PredictDisorder` (the score source) is validated under DISORDER-PRED-001.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M1/M7/S1 lock exact coordinates and scores independently derivable from Campen 2008's TOP-IDP values (a wrong dip detector or wrong score normalization fails the `.Within(1e-6)` assertions). Verified by my standalone Python reimplementation, not by the repo's output.
- **No green-washing:** exact equality on all known values (coordinates with `EqualTo`, scores with `Within(1e-6)`); the only `GreaterThan`/`InRange` assertions are the genuine monotonicity-and-bounds invariants in M7 and the INV property test, which sit alongside the exact-value assertions — not in place of them. No skipped/ignored tests; no widened tolerances.
- **Coverage:** the single public method and its custom-bounds overload are exercised; all Stage-A branches are covered — happy dip (M1), no-flank ordered (M2), no-dip disordered (M3), under-length (M4), over-length (M5), terminal (M6), score monotonicity+bounds (M7), multi-dip independence (S1), case-insensitivity (S2), custom min/max bounds (S3), null/empty (C1), too-short (C2), and INV-1..INV-5 property test.
- **Honest green:** full unfiltered suite **Failed: 0, Passed: 6609**; `dotnet build` 0 errors. The 4 build warnings are pre-existing NUnit-analyzer notices in unrelated test files (ApproximateMatcher, etc.), none in this unit's files.

**Gate result: PASS.** No weak/code-echoing/green-washed/partial tests in this unit.

### Findings / defects (Stage B)

None. One spec-prose nit (N2 above) was corrected; it is documentation only and did not affect code or tests.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (sourced qualitative-flank approximation N1; doc nit N2 fixed; terminology N3 confirmed correct).
- **Stage B:** PASS — code faithfully realises the validated dip-in-disorder description; all independently-derived values match; tests are exact and sourced.
- **End-state:** ✅ CLEAN. No algorithm defect. Algorithm is fully functional.
- **Follow-ups:** none required. The flank-length detail (N1) remains a documented bounded assumption pending retrieval of Oldfield 2005's paywalled Methods; it does not affect the sourced constants and is out of scope for a non-trained heuristic annotator.
