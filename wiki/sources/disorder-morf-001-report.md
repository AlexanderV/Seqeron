---
type: source
title: "Validation report: DISORDER-MORF-001 (MoRF prediction — dip-in-disorder, DisorderPredictor.PredictMoRFs)"
tags: [validation, analysis]
doc_path: docs/Validation/reports/DISORDER-MORF-001.md
sources:
  - docs/Validation/reports/DISORDER-MORF-001.md
source_commit: dc13c70fe90b2fa15b75169c79a58a5cd060d39a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: DISORDER-MORF-001

The two-stage **validation write-up** for test unit **DISORDER-MORF-001** — **MoRF (Molecular
Recognition Feature) prediction**, the "dip within disorder" heuristic that flags a short ordered
segment embedded in a longer intrinsically disordered region, validated 2026-06-16. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code (Stage B),
and the wider campaign is [[validation-and-testing]]. The algorithm, its criterion, oracles and
sub-types are synthesized on the concept [[morf-prediction-dip-in-disorder]]; [[test-unit-registry]]
defines the unit. Distinct from [[disorder-morf-001-evidence]] — the pre-implementation evidence
artifact sourced from `docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: ✅ CLEAN.** No algorithm defect; no code
change required. Full unfiltered suite **6609 passed, 0 failed**; `dotnet build` 0 errors (the 4
build warnings are pre-existing NUnit-analyzer notices in unrelated test files, none in this unit).
Only one spec-prose nit (N2) was corrected — documentation only, no code/test behaviour change.

## Canonical method & source under test

- `DisorderPredictor.PredictMoRFs(string sequence, int minLength = 10, int maxLength = 70)` in
  `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs:615-671`, using the internal
  `PredictDisorder` (`:190`) → `CalculatePerResidueScores` (window 21, `:227`) →
  `CalculateDisorderScore` (normalized TOP-IDP mean, `:255`). Constants `MoRFOrderThreshold = 0.5`
  (`:578`), `MoRFMinLength = 10` (`:584`), `MoRFMaxLength = 70` (`:589`). `PredictDisorder` (the score
  source) is itself validated under DISORDER-PRED-001 ([[intrinsic-disorder-prediction-top-idp]]).

## Stage A — description (algorithm faithfulness)

Confirmed live against **Mohan et al. 2006** (PMID 16935303 — MoRFs are "relatively short (10–70
residues), loosely structured protein regions within longer, largely disordered sequences" that
"undergo disorder-to-order transitions" on binding; α/β/ι subtypes), **Cheng/Oldfield** (PMC2570644 —
the operational "dip" definition and the **0.5** order/disorder threshold; α-MoRF ≈ 20 residues,
candidates ≤ 30), **Oldfield et al. 2005** (PMID 16156658 — MoRE = coupled folding-and-binding within
a longer disordered region), **Campen et al. 2008 TOP-IDP** (PMC2676888 — Table 2 per-residue
propensities: W −0.884, I −0.486, L −0.326, E 0.736, P 0.987; min/max −0.884/0.987, range 1.871), and
**Wikipedia** "Molecular recognition feature".

The description models a MoRF as a maximal interval `[s,e]` with (1) `d(i) < 0.5` for all `i` (ordered
"dip"), (2) `10 ≤ e−s+1 ≤ 70`, (3) `d(s−1) ≥ 0.5` and `d(e+1) ≥ 0.5` (flanked both sides by
disorder). Score = `(0.5 − mean d)/0.5` clamped to `[0,1]`. Every constant traces to a retrieved
source: threshold 0.5 (PMC2570644), length band 10–70 (Mohan 2006 / Wikipedia), per-residue scores =
normalized TOP-IDP (Campen 2008). The score normalization is an honest derivation from the 0.5
threshold (the maximum possible dip depth), not a tuned constant. All four corner cases (fully ordered
→ no flank → no MoRF; fully disordered → no dip → no MoRF; dip outside 10–70 → excluded; terminal dip
→ not flanked both sides → excluded; null/empty → empty) are sourced and defined, not
implementation-defined.

**Independent numeric cross-check.** The validator re-derived the smoothed disorder profile *from the
source TOP-IDP raw values* (standalone Python reimplementation of the window-21 mean, not the repo):
normalized P = 1.000, L = 0.2983, I = 0.2128. For `25P+30L+25P` the profile dips below 0.5 over
residues **[29,50]** (length 22), mean disorder 0.362033, score **0.275934**; for `25P+30I+25P` the
dip is **[28,51]** (length 24), mean 0.300196, score **0.399608**. These independently reproduce the
test's locked values.

### Stage A notes (documented, non-blocking)

- **N1 (PASS-WITH-NOTES):** the exact flank/run-length dip parameters live in **Oldfield 2005's
  paywalled Methods** and could not be retrieved; the unit uses a documented qualitative approximation
  (ordered run < 0.5, flanked by ≥1 disordered residue, 10–70 band). Load-bearing constants (0.5,
  10–70, TOP-IDP) are all source-traceable, so this is a **bounded modeling note, not a correctness
  error** — the same single assumption recorded on [[disorder-morf-001-evidence]].
- **N2 (doc nit, fixed):** TestSpec §4.1/§5.6 M1 row stated coordinates "20–34" for a 15L/20P
  construct that does not match the actual test (`25P+30L+25P` → 29–50). The spec prose was corrected
  to 29–50 / length 22 / score 0.275934. No code/test behaviour change.
- **N3 (terminology, no defect):** Wikipedia lists subtypes as α/β/irregular/complex; the spec uses
  Mohan 2006's "ι (iota)" naming. Mohan 2006 uses "iota-MoRFs" verbatim, so the spec is correct;
  Wikipedia is merely less precise. No change.

## Stage B — implementation

The scan finds maximal runs with `DisorderScore < 0.5`, applies the `[minLength, maxLength]` filter,
requires `≥ 0.5` flanks on both immediate sides (rejecting terminal dips, since `start > 0` /
`end < count−1` are required), and emits `(0.5 − meanDisorder)/0.5` clamped to `[0,1]` — matching the
validated description exactly. Every worked case was recomputed from source TOP-IDP (independent
Python, not the C# output) and matched the locked test expectations:

| Case | Construct | Independent recompute | Test expectation | Match |
|------|-----------|-----------------------|------------------|-------|
| M1 | 25P+30L+25P | (29,50) score 0.2759341 | (29,50) 0.275934 | ✅ |
| M7 (I) | 25P+30I+25P | (28,51) score 0.3996081 | 0.399608 | ✅ |
| M2 | 40L | ∅ | ∅ | ✅ |
| M3 | 40P | ∅ | ∅ | ✅ |
| M4 | 25P+16L+25P | ∅ (dip len 8 < 10) | ∅ | ✅ |
| M5 | 25P+95L+25P | ∅ (dip len 87 > 70) | ∅ | ✅ |
| M6 | 15L+30P | ∅ (terminal) | ∅ | ✅ |
| S1 | two L runs | (29,50)+(89,110) | (29,50)+(89,110) | ✅ |

The default-parameter overload and the custom-bound calls (S3 min/max) are the same method with
different arguments; both exercised and consistent.

**Test-quality audit (HARD gate) — PASS.** M1/M7/S1 lock exact coordinates and scores independently
derivable from Campen 2008's TOP-IDP values (a wrong dip detector or wrong score normalization fails
the `.Within(1e-6)` assertions). Exact equality on all known values (coordinates `EqualTo`, scores
`Within(1e-6)`); the only `GreaterThan`/`InRange` assertions are genuine monotonicity-and-bounds
invariants (M7, INV property test) sitting **alongside** exact-value assertions, not in place of them.
No skipped/ignored tests, no widened tolerances. Coverage spans happy dip (M1), no-flank ordered (M2),
no-dip disordered (M3), under-length (M4), over-length (M5), terminal (M6), score
monotonicity+bounds (M7), multi-dip independence (S1), case-insensitivity (S2), custom bounds (S3),
null/empty (C1), too-short (C2), and INV-1..INV-5.

## Findings & follow-ups

- **No code defect (State CLEAN).** Code faithfully realises the validated dip-in-disorder
  description; all independently-derived values match; tests are exact and sourced.
- **Documentation only:** the M1 coordinate nit (N2) was corrected in the spec prose.
- **Standing bounded assumption (N1, non-blocking):** the exact flank-length detail remains pending
  retrieval of Oldfield 2005's paywalled Methods; it does not affect the sourced constants and is out
  of scope for a non-trained heuristic annotator. No follow-ups required.
</content>
</invoke>
