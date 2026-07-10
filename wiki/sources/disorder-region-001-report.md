---
type: source
title: "Validation report: DISORDER-REGION-001 (disordered-region detection — contiguous-run IDR grouping + MobiDB-lite flavour typing)"
tags: [validation, analysis]
doc_path: docs/Validation/reports/DISORDER-REGION-001.md
sources:
  - docs/Validation/reports/DISORDER-REGION-001.md
source_commit: 8bd6a5e1b48367a44632668b3e6d980ebe3e6d2c
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: DISORDER-REGION-001

The two-stage **validation write-up** for test unit **DISORDER-REGION-001** — the **segment-calling
layer** that collapses the per-residue [[intrinsic-disorder-prediction-top-idp|`PredictDisorder`]]
TOP-IDP disorder profile into **contiguous IDR regions** (start/end/length, minimum-length filter,
`RegionType` label, `Confidence`), plus the opt-in **MobiDB-lite v3 flavour typing**
`ClassifyRegionFlavorMobiDbLite`. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B). The wider campaign is
[[validation-and-testing]]; [[test-unit-registry]] defines the unit. The region-calling rule, the two
classification schemes and their thresholds are synthesized on the concept
[[intrinsic-disorder-prediction-top-idp]] (the shared TOP-IDP anchor this layer aggregates).

Distinct from [[disorder-region-001-evidence]] — the pre-implementation evidence artifact sourced from
`docs/Evidence/` that records the source trace, the 12-row flavour table and the recommended coverage;
**this** page is the independent two-stage re-validation verdict, including a **fidelity defect found
and fully fixed this session**. Sibling reports [[disorder-pred-001-report]] (the windowed predictor),
[[disorder-propensity-001-report]] (the O(1) primitives), [[disorder-lc-001-report]] (SEG
low-complexity) and [[disorder-morf-001-report]] (MoRF) cover different units of the same
protein-disorder family.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN.** Fresh re-validation **2026-06-25**,
supersedes the 2026-06-24 pass — reset to pending because the campaign added the opt-in
`ClassifyRegionFlavorMobiDbLite` (F4). The **region boundary / grouping logic has no defect and needed
no change**; one **fidelity defect in the flavour typing** (histidine omitted from the positive-charge
fraction) was found and **completely fixed in-session** with sourced test coverage. Full unfiltered
`dotnet test Seqeron.sln -c Debug` → **Failed 0** (Genomics 18819 passed, incl. the 2 new His tests);
changed source project builds 0 warnings / 0 errors.

## Canonical methods & source under test

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`:

- Public `PredictDisorder(sequence, windowSize, threshold, minRegionLength)` → **`IdentifyDisorderedRegions`**
  (private, `:358-417`) — the single-pass region scan.
- **`ClassifyDisorderedRegion`** (`:428-462`) / **`CalculateConfidence`** (`:470-475`) — the default
  `RegionType` label and the rescaled per-region confidence.
- Opt-in public **`ClassifyRegionFlavorMobiDbLite(regionSequence)`** — MobiDB-lite 3.0 flavour typing
  (added this campaign; the F4 trigger for re-validation).

## Stage A — description (algorithm faithfulness)

External first sources retrieved live this session: **Ward et al. 2004 (DISOPRED2**, J Mol Biol
337:635-645) confirming long disordered regions = **>30 consecutive residues** (prevalence 2.0%
archaea / 4.2% bacteria / 33% eukaryotes) — the anchor for the code's `sequence.Length > 30` → **"Long
IDR"** label; **Campen et al. 2008 (TOP-IDP, PMC2676888)** for cutoff 0.542 / window 21 (scoring detail
belongs to [[disorder-pred-001-report|DISORDER-PRED-001]], reused here); **van der Lee et al. 2014**
for order↔disorder-transition boundaries and the compositional subtypes; **Das & Pappu 2013 (PNAS,
PMID 23901099)** for FCR/NCPR (FCR > 0.35 strong/weak boundary); and the **MobiDB-lite v3 reference
source** (`BioComputingUP/MobiDB-lite@v3`, fetched verbatim via `gh api`).

Region-calling logic validated as **correct and sourced**:

1. **Consecutive-grouping** — maximal runs of residues with `score ≥ threshold`. Standard.
2. **Threshold 0.542** (Campen 2008) — confirmed.
3. **Minimum region length** — a small **configurable** minimum (impl default 5); external practice
   supports both a short-IDR minimum and the ≥30 long-IDR convention. Within accepted practice.
4. **Coordinate convention** — 0-based, **inclusive** Start and End; `length = End − Start + 1`.
   Internally consistent.

Edge-case semantics all defined & sourced: no disordered residue → 0 regions; all-disordered &
len ≥ minLen → single region `[0, L−1]`; run < minLen excluded; run == minLen **included** (`>=`);
**trailing run to the last residue captured** by an explicit end-of-loop branch (the documented
off-by-one pitfall, spec §1.4.1).

The MobiDB-lite v3 constants were verbatim-matched to source: translation table
`intab='RKDEACFGHILMNPQSTVWY'` / `outab='PPNN____P___________'` ⇒ positive = **{R,K,H}**, negative =
**{D,E}**; charge gate `fcr > 0.35`; enrichment `is_enriched(threshold=0.32)` inclusive `>=`; priority
**charge (PA/PPE/NPE) → C → P → G → low-complexity → Polar {S,T,N,Q}**.

**Notes (PASS-WITH-NOTES).** The **classification heuristics** are disclosed in the spec (§6/§7) as
internal design decisions with **no published source**: default enrichment threshold **0.25**, priority
order **Pro > Acidic > Basic > S/T > Long > Standard**, confidence formula **`(mean − 0.542)/(1 − 0.542)`**,
AA groups `{E,D}`/`{K,R}`/`{S,T}`. These affect only the `RegionType` / `Confidence` **labels, not the
region boundaries or lengths**, and are honestly flagged. A second, cosmetic note: the code's Campen
ranking **comment** writes "…Q,S,K,E,P" whereas Wikipedia lists "…Q,K,S,E,P" (S/K swapped) — a
comment-only annotation; the numeric values (Q 0.318 < S 0.341 < K 0.586 < E 0.736 < P 0.987) are
self-consistent and move no boundary or flavour call.

### Independent cross-check (hand-recomputed this session, Python)

Window scoring + grouping re-implemented independently (window 21, halfWindow 10, threshold 0.542,
minLen 5), reproducing the spec coordinates **exactly**:

| Case | Sequence | Region(s) |
|------|----------|-----------|
| M2 full | P×30 | [0, 29] |
| M6 trail | W10+P20 | [11, 29] |
| S2 lead | P20+W30 | [0, 18] |
| S5 central | W15+P20+W15 | [16, 33] |
| M14 multi | W15+P20+W15+P20+W15 | [(16,33),(51,68)] |

Homopolymer MeanScore / Confidence also reproduced: P30 → 1.0 / 1.0; E30 → 0.866 / 0.707; K30 → 0.786 /
0.532; S30 → 0.655 / 0.246 — all matching the test assertions. FCR/NCPR flavour hand-trace (7 rows
incl. the two His cases) confirmed against code.

## Stage B — implementation (one defect, fixed)

The region scan realises the validated description faithfully:

- **Open region** on first disordered residue (`regionStart = i`, resets `scoreSum`, `:372-377`);
  **close on order residue** with `End = i − 1` (`:381, :389`) — correct inclusive boundary.
- **Trailing region** handled by an explicit post-loop block (`:400-416`, `End = count − 1`) — **no
  off-by-one**, the documented pitfall is covered.
- **Min-length filter** `length >= minLength` (`:382, :403`) — inclusive: run == minLength included,
  minLength−1 excluded.
- **Single-pass** ⇒ regions non-overlapping, sorted by Start; `0 ≤ Start`, `End ≤ count−1 < L`,
  `End − Start + 1 = length ≥ minLength` (INV-1..4).
- **Long IDR** `sequence.Length > 30` (`:458`) matches the >30 Ward 2004 convention (40-mer → Long IDR;
  20 → Standard; 30 itself does **not** qualify).

### Defect found & fixed — histidine omitted from f₊

`ClassifyRegionFlavorMobiDbLite` originally computed `f_plus` over **{R,K} only**, with a docstring
claiming `f₊ = (R+K)/L` "verbatim" from the v3 source. But the v3 `states.py` translation table maps
**H → positive** as well: `f₊ = (R+K+H)/L`. For histidine-containing regions the two diverge —
`HHHHHHHHAA` (f₊ = 0.8) is **PositivePolyelectrolyte** per MobiDB-lite v3, but the old code returned
**WeaklyCharged**. The **16 prior flavour tests all used His-free sequences**, so the gap was invisible
(green-on-blind-spot).

**Fix (this session):** added `H` to the positive count, corrected the docstring
(`f₊ = (R+K+H)/L`, translation table cited), and added two sourced tests — **F5b** (`HHHHHHHHAA` → PPE)
and **F5c** (`HHHHDDDDAA` → Polyampholyte, H balancing D). The validated TOP-IDP **boundaries** and the
default `RegionType` / `Confidence` are untouched.

**No citable standard for CONFIDENCE.** MobiDB-lite / IUPred / PONDR report a per-*residue* disorder
score, not a per-*region* calibrated confidence; no published deterministic region-confidence formula
exists. `Confidence = (mean − 0.542)/(1 − 0.542)` remains a **declared first-principles heuristic**
(boundaries unaffected) — recorded honestly, not fabricated.

### Test-quality audit

`DisorderPredictor_DisorderedRegion_Tests.cs` — **24 tests**, asserting **exact** Start/End/Count and
exact MeanScore/Confidence (not "no throw"), deterministic, covering every Stage-A edge case
(empty, all-ordered, all-disordered, exact-minLen, below-minLen, leading, trailing, central,
multi-region, bounds/length invariants). Boundary fixture → **24 passed / 0 failed**; flavour fixture
(now incl. F5b/F5c) → **16 passed / 0 failed**. Full unfiltered suite → **Failed 0** (Genomics 18819).
No `*Fast`/delegate variants for region grouping; `PredictMoRFs` reuses `PredictDisorder`'s per-residue
scores consistently (out of scope here).

## Runtime enforcement (LimitationPolicy)

The unit's guarded branch — the **uncalibrated per-region `Confidence`** — has **minimum access mode
`Permissive`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default
`LimitationPolicy.DefaultMode = Moderate` it throws `SeqeronLimitationException`; use the boundaries
without a confidence to stay within the validated contract. Additive policy layer; the `CLEAN` verdict
is unchanged.

## Findings & follow-ups

- **Region grouping / threshold / coordinates / edge-cases: correct and sourced** (State CLEAN); no
  off-by-one, all worked examples reproduced.
- **One defect (His omitted from f₊ in the opt-in flavour typing) found and fully fixed** with two new
  sourced tests; the invisible-on-His-free-fixtures gap is closed.
- The default `RegionType` composition thresholds and the rescaled `Confidence` are **disclosed
  first-principles heuristics** (labels only, boundaries unaffected); MobiDB-lite's own 0.35 / 0.32
  thresholds are source-exact. **No follow-ups.**
