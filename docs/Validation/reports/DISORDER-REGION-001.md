# Validation Report: DISORDER-REGION-001 — Disordered Region Detection (IDR grouping)

- **Validated:** 2026-06-25 (fresh re-validation after F4 MobiDB-lite flavour typing was added; supersedes the 2026-06-24 pass)   **Area:** ProteinPred
- **Canonical method(s):** public `PredictDisorder(sequence, windowSize, threshold, minRegionLength)` → `IdentifyDisorderedRegions` (private), `ClassifyDisorderedRegion` / `CalculateConfidence` (private); plus the opt-in public `ClassifyRegionFlavorMobiDbLite(regionSequence)` (MobiDB-lite 3.0 flavour typing)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (one fidelity defect found in `ClassifyRegionFlavorMobiDbLite` — histidine omitted from f₊ — fully fixed this session)

> **Re-validation 2026-06-25.** Reset to ⬜ pending because the campaign added the opt-in
> `ClassifyRegionFlavorMobiDbLite` (F4). Re-validated fresh against external first sources retrieved
> this session (Campen 2008 TOP-IDP, Das & Pappu 2013 PNAS, MobiDB-lite v3 `states.py`/`consensus.py`
> via `gh api`). Found and fixed a charge-class defect (His). End-state CLEAN.

This unit covers the **grouping of per-residue disorder calls into contiguous REGIONS**
(start/end/length, minimum-length filter, classification label, confidence). The upstream
per-residue TOP-IDP window scoring is DISORDER-PRED-001 and is treated here as an input.

## Stage A — Description

### Sources opened & what they confirm
- **Ward et al. (2004), J Mol Biol 337:635-645 (DISOPRED2)** — confirmed via web search:
  long disordered regions are defined as **>30 consecutive residues**, with prevalence
  **2.0% archaea / 4.2% bacteria / 33% eukaryotes**. This anchors the classic ≥/>30-residue
  "long IDR" convention used by the code (`sequence.Length > 30` → "Long IDR").
- **Campen et al. (2008), TOP-IDP (PMC2676888):** contiguous windows predicted disordered
  constitute a disordered region; **cutoff 0.542**, window 21. (Scoring detail belongs to
  DISORDER-PRED-001; reused here.)
- **van der Lee et al. (2014), Chem Rev (PMC4095912):** IDR boundaries are the order↔disorder
  transition; recognized compositional subtypes (proline-rich, acidic, basic, Ser/Thr-rich).
- **Dunker et al. (2001), PMID 11381529:** disorder-/order-promoting AA sets.

### Region-calling logic validated
1. **Consecutive-grouping rule:** maximal runs of residues flagged disordered (`score ≥ threshold`). Standard.
2. **Threshold:** 0.542 (TOP-IDP cutoff, Campen 2008). Confirmed.
3. **Minimum region length:** a small **configurable** minimum (impl default 5). External
   practice supports both ≥30 (long IDR) and a smaller short-IDR minimum; a configurable
   minimum is within accepted practice. The code additionally **labels** runs >30 as
   "Long IDR", correctly anchoring the 30-residue convention to Ward (2004).
4. **Coordinate convention:** 0-based, **inclusive** Start and End; length = End − Start + 1.
   Standard and internally consistent.

### Edge-case semantics (all defined & sourced)
no disordered residue → 0 regions; all disordered & len≥minLen → single region [0, L−1];
run < minLen → excluded; run == minLen → included (`>=`); trailing run to last residue →
captured by explicit end-of-loop branch (the documented off-by-one pitfall, spec §1.4.1).

### Independent cross-check (hand recompute, Python)
Re-implemented the window scoring + grouping independently (window=21, halfWindow=10,
P→1.0, W→0.0, threshold 0.542, minLen 5), reproducing the spec coordinates exactly:

| Case | Sequence | Expected (spec) | Independent recompute |
|------|----------|-----------------|-----------------------|
| M2 full   | P×30          | [0, 29]              | **[0, 29]** ✓ |
| M6 trail  | W10+P20       | [11, 29]             | **[11, 29]** ✓ |
| S2 lead   | P20+W30       | [0, 18]              | **[0, 18]** ✓ |
| S5 central| W15+P20+W15   | [16, 33]             | **[16, 33]** ✓ |
| M14 multi | W15+P20+W15+P20+W15 | two regions    | **[(16,33),(51,68)]** ✓ |

MeanScore / Confidence recomputed: P30 → 1.0 / 1.0; E30 → 0.866 / 0.707; K30 → 0.786 / 0.532;
S30 → 0.655 / 0.246 — all matching the test assertions.

### Findings / divergences
PASS-WITH-NOTES. The grouping rule, threshold, coordinate convention and edge-case semantics
are correct and sourced. The **notes** are the classification heuristics, already disclosed in
the spec (§6/§7) as internal design decisions with no published source: enrichment threshold
0.25, priority order Pro>Acidic>Basic>S/T>Long>Standard, confidence formula
`(mean−0.542)/(1−0.542)`, AA groups {E,D}/{K,R}/{S,T}. These affect only the `RegionType`/
`Confidence` labels, not the region boundaries/lengths, and are honestly flagged ⚠ in the spec.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`:
- `IdentifyDisorderedRegions` (lines 358–417) — single-pass scan.
- `ClassifyDisorderedRegion` (lines 428–462), `CalculateConfidence` (lines 470–475).

### Logic realised correctly? (evidence)
- **Open region:** first disordered residue sets `regionStart = i`, resets `scoreSum` (lines 372–377).
- **Close on order residue:** `length = i − regionStart`; emits `End = i − 1` (lines 381, 389). Correct inclusive boundary.
- **Trailing region:** explicit post-loop block (lines 400–416): `length = count − regionStart`, `End = count − 1`. **No off-by-one** — the documented pitfall is handled.
- **Min-length filter:** `length >= minLength` (lines 382, 403) — inclusive: run == minLength included (S3), minLength−1 excluded (S4). Correct.
- **Single-pass** ⇒ regions non-overlapping and sorted by Start (INV-4); Start ≥ 0, End ≤ count−1 < L (INV-1/2); End−Start+1 = length ≥ minLength (INV-3).
- **Long IDR:** `sequence.Length > 30` (line 458) matches the >30 Ward (2004) convention (40-residue → Long IDR; 20 → Standard IDR; 30 itself would not qualify).

### Cross-verification recomputed vs code
The 5 hand-computed cases above match the exact assertions in M2/M6/S2/S5/M14, which pass
against the actual code. M5 (30×P, minLen 31 → 0 regions) confirms the filter.

### Variant/delegate consistency
No `*Fast`/delegate variants for region grouping. `PredictMoRFs` reuses `PredictDisorder`'s
per-residue scores consistently (out of scope here).

### Test quality audit
`DisorderPredictor_DisorderedRegion_Tests.cs` — 24 tests. Assertions check **exact** Start/End/
Count and exact MeanScore/Confidence values (not "no throw"), are deterministic, and cover every
Stage-A edge case (empty, all-ordered, all-disordered, exact-minLen, below-minLen, leading,
trailing, central, multi-region, bounds/length invariants). High quality. `--filter
~DisorderPredictor_DisorderedRegion_Tests` → **24 passed / 0 failed**.

### Findings / defects
None. The grouping, threshold, minimum-length filter, 0-based inclusive coordinates and the
>30 long-IDR label faithfully realise the validated description.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — grouping/threshold/coordinates/edge-cases correct & sourced;
  classification *labels* rest on disclosed internal heuristics (boundaries unaffected).
- **Stage B: PASS** — code matches; no off-by-one; all worked examples reproduced.
- **State: CLEAN** — no defect found in the boundary logic; the labelling heuristic was revisited (below).
- **Tests:** boundary fixture → 24 passed / 0 failed; new flavor fixture → 16 passed / 0 failed.

## Fresh re-validation (2026-06-25) — TOP-IDP scale, Das & Pappu, MobiDB-lite v3

### External first sources retrieved THIS session
- **Campen et al. (2008) TOP-IDP** (PMC2676888 / eurekaselect.com / PubMed 18991772, Wikipedia
  citing it): order→disorder ranking **W,F,Y,I,M,L,V,N,C,T,A,G,R,D,H,Q,(K/S),E,P**; P most
  disorder-promoting (0.987), W most order-promoting (−0.884). Smoothing **window 21**, cutoff
  **0.542** as used by the code.
- **Das & Pappu (2013) PNAS 110(33):13392-13397, PMID 23901099** — diagram of states: FCR = f₊+f₋,
  NCPR = |f₊−f₋|; **FCR > 0.35** is the strong/weak boundary; strong polyampholytes have NCPR ≤ 0.35,
  strong polyelectrolytes have NCPR > 0.35 (R4 with |NCPR| > 0.3).
- **MobiDB-lite v3 reference source** (`BioComputingUP/MobiDB-lite@v3`, fetched verbatim via `gh api`):
  - `states.py:get_disorder_class` — translation table **intab='RKDEACFGHILMNPQSTVWY'**,
    **outab='PPNN____P___________'** ⇒ positive = **{R,K,H}**, negative = **{D,E}**; gate `fcr > 0.35`;
    PA when `ncpr <= 0.35 or (f_minus > 0.35 and f_plus > 0.35)`; else PPE if `f_plus > 0.35`, NPE if
    `f_minus > 0.35`; else WC.
  - `states.py:is_enriched(subset, threshold=0.32)` — enriched iff `s >= threshold` (inclusive 0.32).
  - `consensus.py:get_region_features` priority: **charge (PA/PPE/NPE) → C → P → G → Low-complexity → Polar{S,T,N,Q}**.

### Stage A — confirmed
TOP-IDP values, window 21, threshold 0.542; Das & Pappu charge thresholds 0.35; MobiDB-lite
enrichment 0.32 (inclusive `>=`), priority order, polar set {S,T,N,Q}. All verbatim-matched to source.
Minor note: the code's comment ranking writes "…Q,S,K,E,P" whereas Wikipedia's Campen ranking lists
"…Q,K,S,E,P" (S/K swapped). This is a comment-only annotation; the code's **numeric** values
(Q 0.318 < S 0.341 < K 0.586 < E 0.736 < P 0.987) are self-consistent and do not affect any
boundary or flavour call. PASS-WITH-NOTES.

### Independent cross-check (hand-computed THIS session)
**Windowed TOP-IDP boundaries** (window 21, halfWindow 10, threshold 0.542), reproduced exactly:
- W10+P20 → region **[11,29]** (pos 11 = 12P/21 = 0.5714 ≥ 0.542; pos 10 = 11P/21 = 0.5238 < 0.542).
- W15+P20+W15 → region **[16,33]** (pos 16/33 = 12P/21 = 0.5714; pos 15/34 = 0.5238 < 0.542).
- Homopolymer mean norms: P 1.000, E 0.866, K 0.786, S 0.655 — match the test assertions.

**FCR/NCPR + flavour** (hand-computed, then confirmed against code):
| Seq | f₊ | f₋ | FCR | NCPR | Expected flavour | Branch |
|-----|----|----|-----|------|------------------|--------|
| RKDERKDE | 0.50 | 0.50 | 1.00 | 0.00 | Polyampholyte | FCR>0.35, NCPR≤0.35 |
| RKRKRKRKRR | 1.00 | 0 | 1.00 | 1.00 | PositivePolyelectrolyte | f₊>0.35 |
| DEDEDEDEDD | 0 | 1.00 | 1.00 | 1.00 | NegativePolyelectrolyte | f₋>0.35 |
| SSTTNNQQAA | 0 | 0 | 0 | 0 | Polar | {S,T,N,Q}=0.8≥0.32 |
| PPPPAAAAAA | 0 | 0 | 0 | 0 | ProlineRich | P=0.4≥0.32 |
| **HHHHHHHHAA** | **0.80** | 0 | 0.80 | 0.80 | **PositivePolyelectrolyte** | **H is positive** |
| **HHHHDDDDAA** | **0.40** | 0.40 | 0.80 | 0.00 | **Polyampholyte** | **H balances D** |

### Stage B — DEFECT found and fixed (histidine in f₊)
`ClassifyRegionFlavorMobiDbLite` originally computed `f_plus` over **{R,K} only** and the docstring
claimed `f₊ = (R+K)/L` "verbatim" from the v3 source. The v3 `states.py` translation table maps
**H → positive** as well (`f₊ = (R+K+H)/L`). For histidine-containing regions the two diverge:
`HHHHHHHHAA` (f₊ = 0.8) is **PositivePolyelectrolyte** per MobiDB-lite v3 but the old code returned
**WeaklyCharged**. The 16 prior flavour tests all used His-free sequences, so the gap was invisible.

**Fix (this session):** added `H` to the positive count, corrected the docstring (`f₊ = (R+K+H)/L`
with the translation table cited), and added two tests — F5b (`HHHHHHHHAA` → PPE) and
F5c (`HHHHDDDDAA` → Polyampholyte, H balancing D). The validated TOP-IDP **boundaries** and the
default `RegionType`/`Confidence` are untouched.

**No citable standard for CONFIDENCE.** MobiDB-lite / IUPred / PONDR report a per-residue disorder
score, not a per-*region* calibrated confidence; no published deterministic region-confidence formula
exists. `Confidence = (mean−0.542)/(1−0.542)` remains a **declared first-principles heuristic**
(boundaries unaffected) — recorded honestly, not fabricated. This is the documented acceptable boundary.

### Test run
Full unfiltered `dotnet test Seqeron.sln -c Debug` → **Failed: 0** (Genomics 18819 passed incl. the 2
new His tests; total across projects all green). Changed source project builds with 0 warnings / 0 errors.

**End-state: CLEAN** — the one defect (His omission) was completely fixed in-session with sourced
test coverage; boundary logic remains correct; the no-confidence-standard is the declared heuristic boundary.

## Runtime enforcement (LimitationPolicy)

Under the default `LimitationPolicy.DefaultMode = Strict`, the uncalibrated per-region `Confidence` (use `PredictDisorderRegions` for the validated TOP-IDP boundaries without a confidence) throws `Seqeron.Genomics.Core.SeqeronLimitationException` (named limitation + workaround; see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement and `LimitationCatalog`). `Permissive` mode returns the historical best-effort value. This is an additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
