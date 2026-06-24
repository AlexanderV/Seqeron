# Validation Report: DISORDER-REGION-001 — Disordered Region Detection (IDR grouping)

- **Validated:** 2026-06-24   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.IdentifyDisorderedRegions(predictions, threshold, minLength)` (private), `ClassifyDisorderedRegion(region)` (private), `CalculateConfidence(meanScore)` (private) — all exercised via public `PredictDisorder(sequence, windowSize, threshold, minRegionLength)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

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
- **State: CLEAN** — no defect found; no code change required.
- **Tests:** unit filter → 24 passed / 0 failed. (No code touched; full-suite run not required.)
