# Validation Report: DISORDER-REGION-001 — Disordered Region Detection (IDR grouping)

- **Validated:** 2026-06-12   **Area:** ProteinPred
- **Canonical method(s):** `DisorderPredictor.IdentifyDisorderedRegions(predictions, threshold, minLength)` (private), `ClassifyDisorderedRegion(region)` (private), exercised via public `PredictDisorder(sequence, windowSize, threshold, minRegionLength)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

This unit covers the **grouping of per-residue disorder calls into contiguous REGIONS** (start/end/length, minimum-length filter, classification). The upstream per-residue scoring (TOP-IDP window) is DISORDER-PRED-001 and was validated separately; here it is treated as an input.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Intrinsically disordered proteins"** and **Ward et al. (2004), DISOPRED2** (via Wikipedia citation): long (>30 residue) disordered **segments** occur in 2.0% archaeal / 4.2% eubacterial / 33.0% eukaryotic proteins. Confirms the classic **≥30 consecutive residues** threshold for a *long* IDR.
- **Web search (DisProt / IUPred2A / SPOT-Disorder usage):** an IDR is operationally defined as a **contiguous run of consecutive residues whose disorder score exceeds a threshold**; a "long disordered region" requires **≥30 consecutive** disordered residues. Short IDRs use a **smaller minimum** (commonly **≥4** consecutive putative disordered residues). This establishes that a *configurable* small minimum length is standard practice, not a fixed value.
- **Campen et al. (2008), TOP-IDP (PMC2676888):** contiguous windows predicted disordered constitute a disordered region; cutoff 0.542, window 21. (Scoring detail — DISORDER-PRED-001.)
- **van der Lee et al. (2014), Chem Rev (PMC4095912):** IDR boundaries are defined by the order↔disorder transition; recognized compositional subtypes (proline-rich, acidic, basic, Ser/Thr-rich).

### Region-calling logic validated
1. **Consecutive-grouping rule:** group maximal runs of residues flagged disordered (`score ≥ threshold`). Confirmed standard.
2. **Threshold:** 0.542 (TOP-IDP cutoff, Campen 2008). Confirmed.
3. **Minimum region length:** a small **configurable** minimum (impl default 5). External sources support both ≥30 (long IDR) and ≥4 (short IDR) — a configurable minimum is within accepted practice. The impl additionally **labels** runs >30 as "Long IDR", correctly anchoring the 30-residue convention to the Ward (2004) source.
4. **Coordinate convention:** 0-based, **inclusive** Start and End; length = End − Start + 1. Standard and internally consistent.

### Edge-case semantics (all have defined, sourced expectations)
no disordered residue → 0 regions; all disordered & len≥minLen → single region spanning whole sequence [0, L−1]; run < minLen → excluded; run == minLen → included (`>=`); trailing run to last residue → captured by explicit end-of-loop branch (the documented off-by-one pitfall).

### Independent cross-check (hand computation)
Re-implemented the grouping independently in Python (window=21, halfWindow=10, P normalized→1.0, W→0.0, threshold 0.542, minLen 5), reproducing the spec's expected coordinates exactly:

| Case | Sequence | Expected (spec) | Independent recompute |
|------|----------|-----------------|-----------------------|
| M6 trailing | W10+P20 | [11, 29] | **[11, 29]** ✓ |
| S2 leading  | P20+W30 | [0, 18]  | **[0, 18]** ✓ |
| S5 central  | W15+P20+W15 | [16, 33] | **[16, 33]** ✓ |

### Findings / divergences
PASS-WITH-NOTES. The grouping rule, threshold, coordinate convention and edge-case semantics are correct and sourced. The notes are the **classification heuristics**, already disclosed in the spec as internal design decisions (no published source): enrichment threshold 0.25, priority order Pro>Acidic>Basic>S/T>Long>Standard, confidence formula `(mean−0.542)/(1−0.542)`, AA groups {E,D}/{K,R}/{S,T}. These affect only the `RegionType`/`Confidence` labels, not the region boundaries/lengths, and are honestly flagged ⚠ in §6/§7 of the spec.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`:
- `IdentifyDisorderedRegions` (lines 367–426) — single-pass scan.
- `ClassifyDisorderedRegion` (lines 437–471), `CalculateConfidence` (lines 479–484).

### Logic realised correctly? (evidence)
- **Open region:** on first disordered residue, `regionStart = i` (line 383).
- **Close region (order residue):** `length = i − regionStart`; residues `regionStart..i−1`; emits `End = i − 1` (lines 390, 399). Correct inclusive boundary.
- **Trailing region:** explicit post-loop block (lines 409–425): `length = count − regionStart`, `End = count − 1`. **No off-by-one** — the documented pitfall (spec §1.4.1) is handled.
- **Min-length filter:** `length >= minLength` (lines 391, 412) — inclusive, so a run of exactly minLength is included (S3) and minLength−1 is excluded (S4). Correct.
- **Single-pass** guarantees regions are non-overlapping and sorted by Start (INV-4).
- Length invariant: emitted regions satisfy `End − Start + 1 = length ≥ minLength` (INV-3). `Start ≥ 0`, `End ≤ count−1 < L` (INV-1/2).

### Cross-verification recomputed vs code
The 3 hand-computed cases above match the exact assertions in M6/S2/S5 tests, which pass against the actual code. M2 (30×P → [0,29]) and M5 (30×P, minLen 31 → 0 regions) confirm full-span and filter behaviour.

### Variant/delegate consistency
No `*Fast`/delegate variants for region grouping. `PredictMoRFs` reuses `PredictDisorder`'s regions consistently (out of scope here).

### Test quality audit
`DisorderPredictor_DisorderedRegion_Tests.cs` — 24 spec rows + invariant tests (29 executed). Assertions check **exact** Start/End/Count and exact MeanScore/Confidence values (not "no throw"), are deterministic, and cover every Stage-A edge case (empty, all-ordered, all-disordered, exact-minLen, below-minLen, leading, trailing, central, multi-region, bounds/length invariants). High quality.

### Findings / defects
None. The grouping, threshold, minimum-length filter and 0-based inclusive coordinates faithfully realise the validated description.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — grouping/threshold/coordinates/edge-cases correct & sourced; classification *labels* rest on disclosed internal heuristics (boundaries unaffected).
- **Stage B: PASS** — code matches; no off-by-one; all worked examples reproduced.
- **State: CLEAN** — no defect found; no code change required.
- **Tests:** `--filter ...~DisorderedRegion` → 29 passed / 0 failed. Full suite: **4486 passed, 0 failed**.
