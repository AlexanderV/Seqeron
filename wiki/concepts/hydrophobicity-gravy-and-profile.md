---
type: concept
title: "Hydrophobicity — Kyte-Doolittle GRAVY index & sliding-window hydropathy profile"
tags: [sequence-statistics, protein, algorithm]
sources:
  - docs/Evidence/SEQ-HYDRO-001-Evidence.md
  - docs/Evidence/SEQ-MW-001-Evidence.md
source_commit: e058738ff312bb90e5022081cf85e0b9da5b67cb
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-hydro-001-evidence
      evidence: "Test Unit ID: SEQ-HYDRO-001 ... Algorithm: Hydrophobicity Analysis (Kyte-Doolittle GRAVY index and sliding-window hydropathy profile)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:transmembrane-helix-prediction
      source: seq-hydro-001-evidence
      evidence: "SEQ-HYDRO-001 uses the same Kyte-Doolittle (1982) scale and the same unweighted sliding-window mean (Biopython protein_scale edge=1.0); its GCAT/Davidson source gives window 19 with 'peaks with scores greater than 1.6' for transmembrane detection — the segment-calling application PROTMOTIF-TM-001 builds on top of this profile."
      confidence: high
      status: current
---

# Hydrophobicity — Kyte-Doolittle GRAVY index & sliding-window hydropathy profile

**Hydrophobicity analysis** scores *how water-avoiding* an amino-acid sequence is, using the
classic **Kyte-Doolittle (1982) hydropathy scale**. The **SEQ-HYDRO-001** unit
([[seq-hydro-001-evidence]]) validates two related outputs: a whole-sequence scalar **GRAVY
index** and a per-position **sliding-window hydropathy profile**. [[test-unit-registry]] tracks
the unit; [[algorithm-validation-evidence]] describes the artifact pattern.

This is a **protein-property member of the SEQ-\* sequence-statistics family** — the
amino-acid hydropathy analogue of the nucleotide [[base-composition]] tally, and a whole-sequence-
scalar sibling of [[molecular-weight]] (SEQ-MW-001), which sums the same 20 residues' Daltons
instead of averaging their hydropathy. It shares the
exact Kyte-Doolittle scale and the unweighted-window-mean mechanism with the *segment-calling*
[[transmembrane-helix-prediction]] (PROTMOTIF-TM-001) and with the `CalculateHydropathy`
utility inside [[intrinsic-disorder-prediction-top-idp]] — **same scale, different output**:
GRAVY is a single number over the whole protein and the profile is the raw per-window signal,
whereas the TM unit thresholds that same profile (W=19, mean ≥ 1.6) to emit membrane-spanning
segments.

## The two outputs

1. **GRAVY (grand average of hydropathy)** = `Σ (kd value of each residue) / length` — the
   sum of per-residue Kyte-Doolittle hydropathy divided by the number of residues (Expasy
   ProtParam; Biopython `gravy()`). **Positive ⇒ hydrophobic, negative ⇒ hydrophilic.**
2. **Sliding-window hydropathy profile** — for a window length `W`, exactly **N − W + 1**
   points, each the **unweighted arithmetic mean** of that window's per-residue values
   (Biopython `protein_scale(edge=1.0)` — every position weighted equally). Empty when
   `W > N`.

## The Kyte-Doolittle scale (all 20 residues)

Verbatim `kd` values (Kyte & Doolittle 1982, via Biopython `ProtParamData.py`) — most
hydrophobic first: **I 4.5, V 4.2, L 3.8, F 2.8, C 2.5, M 1.9, A 1.8**; then **G −0.4, T −0.7,
S −0.8, W −0.9, Y −1.3, P −1.6, H −3.2**; most hydrophilic **N/D/Q/E −3.5, K −3.9, R −4.5**.
(Identical table to [[transmembrane-helix-prediction]]'s `HydropathyScale`.)

## Window-size guidance (source-anchored)

- **W = 9** — "best results" for **surface regions** of globular proteins (GCAT/Davidson).
- **W = 19** — "needed" for **transmembrane** detection; membrane segments are "peaks with
  scores greater than 1.6" — the parameter pair that [[transmembrane-helix-prediction]] adopts
  as its defaults.

## Canonical oracles

Closed-form over the tabulated scale (no library run needed):

| Input | GRAVY | Derivation |
|-------|-------|-----------|
| `A` | **1.8** | 1.8 / 1 |
| `AG` | **0.7** | (1.8 + −0.4) / 2 |
| `FLIV` | **3.825** | (2.8 + 3.8 + 4.5 + 4.2) / 4 |
| `RKDE` | **−3.85** | (−4.5 + −3.9 + −3.5 + −3.5) / 4 |

Window = 3 profile: `FLIV` → 2 windows **[3.7, 4.1666666667]**; `AG` (W > N) → **empty**.

## Contract and the one deviation

- **Case-insensitive** — the implementation uppercases; the scale is keyed on uppercase.
- **`W > N` / empty / null** ⇒ empty profile and GRAVY 0 (`range(N−W+1)` = zero iterations).
- **Unknown-residue handling — the sole deviation.** Biopython `gravy()` raises **`KeyError`**
  on any residue outside the 20 canonical amino acids (B/Z/X, gaps, stop). This repository
  instead **skips** unrecognized residues — GRAVY divides by the count of *recognized*
  residues; the profile treats them as 0. Kyte-Doolittle 1982 and Expasy define values only
  for the 20 standard residues and are silent on ambiguity codes/gaps, so **neither rule is
  mandated**. It is an API-shape/robustness choice, not a scoring-constant change: every value
  produced for in-alphabet residues stays **exactly source-conformant** (algorithm doc §5.4).

## Scope

A **sequence-only composition scalar + raw profile** — it reports hydropathy, it does not
*call* features. Membrane-span segment detection (thresholding the W=19 profile) lives in
[[transmembrane-helix-prediction]]; the charge–hydropathy view of disorder is
[[intrinsic-disorder-prediction-top-idp]]. A [[research-grade-limitations|research-grade]]
implementation of the Kyte-Doolittle method.

## References

Kyte J. & Doolittle R.F. (1982) *J Mol Biol* 157(1):105–132 (hydropathy scale + GRAVY);
Expasy ProtParam (GRAVY = Σ hydropathy / residues); Biopython `Bio.SeqUtils.ProtParam`
(`gravy`, `protein_scale(edge=1.0)` — reference implementations); GCAT/Davidson DGPB
(window 9 surface / 19 transmembrane, peak > 1.6); alakazam (CRAN) `gravy`. Full citations in
[[seq-hydro-001-evidence]] (not duplicated here).
