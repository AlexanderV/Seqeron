---
type: concept
title: "Transmembrane helix prediction (Kyte-Doolittle hydropathy sliding window)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/PROTMOTIF-TM-001-Evidence.md
  - docs/Evidence/SEQ-HYDRO-001-Evidence.md
  - docs/Validation/reports/PROTMOTIF-TM-001.md
source_commit: 3e90ada29a32f385c5a7ffa5227e4471967a9915
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: protmotif-tm-001-evidence
      evidence: "Test Unit ID: PROTMOTIF-TM-001 ... Algorithm: Transmembrane Helix Prediction (Kyte-Doolittle hydropathy sliding window)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:intrinsic-disorder-prediction-top-idp
      source: protmotif-tm-001-evidence
      evidence: "Both use the Kyte & Doolittle (1982) hydropathy scale over a sliding window (this unit's HydropathyScale; the disorder anchor's CalculateHydropathy utility returns a window's mean Kyte-Doolittle hydropathy). Same scale, different purpose — membrane-span detection vs charge-hydropathy disorder view."
      confidence: medium
      status: current
---

# Transmembrane helix prediction (Kyte-Doolittle hydropathy sliding window)

Predicting **transmembrane (TM) helices** — the hydrophobic α-helical segments that span a
lipid bilayer — from an amino-acid sequence alone. Seqeron implements the classic
**Kyte-Doolittle (1982) hydropathy sliding window**: slide a window of length **19** along the
protein, score each window as the **arithmetic mean** of its per-residue Kyte-Doolittle
hydropathy values, and emit the contiguous runs whose window mean stays **above threshold 1.6**
as candidate membrane-spanning segments. Validated under test unit **PROTMOTIF-TM-001**; the
validation record is [[protmotif-tm-001-evidence]] and [[test-unit-registry]] tracks the unit.
See [[algorithm-validation-evidence]] for the artifact pattern.

This is a **distinct unit of the protein-motif / ProteinMotif family** (alongside
[[coiled-coil-prediction]], [[common-protein-motifs]],
[[protein-domain-and-signal-peptide-prediction]] and
[[protein-motif-pattern-search]]). It keys off **hydrophobicity averaged over a window** —
the same Kyte-Doolittle scale used by the `CalculateHydropathy` utility inside
[[intrinsic-disorder-prediction-top-idp]] (which applies it to the charge–hydropathy view of
disorder), but here the signal is the sustained hydrophobic stretch long enough to cross the
bilayer. This unit is the **segment-calling layer over the raw hydropathy profile** validated
by [[hydrophobicity-gravy-and-profile]] (SEQ-HYDRO-001) — that concept owns the same
Kyte-Doolittle scale, the identical `edge=1.0` window mean, and the whole-sequence GRAVY
scalar, and its GCAT/Davidson source is where this unit's `W=19` / threshold `1.6` defaults
come from. It differs from [[coiled-coil-prediction]] (which scores heptad **a/d** periodicity,
not raw window-mean hydropathy) and from the signal-peptide model in
[[protein-domain-and-signal-peptide-prediction]] (a cleaved N-terminal targeting sequence, not
a membrane-spanning segment).

## Method

For an uppercased sequence `S` of length `n`, a window length `W` (default 19) and a threshold
`t` (default 1.6):

- **Hydropathy scale.** Each of the 20 standard residues has a fixed Kyte-Doolittle value
  (`HydropathyScale`): most hydrophobic I 4.5, V 4.2, L 3.8, F 2.8, C 2.5, M 1.9, A 1.8; then
  G −0.4, T −0.7, S −0.8, W −0.9, Y −1.3, P −1.6; most hydrophilic H −3.2, and D/E/N/Q −3.5,
  K −3.9, R −4.5.
- **Profile.** For each window start `i ∈ [0, n − W]`, the profile point is the **arithmetic
  mean** of the window's per-residue values: `score(i) = mean(scale[S[i..i+W−1]])`. The profile
  has `n − W + 1` points. Non-standard residues (X, B, Z, `*`) carry no scale value and are
  **excluded from the mean** (mean of the residues that do have values). Biopython's
  `protein_scale(edge=1.0)` reduces to this same mean.
- **Segment call.** A window is membrane-spanning iff `score(i) ≥ t`. Contiguous passing
  windows `[i₀ … i₁]` map to the residue segment `[i₀, i₁ + W − 1]` — the union of residues
  covered by any passing window — reported with `Score` = the peak profile value in the run.

The 19-residue window and the minimum span equal to the window width are justified by the
biophysics: a single α-helix needs ~18–21 residues (≈0.15 nm each) to cross the ~3–4 nm bilayer.

## Parameters and defaults

| Parameter | Default | Role / source |
|-----------|---------|---------------|
| `windowSize` | **19** | sliding-window length "needed" to locate a TM region (Davidson / Kyte-Doolittle); sequences shorter than this yield nothing |
| `threshold` | **1.6** | TM regions are "peaks with scores greater than 1.6" (Davidson, verbatim) |
| scale | Kyte-Doolittle 1982 | the 20 fixed per-residue values (QIAGEN + Davidson tables, matching exactly) |

Both `windowSize` and `threshold` are caller-overridable (a stated parameter — lowering the
threshold raises the segment count deterministically).

## Invariants

- **INV-01** — the profile has exactly `n − W + 1` points; empty when `n < W`.
- **INV-02** — every reported `Score ≥ threshold` and `Start ≤ End`.
- **INV-03** — a segment's `End = lastPassingProfileIndex + W − 1`, clamped to `n − 1` (see
  reporting convention below).
- **INV-04** — sequences shorter than `windowSize` (and null/empty) produce **no segments**.

## Canonical oracles and corner cases

Window means are closed-form over the tabulated scale, so exact oracles are hand-derivable:

- **Single hydrophobic stretch** — `D×10 + L×20 + D×10` (40 aa), W 19 / t 1.6: profile length
  22; passing windows run from profile index **5** to **16**; last passing window (start 16)
  covers residues 16..34 → segment **(Start 5, End 34)**; peak **3.8** (all-L window mean).
- **All-hydrophilic** — `D×40` (D = −3.5) → every window mean = −3.5 < 1.6 → **no segments**.
- **Exactly one poly-Leu window** — `L×19` → profile length 1, value 3.8 ≥ 1.6 →
  **(Start 0, End 18, Score 3.8)**.
- **Scale lookup** — a 19-residue uniform window reproduces the residue's own published value.

Other cases: null / empty / shorter-than-window input → empty; non-standard residues (X, B, Z,
`*`) are excluded from the window mean.

## Reporting convention (the one assumption)

The published rule defines the threshold-crossing *window run* but not how the window-indexed
run maps to exact residue boundaries. This unit reports
`End = lastPassingProfileIndex + windowSize − 1` (clamped to the last residue index) — the last
residue actually covered by an above-threshold window. A **2026-06-16 validation corrected an
earlier off-by-one** (`lastProfileIndex + windowSize`, which named a residue no passing window
covered). This is the sole documented assumption; there are no deviations from the scale or the
mean-windowing rule. The two-stage validation verdict is recorded in [[protmotif-tm-001-report]] —
**Stage A PASS-WITH-NOTES / Stage B PASS, State ✅ CLEAN** (the off-by-one End was fixed in-session;
full suite 6579/0).

## Scope

A **sequence-only, hydropathy-only heuristic** — the Kyte-Doolittle window is the classic first
approximation, not a topology model. It does not infer helix **orientation / in-out topology**,
does not apply the "positive-inside" rule, and does not use the HMM grammar or per-position
statistics of TMHMM / Phobius. It flags hydrophobic membrane-span candidates; a
[[research-grade-limitations|research-grade]] implementation.

## References

Kyte J. & Doolittle R.F. (1982) *J Mol Biol* 157(1):105–132 (hydropathy scale + sliding-window
method); Davidson College DGPB Kyte-Doolittle background (window 19, threshold 1.6, mean
windowing rule); QIAGEN CLC Genomics Workbench Hydrophobicity scales + Davidson per-amino-acid
scores (the 20 scale values); Biopython `Bio.SeqUtils.ProtParam` (`protein_scale`, `gravy` —
mean-windowing reference implementation). Full citations in [[protmotif-tm-001-evidence]] (do
not duplicate here).
