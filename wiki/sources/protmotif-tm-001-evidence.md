---
type: source
title: "Evidence: PROTMOTIF-TM-001 (Transmembrane helix prediction — Kyte-Doolittle hydropathy sliding window)"
tags: [validation, protein]
doc_path: docs/Evidence/PROTMOTIF-TM-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-TM-001-Evidence.md
source_commit: 12eadbfe614a7835d3d88f0c1d2ea1a8cb209612
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-TM-001

The validation-evidence artifact for test unit **PROTMOTIF-TM-001** — **transmembrane
helix prediction** by the **Kyte-Doolittle (1982) hydropathy sliding window**: slide a
window of length 19 along the protein, score each window as the **arithmetic mean** of its
per-residue Kyte-Doolittle hydropathy values, and report the contiguous runs whose window
mean stays above threshold **1.6** as candidate membrane-spanning segments. It is one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the model, contract, invariants and worked oracles are synthesized in
[[transmembrane-helix-prediction]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Kyte J & Doolittle RF 1982**, *A simple method for displaying the hydropathic
    character of a protein* (rank 1 primary, *J Mol Biol* 157:105–132) — the origin of the
    hydropathy scale and the sliding-window method; cited by every other source here.
  - **Davidson College DGPB — Kyte-Doolittle background** (rank 4 educational, quoting the
    primary) — pins the two operating constants verbatim: **window size 19** "is needed" to
    locate a transmembrane region, and "**transmembrane regions are identified by peaks with
    scores greater than 1.6**". Also states the windowing rule: the profile point is the
    **average** of the window's per-residue scores, sliding **one residue at a time**.
  - **QIAGEN CLC Genomics Workbench — Hydrophobicity scales** (rank 3 vendor reference,
    citing Kyte & Doolittle 1982) — the exact 20 one-letter scale values used in the
    implementation's `HydropathyScale` (I 4.5, V 4.2, L 3.8, F 2.8, C 2.5, M 1.9, A 1.8,
    G −0.4, T −0.7, S −0.8, W −0.9, Y −1.3, P −1.6, H −3.2, and D/E/N/Q all −3.5, K −3.9,
    R −4.5).
  - **Davidson College — per-amino-acid scores page** (rank 4) — independent confirmation of
    the QIAGEN scale table (matches exactly).
  - **Biopython `Bio.SeqUtils.ProtParam`** (rank 3 reference implementation) —
    `protein_scale(..., edge=1.0)` reduces to the **arithmetic mean** of the window's scale
    values (weights all 1.0, normalised by `sum(weights)*2+1`), identical to the Davidson
    windowing rule and to this unit; `gravy()` confirms the Kyte-Doolittle convention
    `Σ scale[aa] / length`.
  - **Transmembrane α-helix length** (rank 4 textbook, via WebSearch) — a single α-helix needs
    ~**18–21 residues** to span the ~3–4 nm lipid bilayer (≈0.15 nm/residue), which justifies
    the 19-residue window and the minimum-span filter equal to the window width.

- **Datasets (documented oracles):** hydropathy is a **closed-form mean** over the tabulated
  scale, so exact oracles are hand-derivable:
  - **Single hydrophobic stretch** — `D×10 + L×20 + D×10` (40 aa), window 19 / threshold 1.6:
    profile length `40−19+1 = 22`; windows with mean ≥ 1.6 run from profile index **5** to
    **16**; the last passing window (start 16) covers residues 16..34 → reported segment
    **(Start 5, End 34)**; peak score **3.8** (any all-L window: mean of 19×3.8 = 3.8).
  - **All-hydrophilic negative control** — `D×40` (D = −3.5) → every window mean = −3.5 < 1.6
    → **no segments**.
  - **Exactly one window of poly-Leu** — `L×19` → profile length 1, value 3.8 ≥ 1.6 →
    segment **(Start 0, End 18, Score 3.8)**.

- **Corner cases / failure modes:** a window mean is undefined when fewer than `windowSize`
  residues exist (null / empty / shorter-than-window input → **no segments**); the scale is
  defined only for the 20 standard amino acids, so non-standard characters (X, B, Z, `*`) are
  **excluded from the window mean** (mean of the residues that do have values), per Biopython
  scale-coverage behaviour.

## Deviations and assumptions

**One assumption — segment `End` reporting convention.** The published rule defines the
threshold-crossing *window run* but does not prescribe how the window-indexed run maps back to
exact residue boundaries. This unit reports `End = lastProfileIndex + windowSize − 1` (clamped
to the last residue index) — the last residue actually covered by any above-threshold window
(the union of all passing windows' residue coverage). A 2026-06-16 validation corrected an
earlier off-by-one (`lastProfileIndex + windowSize`, which named a residue not covered by any
passing window). Expected boundary values in the tests derive from this corrected convention.

Recommended coverage (MUST): single hydrophobic stretch → exactly one segment at the expected
(Start, End) with peak 3.8; all-hydrophilic → no segments; a 19-residue uniform window
reproduces each published per-residue scale value (mean = the residue's own score);
null/empty/shorter-than-window → empty. SHOULD: non-standard residues excluded from the window
mean; custom threshold raises/lowers segment count deterministically. COULD: property — every
reported `Score ≥ threshold` and `Start ≤ End`. No contradictions among sources (the two
Davidson pages, QIAGEN and Biopython all agree on the scale and the mean windowing rule).
