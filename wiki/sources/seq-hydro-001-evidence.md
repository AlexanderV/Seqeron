---
type: source
title: "Evidence: SEQ-HYDRO-001 (hydrophobicity — Kyte-Doolittle GRAVY index + sliding-window hydropathy profile)"
tags: [validation, sequence-statistics, protein]
doc_path: docs/Evidence/SEQ-HYDRO-001-Evidence.md
sources:
  - docs/Evidence/SEQ-HYDRO-001-Evidence.md
source_commit: 3e90ada29a32f385c5a7ffa5227e4471967a9915
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-HYDRO-001

The validation-evidence artifact for test unit **SEQ-HYDRO-001** — **hydrophobicity
analysis**: the whole-sequence **Kyte-Doolittle GRAVY index** (grand average of hydropathy)
and the **sliding-window hydropathy profile** over an amino-acid sequence. It is one instance
of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern;
[[test-unit-registry]] tracks the unit. The formula family, contract, oracles and the one
deviation are synthesized on the concept
[[hydrophobicity-gravy-and-profile]].

## What this file records

- **Online sources:**
  - **Biopython `Bio/SeqUtils/ProtParamData.py`** (`kd` scale, rank 3) — the verbatim
    Kyte-Doolittle (1982) hydropathy values for all 20 residues (A 1.8 … I 4.5, V 4.2 …
    R −4.5, D/E/N/Q −3.5), attributed to Kyte & Doolittle, *J Mol Biol* (1982).
  - **Biopython `Bio/SeqUtils/ProtParam.py`** (`gravy`, `protein_scale`, rank 3) —
    `gravy()` = `total_gravy / self.length` (sum of per-residue scale values ÷ length);
    `protein_scale` loops `range(N − W + 1)` (so exactly **N − W + 1** window values) with
    default `edge=1.0` = an **unweighted mean** over each window (no edge down-weighting).
  - **Expasy ProtParam documentation** (rank 2, SIB) — GRAVY, verbatim: "the sum of hydropathy
    values of all the amino acids, divided by the number of residues"; positive ⇒ hydrophobic,
    negative ⇒ hydrophilic; values from Kyte & Doolittle 1982.
  - **GCAT Davidson — Kyte-Doolittle background** (rank 4) — window **9** "best" for surface
    regions of globular proteins; window **19** "needed" for transmembrane detection, where
    "transmembrane regions are identified by peaks with scores greater than 1.6".
  - **alakazam (CRAN) `gravy`** (rank 3) — default scale Kyte & Doolittle, citing the 1982
    *J Mol Biol* 157:105-32 primary.
- **Datasets (hand-derived, closed-form over the kd scale):**
  - GRAVY: `A` → 1.8/1 = **1.8**; `AG` → (1.8 + −0.4)/2 = **0.7**;
    `FLIV` → (2.8 + 3.8 + 4.5 + 4.2)/4 = **3.825**; `RKDE` → (−4.5 + −3.9 + −3.5 + −3.5)/4 =
    **−3.85**.
  - Window = 3 profile: `FLIV` → 2 windows `[(2.8+3.8+4.5)/3, (3.8+4.5+4.2)/3]` =
    **[3.7, 4.1666666667]**; `AG` window 3 → **empty** (W > N).
- **Corner cases / failure modes:** window **larger than the sequence** (`W > N`) ⇒ empty
  profile (`range(N−W+1)` = zero iterations); **empty/null** input ⇒ GRAVY 0 / empty profile;
  **case-insensitive** (implementation uppercases; scale keyed on uppercase); a transmembrane
  window (W=19) over a hydrophobic stretch exceeds 1.6.

## Deviations and assumptions

**One documented deviation — unknown-residue handling.** Biopython `gravy()` indexes the scale
dict directly and raises **`KeyError`** on any residue outside the 20 canonical amino acids
(ambiguity codes B/Z/X, gaps, stop are undefined by the kd scale). The **repository
implementation instead *skips* unrecognized residues** — GRAVY divides by the count of
*recognized* residues, and the profile treats them as 0 (algorithm doc §5.4). Kyte-Doolittle
1982 and the Expasy doc define values only for the 20 standard residues and are silent on
ambiguity codes/gaps, so **neither rule is mandated** by an authoritative source; this is an
API-shape/robustness choice, not a scoring-constant change — every value produced for
in-alphabet residues remains **exactly source-conformant**. **No source contradictions** —
Biopython, Expasy, GCAT/Davidson and alakazam agree on the scale, the GRAVY formula, and the
unweighted-window-mean rule.

Recommended coverage (from the artifact): MUST — GRAVY of single residue / short hydrophobic
& hydrophilic peptides = exact sum/length; case-insensitivity; sliding-window profile =
exactly N−W+1 unweighted window means; window > length ⇒ empty, empty/null ⇒ 0/empty. SHOULD —
unknown residues skipped (the Biopython-KeyError deviation). COULD — a W=19 window over a
hydrophobic stretch exceeds the 1.6 transmembrane threshold.
