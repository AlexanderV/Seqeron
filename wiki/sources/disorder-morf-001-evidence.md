---
type: source
title: "Evidence: DISORDER-MORF-001 (MoRF prediction — dip-in-disorder heuristic)"
tags: [validation, analysis]
doc_path: docs/Evidence/DISORDER-MORF-001-Evidence.md
sources:
  - docs/Evidence/DISORDER-MORF-001-Evidence.md
source_commit: 765f3b80c80cdce808d04012bc1cfb1e421a4a36
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: DISORDER-MORF-001

The validation-evidence artifact for test unit **DISORDER-MORF-001** — **MoRF (Molecular
Recognition Feature) prediction**, the "dip within disorder" heuristic that flags short ordered
segments embedded in longer intrinsically disordered regions. This is the **second ingested unit
of the protein disorder / features family** (after [[disorder-lc-001-evidence|DISORDER-LC-001]]) and
one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. The algorithm itself is written up on the concept page
[[morf-prediction-dip-in-disorder]]; this file records the source trace and worked oracles. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (four, with authority ranks):**
  - **Mohan et al. 2006** "Analysis of molecular recognition features (MoRFs)" (*J Mol Biol*
    362(5):1043–59, PMID 16935303, rank 1) — MoRFs are "relatively short (**10–70 residues**),
    loosely structured protein regions" **within longer, largely disordered sequences** that
    "undergo disorder-to-order transitions" on partner binding; three structural types by bound
    conformation — **α-MoRFs** (α-helix), **β-MoRFs** (β-strand), **ι-MoRFs** (irregular).
  - **Cheng/Oldfield et al.** "Mining α-helix-forming MoRFs with cross species sequence alignments"
    (*Biochemistry*, PMC2570644, rank 1) — the **operational "dip" definition**: the heuristic
    "identifies short regions of order within longer regions of disorder – or 'dips' – in disorder
    prediction profiles"; an α-MoRF is "a short (**around 20 residues**) structural element",
    candidates retrieved as "regions of **30 residues or less**"; disorder profile uses **the
    threshold of 0.5** (values above 0.5 = disorder, a dip below it = relative predicted order);
    exact dip parameters attributed to Oldfield 2005.
  - **Oldfield et al. 2005** "Coupled folding and binding with α-helix-forming molecular recognition
    elements" (*Biochemistry* 44(37):12454–70, PMID 16156658, rank 1) — the MoRE "consists of a
    short region that undergoes coupled binding and folding within a longer region of disorder"; the
    **exact numeric dip-detection parameters** (precise flank lengths, ordered-run window) live in
    the **paywalled Methods** and could **not** be retrieved (see Assumptions).
  - **Wikipedia** "Molecular recognition feature" (rank 4, citing Mohan 2006) — MoRFs are "small
    (10–70 residues) intrinsically disordered regions" that "undergo a disorder-to-order transition
    upon binding to their partners"; disordered before binding, folded after.
- **Underlying per-residue disorder score:** the repository `PredictDisorder` uses the **TOP-IDP
  scale** (Campen et al. 2008, PMC2676888) normalized to `[0,1]` via `(prop − (−0.884)) / 1.871`
  (Campen Table 2); a residue is **ordered** when score `< 0.5`, **disordered** when `≥ 0.5`
  (PMC2570644 0.5 threshold). Sample normalized values: P 1.000 / E 0.866 (disordered), L 0.298 /
  I 0.213 / W 0.000 (ordered).
- **Dataset (documented oracle):** a **synthetic dip-in-disorder construct** — a short homopolymer
  window of an ordered residue (e.g. L, score 0.298) embedded inside long P/E disordered flanks →
  a single ordered dip flanked by disorder = **one MoRF** at the dip's coordinates. Flanks are made
  long enough that `PredictDisorder`'s window averaging lets interior residues reach the pure
  per-residue score.
- **Documented corner cases / failure modes:** fully-ordered protein → no surrounding disorder → no
  MoRF; fully-disordered protein → no ordered dip → no MoRF; ordered run **outside 10–70 residues**
  (Mohan length band) → not a MoRF; a dip at the sequence **terminus** (not flanked by disorder on
  both sides) → not a MoRF.
- **Recommended coverage:** dip-in-disorder → exactly one MoRF at the dip coordinates; fully
  ordered → none; fully disordered → none; sub-10 / super-70 ordered run → none; terminal dip →
  none; MoRF score ∈ `[0,1]` and increases with dip depth (distance below 0.5); two separate dips
  → two non-overlapping MoRFs; case-insensitive input; null/empty → empty.

## Deviations and assumptions

**One ASSUMPTION**, scoped to the flank-length detail only:

1. **Exact dip flank/length detection parameters.** Oldfield et al. 2005 defines the precise numeric
   dip parameters (flank length, ordered-run window) but the Methods section is **paywalled and
   could not be retrieved**. This unit therefore implements the fully-retrievable **qualitative
   criterion**: an **ordered run** (per-residue disorder score `< 0.5`) of **total length within the
   Mohan 10–70 residue band**, flanked on **both** sides by ≥1 disordered residue inside a predicted
   disordered region. The **threshold (0.5)**, the **length band (10–70)**, and the **"order within
   disorder" shape** are all source-traceable and are **not** assumptions — only the flank-length
   detail is a correctness-affecting modeling choice.

No source contradictions — Mohan 2006, Cheng/Oldfield (PMC2570644), Oldfield 2005, and Wikipedia
agree on the 10–70 residue length, the "short order within longer disorder" shape, the disorder-to-
order transition on binding, and (where retrievable) the 0.5 disorder threshold. The bounded MoRF
score in `[0,1]` increasing with dip depth is a documented derivation, not a departure from any
source.
