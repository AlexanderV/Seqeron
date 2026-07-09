---
type: source
title: "Evidence: DISORDER-LC-001 (Protein low-complexity region detection — SEG algorithm)"
tags: [validation, analysis]
doc_path: docs/Evidence/DISORDER-LC-001-Evidence.md
sources:
  - docs/Evidence/DISORDER-LC-001-Evidence.md
source_commit: 35bcacb8fa9a2080358233de559d2678d7600b14
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: DISORDER-LC-001

The validation-evidence artifact for test unit **DISORDER-LC-001** — low-complexity region detection
in **protein** sequences by the **SEG algorithm** (Wootton & Federhen 1993/1996). This is the **first
ingested unit of the protein disorder / features family** (DISORDER-LC / MORF / PRED / PROPENSITY /
REGION) and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern. The algorithm itself is written up on the concept page
[[protein-low-complexity-seg]]; this file records the source trace and worked oracles. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **NCBI BLAST `blast_seg.c`** (NCBI C++ Toolkit, rank 3 reference implementation) — verbatim
    defaults `kSegWindow = 12`, `kSegLocut = 2.2`, `kSegHicut = 2.5` (= trigger window W, trigger
    complexity K1, extension complexity K2); complexity = Shannon entropy in **bits** via
    `s_Entropy` (composition counts normalized by window length, converted to base-2 with
    `NCBIMATH_LN2`); two-stage scan (`s_EntropyOn` / `s_StateOn`) — trigger at entropy ≤ locut,
    extend while ≤ hicut.
  - **SEG program help** (GCG/Weizmann mirror + `ncbi-seg` Ubuntu manpage + rothlab genhelp, rank 3)
    — verbatim `-WINdow 12`, `-LOWcut (K1) 2.2`, `-HIGhcut (K2) 2.5`; complexity in **bits/residue**;
    stage-1 identifies segments with complexity ≤ K1, stage-2 extends into overlapping segments with
    complexity ≤ K2; max complexity for 20-letter alphabet = `log₂(20) ≈ 4.322` bits; complexity
    "defined by equation (3) of Wootton & Federhen (1993)".
  - **Wootton & Federhen (1993/1996)** (primary literature, rank 1, search-level confirmation via
    Semantic Scholar + Oxford Bioinformatics review) — "local compositional complexity"; three
    parameters W/K1/K2 with defaults **W=12, K1=2.2 bits, K2=2.5 bits**; complexity measure "can be
    described by Shannon's Entropy", consistent with the reference implementation.
- **Complexity formula:** `H = −Σ pᵢ·log₂(pᵢ)` (bits/residue), max `log₂(20) = 4.321928`.
- **Datasets (documented oracles):**
  - *Hand-derived window entropies (L = 12, matching `s_Entropy`):* homopolymer `QQQQQQQQQQQQ` → H
    0.0 (triggers); two-residue `AAAAAALLLLLL` → H 1.0 (triggers); four-residue `AAABBBCCCDDD` → H 2.0
    (triggers at K1=2.2 but not a strict K1=0.5); 12-distinct `ACDEFGHIKLMN` → H 3.584963 (> K2 →
    no trigger/extension).
  - *Maximum amino-acid complexity:* `log₂(20) = 4.321928` bits/residue.
- **Documented corner cases:** window longer than sequence → no window → no segments;
  maximal-complexity window (W distinct → `log₂(12) ≈ 3.585 > K2`) → neither triggers nor extends;
  homopolymeric/near-homopolymeric windows (H 0–~1 bit ≪ K1) = the canonical target case.
- **Recommended coverage:** homopolymer ≥ W → one full-span segment; all-distinct → no segments;
  `AAABBBCCCDDD` triggers at K1=2.2 not 0.5; two-residue block → single merged segment;
  two homopolymer runs split by a high-complexity spacer → two segments; `minLength` filter;
  sequence < W → empty; case-insensitivity.

## Deviations and assumptions

Two **ASSUMPTIONs**, both flagged in the artifact as deviations from Wootton & Federhen, neither
moving segment boundaries on the canonical oracle cases:

1. **Region-type label string (`"X-rich"` / `"X/Y-rich"`).** A presentation extension: the dominant
   residue when its fraction > 0.5 → `"X-rich"`, else the top two → `"X/Y-rich"`. SEG defines only
   *where* segments are, not a textual composition label; only the 50% dominance threshold affects
   the label produced and it does not change boundaries. Documented as a deviation, not source-derived.
2. **Greedy single-residue extension.** The reference extends by merging length-W extension windows
   with complexity ≤ K2; the repo grows the contig one residue at a time while the *whole growing
   segment's* entropy stays ≤ K2. For the homopolymer/dipeptide test cases the boundaries are fixed
   by the trigger spans and are identical — an implementation variant.

No contradictions among sources — the NCBI reference implementation, the GCG/manpage program
documentation, and the Wootton & Federhen primary literature agree on the parameters (W=12,
K1=2.2, K2=2.5), the Shannon-entropy bits/residue complexity measure, and the two-stage
trigger/extend scan.
