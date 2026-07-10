---
type: source
title: "Evidence: SEQ-SECSTRUCT-001 (protein secondary structure — Chou-Fasman Pα/Pβ/Pt sliding-window propensity profile)"
tags: [validation, sequence-statistics, protein]
doc_path: docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md
sources:
  - docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md
source_commit: 60297de1f42b4812ee249dc1773b8d88d89aa0d5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-SECSTRUCT-001

The validation-evidence artifact for test unit **SEQ-SECSTRUCT-001** — **protein secondary
structure prediction** by the **Chou-Fasman conformational propensities** Pα / Pβ / Pt
(α-helix / β-sheet / β-turn), evaluated here as a generic **sliding-window mean-propensity
profile** over an amino-acid sequence. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]] tracks the
unit. The propensity table, the profile contract, the oracles and the three assumptions are
synthesized on the concept [[protein-secondary-structure-chou-fasman]].

## What this file records

- **Online sources:**
  - **Wikipedia — Chou–Fasman method** (rank 4, cites Chou & Fasman 1974) — the method assigns
    each residue Pα/Pβ/Pt derived from relative frequencies in known structures; the original
    nucleation rules ("four of any six contiguous residues nucleate a helix", helix-former
    threshold 1.03; "three of any contiguous five" for a sheet, threshold 1.00; turn
    p(t)=p_t(j)·p_t(j+1)·p_t(j+2)·p_t(j+3), cutoff 7.5e-3).
  - **Kelley bioinfo lecture** (rank 4, PDF) — P(a)/P(b)/P(turn) = observed/expected propensity
    (×100 integer convention); the verbatim helix state machine (STEP 1 find 4/6 with P(a)>100
    in a window of 6; STEP 2 extend until 4 contiguous P(a)<100; STEP 3 mean P(a)/P(b) over the
    region; STEP 4 helix iff length>5 and P(a)>P(b)); worked values A 142/83, E 151/37, S 77/75,
    P 55/55, D 101/54.
  - **CSB|SJU (Jakubowski)** (rank 4) — verbatim helix Pα / sheet Pβ for all 20 residues with
    former/breaker classes; **lists Lys Pα = 1.16** (the contested value).
  - **Przytycka (NCBI/NLM) lecture** (rank 4, PDF) — verbatim Pa Pb Pt including **Lys 1.14 0.74
    1.01**, Ala 1.42 0.83 0.66, Asn 0.67 0.89 1.56, Phe 1.13 1.38 0.60, Tyr 0.69 1.47 1.14.
  - **ravihansa3000/ChouFasman reference implementation** (rank 3) — the full 20-residue Pa/Pb/Pt
    table in the integer ×100 convention (Ala 142/83/66 … Val 106/170/50, **Lys 114/74/101**);
    turn position frequencies f(i)..f(i+3) also tabulated but **not used** by the sliding-window
    mean profile under test (only Pa/Pb/Pt are averaged).
  - **BMC Bioinformatics PMC1780123** ("Improved Chou-Fasman method", rank 1, peer-reviewed) —
    restates the original nucleation/extension rules and the segment thresholds (helix ⟨Pα⟩ >
    1.03 and ⟨Pα⟩ > ⟨Pβ⟩; strand > 1.05), citing Chou & Fasman 1974.
- **Datasets (hand-derived, closed-form over the 1978 propensity table):**
  - Single-residue tuples (Pa/Pb/Pt): A 1.42/0.83/0.66, E 1.51/0.37/0.74, V 1.06/1.70/0.50,
    N 0.67/0.89/1.56, **K 1.14/0.74/1.01**.
  - Window mean (window = sequence length): `"AE"` (W 2) → helix (1.42+1.51)/2 = **1.465**, sheet
    (0.83+0.37)/2 = **0.60**, turn (0.66+0.74)/2 = **0.70**; `"AEV"` (W 3) → helix
    (1.42+1.51+1.06)/3 = **1.330**, sheet **0.9666…**, turn **0.6333…**.
- **Corner cases / failure modes:** the table covers exactly the **20 standard residues** — X, B,
  Z, `*`, gaps have no defined propensity and are **excluded** from the average (a window of only
  unknown residues emits nothing); a **window larger than the sequence** yields no scan positions;
  the sliding window steps by 1 producing **N − W + 1** windows in N-terminus order;
  case-insensitive (input uppercased); null/empty → empty result; non-positive window → empty.

## Deviations and assumptions

Three documented assumptions, no source contradictions on the 19 uncontested residues:

1. **Lysine Pα conflict resolved to 1.14, not 1.16.** CSB|SJU lists 1.16; the Przytycka NCBI
   lecture **and** the ravihansa3000 reference implementation (integer 114) list 1.14. The value
   **1.14** is adopted — two independent retrieved sources (one a reference implementation) versus
   a single source for 1.16, and 114 is consistent with the integer convention used for every
   other residue. This is the **only contested value**; the remaining 19 residues are identical
   across all sources.
2. **Default window size = 7 is an API convenience, not a Chou-Fasman constant.** Chou-Fasman
   defines a 6-residue *helix* nucleation window and a 5-residue *sheet* nucleation window, not a
   single 7-residue averaging window. The unit under test is a **generic sliding-window
   mean-propensity profile** whose window length is caller-supplied; tests pass the window
   explicitly and verify the arithmetic mean rather than the published nucleation/extension state
   machine.
3. **Unknown-residue handling = skip-and-exclude.** No retrieved source specifies behaviour for
   non-standard residues inside a window; the implementation excludes them from the per-window
   count/average (a window of only unknown residues emits nothing) — the documented deterministic
   contract.

**Reliability caveat (Wikipedia):** the original parameters were derived from a small,
non-representative sample (29 proteins) and have limited accuracy (~50–60% Q3); the propensities
are nonetheless the formally defined Chou-Fasman values.

Recommended coverage (from the artifact): MUST — single-residue window = that residue's
(Pa,Pb,Pt); lysine window Pa = **1.14** not 1.16; multi-residue window = exact arithmetic mean
(`AE`, `AEV`); N − W + 1 windows in N-terminus order; case-insensitive; unknown residues excluded
(`AXE` W 3 averages only A and E). SHOULD — null/empty → empty; window > sequence → empty;
non-positive window → empty. COULD — helix-favouring vs sheet-favouring peptides show mean Pa >
mean Pb (and vice versa).
