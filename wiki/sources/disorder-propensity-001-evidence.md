---
type: source
title: "Evidence: DISORDER-PROPENSITY-001 (per-residue disorder propensity — raw TOP-IDP lookup + Dunker classification)"
tags: [validation, analysis]
doc_path: docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md
sources:
  - docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md
source_commit: 1934ccea3f2b16ce6f2ae8fb793ac8f0704a6500
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: DISORDER-PROPENSITY-001

The validation-evidence artifact for test unit **DISORDER-PROPENSITY-001** — **per-residue disorder
propensity**, the low-level lookup layer beneath the sliding-window predictor. This is the **fourth
ingested unit of the protein disorder / features family** (after [[disorder-lc-001-evidence|DISORDER-LC-001]],
[[disorder-morf-001-evidence|DISORDER-MORF-001]] and [[disorder-pred-001-evidence|DISORDER-PRED-001]])
and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. It does **not** introduce a new algorithm: the four in-scope methods are the **raw TOP-IDP
scale lookup + Dunker order/disorder classification primitives** that the shared
[[intrinsic-disorder-prediction-top-idp|`PredictDisorder` anchor]] normalizes and windows. This file
records the source trace and worked oracles; the algorithm write-up lives on that concept page. See
[[test-unit-registry]] for how units are tracked.

## Scope — the four primitive methods

Unlike DISORDER-PRED-001 (which tests the *normalized, windowed* score), this unit tests the
**per-residue primitives** directly:

- `GetDisorderPropensity(residue)` → returns the **raw, un-normalized TOP-IDP Table 2 value**
  (e.g. `W = −0.884`, `P = +0.987`), **not** the `[0,1]` normalized score used inside `PredictDisorder`.
- `IsDisorderPromoting(residue)` → `true` iff the residue is in the Dunker disorder-promoting set.
- `DisorderPromotingAminoAcids` → the set `{A, E, G, K, P, Q, R, S}` (8 members).
- `OrderPromotingAminoAcids` → the set `{C, F, I, L, N, V, W, Y}` (8 members).

The 4 ambiguous residues `{D, H, M, T}` are in **neither** public set (so `IsDisorderPromoting` is
`false` for them).

## What this file records

- **Online sources (three, with authority ranks):**
  - **Campen et al. 2008** "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic
    Disorder" (*Protein Pept Lett* 15(9):956–963, PMC2676888, PMID 18991772, **rank 1 / primary**) —
    the **exact 20-value TOP-IDP scale** (Table 2, order→disorder: W −0.884, F −0.697, Y −0.510,
    I −0.486, M −0.397, L −0.326, V −0.121, N 0.007, C 0.020, T 0.059, A 0.060, G 0.166, R 0.180,
    D 0.192, H 0.303, Q 0.318, S 0.341, K 0.586, E 0.736, P 0.987), the **anchor residues** (W most
    order-promoting, P most disorder-promoting), and — recorded **for completeness only** — the
    prediction cutoff `0.542` with `I_Top-IDP = −(⟨Top-IDP⟩ − 0.542)`, which governs `PredictDisorder`,
    **not** the four methods in this unit's scope.
  - **Wikipedia** "Intrinsically disordered proteins" (**rank 4**; used only for the Dunker et al.
    2001 primary classification it cites) — disorder-promoting `A,R,G,Q,S,P,E,K` (8), order-promoting
    `W,C,F,I,Y,V,L,N` (8), ambiguous `H,M,T,D` (4); the continuum ranking string
    `W,F,Y,I,M,L,V,N,C,T,A,G,R,D,H,Q,K,S,E,P` (agrees with Campen, cross-source consistency); property
    note that disorder-promoting residues are hydrophilic/charged and order-promoting ones
    hydrophobic/uncharged.
  - **Dunker et al. 2001** "Intrinsically disordered protein" (*J Mol Graph Model* 19(1):26–59, PMID
    11381529, **rank 1**) — used as the **primary-citation locator** confirming author/year/journal
    for the classification the Wikipedia source cites (the verbatim residue sets are extracted from
    Wikipedia, cited to this primary, not recalled).
- **Datasets:** (1) the raw TOP-IDP scale — all 20 standard residues with their Table 2 values;
  (2) the Dunker disorder/order/ambiguous partition (8 + 8 + 4 = 20).
- **Documented corner cases / failure modes:** the scale is defined for the **20 standard residues
  only** — no value for non-standard / ambiguity codes (B, J, O, U, X, Z) or gaps; **unknown residue
  → 0.0** (`GetValueOrDefault(..., 0)` implementation contract, *not* a source-defined value);
  **case handling** — input is upper-cased before lookup, so `'p'` and `'P'` return the same value.
- **Recommended coverage:** MUST — `GetDisorderPropensity` returns the exact Table 2 value for all
  20 residues; `IsDisorderPromoting` true for `{A,R,G,Q,S,P,E,K}` and false for the order-promoting
  and ambiguous residues; the two public sets are exactly `{A,E,G,K,P,Q,R,S}` and `{C,F,I,L,N,V,W,Y}`
  (8 each); the three classification sets are **pairwise disjoint and cover all 20 residues** (derived
  consistency check). SHOULD — unknown residue → 0.0; lowercase input equals uppercase. COULD — W is
  the global minimum (−0.884) and P the global maximum (+0.987) of `GetDisorderPropensity`.

## Deviations and assumptions

**Two documented assumptions (both implementation-side, not source contradictions):**

1. **Unknown-residue propensity = 0.0** — Campen (2008) defines values only for the 20 standard
   residues; returning `0.0` for any out-of-scale character is an implementation contract
   (`GetValueOrDefault(..., 0)`), tested as a documented contract, not against the source.
2. **Ranking-vs-value discrepancy for S/K** — the Campen and Wikipedia *rendered ranking strings*
   place `…Q, K, S, E, P`, whereas the per-residue Table 2 **values** give `S = 0.341 < K = 0.586`
   (so by value the order is `Q, S, K, E, P`). The **numeric Table 2 values are authoritative** and
   are what the implementation and tests use; the ranking string is a presentation-order artifact.
   No correctness impact on the four scope methods (which use values and set membership, not the rank
   string). See also the same note on [[disorder-pred-001-evidence|DISORDER-PRED-001]].

No source contradictions.
