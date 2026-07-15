---
type: source
title: "Evidence: SEQ-PI-001 (isoelectric point — pH at which a protein's net charge is zero)"
tags: [validation, sequence-statistics, protein]
doc_path: docs/Evidence/SEQ-PI-001-Evidence.md
sources:
  - docs/Evidence/SEQ-PI-001-Evidence.md
source_commit: 8a4f33ace0d47f5ad2116aa8e775cab5608ccfc2
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-PI-001

The validation-evidence artifact for test unit **SEQ-PI-001** — **isoelectric point (pI)
calculation**: the **pH at which a protein's net charge is zero**, found by evaluating the
Henderson–Hasselbalch net-charge function over the ionizable groups and locating its zero
crossing. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]] tracks the
unit. The charge formula, pKa table, bisection method, contract, oracles and the two assumptions
are synthesized on the concept [[isoelectric-point]].

## What this file records

- **Online sources (all authority rank 3, reference implementations/tools):**
  - **EMBOSS `iep` documentation** — purpose "Calculate the isoelectric point of a protein from
    its amino acid composition assuming that no electrostatic interactions change the propensity
    for ionization"; pI = the pH at which net charge = 0 (positive and negative charges balance).
    Supplies the **EMBOSS `Epk.dat` pKa table** (verbatim): N-terminus 8.6 (+), C-terminus 3.6 (−),
    C 8.5 (−), D 3.9 (−), E 4.1 (−), H 6.5 (+), K 10.8 (+), R 12.5 (+), Y 10.1 (−).
  - **Peptides R package `charge_pI.cpp`** (Osorio et al. 2015) — the **net-charge formula**
    (Henderson–Hasselbalch, Moore 1985): basic/positive groups contribute `+1 / (1 + 10^(pH − pKa))`,
    acidic/negative groups contribute `−1 / (1 + 10^(pKa − pH))`. N-terminus keyed `'n'` (positive),
    C-terminus `'c'` (negative), each added once per chain. Acidic set = {D, E, C, Y, C-term};
    basic set = {R, K, H, N-term}.
  - **Peptides `charge()` documentation** — lists nine selectable pKa scales (Bjellqvist, Dawson,
    EMBOSS, Lehninger, Murray, Rodwell, Sillero, Solomon, Stryer) and a **worked EMBOSS-scale
    example**: sequence `FLPVLAGLTPSIVPKLVCLLTKKC` → net charge 3.037398 @ pH 5, 2.914112 @ pH 7,
    0.7184524 @ pH 9.
  - **seqinr `computePI` documentation** — pI = "the pH at which the protein has a neutral charge";
    uses Bjellqvist et al. pK values, the same algorithm as ExPASy Compute pI. Worked value (Bjellqvist
    scale) pI of `ACDEFGHIKLMNPQRSTVWY` = 6.78454 — recorded **only** to document scale dependence
    (this unit targets the EMBOSS scale, NOT used as an expected value).
  - **ExPASy Compute pI/Mw documentation** — pI computed from Bjellqvist et al. (1993) pK values
    (immobilised-pH-gradient migration); pairs pI with molecular weight in one tool (see
    [[molecular-weight]]). Accuracy caveat: predictions for **highly basic and small proteins** can
    be problematic; poor buffer capacity increases error.

- **Datasets:**
  - **Net-charge reference** (Peptides EMBOSS scale): `FLPVLAGLTPSIVPKLVCLLTKKC` → 3.037398 / 2.914112 /
    0.7184524 at pH 5 / 7 / 9 — validates the charge formula + EMBOSS pKa (the repo charge function
    reproduces these to 6 dp).
  - **Derived pI values** (EMBOSS scale, bisection over [0,14] to ±0.01, in-session, traceable via the
    confirmed charge function): `FLPVLAGLTPSIVPKLVCLLTKKC` → **9.67**; `A` → **6.10** (termini only:
    midpoint of N 8.6 / C 3.6); `AG` → **6.10** (termini only); `D` → **3.75**; `K` → **9.70**;
    `DDDD` → **3.23**; `KKKK` → **11.27**; `ACDEFGHIKLMNPQRSTVWY` → **7.36** (one of each, EMBOSS scale).

- **Corner cases / failure modes:** **no electrostatic interactions** — each ionizable group titrates
  independently, so pI is a function of **amino-acid composition only, not sequence order** (permutations
  give identical pI); **small / highly basic proteins** — predicted pI may be inaccurate (accuracy caveat,
  not a correctness rule).

## Deviations and assumptions

**Two documented assumptions (no source-contradicting deviation of the charge model):**
1. **Empty / null → neutral 7.0.** No authoritative source defines pI for a zero-length protein (a real
   protein always has both termini). The repository input-guard convention (sibling SEQ-\* statistics
   return a neutral/zero sentinel for empty input) returns 7.0 — a non-correctness-affecting guard, not an
   algorithm output.
2. **pKa scale = EMBOSS.** Multiple published scales exist (Bjellqvist, EMBOSS, Lehninger, …) giving
   slightly different pI. The single-pKa-per-residue model here matches the **EMBOSS** scale (not the
   position-dependent Bjellqvist model), so the EMBOSS `Epk.dat` values are the authoritative constant set.
   The seqinr Bjellqvist worked value (6.78454) is therefore NOT an expected value for this
   implementation — recorded only to document scale dependence.

Recommended coverage (from the artifact): MUST — net charge of `FLPVLAGLTPSIVPKLVCLLTKKC` = 3.037398 /
2.914112 / 0.7184524 at pH 5/7/9; pI of the same ≈ 9.67 (basic); pI ∈ [0,14] invariant for any input;
termini-only `A`/`AG` → 6.10 (midpoint of N/C-term pKa). SHOULD — acidic-only `DDDD` ≈ 3.23 / basic-only
`KKKK` ≈ 11.27 (monotonic charge response); empty/null → 7.0. COULD — order-independence (permutation ⇒
identical pI, the "no electrostatic interactions" assumption). **No source contradictions** — EMBOSS,
Peptides, seqinr and ExPASy agree on the charge model; they differ only in the pKa scale, which is a
documented parameter choice.
