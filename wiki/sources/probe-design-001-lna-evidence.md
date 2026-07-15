---
type: source
title: "Evidence: PROBE-DESIGN-001-LNA (LNA-adjusted nearest-neighbour Tm)"
tags: [validation, primer]
doc_path: docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md
sources:
  - docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md
source_commit: 56951fd9ed19b4c55b2183678292d5041eb772e1
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROBE-DESIGN-001-LNA

The validation-evidence artifact for the **LNA (locked nucleic acid) variant** of test unit
**PROBE-DESIGN-001** — an **LNA-adjusted nearest-neighbour melting-temperature** calculation
that raises the Tm (and specificity) of a hybridization probe by locking one or more internal
bases. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesized algorithm is
folded into the PROBE-family anchor [[taqman-probe-design-rules]] (LNA/Tm-adjustment section).
It reuses the repository's SantaLucia-1998-unified DNA nearest-neighbour thermodynamics — the
same NN engine validated under PRIMER-TM-001 (see [[primer-dimer-thermodynamics-tm]]).

## What this file records

- **Online sources:**
  - **McTigue, Peterson & Kahn (2004),** *Biochemistry* 43(18):5388–5405 (rank 1, primary
    peer-reviewed) — hybridization ΔH°/ΔS°/Tm measured from absorbance melting curves for 100
    duplexes carrying **single internal LNA nucleotides**; reported as ΔΔH°/ΔΔS° increments for
    **all 32** LNA+DNA:DNA nearest neighbours. Internal LNA **raises Tm** ("the largest known
    increase in thermal stability of any modified DNA duplex"), sequence-dependently.
  - **MELTING 5 reference implementation** (rank 3) — the verbatim McTigue NN parameter table
    (`McTigue2004lockedmn.xml`, values in **cal/mol** and cal/(mol·K)) and
    `McTigue04LockedAcid.java`. Confirms the **additive-increment application model**: compute
    the plain DNA NN ΔH°/ΔS° first, then add the LNA increment per LNA-containing NN step; and
    the **terminal-LNA rejection** (`isApplicable` returns not-applicable for an LNA at duplex
    position 0 or last — "not established for terminal locked nucleic acids"). Key notation
    `<NN-with-L>/<complement>` where `L` follows the locked base (e.g. `TTL/AA`, `TLG/AC`).
  - **rmelting tutorial worked example** (rank 3) — `melting("CCATTLGCTACC", conc 1e-4,
    dnadna, Na 1, method mct04)` → **Tm 63.61426 °C**, ΔH −81100 cal/mol, ΔS −222.5 cal/(mol·K).
    Confirms the `L`-suffix notation (5th base of `CCATTGCTACC`, 0-based index 4, is locked).
- **Datasets:** the MELTING worked example above (single internal LNA), and the **verbatim
  32-key McTigue increment table** (ΔΔH cal/mol, ΔΔS cal/(mol·K)); increments range from
  strongly stabilizing (e.g. `GLA/CT` +3162 cal) to destabilizing per-step (e.g. `GLG/CC`
  −2844 cal), but the net effect on Tm is stabilizing for internal LNA.
- **Hand-derived oracle (implementation-independent):** duplex `CCATTGCTACC`, LNA at index 4,
  C_T 1e-4 M, [Na⁺] 1 M reference state. Base DNA NN (SantaLucia 1998 unified) ΔH° −80.8 /
  ΔS° −221.7; add LNA steps `TTL/AA` (+2326, +8.1) and `TLG/AC` (−1540, −3.0) → **ΔH° −80.014
  kcal/mol, ΔS° −216.6 cal/(mol·K)**. Tm = ΔH°·1000/(ΔS° + R·ln(C_T/4)) − 273.15, R = 1.9872:
  all-DNA **59.692 °C**, LNA-adjusted **63.528 °C**, **ΔTm +3.84 °C** (LNA raises Tm); agrees
  with MELTING `mct04` (63.614 °C) to **0.086 °C**.
- **Corner cases / failure modes:** terminal LNA (index 0 or last) has no McTigue increment →
  **not computable** (must not be treated as internal); per-LNA ΔΔ is strongly
  context-dependent (some ΔΔH negative) so tests must use the actual table, not one average;
  non-ACGT base → the underlying DNA NN lookup fails → not computable.

## Deviations and assumptions

- **ASSUMPTION: base DNA NN model = SantaLucia 1998 unified.** The library adds the McTigue
  increments to its own SantaLucia-1998-unified DNA NN ΔH°/ΔS° (`CalculateNearestNeighbor…`)
  rather than to McTigue's own reference DNA set, keeping the LNA Tm consistent with the rest
  of the library's NN Tm. The ~0.09 °C residual vs MELTING is attributable to this base-model
  choice (documented); the **increment values themselves are applied exactly as published**.

No source contradictions — McTigue 2004, MELTING, and the rmelting worked example agree; the
only residual (0.086 °C) is the documented base-NN-set difference. Recommended coverage —
**MUST:** ΔH°/ΔS° of `CCATTGCTACC` (LNA@4) = −80.014 / −216.6; LNA Tm = 63.528 °C and within
~0.1 °C of MELTING 63.614; LNA **raises** Tm vs all-DNA (63.53 > 59.69); a **terminal LNA is
rejected** (not-computable). **SHOULD:** all 32 increment keys present, negative-ΔΔH key
(`GLG/CC` −2844) applied with correct sign; null/empty/short/non-ACGT/out-of-range index →
not-computable. **COULD:** MGB shorter-probe window (13–20 nt) and 3'-MGB placement guidance.
