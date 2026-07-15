---
type: source
title: "Evidence: PRIMER-TM-001-HAIRPIN (Intramolecular hairpin self-folding MFE + unimolecular hairpin Tm)"
tags: [validation, primer]
doc_path: docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md
sources:
  - docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md
source_commit: 6c16153e119b7de7b8958cdb6c9dfc7fb2d092a8
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PRIMER-TM-001-HAIRPIN

The validation-evidence artifact for test unit **PRIMER-TM-001** — the **hairpin /
secondary-structure** extension: intramolecular self-folding (stem + loop) MFE detection and
the **unimolecular hairpin Tm** of a single primer. It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The shared
SantaLucia unified nearest-neighbour (NN) DNA thermodynamics, the hairpin loop table, and the
concentration-independent hairpin Tm are synthesized alongside the dimer engine in
[[primer-dimer-thermodynamics-tm]] (same `thal.c`/ntthal hairpin mode, same `NnUnifiedParams`
stem stacks). See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **SantaLucia & Hicks (2004)**, *The Thermodynamics of DNA Structural Motifs* (rank 1,
    Annu Rev Biophys Biomol Struct 33:415–440) — the primary source for the hairpin model:
    - **Table 1** unified Watson-Crick NN ΔH°/ΔS° at 1 M NaCl — reused **verbatim** for the
      stem stacks (the repo's `NnUnifiedParams`, already validated under PRIMER-TM-001-NN); no
      new NN values are introduced by the hairpin folder.
    - **Table 4 hairpin-loop ΔG°37 increments** (1 M NaCl), verbatim by loop size:
      3→3.5, 4→3.5, 5→3.3, 6→4.0, 7→4.2, 8→4.3, 9→4.5, 10→4.6, 12→5.0, 14→5.1, 16→5.3,
      18→5.5, 20→5.7, 25→6.1, 30→6.3 kcal/mol.
    - **Loop enthalpy/entropy rule** (Table 4 footnote a): all loop ΔH° = 0; loop ΔS° derived
      as `ΔS° = −ΔG°37·1000/310.15` (destabilising loop ⇒ ΔG°37 > 0).
    - **Complete hairpin** = salt-corrected stem NN stacks (Table 1, Eq. 3) **+** loop energy
      (Eqs. 8–10); for loops ≥ 5 the total adds a terminal-mismatch increment (Eq. 10).
    - **Unimolecular hairpin Tm (Eq. 11):** `Tm = ΔH°·1000/ΔS° − 273.15` — **no**
      `R·ln(C_T/x)` strand-concentration term (contrast the bimolecular duplex Eq. 3).
    - **Steric minimum:** loops shorter than 3 nt are sterically prohibited.
    - **Length-3 / length-4 special handling** (Eqs. 8–9): triloop bonus + closing-A·T penalty
      (+0.5 kcal/mol, only when the loop is closed by an A·T pair) for length 3; tetraloop bonus
      + terminal-mismatch increment for length 4. The triloop/tetraloop bonus values and the
      terminal-mismatch increment live in the Annual-Reviews **supplementary** material (not in
      the article body, not bundled) → exposed as an opt-in caller-supplied increment.
    - **Jacobson-Stockmayer large-loop extrapolation** (Eq. 7):
      `ΔG°37(loop-n) = ΔG°37(loop-x) + 2.44·R·310.15·ln(n/x)`, x = longest tabulated loop ≤ n
      (coefficient 2.44 from DNA kinetics, preferred over the older 1.75).
  - **SantaLucia (1998)** PNAS 95:1460 unified NN set (rank 1) — cross-check; the Table 1
    values are identical to 2004 and already embedded as `NnUnifiedParams`.
  - **Vallone & Benight (1999)** (rank 1, PMID 10423551) — short-DNA-hairpin melting studies:
    hairpin Tm is **concentration-independent** over 0.5–260 µM (unimolecular transition),
    independently confirming Eq. 11's lack of a concentration term.
- **Datasets (documented hand-derived oracles):**
  - *Canonical `GGGCTTTTGCCC`* (4-bp stem GGGC/GCCC, 4-nt TTTT loop): stem NN steps GG,GG,GC →
    ΔH° = −25.8 kcal/mol, ΔS° = −64.2 cal/(K·mol); loop-4 ΔG°37 = 3.5, loop ΔH° = 0,
    loop ΔS° = −3.5·1000/310.15 = −11.28486216346929; **totals** ΔH° = −25.8,
    ΔS° = −75.48486216346927, ΔG°37 = −2.3883700000000054, **Tm = 68.6403836682880 °C** (Eq. 11).
  - *`GGGCAAAAAGCCC`* (same 4-bp stem, 5-nt AAAAA loop): loop-5 ΔG°37 = 3.3 → total
    ΔS° ≈ −74.83968, ΔG°37 ≈ −2.58837, **Tm ≈ 71.585 °C** (exercises a second loop size).
  - *Non-hairpin poly-A `AAAAAAAAAAAA`*: no Watson-Crick stem possible → no hairpin (null),
    Tm = NaN.
- **Corner cases / failure modes:**
  - **Loop < 3 nt** sterically prohibited — no hairpin closes a 0–2-nt loop.
  - **Length-3 / length-4 loops** require the supplementary triloop/tetraloop bonus + (length-4)
    terminal mismatch; without them the computed hairpin is the exact, sourced stem-stack +
    loop-initiation core, and the bonus is an opt-in caller increment.
  - **No complementary stem** (homopolymer / oligo with no ≥ 2-bp WC stem closing a ≥ 3-nt loop)
    → no structure / not computable.
  - **Unimolecular vs bimolecular:** the bimolecular duplex-initiation term (+0.2/−5.7)
    nucleates two separate strands and is excluded; the loop-initiation term is the unimolecular
    nucleation cost.

## Deviations and assumptions

- **ASSUMPTION — bimolecular duplex-initiation term excluded.** The hairpin = base-pair NN
  contributions + loop energy; the +0.2/−5.7 bimolecular initiation is a two-strand nucleation
  cost and is omitted (loop initiation is the unimolecular nucleation cost). Consistent with the
  UNAFold convention — the standard interpretation, not an invented value.
- **ASSUMPTION — terminal-A·T penalty not applied at the open stem end.** The hairpin equations
  (8–10) add only stem NN stacks + loop; the terminal-A·T penalty is a duplex-end term of Eq. 3.
  Omitted to keep the hairpin core exact and sourced; applying it would need the supplementary
  terminal handling.
- No contradictions among sources.

Recommended coverage (MUST): `GGGCTTTTGCCC` returns ΔH° = −25.8, ΔS° = −75.48486216346927,
ΔG°37 = −2.3883700000000054, Tm = 68.6403836682880 °C exactly; the folder **finds** the correct
hairpin (stem 4, loop 4, span covering the whole oligo), not a worse partial structure; poly-A →
no hairpin (null) + NaN Tm; loop ΔS° follows `ΔS° = −ΔG°37·1000/310.15` (a sign / T error must
fail); the hairpin Tm carries **no** concentration term (Eq. 11 regardless of strand conc).
SHOULD: 5-nt loop `GGGCAAAAAGCCC` uses the Table-4 size-5 increment (3.3); null/empty/non-ACGT →
null/NaN. COULD: Jacobson-Stockmayer extrapolation returns Table-4 values exactly at tabulated
sizes with monotonically increasing ΔG°37 between them.
