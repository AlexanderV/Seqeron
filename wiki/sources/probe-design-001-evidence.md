---
type: source
title: "Evidence: PROBE-DESIGN-001 (TaqMan hydrolysis-probe design rules)"
tags: [validation, primer]
doc_path: docs/Evidence/PROBE-DESIGN-001-Evidence.md
sources:
  - docs/Evidence/PROBE-DESIGN-001-Evidence.md
source_commit: b44e3684e0fc29f5d1f2ec7c66342736ac3fd842
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROBE-DESIGN-001

The validation-evidence artifact for test unit **PROBE-DESIGN-001** — the **TaqMan
(5'-nuclease hydrolysis probe) design rules** added as an **opt-in mode** on top of the
generic hybridization-probe designer. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for
how units are tracked. The synthesized algorithm lives on the concept page
[[taqman-probe-design-rules]]; it is a sibling of the PCR-primer units
[[primer-dimer-thermodynamics-tm]] and [[primer3-weighted-penalty-objective]] in the
PRIMER/MolTools reagent-design family.

The generic probe designer (length/Tm/GC windows, homopolymer, self-complementarity, hairpin,
specificity) is the **unchanged default** and its Tm formulas were validated under
PRIMER-TM-001; this file records only the TaqMan-specific additions.

## What this file records

- **Online sources:**
  - **Applied Biosystems / Thermo Fisher — "Designing a TaqMan Gene Expression Assay"**
    (rank 2, chemistry manufacturer / canonical guidance) — **no G at the 5' end** (interferes
    with reporter-dye fluorescence), **probe Tm ~10 °C above primer Tm** (probe binds template
    before Taq polymerase reaches it), and the **antisense-strand fallback** when a 5'-G is
    unavoidable on the sense strand.
  - **PREMIER Biosoft — "TaqMan probe design tips"** (rank 3, restating ABI guidance) —
    **length 18–22 bp**, **G+C 30–80%**, **more Cs than Gs and not a G at the 5' end**, and
    **no run of ≥4 consecutive Gs** (worst-case homopolymer).
  - **ScienceDirect Topics — "TaqMan"** (rank 1–3, peer-reviewed reference work) — the
    **5'-G quenching rationale**: a G adjacent to the reporter dye **quenches fluorescence even
    after cleavage**, which is why the 5'-G rule is treated as a **hard** rule.
  - **IDT — "Designing PCR primers and probes"** (rank 3, corroborating) — probe Tm ≈ 10 °C
    above primer Tm (primer Tm ≈ 58–60 °C); probe length 18–30 nt (the tighter 18–22 from
    ABI/PREMIER used as default); more Cs than Gs; runs of identical nt ≤ 4.
- **Datasets (hand-derived TaqMan rule examples):** Tm via the repository salt-adjusted formula
  `81.5 + 16.6·log₁₀[Na⁺] + 41·GC − 600/N` at [Na⁺] = 0.05 M (PRIMER-TM-001-validated).
  `CCATCACCCTACATCACC` (18 nt, 5'-C, C=10/G=0, GC 0.5556, Tm 49.3473) → passes all when
  primerTm ≤ 39.35; `GCATCACCCTACATCACC` → fails 5'-G; `ACCCCGGGGACCCTACAT` → fails GGGG run;
  `ACGGGAGGTAGGTAGGTA` (C=1/G=9) → fails more-G-than-C; `CCATCACCCTACATCA` (16 nt) → fails
  length < 18; `CCCGCCCCGCCCCGCCCC` (GC 100%) → fails GC > 80%; sense `GTTAGGGTTAGGGTTAGG` →
  RC `CCTAACCCTAACCCTAAC` → strand selection picks the antisense probe.
- **Corner cases / failure modes:** 5'-G unavoidable on sense strand → design on the complement
  (antisense) strand (`SelectTaqManStrand`); a run of ≥4 Gs flagged separately from the generic
  homopolymer cap.

## Deviations and assumptions

- **ASSUMPTION: 18–22 nt default length window** — ABI/PREMIER give 18–22 while IDT/Thermo
  allow up to 30; the tighter 18–22 is the default, exposed via `minLength`/`maxLength`, so no
  un-citable value is hard-coded.
- **ASSUMPTION: probe-Tm gate uses the repository salt-adjusted Tm** — the "+10 °C above primer"
  rule is sourced; the Tm engine is the existing PRIMER-TM-001-validated formula (not a
  TaqMan-specific nearest-neighbour calc), and the caller supplies the primer Tm.

No contradictions among the four sources — they corroborate the rule set point-for-point.
Recommended coverage — **MUST:** a 5'-G probe is flagged (`NoGuanineAt5Prime == false`,
`PassesAll == false`); the more-C-than-G rule; a run of ≥4 Gs; GC outside 30–80% and length
outside 18–22; the probe-Tm gate vs the supplied primer Tm (fails when Tm < primerTm + 10); a
fully compliant probe accepted (`PassesAll == true`); strand selection picks the
no-5'-G / more-C-than-G strand on the known example.
