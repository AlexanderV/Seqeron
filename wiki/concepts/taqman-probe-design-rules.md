---
type: concept
title: "TaqMan hydrolysis-probe design rules"
tags: [primer, algorithm, validation]
sources:
  - docs/Evidence/PROBE-DESIGN-001-Evidence.md
source_commit: b44e3684e0fc29f5d1f2ec7c66342736ac3fd842
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: probe-design-001-evidence
      evidence: "Test Unit ID: PROBE-DESIGN-001 ... Algorithm: Hybridization Probe Design — TaqMan (5'-nuclease hydrolysis probe) rules (opt-in extension)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: probe-design-001-evidence
      evidence: "The probe-Tm gate uses the repository's existing (PRIMER-TM-001-validated) salt-adjusted Tm formula; sibling MolTools reagent-design unit in the PRIMER/MolTools family"
      confidence: high
      status: current
---

# TaqMan hydrolysis-probe design rules

The **TaqMan (5'-nuclease hydrolysis probe) design rules** applied to a qPCR/FISH
hybridization probe (test unit **PROBE-DESIGN-001**). This is an **opt-in extension** to the
generic hybridization-probe designer — a genuinely distinct algorithm in the PCR/MolTools
reagent-design family, sibling to the primer units [[primer-dimer-thermodynamics-tm]] and
[[primer3-weighted-penalty-objective]]. The literature-traced record is
[[probe-design-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Scope: an opt-in layer over the generic probe designer

The **generic** hybridization-probe designer — length / Tm / GC windows, homopolymer cap,
self-complementarity, hairpin, specificity — is the **unchanged default**, and its
thermodynamic Tm was validated under PRIMER-TM-001 (see [[primer-dimer-thermodynamics-tm]]).
PROBE-DESIGN-001 adds the **TaqMan-specific hard rules** as an opt-in mode. What makes it a
distinct algorithm rather than a cross-link is that these constraints are peculiar to a
**5'-nuclease hydrolysis (dual-labelled reporter/quencher) probe** and have no analogue in
primer selection.

## The TaqMan-specific constraints

Sourced from Applied Biosystems / Thermo Fisher (chemistry manufacturer, authority rank 2)
and PREMIER Biosoft (rank 3), corroborated by IDT and a ScienceDirect topic page:

1. **No G at the 5' end (hard rule).** A guanine adjacent to the reporter dye at the 5' end
   **quenches reporter fluorescence even after cleavage**, so it cannot be rescued by the
   5'-nuclease hydrolysis step — hence a hard rule, not a soft penalty. Exposed as a
   `NoGuanineAt5Prime` flag; any violation forces `PassesAll == false`.
2. **More Cs than Gs.** The probe should contain strictly more C than G.
3. **No run of ≥4 consecutive Gs.** Called out separately from the generic homopolymer cap as
   the worst-case homopolymer for a probe.
4. **G+C content 30–80%.**
5. **Length 18–22 nt (default).** Applied Biosystems / PREMIER Biosoft give 18–22; IDT/Thermo
   elsewhere allow up to 30. The tighter 18–22 is the default, exposed via
   `minLength` / `maxLength` so no un-citable value is hard-coded.
6. **Probe Tm ≥ primer Tm + 10 °C.** The probe must melt ~10 °C above the primer so it binds
   the template before Taq polymerase reaches it. The **caller supplies the primer Tm**; the
   probe Tm uses the repository's salt-adjusted formula
   `81.5 + 16.6·log₁₀[Na⁺] + 41·GC − 600/N` ([Na⁺] = 0.05 M), validated under PRIMER-TM-001 —
   **not** a TaqMan-specific nearest-neighbour calculation.

## Antisense strand fallback

If **no** compliant probe (no 5'-G, more C than G) exists on the sense strand, the design is
moved to the **complement (antisense) strand** — implemented by `SelectTaqManStrand`. On the
worked example the sense probe `GTTAGGGTTAGGGTTAGG` (5'-G, G≫C) is rejected and its
reverse-complement `CCTAACCCTAACCCTAAC` (5'-C, C≫G, no G-run) is selected instead.

## Worked oracles

Tm via the salt-adjusted formula at [Na⁺] = 0.05 M (see the Evidence table):

- `CCATCACCCTACATCACC` (18 nt, 5'-C, C=10 G=0, GC 0.5556, Tm 49.35 °C) → **passes all** when
  primer Tm ≤ 39.35 °C.
- `GCATCACCCTACATCACC` → **fails** (5'-G).
- `ACCCCGGGGACCCTACAT` → **fails** (GGGG run).
- `ACGGGAGGTAGGTAGGTA` (C=1, G=9) → **fails** (more G than C).
- `CCATCACCCTACATCA` (16 nt) → **fails** (length < 18).
- `CCCGCCCCGCCCCGCCCC` (GC 100%) → **fails** (GC > 80%).

## Assumptions and contract

- **ASSUMPTION: 18–22 nt default window** — parameterised, so the value is a default, not a
  hard-coded un-citable constant.
- **ASSUMPTION: probe-Tm gate uses the repository salt-adjusted Tm** — the "+10 °C above
  primer" rule is sourced; the Tm engine is the existing PRIMER-TM-001-validated formula, and
  the primer Tm is caller-supplied.

No source contradictions; the four vendor/reference sources corroborate the rule set
point-for-point.
