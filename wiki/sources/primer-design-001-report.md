---
type: source
title: "Validation report: PRIMER-DESIGN-001 (PCR primer-pair design)"
tags: [validation, primer, algorithm]
doc_path: docs/Validation/reports/PRIMER-DESIGN-001.md
sources:
  - docs/Validation/reports/PRIMER-DESIGN-001.md
source_commit: c4df9b3a138915b4803ed7b97fdeab6658963d04
ingested: 2026-07-15
created: 2026-07-15
updated: 2026-07-15
---

# Validation report: PRIMER-DESIGN-001

The two-stage validation report for `PrimerDesigner.DesignPrimers`, candidate
generation, and per-primer evaluation, summarized in [[primer-design]]. It is
governed by the [[validation-protocol]].

## Verdict and supported relationships

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN.** The report
confirms that candidates must pass hard length/GC/Tm/run/structure constraints
before a legacy additive score greedily ranks them; this supports the relation
to [[primer3-weighted-penalty-objective]] as the richer replacement objective.
It also verifies that pair compatibility requires a Tm difference of at most
5 °C and no primer dimer, so [[primer-design]] consumes the signals represented
by [[primer-dimer-thermodynamics-tm]] and the [[primer-structure-qc-screens]].
The filtered design run recorded **250 passed / 0 failed**.
