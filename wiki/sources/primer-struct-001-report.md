---
type: source
title: "Validation report: PRIMER-STRUCT-001 (primer structure analysis)"
tags: [validation, primer, algorithm]
doc_path: docs/Validation/reports/PRIMER-STRUCT-001.md
sources:
  - docs/Validation/reports/PRIMER-STRUCT-001.md
source_commit: c4df9b3a138915b4803ed7b97fdeab6658963d04
ingested: 2026-07-15
created: 2026-07-15
updated: 2026-07-15
---

# Validation report: PRIMER-STRUCT-001

The two-stage validation report for hairpin, primer-dimer, 3′-stability,
homopolymer, and dinucleotide-repeat screens summarized in
[[primer-structure-qc-screens]]. It follows the [[validation-protocol]].

## Verdict and supported relationships

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** The report verifies these
as fast boolean/scalar candidate screens rather than a full thermodynamic
folding model, supporting the `alternative_to` relationship with
[[primer-dimer-thermodynamics-tm]]. It also traces their direct use by
`EvaluatePrimer` and `DesignPrimers`, supporting the relationship to
[[primer-design]]. The SantaLucia nearest-neighbor values and Primer3 anchor
cases are independently re-derived in the source report.
