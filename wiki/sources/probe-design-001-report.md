---
type: source
title: "Validation report: PROBE-DESIGN-001 (hybridization probe design)"
tags: [validation, probe, algorithm]
doc_path: docs/Validation/reports/PROBE-DESIGN-001.md
sources:
  - docs/Validation/reports/PROBE-DESIGN-001.md
source_commit: c4df9b3a138915b4803ed7b97fdeab6658963d04
ingested: 2026-07-15
created: 2026-07-15
updated: 2026-07-15
---

# Validation report: PROBE-DESIGN-001

The two-stage validation report for generic hybridization-probe design,
specificity scoring, and the opt-in TaqMan extension summarized in
[[hybridization-probe-design]]. It follows the [[validation-protocol]].

## Verdict and supported relationships

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN.** The report
confirms that [[taqman-probe-design-rules]] are an opt-in extension that leaves
the generic designer unchanged. It traces the suffix-tree overload through
`CheckSpecificity`, supporting the relationship to
[[probe-offtarget-specificity-scan]], and independently checks the Wallace and
salt-adjusted formulas behind the relationship to [[melting-temperature]]. The
generic `ProbeDesigner` validation recorded **91 passed / 0 failed** before the
fresh TaqMan re-validation documented in the same source.

