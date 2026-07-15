---
type: source
title: "Checklist 03: Fuzzing"
tags: [validation, testing, methodology]
doc_path: docs/checklists/03_FUZZING.md
sources:
  - docs/checklists/03_FUZZING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 03: Fuzzing

The **P2** per-unit checklist for fuzzing — a 258-row table mapping every test unit to the fuzz
strategies and fuzz targets that must be exercised. Synthesized in the concept [[fuzzing]]; part
of the [[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Purpose:** feed random/invalid/boundary input to surface crashes, hangs, and unhandled
  exceptions; critical for file-format parsers and input-validation points.
- **Strategy legend:** RB (Random Bytes), TF (Truncated Fields), MC (Malformed Content), BE
  (Boundary Exploitation: 0, −1, MaxInt, empty), INJ (Injection: special chars, null bytes,
  unicode), OVF (Overflow: extreme lengths/nesting).
- **Starting point (recorded):** "**Zero for Seqeron.Genomics**" — only `SuffixTreeFuzzTests`
  (corruption headers) existed before the campaign.
- **Per-unit table:** all **258 units ☑**; each row lists strategies + concrete targets (empty
  string, non-ACGT chars, null, unicode, extremely long, pattern > sequence length).
- **Summary:** 258 complete; priority split 12 high (parsers + validation), 45 medium (boundary
  inputs), 29 lower (algorithm-specific edge cases).

## Deviations and contradictions

None. The "zero → 258" trajectory is the largest single coverage jump among the ten methodologies.
