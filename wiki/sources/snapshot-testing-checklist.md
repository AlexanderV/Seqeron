---
type: source
title: "Checklist 05: Snapshot / Approval Testing (Verify)"
tags: [validation, testing, methodology]
doc_path: docs/checklists/05_SNAPSHOT_TESTING.md
sources:
  - docs/checklists/05_SNAPSHOT_TESTING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Checklist 05: Snapshot / Approval Testing (Verify)

The **P1** per-unit checklist for snapshot/approval testing — a 255-row table tracking golden
`.verified.txt` coverage. Synthesized in the concept [[snapshot-testing]]; part of the
[[validation-and-testing]] program tracked in the [[test-unit-registry]].

## What this file records

- **Framework:** Verify + VerifyNUnit. Serialise full output → committed `.verified.txt`; diff on
  every run; changes require explicit developer approval.
- **Workflow:** first run creates the golden file; later runs fail on diff; behaviour change →
  review diff and accept new snapshot.
- **Per-unit table:** **37 / 255 ☑ complete, 221 not started** — the least-complete of the ten.
  ~20 snapshot files under `Snapshots/` cover 20 families; 5 new files still needed (Composition,
  Kmer, Translation, Oncology), ~10 existing to extend.

## Deviations and contradictions

None internal. This is the concrete evidence that the testing program is **not** uniformly
complete — snapshot coverage is a real tracked gap, unlike the 258/258 invariant methodologies.
Overlaps [[characterization-testing]] in mechanism (golden master).
