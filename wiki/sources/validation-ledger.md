---
type: source
title: "Validation ledger — live per-unit validation status tracker"
tags: [validation, testing, governance]
doc_path: docs/Validation/VALIDATION_LEDGER.md
source_commit: 255a8836226b75d0cded56b57ffdb14184737664
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: supersedes
      object: source:findings-register
      source: validation-ledger
      evidence: "The 2026-06-24 full re-validation reset in VALIDATION_LEDGER.md moved every unit back to pending, superseding the 2026-06-12 FINDINGS_REGISTER.md dispositions; the register itself states it was 'explicitly superseded on 2026-06-24' by the ledger and directs readers to the ledger as ground truth for current state."
      confidence: high
      status: current
---

# Validation ledger — live per-unit validation status tracker

The **ground-truth tracker of where each test unit's validation actually stands right now**. Where the
[[test-unit-registry]] defines the ID scheme and per-unit *spec/record*, and the [[findings-register]]
holds a *superseded snapshot* of one triage pass, this ledger is the **live status board**: a per-unit
pass/fail matrix, refreshed by reset banners whenever code churn invalidates prior results. Every other
governance page defers to it — the register, the [[validation-findings-disposition]] process, and
[[validation-and-testing]] all name `docs/Validation/VALIDATION_LEDGER.md` as the authority for
"where things stand."

## What a row records

Each unit gets two per-stage marks plus an end state, driven by the two-context
[[validation-protocol]] (one **fresh session per unit**, implementer and validator deliberately different
contexts, **Stage A description before Stage B implementation**, checked against **external primary sources**
with hand/reference cross-checks and mutation checks):

- **Stage A** (algorithm *description*) and **Stage B** (*implementation*): ✅ PASS · 🟡 PASS-WITH-NOTES · ❌ FAIL · ⬜ pending.
- **State** (end of session): ✅ CLEAN (fully functional) · 🔧 LIMITED (see report) · ↩︎ DUPLICATE-OF · ⬜ pending.

PASS-WITH-NOTES are documented by-design scope boundaries (declared heuristics, stale-spec wording), not
defects. The 86 original per-unit reports were committed once for provenance, then consolidated into this
ledger + the [[findings-register]]; any report is recoverable via `git show cb113ce:docs/Validation/reports/{UNIT-ID}.md`.

## Three phases

The ledger is partitioned into three campaigns, each with its own reset banner and progress line:

- **Phase 1 — implemented (☑) units.** A full **re-validation RESET on 2026-06-24** put all units back to
  ⬜ pending (prior Stage/State results SUPERSEDED, retained only as historical evidence); the same day it
  completed **86/86 ✅ CLEAN** (0 FAIL; Stage-A 63 PASS/23 notes, Stage-B 79 PASS/7 notes). One real
  defect was found and fixed in-session — **PARSE-GENBANK-001** (multi-line qualifier reconstruction adding
  spurious spaces to wrapped `/translation`). A **post-completion re-reset on 2026-06-25** then knocked
  **19** previously-CLEAN units back to pending after later limitation-elimination tiers (G–N) touched them,
  so the headline is no longer a clean 86/86. A separate section tracks **24 new campaign units**
  (McCaskill accessibility, Plan7 HMMER, ntthal dimer/hairpin Tm, MHCflurry/SMM/BIMAS/ν-SVR, CheckM, TETRA,
  MaxEntScan, context++/PCT, pre-miRNA classifier, α-satellite/HOR, …) pending first validation.
- **Phase 2 — 148 Phase-2 registry units.** Reset 2026-06-24, progress line reads 0/148 pending, but the
  campaign-result note records **13 genuine algorithm defects found and fully fixed** (e.g. SV-DETECT-001
  tandem-dup signature, SV-CNV-001 half-integer rounding, TRANS-SPLICE-001 A5SS/A3SS swap, MIRNA-PAIR-001
  DNA-T→U pairing, RNA-PARTITION-001 McCaskill probabilities), plus many green-washed tests re-anchored to
  externally-sourced literals. SEQ-COMPOSITION-001 is flagged ↩︎ DUPLICATE-OF SEQ-STATS-001.
- **Phase 3 — 12 enhanced units.** After the open-questions program shipped 8 enhancements (Group 1 features
  + Group 2 breaking changes), each touched unit is independently re-validated in a fresh session under the
  rule **"tests follow the source, code obeys the tests, tests are never bent to the code."** This pass
  caught **1 real latent defect** the implementer's own tests missed — **PHYLO-NEWICK-001** silently
  accepting unbalanced parentheses (`((A,B);` → wrong topology), now throws.

## LIMITED units and deferred big fixes

The two units once marked 🔧 LIMITED — **META-CLASS-001** (flat best-hit classifier → real Kraken
taxonomy-tree + k-mer-LCA + RTL) and **SPLICE-PREDICT-001** (spliced-sequence/exon-set consistency) — were
both **fixed to ✅ CLEAN**. The "deferred BIG fixes" backlog (guide-tree progressive MSA, circular
restriction digest, N-ary phylo trees, published CRISPR scoring models, EMBL remote refs) is now mostly
**DONE** as approved breaking/additive API changes, cross-referenced to the [[findings-register]]'s C-items.

## Where this fits

- [[validation-protocol]] — the two-stage, one-session-per-unit methodology whose sessions produce every row of this board.
- [[test-unit-registry]] — defines the units this ledger tracks the live state of (the ID scheme / per-unit spec vs. this pass/fail board).
- [[findings-register]] — the 2026-06-12 disposition snapshot this ledger **supersedes** via the 2026-06-24 reset.
- [[validation-findings-disposition]] — the A/B/C/D triage process; treats this ledger as ground truth for current state.
- [[validation-and-testing]] — the overall correctness strategy this ledger operationalizes unit-by-unit.
