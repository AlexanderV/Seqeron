# Review campaign queue (orchestrator state)

The orchestrator (cloud session until 2026-10-02, then a local session — see LOCAL_RUN.md) keeps this file current
and commits it after every change. Batch definitions (units, owned files) are in README.md → "Batches".

## Standing user rules
- **Max 2 batches concurrently.** Launch the next pending batch only when a slot is free.
- **Pause on request.** "Пауза / стоп / только Bxx" from the user overrides the queue until the user says otherwise.
- **Rate limits:** ignore `allowed_warning`; when a session ends with a usage/session limit, relaunch it as a resume
  right after the reset time (never in a loop before the reset).
- **Finish everything doable:** a batch is done only when its `## Leftovers` has zero DOABLE items
  (LIMITATIONS proposals count as leftovers). Otherwise run a FINISHER for that batch before moving on.
- **Order inside a batch:** units → completeness-audit loop → duplication sweep → final audit → heavy tier once → report.
- **Test tiers:** full fast tier per unit; targeted tier per audit WP; full fast tier once per audit round; heavy tier once.
- **No lost work:** small WPs pushed immediately; auditor worklist persisted in Bxx.md; WIP branch
  `claude/stoic-maxwell-0olr2z-wip-<BATCH>` every ~20 min.
- **Queue order:** new batches first (B08 … B23, B25, B26), THEN finishers FIN-B01/B02/B03/B24
  (re-audit their reports first — later batches may have resolved items), THEN Phase 2 (cross-batch dedup),
  THEN final consolidation (VALIDATION_LEDGER, FINDINGS_REGISTER, wiki ingest, final summary).
- **Suffix-tree upgrade plan (ST-UPGRADE, stages 0–8):** proposed, NOT approved — do not start without the user's yes.

## State (2026-10-02 ~11:30 UTC)
- **Running:** B07 (cloud session_01BqQH1AyNFGsM2pKjjKax9Z, audit round 3 of the completeness loop; worklist in B07.md).
  To be moved to a local run (tools/review/prompts/B07.md) when the local setup is ready.
- **Paused by user (10-01 22:05):** B08 — 2/6 units done (CRISPR-PAM-001, CRISPR-GUIDE-001). Resume only on the user's word
  (prompt: tools/review/prompts/B08.md).
- **User rule since 10-02:** ONLY B07 runs until the user says otherwise.
- **Pending (in order):** B08 (paused) B09 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B25 B26,
  then FIN-B01 FIN-B02 FIN-B03 FIN-B24, then Phase 2, then consolidation.
- **Done:** B01, B02, B03, B24 (finishers pending, see leftovers notes in their reports), B04 (F1–F65), B05 (F1–F35),
  B06 (F1–F38, 43841d1, Leftovers none).

## Known cross-batch items to carry into Phase 2
- B04: `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` generator can draw w = k = 5 → draw w from [k+1, 60] (flaky ~1/30).
- B06 report "Cross-batch requests" (B01, B04, B05, B07, B09, B10, B17, B18, MCP owner) — k-mer counting routed to canonical KmerAnalyzer.
- Stray branch `claude/stoic-maxwell-0olr2z-wip-TEST` (proxy cannot delete it) — ignore or delete in the GitHub UI.

## Log
- 10-02 06:54–07:00 README: no-lost-work rule (90d60c7), WIP branch (a0272b2, de7d295); B07 relaunched with them.
- 10-02 ~11:30 README: targeted test tier (9dc5a87); local-run kit added (tools/review/, LOCAL_RUN.md).
