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

## State (2026-10-09 03:55 UTC) — B07 DONE; EVERYTHING PAUSED by user
- **B07 done** (user 10-08 ~21:20: "B07 доходим до конца и останавливаемся"). resume5 (cloud
  session_01XittAhqmUUdF5idpt2X6rU, archived 10-09 03:53) ran audit rounds 3–10 (round 10 = zero DOABLE),
  the duplication sweep (DUP-1), the final audit (AF-4/5/7) and the heavy tier once, then wrote the final report
  (ecb6a84d). All 10 units FIXED (F1–F75); `## Leftovers`: none DOABLE; BLOCKED only (MGB ΔTm, LNA on both
  strands / terminal LNA, genome-wide primer specificity, extreme-score K runtime). Full `Seqeron.Genomics.Tests`
  24397/0, Mcp.MolTools 214/0. WIP branch -wip-B07 = main (clean). resume5 cost ≈ $158.9.
- **Everything paused by user:** nothing launched after B07 — no finisher, no B08, no other batch. B08
  (session_01NvVQb3FM1E4YGuhxEtXqjH) stays paused, not archived. Check-in trigger trig_01MvAagyZjUaqveNm9GXWGtm
  disabled. Old B07 session_01AJjD6znMPKvtM3iddFir9r idle, not archived. Branch carries 54 commits not yet in master
  (no PR yet — only on the user's request).
- Batch status: fully done B04 B05 B06 B07; units done but FINISHER needed B01 B02 (1 open item) B03 B24 (2 LIMITED);
  B08 2/6 paused; B09–B23, B25, B26 not run under the full process (first-pass 09-28: 0–3 units each);
  Phase 2 and consolidation not started.
- Next when the user resumes: queue order as above (B08 resume first, then new batches, then finishers).

### Earlier state (2026-10-08 11:48 UTC)
- **Running:** B07 resume5 (cloud session_01XittAhqmUUdF5idpt2X6rU, launched 10-08 11:48 on the user's word
  "продолжаем B07"): audit round 3 (A3-8…A3-14, A3-17…A3-19) → further rounds → sweep → final audit → heavy → report.
  Everything else stays paused (B08 included).

### Earlier state (2026-10-02 19:40 UTC) — EVERYTHING PAUSED by user (save weekly tokens)
- **Paused by user (10-02 ~19:30):** B07 (cloud session_01AJjD6znMPKvtM3iddFir9r, idle, not archived). Stopped cleanly after WP3-9:
  last main commit cfe64a7 (F47); WIP branch -wip-B07 = cfe64a7 (clean, nothing unpushed). Audit round 3: 10 items open
  (A3-8…A3-14, A3-17…A3-19, unticked in B07.md worklist); then full fast tier, further audit rounds, sweep, final audit, heavy, report.
  Resume: tools/review/prompts/B07.md (local) or the cloud resume prompt; state from B07.md worklist + git log. Check-in trigger disabled.
- **2026-10-08:** campaign work so far (B01–B06, B24 done; B07 partial; B08 2/6) merged to master via
  AlexanderV/Seqeron#12 (merge commit 7ee2db6) after a CI fix (adbfcaf: B07 heavy tests aligned with the
  documented contracts, skills catalog regenerated). `claude/stoic-maxwell-0olr2z` and `-wip-B07` restarted
  from master 7ee2db6 (same names); the next PR carries only new work.
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
