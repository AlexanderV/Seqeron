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

## State (2026-10-09 20:45 UTC) — FIN-B24 RUNNING (only it)
- User 10-09 ~20:40 "Давай доробимо B24": session_01SLG1sxNgXVQkrwHr1EHya2 = FIN-B24 (completeness-audit loop over the
  old-process B24 report → all DOABLE incl. LIMITED ONCO-CNA-002 / ONCO-PURITY-001 and the LIMITATIONS proposals →
  dedup sweep → final audit → heavy tier once → `## Leftovers`). WIP branch `claude/stoic-maxwell-0olr2z-wip-B24`.
- Orchestrator permission: FIN-B24 may edit ONLY the ONCO-PURITY-001 methods (`EstimatePurity*`) in B22-owned
  `OncologyAnalyzer.SomaticCalling.cs` (B24 F7/F8 patches), logged in B22.md "## Changes made by FIN-B24 in this file";
  every other B22 item stays a request. B25's ONCO-CHIP-001 methods in Clonality.cs are excluded.
- Everything else paused; B08 paused (not archived); check-in trigger re-armed (FIN-B24 check-in).

### Earlier state (2026-10-09 20:25 UTC) — R1 DONE; EVERYTHING PAUSED by user
- R1 done (session_01N4ZfEhLRF1Fh8ZLet3KyaQ, 20:04–20:20, archived, cost ≈ $1.4): test defect, production code unchanged
  (bd780061 + 233381ac). Primer3 `oligo_repeat_library_mispriming` keeps the first entry above the INTEGER running max, so
  only ⌊score⌋ is order-invariant/monotone; the port matches primer3-py 2.3.1 bit-exactly. Property now asserts the integer
  part; regression locks the primer3-py values. Recorded in B07.md "## Post-completion changes".
- Nothing running; B08 paused (not archived); check-in trigger disabled.

### Earlier state (2026-10-09 20:05 UTC) — R1 (B07 post-completion) RUNNING (only it)
- User 10-09 ~20:00 "Запускай": session_01N4ZfEhLRF1Fh8ZLet3KyaQ root-causes register item R1 (B07 FsCheck counterexample
  in LibraryMispriming_OrderInvariant_MonotoneInWeight) under the new "Requests to a FINISHED batch" rule. Nothing else
  runs; after it everything stops again. Check-in trigger re-enabled (+30 min).

### Earlier state (2026-10-09 19:40 UTC) — FIN-B01 DONE; EVERYTHING PAUSED by user
- **FIN-B01 done** (session_018dhx24UekZpd2HKsvqhw3X, 15:57–19:19, archived 19:36, cost ≈ $39.4): audit rounds 1–3 +
  final round (20 DOABLE items), fixes F14–F31, dedup sweep c019eae4, heavy once f2e48444 (+29 tests), final report
  0ca2e247. `## Leftovers`: B01 none (all LIMITATIONS proposals resolved). Full Seqeron.Genomics.Tests 24544/1 — the 1 is
  B07 `B07PrimerDesignMetamorphicTests.LibraryMispriming_OrderInvariant_MonotoneInWeight` (seed-dependent, passes on re-run);
  Mcp.Sequence 112/0, Mcp.Analysis 377/0. WIP branch -wip-B01 clean (= 0ca2e247).
- Orchestrator fix in B01.md: "B07 ProbeDesigner Wallace count loops" was a false positive (Wallace Tm already routed to
  ThermoConstants.CalculateBasicTm by B07; :633 C-vs-G TaqMan rule and :741 G-run rule are not GC counting) → moved to done.
- Everything paused again: nothing launched. B08 paused (not archived). Trigger trig_01MvAagyZjUaqveNm9GXWGtm disabled.
- Batch status: fully done B01 B04 B05 B06 B07; finisher needed B02 (1 open item) B03 B24 (2 LIMITED); B08 2/6 paused;
  B09–B23, B25, B26 not run under the full process; Phase 2 and consolidation not started.

### Earlier state (2026-10-09 15:58 UTC) — FIN-B01 RUNNING (only it)
- **User 10-09 ~15:50: "А тепер давай доробимо B01"** → FIN-B01 launched as cloud session_018dhx24UekZpd2HKsvqhw3X
  (prompt: B01 finisher — re-audit B01.md, completeness-audit loop, all DOABLE incl. LIMITATIONS proposals and
  cross-batch requests addressed to B01 (B06 CountKmersSpan, B03 R18, StatisticsHelper.ShannonEntropy), dedup sweep,
  final audit, heavy once, final B01.md with `## Leftovers`). WIP branch `claude/stoic-maxwell-0olr2z-wip-B01`.
  (A first launch, session_01ESJRAKGAm3Af18dZknYwM9, got an empty prompt and was archived at once.)
- Only FIN-B01 runs; after it finishes everything stops again (launch nothing). B08 stays paused, not archived.
  Check-in trigger trig_01MvAagyZjUaqveNm9GXWGtm re-enabled (every ~45 min, disables itself after FIN-B01).

### Earlier state (2026-10-09 03:55 UTC) — B07 DONE; EVERYTHING PAUSED by user
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

## Finished batches (rule "Requests to a FINISHED batch", README → Ownership)
B01 (FIN-B01 2026-10-09), B04, B05, B06, B07. Add a batch here when its final report with `## Leftovers` is pushed.

## Requests to finished batches (register — orchestrator checks it on every check-in)
| # | Target | File:line / test | Evidence | From | Status |
|---|---|---|---|---|---|
| R1 | B07 | `B07PrimerDesignMetamorphicTests.LibraryMispriming_OrderInvariant_MonotoneInWeight` | FsCheck counterexample in the FIN-B01 full heavy run (passed on re-run) — root-cause in `PrimerDesigner` library mispriming or the generator | FIN-B01 | done (bd78006) |

Pre-rule candidates (recorded before 2026-10-09 as "open" in reports; VERIFY each against current code before
acting — several may already be resolved): B06.md → B04 `SequenceComplexity.ShannonEntropyBits` visibility (B01 F-entry
may have resolved via `StatisticsHelper`), B04 `Bbduk_MonotoneInCutoff` generator (w = k = 5); B06.md → B05
`MotifFinder` Core `CountKmersSpan` callers → `KmerAnalyzer.CountKmers`, B05.md D3 stale rationale, `FindExactMotif`
vs B09 `FindMotif`; B05.md → B04 `RepeatFinder.EnumerateForwardMaximalPairs` → `FindMaximalRepeatedPairs`; B05.md → B07
`ProbeDesigner.FindApproximateMatches` → `ApproximateMatcher.FindWithMismatches`; B04.md → B06
`KmerAnalyzer.CalculateKmerEntropy` delegation + `kmer_entropy` MCP doc example; B01.md → B06 `K-mer_Counting.md` stale
sentence, B05 D3 rationale.

## Known cross-batch items to carry into Phase 2
- B04: `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` generator can draw w = k = 5 → draw w from [k+1, 60] (flaky ~1/30).
- B06 report "Cross-batch requests" (B01, B04, B05, B07, B09, B10, B17, B18, MCP owner) — k-mer counting routed to canonical KmerAnalyzer.
- Stray branch `claude/stoic-maxwell-0olr2z-wip-TEST` (proxy cannot delete it) — ignore or delete in the GitHub UI.

## Log
- 10-02 06:54–07:00 README: no-lost-work rule (90d60c7), WIP branch (a0272b2, de7d295); B07 relaunched with them.
- 10-02 ~11:30 README: targeted test tier (9dc5a87); local-run kit added (tools/review/, LOCAL_RUN.md).
