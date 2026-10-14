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
- **Model — claude-opus-5-5 for EVERY session (user 2026-10-14):** these are algorithms; output from another model
  does real damage, so the model is checked, never assumed. `create_session`'s `model` is not always honoured.
  Verify right after launch AND at every check-in with `get_session`: BOTH `session_context.model` and
  `external_metadata.last_served_model` must read `claude-opus-5-5` (a turn-level fallback shows up only in the
  latter). If 5-5 cannot be served: **let the session run as it is, but its batch is NOT finished** — do not add it
  to "Finished batches", and record it under "Model-flagged batches" below with what was produced on which model, so
  it can be re-verified when 5-5 becomes available again. This applies to every session of the campaign, including
  the orchestrator's own and any future resume/finisher, not only to the session where it was first noticed.

## State (2026-10-14 10:37 UTC) — FIN-B02 + FIN-B03 + B08 RUNNING (3 in parallel, user's word)
- User 10-14 ~10:30 "Запусти ще B03" + "І резумни b08" (explicitly above the usual max-2 rule):
  - FIN-B03 session_019yqGyk8v1PunaUYAP83i7B (re-audit B03.md + requests to B03, ~22 known leftover items incl. B02
    ProteinPhysicochemistry routing and R16 codon splitter, MCP R-items; F from F22; WIP `-wip-B03`).
  - B08 resume session_01E7mHtXGjdaPccegYRGjrba (CRISPR-PAM/GUIDE done F1–F8; remaining CRISPR-OFF, RESTR-FIND,
    RESTR-DIGEST, RESTR-FILTER, then audit loop/dedup/heavy; F from F9; WIP `-wip-B08`; adds the RestrictionAnalyzer
    raw-pattern API requested by B02 — the CodonOptimizer caller switch stays with B02). Old paused B08 session
    session_01NvVQb3FM1E4YGuhxEtXqjH archived (superseded; it cost ≈ $28.9).
  - **Model rule (user 10-14): every campaign session MUST run claude-opus-5-5** ("це алгоритми; шкода від кривих
    значна"). The first resume session_01E7mHtXGjdaPccegYRGjrba was served claude-opus-5 (its stored session model
    switched at creation although configured 5-5) → INTERRUPTED at 10:53 after ≈ $10, CRISPR-OFF-001 work kept only on
    `-wip-B08` a80bdfa7. Relaunched as **B08 resume-2 session_01JEZHM6XkDpcN8pugzHb9ns** (verified: session model and
    last served model = claude-opus-5-5), whose prompt cherry-picks a80bdfa7 and RE-VERIFIES it against CRISPOR/Hsu/
    Doench before committing. Every check-in now also verifies each session's model and interrupts on a fallback.
  - All three share the branch: no cross-edits of each other's files, `pull --rebase` before every push.
- FIN-B02 at 10:29: round 1 16/19 done (F31+), getorf modes in progress.
- Check-in trigger covers all three (~45 min).

### Earlier state (2026-10-14 07:25 UTC) — FIN-B02 RUNNING (only it)
- User 10-14 "давай допрацюємо B02": session_016tDCWJHjR3NXKjhHLFynHz = FIN-B02 (re-audit B02.md + cross-batch requests to
  B02, completeness-audit loop over 14 known leftover items incl. CODON-RARE-001 Sherlocc `<` vs `≤`, dedup sweep, final
  audit, heavy tier once, `## Leftovers`). F-entries from F31. WIP branch `claude/stoic-maxwell-0olr2z-wip-B02`.
  Check-in trigger re-enabled (~45 min). B08 still paused (not archived); nothing else running.

### Earlier state (2026-10-10 23:40 UTC) — FIN-B24 DONE; EVERYTHING PAUSED
- FIN-B24 done (session_01SLG1sxNgXVQkrwHr1EHya2, 10-09 20:44 – 10-10 23:01, archived, cost ≈ $204.6): audit rounds 1–3 +
  final round, F24–F66 (43 fixes), dedup sweep, heavy tier once (39 new tests, B24 heavy classes 930/930 ×15), full
  unfiltered `Seqeron.Genomics.Tests` 25159/0 (84cbf060). LIMITED ONCO-CNA-002 / ONCO-PURITY-001 → FIXED. Report tip
  03a06ee5; WIP branch `-wip-B24` clean (= main).
- B24.md `## Leftovers`: DOABLE none. BLOCKED with proof: built-in Battenberg IMPUTE2/Beagle5 imputation (binary-only +
  multi-GB 1000G panel, proxy 403; phased path F40 ported), LICHeE post-`fixNetwork` identity-hash group order
  (lichee.jar nondeterministic, 156/4000), MATLAB `unique` tie semantics (MATLAB proprietary; Octave-locked).
  Open cross-batch requests → B17, B19, B21, B22, B23, B25, Phase 2 (none finished — no register entries).
- FIN-B24 also fixed B04's `Bbduk_MonotoneInCutoff` generator (w > k) test-only under the FINISHED-batch rule (84cbf060,
  R2 below).
- Nothing running; B08 paused (not archived); check-in trigger disabled.
- Billing note (user question 10-12): cloud sessions always bill the claude.ai subscription (an API key in the cloud env is
  ignored). API (Console) money can pay only for a LOCAL run — LOCAL_RUN.md + `tools/review/run-batch.sh` with
  `ANTHROPIC_API_KEY` set (or `/login` with a Console account). Reported session costs are API-rate equivalents.
  Moving to a local run was offered; no decision yet.

### Earlier state (2026-10-09 20:45 UTC) — FIN-B24 RUNNING (only it)
- User 10-09 ~20:40 "Давай доробимо B24": session_01SLG1sxNgXVQkrwHr1EHya2 = FIN-B24 (completeness-audit loop over the
  old-process B24 report → all DOABLE incl. LIMITED ONCO-CNA-002 / ONCO-PURITY-001 and the LIMITATIONS proposals →
  dedup sweep → final audit → heavy tier once → `## Leftovers`). WIP branch `claude/stoic-maxwell-0olr2z-wip-B24`.
- Orchestrator permission: FIN-B24 may edit ONLY the ONCO-PURITY-001 methods (`EstimatePurity*`) in B22-owned
  `OncologyAnalyzer.SomaticCalling.cs` (B24 F7/F8 patches), logged in B22.md "## Changes made by FIN-B24 in this file";
  every other B22 item stays a request. B25's ONCO-CHIP-001 methods in Clonality.cs are excluded.
  **User approved this exception 2026-10-10 ~13:20 UTC** ("цього разу можеш правити. Це критична функціональність").
  Done so far: WP1 4f6fd71c (F7/F8 → FIN-B24 F24/F25), logged in B22.md "## Changes made by FIN-B24 in this file".
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

## Model-flagged batches (work produced on a model other than claude-opus-5-5 — NOT finished; re-verify when 5-5 is back)
Per the model rule above. A batch listed here stays out of "Finished batches" even when its report ends with a clean
`## Leftovers`; what it needs is a re-verification pass on claude-opus-5-5, not more work.
| Batch | What was produced off-model | Session | Model actually served | Status |
|---|---|---|---|---|
| B08 | CRISPR-OFF-001 start: F9 MIT silent-100 aggregate fix + CFD aggregate + dedup, WIP only (`-wip-B08` a80bdfa7) | session_01E7mHtXGjdaPccegYRGjrba (interrupted 10:53, ≈ $10, archived) | claude-opus-5 (session model 5-5, served 5) | being re-verified against CRISPOR/Hsu/Doench by B08 resume-2 session_01JEZHM6XkDpcN8pugzHb9ns on verified 5-5 — clears if it holds |
| B08 | CRISPR-OFF-001 work, still WIP-only (`-wip-B08` cda5b4c8) | session_01JEZHM6XkDpcN8pugzHb9ns (resume-2; verified 5-5 at launch 10:54, fell back mid-session; stopped by the five-hour limit ≈ 12:05 after ≈ $199.3; archived) | claude-opus-5 | OPEN — re-verified from scratch by resume-3 session_01WLgw5gNJSSJdHmg5UKDkPY |
| B08 | CRISPR-OFF-001 as re-verified and committed by resume-3, and everything it does next | session_01WLgw5gNJSSJdHmg5UKDkPY (resume-3, RUNNING) | claude-opus-5 — off-model from its first turn, the fourth B08 session in a row | OPEN — left running per the model rule; B08 cannot be marked finished |
| B08 | CRISPR-PAM-001 F1–F2 (commit 73fe1972) + the B08.md skeleton | session_01YRNyThk8LUUghgQswiDmsF (first pass, archived, ≈ $19.6) | claude-opus-5 (session model AND last served) | OPEN — needs a re-verification pass on 5-5 |
| B08 | CRISPR-GUIDE-001 F3–F8 (commits 2fc133c2, e8e6c72) | session_01NvVQb3FM1E4YGuhxEtXqjH (2026-10-01, archived) | claude-opus-5 (configured + session model 5-5, served 5) | OPEN — needs a re-verification pass on 5-5 before B08 may be called finished |

**Provenance audit done 2026-10-14** (all 36 campaign sessions queried; table in `scratchpad/model-provenance.md`):
every session was *configured* claude-opus-5-5, and every batch except B08 was also *served* 5-5 — B01–B07, B09–B26,
B24, the B07 resumes, FIN-B01, FIN-B24 and R1 all clean, so the six entries in "Finished batches" stand. Only B08 is
affected — and it keeps being affected: all five B08 sessions (first pass, 10-01 resume, today's interrupted resume,
resume-2 and resume-3) were served claude-opus-5 — resume-2 after passing the 5-5 check at launch, resume-3 from its
first turn. Relaunching does not fix it, so B08 is simply run as is and flagged. So every B08 fix, F1–F8 and the
CRISPR-OFF work in flight, is off-model, as the rows above record. The one record
with no data is the FIN-B01 launch that was aborted over the literal `$(cat …)` prompt
(session_01ESJRAKGAm3Af18dZknYwM9) — it produced no commits, so nothing depends on it.

## Finished batches (rule "Requests to a FINISHED batch", README → Ownership)
B01 (FIN-B01 2026-10-09), B04, B05, B06, B07, B24 (FIN-B24 2026-10-10). Add a batch here when its final report with
`## Leftovers` is pushed **and** every session that produced its work was served claude-opus-5-5 (see the model rule
and "Model-flagged batches" above); the batches already listed predate that check and the provenance audit may yet
move one of them.

## Requests to finished batches (register — orchestrator checks it on every check-in)
| # | Target | File:line / test | Evidence | From | Status |
|---|---|---|---|---|---|
| R1 | B07 | `B07PrimerDesignMetamorphicTests.LibraryMispriming_OrderInvariant_MonotoneInWeight` | FsCheck counterexample in the FIN-B01 full heavy run (passed on re-run) — root-cause in `PrimerDesigner` library mispriming or the generator | FIN-B01 | done (bd78006) |
| R2 | B04 | `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` | generator drew w = k, violating BBTools `EntropyTracker` k < windowBases (FIN-B24 full heavy run) | FIN-B24 | done (84cbf060, test-only) |

Pre-rule candidates (recorded before 2026-10-09 as "open" in reports; VERIFY each against current code before
acting — several may already be resolved): B06.md → B04 `SequenceComplexity.ShannonEntropyBits` visibility (B01 F-entry
may have resolved via `StatisticsHelper`), ~~B04 `Bbduk_MonotoneInCutoff` generator (w = k = 5)~~ done R2; B06.md → B05
`MotifFinder` Core `CountKmersSpan` callers → `KmerAnalyzer.CountKmers`, B05.md D3 stale rationale, `FindExactMotif`
vs B09 `FindMotif`; B05.md → B04 `RepeatFinder.EnumerateForwardMaximalPairs` → `FindMaximalRepeatedPairs`; B05.md → B07
`ProbeDesigner.FindApproximateMatches` → `ApproximateMatcher.FindWithMismatches`; B04.md → B06
`KmerAnalyzer.CalculateKmerEntropy` delegation + `kmer_entropy` MCP doc example; B01.md → B06 `K-mer_Counting.md` stale
sentence, B05 D3 rationale.

## Known cross-batch items to carry into Phase 2
- ~~B04: `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` generator can draw w = k = 5 → draw w from [k+1, 60] (flaky ~1/30).~~ Done R2 (84cbf060).
- B24.md "## Leftovers" → open requests to B17/B19/B21/B22/B23/B25 + Phase 2 (Median, MAD, KMeans advisory); orchestrator notes for consolidation: ALGORITHMS_CHECKLIST_V2 rows (ONCO-ASCAT-001 method table, l.4784 ONCO-CCF-001 Ckmeans, l.4808 ONCO-HETERO-001 citation), stale Oncology wiki pages, non-B24 Oncology docs linking the pre-split `OncologyAnalyzer.cs`.
- B06 report "Cross-batch requests" (B01, B04, B05, B07, B09, B10, B17, B18, MCP owner) — k-mer counting routed to canonical KmerAnalyzer.
- Stray branch `claude/stoic-maxwell-0olr2z-wip-TEST` (proxy cannot delete it) — ignore or delete in the GitHub UI.

## Log
- 10-02 06:54–07:00 README: no-lost-work rule (90d60c7), WIP branch (a0272b2, de7d295); B07 relaunched with them.
- 10-02 ~11:30 README: targeted test tier (9dc5a87); local-run kit added (tools/review/, LOCAL_RUN.md).
