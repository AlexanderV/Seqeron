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
- **Decide before spending (user 2026-10-14, after the B08 restarts):** before any interrupt, relaunch or new session
  that costs real money, say the cost and the expected gain first, and do not retry while the cause of the failure is
  still true (e.g. the model fallback is account-wide at that minute — then report and wait, do not relaunch). A batch
  already flagged unfinishable gets no further spend without the user's decision. Never create a session for work the
  user may be restarting themselves — ask who does it. Waste this rule exists to prevent: the 10:53 B08 relaunch
  ($10, the fallback was account-wide and the new session fell back too), ~$199 spent on resume-2 without first asking
  whether an off-model B08 was worth full price, and the duplicate FIN-B03 at 14:03 ($0.22).
- **Scope of a model re-verification pass (user 2026-10-14) — this is the SPEC of the pass that is now DEFERRED, not cancelled (see the 5-5-debt rule below): when the debt is eventually paid, this is what the pass does.** Off-model work does NOT have to be redone. What is
  model-independent stays valid — oracle numbers recorded in Evidence (CRISPOR, Biopython Restriction, matrices read
  from the papers), the tests that lock them bit-exactly, and the cited source line or formula. What needs a pass on
  claude-opus-5-5 is the judgement: which convention is canonical where sources disagree, whether a divergence is a
  defect or a documented deviation, whether a BLOCKED proof really proves it, whether an unsourced heuristic may stay.
  So the owed B08 pass — and now B02's as well — is a REVIEW of those decisions against the recorded evidence, not a
  re-run of the batch, and it must check each test expectation against the oracle number in Evidence: green tests
  written off-model can lock in a wrong expectation. Every off-model F-entry must therefore carry its oracle number in
  Evidence. This pass is **owed, not scheduled** — nothing is launched for it now, and it does not block any batch.
- **Model — claude-opus-5-5 for EVERY session (user 2026-10-14):** these are algorithms; output from another model
  does real damage, so the model is checked, never assumed. `create_session`'s `model` is not always honoured.
  Verify right after launch AND at every check-in with `get_session`: BOTH `session_context.model` and
  `external_metadata.last_served_model` must read `claude-opus-5-5` (a turn-level fallback shows up only in the
  latter). If 5-5 cannot be served: **let the session run as it is** and record what was produced on which model
  under "5-5 debt" below. This applies to every session of the campaign, including the orchestrator's own and any
  future resume/finisher, not only to the session where it was first noticed.
- **SUPERSEDED — the 5-5 debt is a note, not a blocker (user 2026-10-14 23:15):** "На 5.5 не запускається — так і не
  запуститься. Що ти хочеш від нього? На 5-му перевірив? Ок? Помітив, що борг 5.5 і все." Subagents are served
  `claude-opus-5` and will keep being served it, so withholding completion over the model blocks the campaign on
  something no session can fix. From now on: work verified on claude-opus-5 **counts as verified**; the batch is
  finished when its `## Leftovers` has zero DOABLE items, and **no session is launched for the sole purpose of clearing
  the debt**. But the debt is **not written off** (user 23:30: "Не треба прибирати всі згадки про відкладену перевірку
  на 5-5. Вони потім, наступного разу, не скоро будуть валідуватись. Ми повинні помітити і залишити до наступного
  разу"): every off-model entry stays listed under "5-5 debt" as an **owed, deferred pass**, with the F-entries named,
  so that whenever a model that is actually served can do it, the pass knows exactly what to re-read. The campaign's
  final consolidation carries that register forward instead of closing it. The rule above still stands for what it is
  good for: always check and record which model served a session, and never silently pretend a session ran on 5-5.
  (VER-B02, launched at 23:09 to clear B02's debt, was stopped under this rule after ~6 minutes — the pass it was to
  run is deferred, not cancelled.)

**Subagents fall back off-model even when the lead does not (2026-10-14).** FIN-B02's lead was served claude-opus-5-5
for its whole run while its subagents were served `claude-opus-5` — including both of its own re-verification reviewers.
So the launch-time model check on the lead is NOT sufficient: for any pass whose deliverable *is* the model (a 5-5
re-verification), every judgement must be made in the lead's own context and subagents used only for verdict-free
retrieval (install a tool, run a command, grep, list file:line). And such a pass that is itself served off-model
produces nothing. Under the superseding rule above this no longer blocks anything: it is the reason the 5-5 debt
is simply recorded instead of chased.

## State (2026-10-14 23:20 UTC) — B02 FINISHED (5-5 debt noted); SRC + B08 resume-4 running
- **FIN-B03 DONE** (session_019yqGyk8v1PunaUYAP83i7B archived, 10:35 – 17:57, ≈ $75.1, served claude-opus-5-5
  throughout, as was the B03 first pass — so B03 goes into "Finished batches"): audit rounds 1–3 + final round,
  F22 onwards, dedup sweep, heavy tier once (da6bd804), final report 5bcceeca. 73 changed docs listed for wiki
  re-ingest (orchestrator task).
- B03.md `## Leftovers`: DOABLE none. BLOCKED with proof: ProtScale exponential weight-variation model (no source
  gives the equation; mirrors 403/ENOTFOUND; the one sample fits a ∈ [0.81, 0.94] with no natural constant),
  pI with PTM charges (one-letter API, no reference models PTMs), pI with coupled electrostatics (needs 3-D
  coordinates). Open cross-batch requests: **R26, R27, C3 → B02 (RELAY, B02 is in flight, not the register)** —
  monoisotopic amino-acid table + water 18.010565 and a parameterised pK table into
  `Infrastructure/ProteinPhysicochemistry.cs`, and `SequenceExtensions.ToUpperAscii` made public; A33 → B14,
  A34 → B15, A35 → B11 (targets not finished); QUEUE R3 → B01 and R4 → B07 already in the register below.
- Still running: FIN-B02 (round 3, WP22 landed b2148fa3, WP23 CodonW Fop/CBI next) and B08 resume-3
  (RESTR-DIGEST-001, Biopython cross-check 0 mismatches; 4/6 units on main; off-model, see below).
- **Cost watch:** B08 resume-3 passed the $120 threshold (≈ $125.5 at 18:32, context 257k) without closing its unit.

### Earlier state (2026-10-14 10:37 UTC) — FIN-B02 + FIN-B03 + B08 RUNNING (3 in parallel, user's word)
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

## 5-5 debt (work served by `claude-opus-5` instead of `claude-opus-5-5` — a recorded note, NOT a blocker)
**Campaign-wide (user 2026-10-14 23:15, "І так для всіх взагалі").** This applies to EVERY batch — past, present and
future, B08 included — not to B02 alone. A batch is finished when its `## Leftovers` holds zero DOABLE items; the model
that served the work is **recorded, never required**, and **no session is launched to clear a row here**.

**Every row is an OWED, DEFERRED pass — not a closed one (user 23:30).** The validation will happen later, not soon, so
each row names the F-entries the pass must re-read and the spec it follows ("Scope of a model re-verification pass" in
Standing user rules): a REVIEW of the judgement calls against the oracle numbers already in Evidence, including a check
that every test expectation matches its oracle number. Nothing here is deleted when a batch is marked finished, and the
final consolidation carries this register forward rather than closing it. What remains mandatory meanwhile: check
`session_context.model` and `external_metadata.last_served_model` at launch and at every check-in, and never state that
a session ran on 5-5 when it did not.

| Batch | What was produced off-model | Session | Model actually served | Owed pass (deferred) |
|---|---|---|---|---|
| B08 | CRISPR-PAM-001 F1–F2 (`73fe1972`) + the B08.md skeleton | session_01YRNyThk8LUUghgQswiDmsF (first pass, archived, ≈ $19.6) | claude-opus-5 (session model AND last served) | **OWED** — 5-5 review of F1–F2 against CRISPOR/Doench/Hsu numbers in Evidence |
| B08 | CRISPR-GUIDE-001 F3–F8 (`2fc133c2`, `e8e6c72`) | session_01NvVQb3FM1E4YGuhxEtXqjH (2026-10-01, archived) | claude-opus-5 (configured + session model 5-5, served 5) | **OWED** — 5-5 review of F3–F8 (CRISPR-GUIDE-001 judgement calls + test expectations vs Evidence) |
| B08 | CRISPR-OFF-001 WIP that was never merged (`-wip-B08` a80bdfa7, cda5b4c8) | session_01E7mHtXGjdaPccegYRGjrba (≈ $10) and session_01JEZHM6XkDpcN8pugzHb9ns (resume-2, ≈ $199.3; both archived) | claude-opus-5 (resume-2 passed the 5-5 check at launch, then fell back mid-session) | nothing owed — never merged; resume-3 re-did this work from scratch (its own row carries the debt) |
| B08 | F9–F39: CRISPR-OFF-001 (re-verified claim-by-claim, 5 claims overturned), RESTR-FIND/DIGEST/FILTER-001, audit round 1, WP1/WP2/WP-SRC | session_01WLgw5gNJSSJdHmg5UKDkPY (resume-3, ≈ $206.5, archived) | claude-opus-5 from its first turn, the fourth B08 session in a row | **OWED** — 5-5 review of F9–F39: the overturned CRISPR-OFF claims, the four RESTR units' convention calls and the audit-round-1 BLOCKED proofs. Its oracle runs are on the record (144→0 of 34 906; 313→0 of 1452; 1521/1521 + 75 466/75 466), so the pass reads decisions, not numbers |
| B08 | whatever resume-4 produces | session_01Gei6HyCuq5texhxCCGEvji (RUNNING; 5-5 on its first turn — the first B08 session to get it) | record at each check-in | **OWED only for whatever turns out off-model** — 5-5 on its first turn, so re-check at each check-in; B08 finishes on its leftovers either way |
| B02 | F41, F42, F43, F48, F49, F54, F55, F60, F63, F64 + the WP27a/b review verdicts (B02.md "## Model provenance") | FIN-B02 `session_016tDCWJHjR3NXKjhHLFynHz` — lead 5-5 throughout, its **subagents** served opus-5 for WP9/WP10/WP14/WP15/WP18/WP22 and for BOTH re-verification reviewers (WP27a 114 of 124 turns, WP27b 86 of 106) | claude-opus-5 (runtime fallback; per-turn model fields in the subagent transcripts) | **OWED** — 5-5 review of those ten F-entries and the WP27a/b verdicts. Accepted as verified for completion (B02 is in "Finished batches"): the off-model pass re-ran every oracle (ViennaRNA 2.7.2, dnachisel 3.2.16, Biopython 1.88 across all 27 supported tables, an independent `sherlocc.pl` port) and found two real defects, F63 and F64 — which is also why the deferred pass is worth running |

**Provenance audit done 2026-10-14** (all 36 campaign sessions queried; table in `scratchpad/model-provenance.md`):
every session was *configured* claude-opus-5-5, and every batch except B08 was also *served* 5-5 — B01–B07, B09–B26,
the B07 resumes, FIN-B01, FIN-B24 and R1 all clean. B08 kept being affected across all five of its sessions, and
FIN-B02 showed the subtler form: a 5-5 lead with opus-5 subagents. The one record with no data is the FIN-B01 launch
aborted over the literal `$(cat …)` prompt (session_01ESJRAKGAm3Af18dZknYwM9) — no commits, nothing depends on it.

**Why chasing 5-5 was abandoned.** VER-B02 (session_015916Bdb7VvK2GSYCTiNjbG), launched 23:09 to clear B02's debt,
recorded the decisive evidence itself: it verified 5-5 on its first turn and then fell back to `claude-opus-5`
mid-session, in **both** `session_context.model` and `last_served_model`. So "launch it again and check at launch" does
not hold a session on 5-5 — resume-2 showed the same thing in the morning. It stopped with nothing verified and no code
touched (its one durable product: the oracle environment is pinned and reproducible, B02.md "## Model provenance" →
VER-B02), and was archived. The policy above stops the chasing; it does not retire the pass, which stays owed in the
rows above.

## Consolidation task — re-anchor citations to reachable sources (user 2026-10-14)
The agent proxy answers **403 to CONNECT** for `en.wikipedia.org:443`, `www.addgene.org:443` and
`rebase.neb.com:443`; earlier batches also recorded `rosalind.info`, Kazusa (`www.kazusa.or.jp`), ExPASy
(`web.expasy.org`) and the 1000 Genomes panel as unreachable. No number in the review rests on them — the method is
oracle-first (install the reference tool, compare numbers), and e.g. RESTR-DIGEST's circular `k` / linear `k+1` and
fragment-sum rules were re-confirmed against Biopython `catalyse` on 1452 cases — but two things need work:
1. **275 files** under `docs/Validation/reports/` and `tests/TestSpecs/` cite Wikipedia or Addgene, mostly inherited
   from the pre-campaign pass and no longer re-fetchable here. Where such a citation backs an actual NUMBER
   (e.g. `PRIMER-DESIGN-001`: primer length 18–24, GC 40–60 %, pair Tm ≤ 5 °C partly from Addgene), re-anchor it to a
   reachable source (Primer3 source, EMBOSS, the paper) or mark it "inherited, not re-verified in this environment".
   Batch reports that cite them: B01, B04, B05, B06, B07, B08, README.
2. ~~**REBASE is only available as the snapshot bundled in Biopython** (`Bio.Restriction`), i.e. enzyme definitions as
   of that Biopython release, not live REBASE.~~ **Premise corrected by B08 WP-SRC (resume-3):** REBASE `allenz.104`
   **with its verbatim header** (`REBASE version 104`, `Copyright (c) Dr. Richard J. Roberts, 2021`,
   `Rich Roberts Mar 31 2021`) is served from the allowed host `raw.githubusercontent.com`, vendored in
   `crisprVerse/crisprBase`, and B08 already uses it. What stands as a limitation: **no REBASE release newer than 404
   (Apr 2024) is reachable**, and recognition sites / cut positions still come from the Biopython snapshot (version
   stated in Evidence). Record that, not the old claim.
Hosts to allow in the environment's network settings so later sessions do not hit this: `en.wikipedia.org`,
`www.addgene.org`, `rebase.neb.com`, `rosalind.info`, `www.kazusa.or.jp`, `web.expasy.org`. Kazusa and ExPASy would
also unblock two items currently BLOCKED with proof (B02 E. coli per-thousand table; B03 ProtScale exponential
weighting), so those should be re-opened once access exists. **→ access opened 2026-10-14, see the next section.**

## Network access opened (user 2026-10-14 ~22:00 UTC) — two BLOCKED items re-opened
The user opened the environment's network settings. Verified from the orchestrator container with
`curl -sS -o /dev/null -w '%{http_code}'` **plus** a payload check (not just a status code):
| Host | Result | Payload seen |
|---|---|---|
| `www.kazusa.or.jp` | 200 | Codon Usage Database page for *E. coli* W3110 incl. the per-thousand (`/1000`) columns |
| `web.expasy.org` | 200 | ProtScale page listing "Weight variation" with the **exponential** and **linear** models |
| `en.wikipedia.org` | 200 | article HTML |
| `www.addgene.org` | 200 | protocol HTML |
| `rosalind.info` | 200 | problem HTML |
| `raw.githubusercontent.com`, PyPI | 200 | unchanged, still fine |
| `rebase.neb.com` | **curl exit 60** (TLS certificate verification fails; `--cacert /root/.ccr/ca-bundle.crt` does not help although the proxy reports `bundleCoversEveryHost: true`) | — |

**REBASE decision (user: "брати дані REBASE з дзеркала allenz або ігнорити помилки сертифікату").** We take the
mirror, never the certificate shortcut: **disabling TLS verification and unsetting `HTTPS_PROXY` are forbidden in this
environment** (container rules), so `curl -k` is not an option for any session. The mirror is already in use — B08
sourced the REBASE organism and methylation columns from **REBASE `allenz.104` as vendored in the Bioconductor package
`crisprBase` (`crisprVerse/crisprBase`, reachable over `raw.githubusercontent.com`)**: see
`docs/Evidence/B08-Source-Provenance.md` §3/§3a and `RESTR-FIND-001-Evidence.md:62,122,236`,
`RESTR-FILTER-001-Evidence.md:168`, `RESTR-DIGEST-001-Evidence.md:258`. Rule for every session: enzyme recognition
sites and cut positions come from the Biopython `Bio.Restriction` snapshot (version stated in Evidence); organism /
methylation / prototype columns come from `allenz.104` via `crisprBase`; live `rebase.neb.com` stays unavailable and
that is recorded as an environment limitation, not worked around.

**Re-opened items** (were BLOCKED-with-proof solely because the host answered 403; the proof no longer holds):
| # | Batch | Item | Where it is recorded | Owner of the fix |
|---|---|---|---|---|
| RO1 | B02 | CODON-USAGE-001 organism presets hold per-AA fractions only — add the per-thousand columns from the Kazusa Codon Usage Database (`MolTools/CodonUsageAnalyzer.cs`). **Supersedes B02.md leftover (b)** ("E. coli 316407 per-thousand usage — kazusa.or.jp blocked"), whose proof was the 403 | B02.md `## Leftovers` (b) | **SRC** — FIN-B02 landed at 23:04, so SRC's task C condition is met; commit "(SRC for B02)" |
| RO2 | B03 | PROT-SCALE / ProtScale **exponential** weight-variation model — the equation was unobtainable (mirrors 403/ENOTFOUND; one sample fitted a ∈ [0.81, 0.94]) | B03.md `## Leftovers` (BLOCKED) | SRC session, commit "(SRC for B03)" — B03 is FINISHED |

Both are now DOABLE unless the fetched page still fails to define the quantity; if so, the new proof must quote the
fetched page (not a 403).

Also newly verifiable, but **B08-owned** (handled by the B08 session, not by SRC): the Addgene 8–10 nt seed-range
citation at `src/.../MolTools/CrisprDesigner.cs:258`, which B08 WP-SRC had to leave unverifiable
(nothing depends on it numerically; the comment states its own derivation).

## Running sessions (2026-10-14 23:00 UTC)
| Task | Session | Configured / served | Scope |
|---|---|---|---|
| ~~FIN-B02~~ | `session_016tDCWJHjR3NXKjhHLFynHz` (archived, $188.14) | lead 5-5 / **subagents mostly opus-5** | **DONE 2026-10-14 23:04**, final report `b36eecb3`: 4 audit rounds to 0 DOABLE, F31–F64, dedup sweep, heavy tier, full unfiltered **26280/0**. **FINISHED** (0 DOABLE); its off-model entries are a note in "5-5 debt" |
| SRC (sources & citations) | `session_01Wa9VJs4Gp5H1yEXuySyxf8` | 5-5 / verify at check-in | A: RO2 ProtScale exponential (SRC for B03); B: the citation audit — 428 files, 790 occurrences, **206 distinct URLs** (5 Addgene), classify NUMBER / DEFINITION / BIBLIOGRAPHY, fetch and verify every NUMBER, re-anchor code- and test-driving numbers to primary sources, mark the rest "inherited, not re-verified"; C: RO1 Kazusa per-thousand **only if FIN-B02 has landed**, else fetch + register. Report `SRC.md`, WIP branch `…-wip-SRC` |
| ~~VER-B02~~ | `session_015916Bdb7VvK2GSYCTiNjbG` (stopped 23:15, archived) | 5-5 | Launched to clear B02's 5-5 debt; **stopped by the user's ruling** that the debt is a note and not a blocker. Its design lesson is kept above: subagents fall back off-model even when the lead does not |
| B08 resume-4 | `session_01Gei6HyCuq5texhxCCGEvji` | 5-5 / verify at check-in | the `## Resume state` list in B08.md: D3 → D5+D12 → D14+D17 → D15+D16 → D9 → D10 → D7 → D8 → D11 → D6 → D4 → D18 → D19, then round-1 fast tier → round 2 → dedup → final audit → heavy tier → report. Also the now-reachable Addgene citation at `CrisprDesigner.cs:258` |

**B08's launch prompt predates the 23:15 ruling** — it was told its report must keep B08 out of "Finished batches"
over the model. When it lands, the orchestrator corrects that: B08 goes into "Finished batches" if its `## Leftovers`
has zero DOABLE items, with its off-model entries moved to the "5-5 debt" register. Cloud sessions cannot be messaged,
so this is a check-in task, not a relaunch.

Ownership while these run: SRC must not touch B02-owned or B08-owned files (its prompt names them); B08 keeps its own
`CRISPR-*`/`RESTR-*` citations; nobody force-pushes the shared branch. Two sessions launched at 22:56 without a repo
source (`config:no-git-repo`, my omission) were interrupted after ~20 s and archived — $0.36 total, no commits; the
replacements above carry `source_url` + `source_revision`.

## Finished batches (rule "Requests to a FINISHED batch", README → Ownership)
B01 (FIN-B01 2026-10-09), **B02 (FIN-B02 2026-10-14, 5-5 debt recorded)**, B03 (FIN-B03 2026-10-14), B04, B05, B06,
B07, B24 (FIN-B24 2026-10-10). Add a batch here when its final report with `## Leftovers` is pushed and has zero
DOABLE items. The model is recorded, not required: off-model entries go to the "5-5 debt" register above and do not
hold the batch back (user 2026-10-14 23:15).

## Requests to finished batches (register — orchestrator checks it on every check-in)
| # | Target | File:line / test | Evidence | From | Status |
|---|---|---|---|---|---|
| R1 | B07 | `B07PrimerDesignMetamorphicTests.LibraryMispriming_OrderInvariant_MonotoneInWeight` | FsCheck counterexample in the FIN-B01 full heavy run (passed on re-run) — root-cause in `PrimerDesigner` library mispriming or the generator | FIN-B01 | done (bd78006) |
| R2 | B04 | `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` | generator drew w = k, violating BBTools `EntropyTracker` k < windowBases (FIN-B24 full heavy run) | FIN-B24 | done (84cbf060, test-only) |
| R3 | B01 | `GcSkewCalculator.CalculateWindowedGcContent` (both overloads, GcSkewCalculator.cs:1136/1162) | additive `includePartialWindow` (private `EnumerateWindows` :199 already has it); then SS `CalculateGcContentProfile(…, includePartialWindow)` projects it instead of its own `PartialWindowStarts` copy of the placement (B03.md "Duplication sweep (finisher)") | FIN-B03 | open |
| R4 | B07 | `PrimerProbeProperties.LnaTm_AddingLnaCPosition_NeverLowersTm` (tests/.../Properties/PrimerProbeProperties.cs:2370); `PrimerDesigner.CalculateMeltingTemperatureNNLna` (MolTools/PrimerDesigner.cs:2966) | FsCheck counterexample in the FIN-B03 full heavy run (passed 3/3 on re-run). Brute force over the test's generator: 5/6046 violations, e.g. GCTACGCCCGGGGCTAGTGACAAG, LNA {1,4,7,8,9,10,11,12,14,15,18,19,20,22} + C@13: 94.1352 → 93.9034 °C. The invariant is not exact under the default Owczarzy 2011 model (B07 F27). Check against MELTING 5 whether the model is correct, then narrow the invariant with a source, or fix the model | FIN-B03 | open |
| R6 | B03 | `Analysis/SequenceStatistics.cs` `HydrophobicityProfileOptions.WeightModel` / `HydropathyWeightModel.Exponential` (QUEUE RO2) | ProtScale exponential model derived from the live web tool (34 runs) and locked to its output, `docs/Evidence/SEQ-HYDRO-001-Evidence.md` "SRC RO2" | SRC | done (SRC RO2 commit) |
| R5 | B03 | `Analysis/SequenceStatistics.cs:1184-1191` private `ToUpperAscii(char)` / `ToUpperAscii(string)` | verbatim copy of the new public canonical `Infrastructure/SequenceCase.ToUpperAscii` (FIN-B02 WP24; same body `IsAsciiLetterLower ? c - 32 : c`); caller switch, no behaviour change | FIN-B02 | done (d1e6111) |

Pre-rule candidates (recorded before 2026-10-09 as "open" in reports; VERIFY each against current code before
acting — several may already be resolved): B06.md → B04 `SequenceComplexity.ShannonEntropyBits` visibility (B01 F-entry
may have resolved via `StatisticsHelper`), ~~B04 `Bbduk_MonotoneInCutoff` generator (w = k = 5)~~ done R2; B06.md → B05
`MotifFinder` Core `CountKmersSpan` callers → `KmerAnalyzer.CountKmers`, B05.md D3 stale rationale, `FindExactMotif`
vs B09 `FindMotif`; B05.md → B04 `RepeatFinder.EnumerateForwardMaximalPairs` → `FindMaximalRepeatedPairs`; B05.md → B07
`ProbeDesigner.FindApproximateMatches` → `ApproximateMatcher.FindWithMismatches`; B04.md → B06
`KmerAnalyzer.CalculateKmerEntropy` delegation + `kmer_entropy` MCP doc example; B01.md → B06 `K-mer_Counting.md` stale
sentence, B05 D3 rationale.

## Known cross-batch items to carry into Phase 2
- **The 5-5 debt register must be carried into the final consolidation, not closed.** Every row in "5-5 debt" is an
  owed review pass (spec: Standing user rules → "Scope of a model re-verification pass"). Validation will come later,
  not soon; the final summary lists what is still owed, per batch and per F-entry.
- B02 (FIN-B02 final report `b36eecb3`) → requests to batches not yet reviewed, each to be honoured when that batch runs:
  B11 `GenomeAnnotator.GetCodonUsage` (:1175, :1232) → `Translator.CountCodons` + `CodonUsageAnalyzer.CalculateRscu`;
  B11/B09 ORF-engine unification on `Translator.FindOrfs`; B20 `GenBankParser` private translate (:686-697) →
  `Translator.Translate`; B09 `GenomicAnalyzer` hard-coded stops (:413, :506) → `GeneticCode.Standard.IsStopCodon`;
  B23 `OncologyAnalyzer.Fusions` to-first-stop loop (:546-557) → `Translator.Translate(…, toFirstStop: true)`;
  B12 MFE vs ViennaRNA int21/int22 residual (max |Δ| 2.24 kcal/mol, RNA-STRUCT-001 D3/D4).
  Its three requests to B08 are already covered by B08's own worklist (raw-pattern search done as B08 F22
  `FindPatternOccurrences`; `CrisprDesigner.cs:662` ToUpperAscii = D12; the two MCP doc pairs = D14+D17).
- ~~B04: `B04AuditProperties.Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU` generator can draw w = k = 5 → draw w from [k+1, 60] (flaky ~1/30).~~ Done R2 (84cbf060).
- B24.md "## Leftovers" → open requests to B17/B19/B21/B22/B23/B25 + Phase 2 (Median, MAD, KMeans advisory); orchestrator notes for consolidation: ALGORITHMS_CHECKLIST_V2 rows (ONCO-ASCAT-001 method table, l.4784 ONCO-CCF-001 Ckmeans, l.4808 ONCO-HETERO-001 citation), stale Oncology wiki pages, non-B24 Oncology docs linking the pre-split `OncologyAnalyzer.cs`.
- B06 report "Cross-batch requests" (B01, B04, B05, B07, B09, B10, B17, B18, MCP owner) — k-mer counting routed to canonical KmerAnalyzer.
- Stray branch `claude/stoic-maxwell-0olr2z-wip-TEST` (proxy cannot delete it) — ignore or delete in the GitHub UI.

## Log
- 10-02 06:54–07:00 README: no-lost-work rule (90d60c7), WIP branch (a0272b2, de7d295); B07 relaunched with them.
- 10-02 ~11:30 README: targeted test tier (9dc5a87); local-run kit added (tools/review/, LOCAL_RUN.md).
