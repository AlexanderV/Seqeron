---
type: source
title: "Evidence: KMER-ASYNC-001 (Asynchronous k-mer counting — cancelable count wrapper)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-ASYNC-001-Evidence.md
sources:
  - docs/Evidence/KMER-ASYNC-001-Evidence.md
source_commit: b88d1caa7c9b4f3153e097b3e3055db57c3b3551
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-ASYNC-001

The validation-evidence artifact for test unit **KMER-ASYNC-001** — **asynchronous k-mer
counting** (`KmerAnalyzer.CountKmersAsync`), the cooperatively cancelable, progress-reporting
thread-pool wrapper over the synchronous count. This is the first **K-mer** family Evidence file
and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the operation itself is synthesized in [[asynchronous-kmer-counting]]. See
[[test-unit-registry]] for how units are tracked.

It is **not a distinct counting algorithm**: the file's central assumption is that the async
method's numeric result is *identical* to the synchronous reference `CountKmers` (KMER-COUNT-001)
because k-mer counts are fixed by the formula `L − k + 1`. What KMER-ASYNC-001 uniquely validates is
the **.NET cooperative-cancellation and progress contract** layered on top of that count.

## What this file records

- **Online sources:**
  - **Wikipedia — K-mer** (authority rank 4, for the combinatorial definition) — the count
    contract "a sequence of length L will have **L − k + 1** k-mers and there exist **nᵏ** total
    possible k-mers" (n = alphabet size, 4 for DNA), the worked example "**ATGG** has two 3-mers:
    ATG and TGG", and the `GTAGAGCTGT` (L = 10) per-k total table. Corner cases: `k > L` ⇒
    `L − k + 1 ≤ 0` ⇒ no k-mers; empty sequence (`L = 0`) ⇒ no k-mers.
  - **Microsoft Learn — Task Cancellation (.NET)** (authority rank 2) — cooperative cancellation
    ("the requesting code calling `CancellationTokenSource.Cancel` and the user delegate
    terminating in a timely manner"), acknowledging cancellation by throwing
    `OperationCanceledException` via **`ThrowIfCancellationRequested`** (task → **Canceled**
    state), and that awaiting a canceled task throws `TaskCanceledException` (wrapped in
    `AggregateException`; unwrapped on `await`; derives from `OperationCanceledException`).
  - **Microsoft Learn — Task.Run Method** (authority rank 2) — **pre-start cancellation**: "A
    cancellation token allows the work to be cancelled if it has not yet started" (an
    already-signaled token to `Task.Run(func, token)` yields a Canceled task without running the
    delegate), and that the `OperationCanceledException` "is stored into the returned task".

- **Documented corner cases / failure modes:** `k > L` and empty/`L = 0` sequence ⇒ empty result;
  token signaled **before** start ⇒ delegate may not run, task Canceled, awaiting throws
  `OperationCanceledException`; token signaled **during** run ⇒ `ThrowIfCancellationRequested`
  transitions the task to Canceled, awaiting throws `OperationCanceledException`.

- **Datasets (documented oracles):**
  - `GTAGAGCTGT` (L = 10) — k=2 total/distinct **9 / 7** (GT and AG each occur twice); k=3 **8 / 8**;
    k=4 **7 / 7**. The async result must equal the synchronous reference on each.
  - `ATGG`, k = 3 — 3-mers `ATG` (×1) and `TGG` (×1).

- **Test-coverage recommendations:** MUST — async result equals the synchronous reference for the
  `GTAGAGCTGT` (k=2,3,4) and `ATGG` datasets; total-count invariant Σcounts = `L − k + 1`;
  empty/null and `k > L` ⇒ empty dictionary; `k ≤ 0` ⇒ `ArgumentOutOfRangeException` via the awaited
  task; a token signaled before start and one signaled mid-run on a large input ⇒ awaiting throws
  `OperationCanceledException`. SHOULD — large input (≥ `CheckInterval`) returns the synchronous
  result and reports progress to 1.0. COULD — `CancellationToken.None` default and
  case-insensitivity carried through the async path.

## Deviations and assumptions

**One ASSUMPTION**, source-backed and non-correctness-affecting for numeric output:

- **The numeric-correctness contract is identical to the synchronous method.** `CountKmersAsync`
  is a `Task.Run` parallel/async wrapper over the validated synchronous
  `CountKmers(string, int, CancellationToken, IProgress<double>)`; no authoritative source defines
  a different numeric result for an "async k-mer count" (counts are defined purely by `L − k + 1`).
  The async result MUST equal the synchronous reference on identical input — verified by tests
  asserting `CountKmersAsync(seq,k) == CountKmers(seq,k)`. The async layer governs *execution*
  (thread-pool offload, cooperative cancellation, progress), not the count.

No source contradictions — the Wikipedia count definition and the .NET cancellation/`Task.Run`
contracts are orthogonal (one fixes *what* is computed, the others *how* the work runs) and
mutually consistent.
