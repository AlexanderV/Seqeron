---
type: concept
title: "Asynchronous k-mer counting (cancelable, progress-reporting count wrapper)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/KMER-ASYNC-001-Evidence.md
  - docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md
source_commit: b88d1caa7c9b4f3153e097b3e3055db57c3b3551
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: kmer-async-001-evidence
      evidence: "Test Unit ID: KMER-ASYNC-001; Algorithm: Asynchronous K-mer Counting (docs/algorithms/K-mer/Asynchronous_K-mer_Counting.md)"
      confidence: high
      status: current
---

# Asynchronous k-mer counting (cancelable, progress-reporting count wrapper)

The **K-mer** family's thread-pool-offloaded count operation: `KmerAnalyzer.CountKmersAsync`
returns a `Task<Dictionary<string,int>>` mapping each length-`k` substring to its occurrence
count, computed **off the calling thread** so a long input does not block the caller. Its
distinguishing content is *execution*, not arithmetic: the numeric result is defined entirely
by the k-mer count formula and is **byte-for-byte identical** to the synchronous reference
`CountKmers` (test unit KMER-COUNT-001, not yet ingested); what this unit adds and validates is
the **.NET cooperative-cancellation + progress contract** layered over that count. Validated
under test unit **KMER-ASYNC-001** (record [[kmer-async-001-evidence]]);
[[test-unit-registry]] tracks the unit and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## Core model (the count is inherited)

For a sequence `S` of length `L` and window length `k` with `1 ≤ k ≤ L`, the k-mers are the
multiset `{ S[i..i+k−1] : 0 ≤ i ≤ L − k }` — exactly **`L − k + 1`** k-mers (Wikipedia: "a
sequence of length L will have L − k + 1 k-mers and there exist nᵏ total possible k-mers",
n = alphabet size, 4 for DNA). The count of a k-mer `w` is the number of start positions
yielding `w`. When `k > L` or `L = 0` the multiset is empty. This is the *same* count contract
the synchronous method must reproduce; the async method reproduces it exactly because it
delegates to that method.

## The async contract (.NET task cancellation)

The distinctive, separately validated behaviour comes from the .NET task-cancellation model
(Microsoft Learn — *Task Cancellation* / *Task.Run*):

- **Cooperative cancellation.** The synchronous delegate observes the token by calling
  `ThrowIfCancellationRequested` (every `CheckInterval = 1000` windows, to bound polling
  overhead), which throws `OperationCanceledException` and transitions the task to the
  **Canceled** state.
- **Pre-start cancellation.** The token is passed to `Task.Run(func, token)`, so a token
  **already signaled** when the call is made cancels the work *before the delegate runs*.
- **Observing cancellation.** Awaiting a canceled task throws `OperationCanceledException`
  (`TaskCanceledException`, which derives from it, is the wrapped form observed when unwrapped
  from the `AggregateException`).
- **Progress.** An optional `IProgress<double>` reporter is driven 0.0 → 1.0, with the final
  report `1.0` on completion.

## Properties and invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `await CountKmersAsync(S, k)` equals `CountKmers(S, k)` (same keys and counts) | The async method runs the identical synchronous algorithm via `Task.Run`; counts depend only on the k-mer definition |
| INV-02 | Σ counts = `L − k + 1` for `1 ≤ k ≤ L` | Sliding window yields one k-mer per start position |
| INV-03 | A signaled cancellation token ⇒ awaiting the task throws `OperationCanceledException` | Cooperative cancellation: `ThrowIfCancellationRequested` / pre-start `Task.Run` cancellation |
| INV-04 | Empty/null `S` or `k > L` ⇒ empty dictionary; `k ≤ 0` ⇒ `ArgumentOutOfRangeException` (surfaced through the awaited task) | k-mer multiset empty when `L − k + 1 ≤ 0`; synchronous validation preserved through the wrapper |

## Contract and edge cases

- **Inputs:** `sequence` (string; null/empty → empty result; uppercased via `ToUpperInvariant`
  so counting is **case-insensitive**), `k` (int, `> 0`; `k > L` → empty result),
  `cancellationToken` (default `CancellationToken.None`), optional `progress` (`IProgress<double>?`).
- **Alphabet-agnostic:** the alphabet is *not* restricted — non-ACGT characters (e.g. IUPAC `N`)
  are counted literally. Indexing is 0-based; windows are inclusive of length `k`.
- **Edge cases:** empty/null sequence → empty dictionary; `k > L` → empty dictionary; `k ≤ 0` →
  `ArgumentOutOfRangeException` via the awaited task; token signaled before the call *or* during
  the run → awaiting throws `OperationCanceledException`; lowercase/mixed case → same counts as
  uppercase.

## Worked oracles (from the Wikipedia K-mer datasets)

- `CountKmersAsync("ATGG", 3)` → `{ "ATG":1, "TGG":1 }` (Wikipedia: "ATGG has two 3-mers: ATG and TGG"; 4 − 3 + 1 = 2).
- `GTAGAGCTGT` (`L = 10`), total = `L − k + 1` / distinct: **k=2 → 9 / 7** (GT and AG each occur
  twice), **k=3 → 8 / 8** (all windows distinct), **k=4 → 7 / 7** (all windows distinct). The
  async result must equal the synchronous reference on each.

## Deviations and assumptions

- **ASSUMPTION (source-backed, non-correctness-affecting): the numeric contract is identical to
  the synchronous method.** No authoritative source defines a *different* numeric result for an
  "async k-mer count" — the counts are fixed by the `L − k + 1` formula. So the async output MUST
  equal `CountKmers` on identical input; the async layer only governs *how* the work runs, not
  *what* it computes. Verified by tests asserting `CountKmersAsync(seq,k) == CountKmers(seq,k)`.
- **Not implemented:** parallel partitioning of the window scan across cores. The wrapper offloads
  to a *single* thread-pool thread, so it is **not** asymptotically faster than the synchronous
  method (O(n·k) either way); its purpose is non-blocking execution with cancellation and progress.
- **Search-infrastructure decision:** the repository suffix tree was evaluated and **not** used —
  full-spectrum counting is a single linear sliding-window pass, not a fixed-pattern occurrence
  query, so a suffix tree would add O(n) construction without benefit.

## Relation to other units

- **Wrapper of the synchronous count** `KmerAnalyzer.CountKmers` (KMER-COUNT-001, `docs/algorithms/K-mer/K-mer_Counting.md`; not yet ingested) — the async method delegates verbatim to the synchronous cancellation/progress overload. The `L − k + 1` count contract is shared; a future sync-count concept will own the count definition, and this page will link to it.
- The same k-mer sliding-window primitive underlies the K-mer family's other operations (frequency analysis, generation, statistics, both-strand counting) and downstream [[de-bruijn-graph-assembly]] and [[kmer-spectrum-error-correction]].
