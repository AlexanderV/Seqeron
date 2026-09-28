# Asynchronous K-mer Counting

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-ASYNC-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Asynchronous k-mer counting computes, off the calling thread, the number of occurrences
of every length-*k* substring (k-mer) in a sequence. It is the cooperatively cancelable,
progress-reporting variant of synchronous k-mer counting (KMER-COUNT-001): the numeric
result is defined entirely by the k-mer definition [1] and is identical to the
synchronous reference, while the operation runs on the thread pool so long inputs do not
block the caller. The computation is exact (not heuristic): for a sequence of length *L*
it produces exactly *L − k + 1* k-mer occurrences [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a contiguous substring of length *k* over the sequence alphabet (DNA: A,C,G,T).
Counting k-mers underpins assembly, alignment-free comparison, and repeat detection. The
"asynchronous" qualifier concerns *execution* (thread-pool offload + cancellation), not
the mathematical definition of the count.

### 2.2 Core Model

For a sequence *S* of length *L* and window length *k* with 1 ≤ *k* ≤ *L*, the multiset of
k-mers is { S[i .. i+k−1] : 0 ≤ i ≤ L − k }, giving exactly *L − k + 1* k-mers [1]. The
count of a k-mer *w* is the number of start positions *i* yielding *w*. When *k* > *L* or
*L* = 0 the multiset is empty [1].

The asynchronous contract is governed by the .NET task-cancellation model [2][3]:
cancellation is cooperative; the delegate observes the token via
`ThrowIfCancellationRequested`, which throws `OperationCanceledException` and transitions
the task to the Canceled state [2]; awaiting a canceled task throws
`OperationCanceledException` (`TaskCanceledException` is the wrapped form) [2]. A token
passed to `Task.Run` cancels the work if it has not yet started [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `await CountKmersAsync(S, k)` equals `CountKmers(S, k)` (same keys and counts) | The async method runs the identical synchronous algorithm via `Task.Run` [3]; counts depend only on the k-mer definition [1] |
| INV-02 | Σ counts = *L − k + 1* for 1 ≤ *k* ≤ *L* | Sliding window yields one k-mer per start position [1] |
| INV-03 | A signaled cancellation token ⇒ awaiting the task throws `OperationCanceledException` | Cooperative cancellation model: `ThrowIfCancellationRequested` / pre-start cancellation [2][3] |
| INV-04 | Empty/null *S* or *k* > *L* ⇒ empty result; *k* ≤ 0 (non-empty *S*) ⇒ `ArgumentOutOfRangeException` thrown synchronously from the call | k-mer multiset empty when *L − k + 1* ≤ 0 [1]; TAP: usage errors are thrown directly from the call, all other outcomes via the task [4] |
| INV-05 | Progress: one report i/(L−k+1) per 1000-window checkpoint (strictly increasing, in [0,1)), then exactly one final 1.0 on every successful completion (also for empty results); none after cancellation | Progress reported synchronously to `IProgress<T>` [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Sequence to analyze | Null/empty ⇒ empty result; normalized to uppercase |
| k | int | required | K-mer length | Must be > 0; *k* > *L* ⇒ empty result |
| cancellationToken | CancellationToken | `default` | Cooperative cancellation token | Signaled ⇒ task canceled |
| progress | IProgress&lt;double&gt;? | null | Optional 0.0–1.0 progress reporter (reported synchronously [4]) | Checkpoint fractions i/(L−k+1) every 1000 windows, then final report = 1.0 on every successful completion (including empty results) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (result) | Task&lt;Dictionary&lt;string,int&gt;&gt; | Task whose result maps each k-mer to its occurrence count; identical to the synchronous method |

### 3.3 Preconditions and Validation

Null/empty sequence ⇒ empty dictionary (for any *k*). *k* > *L* ⇒ empty dictionary. *k* ≤ 0 on
non-empty input ⇒ `ArgumentOutOfRangeException`, thrown **synchronously from the
`CountKmersAsync` call** (TAP usage error [4]; also takes precedence over an already-canceled
token). Input is uppercased
(`ToUpperInvariant`), so counting is case-insensitive; the alphabet is not restricted
(non-ACGT characters, e.g. IUPAC `N`, are counted as-is). Indexing is 0-based; windows are
inclusive of length *k*. Cancellation: a signaled token ⇒ awaiting throws
`OperationCanceledException` [2][3].

## 4. Algorithm

### 4.1 High-Level Steps

0. `CountKmersAsync` validates the usage error *k* ≤ 0 (non-empty input) synchronously [4].
1. `CountKmersAsync` queues the synchronous `CountKmers(sequence, k, token, progress)` on
   the thread pool via `Task.Run(..., token)` [3].
2. The synchronous method validates inputs, uppercases the sequence, and slides a
   length-*k* window, incrementing a dictionary count per window [1].
3. Periodically (every `CheckInterval` windows) it calls
   `cancellationToken.ThrowIfCancellationRequested()` and reports progress [2].
4. On completion it reports progress 1.0 and the task's result is the count dictionary.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CountKmersAsync | O(n·k) | O(u·k) | n = L − k + 1 windows, each building a length-*k* string key; u = number of distinct k-mers. Thread-pool offload does not change asymptotic cost |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.CountKmersAsync(string, int, CancellationToken, IProgress<double>)`:
  thread-pool wrapper returning a `Task<Dictionary<string,int>>`.
- `KmerAnalyzer.CountKmers(string, int, CancellationToken, IProgress<double>)`:
  synchronous reference invoked by the async method.
- `KmerAnalyzer.CountKmersSpan(ReadOnlySpan<char>, int)`: span variant (deeply tested
  under KMER-COUNT-001).

### 5.2 Current Behavior

`CountKmersAsync` delegates verbatim to the synchronous cancellation/progress overload,
so its numeric output is byte-for-byte identical to `CountKmers`. The token is passed to
`Task.Run`, so a token already signaled when the call is made cancels the work before the
delegate runs [3]; a token signaled while running is observed at the next
`ThrowIfCancellationRequested` checkpoint [2]. Cancellation checks occur every
`CheckInterval` (1000) windows to bound polling overhead.

**Search-infrastructure decision:** the repository suffix tree was evaluated and not
used. K-mer *counting* is a single linear pass tallying every fixed-length window, not an
occurrence-enumeration query of a fixed pattern against a text; an O(n) sliding window is
the appropriate algorithm and matches the synchronous reference. A suffix tree would add
O(n) construction without improving a one-pass full-spectrum count.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exactly *L − k + 1* k-mers per sequence; per-k-mer occurrence counts [1].
- Empty result for *L* = 0 or *k* > *L* [1].
- Cooperative cancellation via `ThrowIfCancellationRequested`; pre-start cancellation via
  the `Task.Run` token; `OperationCanceledException` on await [2][3].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Parallel partitioning of the window scan across cores; **users should rely on:** the
  current single-threaded thread-pool offload, which already satisfies the async contract
  and matches the synchronous result exactly.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / null sequence | Empty dictionary | *L* = 0 ⇒ no k-mers [1] |
| k > L | Empty dictionary | *L − k + 1* ≤ 0 [1] |
| k ≤ 0 (non-empty S) | `ArgumentOutOfRangeException` thrown from the call (before any task) | TAP usage error [4] |
| Token signaled before call | Task `Canceled`, delegate never runs; awaiting throws `OperationCanceledException` whose `CancellationToken` is the caller's token | Task.Run cancels work not yet started [3]; TAP [4] |
| Token signaled during run | Observed at the next 1000-window checkpoint; task `Canceled` with the same token; no final 1.0 progress | ThrowIfCancellationRequested [2]; TAP [4] |
| Token signaled after the last checkpoint | Result is produced (`RanToCompletion`) | TAP: the operation "need not accept the cancellation request" [4] |
| Lowercase / mixed case | Same counts as uppercase | Input uppercased |

### 6.2 Limitations

Counting is alphabet-agnostic (non-ACGT characters are counted literally). The async
variant offloads to a single thread-pool thread; it does not parallelize the scan, so it
is not faster asymptotically than the synchronous method — its purpose is
non-blocking execution with cancellation and progress. TAP guidance [4] advises exposing
purely compute-bound work only synchronously (letting callers choose `Task.Run`); the async
wrapper is retained for public-API compatibility. The method is thread-safe (static, no shared
mutable state); the `IProgress<double>` callback runs on the thread-pool worker.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var counts = await KmerAnalyzer.CountKmersAsync("ATGG", 3);
// counts = { "ATG": 1, "TGG": 1 }  (Wikipedia: "ATGG has two 3-mers: ATG and TGG")
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_CountKmersAsync_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountKmersAsync_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [KMER-ASYNC-001-Evidence.md](../../../docs/Evidence/KMER-ASYNC-001-Evidence.md)
- Related algorithms: [K-mer_Counting](./K-mer_Counting.md)

## 8. References

1. Wikipedia. 2026. K-mer. https://en.wikipedia.org/wiki/K-mer (accessed 2026-06-14).
2. Microsoft. 2025. Task Cancellation — .NET. Microsoft Learn. https://learn.microsoft.com/en-us/dotnet/standard/parallel-programming/task-cancellation (accessed 2026-06-14).
3. Microsoft. 2025. Task.Run Method (System.Threading.Tasks). Microsoft Learn. https://learn.microsoft.com/en-us/dotnet/api/system.threading.tasks.task.run (accessed 2026-06-14).
4. Microsoft. 2026. Task-based asynchronous pattern (TAP) in .NET; Implementing the Task-based Asynchronous Pattern. Microsoft Learn. https://learn.microsoft.com/en-us/dotnet/standard/asynchronous-programming-patterns/task-based-asynchronous-pattern-tap (source opened as dotnet/docs `docs/standard/asynchronous-programming-patterns/*.md` on raw.githubusercontent.com, 2026-09-28).
