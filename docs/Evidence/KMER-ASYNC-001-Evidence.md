# Evidence Artifact: KMER-ASYNC-001

**Test Unit ID:** KMER-ASYNC-001
**Algorithm:** Asynchronous K-mer Counting
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — K-mer

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14 (retrieved via WebFetch on the URL above)
**Authority rank:** 4 (Wikipedia article citing primary sources; used for the
combinatorial L − k + 1 definition, which is the canonical k-mer count formula)

**Key Extracted Points:**

1. **K-mer count formula (verbatim):** "a sequence of length L will have L − k + 1
   k-mers and there exist n^k total possible k-mers" (n = alphabet size, 4 for DNA).
   This is the count contract the asynchronous method must reproduce exactly.
2. **Worked example ATGG:** "The sequence ATGG has two 3-mers: ATG and TGG."
   (4 − 3 + 1 = 2.)
3. **Worked sequence GTAGAGCTGT** (L = 10). The *total* k-mer count is the authoritative
   quantity from the formula, total k-mers = L − k + 1:

   | k | Total k-mers (L−k+1) |
   |---|----------------------|
   | 1 | 10 |
   | 2 | 9  |
   | 3 | 8  |
   | 4 | 7  |
   | 5 | 6  |
   | 6 | 5  |
   | 7 | 4  |
   | 8 | 3  |
   | 9 | 2  |
   | 10| 1  |

   Distinct-k-mer counts are derived directly from the sequence (the Wikipedia
   "examples" listing is not a full distinct count and was not used as authority):
   for k=2 the windows are GT,TA,AG,GA,AG,GC,CT,TG,GT ⇒ 7 distinct (GT and AG each occur
   twice); for k=3 all 8 windows are distinct; for k=4 all 7 windows are distinct.

### Microsoft Learn — Task Cancellation (.NET)

**URL:** https://learn.microsoft.com/en-us/dotnet/standard/parallel-programming/task-cancellation
**Accessed:** 2026-06-14 (retrieved via WebFetch on the URL above)
**Authority rank:** 2 (official .NET platform specification for the cancellation contract)

**Key Extracted Points:**

1. **Cooperative cancellation:** "cancellation involves cooperation between the user
   delegate ... and the code that requested the cancellation. A successful cancellation
   involves the requesting code calling the CancellationTokenSource.Cancel method and
   the user delegate terminating the operation in a timely manner."
2. **Throwing to acknowledge cancellation (verbatim):** "By throwing an
   OperationCanceledException and passing it the token on which cancellation was
   requested. The preferred way to perform is to use the ThrowIfCancellationRequested
   method. A task that's canceled in this way transitions to the Canceled state."
3. **Awaiting a canceled task (verbatim):** "If you're waiting on a Task that
   transitions to the Canceled state, a System.Threading.Tasks.TaskCanceledException
   exception (wrapped in an AggregateException exception) is thrown." (When `await`-ed,
   the unwrapped exception is observed; `TaskCanceledException` derives from
   `OperationCanceledException`.)

### Microsoft Learn — Task.Run Method

**URL:** https://learn.microsoft.com/en-us/dotnet/api/system.threading.tasks.task.run
**Accessed:** 2026-06-14 (retrieved via WebFetch on the URL above)
**Authority rank:** 2 (official .NET API reference)

**Key Extracted Points:**

1. **Pre-start cancellation (verbatim):** "A cancellation token allows the work to be
   cancelled if it has not yet started." So passing an already-signaled token to
   `Task.Run(func, token)` produces a Canceled task without running the delegate.
2. **Stored exception (verbatim from the OperationCanceledException exception entry):**
   "The task has been canceled. This exception is stored into the returned task." —
   i.e. awaiting the task surfaces an `OperationCanceledException`.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia — K-mer

1. **k > L:** L − k + 1 ≤ 0 ⇒ no k-mers exist (empty result).
2. **Empty sequence:** L = 0 ⇒ no k-mers.

### From Microsoft Learn — Task Cancellation / Task.Run

1. **Token signaled before start:** delegate may not run; task is Canceled; awaiting
   throws `OperationCanceledException`.
2. **Token signaled during run:** delegate calls `ThrowIfCancellationRequested`,
   transitions task to Canceled; awaiting throws `OperationCanceledException`.

---

## Test Datasets

### Dataset: GTAGAGCTGT (Wikipedia K-mer table)

**Source:** Wikipedia — K-mer, https://en.wikipedia.org/wiki/K-mer (accessed 2026-06-14)

| Parameter | Value |
|-----------|-------|
| Sequence | GTAGAGCTGT |
| Length L | 10 |
| k = 2 → total / distinct | 9 / 7 |
| k = 3 → total / distinct | 8 / 8 |
| k = 4 → total / distinct | 7 / 7 |

### Dataset: ATGG (Wikipedia worked example)

**Source:** Wikipedia — K-mer (accessed 2026-06-14)

| Parameter | Value |
|-----------|-------|
| Sequence | ATGG |
| k | 3 |
| 3-mers | ATG (×1), TGG (×1) |

---

## Assumptions

1. **ASSUMPTION: Numeric correctness contract is identical to the synchronous method.**
   The asynchronous API `CountKmersAsync` is, per the repository design and the .NET
   `Task.Run` model, a parallel/async wrapper over the validated synchronous
   `CountKmers(string, int, CancellationToken, IProgress<double>)` method
   (KMER-COUNT-001). No authoritative source defines a *different* numeric result for an
   "async k-mer count"; the counts/frequencies are defined purely by the k-mer formula
   (L − k + 1). This assumption is therefore non-correctness-affecting for numeric
   output: the async result MUST equal the synchronous reference on identical input.
   Verified by tests asserting `CountKmersAsync(seq,k) == CountKmers(seq,k)`.

---

## Recommendations for Test Coverage

1. **MUST Test:** Async result equals synchronous reference for the Wikipedia
   GTAGAGCTGT dataset (k=2,3,4) and the ATGG example — Evidence: Wikipedia K-mer table.
2. **MUST Test:** Total-count invariant Σcounts = L − k + 1 holds for the async result —
   Evidence: Wikipedia "a sequence of length L will have L − k + 1 k-mers".
3. **MUST Test:** Empty / null sequence and k > L yield an empty dictionary (async) —
   Evidence: Wikipedia corner cases.
4. **MUST Test:** k ≤ 0 surfaces `ArgumentOutOfRangeException` via the awaited task —
   Evidence: synchronous contract (KMER-COUNT-001) preserved through `Task.Run`.
5. **MUST Test:** A token already signaled (and one signaled mid-run on a large input)
   causes awaiting the task to throw `OperationCanceledException` — Evidence: Microsoft
   Learn Task Cancellation / Task.Run.
6. **SHOULD Test:** Large input (≥ check interval) returns the same result as the
   synchronous method and reports progress to 1.0 — Rationale: async path is intended
   for long-running inputs; determinism on large input must be confirmed.
7. **COULD Test:** Default (`CancellationToken.None`) and case-insensitivity carried
   through the async path — Rationale: API consistency with the canonical method.

---

## References

1. Wikipedia. 2026. K-mer. https://en.wikipedia.org/wiki/K-mer (accessed 2026-06-14).
2. Microsoft. 2025. Task Cancellation — .NET. Microsoft Learn.
   https://learn.microsoft.com/en-us/dotnet/standard/parallel-programming/task-cancellation
   (accessed 2026-06-14).
3. Microsoft. 2025. Task.Run Method (System.Threading.Tasks). Microsoft Learn.
   https://learn.microsoft.com/en-us/dotnet/api/system.threading.tasks.task.run
   (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation for KMER-ASYNC-001.
