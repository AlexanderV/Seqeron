# Test Specification: KMER-ASYNC-001

**Test Unit ID:** KMER-ASYNC-001
**Area:** K-mer Analysis
**Algorithm:** Asynchronous K-mer Counting
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — K-mer (L − k + 1 formula, GTAGAGCTGT table, ATGG example) | 4 | https://en.wikipedia.org/wiki/K-mer | 2026-06-14 |
| 2 | Microsoft Learn — Task Cancellation (.NET) | 2 | https://learn.microsoft.com/en-us/dotnet/standard/parallel-programming/task-cancellation | 2026-06-14 |
| 3 | Microsoft Learn — Task.Run Method | 2 | https://learn.microsoft.com/en-us/dotnet/api/system.threading.tasks.task.run | 2026-06-14 |

### 1.2 Key Evidence Points

1. A sequence of length L has exactly L − k + 1 k-mers — Wikipedia K-mer.
2. The async method is a `Task.Run` wrapper over the synchronous reference; numeric
   counts are defined solely by the k-mer formula, so the async result MUST equal the
   synchronous result on identical input — Wikipedia K-mer + repository design.
3. Throwing `OperationCanceledException` via `ThrowIfCancellationRequested` transitions
   the task to Canceled; awaiting a canceled task throws `OperationCanceledException`
   (a `TaskCanceledException` is the wrapped form) — Microsoft Learn Task Cancellation.
4. `Task.Run(func, token)` cancels the work if it has not yet started; the canceled
   exception is stored in the returned task — Microsoft Learn Task.Run.

### 1.3 Documented Corner Cases

- Empty sequence (L = 0) and k > L ⇒ no k-mers (empty dictionary) — Wikipedia.
- Token signaled before/while running ⇒ task Canceled; awaiting throws — Microsoft Learn.

### 1.4 Known Failure Modes / Pitfalls

1. Race-dependent assertions on async cancellation — mitigated by awaiting completion
   and asserting the thrown exception type, never timing — Microsoft Learn Task Cancellation.
2. k ≤ 0 must surface `ArgumentOutOfRangeException` from the awaited task (preserved
   through `Task.Run`) — synchronous contract (KMER-COUNT-001).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CountKmersAsync(string, int, CancellationToken, IProgress<double>)` | KmerAnalyzer | **Canonical** | Async wrapper; result MUST equal synchronous `CountKmers`; cancellation contract |
| `CountKmers(string, int, CancellationToken, IProgress<double>)` | KmerAnalyzer | **Internal** | Synchronous reference the async method delegates to (deep-tested in KMER-COUNT-001) |
| `CountKmersSpan(ReadOnlySpan<char>, int)` | KmerAnalyzer | **Delegate** | Span variant; deeply tested under KMER-COUNT-001 — smoke only here |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `await CountKmersAsync(seq, k)` returns a dictionary equal (keys + counts) to `CountKmers(seq, k)` for the same input | Yes | Wikipedia K-mer + Task.Run wrapper design |
| INV-2 | Σ(counts) = L − k + 1 for 1 ≤ k ≤ L (async result) | Yes | Wikipedia: "a sequence of length L will have L − k + 1 k-mers" |
| INV-3 | A signaled cancellation token ⇒ awaiting the task throws `OperationCanceledException` | Yes | Microsoft Learn Task Cancellation / Task.Run |
| INV-4 | Empty/null sequence or k > L ⇒ empty dictionary | Yes | Wikipedia corner cases |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Async == sync (GTAGAGCTGT, k=3) | Determinism: async result identical to synchronous reference | Dictionaries equal (8 total, 8 distinct) | Wikipedia L−k+1; INV-1 |
| M2 | Async == sync (GTAGAGCTGT, k=2 and k=4) | Determinism across multiple k | Equal dicts: k=2 → 9 total/7 distinct (GT,AG ×2), k=4 → 7 total/7 distinct | Wikipedia L−k+1 + sequence-derived distinct; INV-1 |
| M3 | ATGG, k=3 exact counts | Worked example | {ATG:1, TGG:1}, 2 entries | Wikipedia "ATGG has two 3-mers: ATG and TGG" |
| M4 | Total-count invariant | Σcounts = L − k + 1 for k = 1..10 on GTAGAGCTGT | Sums = 10,9,8,7,6,5,4,3,2,1 | Wikipedia formula; INV-2 |
| M5 | Empty sequence | Async on "" | Empty dictionary | Wikipedia (L=0); INV-4 |
| M6 | Null sequence | Async on null | Empty dictionary | Synchronous contract; INV-4 |
| M7 | k > length | Async on "ACG", k=4 | Empty dictionary | Wikipedia (L−k+1≤0); INV-4 |
| M8 | k ≤ 0 | Async on valid seq, k=0 | Awaited task throws `ArgumentOutOfRangeException` | Synchronous contract (KMER-COUNT-001) |
| M9 | Pre-signaled token | Token already canceled before call | Awaiting throws `OperationCanceledException` | Task.Run "cancelled if not yet started"; INV-3 |
| M10 | Cancellation on large input | Token canceled; large (> check interval) sequence | Awaiting throws `OperationCanceledException` | Task Cancellation ThrowIfCancellationRequested; INV-3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Large input determinism | 5000 nt random-but-fixed sequence | Async result == sync result | Async path targets long inputs |
| S2 | Case-insensitive (async) | lowercase input | Same dict as uppercase | API consistency with canonical |
| S3 | Default token | `CancellationToken.None` (no token arg) | Completes; equals sync | Default-parameter path |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Progress reporting | `IProgress<double>` invoked, final value 1.0 | Last reported value = 1.0 | Async progress contract |
| C2 | CountKmersSpan smoke | Span delegate equals CountKmers | Equal dict | Delegate (deep-tested elsewhere) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests exist for `CountKmersAsync`. KMER-COUNT-001
  (`KmerAnalyzer_CountKmers_Tests.cs`) deeply tests the synchronous `CountKmers` and
  `CountKmersSpan`; it has only a smoke test for the cancellation overload, none for the
  async wrapper. This unit adds the canonical async fixture.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (async==sync k=3) | ❌ Missing | New unit |
| M2 (async==sync k=2,k=4) | ❌ Missing | New unit |
| M3 (ATGG k=3) | ❌ Missing | New unit |
| M4 (total-count invariant k=1..10) | ❌ Missing | New unit |
| M5 (empty) | ❌ Missing | New unit |
| M6 (null) | ❌ Missing | New unit |
| M7 (k>L) | ❌ Missing | New unit |
| M8 (k≤0 throws) | ❌ Missing | New unit |
| M9 (pre-signaled token) | ❌ Missing | New unit |
| M10 (large input cancel) | ❌ Missing | New unit |
| S1 (large determinism) | ❌ Missing | New unit |
| S2 (case-insensitive) | ❌ Missing | New unit |
| S3 (default token) | ❌ Missing | New unit |
| C1 (progress) | ❌ Missing | New unit |
| C2 (span smoke) | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmersAsync_Tests.cs` — all async, cancellation, determinism, and edge cases for this unit.
- **Remove:** nothing. Deep `CountKmers`/`CountKmersSpan` coverage stays in KMER-COUNT-001; `CountKmersSpan` is smoke-only here (C2) to avoid duplication.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| KmerAnalyzer_CountKmersAsync_Tests.cs | Canonical fixture for KMER-ASYNC-001 | 15 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented async==sync, total=8, distinct=8 | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented k=2/k=4 equality + exact distinct/repeat counts | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented ATGG → {ATG:1,TGG:1} | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented Σcounts=L−k+1, k=1..10 | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented empty → empty | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented null → empty | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented k>L → empty | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented k=0/-1 → ArgumentOutOfRangeException | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented pre-signaled token → OperationCanceledException | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented large+canceled → OperationCanceledException | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented 5000 nt async==sync | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented lowercase==uppercase | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented default-token equals sync | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented progress final value 1.0 | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented CountKmersSpan smoke | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | CountKmersAsync_Gtagagctgt_K3_EqualsSynchronousReference |
| M2 | ✅ Covered | CountKmersAsync_Gtagagctgt_K2AndK4_EqualSynchronousReference |
| M3 | ✅ Covered | CountKmersAsync_Atgg_K3_ReturnsAtgAndTgg |
| M4 | ✅ Covered | CountKmersAsync_Gtagagctgt_TotalCountInvariant_ForAllK |
| M5 | ✅ Covered | CountKmersAsync_EmptySequence_ReturnsEmptyDictionary |
| M6 | ✅ Covered | CountKmersAsync_NullSequence_ReturnsEmptyDictionary |
| M7 | ✅ Covered | CountKmersAsync_KLargerThanSequence_ReturnsEmptyDictionary |
| M8 | ✅ Covered | CountKmersAsync_InvalidK_ThrowsArgumentOutOfRangeException |
| M9 | ✅ Covered | CountKmersAsync_TokenSignaledBeforeCall_ThrowsOperationCanceled |
| M10 | ✅ Covered | CountKmersAsync_CancelledLargeInput_ThrowsOperationCanceled |
| S1 | ✅ Covered | CountKmersAsync_LargeInput_EqualsSynchronousReference |
| S2 | ✅ Covered | CountKmersAsync_LowercaseInput_EqualsUppercaseResult |
| S3 | ✅ Covered | CountKmersAsync_DefaultToken_CompletesAndEqualsSync |
| C1 | ✅ Covered | CountKmersAsync_WithProgress_ReportsCompletion |
| C2 | ✅ Covered | CountKmersSpan_DelegatesToSpanCount_EqualsCountKmers |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Async numeric contract == synchronous reference (no distinct "async" numeric behavior in any source) | INV-1, M1–M4, S1 |

---

## 7. Open Questions / Decisions

1. Decision: `CountKmersSpan` is mapped to KMER-ASYNC-001 in the Registry but is already
   deeply tested under KMER-COUNT-001; here it is a Delegate (smoke only) to avoid
   duplicate canonical coverage. The async method (`CountKmersAsync`) is the genuine
   canonical method for this unit.
