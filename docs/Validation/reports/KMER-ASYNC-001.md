# Validation Report: KMER-ASYNC-001 — Asynchronous K-mer Counting

- **Validated:** 2026-06-16   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.CountKmersAsync(string, int, CancellationToken, IProgress<double>)`
  (thread-pool wrapper over `CountKmers(string, int, CancellationToken, IProgress<double>)`);
  `CountKmersSpan` delegate (smoke only).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened (retrieved this session)

- **Wikipedia — K-mer** (`https://en.wikipedia.org/wiki/K-mer`, WebFetch 2026-06-16). Confirmed verbatim:
  - Formula: *"a sequence of length L will have L − k + 1 k-mers and there exist n^k total possible
    k-mers, where n is number of possible monomers."*
  - ATGG example: *"The sequence ATGG has two 3-mers: ATG and TGG."*
  - GTAGAGCTGT table (total k-mers): k=2→9, k=3→8, k=4→7, k=5→6, k=6→5, k=7→4, k=8→3, k=9→2,
    k=10→1; **k=1 → "4 distinct k-mers"** (the only distinct figure the source gives for this sequence).
- **Microsoft Learn — Task Cancellation (.NET)** (WebFetch 2026-06-16). Confirmed verbatim: throwing
  `OperationCanceledException` via `ThrowIfCancellationRequested` transitions the task to the **Canceled**
  state; *"If you're waiting on a Task that transitions to the Canceled state, a … `TaskCanceledException`
  exception (wrapped in an `AggregateException`) is thrown."* (`TaskCanceledException` derives from
  `OperationCanceledException`; when `await`-ed the unwrapped form is observed).
- **Microsoft Learn — Task.Run** (per Evidence; cancellation-token pre-start semantics consistent with the
  Task Cancellation page already fetched): a token passed to `Task.Run` cancels work that has not yet started.

### Formula check
L − k + 1 overlapping windows; per-k-mer count = number of start positions yielding that k-mer; empty when
L = 0 or k > L (L − k + 1 ≤ 0). Matches the cited source exactly. The "async" qualifier governs **execution**
(thread-pool offload + cooperative cancellation + progress), not the numeric definition — so the async result
MUST equal the synchronous reference on identical input (INV-1). No authoritative source defines a different
numeric result for an "async k-mer count"; assumption A1 is sound.

### Edge-case semantics
Empty/null sequence → empty; k > L → empty; k ≤ 0 → `ArgumentOutOfRangeException`; signaled token →
`OperationCanceledException` on await. All defined and sourced (Wikipedia + Microsoft Learn).

### Independent cross-check (hand-computed on GTAGAGCTGT, L=10)
- k=2 windows: GT,TA,AG,GA,AG,GC,CT,TG,GT → 9 total, **7 distinct** (GT@{0,8}=2, AG@{2,4}=2).
- k=3 windows: GTA,TAG,AGA,GAG,AGC,GCT,CTG,TGT → 8 total, **8 distinct** (none repeat).
- k=4 windows: GTAG,TAGA,AGAG,GAGC,AGCT,GCTG,CTGT → 7 total, **7 distinct**.
- k=1: bases G,T,A,C → **4 distinct** (matches Wikipedia's explicit "4 distinct k-mers"), 10 total.
- ATGG, k=3: ATG×1, TGG×1 (matches Wikipedia).
All independently confirmed; the distinct counts for k≥2 are sequence-derived (Wikipedia tabulates only totals
for k≥2 and distinct only for k=1).

### Findings / divergences
None. Description is biologically and computationally correct and matches the sources. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`:
- `CountKmersAsync` (L93–100): `Task.Run(() => CountKmers(sequence, k, cancellationToken, progress), cancellationToken)`.
  Token passed to both the delegate and `Task.Run` — pre-start cancellation honored; in-run cancellation observed
  at the `ThrowIfCancellationRequested` checkpoints.
- `CountKmers(string,int,CancellationToken,IProgress<double>)` (L52–88): validation (null/empty → empty;
  k≤0 → `ArgumentOutOfRangeException`; k>L → empty), `ToUpperInvariant`, sliding window with `TryAdd`/increment,
  cancellation + progress every 1000 windows, final `progress.Report(1.0)`.

### Formula realised correctly?
Yes. Single linear O(n·k) pass over L−k+1 windows; the async path delegates **verbatim** to the synchronous
method, so its numeric output is byte-for-byte identical (INV-1). The work runs on a single thread-pool thread
(no parallel partitioning), so there are no concurrent writers to the dictionary — **no lost-update / race
risk**, and the result is deterministic. Cancellation/await contract matches the .NET model exactly.

### Cross-verification table recomputed vs code (tests executed this session)
| Input | k | Total (L−k+1) | Distinct | Source |
|-------|---|---------------|----------|--------|
| GTAGAGCTGT | 1 | 10 | 4 | Wikipedia (explicit "4 distinct") |
| GTAGAGCTGT | 2 | 9 | 7 (GT=2, AG=2) | Wikipedia total + hand-derived distinct |
| GTAGAGCTGT | 3 | 8 | 8 | Wikipedia total + hand-derived distinct |
| GTAGAGCTGT | 4 | 7 | 7 | Wikipedia total + hand-derived distinct |
| ATGG | 3 | 2 | 2 ({ATG:1,TGG:1}) | Wikipedia |

All values reproduced by the running async code (fixture green, 16/16).

### Variant/delegate consistency
`CountKmersSpan` and the sync `CountKmers` agree with the async result (C2, S1, S3). Default-token path (S3) and
case-insensitivity (S2) carried through the async wrapper correctly.

### Test quality audit (test-quality gate)
- **Sourced expectations, not code echoes:** the async==sync equality assertions (the INV-1 / async-correctness
  check) are each additionally anchored to **sourced absolute values** — `Sum()==L−k+1` (Wikipedia formula) and
  exact distinct/per-k-mer counts (M1 8/8, M2 9/7 + GT=2/AG=2, M3 {ATG:1,TGG:1}, M4 k=1..10 totals, new M4b k=1
  distinct=4). A deliberately-wrong implementation that broke the count would fail the absolute assertions, not
  merely the equality. PASS.
- **No green-washing:** all assertions use exact `Is.EqualTo`; cancellation tests assert exact thrown type
  (`OperationCanceledException`, matching `TaskCanceledException` subclass); no weakened/widened/skipped tests.
- **Coverage:** all public surface exercised — async overload (4-arg + default-token), sync delegate, span
  delegate, progress; all Stage-A branches/edges (empty, null, k>L, k≤0, pre-signaled token, mid-run cancel on
  >check-interval input, lowercase, large determinism).
- **Honest green:** full unfiltered suite **6607 passed, Failed: 0** (async fixture 15→16 after adding M4b);
  build 0 errors, no new warnings.

### Findings / defects
1. **Stale doc-comment (minor, Stage B):** M1's `<summary>` claimed *"k=3 yields 8 total, 6 unique k-mers"* —
   factually wrong (k=3 has 8 *distinct* 3-mers; the assertion itself already used 8, and the inline comment on
   the next line listed all 8). **Fixed** this session: comment corrected to the verified 8-distinct value with
   the enumerated windows. No code/assertion change needed.
2. **Coverage strengthening:** added **M4b** (`CountKmersAsync_Gtagagctgt_K1_HasFourDistinctBases`) locking
   Wikipedia's only explicit distinct figure (k=1 → 4 distinct, 10 total) through the async path — the one
   externally-tabulated distinct value previously unasserted.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches Wikipedia + Microsoft Learn verbatim; formula, edge cases, and
  cancellation contract are correct and sourced.
- **Stage B: PASS-WITH-NOTES.** Implementation is correct, deterministic, race-free, and identical to the sync
  reference. One stale/incorrect test doc-comment fixed; one sourced distinct-count test added. No algorithmic
  defect.
- **End-state: ✅ CLEAN.** All findings fully resolved this session; build 0 errors, full unfiltered suite 6607
  passed / 0 failed.
