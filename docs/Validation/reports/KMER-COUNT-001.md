# Validation Report: KMER-COUNT-001 — K-mer Counting

- **Validated:** 2026-06-12   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.CountKmers(string, int)` (+ overloads: CancellationToken, `DnaSequence`, async); `SequenceExtensions.CountKmersSpan(ReadOnlySpan<char>, int)`; `KmerAnalyzer.CountKmersBothStrands(DnaSequence, int)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — K-mer** (https://en.wikipedia.org/wiki/K-mer): confirms
  - Formula: "A sequence of length L will have **L − k + 1** k-mers."
  - ATGG example: "The sequence ATGG has two 3-mers: ATG and TGG."
  - Introduction table for `GTAGAGCTGT`: k=2 → `GT, TA, AG, GA, AG, GC, CT, TG, GT`; k=4 → `GTAG, TAGA, AGAG, GAGC, AGCT, GCTG, CTGT`.
  - K-mers are overlapping substrings extracted with a step-1 sliding window; 4^k (n^k) total possible k-mers; "k-mer spectrum" = multiplicity tallying.
- **Rosalind — KMER (k-Mer Composition)** (https://rosalind.info/problems/kmer/): k-mer composition is the array A where A[m] = number of times the m-th lexicographically-ordered k-mer appears in s; all 4^k=256 4-mers reported in lexicographic order (AAAA…TTTT). **Counting is non-canonical** — forward strand only, no reverse-complement collapse.
- **Rosalind — BA1E (Clump Finding)**: confirms window-based k-mer occurrence counting used by `FindClumps` (out of scope for the count itself but consistent).

### Formula / definition check
- Number of windows = L − k + 1, step 1, overlapping. ✔
- Tally = count of occurrences of each distinct k-mer. ✔
- **Canonical/reverse-complement collapse: NOT claimed** by spec for `CountKmers`. The spec and Rosalind both count plain (single-strand) k-mers. The separate `CountKmersBothStrands` variant is an *additive* forward + reverse-complement combination (documented as such), not a canonical min(kmer, revcomp) collapse — this is its stated, distinct purpose, and is internally consistent.
- k > L → 0 k-mers (empty); k ≤ 0 → invalid (spec: throw ArgumentOutOfRangeException). ✔
- Non-ACGT windows (e.g. N) are **counted as-is**, not skipped (matches "k-mer = any length-k substring"). ✔

### Edge-case semantics
Empty/null → empty dict; k>L → empty dict; k=L → exactly 1 k-mer; k=1 → per-nucleotide counts; case normalized to uppercase. All defined and sourced.

### Independent cross-check (exact numbers)
- Hand-computed worked example: `ACGTACGT`, k=2 → windows AC,CG,GT,TA,AC,CG,GT → **{AC:2, CG:2, GT:2, TA:1}**, total 7 = 8−2+1. ✔
- Independently recomputed the Rosalind 415-bp sample (Python `collections.Counter`): length **415**, total **412** (=415−4+1), **209** unique 4-mers, **47** zeros; AAAA=4, AACT=5, AGCT=3, GCTG=5, TATA=2, TTTT=1, CCCC=0. All match the test assertions exactly. (Note: the WebFetch preview truncated the displayed sample to ~309 bp, but the full Rosalind sample / the value baked into the test is 415 bp — verified by direct recount.)

### Findings
No divergences. Description is biologically and mathematically correct.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs:20-42` (canonical `CountKmers`), `:52-88` (CancellationToken overload), `:105-120` (DnaSequence wrappers), `:371-387` (CountKmersBothStrands).
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:225-245` (`CountKmersSpan`).

### Formula realised correctly?
- Loop `for (int i = 0; i <= seq.Length - k; i++)` yields exactly L−k+1 windows of length k, step 1, overlapping. ✔
- `seq.ToUpperInvariant()` normalizes case before tallying. ✔ (all three counting paths: string, CancellationToken span, and `CountKmersSpan` extension each apply `ToUpperInvariant`.)
- Guards: `string.IsNullOrEmpty → empty`; `k <= 0 → ArgumentOutOfRangeException`; `k > Length → empty`. ✔ Matches Stage-A semantics exactly. (`CountKmersSpan` checks `k<=0` then `sequence.Length < k`, equivalent.)
- Non-ACGT (N) characters pass through unmodified, counted as-is. ✔

### Cross-verification table recomputed vs code (via test run)
| Case | Source value | Code result |
|------|-----|-----|
| ATGG, k=3 | {ATG, TGG} | matches |
| GTAGAGCTGT, k=2 | sum 9, 7 unique, GT=2, AG=2 | matches |
| GTAGAGCTGT, k=4 | 7 unique each ×1 | matches |
| GTAGAGCTGT, k=1..10 | L−k+1 series | matches |
| Rosalind 415bp, k=4 | total 412, 209 unique, spot counts | matches |
| ACGTACGT, k=2 (hand) | {AC:2,CG:2,GT:2,TA:1} sum 7 | matches |

### Variant/delegate consistency
- `CountKmersSpan` and CancellationToken overload both normalize case and produce identical results to the canonical (locked by regression tests `CountKmersSpan_MixedCase_MatchesCountKmers`, `CountKmers_CancellationOverload_NormalizesCase`). ✔
- `CountKmers(DnaSequence,…)` delegates to the string version. ✔
- `CountKmersBothStrands` = forward + reverse-complement additive merge; invariant sum = 2(L−k+1) verified. ✔

### Test quality audit
36-spec / 56-runtime tests assert exact externally-sourced values (Wikipedia ATGG/GTAGAGCTGT tables, Rosalind 415-bp sample with 209 unique + spot counts, hand-computed overlaps), not implementation echoes. Edge cases (empty, null, k>L, k=L, k=1, k≤0 throw, lowercase, mixed case, N base) all covered. Tests are deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — L−k+1, overlapping step-1 windows, non-canonical (single-strand) tallying, all confirmed against Wikipedia + Rosalind; worked example and full Rosalind sample independently recomputed.
- **Stage B: PASS** — code faithfully realises the validated description across all variants; edge cases handled as specified; tests lock sourced values.
- **End-state: CLEAN** — no defect. CountKmers filter: 56 passed / 0 failed. Full `Seqeron.Genomics.Tests` suite: 4461 passed / 0 failed. No code changes.
