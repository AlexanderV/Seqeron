# Validation Report: KMER-COUNT-001 — K-mer Counting

- **Validated:** 2026-06-24   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.CountKmers(string, int)` (+ overloads: CancellationToken/IProgress, `DnaSequence`, async `CountKmersAsync`); `SequenceExtensions.CountKmersSpan(ReadOnlySpan<char>, int)` (and `KmerAnalyzer.CountKmersSpan` forwarder); `KmerAnalyzer.CountKmersBothStrands(string, int)` (+ `DnaSequence` overload)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — K-mer**: formula "a sequence of length L will have **L − k + 1** k-mers"; the ATGG example "two 3-mers: ATG and TGG"; the Introduction table for `GTAGAGCTGT` (k=2: GT,TA,AG,GA,AG,GC,CT,TG,GT; k=4: GTAG,TAGA,AGAG,GAGC,AGCT,GCTG,CTGT). K-mers are overlapping substrings via a step-1 sliding window; 4^k possible k-mers; "spectrum" = multiplicity tally.
- **Rosalind — KMER (k-Mer Composition)**: A[m] = number of times the m-th lexicographically-ordered k-mer appears in s; all 4^k = 256 4-mers reported. Counting is **non-canonical** (forward strand only; no reverse-complement collapse).
- **Rosalind — BA1E (Clump Finding)**: window-based occurrence counting; consistent with `FindClumps` (out of scope for the count itself).

### Formula / definition check
- Windows = L − k + 1, step 1, overlapping. ✔
- Tally = occurrence count of each distinct k-mer. ✔
- Canonical (min(kmer, revcomp)) collapse is **not** claimed for `CountKmers` — single-strand, matching Rosalind. The separate `CountKmersBothStrands` is an **additive** forward + reverse-complement merge: `count[w] = forward[w] + forward[RC(w)]`, documented (kPAL "balance"; generalized second Chargaff rule; cf. Marçais & Kingsford canonical counting). Distinct, internally consistent purpose. Total = 2·(L − k + 1).
- k > L → empty; k ≤ 0 → `ArgumentOutOfRangeException`. ✔
- Non-ACGT windows (e.g. N) counted as-is, not skipped. ✔ Case normalized to uppercase. ✔

### Edge-case semantics
empty/null → empty dict; k>L → empty dict; k=L → exactly 1 k-mer; k=1 → per-nucleotide composition; lower/mixed case folded. All defined and sourced.

### Independent cross-check (exact numbers)
- Hand: `ACGTACGT`, k=2 → AC,CG,GT,TA,AC,CG,GT → {AC:2, CG:2, GT:2, TA:1}, total 7 = 8−2+1. ✔
- Rosalind 415-bp 4-mer sample (recomputed in the prior validation via `collections.Counter`): length 415, total 412 = 415−4+1, **209 unique**, 47 absent; AAAA=4, AACT=5, AGCT=3, GCTG=5, TATA=2, TTTT=1, CCCC=0 — all match the test assertions.

### Findings / divergences
None. Description is biologically and mathematically correct.

## Stage B — Implementation

### Code path reviewed
- `KmerAnalyzer.cs:20-42` canonical `CountKmers(string,int)`; `:52-88` CancellationToken/IProgress overload; `:93-100` `CountKmersAsync`; `:105-120` `DnaSequence` wrappers; `:125-128` `CountKmersSpan` forwarder; `:476-490` `CountKmersBothStrands(string,int)`; `:501-502` `DnaSequence` overload.
- `SequenceExtensions.cs:353-373` `CountKmersSpan(ReadOnlySpan<char>,int)`.

### Formula realised correctly?
- Loop `for (i = 0; i <= seq.Length - k; i++)` → exactly L−k+1 length-k overlapping windows, step 1. ✔
- `ToUpperInvariant()` applied in all three counting paths (string canonical, CancellationToken overload, span extension). ✔
- Guards: `IsNullOrEmpty → empty`; `k <= 0 → throw`; `k > Length → empty`. ✔
- Non-ACGT (N) passes through unmodified, counted as-is. ✔
- `CountKmersBothStrands` = `CountKmers(seq) ⊕ CountKmers(RC(seq))` additively merged → `forward[w] + forward[RC(w)]`, total 2(L−k+1). ✔

### Cross-verification table recomputed vs code (test run)
| Case | Source value | Code result |
|------|-----|-----|
| ATGG, k=3 | {ATG, TGG} | matches |
| GTAGAGCTGT, k=2 | sum 9, 7 unique, GT=2, AG=2 | matches |
| GTAGAGCTGT, k=4 | 7 unique ×1 | matches |
| GTAGAGCTGT, k=1..10 | L−k+1 series | matches |
| Rosalind 415bp, k=4 | total 412, 209 unique, spot counts | matches |
| ACGTACGT, k=2 (hand) | {AC:2,CG:2,GT:2,TA:1} sum 7 | matches |

### Variant/delegate consistency
- `CountKmersSpan` and the CancellationToken overload normalize case and equal the canonical (locked by `CountKmersSpan_MixedCase_MatchesCountKmers`, `CountKmers_CancellationOverload_NormalizesCase`). ✔
- `CountKmers(DnaSequence,…)` delegates to the string version; `CountKmersBothStrands(DnaSequence,…)` forwards to the string overload. ✔

### Numerical robustness
Counts are `int`; ranges far below overflow. No division except the progress fraction `i/total`, where `total = L − k + 1 ≥ 1` whenever the loop runs (guarded by the k>L early-return). ✔

### Test quality audit
Tests assert exact externally-sourced values (Wikipedia ATGG/GTAGAGCTGT tables, Rosalind 415-bp sample with 209 unique + 7 spot counts, hand-computed overlaps), not implementation echoes. Edge cases (empty, null, k>L, k=L, k=1, k≤0 throw, lowercase, mixed case, N base) covered and deterministic. Framework is NUnit.

### Findings / defects
- **Minor (benign):** `CountKmersSpan` orders its guards `k <= 0` (throw) before `Length < k` (return empty), whereas `CountKmers(string,…)` checks `IsNullOrEmpty` (return empty) before `k <= 0`. So `CountKmersSpan(<empty>, 0)` throws while `CountKmers("", 0)` returns empty. Both are degenerate invalid-input cases with no biological meaning and no Stage-A-mandated behaviour; not a defect against the validated description.

## Verdict & follow-ups
- **Stage A: PASS** — L−k+1 overlapping step-1 windows, non-canonical (single-strand) tally, confirmed against Wikipedia + Rosalind; worked example and full Rosalind sample independently recomputed.
- **Stage B: PASS** — code faithfully realises the validated description across all variants; edge cases handled as specified; tests lock sourced values.
- **End-state: CLEAN** — no defect requiring a fix. CountKmers test filter: 120 passed / 0 failed. No code changes (implementation unchanged since the prior CLEAN validation; only a doc note on the harmless guard-ordering asymmetry).
