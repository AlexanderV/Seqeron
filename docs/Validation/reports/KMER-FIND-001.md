# Validation Report: KMER-FIND-001 — Frequent Words / Unique k-mers / Clump Finding

- **Validated:** 2026-06-12   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.FindMostFrequentKmers(seq, k)`, `KmerAnalyzer.FindUniqueKmers(seq, k)`, `KmerAnalyzer.FindClumps(seq, k, windowSize, minOccurrences)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_Find_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Rosalind BA1B — Frequent Words Problem** (https://rosalind.info/problems/ba1b/):
  - "Find the most frequent k-mers in a string." A most frequent k-mer is a Pattern
    that **maximizes Count(Text, Pattern)** among all k-mers.
  - `Count(Text, Pattern)` counts **overlapping** occurrences (standard substring count).
  - Return: **"All most frequent k-mers in Text (in any order)"** → all ties returned.
  - Sample: `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4 → output `CATG GCAT`.
- **Rosalind BA1E — Clump Finding Problem** (https://rosalind.info/problems/ba1e/):
  - A Pattern forms an **(L,t)-clump** inside Genome if there is an interval of Genome
    of **length L** in which Pattern appears **at least t times** (with overlaps).
  - Parameters: `k` = pattern length, `L` = window/interval length, `t` = min occurrences.
  - Return: **all distinct k-mers** forming (L,t)-clumps (a set; order-independent).
  - Sample: Genome `CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC`,
    k=5, L=75, t=4 → output `CGACA GAAGA AATGT`.
- **Wikipedia (k-mer):** a k-mer is a length-k substring; a sequence of length L yields
  `L − k + 1` k-mers (overlapping); if L < k there are no k-mers. Unique k-mers (count = 1)
  are used for genomic fingerprinting.

### Formula / definition check
- Frequent words = argmax over overlapping counts, returning the full tie set. ✓
- (L,t)-clump = ∃ window of length L containing ≥ t overlapping occurrences of the k-mer;
  output is the distinct set. ✓ Parameter meanings (k, L=window, t=threshold) match. ✓
- Position convention: neither problem returns positions; only the k-mer strings. ✓

### Edge-case semantics
- All k-mers distinct → every k-mer is "most frequent" (max count = 1) and every k-mer is unique. ✓
- t larger than any window count → empty clump set. ✓
- Window length L > genome length → no valid window → empty. ✓
- k > L → no k-mer fits in a window → empty. ✓
- L < k (and empty input) → no k-mers. ✓

### Independent cross-check (hand computation)
- **BA1B**, `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4:
  - `CATG` at 0-based positions 6, 13, 19 → 3 occurrences.
  - `GCAT` at 0-based positions 5, 12, 18 → 3 occurrences.
  - No other 4-mer reaches 3 → max = 3 → **{CATG, GCAT}**. Matches Rosalind. ✓
- **BA1E** sample → **{CGACA, GAAGA, AATGT}** per Rosalind. ✓ (verified against code output)

### Findings / divergences
None. The spec's stated semantics match Rosalind BA1B/BA1E and Wikipedia exactly.

## Stage B — Implementation

### Code path reviewed
- `FindMostFrequentKmers` — `KmerAnalyzer.cs:156-169`
- `CountKmers` (overlapping count, `i <= len-k`) — `KmerAnalyzer.cs:20-42`
- `FindUniqueKmers` — `KmerAnalyzer.cs:216-220`
- `FindClumps` — `KmerAnalyzer.cs:297-348`

### Formula realised correctly? (evidence)
- **FindMostFrequentKmers:** counts overlapping k-mers, computes `maxCount = counts.Values.Max()`,
  returns every k-mer with `Value == maxCount` → all ties returned. Empty seq / k > length →
  `CountKmers` returns an empty dict → `counts.Count == 0` → yields nothing. ✓
- **FindUniqueKmers:** returns k-mers with `Value == 1`. ✓
- **FindClumps:** sliding window of length L. Initial window seeds counts for the
  `windowSize - k + 1` k-mers in `[0, windowSize)`. Slide loop `i = 1 .. len - windowSize`:
  removes the k-mer leaving (position `i-1`), adds the k-mer entering (position
  `i + windowSize - k`, the last k-mer fully inside the new window `[i, i+windowSize)`),
  then collects any k-mer with `count >= minOccurrences` into a `HashSet` (distinct set). ✓
  Overlapping occurrences are counted (window populated from contiguous overlapping k-mers). ✓
  Guards: empty / `k<=0` / `windowSize<k` (k>L) / `minOccurrences<=0` / `windowSize>length`
  all yield empty. ✓

### Cross-verification table recomputed vs code
| Case | Input | Expected (source) | Code output | Match |
|------|-------|-------------------|-------------|-------|
| BA1B | `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4 | {CATG, GCAT} | {CATG, GCAT} | ✓ |
| BA1E | sample genome, k=5, L=75, t=4 | {CGACA, GAAGA, AATGT} | {CGACA, GAAGA, AATGT} | ✓ |
| Unique | `ACGTACGT`, k=4 | {CGTA, GTAC, TACG} | same | ✓ |
| Homopolymer unique | `AAAA`, k=2 | {} | {} | ✓ |
| Clump simple | `AAAAA`, k=3, L=5, t=3 | {AAA} | contains AAA | ✓ |
| k > L | `ACGTACGT`, k=5, L=4 | {} | {} | ✓ |
| L > length | `ACGT`, L=10 | {} | {} | ✓ |

(BA1B and BA1E outputs verified to be *exactly* the expected sets — see test changes below.)

### Variant/delegate consistency
- `FindMostFrequentKmers`, `FindUniqueKmers`, and `FindClumps` all build on the same
  overlapping `CountKmers` counting (or an equivalent dictionary), so counting is consistent.

### Test quality audit
- 19 tests covering M1–M12 + S1–S5. Originally M1 (BA1B) and M9 (BA1E) used only
  `Does.Contain` assertions, which would not catch an implementation returning *extra* or
  *more frequent* k-mers. **Strengthened both to `Is.EquivalentTo` the exact sourced sets**
  `{CATG, GCAT}` and `{CGACA, GAAGA, AATGT}` to lock the Rosalind reference values.
- All other tests assert concrete sourced values and edge-case behaviour; deterministic.

### Findings / defects
No correctness defect in the implementation. One test-strength improvement applied (above).

## Verdict & follow-ups
- **Stage A: PASS** — definitions, parameters, overlapping counting, tie/set semantics and
  both sample outputs independently confirmed against Rosalind BA1B/BA1E and Wikipedia.
- **Stage B: PASS** — code faithfully implements all three operations; cross-check values
  reproduced exactly.
- **State: CLEAN.** Code correct; M1/M9 tests strengthened to exact-set assertions.
- **Tests:** `KmerAnalyzer_Find_Tests` 19/19 pass; full `Seqeron.Genomics.Tests` 4461 pass
  (one intermittent flaky perf test, unrelated to k-mer logic, passed on re-run).
- **Files changed:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_Find_Tests.cs`.
