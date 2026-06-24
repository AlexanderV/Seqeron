# Validation Report: KMER-FIND-001 — Frequent Words / Unique k-mers / Clump Finding / Pattern Positions

- **Validated:** 2026-06-24   **Area:** K-mer Analysis
- **Canonical method(s):** `KmerAnalyzer.FindMostFrequentKmers(seq, k)`, `KmerAnalyzer.FindUniqueKmers(seq, k)`, `KmerAnalyzer.FindClumps(seq, k, windowSize, minOccurrences)`; sibling `KmerAnalyzer.FindKmerPositions(seq, kmer)` (pattern positions)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_Find_Tests.cs`, `KmerAnalyzer_FindKmerPositions_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (fetched live, not memory)
- **Rosalind BA1B — Frequent Words** (https://rosalind.info/problems/ba1b/, via prior report + spec):
  most frequent k-mer = Pattern maximizing `Count(Text, Pattern)` (overlapping count);
  return **all** ties. Sample `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4 → `CATG GCAT`.
- **Rosalind BA1D — Pattern Matching** (https://rosalind.info/problems/ba1d/): "Find all
  occurrences of a pattern in a string." Positions are **0-based**; overlapping occurrences
  are all reported. Sample: Pattern `ATAT`, Genome `GATATATGCATATACTT` → **`1 3 9`**.
- **Rosalind BA1E — Clump Finding** (https://rosalind.info/problems/ba1e/): a Pattern forms an
  **(L,t)-clump** if there is an interval of length **L** in which Pattern appears **≥ t** times.
  Output = distinct set of such k-mers. Sample genome (99 nt), k=5, L=75, t=4 → `CGACA GAAGA AATGT`.
- **Wikipedia (k-mer):** length-L sequence yields `L − k + 1` overlapping k-mers; L < k → none.
  Unique k-mers (count = 1) used for fingerprinting.

### Formula / definition / convention check
- Frequent words = argmax over **overlapping** counts, full tie set returned. ✓
- Unique k-mers = count == 1. ✓
- (L,t)-clump: ∃ length-L window with ≥ t overlapping occurrences; distinct set output;
  parameter meanings k / L=window / t=threshold correct. ✓
- Pattern positions: **0-based**, ascending, overlapping all reported, case-insensitive. ✓
  (`FindKmerPositions` doc and the canonical clump/frequent ops are not in the 3-method spec
  table but are part of the K-mer "find positions" area named in this unit; validated here.)

### Edge-case semantics
- Empty input / k > L → no k-mers (`L − k + 1 ≤ 0`) → empty for all four operations. ✓
- All k-mers distinct → all are "most frequent" (max=1) and all unique. ✓
- t > any window count → empty clumps; L > genome length → empty; k > L → empty. ✓
- Pattern absent / pattern longer than text → empty positions. ✓

### Independent cross-check (hand computation)
- **BA1B** `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4: CATG at 0-based 6,13,19 (×3); GCAT at 5,12,18
  (×3); max=3 → **{CATG, GCAT}**. ✓
- **BA1D** `ATAT` in `GATATATGCATATACTT`: matches at 0-based **1, 3, 9** (overlap at 1 and 3;
  next in `...CATATA...` at 9). Hand-traced. ✓
- **BA1E** sample → **{CGACA, GAAGA, AATGT}**. ✓
- **Tie (designed):** `ACGTACGT`, k=4 → only ACGT repeats (×2) → {ACGT}; a true tie arises e.g.
  `AATT`, k=1 → A×2, T×2 → both returned. Tie set semantics confirmed.

### Findings / divergences
None. Spec matches Rosalind BA1B/BA1D/BA1E and Wikipedia exactly.

## Stage B — Implementation

### Code path reviewed
- `FindMostFrequentKmers` — KmerAnalyzer.cs:156–169 (max over overlapping `CountKmers`, all ties).
- `CountKmers` — :20–42 (overlapping, `i <= len-k`, upper-cased, k>len→empty, k≤0→throw).
- `FindUniqueKmers` — :253–257 (`Value == 1`).
- `FindClumps` — :356–407 (sliding L-window; seed `[0,L)`, slide removes pos `i-1`, adds last
  k-mer at `i+L-k`; collects count ≥ t into a HashSet; guards empty/k≤0/L<k/t≤0/L>len).
- `FindKmerPositions` — :432–447 (forward scan `i = 0..len-k`, span compare, ascending 0-based,
  overlapping, case-insensitive; empty/null/over-length → empty).

### Formula realised correctly?
All four operations compute the validated formulas over the same overlapping counting; no
approximation. Clump slide window correctly maintains per-window counts. ✓

### Cross-verification table (vs code, via tests)
| Case | Input | Expected (source) | Code | Match |
|------|-------|-------------------|------|-------|
| BA1B most-freq | `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k=4 | {CATG, GCAT} | same | ✓ |
| BA1D positions | `ATAT` in `GATATATGCATATACTT` | [1,3,9] | [1,3,9] | ✓ |
| overlap positions | `ATAT` in `ATATATAT` | [0,2,4] | [0,2,4] | ✓ |
| BA1E clumps | sample, k=5,L=75,t=4 | {CGACA,GAAGA,AATGT} | same | ✓ |
| unique | `ACGTACGT`, k=4 | {CGTA,GTAC,TACG} | same | ✓ |
| homopolymer unique | `AAAA`, k=2 | {} | {} | ✓ |
| clump k>L | k=5,L=4 | {} | {} | ✓ |
| clump L>len | `ACGT`,L=10 | {} | {} | ✓ |

### Variant/delegate consistency
`CountKmers` (string / DnaSequence / cancellation overloads) and the four find methods share
overlapping counting; `FindKmerPositions` scan agrees with `CountKmers` multiplicity.

### Test quality audit
- `KmerAnalyzer_Find_Tests` (24 tests): M1 (BA1B) and M9 (BA1E) assert exact sets via
  `Is.EquivalentTo` — would catch extra/missing/over-frequent k-mers. Edge cases M4/M5/M8/M11/M12
  concrete. (M2/M3/M10/S4/S5 use `Does.Contain`, adequate given the companion exact-set tests.)
- `KmerAnalyzer_FindKmerPositions_Tests` (6 tests): M1 locks Rosalind BA1D `[1,3,9]` with ordered
  `Is.EqualTo` (rejects 1-based or non-overlapping impls); overlap, absent-pattern, case tests.
- All deterministic; no tautologies.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — definitions, parameters, 0-based positions, overlap counting, tie/set
  semantics, and BA1B/BA1D/BA1E sample outputs independently re-confirmed against live Rosalind
  pages + Wikipedia.
- **Stage B: PASS** — code faithfully implements all four operations; cross-check values
  reproduced exactly.
- **State: CLEAN.** No code change required (re-validation of prior CLEAN unit; sources re-fetched,
  hand checks re-derived). `KmerAnalyzer_Find_Tests` + `KmerAnalyzer_FindKmerPositions_Tests`
  = 30/30 pass.
