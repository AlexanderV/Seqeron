# Validation Report: KMER-POSITIONS-001 — K-mer Positions

- **Validated:** 2026-06-16   **Area:** K-mer
- **Canonical method(s):** `KmerAnalyzer.FindKmerPositions(string sequence, string kmer)` → `IEnumerable<int>`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened (retrieved this session)

1. **Rosalind BA1D — "Find All Occurrences of a Pattern in a String"** — https://rosalind.info/problems/ba1d/ (WebFetch, 2026-06-16). Confirmed verbatim:
   - Output is "all starting positions in *Genome* where *Pattern* appears as a substring. **Use 0-based indexing.**"
   - Overlapping occurrences are all reported — the problem's own example states "ATA occurs three times in CGATATATCC​ATAG" (overlapping matches counted separately).
   - Sample: Pattern `ATAT`, Genome `GATATATGCATATACTT` → Output `1 3 9`.
2. **Wikipedia — "k-mer"** — https://en.wikipedia.org/wiki/K-mer (WebFetch, 2026-06-16). Confirmed:
   - A k-mer is a "substring of length k contained within a biological sequence."
   - Number of overlapping k-mers in a length-L sequence = **L − k + 1**.
   - Worked example `AGAT`: 2-mers are `AG`, `GA`, `AT` (4 − 2 + 1 = 3).

### Formula check

`Occ(P,T) = { i ∈ [0, L−k] : T[i..i+k) = P }`, reported ascending, 0-based. This matches Rosalind BA1D (0-based, all overlapping starts) and the Wikipedia L−k+1 window count exactly. The algorithm doc §2.2 and INV-01..04 reproduce this faithfully.

### Edge-case semantics check

- Pattern longer than text → L−k+1 ≤ 0 → empty (Wikipedia count formula). ✔ sourced.
- Pattern absent → empty (only matching starts reported, BA1D). ✔ sourced.
- Pattern equals whole sequence → `[0]` (one occurrence). ✔ derived from BA1D predicate.
- Self-overlap (`AA` in `AAAA` → `[0,1,2]`) → BA1D overlapping rule. ✔ sourced.
- Null/empty input → empty; case-insensitive matching — **not source-mandated**; recorded honestly as repository-convention assumptions (Evidence §Assumptions 2 & 3). Acceptable documented divergence.

### Independent cross-check (numbers retrieved/computed this session)

Independent overlapping-search reference implementation in Python (`str.find` loop with `i+1` advance), run this session:

| Case | Pattern | Text | Python reference | External-source value |
|------|---------|------|------------------|-----------------------|
| M1 | ATAT | GATATATGCATATACTT | `[1, 3, 9]` | Rosalind BA1D `1 3 9` ✔ |
| M2 | AA | AAAA | `[0, 1, 2]` | BA1D overlapping rule ✔ |
| M3 | ana | banana | `[1, 3]` | hand-computed ✔ |
| M4 | AG/GA/AT | AGAT | `[0]`/`[1]`/`[2]` | Wikipedia k-mer ✔ |
| M5 | GG | ATATAT | `[]` | BA1D (no match) ✔ |
| M6 | ATAT | ATATATAT | `[0, 2, 4]` | hand-computed (L−k=4) ✔ |
| S1 | ACGT | AC | `[]` | L−k+1 ≤ 0 ✔ |
| S2 | ACGT | ACGT | `[0]` | whole-sequence ✔ |
| — | ATA | CGATATATCCATAG | `[2, 4, 10]` | Rosalind's own "occurs three times" example ✔ |

The Rosalind in-problem overlapping example (`ATA` in `CGATATATCC​ATAG` → three occurrences) was reproduced exactly by the independent implementation, confirming the overlapping convention beyond the published sample.

### Findings / divergences

None affecting correctness. Two non-sourced conventions (case-insensitive matching, null/empty → empty) are explicitly recorded as assumptions in the Evidence doc; they do not alter any all-uppercase evidence example. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs:432-447`.

```csharp
if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(kmer)) yield break;
var seq = sequence.ToUpperInvariant();
var km  = kmer.ToUpperInvariant();
for (int i = 0; i <= seq.Length - km.Length; i++)
    if (seq.AsSpan(i, km.Length).SequenceEqual(km.AsSpan()))
        yield return i;
```

### Formula realised correctly?

Yes. The loop bound `i <= L − k` is exactly the L−k+1 candidate window range; every candidate is tested (overlapping starts all emitted); positions are yielded in ascending order by construction (no sort needed); `SequenceEqual` over `ReadOnlySpan<char>` is an exact, allocation-light substring compare. 0-based by C# string indexing. Matches the validated description verbatim.

### Cross-verification table recomputed vs code

The full test suite exercises every row of the Stage-A table against the actual code; all pass (see below). Values match the external-source column above exactly — no code-echo: M1 and M2 assertions are written so a 1-based or non-overlapping implementation would FAIL (`[1,3,9]` not `[2,4,10]`; `[0,1,2]` not `[0,2]`).

### Variant/delegate consistency

Single canonical method; no overloads (`grep` over `src/` confirms only the one definition plus an MCP wrapper `AnalysisTools.cs:128` that simply forwards). Nothing else to reconcile.

### Numerical robustness

No arithmetic precision concerns (integer index loop). `seq.Length - km.Length` can be negative when k > L, making the loop body never execute → empty, matching S1. Null/empty guarded before any `.Length` access. No overflow/div-by-zero paths.

### Test quality audit (HARD gate)

- **Sourced, not echoed:** every assertion is `Is.EqualTo(exact)` traced to Rosalind/Wikipedia/hand-computation; M1 & M2 comments document that a wrong (1-based / non-overlapping) impl is rejected — confirmed these would fail such an impl.
- **No green-washing:** no Greater/AtLeast/Contains/range where an exact value is known; no widened tolerance; no skip/ignore/comment-out; no expected value tuned to actual output.
- **Coverage:** all 11 cases (M1–M6, S1–S3, C1–C2) cover the single public method, every Stage-A branch (happy/overlapping/absent/longer-than-text/equals-whole/case-insensitive/null seq/empty seq/null kmer/empty kmer), the L−k+1 boundary, and the ascending-order invariant. No untested public surface.
- **Honest green:** full unfiltered suite `dotnet test` → **Failed: 0, Passed: 6607, Skipped: 0**; `dotnet build` → 0 errors (4 pre-existing warnings, all in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, none in the KMER-POSITIONS file).

**Test-quality gate: PASS.**

### Findings / defects

None. **Stage B: PASS.**

## Verdict & follow-ups

- **Stage A: PASS** — description matches Rosalind BA1D and Wikipedia k-mer exactly; two non-sourced conventions are honestly recorded as assumptions.
- **Stage B: PASS** — code realises the validated formula verbatim; tests assert exact externally-sourced values and would reject a deliberately-wrong implementation.
- **End-state: CLEAN** — no defect found; no code or test change required.
- No defects logged.
