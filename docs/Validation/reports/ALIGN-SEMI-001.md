# Validation Report: ALIGN-SEMI-001 — Semi-Global (Fitting / Query-in-Reference) Alignment

- **Validated:** 2026-06-12   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.SemiGlobalAlign(DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null)` — core in `SemiGlobalAlignCore`, shared `Traceback`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Exact variant (which ends are free)

The unit implements the **fitting alignment** variant of the semi-global family:

- **seq1 = query** is **fully (globally) aligned** — `RemoveGaps(aligned1) == query` (INV-3). First column initialized to `i * gap` (penalized).
- **seq2 = reference** has **both end gaps free** — leading (first row = 0) and trailing (optimum scanned over the whole last row). `RemoveGaps(aligned2)` is a substring of the reference (INV-4).

This is the Rosalind SIMS definition. Note the s/t naming is mirrored between Rosalind and Seqeron: Rosalind says "a substring of **s** against all of **t**" (s = long, t = short/motif). Seqeron names the short fully-aligned sequence `seq1`/query and the long free-end sequence `seq2`/reference. The semantics — *query fully consumed, reference both ends free* — are identical. This is **not** the four-ends-free "overlap"/SMGB variant, and the spec/Evidence explicitly and correctly state the fitting choice.

## Stage A — Description

### Sources opened & what they confirm
- **Rosalind SIMS** (https://rosalind.info/problems/sims/) — fitting alignment = "an alignment of a substring of s against all of t". DP init: first row = 0 and gaps in the short string's ends free; optimum scanned over the entire last row; recurrence is NW-style (diag/up/left max), **no zero floor**. Confirms: free gaps at both ends of the long string only; the short string is fully aligned.
- **Wikipedia: Sequence alignment** (semi-global/glocal section) — confirms semi-global does NOT penalize end gaps; explicit gene-vs-chromosome use case: "the short sequence should be globally (fully) aligned but only a local (partial) alignment is desired for the long sequence." Matches the implemented variant.
- **Wikipedia: Needleman–Wunsch** — recurrence `F(i,j)=max(F(i-1,j-1)+s, F(i-1,j)+d, F(i,j-1)+d)`, O(mn) time/space; semi-global differs only in (a) zeroing the free-start row/col and (b) taking the optimum from the free-end row/col instead of the corner.

### Formula check
Recurrence in Evidence §3.2 matches NW exactly with linear gap `d`, no `max(0,·)` floor (distinguishes it from Smith–Waterman). Init table (Evidence §3.1): fitting row0 = 0, col0 = `d·i`, optimum = `max_j F(m,j)`. All consistent with sources.

### Edge-case semantics
Query embedded / at start / at end / identical / all-mismatch (negative score, no floor) / gap-forced / null → all have defined, sourced expected behaviour. Consistent with NW recurrence and fitting definition.

### Independent cross-check (numbers)
An independent Python reimplementation of the fitting DP (init row0=0, col0=`i·gap`, NW recurrence, optimum = max of last row) reproduced **every** spec score:

| Case | query | ref | scoring (M/Mm/Gap) | Expected | Recomputed |
|------|-------|-----|--------------------|----------|------------|
| M1/M5 | ATGC | AAAATGCAAA | 1/-1/-1 | 4 | 4 |
| M4 | ATGC | ATGCAAAA | 1/-1/-1 | 4 | 4 |
| S1 | ATGCATGC | ATGCATGC | 1/-1/-1 | 8 | 8 |
| S2/MAX | ATG | ATGCCC(/CCCCC) | 1/-1/-1 | 3 | 3 |
| S3 | CCC | ATGCCC | 1/-1/-1 | 3 | 3 |
| S4 | ATGC | AATGCCC | 5/-3/-2 | 20 | 20 |
| NEG | AAAA | CCCC | 1/-1/-1 | -4 | -4 |
| OFS | ACG | AACGG | 1/-1/-1 | 3 | 3 |
| INV | GCATGCG | AAAGCATGCGAAA | 1/-1/-1 | 7 | 7 |
| MIX | AGT | AAACTAAA | 1/-1/-1 | 1 | 1 |
| GAP | ACGT | AGT | 1/-1/-1 | 2 | 2 |

**Distinguishing-from-global check (M1):** fitting score = **4**, while full global NW corner score for the same pair = **-2** (the 6 unmatched flanking reference bases penalized as -6). The +6 difference is exactly the freed end gaps, confirming the variant truly does not penalize reference end gaps.

**No-zero-floor check (NEG/MIX):** NEG = -4 (negative, so not local/SW); MIX = 1 with a mismatch in the path, confirming the NW recurrence with no `max(0,·)`.

### Findings / divergences
None. Description is mathematically correct and the exact variant (fitting) is explicitly and accurately stated.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs`:
- `SemiGlobalAlign` (412–421): null-guards both args (`ArgumentNullException.ThrowIfNull`); defaults to `SimpleDna` (Match=1, Mismatch=-1, GapExtend=-1).
- `SemiGlobalAlignCore` (423–462): `score[0,j]` left at 0 (free leading ref gaps); `score[i,0]=i*GapExtend` (query penalized → fully aligned); NW recurrence with linear `GapExtend`, no zero floor; optimum scanned over the whole last row `max_j score[m,j]` → `maxJ`.
- `Traceback` (468–534): for SemiGlobal, appends trailing reference bases (`j+1..len`) as gaps in query (free trailing gaps); walks back to top row; when `i==0, j>0` the else-branch emits leading reference gaps for free. Score returned is `score[endI, endJ] = score[m, maxJ]` (the fitting optimum, not the corner).

### Formula realised correctly?
Yes. Init (a) and optimum (b) match the validated description exactly; recurrence is plain NW with no floor. The 12 hand-computed scores all matched the running code (filtered test run), so the code reproduces the independently recomputed DP.

### Cross-verification table recomputed vs code
All 20 SemiGlobal-matching tests pass; the in-test inline DP matrices (e.g. OFS, MIX, GAP, S2, S3, S4) match my independent Python matrices. Score-from-last-row (MAX) and negative-score (NEG) are explicitly asserted to exact values.

### Variant/delegate consistency
Single canonical method; no `*Fast` or delegate variants. `Traceback` is shared with Global but branches correctly on `AlignmentType.SemiGlobal`.

### Test quality audit
Tests assert exact sourced scores (`Is.EqualTo(...)`) plus structural invariants (INV-1..5), not "no-throw" tautologies. Deterministic. Edge cases covered: embedded, start, end, identical, all-mismatch (negative), offset, mixed match/mismatch, gap-forced, null×2, score-is-last-row-not-corner.

### Findings / defects
None. No end-gap penalty leakage (not global), optimum taken from last row (not corner), and code matches the stated fitting variant (not mislabelled overlap).

## Verdict & follow-ups
- **Stage A: PASS.** Description correct; fitting variant explicitly and accurately stated; differs from global as expected.
- **Stage B: PASS.** Code faithfully realises the fitting DP; all 12 spec scores recomputed and reproduced; tests are exact and cover all Stage-A edge cases.
- **State: CLEAN** — no defect found. No code or test changes. Full suite 4461 passed / 0 failed (baseline).
- Follow-ups: none.
