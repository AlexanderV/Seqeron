# Validation Report: ALIGN-SEMI-001 — Semi-Global Alignment (Fitting / Query-in-Reference)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` → `SemiGlobalAlignCore` (`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:412-462`), shared `Traceback` (`:468-538`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Rosalind SIMS** (https://rosalind.info/problems/sims/): fitting alignment = "an alignment of a substring of *s* against all of *t*", i.e. the short sequence is **fully** aligned and the long sequence has **free leading and trailing end gaps**. DP init: first row = 0 (free leading gaps in the long seq), first column = decreasing gap penalties (short seq fully aligned), optimum read from the extreme of the row/column that fully consumes the short seq, traceback back to the first row. Match +1 / mismatch −1.
- **Wikipedia: Sequence alignment** — "semi-global or 'glocal'… search for the best possible partial alignment"; the canonical use case is "one sequence is short (a gene) and the other very long (a chromosome)… the short sequence should be globally (fully) aligned but only a local (partial) alignment is desired for the long sequence." Confirms free end gaps on the long sequence only.
- **Wikipedia: Needleman–Wunsch** — recurrence `F(i,j)=max(F(i−1,j−1)+s, F(i−1,j)+d, F(i,j−1)+d)`, **no zero floor** (so scores may be negative), O(mn).

### Formula / convention check
The TestSpec/Evidence model (§3.1–3.3) is exactly the fitting variant:
- First row `F(0,j)=0` → free leading reference gaps.
- First column `F(i,0)=i·d` → query (short seq) fully aligned.
- Recurrence = NW, no zero floor.
- Optimum = `max_j F(m,j)` (free trailing reference gaps).

Role-labelling note (not a defect): Rosalind labels the **long** string *s* and the short motif *t*; Seqeron labels seq1 = query (short, fully aligned) and seq2 = reference (long, free end gaps). The init mapping (first **row** = 0 free on the long seq, first **column** penalized on the short seq, max over the **last row**) is the transpose-consistent realisation of SIMS — semantics identical: short fully aligned, long ends free. Verified consistent.

### Edge-case semantics
All §1.3/§4 corner cases have sourced, defined behaviour: query=ref → global=fitting (m=n); all mismatches → exact negative score (no zero floor); query longer than ref → forced reference gap (NW up-move); null → `ArgumentNullException` (.NET convention). All standard and correct.

### Independent cross-check (hand computation)
- **M1** query=ATGC in ref=AAAATGCAAA: query appears exactly at ref[3..6]; 4 matches, leading/trailing AAA free → **score 4**. ✓
- **GAP** query=ACGT (m=4) vs ref=AGT (n=3): ACGT / A-GT = 3 match +1 ref-gap = 3−1 = **2**. ✓
- **MIX** query=AGT vs ref=AAACTAAA: best window ACT = A(+1) G/C(−1) T(+1) = **1**. ✓
- **MAX** query=ATG vs ref=ATGCCC: fitting `max_j F(3,j)=F(3,3)=3`, whereas global `F(3,6)=0` — confirms optimum must be max of last row, not bottom-right. ✓

### Findings / divergences
None. Description is mathematically and biologically correct; matches Rosalind SIMS fitting alignment and the Wikipedia semi-global use case. The deliberate scope decision (single fitting member of the semi-global family; linear gap, GapOpen unused) is documented and sourced.

## Stage B — Implementation

### Code path reviewed
`SemiGlobalAlignCore` (`SequenceAligner.cs:423-462`) and `Traceback` (`:468-538`).

### Formula realised correctly? (evidence)
- First column `score[i,0]=i*GapExtend` (`:431-432`) → query fully aligned. ✓
- First row left at default 0 (`:428`) → free leading reference gaps. ✓
- Recurrence `Max(diag, up, left)`, no zero floor (`:441-445`). ✓
- `maxJ = argmax_j score[m,j]` (`:449-459`), returned score = `score[m, maxJ]` (`Traceback` returns `score[endI,endJ]`, `:532`) = `max_j F(m,j)` → INV-5. ✓
- Trailing reference suffix `seq2[maxJ..n-1]` appended as gaps in seq1, in reverse so the post-reverse output reproduces the reference (`:485-492`); leading reference handled by the `i==0, j>0` left-move branch (`:509-520`). ✓

### Cross-verification table recomputed vs code (tests run)
| Case | Input | Hand value | Code | Match |
|------|-------|-----------|------|-------|
| M1 | ATGC / AAAATGCAAA | 4 | 4 | ✓ |
| S1 | ATGCATGC / ATGCATGC | 8 | 8 | ✓ |
| S2 | ATG / ATGCCCCC | 3 | 3 | ✓ |
| S3 | CCC / ATGCCC | 3 | 3 | ✓ |
| S4 | ATGC / AATGCCC (5,−3,−2) | 20 | 20 | ✓ |
| NEG | AAAA / CCCC | −4 | −4 | ✓ |
| MAX | ATG / ATGCCC | 3 (not 0) | 3 | ✓ |
| OFS | ACG / AACGG | 3 | 3 | ✓ |
| MIX | AGT / AAACTAAA | 1 | 1 | ✓ |
| GAP | ACGT / AGT | 2 | 2 | ✓ |
| INV | GCATGCG / AAAGCATGCGAAA | 7 | 7 | ✓ |

### Variant/delegate consistency
Single canonical method. MCP wrapper `AlignmentTools.SemiGlobalAlign` delegates directly to `SequenceAligner.SemiGlobalAlign` (no separate logic). Consistent.

### Test quality audit
17 canonical tests (`SequenceAligner_SemiGlobalAlign_Tests.cs`) + 1 property test (`Properties/AlignmentProperties.cs`). All assert **exact** sourced scores with hand-computed DP last rows documented inline; invariants INV-1..5 validated; edge paths (negative score, gap, mixed match/mismatch, offset, max-of-last-row vs bottom-right, null) covered. No "no-throw"/tautological assertions. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS** → **State: CLEAN**. No code changed.
- All 17 canonical + property tests pass; full `SequenceAligner_SemiGlobalAlign_Tests` filter = 17/0 and `AlignmentProperties` filter = 22/0.
- Optional (non-blocking) documentation nicety: the unit is precisely the **fitting** member of the semi-global family (both reference ends free; the short/query seq fully aligned), not the four-ends-free SMGB variant — already correctly stated in TestSpec §6 and Evidence §2.3. No action required.
