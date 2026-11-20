# Validation Report: ALIGN-SEMI-001 — Semi-Global Alignment (Fitting / Query-in-Reference)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.SemiGlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` → `SemiGlobalAlignCore` (`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:412-462`), shared `Traceback` (`:468-538`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS. **State:** FIXED (not CLEAN: the 2026-06 verdict below is superseded).
- **Fixes (B09.md, F-numbers):** F3 `SemiGlobalAlignAffine` (the unit had no affine model; TestSpec NEG / INV-4 corrected), F39 contract of the string overload, F40 `Start/EndPosition2` = the fitted interval (were the full-input bounds), F41 `MatrixScoring`, F44 overlap and the 16 ends-free policies (`EndGaps`), F45 general gap costs, F46 exact linear-space fitting / ends-free.
- **Oracle numbers:** Biopython 1.88 `align().score` 3 590 / 3 590 for the 16 policies, co-optimal membership 2 050 / 2 050, parasail exact variants equal (end-flag variants lower on 23 rows: parasail skips border end cells), needle default `-endweight N` 277 / 287 (never above), F40 823 / 823 enumerated Biopython sets contain the fitted interval + brute-force oracle on 600 pairs.
- **Open items:** none doable. 'No code changed' of the 2026-06 pass is obsolete.
- The unit's TestSpec (`tests/TestSpecs/ALIGN-SEMI-001.md`), Evidence (`docs/Evidence/ALIGN-SEMI-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

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

## 2026-11 audit round 1 correction (A1-1, A1-4, B09 F39/F40)

"No code changed" is superseded by the audit-round changes: (1) `SemiGlobalAlign(string, string, ...)` was added (null -> `ArgumentNullException`, empty not special-cased); (2) `StartPosition2`/`EndPosition2` are now the **fitted interval** of the reference (first/last reference index consumed between the free end gaps; empty interval: `End2 = Start2 - 1`) instead of the full-input bounds, for `SemiGlobalAlign` and `SemiGlobalAlignAffine`. Oracle: 4090 random/edge cases (1382 semi-global rows), interval in the Biopython 1.88 co-optimal interval set in 823/823 fully enumerated rows, parasail 1.3.4 score equal, brute-force fitting definition equal (600 pairs).

## 2026-11 audit round 1 (A1-2, B09 F41)

Matrix-scored `SemiGlobalAlign`/`SemiGlobalAlignAffine`: see ALIGN-GLOBAL-001 report and B09.md F41.

## 2026-11 audit round 1 (A1-3, B09 F44; A1-5, F45)

The overlap / all-ends-free members of the semi-global family are implemented (`EndGaps`, 16 policies, linear + affine, scalar + matrix scoring; Biopython `left/right_insertion/deletion_score`, parasail `sg_*`, EMBOSS needle default `-endweight N`; 3 590 oracle rows, `SequenceAligner_EndGaps_Tests` 14 tests + MCP `endGaps`), so the "simplified subset / not implemented" statements of Stage A/B above are superseded; `Semi_Global_Alignment.md` status Simplified -> Production. General gap-cost functions: `SemiGlobalAlignGeneralGap` (F45). See B09.md F44/F45.

## 2026-11 audit round 1 (A1-6, B09 F46)

Exact linear-space alignment (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991 scheme): every pairwise aligner switches automatically above 25 000 000 cells (identical score, rows, coordinates, ties), explicit `*LinearSpace` methods for any size, MCP automatic. Oracles: Biopython 1.88 / parasail 1.3.4 (1 280 pairs), EMBOSS `stretcher` / `matcher` 6.6.0 (160 + 160 rows; large-pair scores equal up to 50 kb), full-matrix differential (4 750 forced-recursion cases + kb pairs + the automatic switch). 20 kb x 20 kb affine: 111 MB peak (4.8 GB full); 100 kb x 100 kb: 185 MB. Memory guards for the non-affine general-gap aligners (4 M cells) and the profile MSA DP (50 M cells). Algorithm document: `Linear_Space_Alignment.md`.
- **Audit r3 (A3-10, F137):** per-direction / per-end gap scores (`GapScoreProfile`, Biopython's twelve properties) on the existing cores; the "out of scope" statement of the algorithm doc is removed (see Semi_Global_Alignment.md §5.5); 998 / 1000 Biopython rows equal, the 2 others a documented Biopython linear-NW corner.
