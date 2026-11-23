# Validation Report: ALIGN-GLOBAL-001 — Global Alignment (Needleman–Wunsch)

- **Validated:** 2026-06-24   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.GlobalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` → `GlobalAlignCore` + `Traceback`; delegates: `GlobalAlign(string,string,ScoringMatrix?)`, `GlobalAlign(string|DnaSequence, …, CancellationToken, IProgress?)`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES. **Stage B:** PASS-WITH-NOTES. **State:** FIXED (not CLEAN: the 2026-06 verdict below is superseded).
- **Fixes (B09.md, F-numbers):** F1 `GlobalAlignAffine` (Gotoh + Flouri et al. 2015 border initialisation; the unit had no affine model and `GapOpen` was ignored), F39 one empty/null contract for every overload (empty input = all-gap NW alignment; the 2026-06 'optional tidy' became the contract), F41 `MatrixScoring` (BLOSUM62 / PAM250 / EDNAFULL, fractional gaps), F44 ends-free policies via `EndGaps`, F45 arbitrary gap-cost functions (Waterman-Smith-Beyer), F46 exact linear-space alignment (automatic above 25 M cells), F59 co-optimal global / ends-free enumeration.
- **Oracle numbers:** independent re-verification 3 x 3 936 rows vs Biopython 1.88 `PairwiseAligner` / parasail 1.3.4 (0 mismatches); F41 Biopython 1 986 / 1 986, parasail 616 / 616, EMBOSS needle / water 1 103 / 1 121 (end-gap quirk of the program); F46 4 750 forced-recursion cases identical to the full matrix; F59 Biopython `align()` sets 1 859 / 1 859 comparable rows + exhaustive enumeration.
- **Open items:** none doable. The class test count of the 2026-06 pass (13/13) is obsolete: see TestSpec ALIGN-GLOBAL-001 and B09.md for the current tests (`SequenceAligner_GlobalAlign_Tests`, `_MatrixScoring_Tests`, `_LinearSpace_Tests`, `_EndGaps_Tests`, `_GlobalCoOptimal_Tests`, `_GeneralGap_Tests`).
- The unit's TestSpec (`tests/TestSpecs/ALIGN-GLOBAL-001.md`), Evidence (`docs/Evidence/ALIGN-GLOBAL-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Needleman–Wunsch algorithm** (fetched live this session). Confirms verbatim:
  - Initialization: `F(0,j) = d·j`, `F(i,0) = d·i`.
  - Recurrence: `F(i,j) = max(F(i−1,j−1)+S(Aᵢ,Bⱼ), F(i,j−1)+d, F(i−1,j)+d)`.
  - Optimal score location: "The entry F(n,m) gives the maximum score among all possible alignments."
  - Worked example: GCATGCG vs GATTACA, match +1 / mismatch −1 / gap −1 → **optimal score 0**; one optimal alignment `GCAT-GCG / G-ATTACA`.
- Evidence doc (`docs/Evidence/ALIGN-GLOBAL-001-Evidence.md`) and TestSpec reproduce these formulas accurately, including the border row `0,−1,−2,…,−7`.

### Formula check
Single linear gap penalty `d` (true global, end gaps penalized). The TestSpec/Evidence map `d = ScoringMatrix.GapExtend`; `GapOpen` is unused by NW. Matches the Wikipedia "Advanced presentation of algorithm" equations exactly. NW (1970) and the standard linear-gap formulation are correctly described.

### Edge-case semantics
- Identical sequences → diagonal, score `n·Match` (sourced, M6).
- Completely different equal-length → all mismatches, negative score (sourced corner case, M5).
- Unequal lengths → border `d·i`/`d·j` drives the optimal (M3/M4).
- Multiple optimal alignments → one returned deterministically; explicitly allowed by source. Note: Wikipedia shows `GCAT-GCG/G-ATTACA` while TestSpec quotes `GCATG-CG/G-ATTACA` — both score 0, both valid optima.

### Independent cross-check (numbers)
Hand-recomputed the DP matrix for GCATGCG/GATTACA. Row 0 = `0 −1 −2 −3 −4 −5 −6 −7`; Row 1 (G) = `−1 1 0 −1 −2 −3 −4 −5`, reproducing the canonical Wikipedia table from the recurrence. Source-stated optimal F(7,7) = 0, confirmed by code (test M1 passes with Score == 0).

### Findings / divergences
None. Description is mathematically faithful to the primary/reference sources.

## Stage B — Implementation

### Code path reviewed
`SequenceAligner.cs`: `GlobalAlignCore` (L220–273) + `Traceback` (L468–538) — canonical; pooled flat DP buffer copied to 2D for traceback. Cancellation overload (L96–202) carries its own inline fill+traceback. `GlobalAlign` DnaSequence/string delegates L58–84, 207–218.

### Formula realised correctly?
Yes. Init `buf[i*cols]=i·GapExtend`, `buf[j]=j·GapExtend` (L239–242); recurrence `max(diag+match/mismatch, up+GapExtend, left+GapExtend)` (L253–259); returned `Score = score[m,n]` (`Traceback` `score[endI,endJ]` with endI=m, endJ=n). The cancellation overload (L119–202) computes the identical matrix and returns `matrix[m,n]`. Both paths are the validated NW linear-gap model with end gaps penalized.

### Cross-verification table recomputed vs code (test run, 13/13 pass)
| ID | Input | Expected | Result |
|----|-------|----------|--------|
| M1 | GCATGCG/GATTACA (+1/−1/−1) | 0 | PASS |
| M3 | T/ACGT, gap −1 | −2 | PASS |
| M4 | ACGT/T, gap −1 | −2 | PASS |
| M5 | AAAA/TTTT | −4 | PASS |
| M6 | ACGTACGT (identical) | 8, no gaps | PASS |
| M7 | ACGT/AGT | 2 | PASS |
| M8 | ACGT/ACGT match +5 | 20 | PASS |
| M10 | ACGT/AGT stats | M3,Mm0,G1,Id75%,Gap25% | PASS |
| M11 | symmetry score(A,B)=score(B,A) | equal | PASS |

### Variant/delegate consistency
- String overload == DnaSequence overload (M9, PASS).
- Cancellation overload == standard overload (`AlignmentProperties.GlobalAlign_CancellationOverload_SameResultAsStandard`, PASS), and same matrix/score logic by inspection.
- Property tests (FsCheck): symmetry, equal aligned length, aligned length ≥ max(input), identity→max score, determinism — all present in `Properties/AlignmentProperties.cs`.

### Numerical robustness
Integer DP; no division in the core. Scores bounded by `len·max(|Match|,|GapExtend|)`; no overflow on stated DNA-length ranges.

### Test quality audit
Assertions check exact sourced scores/alignments and statistics (not just "no throw"); INV-1/2/3 verified by reconstructing originals and recomputing the per-position score. Deterministic. Edge cases (unequal lengths, all-mismatch, identical, single deletion, empty/null) covered.

### Findings / defects
- **Minor (documented, not a defect):** the `DnaSequence` empty-input case is not guarded — `GlobalAlign(new DnaSequence(""), …)` runs the core and returns a score-0 empty-coordinate alignment, whereas the **string** overload returns `AlignmentResult.Empty`. Both are defensible (NW of empty vs empty = 0). The TestSpec frames empty/null handling as API-contract, not algorithm spec. No correctness impact on the NW model.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- No code changes. Full GlobalAlign test class 13/13 green; corroborating property/cancellation tests green.
- No defects logged by the first pass (its "optional future tidy" about empty inputs was REFUTED by audit A1-1 / F39: see the correction below).

## 2026-11 audit round 1 correction (A1-1, B09 F39)

The "Minor (documented, not a defect)" empty-input note and the "optional future tidy" above are **superseded**: the string overloads' silent `AlignmentResult.Empty` for an empty side was the defect (the Needleman-Wunsch optimum of `""` vs `ACGT` is the all-gap alignment, border `F(0,j) = d·j`), the opposite of the suggested tidy. All `GlobalAlign`/`GlobalAlignAffine` overloads now follow one contract: null -> `ArgumentNullException`, empty side -> all-gap alignment, both empty -> score 0 / empty rows / coordinates 0, 0, -1, -1. Biopython 1.88, parasail 1.3.4, EMBOSS needle 6.6 and scikit-bio reject zero-length input, so the empty case has no tool oracle (see Evidence). Tests: `EmptyInput_AllOverloads_AreTheAllGapNeedlemanWunschAlignment`, `EmptyInput_LinearScore_IsGapTimesLength`, `NullString_ThrowsArgumentNullException_OnEveryStringOverload`.

## 2026-11 audit round 1 (A1-2, B09 F41)

Substitution-matrix scoring and fractional gap scores for global/local/semi-global, linear and affine: additive `MatrixScoring` overloads (string, ProteinSequence, DnaSequence) on the existing DP cores (no second DP), `AlignmentResult.ExactScore`, MCP `matrix` parameter on the six pairwise tools. Oracle: Biopython 1.88 1 986/1 986, parasail 1.3.4 616/616, EMBOSS needle/water 6.6.0 1 103/1 121 (EMBOSS-internal end-gap quirk). Details: B09.md F41, algorithm doc §5.4.

## 2026-11 audit round 1 (A1-3, A1-5, B09 F44/F45)

`SemiGlobalAlign(…, EndGaps.All)` is EMBOSS needle's default ends-free alignment and `EndGaps.None` equals this unit's global alignment (F44). General gap-cost functions (Waterman-Smith-Beyer, any w(k)): `GlobalAlignGeneralGap` (F45; Biopython callable `gap_score` oracle, 960 global/local rows + 480 ends-free rows); `Global_Alignment_Needleman_Wunsch.md` §5.5.

## 2026-11 audit round 1 (A1-6, B09 F46)

Exact linear-space alignment (Hirschberg 1975 / Myers-Miller 1988 / Huang-Miller 1991 scheme): every pairwise aligner switches automatically above 25 000 000 cells (identical score, rows, coordinates, ties), explicit `*LinearSpace` methods for any size, MCP automatic. Oracles: Biopython 1.88 / parasail 1.3.4 (1 280 pairs), EMBOSS `stretcher` / `matcher` 6.6.0 (160 + 160 rows; large-pair scores equal up to 50 kb), full-matrix differential (4 750 forced-recursion cases + kb pairs + the automatic switch). 20 kb x 20 kb affine: 111 MB peak (4.8 GB full); 100 kb x 100 kb: 185 MB. Memory guards for the non-affine general-gap aligners (4 M cells) and the profile MSA DP (50 M cells). Algorithm document: `Linear_Space_Alignment.md`.

## 2026-11 audit round 1 (A1-7b, B09 F59)

The "one optimal alignment only" limitation of the global problem is lifted: `GlobalAlignCoOptimal` / `GlobalAlignAffineCoOptimal` / `SemiGlobalAlign*CoOptimal` (+ MCP `global_align_co_optimal`) enumerate every co-optimal alignment (capped) from the shared fill / border helpers of the single aligners - Biopython 1.88 set equal on 1 859 / 1 859 comparable oracle rows (all 16 end-gap policies, linear/affine, scalar and BLOSUM62), exhaustive enumeration equal on tiny pairs. See the algorithm doc §5.6.

## Audit r2 addendum (A2-37, A2-43, B09 F80 / F85): one scoring contract, `null` scoring overloads

Scalar `ScoringMatrix` and `MatrixScoring` now share one contract: gap scores <= 0 (and > `int.MinValue`), path score below 2^28 (`ArgumentException` instead of a silently wrapped score; e.g. `GlobalAlign("A"x50, "A"x50, new ScoringMatrix(1e9, -1, -2, -1))` returned 1 999 999 904 for the true 5e10); sum-of-pairs totals accumulate in 64 bits. `null` for the scoring argument compiles on every overload pair (`[OverloadResolutionPriority(1)]`). The MCP pairwise tools delegate empty input to the library (all-gap alignment). Tests: `SequenceAligner_ScoringContract_Tests`, `SequenceAligner_InputContract_Tests`, MCP `EmptyInputContractTests`.

## Audit r2 addendum (A2-28, B09 F121 / F122): stretcher and needle options

`SequenceAligner.StretcherAlign` (Myers-Miller as the EMBOSS program, 689 / 689 pairs equal) and `SequenceAligner.NeedleAlign` (the program's single-precision recurrences with `-endweight` / `-endopen` / `-endextend`, 887 / 887 pairs equal); MCP `stretcher_align`, `needle_align` (audit A3-17, B09 F141). Audit A3-14b (F143): the "at least 90 % reaches the plain optimum" tolerance is replaced by an independent score-only program of the needle recurrence (492 / 492 dyadic rows; the one sampled row below the plain optimum is id 234) and stretcher equals the plain optimum on 300 / 300.

## Audit r2 addendum (A2-29, B09 F123): dual affine gaps

`GlobalAlignDualAffine` / `LocalAlignDualAffine` (ksw2-style two-piece affine gaps, optional band), MCP `dual_affine_align`; 1 600 oracle cases equal.

## Audit r2 addendum (A2-28c / A2-29b, B09 F124 / F125): needleall, ksw2

`SequenceAligner.NeedleAll` (+ `SimilaritySubHeader`, the needle `-nobrief` lines) and `Ksw2Aligner` (ksw_extz / extz2_sse / extd / extd2_sse: band, Z-drop, end bonus, N scoring, CIGAR); MCP `needleall_align`, `ksw2_extend`. needleall 260 / 260 files equal; ksw2 exact cases 1 320 / 1 320 equal an independent program; the real-program oracle is a Docker package (A2-29b stays DOABLE until its output is integrated).

- **Audit A4-3 (F153):** the EMBOSS reader drops gap characters and every other non-letter (except `*`) from plain input sequences; one shared `EmbossSequenceReader` (Infrastructure) serves needle, stretcher, needleall: 240 oracle rows each equal (`emboss_gap_oracle.tsv`).

- **Audit r5 / r6 (A5-7, A6-2, A6-4; F172, F177, F179):** `needle` / `stretcher` / `merger` / `needleall` read the sequence-input qualifiers (`NeedleProgram`, `StretcherProgram`, `MergerProgram`, `NeedleAll`; 180 oracle rows each: ranges, `-sreverse`, case, type overrides, the refused combinations); the qualifier text resolves like the ACD library (`needle a.fa b.fa -sbegin 3` = `-sbegin2 3`, `-sbeg1`, `-sbegin1=5`: 2 080 + 640 rows); a null member / `Name` / `Residues` of a set is an `ArgumentNullException` naming the parameter and the index.

- **Audit R14-1 (F288):** `AlignmentResult.End*` is 0-based inclusive in every producer (Needle / Water / Stretcher / dual affine were exclusive, the EMBOSS formats of the MCP core aligners were wrong: markx10 lengths, negative offsets, double-counted fitting flanks, local results printed as needle reports); `EmbossAlignedSequence.Flanks` derives offset / offend from the rows; contract test over 19 producers; markx10 / pair / sam MCP tests equal `water`; the program oracles (needle / water / stretcher / matcher / dual affine) stay equal.

- **Audit R14-15 (F294):** every public pairwise aligner has an explicit trailing `CancellationToken` overload (163 generated by `generate_cancellable_overloads.py`; they enter the ambient scope the dynamic programs poll); `IProgress<double>` stays on the linear `GlobalAlign`.

- **Audit R14-8 (F301):** all 30 named matrices of Biopython 1.88 (fractional ones exact), SCHNEIDER codon matrix, `SubstitutionMatrix.Parse` / `ToText`, matrix text in the MCP `matrix` parameter, U / O / J = the `BLASTP` matrix (blastp oracle); 29 matrices cell-equal and 812 scores equal to Biopython.

- **Audit R14-9 (F302):** `AlignmentResultCoordinates` (`ToBioCoordAlignment`, `ToCigar` standard / extended / with soft clips, `ToMd`, `ToCs`, `Format`) and the `format` names bed / psl / chain / maf / sam / exonerate / cigar / md / cs of the pairwise MCP tools; Biopython 2 100 texts, parasail 667 and minimap2 120 records equal.

- **Audit R14-10 (F303):** `PairwiseAligner` (FOGSAA port, wildcard, strand, presets, the twelve gap scores), `AlignmentOperations` + `AlignmentCounts` (every `Bio.Align.Alignment` operation of the audit and the rest of its public members), the co-optimal enumerators with direction-specific gap scores, MCP `pairwise_aligner` / `alignment_counts` / `alignment_operation`; 1 873 + 700 Biopython 1.88 oracle tests; remainder R14-10b registered.

- **Audit R14-11 (F304):** the parasail 1.3.4 score API (`ParasailAligner`, `ParasailProfile`, MCP `parasail_score`): score-only with end coordinates in O(min(m, n)) memory, statistics without a traceback, tables, last row / column, one query against many references, for nw / sw / 11 sg variants; 3 120 parasail results equal.

- **Audit R14-12 (F305):** banded (Ukkonen / Chao-Pearson-Miller) and X-drop alignment: `BandedAlign`, `GlobalAlign[Affine]Banded`, `LocalAlign[Affine]Banded`, `SemiGlobalAlign[Affine]Banded`, `GlobalAlignBandedExact` (proven optimum), MCP `banded_align`; equal to the full programs when the full alignment lies in the band; 30 000 x 30 000 affine global 98 s -> 0.9 s.

- **Audit R14-10b (F307):** `PairwiseAligner` Score / Align with gap-score functions (Waterman-Smith-Beyer) and gap scores above zero, `CodonAligner`, `MultipleSeqAlignment` / `SeqRecord`, MCP `codon_aligner` / `multiple_seq_alignment` / `gapFunctionsJson`; 1 501 + 1 500 + 1 200 Biopython 1.88 oracle cases (7 030 operations) equal; 100 inputs on which the reference reads uninitialised memory are refused.

