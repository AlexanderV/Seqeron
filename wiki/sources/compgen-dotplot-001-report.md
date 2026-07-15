---
type: source
title: "Validation report: COMPGEN-DOTPLOT-001 (dot plot — word-match / k-tuple dot matrix, ComparativeGenomics.GenerateDotPlot)"
tags: [validation, comparative-genomics, governance]
doc_path: docs/Validation/reports/COMPGEN-DOTPLOT-001.md
sources:
  - docs/Validation/reports/COMPGEN-DOTPLOT-001.md
source_commit: 37c54d6df345672ac216015eda5a1639544b4b01
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: COMPGEN-DOTPLOT-001

The two-stage **validation write-up** for test unit **COMPGEN-DOTPLOT-001** — **dot plot
generation** (word-match / k-tuple dot matrix): emit a dot at `(i, j)` wherever a
length-`wordSize` word starting at position i in sequence 1 exactly matches the word starting at
position j in sequence 2. Validated 2026-06-16. This is the *report* artifact that feeds one row
of the [[validation-ledger]]; it records the validator's independent **verdict** on both the
algorithm description (Stage A) and the shipped code (Stage B), within the wider
[[validation-and-testing]] campaign. The algorithm, its parameters, invariants, worked oracles
and corner cases are synthesized in the concept [[dot-plot-word-match]]; [[test-unit-registry]]
tracks the unit. Distinct from [[compgen-dotplot-001-evidence]] — the pre-implementation evidence
artifact sourced from `docs/Evidence/` — this page is the independent two-stage re-validation
verdict.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End state: CLEAN.** No implementation defect and no
code change required. Two **test-quality** issues were found and fixed **in-session**: one weak
assertion strengthened and one coverage gap closed (see below). Honest green: the FULL unfiltered
suite `dotnet test` = **6606 passed, 0 failed, 0 skipped** (the one Skipped is the pre-existing
benchmark guard, not a dot-plot test); `dotnet build` 0 errors.

## Canonical method & source under test

- `ComparativeGenomics.GenerateDotPlot(string sequence1, string sequence2, int wordSize = 10, int stepSize = 1)`
  in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs:1169-1207`.
- Eager validation: `wordSize <= 0` / `stepSize <= 0` → `ArgumentOutOfRangeException`, thrown
  **before** the iterator (surfaces immediately, not on first `MoveNext`).
- Iterator: null/empty or shorter-than-word ⇒ `yield break`. Builds a `SuffixTree` over
  `ToUpperInvariant(seq2)`, slides the word `seq1.Substring(i,w).ToUpperInvariant()` for
  i = 0, step, … ≤ n−w, and yields `(i, j)` for every occurrence `j` from `FindAllOccurrences`.
- `SuffixTree.FindAllOccurrences` (`SuffixTree.Search.cs:81`) collects **all** leaves under the
  matched node ⇒ **all overlapping** start positions (verified: every occurrence, not just first).
- Default `wordSize` constant = 10 (`DefaultDotPlotWordSize`). Single public method — no
  delegates / `*Fast` variants.

## Stage A — description (algorithm faithfulness)

Confirmed against sources **independently retrieved this session**: **Huttley TIB Dotplot**
(match rule `X[i]==Y[j]` at k=1, worked example `AGCGT` vs `AT` → 1-based (1,1),(5,2) ⇒ 0-based
{(0,0),(4,1)}; matrix indexed [seq1-index, seq2-index] confirms x=seq1/y=seq2), **EMBOSS `dottup`
manual + manpage** (exact word/tuple match, no scoring matrix; wordsize default = 10; longer→less
noise/less sensitive, shorter→more sensitive/more noise), and **Wikipedia — Dot plot
(bioinformatics)** (dot on same-location match; main diagonal = self-alignment; off-diagonal lines
= similar/repetitive patterns; tuple noise reduction; cites **Gibbs & McIntyre 1970**).

The implementation realises the match relation
`D = { (i, j) : A[i..i+w−1] = B[j..j+w−1] }, 0 ≤ i ≤ n−w, 0 ≤ j ≤ m−w` (case-insensitive) —
exactly the `dottup` exact-word rule and, at w=1, the single-residue rule `X[i]==Y[j]`.
**Independent cross-checks** hand-derived from the sourced rule: `AGCGT` vs `AT`, k=1 →
{(0,0),(4,1)}; `ACGTACGT` self, w=4 → {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} (all overlapping
occurrences); `ACGT` self, w=1 → exactly {(0,0),(1,1),(2,2),(3,3)} (full main diagonal, nothing
else). Edge semantics all confirmed: word > sequence ⇒ empty; null/empty ⇒ empty; self ⇒ full
main diagonal; disjoint alphabets ⇒ empty; default wordSize = 10.

**Stage A = PASS.** Three documented presentation/contract choices are reasonable and labelled as
decisions, not sourced biology: axis orientation (x=seq1, y=seq2); case-insensitive folding via
`ToUpperInvariant`; and non-positive `wordSize`/`stepSize` throwing `ArgumentOutOfRangeException`
(INV-5, a sibling-convention contract — an undefined window, not a biology claim).

## Stage B — implementation & test-quality audit (HARD gate)

Formula realised correctly (direct enumeration of `D`; case folding; stepSize samples seq1 only;
default constant = 10). Ten cross-verification cases recomputed vs code all matched, including
`AAAA` vs `CCCC` → ∅; `ACG` vs `ACGT`, w=4 → ∅; w=0/−1 & step=0/−1 → throws AOORE;
`ACGTACGT` self w=4 step=4 → {(0,0),(0,4),(4,0),(4,4)}; `acgtacgt` vs `ACGTACGT` w=4 →
case-insensitive same set; `ACGTACGTAC` self default w=10 → {(0,0)} (locks default=10).

Two **test-quality** fixes applied in-session (no product-code change):

1. **Weak assertion strengthened (M3).** `GenerateDotPlot_SelfComparison_ContainsMainDiagonal`
   used `Is.SupersetOf` where the **exact** match set is known (`ACGT` has distinct residues ⇒
   the set is *precisely* the main diagonal). A superset assertion would green-wash an
   implementation emitting spurious off-diagonal dots. Tightened to
   `Is.EquivalentTo {(0,0),(1,1),(2,2),(3,3)}` (the sourced value).
2. **Coverage gap closed (S3).** New `GenerateDotPlot_DefaultWordSize_IsTenAndMatchesFullWord`
   exercises the default-parameter path and locks the documented EMBOSS `dottup` default
   wordSize=10 with the hand-verified exact set {(0,0)} for a length-10 self-comparison.
   Previously no test exercised the default value.

Every Stage-A branch is now covered (match rule, self-diagonal, disjoint, word>seq, null/empty,
invalid wordSize, invalid stepSize, stepSize sampling, case-insensitivity, default wordSize, and
an O(n×m)/multiples-of-step property test C1). Remaining assertions use exact
`Is.EquivalentTo` / `Is.Empty` / exact exception type — sourced values, not code echoes.

## Findings & follow-ups

- **No implementation defect; State CLEAN.** One weak test assertion strengthened (M3) and one
  default-path coverage gap closed (S3), both in-session; no product-code change.
- **No open follow-ups.** Divergences from a bare `dottup` are documented decisions (axis
  orientation, case-insensitive folding, non-positive-window throws), not defects.
