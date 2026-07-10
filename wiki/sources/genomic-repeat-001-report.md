---
type: source
title: "Validation report: GENOMIC-REPEAT-001 (longest repeated substring + all repeats via suffix tree)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/GENOMIC-REPEAT-001.md
sources:
  - docs/Validation/reports/GENOMIC-REPEAT-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: GENOMIC-REPEAT-001

The two-stage **validation write-up** for test unit **GENOMIC-REPEAT-001** (Repeat
Detection — Longest Repeated Substring + all repeats), validated 2026-06-15. This is
the *report* artifact that feeds one row of the [[validation-ledger]]; it records the
validator's **verdict** on both the algorithm description and the shipped code. The
repeat model is summarized in [[repetitive-element-detection]]; the two-stage
methodology is the [[validation-protocol]]. Distinct from the pre-implementation
[[genomic-repeat-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: FAIL → FIXED · End-state: CLEAN.** A real completeness
defect was found in `FindRepeats` and fully corrected this session; description was
already biologically/mathematically correct. Full unfiltered suite **6571 passed / 0
failed** (was 6570 + one new regression guard), changed files build warning-free. Test
quality gate PASS.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `GenomicAnalyzer.FindLongestRepeat(DnaSequence)` and
  `GenomicAnalyzer.FindRepeats(DnaSequence, int minLength)`.
- Sources opened: **CMU 15-451 Lecture #10 §2.1** (suffix-tree longest-repeat = deepest
  internal node with ≥ 2 leaves), **Wikipedia — Longest repeated substring problem**
  (Θ(n) deepest internal node, `$` sentinel; worked `ATCGATCGA$` → `ATCGA`), and
  **GeeksforGeeks — Suffix Tree Application 3** (corroborates the deepest-internal-node
  rule and examples).
- Formula check: LRS = deepest internal-node string-depth (matches CMU §2.1 + Wikipedia).
  The "all repeats" contract (Repeat_Detection.md §1: "every distinct substring occurring
  at least twice with length ≥ a minimum") is a sound definition grounded in CMU §2.1.
  INV-01..INV-06 are genuine properties.
- Edge cases: empty → None; no-repeat (`ACGT`) → None; overlapping occurrences counted;
  `minLength ≤ 0` clamps to 1. All sourced.
- Independent cross-check: brute-force LRS reproduced every cited value — `ATCGATCGA`→
  `ATCGA`@{0,4}; `AAAAAAAAAA`→`AAAAAAAAA`@{0,1}; `ATATATA`→`ATATA`@{0,2}; `ACGT`→none;
  `banana`→`ana`@{1,3}.
- No divergences. **Stage A = PASS.**

## Stage B — implementation (code review + cross-check)

- Code path: `FindLongestRepeat` — `GenomicAnalyzer.cs:25-37` (delegates to
  `SuffixTree.LongestRepeatedSubstring()` + `FindAllOccurrences`, correct). `FindRepeats`
  — `GenomicAnalyzer.cs:48-95` (sorts suffixes, adjacent-pair LCP).
- **Defect (fixed this session):** the original `FindRepeats` emitted only the *full*
  adjacent-pair LCP of each sorted-suffix pair, dropping every shorter repeated prefix.
  Reproduced against the built assembly: `FindRepeats("ACGTACGTTTTTACGT", 3)` returned
  only **5** substrings `{ACGT, CGT, TACGT, TTT, TTTT}` where brute-force ground truth is
  **8** (missing `ACG`, `TAC`, `TACG`) — a violation of the §1 contract. The repo's
  TestSpec/Evidence M6 row asserted the buggy 5-item set — a **code echo** (green against
  the defective code).
- **Fix:** replaced single-LCP emission with a loop over every prefix length
  `effectiveMinLength..lcpLen` of each adjacent-pair LCP, deduplicated via `HashSet`,
  positions resolved by `FindAllOccurrences`. Docstring + Repeat_Detection.md
  §4.1/§5.2/§7 updated. Re-probed: now returns the full 8-item set.
- Cross-verification (post-fix): `ACGTACGTTTTTACGT`@3 → 8 ✓; `ACGTACGT`@2 → 6 ✓;
  `ACGTACGT`@4 → `ACGT`@{0,4} ✓; `ACGTACGT`@5 → empty ✓; all LRS values match.
- `SuffixTreeGenomicsTools.FindLongestRepeat` delegate unaffected and still green.

## Findings

- **F (fixed):** `FindRepeats` completeness defect — corrected; logged in
  FINDINGS_REGISTER as a completeness defect (now resolved).
- Test-quality audit (HARD gate) PASS: M6 corrected from the code-echoed 5-item set to
  the brute-force-sourced 8-item set with exact positions; added **M8**
  (`FindRepeats("ACGTACGT", 2)` → exact 6-item set) as a completeness regression guard
  that **fails against the old implementation**; M1–M5 LRS expectations re-derived by
  brute force (exact values, not ranges); M7/S3 legitimate invariant property guards; no
  skips, no widened tolerances. Honest green: full suite 6571/0.
