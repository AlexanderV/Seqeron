---
type: source
title: "Validation report: ANNOT-REPEAT-001 (repetitive element detection & classification — GenomeAnnotator.FindRepetitiveElements / ClassifyRepeat)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-REPEAT-001.md
sources:
  - docs/Validation/reports/ANNOT-REPEAT-001.md
source_commit: f11e8bc7feeba4d997051309d93f661db1f53382
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-REPEAT-001

The two-stage **validation write-up** for test unit **ANNOT-REPEAT-001** (repetitive element
detection and classification — tandem repeats, inverted repeats, repeat-class assignment),
validated 2026-06-15. This is the *report* artifact that feeds one row of the
[[validation-ledger]] (row **58**); it records the validator's **verdict** on both the algorithm
description and the shipped code. The algorithm and its three sub-problems are summarized in
[[repetitive-element-detection]]; the two-stage methodology is the [[validation-protocol]] under
[[validation-and-testing]]. Distinct from the pre-implementation
[[annot-repeat-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN — one real defect found and fully fixed
this session** (code + tests + Evidence + TestSpec). Full **unfiltered** suite
**6566 passed / 0 failed** (1 pre-existing unrelated benchmark skipped); `dotnet build` 0 errors,
changed files build warning-free. Logged as finding **A17** in `FINDINGS_REGISTER.md`.

## Canonical methods validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`:

- `FindRepetitiveElements(string, int minRepeatLength, int minCopies)` (`:652-792`) — tandem via
  `FindTandemRepeats` + `IsPrimitive`; inverted via `FindInvertedRepeats`.
- `ClassifyRepeat(string, IReadOnlyDictionary<string,string>)` + `IsSimpleRepeat` (`:808-869`).

## Stage A — description (algorithm faithfulness)

Confirmed against five sources opened live: **Wikipedia "Tandem repeat"** (head-to-tail adjacency,
worked example `ATTCG ATTCG ATTCG`, STR < 10 nt, minisatellite 10–60 nt), **NIH/genome.gov
glossary** (the authoritative "two or more" head-to-tail minimum → `minCopies >= 2`),
**Wikipedia "Inverted repeat"** (`5'-TTACGnnnnnnCGTAA-3'`, spacer any length incl. zero, gap-0 ⇒
palindrome), **IUPACpal / Hampson et al. 2021** (BMC Bioinformatics 22:51; formal `W W̄ᴿ` / gapped
`W G W̄ᴿ` with `|G| ≥ 0`, imperfect within `δ_H ≤ k`), and **RepeatMasker** (class vocabulary
SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Low_complexity/RNA/Unknown; classify by best homology).

Hand-checked reverse complements confirmed: revcomp(`GAA`)=`TTC` (so `GAATTC`/EcoRI is a gap-0 IR
palindrome), revcomp(`TTACG`)=`CGTAA` (so `TTACGAAAAAACGTAA` is a gap-6 IR), revcomp(`GGATCC`)=
`GGATCC` (BamHI). Edge-case semantics all matched sources: single copy is not a tandem repeat,
non-primitive units collapse to the primitive period, gap-0 IR = palindrome, no library match →
Unclassified/Unknown with a documented STR-1–6 bp Simple_repeat fallback. Spans 0-based
end-exclusive. **Stage A findings: none.** The one documented simplification (exact substring
containment instead of Smith-Waterman/Repbase homology) is a Framework limitation flagged in the
evidence Assumption, not an error — and was tightened in Stage B.

## Stage B — implementation

Tandem detection scans unit lengths `1..min(len/minCopies, 60)`, skips non-primitive units,
counts adjacent equal blocks, enforces `minCopies`/`minRepeatLength`, and reports left-maximal
arrays de-duplicated by span (legitimate overlapping arrays such as `TTCGA`×2 are genuine, not
false positives; tests filter by exact span). Inverted detection computes
`GetReverseComplementString` of each left arm and searches the right arm within gap 0..50. Input
contract verified: null → `ArgumentNullException`, `minCopies<2` / `minRepeatLength<1` → throws,
empty → empty. Case-insensitivity verified (`gaattc` → `GAATTC` IR).

### Defect found and fixed — `ClassifyRepeat` bidirectional containment

Matching used `query.Contains(element) || element.Contains(query)`. The **reverse direction**
(`element.Contains(query)`) let a trivially short query inherit a class merely because a longer
consensus contained those letters — e.g. `ClassifyRepeat("A", {Alu, Satellite})` → `"SINE/Alu"`
and `ClassifyRepeat("GGCCGGG", {Alu})` → `"SINE/Alu"` (both **WRONG**; a 1-bp query matches almost
any library). This contradicts RepeatMasker, which **screens the query for occurrences of known
library elements** (element found *within* the query) and reports the best match. **Fix:** matching
made one-directional (`element ⊆ query`); XML doc updated to the screen-query-for-elements model.
Post-fix `"A"` and `"GGCCGGG"` → `"Unknown"`, while the real Alu-containing M5 query still →
`"SINE/Alu"`. Evidence Assumption 1 and TestSpec Assumption 1 corrected to drop the "either
direction / is contained by" language.

**Why it was uncaught:** no MUST test exercised the reverse direction or a short/sub-consensus
query, so the defect was invisible and would pass a deliberately-wrong impl (code-echo blind spot).

## Findings

- **One real defect (fixed):** `ClassifyRepeat` bidirectional containment misclassified short
  queries; corrected across code + tests + Evidence + TestSpec. End state CLEAN.
- **Test-quality gate PASS (HARD gate):** before, 13 tests, none reaching the reverse-containment
  path; after, added M7b (`"A"`→Unknown), M7c (`GGCCGGG` Alu-fragment→Unknown), M5b (longest-match
  tiebreak), M6b, and an INV-3 reverse-complement-arm test. All expectations trace to external
  sources; exact-value asserts, no weakened assertions or widened tolerances.
- **Scope note:** `ClassifyRepeat` uses exact-substring containment, not Smith-Waterman/Repbase
  homology — a documented Framework/Simplified limitation, not a defect. See
  [[repetitive-element-detection]] for the full deviation writeup.
