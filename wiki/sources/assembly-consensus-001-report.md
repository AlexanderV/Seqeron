---
type: source
title: "Validation report: ASSEMBLY-CONSENSUS-001 (consensus computation — column-wise plurality/threshold)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-CONSENSUS-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-CONSENSUS-001.md
source_commit: ac0a26ca9923867f41c95e2d4a7046a9712ccca6
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-CONSENSUS-001

The two-stage **validation write-up** for test unit **ASSEMBLY-CONSENSUS-001** (Consensus
Computation — the **C** of Overlap-Layout-[[overlap-layout-consensus-assembly|Consensus]]),
validated 2026-06-15 in a fresh context. This is the *report* artifact that feeds one row of
the [[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The algorithm itself is summarized in the concept
[[consensus-sequence]]; the two-stage methodology is the [[validation-protocol]]. Distinct
from [[assembly-consensus-001-evidence]] (the pre-implementation evidence artifact) — this is
the independent re-validation verdict.

Canonical method under test:
`SequenceAssembler.ComputeConsensus(IReadOnlyList<string> alignedReads, double threshold = 0.5, char ambiguous = 'N')`.

## Verdict

**Stage A: 🟡 PASS-WITH-NOTES · Stage B: ✅ PASS · State: ✅ CLEAN.** No defect found; no code
or test change required. The PASS-WITH-NOTES is *not* a defect — it flags two documented,
parameter-reachable **default-value divergences** from Biopython (below). Full unfiltered suite
**6532 passed, 0 failed**; `dotnet build` 0 errors (the 4 build warnings live in an unrelated
`ApproximateMatcher_EditDistance_Tests.cs`, untouched here).

## Stage A — description (algorithm faithfulness)

Theory checked against primary sources opened this session, not the repo's own assertions:

- **Biopython `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus`** (raw source fetched at tag
  `biopython-179`) — confirmed **verbatim**: signature `dumb_consensus(self, threshold=0.7,
  ambiguous="X", require_multiple=False)`; consensus spans the full alignment length
  (`con_len = get_alignment_length()`); per-column tally **skips gaps** (`!= "-" and != "."`),
  `num_atoms` per non-gap residue; bounds guard `if n < len(record.seq)`; max-set rule
  (`> max_size` resets `max_atoms`, `== max_size` appends); decision rule
  `(len(max_atoms)==1) and (max_size/num_atoms >= threshold)` → residue, else `ambiguous`.
- **EMBOSS `cons`** — plurality cut-off "below which there is no consensus" corroborates the
  threshold→ambiguous semantics (cited, not re-fetched; Biopython already pins the rule exactly).
- **Wikipedia "Consensus sequence"** — the definition ("most frequent residues … at each
  position in a sequence alignment") and IUPAC `N` = any base (default-symbol rationale).

**Formula check.** Commit predicate `unique-max ∧ num_atoms>0 ∧ (max_size/num_atoms ≥ threshold)`,
else ambiguous — matches Biopython with the all-gap short-circuit; inclusive `≥`, gap set
`{'-','.'}`, and length = longest read (`con_len`) all match the source. Edge cases each have a
sourced expected behaviour: empty list → `""`; null → `ArgumentNullException`; all-gap column
(`num_atoms==0`) → ambiguous with **no div-by-zero**; tie → ambiguous; sub-threshold → ambiguous;
ragged reads span the longest read. INV-01..INV-05 are genuine properties.

**Independent cross-check (numbers).** The validator ran the **Biopython 1.85** reference
`dumb_consensus` on the test datasets (gap-padded to equal length so `-` is skipped). All ten
outputs equal the TestSpec expected values — every expected value is externally sourced, not a
code echo:

| Case | Reads | thr | Biopython output |
|------|-------|-----|------------------|
| M1 | ACGT,ACGT,ACGT | 0.5 | `ACGT` |
| M2 | ACGT,ACGT,ACGT,TCGT | 0.7 | `ACGT` |
| M3 | AC,AC,TC | 0.7 | `NC` |
| M4 | A,G | 0.5 | `N` |
| M5 | A-GT,A.GT,ACGT | 0.5 | `ACGT` |
| M6 | ACGT,ACG | 0.5 | `ACGT` |
| M7 | A-T,A-T | 0.5 | `ANT` |
| M9a | A,A,A,A,T | 0.7 | `A` |
| M9b | A,A,T | 0.7 | `N` |
| C1 | A,G (ambiguous=X) | 0.5 | `X` |

Biopython's own deprecation-notice example (`ACGT/ATGT/ATGT` → `ANGT` at default 0.7; col1
T=2/3≈0.667 < 0.7 → N) independently confirms the rule.

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:861-926`:
null guard `ArgumentNullException.ThrowIfNull` (L866), empty → `""` (L868); `length` = max read
length = `con_len` (L871-875); per-column reset of `counts`/`numAtoms` (L881-882) + bounds
`pos >= read.Length → continue` (L886, = Biopython `n < len(record.seq)`); gap skip `'-'||'.'`
after `ToUpperInvariant` (L888-889); tie via `maxCount` (L896-912, = `len(max_atoms)`); commit
predicate `maxCount==1 && numAtoms>0 && (double)maxSize/numAtoms >= threshold` (L918-920) else
`ambiguous` (L922), with `numAtoms>0` short-circuiting before the division. Every clause maps
1:1 to the Biopython rule. Hand-recomputation of each table case matches both code and
Biopython 1.85.

**Variant/delegate consistency.** Only one canonical method. The MCP wrapper
`AlignmentTools.ComputeConsensus` (`src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs:348`)
delegates to it with defaults plus an equal-length precondition; no divergent consensus logic.

**Test-quality audit (HARD gate: PASS).** All 12 unit tests assert exact `Is.EqualTo`
full-string values re-derived from the Biopython rule (no Greater/AtLeast/Contains, no widened
tolerances, no skipped tests). A deliberately-wrong impl (`MaxBy` tie-break, or `>` instead of
`≥`) would fail M3/M4/M9. Coverage exercises every decision branch and both defaults: unanimous
(M1), threshold-above (M2/M9a), threshold-below (M3/M9b), tie (M4), gap-skip for both `-` and
`.` (M5), ragged (M6), all-gap column (M7), empty (M8), null (S1), lowercase (S2), custom
ambiguous symbol (C1).

## Findings

- **None.** No defects, no code or test change. State ✅ CLEAN.
- **Two documented default divergences (carried, not defects)** — the reason Stage A is
  PASS-WITH-NOTES rather than PASS; both are presentation-only and fully parameter-reachable, so
  they never change the decision rule (detailed on the concept [[consensus-sequence]]):
  - `threshold` default **0.5** (simple-majority / EMBOSS plurality) vs Biopython's documented
    **0.7** — pass `threshold: 0.7` to reproduce Biopython exactly (verified by M2/M3/M9).
  - `ambiguous` default **`'N'`** (DNA IUPAC any-base) vs Biopython's **`'X'`** (protein) — pass
    `ambiguous: 'X'` for parity (verified by C1).
- The Biopython `require_multiple` flag is intentionally **not implemented** (documented; default
  `False` is the rule under test) — out of scope, not a defect.
- **No follow-ups.**
