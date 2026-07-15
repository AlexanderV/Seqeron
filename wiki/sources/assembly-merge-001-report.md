---
type: source
title: "Validation report: ASSEMBLY-MERGE-001 (contig merging — suffix–prefix overlap collapse)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-MERGE-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-MERGE-001.md
source_commit: 6abf4edca8f18ac8c0d17c25f3949d7c1dea135d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-MERGE-001

The two-stage **validation write-up** for test unit **ASSEMBLY-MERGE-001** (contig merging — the
suffix–prefix overlap collapse that stitches two overlapping contigs into one longer superstring),
validated 2026-06-15 in a fresh context. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm description and
the shipped code. The algorithm itself is summarized in the concept
[[contig-merge-overlap-collapse]] (anchor of the assembly MERGE family), and the wider campaign is
[[validation-and-testing]]. Distinct from [[assembly-merge-001-evidence]] (the pre-implementation
evidence artifact, sourced from `docs/Evidence/`) — this is the independent re-validation verdict.

Canonical method under test:
`SequenceAssembler.MergeContigs(string contig1, string contig2, int overlapLength)` — a single static
collapse primitive (no `*Fast`/instance/overload variants), intentionally decoupled from overlap
*discovery* (`FindOverlap`).

## Verdict

**Stage A: ✅ PASS · Stage B: ✅ PASS · State: ✅ CLEAN.** No defect found; the algorithm is fully
functional and no code was changed. All 12 `MergeContigs` NUnit tests pass; the full unfiltered suite
is **6529 passed, Failed: 0**, and the (report/ledger-only) change set builds warning-free (pre-existing
NUnit2007 warnings live in unrelated test files).

## Stage A — description (algorithm faithfulness)

Theory checked against two primary sources opened this session, independent of repo artifacts. Both
Langmead JHU lecture-note PDFs are image/encoded (WebFetch could not read them), so they were
downloaded and text-extracted locally with `pdftotext -layout`:

- **Langmead (JHU) — "Assembly & shortest common superstring"** (assembly_scs.pdf) — the overlap
  definition verbatim (*"a length-l suffix of X matches a length-l prefix of Y, where l is given"*),
  the `suffixPrefixMatch(x, y, k)` primitive (longest suffix ≥ k matching a prefix, else 0; guard
  `if len(x) < k or len(y) < k: return 0`, so overlap ≤ `min(|X|,|Y|)`), superstring-by-collapse
  (*"without requirement of 'shortest,' it's easy: just concatenate them"*), and the greedy traces
  (`BAA`+`AAB` at overlap 2 → `BAAB`; the `{AAA,AAB,ABB,BBB,BBA}` chain → SCS `AAABBBA`, length 7).
- **Langmead (JHU) — "Overlap Layout Consensus assembly"** (assembly_olc.pdf) — contig terminology
  (*"fragments are contigs"*), the OLC stages (*"Layout — bundle stretches of the overlap graph into
  contigs"*; *"Consensus — pick most likely nucleotide sequence for each contig"*), and *"report only
  the longest suffix/prefix match"*.

**Formula check.** The merge primitive is `merge(X, Y, l) = X + Y[l:]` with `|merge| = |X| + |Y| − l`
— the standard superstring collapse (drop the length-`l` prefix of the right string, keep one copy of
the overlap), symbol-exact against the Langmead overlap definition and `suffixPrefixMatch` (the
returned length `ln` is exactly the prefix of `y` removed). No log/normalisation/units involved.

**Edge-case semantics (sourced).** `l = 0` → `X + Y` (sourced: *"just concatenate them"*);
`l ≤ 0` (negative) → `X + Y` (a non-positive length is "no overlap"); `l > min(|X|,|Y|)` → `X + Y`
(sourced from the `suffixPrefixMatch` guard — an overlap cannot exceed the shorter string);
`l = min(|X|,|Y|)` boundary → full collapse of the shorter prefix; null operands →
`ArgumentNullException` (library input-validation convention, not a biology claim).

**Independent cross-check (numbers).** A Python re-implementation of the sourced formula (separate
from repo code) reproduced every spec/test value exactly:

| Inputs | l | Result | Len | Source basis |
|--------|---|--------|-----|--------------|
| `BAA` + `AAB` | 2 | `BAAB` | 4 | Langmead greedy trace (M1) — exact printed string |
| chain `AAA·AAB·ABB·BBB·BBA` | 2 each | `AAABBBA` | 7 | Langmead SCS example (M2) — exact printed string |
| `BAA` + `AAB` | 0 | `BAAAAB` | 6 | "just concatenate" (M3) |
| `ACGT` + `CGTAA` | 3 | `ACGTAA` | 6 | verified overlap (suffix `CGT` = prefix `CGT`) (M4) |
| `GATTACA` + `ACATGAA` | 3 | `GATTACATGAA` | 11 | verified overlap (suffix `ACA` = prefix `ACA`) (M5) |
| `AC` + `GTAA` | 3 (>min 2) | `ACGTAA` | 6 | guard → concat (S1) |
| `BAA` + `AAB` | −2 | `BAAAAB` | 6 | non-positive → concat (S2) |
| `""` + `AAB` | 0 | `AAB` | 3 | identity (C1) |
| `BAA` + `""` | 0 | `BAA` | 3 | identity (C2) |

M1/M2 trace to the **exact strings printed in the primary source**; M4/M5 use *genuine* verified
suffix/prefix overlaps, so the asserted outputs are biologically meaningful, not arbitrary. Stage A
= **PASS**, no divergences.

## Stage B — implementation (code review)

Code path `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:629–641`:

```csharp
ArgumentNullException.ThrowIfNull(contig1);
ArgumentNullException.ThrowIfNull(contig2);
const int NoOverlap = 0;
if (overlapLength <= NoOverlap || overlapLength > Math.Min(contig1.Length, contig2.Length))
    return contig1 + contig2;
return contig1 + contig2.Substring(overlapLength);
```

**Formula realised correctly.** The merge path `contig1 + contig2.Substring(overlapLength)` is
exactly `X + Y[l:]`, giving length `|X| + |Y| − l`. The single fallback covers both `l ≤ 0` (no
overlap → concat) and `l > min(|c1|,|c2|)` (exceeds shorter string → concat), matching the
`suffixPrefixMatch` guard; null operands throw before any indexing. O(n) single substring +
concatenation, as documented.

**Cross-verification.** The full NUnit suite was run (`dotnet test --no-build`); all 12 `MergeContigs`
tests pass and their asserted values are identical to the independent reference table above.

**Test-quality audit (HARD gate: PASS).** M1 (`BAAB`) and M2 (`AAABBBA`) assert the *exact strings
printed in the Langmead source* — a deliberately-wrong merge (trimming `contig1` instead, or
off-by-one) fails them. Every MUST/SHOULD uses `Is.EqualTo(exact)` — no `Greater`/`AtLeast`/`Contains`/
range on a known value; S3/S4 assert the exact exception type; nothing skipped or weakened. All three
code paths are exercised: null-throw (S3, S4), concat-fallback via `l ≤ 0` (M3, S2) and via `l > min`
(S1), and substring-merge (M1, M2, M4, M5, C3), plus the `l = min` boundary (M4), length invariant
(M5), and empty-operand identity (C1, C2). The full unfiltered suite passes **6529/0**. (C3's
`StartsWith`/`EndsWith` property test is supplementary to the exact-value M1–M5 assertions, not the
load-bearing one — noted as acceptable.)

## Findings

- **No algorithm defect. No code changed. End-state ✅ CLEAN.** Description, formula, invariants
  (INV-01…INV-04), and every edge case are independently confirmed against the Langmead SCS and OLC
  primary sources; tests assert exact source-traced values and cover every branch.
- **No contradictions** among the sources — the SCS notes, the OLC notes, and (per the evidence
  artifact) MIT 7.91J Lecture 6 all state the identical suffix-of-X / prefix-of-Y overlap definition.
- **No follow-ups.**
