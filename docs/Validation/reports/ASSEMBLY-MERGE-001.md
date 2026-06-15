# Validation Report: ASSEMBLY-MERGE-001 — Contig Merging (suffix–prefix overlap collapse)

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.MergeContigs(string contig1, string contig2, int overlapLength)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm

Both Langmead JHU lecture-note PDFs are image/encoded and could not be read by WebFetch's text
model, so they were downloaded and text was extracted locally with `pdftotext -layout` and read
directly this session.

1. **Langmead, "Assembly & shortest common superstring"** (assembly_scs.pdf), extracted locally:
   - Overlap definition (verbatim): *"Overlap: length-l suffix of X matches length-l prefix of Y,
     where l is given."* (and *"an 'overlap' is when a suffix of X of length ≥ l exactly matches a
     prefix of Y"*). Confirms the suffix/prefix model used by the unit.
   - `suffixPrefixMatch(x, y, k)` pseudocode (verbatim): docstring *"Return length of longest suffix
     of x of length at least k that matches a prefix of y. Return 0 if there no suffix/prefix match
     has length at least k."* and guard *`if len(x) < k or len(y) < k: return 0`*. Confirms (a) the
     longest match is the reported overlap and (b) the overlap is bounded by `min(|x|,|y|)`.
   - No-overlap = concatenation (verbatim): *"Without requirement of 'shortest,' it's easy: just
     concatenate them."* Confirms INV-02 / overlap-0 behaviour.
   - Greedy-SCS trace (l = 1), first round (verbatim line): `2  BAAB  ABA  ABB  AAA  BBB  BBA  BAB`
     with *"Number in first column = length of overlap merged before that round."* → `BAA` + `AAB`
     at overlap **2** collapses to `BAAB` (M1). Confirmed.
   - Small SCS example (verbatim): `S: AAA AAB ABB BBB BBA` → `SCS(S): AAABBBA`, with the alignment
     `AAA / AAB / ABB / BBB / BBA` each offset by one (overlap 2), length 7 (M2). Confirmed.
2. **Langmead, "Overlap Layout Consensus assembly"** (assembly_olc.pdf), extracted locally:
   - *"Fragments are contigs (short for contiguous)."*
   - Stages (verbatim): *"Layout — Bundle stretches of the overlap graph into contigs"*;
     *"Consensus — Pick most likely nucleotide sequence for each contig."*
   - *"Assume for given string pair we report only the longest suffix/prefix match."*
   These corroborate the contig/overlap collapse model and the "longest match" convention.

### Formula check

The cited merge primitive is `merge(X, Y, l) = X + Y[l:]` with `|merge| = |X| + |Y| − l`. This is the
standard superstring collapse (drop the length-`l` prefix of the right string, keep one copy of the
overlap) consistent with the Langmead overlap definition and `suffixPrefixMatch` (the returned length
`ln` is exactly the prefix of `y` removed). No log/normalisation/units involved; symbol-exact.

### Edge-case semantics check

- `l = 0` → `X + Y` (sourced: "just concatenate them"). Defined and sourced.
- `l ≤ 0` (negative) → no usable overlap → `X + Y`. Derived from the overlap definition (an overlap
  length is a positive suffix/prefix length); a non-positive length is "no overlap".
- `l > min(|X|,|Y|)` → `X + Y`. Sourced directly from the `suffixPrefixMatch` guard
  `if len(x) < k or len(y) < k: return 0` (an overlap cannot be longer than the shorter string).
- `l = min(|X|,|Y|)` (boundary) → full collapse of the shorter prefix. Valid overlap boundary.
- null operands → `ArgumentNullException` (library input-validation convention; not a biology claim).

### Independent cross-check (numbers)

Reference re-implementation of the sourced formula `merge(x,y,l)= x+y[l:] (else x+y)` in Python,
run this session, reproduces every spec/test value exactly:

| Inputs | l | Result | Len | Source basis |
|--------|---|--------|-----|--------------|
| BAA + AAB | 2 | `BAAB` | 4 | Langmead greedy trace (M1) — exact |
| chain AAA·AAB·ABB·BBB·BBA | 2 each | `AAABBBA` | 7 | Langmead SCS example (M2) — exact |
| BAA + AAB | 0 | `BAAAAB` | 6 | "just concatenate" (M3) |
| ACGT + CGTAA | 3 | `ACGTAA` | 6 | formula; verified overlap (suffix "CGT" = prefix "CGT") (M4) |
| GATTACA + ACATGAA | 3 | `GATTACATGAA` | 11 | formula; verified overlap (suffix "ACA" = prefix "ACA") (M5) |
| AC + GTAA | 3 (>min 2) | `ACGTAA` | 6 | guard → concat (S1) |
| BAA + AAB | −2 | `BAAAAB` | 6 | non-positive → concat (S2) |
| "" + AAB | 0 | `AAB` | 3 | identity (C1) |
| BAA + "" | 0 | `BAA` | 3 | identity (C2) |

M1/M2 trace to **exact strings printed in the primary source** (BAAB; AAABBBA). M3–C2 are direct
applications of the sourced formula; M4/M5 additionally use *genuine* verified suffix/prefix overlaps
(checked by hand this session), so the asserted outputs are biologically meaningful, not arbitrary.

### Findings / divergences

None. Description, formula, invariants (INV-01…INV-04), and every edge case are confirmed against
the two primary sources. Stage A: **PASS**.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:629–641`.

```csharp
ArgumentNullException.ThrowIfNull(contig1);
ArgumentNullException.ThrowIfNull(contig2);
const int NoOverlap = 0;
if (overlapLength <= NoOverlap || overlapLength > Math.Min(contig1.Length, contig2.Length))
    return contig1 + contig2;
return contig1 + contig2.Substring(overlapLength);
```

### Formula realised correctly?

Yes. The merge path is `contig1 + contig2.Substring(overlapLength)` = `X + Y[l:]`, the exact sourced
primitive, giving length `|X| + |Y| − l`. The fallback covers both `l ≤ 0` (no overlap → concat) and
`l > min(|c1|,|c2|)` (exceeds shorter string → concat), matching the `suffixPrefixMatch` guard. Null
operands throw before any indexing. O(n) single substring + concatenation, as documented.

### Cross-verification table recomputed vs code

The full NUnit suite was run (`dotnet test --no-build`); all 12 `MergeContigs` tests pass, and their
asserted values are identical to the independent reference table above. Every value matches the
primary source.

### Variant/delegate consistency

`MergeContigs` is a single static primitive with no `*Fast`/instance/overload variants; nothing to
reconcile. It is intentionally decoupled from overlap *discovery* (`FindOverlap`), consistent with
the documented contract.

### Test quality audit (HARD gate)

- **Sourced, not code echoes:** M1 (`BAAB`) and M2 (`AAABBBA`) assert the *exact strings printed in
  the Langmead primary source* — these would fail against a wrong implementation. M3–C2 assert the
  exact sourced-formula outputs. A deliberately-wrong merge (e.g. trimming `contig1` instead, or
  off-by-one) is caught by the exact-string MUST assertions.
- **No green-washing:** every MUST/SHOULD uses `Is.EqualTo(exact)`; no Greater/AtLeast/Contains/range
  on a known value; S3/S4 assert the exact exception type; no test skipped/ignored/weakened.
- **Coverage of all logic:** all three code paths exercised — null-throw (S3, S4), concat-fallback via
  `l ≤ 0` (M3, S2) and via `l > min` (S1), and substring-merge (M1, M2, M4, M5, C3). Boundary
  `l = shared-prefix` (M4), length invariant (M5), identity with empty operands (C1, C2), and the
  prefix/suffix containment invariant INV-04 (C3) are all covered.
- **Honest green:** the FULL unfiltered suite passes — **Failed: 0, Passed: 6529** — and the changed
  files (none; report/ledger only) build warning-free. (Pre-existing NUnit2007 warnings live in
  unrelated test files and were not introduced here.)

Note: C3 (`StartsWith`/`EndsWith`) is a property test that alone would be weak, but it is supplementary
to the exact-value M1–M5 tests; it is not the load-bearing assertion for the merged output. Acceptable.

### Findings / defects

None.

## Verdict & follow-ups

- **Stage A: PASS.** Description, formula, invariants, and edge cases independently confirmed against
  the Langmead SCS and OLC primary sources (extracted locally this session).
- **Stage B: PASS.** Code realises `X + Y[l:]` exactly with the sourced concat fallback; tests assert
  exact source-traced values and cover every branch and documented edge case.
- **Test-quality gate: PASS** (sourced values, no green-washing, full-branch coverage, honest green
  6529/0).
- **End-state: CLEAN.** No defect found; algorithm fully functional. No code changed.
