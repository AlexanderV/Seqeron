# Validation Report: SEQ-GC-001 — GC Content Calculation

- **Validated:** 2026-06-12   **Area:** Composition
- **Canonical method(s):** `SequenceExtensions.CalculateGcContent(ReadOnlySpan<char>)` (percentage),
  `CalculateGcFraction(ReadOnlySpan<char>)` (fraction); delegates `*Fast(string)`, `DnaSequence.GcContent()`.
- **Stage A verdict:** 🟡 PASS-WITH-NOTES
- **Stage B verdict:** ✅ PASS

## Stage A — Description

**Sources opened.** Wikipedia *GC-content* gives the formula
(G+C)/(A+T+G+C) × 100%. Biopython `gc_fraction` (default `"remove"` mode) counts only valid
nucleotides in the denominator. Both are correctly cited in `TestSpecs/SEQ-GC-001.md`.

**Formula check.** Spec formula matches the source exactly: numerator = G+C, denominator =
A+T+G+C (+U for RNA), ×100 for percentage, no ×100 for fraction. INV-3 (content = fraction×100)
is a true identity.

**Edge-case semantics.** All edge cases have a *defined, sourced* expected value: empty → 0
(Biopython "returns zero for an empty sequence"), no valid nucleotides → 0, ambiguous chars
excluded from both numerator and denominator ("remove" mode), U treated as a valid non-GC
nucleotide. No "implementation-defined" gaps.

**Independent cross-check.** Hand-computed against the Biopython `gc_fraction(seq,"remove")`
values: `""`→0, `ACTG`→0.50, `ACTGN`→0.50, `CCTGNN`→0.75, `GDVV`→1.00, `GCAU`→0.50,
`GGAUCUUCGGAUCU`→0.50, `NNNNN`→0. All consistent with the formula.

**Finding (note, not a defect).** Biopython additionally counts the ambiguity code **S** (G|C)
as GC and **W** (A|T) as AT in its denominator; Seqeron excludes S and W entirely (treats them
like any other non-ACGTU character). The spec documents this divergence (§1.3). For standard
A/C/G/T/U sequences the results are identical; they only differ on sequences containing literal
S/W ambiguity codes. → recorded as PASS-WITH-NOTES, not a failure.

## Stage B — Implementation

**Code path.** [SequenceExtensions.cs:28-58](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs#L28-L58)
(`CalculateGcContent`) and [:81-111](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs#L81-L111)
(`CalculateGcFraction`).

**Formula realised correctly.** The loop counts G/C/g/c into both `gcCount` and `validCount`,
and A/T/U (both cases) into `validCount` only; everything else is skipped. Returns
`gcCount/validCount*100` (or `/validCount` for fraction), with `validCount==0 → 0` and
`IsEmpty → 0`. This is exactly the validated formula and the "remove"-mode denominator. The
S/W exclusion noted in Stage A is visible here (no `S`/`W` cases in the switch).

**Cross-verification recomputed against the code** — traced each value through the switch:

| Input | Source (Biopython) | Code result | Match |
|-------|--------------------|-------------|-------|
| `""` | 0 | 0 (IsEmpty) | ✓ |
| `ACTG` | 0.50 | 2/4 = 0.50 | ✓ |
| `ACTGN` | 0.50 | 2/4 (N skipped) | ✓ |
| `CCTGNN` | 0.75 | 3/4 (NN skipped) | ✓ |
| `GDVV` | 1.00 | 1/1 (D,V skipped) | ✓ |
| `GCAU` | 0.50 | 2/4 | ✓ |
| `GGAUCUUCGGAUCU` | 0.50 | 7/14 | ✓ |
| `NNNNN` | 0 | validCount 0 → 0 | ✓ |

**Variant/delegate consistency.** `CalculateGcContentFast`/`CalculateGcFractionFast` and
`DnaSequence.GcContent()` all forward to the canonical Span methods — no independent logic to
diverge.

**Numerical robustness.** Integer counts, double division; div-by-zero guarded by the
`validCount==0` check. No overflow concern for realistic lengths.

**Test quality.** `SequenceExtensions_CalculateGcContent_Tests.cs` asserts the exact Biopython
cross-check values (CCTGNN=75, GDVV=100, GCAU=50, NNNNN=0, RNA=50) with `Is.EqualTo`, plus
property-based invariants in `Properties/GcContentProperties.cs`. Assertions are exact-value,
not "no-throw" tautologies, and cover the Stage-A edge cases.

## Verdict & follow-ups

- **Stage A: 🟡 PASS-WITH-NOTES** — description correct; documented S/W ambiguity-code divergence
  from Biopython (immaterial for standard sequences).
- **Stage B: ✅ PASS** — implementation faithfully realises the validated formula; all
  cross-check values reproduce; delegates consistent; tests assert exact sourced values.
- **No defects logged.** Optional future enhancement (not required): if exact Biopython parity
  on S/W is ever desired, count S as GC and W as AT — but this would diverge from the strict
  Wikipedia ACGT(U) formula, so leave as-is unless a consumer needs it.
