---
type: concept
title: "Nucleotide composition skew (AT skew / GC skew)"
tags: [sequence-statistics, composition, chromosome]
mcp_tools:
  - at_skew
  - nucleotide_composition
sources:
  - docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md
  - docs/Evidence/SEQ-STATS-001-Evidence.md
  - docs/Evidence/SEQ-ATSKEW-001-Evidence.md
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
  - docs/Evidence/SEQ-REPLICATION-001-Evidence.md
source_commit: 4beb586f662e59a317383820445c2dc5b176a08e
created: 2026-07-10
updated: 2026-07-13
---

# Nucleotide composition skew (AT skew / GC skew)

**Strand compositional skew** measures how far a single DNA strand departs from
intra-strand equifrequency of complementary bases. The family has two members with
identical shape:

- **AT skew** = `(A − T) / (A + T)`
- **GC skew** = `(G − C) / (G + C)`

Both were introduced by **Lobry (1996)** — the founding observation was a *"departure
from intrastrand equifrequency between A and T or between C and G, showing that the
substitution patterns of the two strands of DNA were asymmetric."* Both skew members were
first delivered *together* under the original **SEQ-STATS-001** sequence-statistics umbrella
([[seq-stats-001-evidence]]), alongside the [[base-composition]] tally. The **SEQ-ATSKEW-001**
test unit ([[seq-atskew-001-evidence]]) later validates the AT-skew member as its own registry
entry; GC skew is the sibling member. Both skews are also computed together (as
`OverallGcSkew`/`OverallAtSkew`, plus a **windowed** GC-skew profile and its population variance)
by the composite **SEQ-GC-ANALYSIS-001** GC-analysis unit — see
[[windowed-gc-profile-and-variance]] and [[seq-gc-analysis-001-evidence]]; a standalone `gc-skew`
unit is separately flagged, not yet ingested. [[test-unit-registry]] tracks the units and
[[algorithm-validation-evidence]] describes the artifact pattern.

## Definition and range

- **Bounded** to `[−1, +1]`. For AT skew: `+1 ⇔ T = 0` (all A among A/T), `−1 ⇔ A = 0`
  (all T among A/T); for GC skew symmetrically with G/C.
- **Zero-denominator convention:** when the relevant pair is absent (`A + T = 0`, resp.
  `G + C = 0`) the skew is defined as **`0.0`**, not NaN or an exception. This is the
  Biopython `Bio.SeqUtils.GC_skew` behaviour (`ZeroDivisionError → 0.0`) that the library
  follows.
- **Case-insensitive counting**; **non-canonical symbols ignored** — gaps, `N`, and any
  IUPAC ambiguity code contribute to neither numerator nor denominator ("does NOT look at
  any ambiguous nucleotides"). The library normalizes via `ToUpperInvariant` and counts
  only the two relevant bases.

Worked values (arithmetic consequences of the formula, no library run needed):
`AAAA → +1.0`, `TTTT → −1.0`, `ATAT → 0.0`, `AAAT → +0.5`, `GGCC → 0.0` (no A/T),
`AAATGGGCCC → +0.5` (G/C ignored), `aaat → +0.5` (case-insensitive).

## Why it matters — replication-strand asymmetry

Skew is not random: on a bacterial chromosome the sign of GC (and AT) skew tends to be
constant along a replichore and **flips at the replication origin and terminus**, because
the leading and lagging strands accumulate different mutational/substitution biases. A
**cumulative skew** plotted along the sequence therefore locates the origin/terminus as
its extrema — the practical use that motivated Lobry's work. That locator is a distinct
algorithm in its own right, synthesized on [[replication-origin-cumulative-skew]]
(SEQ-REPLICATION-001): it integrates an *integer* running skew (G:+1, C:−1, A/T:0) and reads
the origin off the global **minimum** and the terminus off the **maximum** (Grigoriev 1998).
AT skew is the weaker,
sometimes atypical signal: **Charneski et al. (2011)** showed Firmicute AT skew arises
from *selection*, not mutation, so the two skews need not co-vary.

Skew is the *asymmetry* view of the same per-base tally that [[base-composition]] counts as
magnitudes/fractions (A/T/G/C/U counts, GC content) — the two are siblings, both surfaced by
the SEQ-COMPOSITION-001 doc. It is also the compositional-asymmetry counterpart to two other
composition statistics in the wiki: the CpG observed/expected density of
[[cpg-island-detection]] (a dinucleotide composition ratio) and the GC-variability heuristic
used inside [[centromere-analysis]] (which flags GC-skew as a chromosome-analysis unit
warranting its own concept — this page).

## Scope / assumptions

The formula itself is fully sourced (Charneski 2011; Lobry 1996; corroborated by the
Wikipedia "GC skew" entry). Only the **symbol-handling convention** for the AT-skew member
(case-insensitive counting, ignore-everything-not-A/T) is taken by analogy from the shipped
`GC_skew` reference implementation, because Biopython ships `GC_skew` but not an AT-skew
line — a documented assumption in [[seq-atskew-001-evidence]], matching the repository
implementation.

## Implementation

Both scalar skews live in `GcSkewCalculator` (`Seqeron.Genomics.Analysis`; primary spec
`docs/algorithms/Extended_GC_Skew_Analysis/AT_Skew.md`, test unit SEQ-ATSKEW-001). AT skew
has two overloads over a shared private core: `CalculateAtSkew(string)` is the canonical entry
(upper-cases, counts `A`/`T` via `string.Count`, returns the skew — `null`/empty ⇒ `0`), and
`CalculateAtSkew(DnaSequence)` forwards the already-normalized value object (`null` ⇒
`ArgumentNullException`). Cost is **O(n) time / O(1) space** — two linear symbol counts. The
repository suffix tree does **not** apply: AT skew is a two-symbol count, not a substring-search
or occurrence-enumeration problem.

This unit returns only the **single global scalar** — windowed/cumulative *AT*-skew profiles
and AT-skew-based origin location are deliberately *not* implemented. For localization along a
sequence the same class exposes the GC-skew-based `CalculateWindowedGcSkew` /
`CalculateCumulativeGcSkew` / `PredictReplicationOrigin` (see
[[windowed-gc-profile-and-variance]] and [[replication-origin-cumulative-skew]]).
