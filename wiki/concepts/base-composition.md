---
type: concept
title: "Sequence composition (base/residue counts, fractions, GC content)"
tags: [sequence-statistics, composition]
mcp_tools:
  - analyze_gc_content
sources:
  - docs/algorithms/Sequence_Composition/Sequence_Composition.md
  - docs/algorithms/Sequence_Composition/Sequence_Composition_Statistics.md
  - docs/algorithms/Statistics/Sequence_Summary.md
  - docs/Evidence/SEQ-STATS-001-Evidence.md
  - docs/Evidence/SEQ-COMPOSITION-001-Evidence.md
  - docs/Evidence/SEQ-DINUC-001-Evidence.md
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
  - docs/Evidence/SEQ-MW-001-Evidence.md
  - docs/Evidence/SEQ-SUMMARY-001-Evidence.md
source_commit: f3f84d5b4393c99d3e3cb25ee1cae338be09e293
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-composition-001-evidence
      evidence: "Test Unit ID: SEQ-COMPOSITION-001 ... Algorithm: Sequence Composition (nucleotide composition + amino-acid composition)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-stats-001-evidence
      evidence: "Test Unit ID: SEQ-STATS-001 ... Algorithm: Sequence Composition Statistics (nucleotide composition, GC content, GC/AT skew) — the original umbrella that first delivered the nucleotide-composition + GC-content methods."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-summary-001-evidence
      evidence: "Test Unit ID: SEQ-SUMMARY-001 ... Algorithm: Sequence Summary — SummarizeNucleotideSequence aggregates the nucleotide composition + GC content tally (plus Shannon entropy, linguistic complexity, Tm) into one per-sequence record."
      confidence: high
      status: current
---

# Sequence composition (base/residue counts, fractions, GC content)

**Sequence composition** is the foundational tally of *what a sequence is made of* — the
counts and fractions of each symbol. It underlies almost every downstream sequence statistic.
The **SEQ-COMPOSITION-001** unit ([[seq-composition-001-evidence]]) validates three related
outputs over the standard {A,T,G,C,U} + amino-acid alphabets; [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the artifact pattern.

**Original home.** The nucleotide-composition and GC-content methods were first delivered under
the **SEQ-STATS-001** sequence-statistics *umbrella* ([[seq-stats-001-evidence]]);
SEQ-COMPOSITION-001 is its later duplicate/consolidated registry entry over those same two
methods, and `SummarizeNucleotideSequence` (validated separately as **SEQ-SUMMARY-001** —
[[seq-summary-001-evidence]]) is the umbrella's aggregation wrapper: a per-sequence summary record
that re-exposes this GC content / count tally **and bundles** the sequence's Shannon entropy,
linguistic complexity, and melting temperature alongside it (a pure field-by-field aggregation of
the already-validated per-metric methods). The skew members of that umbrella live on the sibling
[[nucleotide-composition-skew]].

## The three outputs

1. **Nucleotide composition** — an exact partition of **Length** into per-base counts
   `A, T, G, C, U` plus `CountN` (the `N` "any base" code) and `CountOther` (everything else:
   degenerate IUPAC codes, gaps, `X`, …). The seven counts sum to Length.
2. **GC content** = `(G + C) / (A + T + G + C + U)`, a fraction in **[0, 1]** (Biopython
   `Bio.SeqUtils.gc_fraction`). **Empty sequence ⇒ 0** (the zero-length denominator is handled,
   not an exception). Its complement **AT content** = `(A + T + U) / (A + T + G + C + U)` — note
   uracil is folded into the AT-fraction *numerator* (RNA `U` pairs with `A` like `T`), so for a
   pure-RNA input AtContent counts U. This is a deliberate DNA/RNA asymmetry from the SEQ-STATS-001
   primary spec (`Sequence_Composition_Statistics.md`): the **AT skew** on the sibling
   [[nucleotide-composition-skew]] uses the DNA-specific `(A − T)/(A + T)` — U is *not* in the skew,
   matching the Lobry/Wikipedia definition, whereas the AT *content* fraction does include it.
3. **Amino-acid composition** — the analogous exact per-residue count over the **20** standard
   IUPAC single-letter codes (A C D E F G H I K L M N P Q R S T V W Y) plus Length.

Worked values (arithmetic consequences of the definitions): `ATGC` → GC content **0.5**;
`GGGC` → **1.0**; `AAUUGGCC` → A/T/G/C/U = 2/0/2/2/2, GC content **0.5**; `MKVLWA` → six
residues each count 1, Length 6.

## Counting conventions

- **Case-insensitive** — Biopython's GC counting explicitly includes lowercase (`"CGScgs"`),
  so composition ignores case; the library normalizes before counting.
- **Canonical vs non-canonical** — only A/T/G/C/U are canonical bases; `N` and degenerate
  codes are **tracked separately** (`CountN`/`CountOther`), never folded into the four/five
  canonical bases.
- **Degenerate codes excluded from GC/AT totals** — a documented **assumption**: Biopython's
  `gc_fraction` counts `S` toward GC and `W` toward the denominator, whereas this repository
  counts only A/T/G/C/U. Over the {A,T,G,C,U} alphabet (this unit's scope) the two agree
  exactly; the divergence is confined to degenerate symbols. See [[seq-composition-001-evidence]].

## Relationship to neighbouring composition statistics

Composition is the base layer that several other wiki concepts build on:

- **Strand skew** — [[nucleotide-composition-skew]] derives `(G−C)/(G+C)` and `(A−T)/(A+T)`
  from these same base counts (with the zero-denominator ⇒ `0.0` convention). Skew is the
  *asymmetry* view of the same tally; this page is the *magnitude/fraction* view. The
  SEQ-COMPOSITION-001 doc mentions both GC/AT skew alongside GC content, which is why the two
  concepts are siblings.
- **Dinucleotide composition** — [[dinucleotide-relative-abundance]] counts adjacent base pairs and
  scores each pair's Karlin odds ratio `ρ_XY = f_XY/(f_X·f_Y)` against these single-base frequencies,
  the dinucleotide-frequency generalization of single-base composition; [[cpg-island-detection]] is
  its `CG`-specialized CpG observed/expected ratio.
- **Windowed composition entropy** — [[windowed-sequence-complexity-profile]] computes a
  Shannon entropy of *base composition* per window; composition is its per-window input.
- **Molecular weight** — [[molecular-weight]] (SEQ-MW-001) is the **mass-weighted** view of this
  same per-monomer tally: `Σ (monomer mass) − (len−1)·water` over the identical {A,T,G,C,U}+
  amino-acid alphabets, sharing the case-fold and skip-unknown contract. Composition counts the
  monomers; MW sums their Daltons.
- **Windowed GC profile & variability** — [[windowed-gc-profile-and-variance]] slides a window
  along the sequence, emitting a per-window GC%/GC-skew profile and the population variance of
  each series (the composite `GcAnalysisResult`, SEQ-GC-ANALYSIS-001, which re-exposes GC
  content as a **percentage ×100** rather than this page's `[0,1]` fraction). [[centromere-analysis]]
  uses a related GC-content heuristic over windows.

GC content is also a design/QC constraint elsewhere in the library (e.g. the 30–80% probe GC
window in [[taqman-probe-design-rules]] and the 40–60% balanced-GC codon-optimization
strategy in [[codon-optimization]]) — all reading the same underlying composition.

## Implementation surface (SEQ-COMPOSITION-001 primary spec)

The canonical primary spec (`docs/algorithms/Sequence_Composition/Sequence_Composition.md`)
adds method-level detail to the Evidence view. Two entry points live in
`Seqeron.Genomics.Analysis/SequenceStatistics.cs`:

- `SequenceStatistics.CalculateNucleotideComposition(string)` → the `NucleotideComposition`
  record struct (`Length`, per-base `A/T/G/C/U` counts, `CountN`, `CountOther`, `GcContent`,
  `AtContent`, `GcSkew`, `AtSkew`). A single upper-casing linear pass classifies each character
  into exactly one bucket, then closed-form ratios with zero-denominator guards. **O(n)** time,
  **O(1)** space (fixed-size counters).
- `SequenceStatistics.CalculateAminoAcidComposition(string)` → the `AminoAcidComposition`
  record struct (`Length`, per-residue `Counts` over the 20 standard codes, plus derived
  physico-chemical fields owned by the SEQ-MW/PI/HYDRO units). **O(n)** time, **O(σ)** space
  (σ = distinct residues ≤ 26).

**Consolidation, not re-implementation.** The spec's own note records that SEQ-COMPOSITION-001 is
a *duplicate Registry entry* over the two methods first delivered under **SEQ-STATS-001**
([[seq-stats-001-evidence]]); it is resolved by consolidation, and the full GC-content/skew
treatment lives on that SEQ-STATS-001 umbrella (its skew members on [[nucleotide-composition-skew]]). No
search/matching is involved (a pure counter scan), so the repository suffix tree does not apply.

**Invariants** (spec §2.4): `0 ≤ GcContent, AtContent ≤ 1` (INV-01, subset-count over
canonical-base total); `−1 ≤ GcSkew, AtSkew ≤ 1` (INV-02, difference over sum of non-negative
counts); `CountA+T+G+C+U+N+Other = Length` (INV-03, each char in exactly one bucket); amino-acid
`Length = Σ Counts` (INV-04). Zero-denominator skew ⇒ `0.0`, matching Biopython's
`ZeroDivisionError` handling.

**Not implemented.** Weighted ambiguous-base GC (Biopython `gc_fraction(ambiguous="weighted")`)
is out of scope; degenerate-heavy sequences should use Biopython `gc_fraction` directly. Over the
canonical {A,T,G,C,U} alphabet the two agree exactly.

## The `SummarizeNucleotideSequence` aggregator (SEQ-SUMMARY-001 primary spec)

`SummarizeNucleotideSequence` is the top-level per-sequence **summary record** that bundles the
headline descriptors of a DNA/RNA sequence into one object so a caller gets them in a single call.
Its primary spec (`docs/algorithms/Statistics/Sequence_Summary.md`, unit SEQ-SUMMARY-001) is
reconciled here — it is a **pure field-by-field aggregation** that introduces no new computation:
every field reproduces, bit-for-bit, the value its canonical per-metric method returns on the same
input. Its correctness contract is exactly that field-wise equality, so the summary owns no new
concept — each aggregated method already has its home (this page, plus [[shannon-entropy]],
[[linguistic-complexity]], [[melting-temperature]]). The Evidence view is [[seq-summary-001-evidence]].

**Entry point** — `SequenceStatistics.SummarizeNucleotideSequence(string sequence)` in
`Seqeron.Genomics.Analysis/SequenceStatistics.cs`, returning the `readonly record struct`
`SequenceStatistics.SequenceSummary`. **Field → sub-analyzer** map (each row is an invariant,
INV-01…INV-06):

| Summary field | Type | Delegated to | Source of value |
|---------------|------|--------------|-----------------|
| `Length` | `int` | `CalculateNucleotideComposition(S).Length` | raw character count |
| `GcContent` | `double` | `CalculateNucleotideComposition(S).GcContent` | `(G+C)/(A+T+G+C+U)` ∈ [0,1] (this page) |
| `Entropy` | `double` | `CalculateShannonEntropy(S)` | Shannon `H = −Σ p·log₂ p` (bits); the **general-alphabet** `SequenceStatistics` kernel of [[shannon-entropy]] |
| `Complexity` | `double` | `CalculateLinguisticComplexity(S)` | linguistic vocabulary-usage in (0,1) — the **mean**-based [[linguistic-complexity]], *not* the strict Trifonov product |
| `MeltingTemperature` | `double` | `CalculateMeltingTemperature(S, useWallaceRule: S.Length < 14)` | Wallace `2(A+T)+4(G+C)` for len<14 else GC/Marmur-Doty `64.9 + 41·(GC−16.4)/N` ([[melting-temperature]]) |
| `Composition` | `IReadOnlyDictionary<char,int>` | counts from `CalculateNucleotideComposition(S)` | a **6-entry** dict keyed A,T,G,C,U,N |

**The Tm length-dispatch is the summary's only branching decision:** it passes
`useWallaceRule: sequence.Length < 14` — the strict `<` boundary means length exactly **14** takes
the GC/Marmur-Doty branch. The 14-nt threshold is the sibling SEQ-TM-001 convention
(`ThermoConstants.WallaceMaxLength`); the summary is tested only for **equality** with
`CalculateMeltingTemperature`, so the threshold's correctness belongs to that unit, not the
aggregation.

**Cost** — **O(n)** time / **O(n)** space, dominated by the linguistic-complexity term (k-mer sets
up to k=6); the other metrics are O(n) single passes. No string searching/matching is performed, so
the repository suffix tree does not apply (a counting-and-arithmetic aggregation).

**Edge cases** — null/empty input returns the **degenerate zero summary** (Length 0, GcContent 0,
Entropy 0, Complexity 0, MeltingTemperature 0, all composition counts 0) with no exception, matching
each per-metric method's empty handling; input is **case-insensitive** (each delegate uppercases
internally); RNA `U` and `N` are counted (distinct from T). Worked oracle: `"ATGCATGC"` → Length 8,
GcContent 0.5, Entropy 2.0 (= log₂4), MeltingTemperature 24.0 (len<14 → Wallace `2·4 + 4·4`).
