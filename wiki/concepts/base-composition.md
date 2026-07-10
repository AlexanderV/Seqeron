---
type: concept
title: "Sequence composition (base/residue counts, fractions, GC content)"
tags: [sequence-statistics, composition]
sources:
  - docs/Evidence/SEQ-STATS-001-Evidence.md
  - docs/Evidence/SEQ-COMPOSITION-001-Evidence.md
  - docs/Evidence/SEQ-DINUC-001-Evidence.md
  - docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
  - docs/Evidence/SEQ-MW-001-Evidence.md
source_commit: 6e8fde12868aa0db4347950e4cf52449588e0b68
created: 2026-07-10
updated: 2026-07-10
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
methods, and `SummarizeNucleotideSequence` is the umbrella's thin aggregation wrapper re-exposing
this GC content / count tally. The skew members of that umbrella live on the sibling
[[nucleotide-composition-skew]].

## The three outputs

1. **Nucleotide composition** — an exact partition of **Length** into per-base counts
   `A, T, G, C, U` plus `CountN` (the `N` "any base" code) and `CountOther` (everything else:
   degenerate IUPAC codes, gaps, `X`, …). The seven counts sum to Length.
2. **GC content** = `(G + C) / (A + T + G + C + U)`, a fraction in **[0, 1]** (Biopython
   `Bio.SeqUtils.gc_fraction`). **Empty sequence ⇒ 0** (the zero-length denominator is handled,
   not an exception).
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
