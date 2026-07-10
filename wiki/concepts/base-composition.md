---
type: concept
title: "Sequence composition (base/residue counts, fractions, GC content)"
tags: [sequence-statistics, composition]
sources:
  - docs/Evidence/SEQ-COMPOSITION-001-Evidence.md
source_commit: 2fa9affeb77d7240ffffd91ffd809647c4297484
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
---

# Sequence composition (base/residue counts, fractions, GC content)

**Sequence composition** is the foundational tally of *what a sequence is made of* — the
counts and fractions of each symbol. It underlies almost every downstream sequence statistic.
The **SEQ-COMPOSITION-001** unit ([[seq-composition-001-evidence]]) validates three related
outputs over the standard {A,T,G,C,U} + amino-acid alphabets; [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the artifact pattern.

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
- **Dinucleotide composition** — [[cpg-island-detection]] uses a CpG observed/expected ratio,
  the dinucleotide-frequency generalization of single-base composition.
- **Windowed composition entropy** — [[windowed-sequence-complexity-profile]] computes a
  Shannon entropy of *base composition* per window; composition is its per-window input.
- **GC variability** — [[centromere-analysis]] uses a GC-content heuristic over windows.

GC content is also a design/QC constraint elsewhere in the library (e.g. the 30–80% probe GC
window in [[taqman-probe-design-rules]] and the 40–60% balanced-GC codon-optimization
strategy in [[codon-optimization]]) — all reading the same underlying composition.
