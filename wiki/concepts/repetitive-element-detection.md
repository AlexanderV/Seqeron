---
type: concept
title: "Repetitive element detection and classification"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/ANNOT-REPEAT-001-Evidence.md
source_commit: 3ecff734f0954cc6cd463f353d3994e99068f93f
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-repeat-001-evidence
      evidence: "Test Unit ID: ANNOT-REPEAT-001 ... Algorithm: Repetitive Element Detection and Classification (tandem repeats, inverted repeats, repeat-class assignment)"
      confidence: high
      status: current
---

# Repetitive element detection and classification

Finding and typing the repeated substrings that make up a large fraction of genomes. Seqeron's
repeat analyzer covers three distinct sub-problems, each with its own source-backed definition,
validated together as [[annot-repeat-001-evidence|ANNOT-REPEAT-001]]. This page is the shared
anchor for the whole **repeats / tandem family** — sibling units (GENOMIC-REPEAT,
GENOMIC-TANDEM, microsatellite/STR, low-complexity, etc.) should link here rather than
re-deriving the same definitions. See [[test-unit-registry]] for how units are tracked and
[[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The three sub-problems

### 1. Tandem repeats (head-to-tail)

A **tandem repeat** is a motif whose copies are directly adjacent with no intervening sequence —
a *head-to-tail* arrangement (e.g. `ATTCG ATTCG ATTCG`). The minimum is **two** adjacent copies;
a motif occurring once is not a tandem repeat. By repeat-unit length:

- **Microsatellite / STR** — 1–6 bp unit (Simple_repeat).
- **Minisatellite** — 10–60 bp unit.
- **Satellite** — larger units.

**Primitive-unit rule:** report the shortest period. `AAAAAA` is the mononucleotide `A`
repeated, not the dinucleotide `AA` or trinucleotide `AAA` — reporting a non-primitive unit
double-counts.

### 2. Inverted repeats (reverse-complement arms)

An **inverted repeat** (IR) is a left arm `W` followed downstream by its reverse complement.
Following IUPACpal (Hampson et al. 2021), an IR has the form **W W̄ᴿ** (perfect, ungapped), or
**W G W̄ᴿ** for a gapped IR with a spacer/loop `G`, `|G| ≥ 0`. An **imperfect** IR allows up to
`k` mismatches: Hamming distance `δ_H(W, W̄ᴿ) ≤ k`.

- **Zero gap ⇒ palindrome.** When `|G| = 0` the composite is an even-length reverse-complement
  palindrome (e.g. `GAATTC` → arm `GAA`, revcomp arm `TTC`).
- Detection parameters: minimum arm length, maximum arm length, maximum gap, maximum mismatches.

### 3. Repeat-class assignment

Assigning a query to a repeat class from the RepeatMasker/Repbase vocabulary:
**SINE, LINE, LTR, DNA** (DNA transposons), **Satellite, Simple_repeat, Low_complexity,
Small RNA, Unclassified/Unknown**. RepeatMasker itself classifies by *homology* — best
Smith-Waterman-Gotoh match above a score threshold against the Repbase library — and returns
Unknown when nothing matches.

## Deviation: classification is exact-substring, not scored alignment

Seqeron's `ClassifyRepeat(sequence, repeatDb)` **does not** run Smith-Waterman against a curated
Repbase library (out of scope for one unit). Instead it screens the query for library elements
**exactly contained** within it (element ⊆ query, one-directional) and assigns the class of the
**longest** such match; with no match it falls back to motif-size Simple_repeat classification,
else Unknown. Only the *matching relaxation* (exact substring vs. scored homology) is assumed —
the class vocabulary is source-backed. The one-directional containment prevents a trivially short
query from being forced into a class just because a longer consensus happens to contain its
letters. Documented as a Framework/Simplified [[research-grade-limitations|limitation]], not an
invented constant.

## Structural invariants (good test oracles)

- Every reported tandem repeat's `sequence` equals `input[start..end]` and is an integer number
  of unit copies.
- IR arms are exact reverse complements (within `k` mismatches for imperfect IRs).
- Reported spans are half-open `[start, end)` (the worked example `ATTCGATTCGATTCG` → unit
  `ATTCG`, 3 copies, span `[0, 15)`).
