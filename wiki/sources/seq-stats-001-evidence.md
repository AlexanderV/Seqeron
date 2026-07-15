---
type: source
title: "Evidence: SEQ-STATS-001 (sequence composition statistics — the ORIGINAL composition umbrella: nucleotide composition, GC content, GC/AT skew)"
tags: [validation, sequence-statistics, composition]
doc_path: docs/Evidence/SEQ-STATS-001-Evidence.md
sources:
  - docs/Evidence/SEQ-STATS-001-Evidence.md
source_commit: 6e8fde12868aa0db4347950e4cf52449588e0b68
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-STATS-001

The validation-evidence artifact for test unit **SEQ-STATS-001** — **sequence composition
statistics**: nucleotide composition (exact A/T/G/C/U/N/Other partition of Length), **GC
content** `(G+C)/(A+T+G+C+U)`, and the two strand-asymmetry skews **GC skew** `(G−C)/(G+C)`
and **AT skew** `(A−T)/(A+T)`. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]] tracks the
unit.

## The original composition umbrella — connective tissue

**SEQ-STATS-001 is the ORIGINAL home of the sequence-composition methods.** Several later
composition units are documentation/registry consolidations *over* the methods first delivered
here, not fresh implementations. Most explicitly, [[seq-composition-001-evidence|SEQ-COMPOSITION-001]]
records in its Change History that it is a **duplicate/consolidated Registry entry for the two
composition methods already delivered under SEQ-STATS-001** (TestSpec §7). So this artifact is
the connective tissue tying together the piecemeal-ingested composition family:

- **Nucleotide composition + GC content** — synthesized on [[base-composition]] (counts,
  fractions, the A/T/G/C/U/`CountN`/`CountOther` partition, `(G+C)/(A+T+G+C+U)` in `[0,1]`).
- **GC skew + AT skew** — synthesized on the sibling family concept
  [[nucleotide-composition-skew]] (both `(X−Y)/(X+Y)` members, bounded `[−1,+1]`,
  zero-denominator ⇒ `0.0`). SEQ-STATS-001 is where GC skew **and** AT skew were first
  delivered *together*, before AT skew was split into its own [[seq-atskew-001-evidence|SEQ-ATSKEW-001]]
  registry unit and before both were re-bundled by the composite
  [[windowed-gc-profile-and-variance|SEQ-GC-ANALYSIS-001]] (`OverallGcSkew`/`OverallAtSkew` + a
  windowed profile).
- **Aggregation & protein variant** — the coverage recommendations also exercise
  `SummarizeNucleotideSequence` (the summary wrapper re-exposing the same GC content / counts —
  validated in its own right as **SEQ-SUMMARY-001**, [[seq-summary-001-evidence]], where it is
  shown to also bundle the sequence's Shannon entropy, linguistic complexity, and Tm) and
  `CalculateAminoAcidComposition` (exact per-residue counts + Length over the 20
  IUPAC amino-acid letters — the amino-acid composition also covered by [[base-composition]];
  MW/pI/hydrophobicity are explicitly *out of scope* here, belonging to SEQ-MW / SEQ-PI /
  SEQ-HYDRO units — see [[molecular-weight]], [[isoelectric-point]],
  [[hydrophobicity-gravy-and-profile]]).

No new concept is warranted: every method this umbrella exercises is already synthesized on an
existing concept page. This page records only what the artifact adds.

## What this file records

- **Online sources:**
  - **Biopython `Bio.SeqUtils`** (`gc_fraction`, `GC_skew`; rank 3, reference impl) —
    `gc_fraction` core `gc = sum(seq.count(x) for x in "CGScgs")`, `length = gc + Σ ATWUatwu`,
    returns `gc/length`, a float in **[0, 1]**; default `ambiguous="remove"` counts only GCS
    and includes only ACTGSWU in the length; **empty sequence ⇒ 0**; **case-insensitive**
    (`"CGScgs"` includes lowercase). `GC_skew` = `(g−c)/(g+c)` (`g=count("G")+count("g")`,
    `c=count("C")+count("c")`, default `window=100`), catching `ZeroDivisionError → 0.0`.
  - **Wikipedia "GC skew"** (rank 4, tracing to Lobry 1996) — `GC skew = (G−C)/(G+C)`,
    `AT skew = (A−T)/(A+T)`; positive GC skew ⇒ G over C, negative ⇒ C over G; Lobry's original
    notation `(C−G)/(C+G)` is flipped by modern implementations to `(G−C)/(G+C)`. (Skew detail
    lives on [[nucleotide-composition-skew]].)
  - **Lobry (1996)** *Mol Biol Evol* 13(5):660–665 (DOI 10.1093/oxfordjournals.molbev.a025626;
    rank 1, primary) — the founding *"departure from intrastrand equifrequency between A and T
    or between C and G"* that makes G/C and A/T skew meaningful.
- **Datasets (hand-derived worked examples, arithmetic — no library run needed):**

  | Input | A/T/G/C/U | GC content | GC skew | AT skew |
  |-------|-----------|-----------|---------|---------|
  | `ATGC` | 1/1/1/1/0 | 2/4 = **0.5** | 0/2 = 0 | 0/2 = 0 |
  | `GGGC` | 0/0/3/1/0 | 4/4 = **1.0** | 2/4 = 0.5 | **0** (a+t=0) |
  | `AAAT` | 3/1/0/0/0 | 0/4 = **0.0** | **0** (g+c=0) | 2/4 = 0.5 |
  | `GCCC` | 0/0/1/3/0 | 4/4 = **1.0** | −2/4 = −0.5 | **0** (a+t=0) |
  | `AAUUGGCC` | 2/0/2/2/2 | 4/8 = **0.5** | 0/4 = 0 | 2/2 = **1.0** |

- **Corner cases / failure modes:** empty/null ⇒ **all-zero composition** (GC content 0);
  **no G or C** ⇒ GC skew **0.0**, **no A or T** ⇒ AT skew **0.0** (both from Biopython's
  caught `ZeroDivisionError`); **mixed case** ⇒ result is case-insensitive (GC counting
  includes lowercase).

## Deviations and assumptions

**One documented assumption — degenerate IUPAC codes (S, W, R, Y, …) are NOT counted toward
composition totals.** Biopython's `gc_fraction` counts `S` toward GC and `W` toward the length
denominator; the repository counts **only A/T/G/C/U** toward GC/AT totals and routes other
letters to `CountN`/`CountOther`. Over the standard {A,T,G,C,U} alphabet (this unit's scope) the
two agree **exactly**; the difference manifests only on degenerate symbols. Documented as an
intentional simplification, not an invented constant. **No source contradictions** — Biopython,
Wikipedia, and Lobry agree on the formulas and alphabet.

Recommended coverage (from the artifact): MUST — GC content `(G+C)/(A+T+G+C+U)`; GC skew
`(G−C)/(G+C)` incl. a negative case; AT skew `(A−T)/(A+T)`; exact A/T/G/C/U/N/Other counts +
Length; empty/null ⇒ all-zero. SHOULD — case-insensitivity; zero-denominator skews ⇒ 0. COULD —
`SummarizeNucleotideSequence` aggregates the same GC content/counts (plus entropy/complexity/Tm —
[[seq-summary-001-evidence]]); `CalculateAminoAcidComposition` returns exact residue counts + length.
