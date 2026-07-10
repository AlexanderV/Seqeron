---
type: concept
title: "Known motif search (multi-pattern exact substring matching)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
  - docs/algorithms/Motif_Analysis/Known_Motif_Search.md
source_commit: 58f37bc5de666c59a60b8e7997c0894c16768c96
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: genomic-motifs-001-evidence
      evidence: "Test Unit ID: GENOMIC-MOTIFS-001 ... Algorithm: Known Motif Search (multi-pattern exact substring matching)"
      confidence: high
      status: current
---

# Known motif search (multi-pattern exact substring matching)

**Known motif search** locates a *set of already-known query motifs* in a subject
sequence by **exact** substring matching, reporting **all** start positions of each
motif. It is the "you know what you are looking for" end of motif analysis — distinct
from *discovery* (finding over-represented / unknown motifs) and from *degenerate*
matching (IUPAC / PROSITE / PWM). Its *inverse* is de novo
[[overrepresented-kmer-discovery|motif discovery]] — which searches for **unknown**
over-represented motifs rather than matching a supplied set. Seqeron exposes it as `GenomicAnalyzer.FindMotif`.
Validated under test unit **GENOMIC-MOTIFS-001**; the validation record is
[[genomic-motifs-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The exact-matching problem

Given a text `T` (length m) and a pattern `P` (length n), the exact-matching problem
is *"find all occurrences of P in T"* — the **set of all start positions** where `P`
aligns in `T` (Gusfield, via the Tufts COMP 150GEN notes). "Known motif search"
generalizes this to a **set** of query motifs: each motif independently maps to its
list of positions in the subject (Biopython `Seq.search` semantics — a list of query
terms yields per-motif hits).

## Overlaps are all reported (THE correctness rule)

The one rule a correct implementation must honour: **overlapping occurrences are all
counted**. `P=aaa`, `T=aaaaa` → **three** occurrences at positions **0, 1, 2**
(Gusfield). A greedy/non-overlapping scan that skips past each match under-reports and
is a defect for this problem. Biopython confirms the same distinction:
`Seq("AAAA").count_overlap("AA")` returns **3**, whereas non-overlapping
`count("AA")` returns 2.

## API contract and invariants

| Aspect | Behaviour |
|--------|-----------|
| Positions | **0-based** start indices, **sorted ascending** (deterministic public contract) |
| Multi-motif set | per-motif position list; a motif with **no** occurrence is **omitted** from the result (not an empty entry) |
| Case | motifs **upper-cased** before search; the result is keyed by the upper-cased motif |
| Empty / whitespace motif | contributes **no** entry — `FindMotif` returns `Array.Empty<int>()` (an empty pattern has no meaningful occurrence set) |
| Absent pattern | the empty set of positions |

## Worked oracles

- **Overlapping homopolymer** (Gusfield): T=`AAAAA`, P=`AAA` → **{0, 1, 2}** (count 3).
- **EcoRI motif** (real biology): T=`GAATTCAAAGAATTC`, P=`GAATTC` (the EcoRI
  palindrome, Wikipedia "Restriction site") → **{0, 9}** (count 2, non-overlapping by
  construction).
- **Multi-motif set** (Biopython): T=`ACGTACGTAA`, set `{ACGT, AA, TTT}` →
  `ACGT`→{0, 4}, `AA`→{8}, `TTT`→(absent, omitted).

## Scope

This is the **exact** DNA-motif finder. The degenerate/consensus family lives in
separate units (IUPAC-degenerate matching, PROSITE patterns, position-weight-matrix
scanning) that relax the exact-equality test; known-motif search is the exact-equality
baseline. A worked biological instance of the position-weight-matrix branch is
[[splice-acceptor-site-prediction]] (the AcceptorPwm 3' splice-site scorer). The
**canonical-biology** specialization is [[regulatory-element-detection]]
(MOTIF-REGULATORY-001) — the same all-occurrences scan applied to a *fixed, cited catalog*
of regulatory consensus strings rather than a caller-supplied set. It shares the *exact all-occurrences* semantics with the exact-match engine
behind [[longest-common-substring]] and [[dot-plot-word-match]], but operates on a
caller-supplied set of literal query motifs rather than deriving shared substrings.
The **single-pattern** K-mer-family counterpart is [[k-mer-positions]]
(`KmerAnalyzer.FindKmerPositions`) — the same 0-based ascending all-overlapping-occurrences
semantics for one k-mer against one sequence (a single position list rather than a per-motif map).

## Reference sources

**Tufts COMP 150GEN** (Gusfield-derived exact-matching definition + overlapping
`aaa`/`aaaaa` rule), **Biopython `Bio.Seq`** (`search` multi-motif semantics,
`count_overlap` overlap confirmation, upper-casing), and **Wikipedia "Restriction
site"** (EcoRI `GAATTC` biological motif). **No deviations**; the only assumptions are
two API-shape policies — empty-motif → no entry, and upper-case result keys — neither
of which any source defines as an algorithm-correctness parameter.
