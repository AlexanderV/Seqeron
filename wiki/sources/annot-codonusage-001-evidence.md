---
type: source
title: "Evidence: ANNOT-CODONUSAGE-001 (Relative Synonymous Codon Usage — RSCU)"
tags: [validation, annotation]
doc_path: docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
sources:
  - docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
source_commit: 1a7037743423e6db365bbc5460f2d2d04f9384a5
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ANNOT-CODONUSAGE-001

The validation-evidence artifact for test unit **ANNOT-CODONUSAGE-001** (Relative Synonymous
Codon Usage — RSCU). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is
summarized in [[relative-synonymous-codon-usage]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13) — the LIRMM (Rivals et al., Université de
  Montpellier) RSCU methods page for the verbatim `RSCU_{i,j} = (n_i · x_{i,j}) / Σ x_{i,j}`
  formula, the [0, n_i] range, and the no-bias = 1.0 property; PMC2528880 (begomovirus codon
  usage) for the observed-over-expected definition and the >1 preferred / <1 under-represented
  reading; **Sharp & Li (1986)**, *Nucleic Acids Research* 14(19):7737–7749, as the original
  primary citation; the **CodonU** reference implementation (`internal_comp.py::rscu`) for the
  exact aggregation and sense-codon handling; and NCBI genetic-code table 1 for the
  synonymous-family map.
- **Algorithm spec** — for codon `j` of amino acid `i`, `RSCU = n_i · x_{i,j} / Σ x_{i,j}`
  where `n_i` is the synonymous-family size and the denominator sums occurrences over all
  synonymous codons of `i`. Uniform usage ⇒ 1.0; the family's RSCU values sum to `n_i`.
- **Datasets** — Leucine-only worked example `CTTCTTCTGTTA` ⇒ RSCU 3.0 / 1.5 / 1.5 / 0 / 0 / 0
  (Σ = 6 = n_i); uniform Phe `TTTTTC` ⇒ 1.0 / 1.0; single-codon Met `ATGATG` ⇒ 1.0.
- **Corner cases** — counts pooled across all reference sequences before computing (RSCU is on
  the aggregate, not per-sequence); stop codons excluded (sense codons only, Biopython
  `forward_table`); single-codon amino acids (Met, Trp; n_i = 1) are always exactly 1.0; a
  fully-unobserved family (Σ = 0) is undefined in the base definition.
- **Recommended coverage** — MUST tests for the Leu values, uniform = 1.0, single-codon = 1.0,
  multi-sequence pooling, stop-codon exclusion, and the Σ-over-family = n_i invariant; SHOULD
  tests for case-insensitivity and partial trailing-codon skip; COULD test for null/empty.

## Assumptions (from the artifact)

1. **Genetic code defaults to Standard (NCBI table 1).** An overload accepts a `GeneticCode`
   for non-standard tables; the default is an API convention, not a correctness gap.
2. **Zero-count family ⇒ RSCU 0 for every codon in it.** The base definition (LIRMM, raw
   counts) leaves RSCU undefined when Σ = 0; Seqeron reports 0.0 per member to avoid division
   by zero. The **CAI 0.5 pseudocount** (Sharp & Li 1987, for `log(0)`) is deliberately **not**
   applied — it is a CAI-specific convention, distinct from plain RSCU, and would belong to the
   sibling CAI unit, not here.

No contradictions between sources; the LIRMM formula, PMC2528880 definition, and CodonU code
are algebraically identical. The only deviations are the two API-default assumptions above,
neither of which affects a real (observed-family) RSCU value.
