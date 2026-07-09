---
type: concept
title: "SBS-96 trinucleotide context catalog — pyrimidine-strand folding"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-SIG-001-Evidence.md
source_commit: 6fdbd84d46e8ac221dadd222b315412645d44051
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-sig-001-evidence
      evidence: "Test Unit ID: ONCO-SIG-001 ... Algorithm: SBS-96 Single-Base-Substitution Trinucleotide Context Catalog (pyrimidine-strand folding)"
      confidence: high
      status: current
---

# SBS-96 trinucleotide context catalog

The Oncology family's **mutational-signature catalog** unit (**ONCO-SIG-001**): the first step of
single-base-substitution (SBS) mutational-signature analysis — classifying each somatic single-base
substitution into one of the **96 canonical COSMIC channels** and building the per-channel count vector
(the SBS-96 "spectrum") that downstream signature-exposure fitting decomposes. The literature-traced
record is [[onco-sig-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

This unit is **only the catalog / classification**. Fitting the spectrum against COSMIC reference
signatures (NMF extraction / NNLS exposure estimation) is a separate downstream concern.

## The 96 channels: 6 × 4 × 4

A channel is a substitution together with its immediate 5′ and 3′ flanking bases (the **trinucleotide
context**, with the mutated base **centred**):

- **6 substitution subtypes**, each named by the **pyrimidine** of the mutated Watson–Crick base pair:
  `C>A, C>G, C>T, T>A, T>C, T>G` (the SBS-6 classification).
- **× 4** choices of 5′ base **× 4** choices of 3′ base = 16 trinucleotides per subtype.
- **6 × 16 = 96** channels, labelled `5'[REF>ALT]3'` (e.g. `A[C>A]A`, `T[C>T]G`).

## Pyrimidine-strand folding (the key rule)

Because a substitution on one strand is equivalent to its reverse complement on the other, every SBS is
reported on the **pyrimidine strand**. When the mutated (reference) base is a **purine (A or G)**, it is
not one of the six canonical subtypes and MUST be folded:

1. Reverse-complement the **trinucleotide context** (complement map `A↔T, C↔G`, then reverse).
2. Complement the **substitution** itself (ref and alt).

The centre base then becomes a pyrimidine (C or T) and the mutation lands in one of the six subtypes.
Pyrimidine-reference mutations (C or T) are already canonical and classify **to themselves with the centre
base unchanged**.

**Worked fold:** plus-strand context `5'-TGA-3'`, mutation `G→T`. Reference `G` is a purine → fold. Complement
each base `TGA → ACT`, reverse → `TCA`. Centre `G→C` (pyrimidine), alt `T→A` → substitution `C>A`. Result:
5′ = T, sub = C>A, 3′ = A → **T[C>A]A**.

## Worked oracles

| 5′ | Ref | Alt | 3′ | Fold | Channel |
|----|-----|-----|----|------|---------|
| A | C | A | A | none (pyrimidine) | A[C>A]A |
| T | C | T | G | none (pyrimidine) | T[C>T]G |
| G | T | C | A | none (pyrimidine) | G[T>C]A |
| T | G | T | A | revcomp TGA→TCA, G>T→C>A | T[C>A]A |
| C | A | G | T | revcomp CAT→ATG, A>G→T>C | A[T>C]G |
| G | G | C | C | revcomp GGC→GCC, G>C→C>G | G[C>G]C |
| A | A | T | A | revcomp AAA→TTT, A>T→T>A | T[T>A]T |

## Catalog invariants

Building the 96-channel catalog from a multiset of variants is a **partition**:

- **Exactly 96 channels** (6 × 4 × 4); the enumerated key set is exactly the 96 canonical pyrimidine labels.
- **Σ per-channel counts = number of classifiable SBS variants**.
- Every input classifies to **one** of the 96 canonical labels.

## Corner cases and assumptions

- **Purine reference → fold** (never counted directly). A/G-reference mutations are reverse-complemented
  before counting.
- **Non-SBS variants excluded.** Only single-base substitutions are SBS-96 events; indels and doublet/
  multi-base substitutions belong to the separate ID and DBS catalogues.
- **Non-ACGT context** (e.g. `N` flank) has no defined trinucleotide context → not classifiable.
- **`ref == alt`** is not a mutation → out of scope.
- **ASSUMPTION — label rendering is cosmetic.** Bracket form `A[C>A]A` vs underline `ACA>AAA` is a display
  choice; the partition into the 96 pyrimidine-keyed classes is identical either way.

## Relation to the oncology family

The catalog is the **substrate** for mutational-signature analysis — a somatic-mutation-spectrum summary
independent of the allele-specific copy-number layer. It sits alongside the other tumor mutational-burden /
mutation-pattern biomarkers and is orthogonal to the copy-number-scar
[[homologous-recombination-deficiency-score]] and the immunotherapy biomarker
[[microsatellite-instability-detection]]; a high-`C>T`-at-CpG or SBS-specific spectrum is the kind of signal
those interpretation layers ([[cancer-variant-tier-classification-amp-asco-cap]],
[[clinical-actionability-oncokb-levels]]) can act on once exposures are fit.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the SBS-96 classification and catalog
construction. The 6×4×4 = 96 definition and pyrimidine convention are COSMIC SBS96, SigProfilerMatrixGenerator
(Bergstrom 2019), and Alexandrov 2013; the reverse-complement fold is SigProfiler's verbatim purine rule with
the Watson–Crick complement map. **Not for clinical or diagnostic use.** No source contradictions.
