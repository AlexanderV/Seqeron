---
type: source
title: "Evidence: ONCO-SIG-001 (SBS-96 trinucleotide context catalog — pyrimidine-strand folding)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SIG-001-Evidence.md
sources:
  - docs/Evidence/ONCO-SIG-001-Evidence.md
source_commit: 6fdbd84d46e8ac221dadd222b315412645d44051
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SIG-001

The validation-evidence artifact for test unit **ONCO-SIG-001** — the **SBS-96 single-base-substitution
trinucleotide context catalog**: classification of somatic single-base substitutions into the 96 canonical
COSMIC channels by folding purine-reference mutations onto the pyrimidine strand. The **twenty-ninth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its own
concept, [[sbs96-mutational-signature-catalog]]; [[test-unit-registry]] tracks the unit.

Note: this unit is **only the 96-channel catalog / classification** step. The downstream NMF/NNLS
signature-**exposure** fitting against COSMIC reference signatures is a separate concern, not covered here.

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **COSMIC — SBS96 Mutational Signatures** (Wellcome Sanger Institute; rank 5 database / rank 2 de-facto
    standard catalogue) — six substitution subtypes **C>A, C>G, C>T, T>A, T>C, T>G**; the 96 mutation types =
    "6 types of substitution × 4 types of 5′ base × 4 types of 3′ base"; each substitution is "referred to by
    the **pyrimidine** of the mutated Watson–Crick base pair"; context = the bases "immediately 5′ and 3′ to
    each mutated base".
  - **SigProfilerMatrixGenerator — Bergstrom et al. (2019), BMC Genomics 20:685** (PMC6717374; rank 1
    peer-reviewed + rank 3 reference implementation that produces COSMIC matrices) — SBS-6 = the six classes;
    SBS-96 elaborates each with one 5′ and one 3′ adjacent base (4×4 = 16 trinucleotides per class → 6×16 = 96);
    the mutated base is **centred** in the trinucleotide (e.g. ACA>AAA); **reverse-complement folding**
    (verbatim): "using the purine base of the Watson–Crick base-pair for classifying mutation types will
    require taking the reverse complement sequence of each of the classes of SBS-96".
  - **Alexandrov et al. (2013), Nature 500:415-421** (rank 1 primary; body behind an auth gate, verbatim
    definitions supplied by the search-index summary + the same group's COSMIC/SigProfiler pages) — the 96
    types = 6 substitution × 4 × 4; the six subtypes; the study analysed 4,938,362 mutations from 7,042 cancers.
  - **Complementarity (molecular biology) — Wikipedia** (rank 4) — the complement map A↔T, C↔G used for the
    reverse-complement fold.

- **Documented corner cases / failure modes:** purine reference base (A/G) is **not** one of the six canonical
  pyrimidine substitutions and MUST be reverse-complemented before counting; non-SBS variants (indels,
  doublet/multi-base DBS) belong to other catalogues (ID, DBS), not the 96-channel SBS spectrum; a non-ACGT
  flanking base (e.g. N) has no defined context and cannot be classified; `ref == alt` is not a mutation
  (out of scope).

- **Datasets (deterministic worked oracles):** seven `5′,ref,alt,3′` rows → expected channel, each computed
  independently of any implementation. The three pyrimidine-reference rows classify to themselves
  (A,C,A,A → **A[C>A]A**; T,C,T,G → **T[C>T]G**; G,T,C,A → **G[T>C]A**); the four purine-reference rows fold
  (T,G,T,A → **T[C>A]A**; C,A,G,T → **A[T>C]G**; G,G,C,C → **G[C>G]C**; A,A,T,A → **T[T>A]T**). Worked fold of
  row 4: plus-strand TGA, G→T; ref G is purine → reverse-complement context TGA→TCA and substitution G>T→C>A →
  **T[C>A]A**. Catalog-count invariant: exactly **96 channels (6×4×4)**; Σ channel counts = number of
  classifiable SBS variants; every key is one of the 96 canonical pyrimidine labels.

- **Coverage recommendations:** MUST test the six pyrimidine substitutions classify to themselves
  (centre base unchanged); MUST test purine-reference mutations fold by reverse-complement (the seven rows);
  MUST test the catalog partition (Σ counts = classifiable variants, keys ⊆ 96 labels); MUST test the enumerated
  channel set is exactly the 96 canonical labels; SHOULD test null/empty, non-ACGT context, `ref==alt`,
  non-single-base variant; COULD test case-insensitive (lower-case) input.

## Deviations and assumptions

- **ASSUMPTION — channel label format `5'[REF>ALT]3'`.** COSMIC/SigProfiler render the trinucleotide with the
  mutated base centred; the exact textual form (bracket `A[C>A]A` vs underline `ACA>AAA`) is a display choice
  that does **not** change which variants fall in which class — the partition into the 96 pyrimidine-keyed
  classes is identical either way. The bracket form is adopted for the keys.

No source contradictions — COSMIC, SigProfilerMatrixGenerator, and Alexandrov 2013 give the identical
6×4×4 = 96 definition and pyrimidine convention; the Wikipedia complement map supplies the folding chemistry.
