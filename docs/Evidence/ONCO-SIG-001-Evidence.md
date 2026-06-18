# Evidence Artifact: ONCO-SIG-001

**Test Unit ID:** ONCO-SIG-001
**Algorithm:** SBS-96 Single-Base-Substitution Trinucleotide Context Catalog (pyrimidine-strand folding)
**Date Collected:** 2026-06-14

---

## Online Sources

### COSMIC — SBS96 Mutational Signatures (Wellcome Sanger Institute)

**URL:** https://cancer.sanger.ac.uk/signatures/sbs/sbs96/
**Retrieved:** WebFetch of the URL on 2026-06-14 (search query that surfaced it:
"COSMIC SBS mutational signatures 96 channels trinucleotide pyrimidine reverse complement classification").
**Accessed:** 2026-06-14
**Authority rank:** 5 (well-maintained bioinformatics database) / 2 (the de-facto standard catalogue)

**Key Extracted Points:**

1. **Six substitution subtypes (verbatim):** "C>A, C>G, C>T, T>A, T>C, and T>G." The page states
   the 96 mutation types are formed by "6 types of substitution x 4 types of 5' base x 4 types of 3' base".
2. **Pyrimidine convention (verbatim):** "Each of the substitutions is referred to by the pyrimidine of the
   mutated Watson—Crick base pair." Considering the pyrimidines of the Watson-Crick base pairs, there are
   only six possible substitutions (the six above).
3. **Context:** mutation types incorporate the bases "immediately 5' and 3' to each mutated base", generating
   96 possible mutation types.

### SigProfilerMatrixGenerator — Bergstrom et al. (2019), BMC Genomics 20:685 (PMC6717374)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6717374/
**Retrieved:** WebFetch of the URL on 2026-06-14 (search query: "SigProfilerMatrixGenerator 96 SBS reverse
complement pyrimidine reference base GitHub source").
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed) + 3 (reference implementation, the tool that produces COSMIC matrices)

**Key Extracted Points:**

1. **SBS-6 (verbatim):** "the most commonly used SBS-6 classification of single base substitutions can be
   written as: C > A, C > G, C > T, T > A, T > C, and T > G."
2. **SBS-96 construction (verbatim):** "A commonly used classification for analysis of mutational signatures
   is SBS-96, where each of the classes in SBS-6 is further elaborated using one base adjacent at the 5' of
   the mutation and one base adjacent at the 3' of the mutation. Thus, for a C > A mutation, there are
   sixteen possible trinucleotide (4 types of 5' base * 4 types of 3' base)." → 6 × 16 = 96 channels.
3. **Reverse-complement folding (verbatim):** "Please note that using the purine base of the Watson-Crick
   base-pair for classifying mutation types will require taking the reverse complement sequence of each of
   the classes of SBS-96." → mutations whose reference (mutated) base is a purine (A or G) are folded to the
   pyrimidine strand by reverse-complementing the trinucleotide context AND the substitution.
4. **Mutated base centred:** the mutated base sits in the middle of the trinucleotide (e.g. ACA > AAA, centre
   base mutated).

### Alexandrov et al. (2013) — "Signatures of mutational processes in human cancer", Nature 500:415-421

**URL:** https://www.nature.com/articles/nature12477
**Retrieved:** WebFetch of the URL on 2026-06-14 (search query: "Alexandrov 2013 Nature signatures mutational
processes human cancer 96 trinucleotide context six substitution types"). The article body redirects to an
IdP auth gate; the abstract/summary returned by the search index and the COSMIC/SigProfiler pages (above,
authored by the same group) supplied the verbatim definitions used here.
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary)

**Key Extracted Points:**

1. **96 types (verbatim from the retrieved summary):** "Each of the substitutions is examined by
   incorporating information on the bases immediately 5' and 3' to each mutated base generating 96 possible
   mutation types (6 types of substitution * 4 types of 5' base * 4 types of 3' base)."
2. **Six subtypes (verbatim from the retrieved summary):** "the six substitution subtypes: C>A, C>G, C>T,
   T>A, T>C, and T>G."
3. **Dataset scale (verbatim):** the study analysed "4,938,362 mutations from 7,042 cancers".

### Complementarity (molecular biology) — Wikipedia (Watson-Crick pairing)

**URL:** https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)
**Retrieved:** WebFetch of the URL on 2026-06-14.
**Accessed:** 2026-06-14
**Authority rank:** 4 (cites primary base-pairing chemistry)

**Key Extracted Points:**

1. **Base complements (verbatim table):** DNA bases adenine(A), thymine(T), guanine(G), cytosine(C) with
   complements "A = T, G ≡ C" — i.e. A↔T and C↔G. This is the complement map used for reverse-complement
   folding of purine-reference mutations onto the pyrimidine strand.

---

## Documented Corner Cases and Failure Modes

### From SigProfilerMatrixGenerator (PMC6717374) / COSMIC

1. **Purine reference base:** a substitution whose mutated (reference) base is a purine (A or G) is NOT one
   of the six canonical pyrimidine substitutions; it MUST be reverse-complemented to its pyrimidine equivalent
   before counting ("using the purine base ... will require taking the reverse complement sequence").
2. **Non-SBS variants:** only single-base substitutions are SBS-96 events. Indels, doublet/multi-base
   substitutions (DBS), and non-substitution variants belong to other catalogues (ID, DBS) and are not
   counted in the 96-channel SBS spectrum.

### Derived corner cases (combinatorial, not requiring separate literature)

3. **Non-ACGT context base:** a flanking base that is not A/C/G/T (e.g. N) has no defined trinucleotide
   context and cannot be classified.
4. **ref == alt:** not a substitution (no mutation); out of scope of the SBS catalogue.

---

## Test Datasets

### Dataset: Worked SBS-96 classifications (derived from the verbatim folding rules above)

**Source:** COSMIC SBS96 + SigProfilerMatrixGenerator (Bergstrom 2019) folding rule; complement map from
Watson-Crick pairing. Each row is computed independently of any implementation.

| 5' | Ref | Alt | 3' | Pyrimidine? | Folding step | Expected channel |
|----|-----|-----|----|-------------|--------------|------------------|
| A | C | A | A | yes (C) | none | A[C>A]A |
| T | C | T | G | yes (C) | none | T[C>T]G |
| G | T | C | A | yes (T) | none | G[T>C]A |
| T | G | T | A | no (G, purine) | revcomp context TGA→TCA, sub G>T → C>A | T[C>A]A |
| C | A | G | T | no (A, purine) | revcomp context CAT→ATG, sub A>G → T>C | A[T>C]G |
| G | G | C | C | no (G, purine) | revcomp context GGC→GCC, sub G>C → C>G | G[C>G]C |
| A | A | T | A | no (A, purine) | revcomp context AAA→TTT, sub A>T → T>A | T[T>A]T |

**Folding worked example (row 4), step by step:**
- Plus-strand trinucleotide 5'-T G A-3', mutation G→T.
- Reference base G is a purine → fold to pyrimidine strand by reverse-complement.
- Complement each base (A↔T, C↔G): T G A → A C T; reverse: T C A. So context becomes 5'-T C A-3'.
- Centre base G complements to C (pyrimidine); alt T complements to A → substitution C>A.
- Result: 5'=T, sub=C>A, 3'=A → **T[C>A]A**.

### Dataset: Catalog count invariant

**Source:** definition (the catalog is a partition of the input SBS variants into 96 channels).

| Property | Value |
|----------|-------|
| Number of channels | 96 (6 × 4 × 4) |
| Sum of channel counts | equals number of classifiable SBS variants |
| Channel for any input | one of the 96 canonical pyrimidine labels |

---

## Assumptions

1. **ASSUMPTION: Channel label format `5'[REF>ALT]3'`** — COSMIC/SigProfiler render the trinucleotide with
   the mutated base in the centre; the exact textual rendering (bracket form `A[C>A]A` vs underline `ACA>AAA`)
   is a display/labeling choice (non-correctness-affecting per the assumption test: the partition of variants
   into the 96 pyrimidine-keyed classes is identical either way). The bracket form `5'[REF>ALT]3'` is adopted
   for the keys; it does not change which variants fall in which class.

---

## Recommendations for Test Coverage

1. **MUST Test:** the six pyrimidine substitutions classify to themselves with the centre base unchanged
   (e.g. A[C>A]A, T[C>T]G, G[T>C]A) — Evidence: COSMIC SBS96, SigProfiler SBS-96.
2. **MUST Test:** purine-reference mutations fold by reverse-complement (T,G,T,A → T[C>A]A; the seven worked
   rows above) — Evidence: SigProfiler "using the purine base ... will require taking the reverse complement".
3. **MUST Test:** building the 96-channel catalog from a multiset of variants yields per-channel counts whose
   sum equals the number of classifiable SBS variants and whose keys are among the 96 canonical labels —
   Evidence: definition (96 = 6×4×4) + partition property.
4. **MUST Test:** the enumerated channel set is exactly the 96 canonical pyrimidine labels — Evidence: 6×4×4.
5. **SHOULD Test:** null / empty input; non-ACGT context base; ref == alt; non-single-base variant — Rationale:
   documented corner cases (purine fold, non-SBS exclusion) and standard input validation.
6. **COULD Test:** case-insensitive bases (lower-case input) classify identically — Rationale: robustness; not
   a literature requirement.

---

## References

1. Alexandrov, L.B., Nik-Zainal, S., Wedge, D.C., et al. (2013). Signatures of mutational processes in human
   cancer. Nature 500(7463):415-421. https://www.nature.com/articles/nature12477 (DOI: 10.1038/nature12477)
2. COSMIC Mutational Signatures — SBS96. Wellcome Sanger Institute.
   https://cancer.sanger.ac.uk/signatures/sbs/sbs96/ (accessed 2026-06-14)
3. Bergstrom, E.N., Huang, M.N., Mahto, U., et al. (2019). SigProfilerMatrixGenerator: a tool for visualizing
   and exploring patterns of small mutational events. BMC Genomics 20:685.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC6717374/ (DOI: 10.1186/s12864-019-6041-2)
4. Complementarity (molecular biology). Wikipedia.
   https://en.wikipedia.org/wiki/Complementarity_(molecular_biology) (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation (ONCO-SIG-001 SBS-96 trinucleotide context catalog).
