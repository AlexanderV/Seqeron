# Evidence Artifact: KMER-GENERATE-001

**Test Unit ID:** KMER-GENERATE-001
**Algorithm:** K-mer Generation (enumerate all possible k-mers over an alphabet)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — K-mer

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14
**Authority rank:** 4 (encyclopedia article; the definition and the n^k count are standard and uncontested; used here for the foundational definition and the size of the k-mer universe)

**Retrieval:** WebSearch query `k-mer definition all possible k-mers 4^k DNA alphabet ACGT bioinformatics`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Definition:** "In bioinformatics, k-mers are substrings of length k contained within a biological sequence." (verbatim from fetched text) — a k-mer is a length-k string over the sequence alphabet.
2. **Size of the k-mer universe (formula):** "there exist n^k total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)." (verbatim) — for the DNA alphabet {A,C,G,T} the count of *all possible* k-mers is **4^k**.
3. **Worked example (AGAT):** "such that the sequence AGAT would have four monomers (A, G, A, and T), three 2-mers (AG, GA, AT), two 3-mers (AGA and GAT) and one 4-mer (AGAT)." (verbatim) — confirms a k-mer is a contiguous length-k string over {A,C,G,T}.

### BioInfoLogics — k-mer counting, part I: Introduction

**URL:** https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
**Accessed:** 2026-06-14
**Authority rank:** 3 (technical reference accompanying an established k-mer-counting toolchain; author Bernardo J. Clavijo, 2018; used for an explicit statement of the 4^k universe size)

**Retrieval:** WebSearch query `k-mer definition all possible k-mers 4^k DNA alphabet ACGT bioinformatics`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Definition:** "A k-mer is just a sequence of k characters in a string (or nucleotides in a DNA sequence)." (verbatim).
2. **All possible k-mers = 4^k:** "Since each of the k nucleotides in a k-mer can take any of the A, C, G or T values, the possible combinations of k positions are computed as 4^k." (verbatim) — independently confirms point 2 above and gives the constructive reason (each of the k positions is independently chosen from the 4-letter alphabet, i.e. the k-fold Cartesian product of the alphabet).

### Python Standard Library — itertools.product

**URL:** https://docs.python.org/3/library/itertools.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (official CPython reference documentation; the canonical reference-implementation behaviour for Cartesian-product enumeration that the generation algorithm realises)

**Retrieval:** WebSearch query `Cartesian product enumerate all strings length k over alphabet lexicographic order itertools.product Biopython k-mer`, then WebFetch of the documentation URL above.

**Key Extracted Points:**

1. **Cartesian product:** `product(A, repeat=k)` computes the k-fold Cartesian product of the alphabet `A` — i.e. all length-k tuples — which, joined into strings, is exactly the set of all possible k-mers. "To compute the product of an iterable with itself, specify the number of repetitions with the optional repeat keyword argument. For example, `product(A, repeat=4)` means the same as `product(A, A, A, A)`." (verbatim).
2. **Lexicographic ordering guarantee:** "The nested loops cycle like an odometer with the rightmost element advancing on every iteration. This pattern creates a lexicographic ordering so that if the input's iterables are sorted, the product tuples are emitted in sorted order." (verbatim) — the rightmost (last) position varies fastest; if the alphabet is supplied in sorted order, the enumerated k-mers come out in lexicographic order. {A,C,G,T} is already sorted, so DNA k-mers are emitted as AAA, AAC, AAG, AAT, ACA, … , TTT.
3. **Worked enumeration example:** "product(range(2), repeat=3) → 000 001 010 011 100 101 110 111" (verbatim) — the canonical odometer ordering over a 2-letter sorted alphabet, structurally identical to k-mer generation over {A,C,G,T}.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia / BioInfoLogics

1. **k-mer length must be positive:** a k-mer is a string "of length k"; k ≤ 0 has no meaning as a k-mer length (no length-0-or-negative substring is a k-mer in the definition).

### From itertools.product

1. **Single-letter alphabet:** with a 1-letter alphabet, the k-fold product yields exactly one k-mer (the homopolymer), since 1^k = 1.
2. **Ordering depends on alphabet order:** the lexicographic-emission guarantee holds *only if the input alphabet is sorted*; an unsorted alphabet still yields all 4^k strings but in the alphabet's own positional order, not lexicographic.

---

## Test Datasets

### Dataset: DNA k=1 (Wikipedia monomers)

**Source:** Wikipedia — K-mer (n^k formula, DNA monomers A,C,G,T)

| Parameter | Value |
|-----------|-------|
| alphabet | ACGT |
| k | 1 |
| count (4^1) | 4 |
| all k-mers (lexicographic) | A, C, G, T |

### Dataset: DNA k=2 (4^2 = 16, lexicographic Cartesian product)

**Source:** Wikipedia n^k + itertools.product odometer ordering (sorted alphabet)

| Parameter | Value |
|-----------|-------|
| alphabet | ACGT |
| k | 2 |
| count (4^2) | 16 |
| all k-mers (lexicographic, in order) | AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT |
| first | AA |
| last | TT |

### Dataset: DNA k=3 (4^3 = 64)

**Source:** Wikipedia n^k + itertools.product odometer ordering

| Parameter | Value |
|-----------|-------|
| alphabet | ACGT |
| k | 3 |
| count (4^3) | 64 |
| first (lexicographic) | AAA |
| second | AAC |
| last | TTT |
| all distinct | yes (no duplicates; |set| = 64) |

### Dataset: Protein alphabet (20^k)

**Source:** Wikipedia n^k formula generalised to an n-letter alphabet (n = 20 standard amino acids)

| Parameter | Value |
|-----------|-------|
| alphabet | ACDEFGHIKLMNPQRSTVWY (20 letters, sorted) |
| k | 2 |
| count (20^2) | 400 |

### Dataset: Single-letter alphabet (1^k)

**Source:** itertools.product (k-fold product of a 1-element pool)

| Parameter | Value |
|-----------|-------|
| alphabet | A |
| k | 4 |
| count (1^4) | 1 |
| all k-mers | AAAA |

---

## Assumptions

1. **ASSUMPTION: default alphabet = "ACGT" and case.** Sources define k-mers over the DNA alphabet {A,C,G,T} in upper case; the default `alphabet` parameter is `"ACGT"` (matching sibling `KmerAnalyzer` methods). The alphabet is used verbatim; lexicographic ordering of the output is guaranteed only when the supplied alphabet is itself in sorted order (per itertools.product). This is a documented property of the algorithm, not an invented value.

---

## Recommendations for Test Coverage

1. **MUST Test:** `GenerateAllKmers(1)` returns exactly {A, C, G, T} (4^1 = 4). — Evidence: Wikipedia n^k, DNA monomers.
2. **MUST Test:** `GenerateAllKmers(2)` returns all 16 (4^2) 2-mers in lexicographic order AA..TT. — Evidence: Wikipedia n^k + itertools.product odometer ordering on sorted alphabet.
3. **MUST Test:** `GenerateAllKmers(3)` yields 64 distinct k-mers, first AAA, last TTT. — Evidence: Wikipedia n^k; itertools.product ordering.
4. **MUST Test:** count equals 4^k for several k (1..6) with the default DNA alphabet. — Evidence: Wikipedia/BioInfoLogics 4^k formula.
5. **MUST Test:** custom alphabet count equals (alphabet length)^k (e.g. 20-letter protein alphabet, k=2 → 400). — Evidence: generalised n^k formula.
6. **SHOULD Test:** single-letter alphabet, k=4 → exactly {AAAA} (1^4 = 1). — Rationale: degenerate 1^k case (itertools.product).
7. **SHOULD Test:** output has no duplicates and is exactly the Cartesian-product set. — Rationale: each k-mer is a unique length-k tuple over the alphabet.
8. **SHOULD Test:** `k <= 0` throws ArgumentOutOfRangeException. — Rationale: k-mer length must be positive (definition).
9. **SHOULD Test:** empty/null alphabet throws ArgumentException. — Rationale: no alphabet ⇒ no k-mers can be formed.

---

## References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo BJ. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Python Software Foundation. 2026. *itertools — Functions creating iterators for efficient looping* (itertools.product). Python 3 Standard Library documentation. https://docs.python.org/3/library/itertools.html

---

## Change History

- **2026-06-14**: Initial documentation.
