# Evidence Artifact: KMER-UNIQUE-001

**Test Unit ID:** KMER-UNIQUE-001
**Algorithm:** Unique K-mers / K-mers with Minimum Count (k-mer frequency filtering)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — K-mer

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14
**Authority rank:** 4 (encyclopedia article; used for the foundational definition and worked example, which are standard and uncontested)

**Retrieval:** WebSearch query `k-mer definition substring length k overlapping bioinformatics number of k-mers n-k+1`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Definition:** "In bioinformatics, k-mers are substrings of length k contained within a biological sequence." (verbatim from fetched text)
2. **Total-count formula:** A sequence of length L contains **L − k + 1** total k-mers (overlapping, step size +1). The number of *possible* k-mers over an alphabet of n monomers is **n^k** (n = 4 for DNA).
3. **Worked example (AGAT):** monomers (k=1): A, G, A, T → 4; 2-mers: AG, GA, AT → 3; 3-mers: AGA, GAT → 2; 4-mers: AGAT → 1. Confirms overlapping extraction with single-position advancement.
4. **Overlap:** k-mers are extracted by a sliding window of step 1, so adjacent k-mers overlap by k−1 characters.

### BioInfoLogics — k-mer counting, part I: Introduction

**URL:** https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
**Accessed:** 2026-06-14
**Authority rank:** 3 (technical reference accompanying an established k-mer-counting toolchain; used for the precise total/distinct/unique terminology and a worked frequency table)

**Retrieval:** WebSearch query `singleton unique k-mers appearing exactly once sequence definition`, then WebFetch of the article URL above.

**Key Extracted Points:**

1. **Total k-mers:** the sum of all k-mers extracted from a sequence counting duplicates; for length L there are L − k + 1 of them.
2. **Distinct k-mers:** "_Distinct k-mers_ are counted only once, even if they appear more times." (verbatim) — i.e. the number of different k-mer strings present.
3. **Unique k-mers:** "_Unique k-mers_ are those that appear only once." (verbatim) — k-mers whose frequency is exactly 1. This is the canonical definition `FindUniqueKmers` must implement.
4. **Worked example (ATCGATCAC, k=3, non-canonical):** 7 total 3-mers, 6 distinct, 5 unique. The 5 unique 3-mers occurring exactly once are TCG, CGA, GAT, TCA, CAC; ATC occurs twice (positions 0 and 4) and is therefore NOT unique.

### Compeau & Pevzner — Bioinformatics Algorithms: An Active Learning Approach (k-mer counting / Frequent Words)

**URL:** https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637 (catalogue/description page surfaced for the textbook)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed textbook)

**Retrieval:** WebSearch query `Compeau Pevzner Bioinformatics Algorithms k-mer Count Pattern Text frequency definition L-k+1`.

**Key Extracted Points:**

1. **Count(Text, Pattern):** the number of times a k-mer `Pattern` appears as a substring of `Text` (overlapping occurrences counted). This is the per-k-mer frequency on which both "unique" (Count = 1) and "min-count" (Count ≥ t) filters operate.
2. **Most-frequent / recurrent k-mers:** a k-mer is a most frequent k-mer if it maximises Count(Text, Pattern); selecting k-mers whose Count ≥ a threshold t is the standard way to isolate recurrent/over-represented k-mers (the basis for `FindKmersWithMinCount`).

---

## Documented Corner Cases and Failure Modes

### From Wikipedia — K-mer

1. **k > L:** when the k-mer length exceeds the sequence length, L − k + 1 ≤ 0, so the sequence contains zero k-mers.

### From BioInfoLogics

1. **Repeated k-mers excluded from "unique":** a k-mer that appears two or more times (e.g. ATC in ATCGATCAC) is distinct but NOT unique; the unique set is strictly the frequency-1 subset.

---

## Test Datasets

### Dataset: ATCGATCAC (BioInfoLogics worked example)

**Source:** BioInfoLogics — k-mer counting, part I (worked table, non-canonical k=3)

| Parameter | Value |
|-----------|-------|
| Sequence | ATCGATCAC |
| k | 3 |
| Total 3-mers | 7 (ATC, TCG, CGA, GAT, ATC, TCA, CAC) |
| Distinct 3-mers | 6 |
| Unique 3-mers (Count = 1) | 5 → {TCG, CGA, GAT, TCA, CAC} |
| Non-unique 3-mers | ATC (Count = 2) |

### Dataset: ACGTACGT (derived from the k-mer definition for min-count filtering)

**Source:** k-mer definition (Wikipedia L−k+1 + Compeau & Pevzner Count); occurrences enumerated directly.

| Parameter | Value |
|-----------|-------|
| Sequence | ACGTACGT |
| k | 4 |
| Total 4-mers | 5 (ACGT, CGTA, GTAC, TACG, ACGT) |
| Counts | ACGT=2, CGTA=1, GTAC=1, TACG=1 |
| FindKmersWithMinCount(..,4,2) | {(ACGT, 2)} |
| FindKmersWithMinCount(..,4,1) | all 4 distinct, ordered by Count desc (ACGT first, Count 2) |
| FindUniqueKmers(..,4) | {CGTA, GTAC, TACG} (the Count=1 set) |

### Dataset: AGAT (Wikipedia worked example, all distinct)

**Source:** Wikipedia — K-mer (AGAT example)

| Parameter | Value |
|-----------|-------|
| Sequence | AGAT |
| k | 2 |
| 2-mers | AG, GA, AT (all Count = 1) |
| FindUniqueKmers(AGAT, 2) | {AG, GA, AT} (all three unique) |

---

## Assumptions

1. **ASSUMPTION: minCount ≤ 0 behaviour** — Authoritative sources define min-count filtering only for a meaningful threshold (t ≥ 1, recurrent k-mers). For minCount ≤ 1 the predicate `Count ≥ minCount` is satisfied by every observed k-mer; the implementation returns all distinct k-mers ordered by count. This is the mathematically consistent extension of `Count ≥ t`, not an invented value, so it is treated as defined behaviour rather than a correctness-affecting unknown.
2. **ASSUMPTION: case normalisation** — sources use upper-case DNA; the implementation upper-cases input (consistent with sibling `KmerAnalyzer` methods) so that case variants count as the same k-mer. No source contradicts this.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindUniqueKmers(ATCGATCAC, 3)` returns exactly {TCG, CGA, GAT, TCA, CAC} (Count = 1 set; ATC excluded). — Evidence: BioInfoLogics worked table.
2. **MUST Test:** `FindUniqueKmers(AGAT, 2)` returns all three 2-mers {AG, GA, AT}. — Evidence: Wikipedia AGAT example.
3. **MUST Test:** `FindKmersWithMinCount(ACGTACGT, 4, 2)` returns exactly {(ACGT, 2)}. — Evidence: definition + enumerated occurrences.
4. **MUST Test:** `FindKmersWithMinCount(ACGTACGT, 4, 1)` returns all 4 distinct 4-mers ordered by Count descending, ACGT (Count 2) first. — Evidence: Compeau & Pevzner Count ≥ t + ordering.
5. **MUST Test:** homopolymer (e.g. AAAAA, k=3) has zero unique k-mers (the single distinct 3-mer AAA has Count = 3). — Evidence: definition (Count > 1 ⇒ not unique).
6. **SHOULD Test:** empty sequence and k > length return empty for both methods. — Rationale: L − k + 1 ≤ 0 ⇒ no k-mers (Wikipedia).
7. **SHOULD Test:** k ≤ 0 throws ArgumentOutOfRangeException. — Rationale: k-mer length must be positive (definition: substrings of length k).
8. **COULD Test:** case-insensitivity (lower-case input yields same unique set). — Rationale: documented normalisation assumption.

---

## References

1. Wikipedia contributors. 2026. *K-mer*. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Bernardo Clavijo et al. 2018. *k-mer counting, part I: Introduction*. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach* (2nd ed.). Active Learning Publishers. https://www.amazon.com/BIOINFORMATICS-ALGORITHMS-Phillip-Compeau/dp/0990374637

---

## Change History

- **2026-06-14**: Initial documentation.
