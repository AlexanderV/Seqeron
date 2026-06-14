# Evidence Artifact: MOTIF-DISCOVER-001

**Test Unit ID:** MOTIF-DISCOVER-001
**Algorithm:** Motif Discovery via Overrepresented k-mers (observed/expected enrichment)
**Date Collected:** 2026-06-14

---

## Online Sources

### Compeau & Pevzner, *Bioinformatics Algorithms* — "Expected number of occurrences of a k-mer in a DNA string"

**URL:** https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string
**Accessed:** 2026-06-14 (retrieved with WebFetch; this wiki page reproduces the formula and worked example from Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*, "Finding Regulatory Motifs in DNA Sequences")
**Authority rank:** 1 (peer-reviewed textbook content, reproduced verbatim)

**Key Extracted Points:**

1. **Background model:** Retrieved text states "the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25)." This is the zero-order i.i.d. uniform background.
2. **Probability formula:** Verbatim from the page: `Pr(N,A,Pattern,t) ≈ ( N−t·(k−1) | t )/A^(t·k)`, where N = sequence length, A = alphabet size (=4 for DNA), k = pattern length, t = number of occurrences, and `( m | n )` is a binomial coefficient.
3. **Expected count (t = 1):** For a single expected occurrence the binomial `( N−(k−1) | 1 ) = N−k+1`, so the expected number of occurrences of a specific k-mer is `(N − k + 1) / 4^k`: the number of length-k windows `N − k + 1` divided by the number `4^k` of distinct DNA k-mers.
4. **Worked example:** Retrieved text: `Pr(1000, 4, 9, 1)·500 ≈ 1.9` — a random 9-mer is expected ≈ 2 times across 500 sequences of length 1000.

### O/E ratio confirmation (search corroboration)

**URL:** WebSearch query "expected number of occurrences k-mer in random sequence overrepresented motif observed expected ratio formula" (2026-06-14); top results include monaLisa `getKmerFreq` and the PeerJ supplemental "Ratio of Observed and Expected k-mer counts".
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation / methods documentation)

**Key Extracted Points:**

1. **O/E ratio:** The retrieved search summary states the standard overrepresentation statistic is the "ratio of observed to expected frequency of a k-mer (O/E ratio), where O is the observed frequency." Enrichment > 1 ⇒ the k-mer is overrepresented.
2. **Expected = length/4^k:** The retrieved summary states "The expected number of occurrences of a k-mer is N/4^k" for the zero-order uniform model (consistent with `(N−k+1)/4^k` once the exact window count is used).

---

## Documented Corner Cases and Failure Modes

### From Compeau & Pevzner

1. **Self-overlap approximation:** The probability formula "is only an approximation because it assumes that the pattern cannot overlap with itself, which is not always the case." This affects the *probability* statistic, not the deterministic observed count or the expected-count denominator used for the O/E ratio.
2. **k > N (no windows):** When `k > N` there are zero length-k windows (`N − k + 1 ≤ 0`); no k-mer can be counted, so the discovery returns nothing.

---

## Test Datasets

### Dataset: Tandem-repeat string "ATGCATGCATGC"

**Source:** Derived directly from the Compeau & Pevzner expected-count formula `(N−k+1)/4^k`.

| Parameter | Value |
|-----------|-------|
| Sequence | ATGCATGCATGC |
| N (length) | 12 |
| k | 4 |
| Window count N−k+1 | 9 |
| Expected count E = 9/4^4 | 9/256 = 0.03515625 |
| Observed count of "ATGC" | 3 (positions 0, 4, 8) |
| Enrichment = 3 / (9/256) | 768/9 = 85.3333… |

### Dataset: Homopolymer "AAAAAAAAAA"

**Source:** Derived from `(N−k+1)/4^k`.

| Parameter | Value |
|-----------|-------|
| Sequence | AAAAAAAAAA |
| N (length) | 10 |
| k | 3 |
| Window count N−k+1 | 8 |
| Expected count E = 8/4^3 | 8/64 = 0.125 |
| Observed count of "AAA" | 8 (positions 0..7) |
| Enrichment = 8 / 0.125 | 64.0 |

---

## Assumptions

1. **ASSUMPTION: minCount filter** — The `minCount` parameter (default 2) that filters which k-mers are returned is a presentation threshold, not part of the published statistic; the formal O/E value is defined for every k-mer regardless of the cutoff. Changing it only includes/excludes rows, never alters a returned k-mer's Count, Positions, or Enrichment, so it is not correctness-affecting for any returned record.

---

## Recommendations for Test Coverage

1. **MUST Test:** Exact O/E enrichment for "ATGC" in "ATGCATGCATGC" (k=4) equals 768/9 — Evidence: Compeau & Pevzner expected-count formula.
2. **MUST Test:** Exact O/E enrichment for "AAA" in "AAAAAAAAAA" (k=3) equals 64.0 — Evidence: same formula.
3. **MUST Test:** Observed count and positions of a repeated k-mer (e.g. "ATGC" at 0/4/8) — Evidence: deterministic window enumeration.
4. **SHOULD Test:** minCount filter excludes k-mers below the threshold — Rationale: documented filter semantics.
5. **SHOULD Test:** null sequence → ArgumentNullException; k < 1 → ArgumentOutOfRangeException — Rationale: documented validation.
6. **COULD Test:** k > N returns no motifs — Rationale: no windows exist (corner case).

---

## References

1. Compeau P, Pevzner P (2015). *Bioinformatics Algorithms: An Active Learning Approach*, 2nd ed., Chapter 2 ("Which DNA Patterns Play the Role of Molecular Clocks?"). Active Learning Publishers. Formula reproduced at https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string (accessed 2026-06-14).
2. monaLisa (fmicompbio), `getKmerFreq` — observed vs expected k-mer frequencies and log2 enrichment. https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
