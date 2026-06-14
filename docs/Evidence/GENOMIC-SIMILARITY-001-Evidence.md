# Evidence Artifact: GENOMIC-SIMILARITY-001

**Test Unit ID:** GENOMIC-SIMILARITY-001
**Algorithm:** Sequence Similarity (k-mer Jaccard index)
**Date Collected:** 2026-06-14

---

## Online Sources

### Jaccard index (Wikipedia, citing Jaccard 1901 primary)

**URL:** https://en.wikipedia.org/wiki/Jaccard_index
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the primary source Jaccard 1901, rank 1)

**Retrieval:** Web search query `Jaccard index definition intersection over union 1901 Jaccard coefficient formula`, then fetched the article URL.

**Key Extracted Points:**

1. **Formula (verbatim):** "J ( A , B ) = | A ∩ B | | A ∪ B | = | A ∩ B | | A | + | B | − | A ∩ B |" — the Jaccard index is the size of the intersection divided by the size of the union of the two sample sets.
2. **Range (verbatim):** "By definition, 0 ≤ J ( A , B ) ≤ 1."
3. **Scope (verbatim):** "The Jaccard index measures similarity between finite **non-empty** sample sets." The article gives no value for the empty-union case; for measure spaces it notes: "The definition is not well-defined when μ(A ∪ B) = 0."
4. **Jaccard distance (verbatim):** "d_J(A,B) = 1 − J(A,B) = (|A ∪ B| − |A ∩ B|) / |A ∪ B|".
5. **Primary citation:** Jaccard, Paul (1901). "Étude comparative de la distribution florale dans une portion des Alpes et des Jura." *Bulletin de la Société vaudoise des sciences naturelles*, 37(142):547–579.

### Mash: fast genome and metagenome distance estimation using MinHash (Ondov et al. 2016)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Genome Biology*) — k-mer-set application of the Jaccard index

**Retrieval:** Web search `k-mer Jaccard similarity sequence MinHash Mash genome distance definition` and `Ondov Mash 2016 "Jaccard index" "fraction of shared k-mers"`, then fetched the PMC full text.

**Key Extracted Points:**

1. **k-mer set Jaccard (verbatim):** "the Jaccard index is simply the fraction of shared hashes ... out of all distinct hashes in _A_ and _B_". For exact k-mer sets this is the fraction of shared k-mers out of all distinct k-mers in the two sequences.
2. **Formal definition (verbatim):** "J(A,B) = |A∩B|/|A∪B|" over k-mer sets A and B; the MinHash sketch estimate is "j(A_s,B_s) = |A_s ∩ B_s| / s".
3. **Full citation:** Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biology*. 2016;17:132. DOI: 10.1186/s13059-016-0997-x.

### Mash documentation — Distance Estimation

**URL:** https://github.com/marbl/Mash/blob/master/doc/sphinx/distances.rst (and https://mash.readthedocs.io/en/latest/distances.html)
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference-implementation documentation)

**Retrieval:** Fetched both the GitHub source `.rst` and the rendered readthedocs page.

**Key Extracted Points:**

1. **Sketch Jaccard (verbatim):** "j(A_s,B_s) = |A_s ∩ B_s| / s" where A_s, B_s are k-mer subsets whose union equals the sketch size s — confirms Jaccard as shared-k-mers / distinct-k-mers.

---

## Documented Corner Cases and Failure Modes

### From Jaccard index (Wikipedia / Jaccard 1901)

1. **Empty union (both sets empty):** The Jaccard index is defined for non-empty sets only; it is "not well-defined when μ(A ∪ B) = 0". No authoritative value is assigned. The repository implementation returns 0 in this case (see Assumptions).
2. **Identical sets:** J = 1 (intersection = union), the maximum.
3. **Disjoint sets:** J = 0 (empty intersection, non-empty union).

### From Mash (Ondov et al. 2016)

1. **k-mer set decomposition:** distinct k-mers (a set, not a multiset) are compared; repeated k-mers within a sequence count once.
2. **Choice of k:** changes resolution of the comparison; the Jaccard formula is unchanged. Mash uses k=21 for whole genomes.

---

## Test Datasets

### Dataset: Hand-derived k-mer Jaccard worked examples (k=3)

**Source:** Derived directly from J(A,B)=|A∩B|/|A∪B| (Jaccard 1901; Ondov et al. 2016). Each k-mer set enumerated by hand and independently confirmed by a Python set computation in this session.

| seq1 | seq2 | k | A = distinct k-mers(seq1) | B = distinct k-mers(seq2) | \|A∩B\| | \|A∪B\| | J×100 |
|------|------|---|---------------------------|---------------------------|---------|---------|-------|
| ACGTACGT | ACGTACGA | 3 | {ACG,CGT,GTA,TAC} | {ACG,CGT,GTA,TAC,CGA} | 4 | 5 | 80.0 |
| ACGT | ACGA | 3 | {ACG,CGT} | {ACG,CGA} | 1 | 3 | 33.333…(100/3) |
| ACGTACGT | ACGTACGT | 3 | {ACG,CGT,GTA,TAC} | same | 4 | 4 | 100.0 |
| AAAAA | CCCCC | 3 | {AAA} | {CCC} | 0 | 2 | 0.0 |

---

## Assumptions

1. **ASSUMPTION: Empty-union return value** — When both k-mer sets are empty (both sequences empty, or both shorter than k), the union is empty and the Jaccard index is mathematically undefined (Jaccard 1901; Wikipedia: "not well-defined when μ(A ∪ B) = 0"). The implementation returns 0.0 (interpreted as "no shared content / no similarity"). This is an implementation convention, not a source-mandated value; either 0 or 1 appears in practice. Documented and tested as the implementation contract, not asserted as the literature value.
2. **ASSUMPTION: Percentage scaling (×100)** — The formal Jaccard index is in [0,1]. This method multiplies by 100 to report a percentage in [0,100]. The factor is a presentation convention, not part of the coefficient; it does not change relative ordering.
3. **ASSUMPTION: Default k = 5** — No authoritative source mandates a default k for short-DNA similarity; Mash uses k=21 for genomes. k=5 is a project default; it only sets comparison resolution, not the formula. All evidence-based tests pass k explicitly.

---

## Recommendations for Test Coverage

1. **MUST Test:** Partial overlap with exact fractional Jaccard (ACGTACGT vs ACGTACGA, k=3 → 80.0). — Evidence: Jaccard 1901; Ondov et al. 2016.
2. **MUST Test:** Identical sequences → 100.0 (J=1). — Evidence: Jaccard 1901 (max = 1).
3. **MUST Test:** Disjoint k-mer sets → 0.0 (J=0). — Evidence: Jaccard 1901 (min = 0).
4. **MUST Test:** Non-integer fraction case (ACGT vs ACGA, k=3 → 100/3). — Evidence: Jaccard formula.
5. **MUST Test:** Distinct-k-mer (set) semantics — repeated k-mers counted once (e.g. AAAAAA vs AAAA, k=3 → 100). — Evidence: Ondov et al. 2016 (distinct hashes).
6. **MUST Test:** Symmetry J(A,B)=J(B,A). — Evidence: Jaccard formula is symmetric.
7. **SHOULD Test:** Empty-union convention returns 0 (both empty; both shorter than k). — Rationale: documented implementation contract for an undefined case.
8. **SHOULD Test:** Input validation — null sequence → ArgumentNullException; kmerSize < 1 → ArgumentOutOfRangeException. — Rationale: documented failure modes.
9. **COULD Test:** Range invariant 0 ≤ result ≤ 100 over varied inputs. — Rationale: bound from Jaccard range scaled by 100.

---

## References

1. Jaccard, Paul (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles, 37(142):547–579. https://en.wikipedia.org/wiki/Jaccard_index (primary cited therein)
2. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM (2016). Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology 17:132. https://doi.org/10.1186/s13059-016-0997-x (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/)
3. Mash documentation — Distance Estimation. marbl/Mash. https://github.com/marbl/Mash/blob/master/doc/sphinx/distances.rst

---

## Change History

- **2026-06-14**: Initial documentation.
