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
- **2026-11-18** (B09 re-validation): see "B09 re-validation oracles" below.

---

## B09 re-validation oracles (2026-11-18)

Sources actually opened this session (en.wikipedia.org / ncbi are blocked; the Jaccard 1901 definition and the Ondov et al.
text above were taken from the earlier fetch):

- Mash `doc/sphinx/sketches.rst` and `src/mash/Sketch.cpp` (raw.githubusercontent.com/marbl/Mash/master, neutral User-Agent):
  "will ignore strandedness by using canonical k-mers ... Strandedness can be preserved with -n (non-canonical)"; `Sketch.cpp`
  `kmer = (noncanonical || memcmp(kmer_fwd, kmer_rev, k) <= 0) ? kmer_fwd : kmer_rev`. Canonical k-mers are therefore the Mash default.
- scipy 1.18.1 `scipy.spatial.distance.jaccard` (installed package) and sourmash 4.9.4 `MinHash.jaccard` (PyPI, installed).

Oracle `scripts/analysis/generate_similarity_oracle.py` -> `TestData/Analysis/similarity_oracle.tsv` (2086 rows = 164 pairs
x 7 k). Per row it asserts (a) Python set J == scipy `1 - jaccard(presence vectors)` and (b) Python set over canonical
k-mers (Biopython `Seq.reverse_complement`) == sourmash exact Jaccard (`scaled=1`); the TSV stores literal and canonical percent.
924 of 2086 rows differ between literal and canonical, so the file discriminates the two modes.

Spot values (k = 5): ATGCATGC vs ATGCATTT -> 100*2/6 = 33.333333333333336 (the MCP doc's former "0.75" was wrong); AACCCGGG vs reverse
complement CCCGGGTT, k = 3: literal 50.0, canonical 100.0.

Empty-union conventions of the references: scipy `jaccard([], [])` = 0.0 (distance, i.e. J = 1), sourmash `jaccard` of two empty
sketches = 0.0, Mash `shared == union` -> distance 0 (J = 1). No source fixes J for an empty union; the library keeps 0 (ASM-1).

## 2026-11-19 audit round 2 (A2-31 / A2-35, B09 F95): sketch estimators, containment, Jaccard distance, ANI

- **Source definitions.** Ondov et al. 2016 (Mash, Genome Biol 17:132): bottom-s MinHash, J estimate x/s, `D = -ln(2J/(1+J))/k`; Irber et al. 2022 / Pierce et al. 2019 (sourmash FracMinHash `scaled`) and Koslicki & Zabeti 2019 (containment index); Hera, Pierce, Koslicki et al. 2023 (Genome Res, "Debiasing FracMinHash and deriving confidence intervals for mutation rates...") for the containment / Jaccard ANI point estimates `C^(1/k)` and `(2J/(1+J))^(1/k)` as implemented by sourmash 4.9.4 `distance_utils.containment_to_distance` / `jaccard_to_distance` (the installed source was read: `point_estimate = 1 - containment ** (1/ksize)`; 0 -> distance 1, 1 -> distance 0).
- **Oracle numbers.** `scripts/analysis/generate_similarity_sketch_oracle.py`: 4290 rows — 990 exact (Python sets of canonical k-mers, asserted equal to sourmash `MinHash(scaled=1)` Jaccard and both containments; `containment_to_distance` ANI both directions; Jaccard ANI formula asserted against `jaccard_to_distance`), 1980 FracMinHash (scaled 2/3/5/10, sourmash's own jaccard / contained_by / max_containment, k = 3 … 31), 1320 bottom-s (`MinHash(n=5..500)`, k = 17 … 31; Mash hashes 64 bits there, 32 bits for k <= 16 where only the Mash definition applies and the canonical `KmerAnalyzer` sketch is the reference). All equal within 1e-12.
- **Implementation.** `GenomicAnalyzer.Similarity.cs` selects the estimator and converts; sketches and distances are `KmerAnalyzer` (B06) calls (no copy of the hash, sketch or Mash formula); only the two one-line ANI conversions are new.
- **Tests.** `GenomicAnalyzer_SimilaritySketch_Oracle_Tests` (6), MCP `calculate_similarity` (`sketchSize`, `scaled`, the added fields; SuffixTree.Mcp.Core tests).

## 2026-11 audit round 3 (A3-11, B09 F138): ANI confidence intervals and the `mash dist` comparison

Source: sourmash 4.9.4 `sourmash/distance_utils.py` (installed Python package, read in full: `containment_to_distance`, `jaccard_to_distance`, `var_n_mutated`, `exp_n_mutated`, `get_exp_probability_nothing_common`, `set_size_exact_prob`) and `minhash.py` (`containment_ani`, `max_containment_ani`, `jaccard_ani`, `size_is_accurate`: n_unique_kmers = len(sketch) x scaled; jaccard: round((len A + len B) / 2 x scaled)); equations from github.com/KoslickiLab/mutation-rate-ci-calculator, Hera et al. 2023 (Genome Res). Ported with scipy's own algorithms: `brentq` (scipy/optimize/Zeros/brentq.c, bracket [1e-7, 1 - 1e-7], xtol 2e-12, rtol 4 eps), `norm.ppf` (Wichura AS 241), `binom.cdf` (regularized incomplete beta, existing `StatisticsHelper.RegularizedIncompleteBeta`). Oracle `generate_ani_ci_oracle.py`: 280 random pairs (400-4000 nt, 0.5-25 % mutation, reverse complements, unrelated, N) x k in 5..31 x scaled in 1..20 x confidence 0.90 / 0.95 / 0.99 x {containment A in B, B in A, max containment, Jaccard} = **1120 rows**: distance, interval, probability of nothing in common (and its flag), size-inaccuracy flag, Jaccard error (and its flag) equal sourmash on all rows (interval 560+ rows with a CI; 252 size-inaccurate; probability-threshold flags on both values). scipy values of probit / brentq / binom.cdf locked in the unit tests.

`mash dist` 2.3 (Debian mash 2.3+dfsg-6, installed): `generate_mash_dist_oracle.py` -> **400 rows** (200-6000 nt pairs, k in {3, 5, ..., 31}, s in {10, ..., 1000}, canonical and `-n`, N runs, related at 0-30 % and unrelated): shared hashes and denominator equal exactly, distance and p-value equal to the 6 printed significant digits (p-values down to 1e-300). The earlier "Not provided: mash p-value" note was a missing oracle, not a missing formula: `MashPValue` = Mash `pValue` already agreed.
