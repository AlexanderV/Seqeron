# Evidence Artifact: META-PATHWAY-001

**Test Unit ID:** META-PATHWAY-001
**Algorithm:** Metabolic Pathway Enrichment (Over-Representation Analysis via the hypergeometric test)
**Date Collected:** 2026-06-13

---

## Online Sources

### Boyle et al. (2004) — GO::TermFinder (Bioinformatics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3037731/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC article page)
**Authority rank:** 1 (peer-reviewed paper) / 3 (the reference implementation it documents)

**Key Extracted Points:**

1. **Exact p-value formula (verbatim):** The over-representation p-value is the right tail of the
   hypergeometric distribution: `P = 1 − Σ_{i=0}^{k−1} C(M,i)·C(N−M, n−i) / C(N, n)`.
2. **Symbol definitions (verbatim):** *N* = "total number of genes in the background distribution";
   *M* = "number of genes within that distribution that are annotated (either directly or indirectly)
   to the node of interest" (the pathway/term size); *n* = "size of the list of genes of interest"
   (the query); *k* = "number of genes within that list which are annotated to the node" (the overlap).
3. **One-sided / right-tail semantics:** the P-value is "the probability of *x* or more out of *n*
   genes having a given annotation" — i.e. P(X ≥ k), an over-representation (upper-tail) test.
4. **Distribution choice:** the hypergeometric distribution (sampling *without* replacement) is used
   because it is more accurate than the binomial (sampling with replacement).
5. **Worked intuition example:** a bag of 500 red + 500 green beads; drawing 20 without replacement
   and observing 17 green → P = probability of picking 17 or more green from 20. (No numeric result
   is given in the paper; used here only to confirm the right-tail "k or more" semantics.)

### PNNL — Proteomics Data Analysis in R/Bioconductor, §8.2 Over-Representation Analysis

**URL:** https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
**Accessed:** 2026-06-13 (retrieved via WebFetch of the tutorial page)
**Authority rank:** 3 (documents the canonical R/Bioconductor `phyper` reference computation)

**Key Extracted Points:**

1. **Formula (verbatim):** `P(X ≥ x) = 1 − P(X ≤ x−1) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n)`
   — identical to the GO::TermFinder formula.
2. **Symbols:** *N* = number of background genes; *n* = number of "interesting" genes (query);
   *M* = number of genes annotated to gene set *S* (pathway size); *x* = number of interesting genes
   annotated to *S* (overlap).
3. **Reference R call:** `phyper(q = x − 1, m = M, n = N − M, k = n, lower.tail = FALSE)` computes the
   upper tail P(X ≥ x). The `q = x − 1` shift is what makes `lower.tail = FALSE` give P(X ≥ x) instead
   of P(X > x).
4. **Worked numeric example (verbatim values):** background N = 8000; query/cluster *C* of n = 400;
   gene set *S* of M = 100 annotated genes; overlap x = 20. Reported result **P(X ≥ 20) = 7.88 × 10⁻⁸**.

---

## Documented Corner Cases and Failure Modes

### From PNNL §8.2 / GO::TermFinder

1. **x = 0 (no overlap):** the sum Σ_{i=0}^{x−1} is empty, so P = 1 − 0 = 1. No over-representation
   is possible, so the largest p-value (1.0) is correct.
2. **Empty / degenerate population (N = 0, M = 0, or n = 0):** no successes can be drawn; the
   probability of "at least one" success is 0, hence P(X ≥ x) = 1 for any x.
3. **Sampling-without-replacement constraint:** terms with i > M or (n − i) > (N − M) are infeasible
   (C(·,·) = 0) and contribute 0 to the cumulative sum.

---

## Test Datasets

### Dataset: PNNL §8.2 worked example

**Source:** PNNL Proteomics Data Analysis in R/Bioconductor §8.2 (URL above).

| Parameter | Value |
|-----------|-------|
| Background size N | 8000 |
| Pathway / gene-set size M | 400 (symmetric with the query/cluster split; see note) |
| Query size n | 100 |
| Overlap x | 20 |
| Expected P(X ≥ 20) | 7.88 × 10⁻⁸ |

> Note: the hypergeometric is symmetric under swapping the roles of the two marked groups, so the
> PNNL labelling (cluster *C* = 400, gene set *S* = 100) and the GO::TermFinder labelling
> (pathway *M*, query *n*) give the same p-value whether (M, n) = (400, 100) or (100, 400).
> Both orientations were independently computed and yield 7.884747×10⁻⁸.

### Dataset: Small exact hand-derived cases (rational arithmetic)

**Source:** Direct evaluation of the formula above with exact binomial coefficients.

| Case | N | M | n | x | P(X ≥ x) (exact) |
|------|---|---|---|---|------------------|
| All query in pathway | 10 | 5 | 5 | 5 | 1/252 = 0.003968253968… |
| Partial overlap | 4 | 2 | 2 | 1 | 5/6 = 0.833333333… |
| No overlap (x = 0) | 10 | 5 | 5 | 0 | 1 |
| At least one | 10 | 5 | 5 | 1 | 251/252 = 0.996031746… |

Derivation of 1/252: P(X = 5) = C(5,5)·C(5,0)/C(10,5) = 1·1/252. P(X ≥ 5) = P(X = 5) = 1/252.
Derivation of 5/6: P(X = 0) = C(2,0)·C(2,2)/C(4,2) = 1/6 ⇒ P(X ≥ 1) = 1 − 1/6 = 5/6.

---

## Assumptions

1. **ASSUMPTION: Background defaulting** — When the caller supplies no background gene universe, the
   implementation uses the union of all pathway members (plus the query) as the background. The cited
   sources require a *caller-defined* background ("total number of genes in the background distribution");
   they do not prescribe a default. This default is a convenience that does not change the formula and
   is documented as caller-overridable. Non-correctness-affecting for the formula itself when a
   background is supplied; correctness-affecting only for the no-background convenience path, which is
   explicitly tested against a hand-derived expectation.

---

## Recommendations for Test Coverage

1. **MUST Test:** PNNL §8.2 worked example reproduces P(X ≥ 20) = 7.88 × 10⁻⁸ — Evidence: PNNL §8.2.
2. **MUST Test:** Exact small case all-overlap → 1/252 — Evidence: formula (Boyle 2004 / PNNL §8.2).
3. **MUST Test:** Exact small case partial overlap → 5/6 — Evidence: formula (Boyle 2004 / PNNL §8.2).
4. **MUST Test:** x = 0 → p = 1 (empty upper-sum corner case) — Evidence: PNNL §8.2 corner case 1.
5. **MUST Test:** Right-tail monotonicity: result fields (overlap/sizes) match inputs; results sorted
   ascending by p-value — Evidence: Boyle 2004 (P-value ranking of terms).
6. **SHOULD Test:** Symmetry — (M, n) and (n, M) give the same p-value — Rationale: hypergeometric symmetry.
7. **SHOULD Test:** Caller-supplied vs defaulted background changes N and hence p — Rationale: ASSUMPTION path.
8. **COULD Test:** Many pathways processed in one call, each independent — Rationale: API shape.

---

## References

1. Boyle EI, Weng S, Gollub J, Jin H, Botstein D, Cherry JM, Sherlock G (2004). GO::TermFinder — open
   source software for accessing Gene Ontology information and finding significantly enriched Gene
   Ontology terms associated with a list of genes. *Bioinformatics* 20(18):3710–3715.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3037731/ (DOI: 10.1093/bioinformatics/bth456)
2. PNNL Computational Mass Spectrometry. Proteomics Data Analysis in R/Bioconductor, §8.2
   Over-Representation Analysis. https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
   (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
