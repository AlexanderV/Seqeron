# Evidence Artifact: PHYLO-BOOT-001

**Test Unit ID:** PHYLO-BOOT-001
**Algorithm:** Phylogenetic Bootstrap Analysis (Felsenstein's Bootstrap Proportions)
**Date Collected:** 2026-06-13

---

## Online Sources

### Felsenstein, J. (1985) — "Confidence Limits on Phylogenies: An Approach Using the Bootstrap" (Evolution 39(4):783–791)

**URL:** https://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.1985.tb00420.x (DOI 10.1111/j.1558-5646.1985.tb00420.x); abstract/text also via https://www.osti.gov/biblio/6044842
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper — origin of the method)

**Retrieval method:** WebSearch query `Felsenstein 1985 Confidence limits on phylogenies bootstrap approach Evolution`; the Wiley DOI page returned HTTP 403, so the abstract and procedure statements were fetched from the OSTI bibliographic record (WebFetch of https://www.osti.gov/biblio/6044842) and cross-checked against the WebSearch snippet of the Wiley/ZJU PDF.

**Key Extracted Points:**

1. **Bootstrap definition (verbatim from OSTI abstract):** "It involves resampling points from one's own data, with replacement, to create a series of bootstrap samples of the same size as the original data." — sample size of each replicate equals the original data size.
2. **What is resampled (verbatim):** "keep all of the original species while sampling characters with replacement, under the assumption that the characters have been independently drawn by the systematist and have evolved independently." — taxa (rows) are kept; characters/columns (sites) are resampled with replacement.
3. **Support proportion P:** Each replicate is analyzed and the variation among the resulting estimates indicates the error. Monophyletic groups appearing in a majority of the bootstrap samples are shown; P for a group is the fraction of bootstrap samples in which that group appears.
4. **Significance threshold (verbatim):** "If a group shows up 95% of the time or more, the evidence for it is taken to be statistically significant." — descriptive threshold; not a parameter of the computation.

### Lemoine et al. (2018) — "Renewing Felsenstein's Phylogenetic Bootstrap in the Era of Big Data" (Nature; PMC review)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6030568/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed; restates the formal FBP procedure)

**Retrieval method:** WebSearch query `Felsenstein 1985 bootstrap phylogenies pdf ...`; WebFetch of the PMC article URL.

**Key Extracted Points:**

1. **Resampling (verbatim):** "resample, with replacement, the sites of the alignment to obtain pseudo-alignments of the same length" — pseudo-alignment length equals the reference alignment length.
2. **Support computation (verbatim):** "measure the support of every branch in the reference tree as the proportion of pseudo-trees containing that branch." Support is a binary per-replicate assessment: a branch either appears in a bootstrap tree (counted) or does not.
3. **Branches as bipartitions (verbatim):** "Any branch of an X-tree defines a bipartition of X." A clade is the cluster of taxa on one side. The reference tree's internal branches are the entities scored.
4. **Term:** the values are "Felsenstein's bootstrap proportions (FBPs)" / bootstrap support values; practitioners use thresholds such as 70%.

### Biopython `Bio.Phylo.Consensus` — reference implementation (`bootstrap`, `bootstrap_trees`, `get_support`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/Phylo/Consensus.py (source); docs https://biopython.org/docs/latest/api/Bio.Phylo.Consensus.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established library, Biopython master)

**Retrieval method:** WebSearch query `Biopython Phylo Consensus bootstrap get_support ...`; WebFetch of the raw GitHub source file.

**Key Extracted Points:**

1. **Column resampling (verbatim code):** `for j in range(length): col = random.randint(0, length - 1)` then the column `alignment[:, col:col+1]` is appended. The number of resampled columns equals the original alignment length.
2. **Support counting (verbatim code):** for each non-terminal clade of the target tree a bitstring of its terminals is stored with count 0; for each bootstrap tree, every non-terminal clade whose terminal set matches increments the count: `c.confidence = (t + 1) * 100.0 / size`. Support = (number of bootstrap trees containing a clade with the identical terminal set) / (number of trees), here expressed as a percentage.
3. **Clade identity:** clades are compared by their set of terminal (leaf) names, not by branch length or internal labels. Only non-terminal clades (`find_clades(terminal=False)`) are scored — trivial single-leaf clades are not bootstrap entities.

---

## Documented Corner Cases and Failure Modes

### From Felsenstein (1985) / Lemoine et al. (2018)

1. **Reference tree provides the entity set:** support is measured only for clades present in the reference (original-data) tree; clades that appear only in bootstrap replicates are not reported. (Lemoine §FBP, point 2.)
2. **Binary per-replicate scoring:** a clade in a replicate either matches exactly (counted) or does not; there is no partial credit. (Lemoine §FBP, point 2.)
3. **Same-length resampling:** each pseudo-alignment must have the same number of columns as the original; otherwise the resampling distribution is wrong. (Felsenstein abstract; Lemoine.)

### From Biopython `Consensus.py`

1. **Clade equality by terminal set:** two clades match iff their leaf-name sets are identical; ordering and branch lengths are irrelevant.
2. **Non-trivial clades only:** terminal (leaf) clades are excluded from scoring.

---

## Test Datasets

### Dataset: Two-group deterministic alignment (reproducible support via fixed seed)

**Source:** Constructed to exercise the procedure of Felsenstein (1985); resampling made reproducible with a fixed RNG seed (mirrors Biopython's column-index resampling). Expected support values are derived analytically (see below), not copied from a code run.

| Parameter | Value |
|-----------|-------|
| Taxa | A, B, C, D |
| A | `AAAAAAAAAA` |
| B | `AAAAAAAAAA` |
| C | `GGGGGGGGGG` |
| D | `GGGGGGGGGG` |
| Alignment length | 10 |
| Tree method | UPGMA |
| Distance method | JukesCantor |
| Replicates | 100 |
| Seed | 42 (documented, fixed) |

**Derivation of expected support for clade {A,B} (and {C,D}):** Every column is either an all-`A`/all-`A` vs all-`G`/all-`G` site. After resampling *any* multiset of these 10 columns with replacement, the pairwise distances are unchanged: d(A,B)=0, d(C,D)=0, d(A,C)=d(A,D)=d(B,C)=d(B,D)=JC(p=1)=+∞ (saturated). Every bootstrap replicate therefore yields the same UPGMA topology grouping {A,B} and {C,D}. Hence support({A,B}) = support({C,D}) = 100/100 = **1.0** for all replicates and all seeds. (Felsenstein 1985: groups that appear in every replicate have P = 1.0.)

### Dataset: All-identical alignment (no non-trivial internal clades on one side / degenerate distances)

**Source:** Edge case derived from the procedure; all distances are 0.

| Parameter | Value |
|-----------|-------|
| Taxa | A, B, C |
| Each sequence | `ACGTACGT` |
| Expected | Every reported clade has support 1.0 (every replicate produces an identical zero-distance matrix → identical reference and bootstrap topology) |

---

## Assumptions

1. **ASSUMPTION: Rooted-clade scoring rather than unrooted bipartitions.** Felsenstein (1985) and Lemoine et al. (2018) describe bipartitions of unrooted trees, while the repository scores rooted clades (subtree leaf-sets) of UPGMA/NJ trees, matching Biopython's `get_support`, which compares clades by terminal set. For UPGMA (rooted, ultrametric) trees this is the conventional and consistent representation; it is the same entity Biopython scores. Documented as a modeling choice, not an unresolved correctness gap.
2. **ASSUMPTION: Support reported as a proportion in [0,1] rather than a percentage.** Felsenstein/Lemoine express P as a percentage and Biopython multiplies by 100; the repository returns the raw proportion `count/replicates ∈ [0,1]`. This is a units/labeling choice (multiply by 100 to obtain the published percentage) and does not change which clades are reported or their relative ranking.

---

## Recommendations for Test Coverage

1. **MUST Test:** Two-group alignment with a fixed seed yields support 1.0 for clades {A,B} and {C,D} (every replicate reproduces the topology). — Evidence: Felsenstein (1985) "group shows up 100% → P=1"; derivation above.
2. **MUST Test:** All reported support values lie in [0,1]. — Evidence: support = count/replicates, count ∈ [0,replicates] (Biopython `get_support`; Lemoine "proportion").
3. **MUST Test:** Determinism — same seed and inputs give identical results across runs (resampling is RNG-driven). — Evidence: bootstrap is randomized (Felsenstein 1985); reproducibility requires a fixed seed.
4. **MUST Test:** Keys of the returned dictionary equal exactly the non-trivial clades of the reference tree built from the original data. — Evidence: Biopython scores `target_tree.find_clades(terminal=False)`; Lemoine "support of every branch in the reference tree".
5. **MUST Test:** Increasing/decreasing replicates changes the denominator; a support value is always `k/replicates` for integer k. — Evidence: Biopython `(t+1)*100/size`; proportion definition.
6. **SHOULD Test:** All-identical sequences → all reported clades have support 1.0. — Rationale: degenerate zero-distance matrices reproduce one topology every replicate.
7. **SHOULD Test:** Input validation — null sequences, fewer than 2 sequences, replicates < 1 throw `ArgumentException`/`ArgumentNullException`. — Rationale: `BuildTree` requires ≥2 sequences; a bootstrap with <1 replicate is undefined.
8. **COULD Test:** Resampled alignment length equals the original alignment length (each replicate uses exactly `alignmentLength` columns). — Rationale: Felsenstein/Biopython same-length requirement.

---

## References

1. Felsenstein, J. (1985). Confidence Limits on Phylogenies: An Approach Using the Bootstrap. Evolution 39(4):783–791. https://doi.org/10.1111/j.1558-5646.1985.tb00420.x (abstract/text retrieved via https://www.osti.gov/biblio/6044842).
2. Lemoine, F., Domelevo Entfellner, J.-B., Wilkinson, E., Correia, D., Dávila Felipe, M., De Oliveira, T., Gascuel, O. (2018). Renewing Felsenstein's phylogenetic bootstrap in the era of big data. Nature 556:452–456. https://pmc.ncbi.nlm.nih.gov/articles/PMC6030568/
3. Biopython contributors. Bio.Phylo.Consensus module (`bootstrap`, `bootstrap_trees`, `get_support`). Biopython, master branch. https://raw.githubusercontent.com/biopython/biopython/master/Bio/Phylo/Consensus.py (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
