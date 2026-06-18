# Evidence Artifact: PANGEN-CORE-001

**Test Unit ID:** PANGEN-CORE-001
**Algorithm:** Core / Accessory / Unique genome construction (pan-genome partitioning), genome fluidity, open/closed classification
**Date Collected:** 2026-06-13

---

## Online Sources

### Tettelin et al. (2005) — original "pan-genome" definition (PNAS)

**URL:** https://www.pnas.org/doi/10.1073/pnas.0506758102 (DOI 10.1073/pnas.0506758102)
**Retrieved via:** WebSearch query `Tettelin 2005 pan-genome core genome accessory genome definition PNAS Streptococcus` (the publisher page returned HTTP 403 to WebFetch; the definitions below were extracted from the indexed abstract/summary returned by the search and re-confirmed against Tettelin 2008 and the Wikipedia Pan-genome article, both opened below).
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Pan-genome:** The set of all unique genes across a set of genomes — the union of all genes (global gene repertoire).
2. **Core genome:** The set of genes found in *every* genome — the intersection.
3. **Accessory / dispensable genome:** Genes present in some but not all genomes; genes unique to a single genome are strain-specific.

---

### Tettelin, Riley, Cattuto, Medini (2008) — "Comparative genomics: the bacterial pan-genome" (Curr Opin Microbiol)

**URL:** https://www.sciencedirect.com/science/article/abs/pii/S1369527408001239 (DOI 10.1016/j.mib.2008.09.006; PubMed 19086349)
**Retrieved via:** WebSearch query `Tettelin 2008 "comparative genomics" pan-genome definition core dispensable genome "open pan-genome" "closed pan-genome" Curr Opin Microbiol`
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed review)

**Key Extracted Points:**

1. **Core + dispensable decomposition:** A microbial pan-genome is the combination of a *core* genome (genes present in all strains) and a *dispensable* / flexible / accessory genome (genes absent from one or more strains).
2. **Open vs closed:** The pan-genome is classified as *open* when the Heaps'-law power-law exponent indicates unbounded growth (alpha < 1); many bacterial species have open pan-genomes — a small conserved core plus a large dynamic accessory pool.

---

### Wikipedia — "Pan-genome" (citing primaries Tettelin 2005/2008)

**URL:** https://en.wikipedia.org/wiki/Pan-genome
**Retrieved via:** WebFetch of the article URL
**Accessed:** 2026-06-13
**Authority rank:** 4 (encyclopedic, used only to corroborate the primary definitions above)

**Key Extracted Points:**

1. **Components:** core genome = genes present in all individuals; shell/accessory = genes shared by two or more strains; cloud / strain-specific = genes in a single strain (also "dispensable").
2. **Heaps' law classification:** pan-genomes are classified using Heaps' law `N = k·n^(-alpha)`. **Open** pan-genome when **alpha ≤ 1** (genes keep accumulating, no asymptote); **closed** when **alpha > 1** (few new genes per genome; size approaches a limit).
3. **Examples:** *E. coli* — open (~89,000 gene families across ~2,000 genomes); *S. pneumoniae* — closed (negligible new genes after ~50 genomes).

---

### Kislyuk, Haegeman, Bergman, Weitz (2011) — "Genomic fluidity" (BMC Genomics 12:32)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3030549/ (DOI 10.1186/1471-2164-12-32; PubMed 21232151)
**Retrieved via:** WebSearch `Kislyuk 2011 genome fluidity formula unique gene families pan-genome BMC Genomics` → WebFetch of the PMC open-access full text.
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Fluidity formula (verbatim symbols):** `φ = [2 / (N(N−1))] · Σ_{k<l} (U_k + U_l) / (M_k + M_l)`, where `N` = number of genomes, `U_k, U_l` = number of gene families found only in genome k and only in genome l respectively, `M_k, M_l` = total number of gene families in k and l respectively.
2. **Range / interpretation:** fluidity ranges 0..1. `φ = 0` ⇒ genomes identical at the gene level; `φ = 1` ⇒ complete dissimilarity. A fluidity of 0.1 means a pair of genomes have on average 10% unique genes and share 90%.
3. **Variance (jackknife):** `σ² = ((N−1)/N) · Σ_i (φ̂_(i) − φ̂)²`, where `φ̂_(i)` is fluidity estimated from genome pairs not including genome i.

---

### Page et al. (2015) — Roary pan-genome pipeline (Bioinformatics 31(22):3691–3693)

**URL:** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4817141/ (DOI 10.1093/bioinformatics/btv421); README https://raw.githubusercontent.com/sanger-pathogens/Roary/master/README.md
**Retrieved via:** WebSearch `Roary pan-genome core genes 99% definition soft core shell cloud genes Page 2015 Bioinformatics` → WebFetch of the PMC paper and the GitHub raw README.
**Accessed:** 2026-06-13
**Authority rank:** 3 (established reference implementation)

**Key Extracted Points:**

1. **Core definition:** "Core is defined as a gene being in at least 99% of samples, which allows for some assembly errors in very large datasets." Default core threshold flag `-cd [99]`.
2. **Clustering identity:** all-against-all BLASTP with "a user defined percentage sequence identity (default 95%)"; README usage `-i minimum percentage identity for blastp [95]`.
3. **Four-tier categories (from the search-indexed Page et al. classification):** hard core >99% of genomes; soft core 95–99%; shell 15–95%; cloud <15%.

---

### micropan (CRAN) — `heaps()` reference implementation of Heaps-law openness

**URL:** https://rdrr.io/cran/micropan/man/heaps.html and https://search.r-project.org/CRAN/refmans/micropan/html/heaps.html
**Retrieved via:** WebSearch `Kislyuk ... micropan fluidity` and `Tettelin Heaps law ... micropan` → WebFetch of both CRAN doc mirrors.
**Accessed:** 2026-06-13
**Authority rank:** 3 (established reference implementation, cites Tettelin 2008)
**Key Extracted Points:**

1. **Model:** "The Heaps law model is fitted to the number of *new* gene clusters observed when genomes are ordered in a random way." Returns the Intercept and the **decay parameter alpha**.
2. **Openness criterion (verbatim):** "If `alpha<1.0` the pan-genome is open, if `alpha>1.0` it is closed."
3. **Fluidity (companion `fluidity()` doc):** genomic fluidity between two genomes = number of unique gene families divided by total number of gene families, averaged over random pairs (corroborates Kislyuk formula).

---

## Documented Corner Cases and Failure Modes

### From Tettelin (2008) / micropan

1. **Too few genomes for openness fit:** the Heaps-law decay exponent is only meaningful once several genomes are accumulated; with < 3 genomes the new-gene curve is degenerate (micropan fits on the random-order new-cluster curve).

### From Kislyuk (2011)

1. **Single genome / no pairs:** fluidity is defined as an average over genome *pairs*; with N < 2 there are no pairs and fluidity is undefined (taken as 0 by convention here).
2. **Empty gene sets:** if `M_k + M_l = 0` for a pair the term is undefined; such pairs contribute 0.

### From Page et al. (2015)

1. **Soft core vs hard core:** the 99% core threshold tolerates assembly error. Roary's rule is fractional — "a gene being in at least 99% of samples", i.e. `occupancy / N ≥ coreFraction` — **not** `floor(coreFraction · N)`. With few genomes the fractional rule is strict (N=3, 0.99 ⇒ only 3/3 is core; 2/3 = 66.7% is shell/accessory). (Corrected 2026-06-15 during PANGEN-CORE-001 validation; see report.)

---

## Test Datasets

### Dataset: Kislyuk worked interpretation (fluidity)

**Source:** Kislyuk et al. (2011), BMC Genomics 12:32.

| Parameter | Value |
|-----------|-------|
| Two genomes, 10% unique each, 90% shared | fluidity = 0.1 |
| Identical gene content | fluidity = 0 |
| Disjoint gene content | fluidity = 1 |

**Hand-derived three-genome example (from the fluidity equation, used in tests):**

Genome A clusters = {c1, c2, c3}; B = {c1, c2, c4}; C = {c1, c5, c6}, where c1 is shared by all, c2 by A and B.

| Pair | U_k+U_l (symmetric diff size) | M_k+M_l (total) | term |
|------|------|------|------|
| A,B | 2 (c3, c4) | 6 | 2/6 = 0.333333… |
| A,C | 4 (c2,c3,c5,c6) | 6 | 4/6 = 0.666667… |
| B,C | 4 (c2,c4,c5,c6) | 6 | 4/6 = 0.666667… |

φ = (1/3)·(2/6 + 4/6 + 4/6) = (1/3)·(10/6) = 10/18 = **0.5555555555…**

### Dataset: Core/accessory/unique partition (Tettelin/Roary semantics)

**Source:** Tettelin (2005, 2008); Page et al. (2015).

3 genomes; cluster occupancy: c1 in all 3 (core), c2 in 2 (accessory/shell), c3,c4,c5 each in 1 (unique/cloud).
With coreFraction = 1.0 → coreThreshold = 3: core = {c1}; unique (occupancy 1) = {c3,c4,c5}; accessory = {c2}.

---

## Assumptions

1. **ASSUMPTION: clustering identity metric.** Roary/Tettelin use BLASTP-based percentage identity for clustering; this repository's `ConstructPanGenome` delegates to the in-repo `ClusterGenes` (k-mer Jaccard heuristic, threshold default 0.9). The *partitioning* logic under test (core/accessory/unique by cluster occupancy, fluidity, openness) is independent of the upstream identity metric, so test inputs use identical or fully-disjoint sequences where the occupancy is unambiguous. The clustering metric itself is the subject of PANGEN-CLUSTER-001, not this unit.
2. **ASSUMPTION: empty-pair convention for fluidity.** Pairs whose `M_k + M_l = 0` contribute 0 (the equation is undefined; 0 is the neutral element). Not stated explicitly by Kislyuk, but only arises for empty genomes.

---

## Recommendations for Test Coverage

1. **MUST Test:** core/accessory/unique partition by cluster occupancy with the fractional rule `occupancy / N ≥ coreFraction` (Roary "at least 99% of samples", not `floor(coreFraction · N)`) — Evidence: Tettelin (2005, 2008), Page et al. (2015).
2. **MUST Test:** genome fluidity equals the closed-form `φ = [2/(N(N−1))]·Σ_{k<l}(U_k+U_l)/(M_k+M_l)` on the hand-derived 3-genome example (= 0.5̄) — Evidence: Kislyuk (2011).
3. **MUST Test:** fluidity bounds — identical gene content → 0; disjoint gene content → 1 — Evidence: Kislyuk (2011).
4. **MUST Test:** open vs closed classification by Heaps'-law decay exponent alpha (open ⟺ alpha < 1) — Evidence: Tettelin (2008), micropan.
5. **SHOULD Test:** empty input → empty result; single genome → no pairs (fluidity 0), occupancy-1 clusters are unique — Rationale: documented degenerate cases.
6. **COULD Test:** core fraction = CoreGeneCount / TotalClusters invariant — Rationale: definitional consistency.

---

## References

1. Tettelin H, et al. (2005). Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome". *PNAS* 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102
2. Tettelin H, Riley D, Cattuto C, Medini D (2008). Comparative genomics: the bacterial pan-genome. *Curr Opin Microbiol* 11(5):472–477. https://doi.org/10.1016/j.mib.2008.09.006
3. Kislyuk AO, Haegeman B, Bergman NH, Weitz JS (2011). Genomic fluidity: an integrative view of gene diversity within microbial populations. *BMC Genomics* 12:32. https://doi.org/10.1186/1471-2164-12-32
4. Page AJ, et al. (2015). Roary: rapid large-scale prokaryote pan genome analysis. *Bioinformatics* 31(22):3691–3693. https://doi.org/10.1093/bioinformatics/btv421
5. Lagesen K, et al. micropan: Microbial Pan-Genome Analysis — `heaps()`/`fluidity()` reference docs. CRAN. https://rdrr.io/cran/micropan/man/heaps.html
6. Wikipedia contributors. Pan-genome. https://en.wikipedia.org/wiki/Pan-genome (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation (PANGEN-CORE-001).
