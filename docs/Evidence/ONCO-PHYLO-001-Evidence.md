# Evidence Artifact: ONCO-PHYLO-001

**Test Unit ID:** ONCO-PHYLO-001
**Algorithm:** Tumor Phylogeny Reconstruction — clonal tree from CCF clusters (sum rule + lineage-precedence rule)
**Date Collected:** 2026-06-15

---

## Online Sources

### Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S (2015) — LICHeE

**URL:** https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0647-8 (DOI 10.1186/s13059-015-0647-8); full text PDF retrieved at https://arxiv.org/pdf/1412.8574
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper, *Genome Biology*)

**How retrieved:** WebSearch query "LICHeE Popic 2015 Genome Biology fast scalable inference multi-sample cancer lineages VAF perfect phylogeny constraints" → fetched the arXiv preprint PDF (arxiv.org/pdf/1412.8574); the PDF was parsed (pypdf) and the constraint passages quoted below were extracted from the rendered text.

**Key Extracted Points:**

1. **Perfect-phylogeny ordering constraints (verbatim):** "Firstly, (1) a mutation present in a given set of samples cannot be a successor of a mutation that is present in a smaller subset of these samples … (2) a given mutation cannot have a VAF higher than that of its predecessor mutation (except due to CNVs), since all cells containing this mutation will also contain the predecessor. Finally, (3) the sum of the VAFs of mutations disjointly present in distinct subclones cannot exceed the VAF of a common predecessor mutation present in these subclones, since the subclones with the descendent mutations must contain the parent mutations."
2. **Ancestor ≥ descendant edge rule (verbatim, Eq. 2):** an edge `(u → v)` is added only if for all samples i: "(1) `u.VAFi ≥ v.VAFi − ϵuv` and (2) if `u.VAFi = 0`, `v.VAFi = 0`", where `ϵuv` is the VAF noise error margin.
3. **Sum rule (verbatim, Eq. 5):** a valid lineage tree T "must be a spanning tree of the network that satisfies the following requirement `∀ nodes u∈T : ∀i∈ samples : Σ_{v s.t. (u→v)∈T} v.VAFi ≤ u.VAFi + ϵ`. That is, the sum of the VAF centroids of all the children must not exceed the centroid of the parent."
4. **Why inequality, not equality:** "we use inequality here since our method does not require all the true lineage branches to have been observed."

### Zheng L, Dang H, Niknafs N, et al. (2022) — PICTograph

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/ (DOI 10.1093/bioinformatics/btac367, *Bioinformatics* 38(15):3677–3683)
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed paper, *Bioinformatics*)

**How retrieved:** WebSearch query "CITUP tumor phylogeny CCF sum rule pigeonhole ancestor descendant cancer cell fraction clonal tree" → WebFetch of the PMC full-text article with a prompt extracting the sum condition and lineage-precedence rule.

**Key Extracted Points:**

1. **Sum Condition (CCF form, verbatim):** "the CCF of an ancestral clone must be greater than or equal to the sum of CCFs of its descendants." Trees where "the sum of the CCFs of the descendants of a parent node exceeded the parent node's CCF by more than ε₂ were excluded" (default ε₂ = 0.2).
2. **Lineage Precedence (CCF form, verbatim):** "the CCF of any mutation cannot exceed the CCF of its ancestor." Operationally the descendant cluster must be present in a subset of the samples in which the ancestral cluster is present, and the descendant CCF is at most ε₁ greater than the ancestral CCF in each sample (default ε₁ = 0.1).
3. **Sum Condition generalizes the pigeonhole principle:** the constraint that the summed children CCF cannot exceed the parent CCF is the cell-fraction analogue of the pigeonhole principle applied per node.

---

## Documented Corner Cases and Failure Modes

### From Popic et al. (2015) — LICHeE

1. **Branch with insufficient parent budget:** when three candidate sibling groups have CCFs whose sum exceeds the parent CCF, "no more than one of these groups can be a descendant of the parent group, without violating the VAF phylogenetic constraint"; conflicting branches must be removed/re-placed.
2. **Under-determined private mutations:** private (single-sample) mutation groups are "under-constrained (i.e. multiple tree nodes can serve as ancestors)"; their exact placement is ambiguous, so a deterministic tie-break is required.
3. **Noise margin:** equality of the sum rule is relaxed to `≤ u + ϵ` to tolerate VAF/CCF measurement noise.

### From Zheng et al. (2022) — PICTograph

1. **Presence-pattern constraint:** a descendant cluster must be present only in (a subset of) the samples where the ancestor is present; this is the CCF analogue of LICHeE constraint (1).

---

## Test Datasets

### Dataset: Linear (chain) clonal evolution — derived from the lineage-precedence + sum rules

**Source:** Popic et al. (2015) Eq. 2 (ancestor ≥ descendant) and Eq. 5 (sum rule); single sample.

| Cluster | CCF (sample 1) | Expected placement (deepest valid ancestor) |
|---------|----------------|---------------------------------------------|
| Normal (root) | 1.0 | — |
| A | 1.0 | child of Normal (A.CCF ≤ Normal.CCF) |
| B | 0.6 | child of A (0.6 ≤ 1.0; A budget 1.0 ≥ 0.6) |
| C | 0.3 | child of B (0.3 ≤ 0.6; B budget 0.6 ≥ 0.3) |

Edges: Normal→A, A→B, B→C. Trunk = {A}. Branches = {B, C}.

### Dataset: Branching clonal evolution (two samples) — derived from constraints (1),(2),(3)

**Source:** Popic et al. (2015) constraints (1) presence pattern, (2) ancestor ≥ descendant per sample, (3) sum rule per sample.

| Cluster | CCF sample 1 | CCF sample 2 | Expected placement |
|---------|--------------|--------------|--------------------|
| Normal (root) | 1.0 | 1.0 | — |
| A (trunk) | 1.0 | 1.0 | child of Normal |
| B | 0.6 | 0.0 | child of A (B private to s1; not a descendant of C because C.s1=0 < B.s1=0.6) |
| C | 0.0 | 0.7 | child of A (C private to s2; not a descendant of B because B.s2=0 < C.s2=0.7) |

Edges: Normal→A, A→B, A→C. A is the trunk (parent of two sibling branches B, C). Sum rule under A holds per sample: s1 = 0.6+0.0 = 0.6 ≤ 1.0; s2 = 0.0+0.7 = 0.7 ≤ 1.0. B and C cannot be ancestor/descendant of each other (each is absent in the other's sample → constraint (1)/(2) violated both directions).

### Dataset: Sum-rule violation forces branching instead of nesting

**Source:** Popic et al. (2015) Eq. 5 — children CCF sum may not exceed parent CCF.

| Cluster | CCF (sample 1) | Expected placement |
|---------|----------------|--------------------|
| Normal | 1.0 | — |
| A | 1.0 | child of Normal |
| B | 0.6 | child of A (A budget 1.0 ≥ 0.6) |
| C | 0.6 | child of A — NOT child of B, because B already has nothing but B.CCF=0.6 ≥ C.CCF=0.6 would be allowed by lineage rule; however attaching C under A keeps A's children sum 0.6+0.6=1.2 > 1.0 → violates sum rule |

Expected: with B and C both 0.6 they cannot both be children of the same parent (sum 1.2 > 1.0). Deterministic rule attaches B under A first, then C must nest under B (0.6 ≤ 0.6, B budget 0.6 ≥ 0.6) → chain Normal→A→B→C. This shows the sum rule converting a would-be sibling pair into a chain.

---

## Assumptions

1. **ASSUMPTION: Deterministic tie-break for under-constrained placement** — Popic et al. note private/under-constrained clusters admit multiple valid ancestors. To make the output deterministic this implementation attaches each cluster to the *deepest valid ancestor* (the candidate parent with the smallest total CCF whose per-sample sum-rule budget still admits the child); ties broken by ascending cluster id. This selects the most-recent common ancestor consistent with all cited constraints; it does not change which trees are *valid*, only which single valid tree is returned.
2. **ASSUMPTION: Noise margin ε = 0** — the cited sources relax the inequalities by a configurable ε (LICHeE ϵ; PICTograph ε₁=0.1, ε₂=0.2). Because the unit consumes already-clustered CCF point estimates (clustering, with its noise model, is ONCO-CCF-001), the default comparison uses ε = 0 (strict inequalities), exposed as an optional tolerance parameter so callers can supply the source defaults. Setting ε > 0 only widens admissibility; it never changes a strictly-satisfied relationship.

---

## Recommendations for Test Coverage

1. **MUST Test:** Linear chain from descending single-sample CCFs (Normal→A→B→C). — Evidence: Popic 2015 Eq. 2; Zheng 2022 lineage precedence.
2. **MUST Test:** Branching tree from two private single-sample clusters (A→B, A→C siblings). — Evidence: Popic 2015 constraints (1)(2)(3).
3. **MUST Test:** Sum rule rejects two equal-CCF siblings under one parent and forces a chain. — Evidence: Popic 2015 Eq. 5.
4. **MUST Test:** Trunk identification = clusters on the path from root with CCF ≈ 1 in all samples / the unique root child; branch identification = all non-trunk clusters. — Evidence: Popic 2015 (trunk = common predecessor present across all samples).
5. **MUST Test:** Ancestor CCF ≥ descendant CCF holds on every reconstructed edge (invariant). — Evidence: Popic 2015 Eq. 2.
6. **MUST Test:** Per-node sum rule holds on the reconstructed tree (invariant). — Evidence: Popic 2015 Eq. 5; Zheng 2022 sum condition.
7. **SHOULD Test:** Empty input → tree with only the root, no trunk/branch mutations. — Rationale: boundary.
8. **SHOULD Test:** Single cluster → child of root, is the trunk. — Rationale: boundary.
9. **SHOULD Test:** Null input / CCF out of [0,1] / NaN / inconsistent sample count → exceptions. — Rationale: documented validation.
10. **COULD Test:** Tolerance ε > 0 admits a near-violating edge that ε = 0 would reject. — Rationale: optional parameter.

---

## References

1. Popic V, Salari R, Hajirasouliha I, Kashef-Haghighi D, West RB, Batzoglou S. (2015). Fast and scalable inference of multi-sample cancer lineages. *Genome Biology* 16:91. https://doi.org/10.1186/s13059-015-0647-8 (full-text PDF: https://arxiv.org/pdf/1412.8574)
2. Zheng L, Dang H, Niknafs N, et al. (2022). Estimation of cancer cell fractions and clone trees from multi-region sequencing of tumors (PICTograph). *Bioinformatics* 38(15):3677–3683. https://doi.org/10.1093/bioinformatics/btac367 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC9344857/)

---

## Change History

- **2026-06-15**: Initial documentation.
