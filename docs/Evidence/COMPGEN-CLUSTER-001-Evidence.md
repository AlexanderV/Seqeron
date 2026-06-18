# Evidence Artifact: COMPGEN-CLUSTER-001

**Test Unit ID:** COMPGEN-CLUSTER-001
**Algorithm:** Conserved Gene Clusters (common intervals of permutations)
**Date Collected:** 2026-06-14

---

## Online Sources

### Bui-Xuan, Habib, Paul — "MinMax-Profiles: A Unifying View of Common Intervals, Nested Common Intervals and Conserved Intervals of K Permutations" (arXiv:1304.5140)

**URL:** https://arxiv.org/pdf/1304.5140
**Accessed:** 2026-06-14 (PDF fetched via WebFetch; saved locally and extracted with `pdftotext`)
**Authority rank:** 1 (peer-reviewed combinatorics/algorithms, builds directly on Heber & Stoye and Uno & Yagiura)

**Key Extracted Points:**

1. **Interval of a permutation:** "The interval [i, j] of Pk, defined only for 1 ≤ i < j ≤ n, is the set of elements located between position i (included) and position j (included) in Pk." (Generalities §2). When `Pk` is the identity `Idn`, the interval is denoted `(i..j) = {i, i+1, …, j}`.
2. **Common interval (Definition 1, attributed to ref [23] = Uno & Yagiura):** "A common interval of P is a set of integers that is an interval of each Pk, k ∈ [K]." So a common interval is a set of elements that is contiguous in *every* permutation simultaneously.
3. **Worked Example 1 (verbatim):** "Let P1 = Id7 and P2 = (7 2 1 3 6 4 5). Then the common intervals of P = {P1, P2} are (1..2), (1..3), (1..6), (1..7), (3..6), (4..5) and (4..6)." (i.e. the element sets {1,2}, {1,2,3}, {1,…,6}, {1,…,7}, {3,4,5,6}, {4,5}, {4,5,6}).
4. **Interval size:** the interval `[i, j]` is "defined only for 1 ≤ i < j ≤ n" — i.e. intervals (and therefore common intervals) in this formulation have size ≥ 2; singletons are excluded by the `i < j` constraint. The whole set `(1..n)` is always a common interval.

### Uno, Yagiura — "Fast Algorithms to Enumerate All Common Intervals of Two Permutations", Algorithmica 26(2):290–309 (2000)

**URL:** https://doi.org/10.1007/s004539910014 (bibliographic record retrieved from dblp https://dblp.uni-trier.de/rec/journals/algorithmica/UnoY00.html via WebFetch; full-text article body was behind a Springer auth redirect and not opened, so only the citation and the abstract facts below are used)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed; the originating paper of the common-interval model — Definition 1 of the MinMax-Profiles paper cites this work)

**Key Extracted Points:**

1. **Definition (from search abstract / dblp record):** "Given two permutations of n elements, a pair of intervals of these permutations consisting of the same set of elements is called a common interval." This is the model implemented here.
2. **Complexity:** the paper gives an O(n²) algorithm (LHP) and an O(n + K) output-sensitive algorithm (RC), where K is the number of common intervals. The repository implements the simple O(n² · K_genomes) check, adequate for the small gene-cluster inputs in scope.

### Schmidt, Stoye — "Quadratic Time Algorithms for Finding Common Intervals in Two and More Sequences" / Didier et al. — "Extending Common Intervals Searching from Permutations to Sequences" (arXiv:1310.4290)

**URL:** https://arxiv.org/pdf/1310.4290
**Accessed:** 2026-06-14 (PDF fetched via WebFetch; saved locally and extracted with `pdftotext`)
**Authority rank:** 1 (peer-reviewed; extends the model from permutations to sequences with duplicates)

**Key Extracted Points:**

1. **Interval over a sequence (Preliminaries):** "An interval of T is any set I of integers from Σ such that there exist i, j with 1 ≤ i ≤ j ≤ n and I = Set(T[i..j])." A *location* of I on T is such a `[i, j]`.
2. **Common interval of sequences (Definition 1, attributed to refs [8,17]):** "A common interval of two sequences T and S over Σ is a set I of integers that is an interval of both T and S." This is the generalisation handling repeated genes; the gene-cluster setting maps each gene to its ortholog-group label and asks for sets of labels contiguous in every genome.
3. **Worked Example 1 (verbatim):** "Let T = 1 2 5 2 1 4 3 1 2 6 5 and S = 5 6 4 2 3 4 1 5. Then {1, 2} is an interval of T … but is not a common interval of T and S. An example of common interval is {1, 2, 3, 4} which has five locations on T … and two locations on S". Confirms that a set is a common interval iff it appears as a contiguous block (some location) in *every* sequence.

### Heber, Stoye — "Finding All Common Intervals of k Permutations", CPM 2001, LNCS 2089:207–218

**URL:** https://doi.org/10.1007/3-540-48194-X_19 (bibliographic record retrieved via WebSearch result on link.springer.com; used only for citation provenance of the k-permutation generalisation)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed; the canonical k-permutation common-interval paper)

**Key Extracted Points:**

1. **k-permutation generalisation:** common intervals extend to a family of k permutations; an algorithm "finds in a family of k permutations of n elements all z common intervals in optimal O(kn + z) time and O(n) additional space." Confirms the multi-genome (k ≥ 2) cluster definition: a set contiguous in *all* k genomes.

---

## Documented Corner Cases and Failure Modes

### From the MinMax-Profiles paper (Definition 1, §2)

1. **Trivial common intervals:** the full set `(1..n)` is always a common interval of any family; per the `i < j` interval constraint, singletons are not counted as (non-trivial) common intervals. A minimum-size threshold filters short clusters.
2. **Fewer than two permutations:** a common interval is defined over a *family* of permutations; with a single permutation every interval is trivially "common", so the conserved-cluster question is meaningful only for K ≥ 2 genomes.

### From the sequence-extension paper (Didier et al., Preliminaries / Example 1)

1. **Repeated labels:** when genes map to ortholog groups, the same group label can occur multiple times in a genome (paralogs/duplications). A set is a common interval iff *some* contiguous window in each genome has exactly that set of labels (any location suffices).

---

## Test Datasets

### Dataset: MinMax-Profiles Example 1 (two-permutation golden vector)

**Source:** Bui-Xuan, Habib, Paul (arXiv:1304.5140), Example 1.

| Parameter | Value |
|-----------|-------|
| Genome 1 (labels in order) | 1 2 3 4 5 6 7 |
| Genome 2 (labels in order) | 7 2 1 3 6 4 5 |
| All non-trivial common intervals (sets) | {1,2}, {1,2,3}, {3,4,5,6}, {4,5}, {4,5,6}, {1,2,3,4,5,6} |
| Whole-set (trivial) common interval | {1,2,3,4,5,6,7} |

(Independently recomputed by brute force over all element subsets of size ≥ 2: result identical to the paper.)

### Dataset: Sequence-with-duplicates common interval

**Source:** Didier et al. (arXiv:1310.4290), Example 1.

| Parameter | Value |
|-----------|-------|
| Sequence T | 1 2 5 2 1 4 3 1 2 6 5 |
| Sequence S | 5 6 4 2 3 4 1 5 |
| {1,2} is a common interval? | No (interval of T, not of S) |
| {1,2,3,4} is a common interval? | Yes (a contiguous window with this label set exists in both) |

---

## Assumptions

1. **ASSUMPTION: maxGap parameter** — The public method retains a `maxGap` parameter (API/MCP backward compatibility). The strict common-interval model (Uno & Yagiura 2000; Heber & Stoye 2001) is **gap-free**: a cluster must occupy a *contiguous* window in every genome. The validated, tested behaviour is therefore the strict model (the cluster's group set equals the set of all groups in some window, no foreign groups inside). `maxGap` does not relax this strict definition in the validated path and the gene-teams gapped extension (Bergeron, Corteel & Raffinot 2002) is **not** implemented because its source was not retrievable in this session. This is API-shape only; the correctness contract is the strict common-interval definition.

---

## Recommendations for Test Coverage

1. **MUST Test:** Two-genome golden vector P1=Id7, P2=(7 2 1 3 6 4 5) — every non-trivial common interval set must be returned and nothing else. — Evidence: MinMax-Profiles Example 1.
2. **MUST Test:** A set contiguous in genome 1 but split in genome 2 is NOT a cluster (e.g. {1,2} in the golden vector is contiguous in P1 but 1 and 2 are non-adjacent in P2). — Evidence: MinMax-Profiles Example 1 (note {1,2} adjacency in P2: positions of 1 and 2 are 3 and 2 → adjacent, so it IS common; use a set that is genuinely split, e.g. {2,3}: positions in P2 are 2 and 4 → not adjacent → not common).
3. **MUST Test:** Repeated ortholog-group labels (paralogs) — a set that is a contiguous window in each genome despite duplicates is reported. — Evidence: Didier et al. Example 1.
4. **MUST Test:** `minClusterSize` filters out clusters smaller than the threshold. — Evidence: trivial intervals / size threshold (MinMax §2).
5. **SHOULD Test:** Fewer than two genomes returns no clusters. — Rationale: common interval is undefined for K < 2 (MinMax §2 family definition).
6. **SHOULD Test:** Identical gene order across all genomes → every window of size ≥ minClusterSize is conserved. — Rationale: identity vs identity, all intervals common.
7. **COULD Test:** Determinism / order-independence of the returned cluster list. — Rationale: set semantics imply a canonical, reproducible result.

---

## References

1. Bui-Xuan B-M, Habib M, Paul C. 2013. MinMax-Profiles: A Unifying View of Common Intervals, Nested Common Intervals and Conserved Intervals of K Permutations. arXiv:1304.5140. https://arxiv.org/abs/1304.5140
2. Uno T, Yagiura M. 2000. Fast Algorithms to Enumerate All Common Intervals of Two Permutations. Algorithmica 26(2):290–309. https://doi.org/10.1007/s004539910014
3. Didier G, Schmidt T, Stoye J, Tsur D. 2013. Extending Common Intervals Searching from Permutations to Sequences. arXiv:1310.4290. https://arxiv.org/abs/1310.4290
4. Heber S, Stoye J. 2001. Finding All Common Intervals of k Permutations. In: Combinatorial Pattern Matching (CPM 2001), LNCS 2089:207–218. https://doi.org/10.1007/3-540-48194-X_19

---

## Change History

- **2026-06-14**: Initial documentation.
