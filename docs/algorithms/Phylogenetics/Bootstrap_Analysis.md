# Phylogenetic Bootstrap Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Phylogenetics |
| Test Unit ID | PHYLO-BOOT-001 |
| Related Projects | Seqeron.Genomics.Phylogenetics |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Phylogenetic bootstrap analysis estimates how strongly the data support each grouping (clade) in an inferred tree. It is a probabilistic resampling procedure: the alignment columns are resampled with replacement many times, a tree is rebuilt on each pseudo-alignment, and the support of a clade is the proportion of replicate trees that recover it [1][2]. It is the standard way to attach confidence values to a phylogeny and should be used whenever a tree's branches need a quantitative support measure. The result is not exact but converges with the number of replicates; for fixed inputs and a fixed RNG seed it is reproducible.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A multiple sequence alignment has rows (taxa) and columns (homologous sites/characters). A distance-based tree (UPGMA or Neighbor-Joining) summarizes the pairwise distances into a branching topology. Each internal branch of a tree separates the taxa into two sets; the set of taxa below an internal node is a *clade* [2]. Bootstrap analysis asks: if we had sampled a slightly different set of characters, would we still recover the same clades?

### 2.2 Core Model

Given an alignment of length L over n taxa and a tree-building function T:

1. Build the reference tree T₀ from the original alignment and record its set of non-trivial clades C(T₀) [2].
2. For each replicate r = 1..B: draw L column indices uniformly with replacement from {0..L−1}, assemble a pseudo-alignment of the **same length L** keeping all taxa, and build a replicate tree Tᵣ [1][3]. Felsenstein (1985): the method "involves resampling points from one's own data, with replacement, to create a series of bootstrap samples of the same size as the original data" and one should "keep all of the original species while sampling characters with replacement" [1].
3. For each clade c ∈ C(T₀), its bootstrap support is

   support(c) = ( #{ r : c ∈ C(Tᵣ) } ) / B

   i.e. the fraction of replicate trees that contain a clade with the identical leaf-name set [2][3]. Lemoine et al.: "measure the support of every branch in the reference tree as the proportion of pseudo-trees containing that branch" [2]. Biopython's `get_support` computes `(count) * 100 / size` over `find_clades(terminal=False)` comparing clades by their terminal bitstring [3].

A group appearing in 100% of replicates has support 1.0 [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Characters (columns) are independent and identically distributed draws | Resampling distribution misestimates error; support values are biased [1] |
| ASM-02 | The reference tree's clades are the entities of interest; replicate-only clades are not scored | Support is undefined for groupings absent from the original-data tree [2] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 0 ≤ support(c) ≤ 1 for every reported clade c | support = count/B with 0 ≤ count ≤ B [2][3] |
| INV-02 | support(c) = k/B for some integer k | count is an integer number of matching replicate trees [3] |
| INV-03 | The reported clades equal exactly the non-trivial clades of the reference tree | only `C(T₀)` is scored [2][3] |
| INV-04 | For fixed (alignment, B, methods, seed), the result is identical across runs | resampling is the only randomness and is seeded |
| INV-05 | A clade recovered in every replicate has support 1.0 | count = B ⇒ B/B = 1 [1] |

### 2.5 Comparison with Related Methods

| Aspect | Felsenstein bootstrap (this) | Transfer bootstrap (TBE) |
|--------|------------------------------|--------------------------|
| Branch match | exact leaf-set identity (binary) [2] | gradual transfer distance [2] |
| Behavior on large trees | support drops with a single misplaced taxon | more robust to rogue taxa |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequences` | `IReadOnlyDictionary<string,string>` | required | Named aligned sequences | non-null, ≥2 entries, equal length |
| `replicates` | `int` | 100 | Number of bootstrap replicates B | ≥ 1 |
| `distanceMethod` | `DistanceMethod` | `JukesCantor` | Pairwise distance model | — |
| `treeMethod` | `TreeMethod` | `UPGMA` | Tree-construction method | UPGMA or NeighborJoining |
| `seed` | `int` | 42 | RNG seed for column resampling | any int; fixed ⇒ reproducible |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IReadOnlyDictionary<string,double>` | Map from clade key (sorted, `\|`-joined leaf names) to support proportion in [0,1] |

### 3.3 Preconditions and Validation

- `sequences == null` → `ArgumentNullException`.
- `sequences.Count < 2` → `ArgumentException` (a tree needs ≥2 taxa).
- `replicates < 1` → `ArgumentException` (the denominator must be ≥1).
- Unequal sequence lengths surface as `ArgumentException` from `BuildTree` (alignment requirement).
- Distance computation uppercases bases and ignores gaps/ambiguous characters (per `CalculatePairwiseDistance`); clade comparison is by leaf-name set, case-sensitive on taxon names.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; build the reference tree T₀ and collect its non-trivial clades.
2. Initialize a support counter to 0 for each reference clade.
3. For B replicates: resample L column indices with replacement (seeded RNG), build the pseudo-alignment, build the replicate tree, and increment the counter of every reference clade present in it.
4. Return counter/B for each clade.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Clade identity key: the subtree's leaf names, sorted ascending and joined with `|` (see `GetClades`/`CollectClades`). Two clades match iff their key strings are equal — i.e. identical leaf-name sets, independent of branch length or internal labels [3].
- Only non-trivial clades (more than one taxon and fewer than all) are scored, matching Biopython's `find_clades(terminal=False)` [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Bootstrap(B replicates, n taxa, L sites) | O(B · (n·L + n³)) | O(n² + n·L) | Each replicate resamples L columns (n·L) and rebuilds an O(n³) distance tree; clade matching is O(n) per internal node |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PhylogeneticAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Phylogenetics/PhylogeneticAnalyzer.cs)

- `PhylogeneticAnalyzer.Bootstrap(sequences, replicates, distanceMethod, treeMethod, seed)`: runs the resampling, returns clade→support proportions.
- Reuses `BuildTree` (reference and replicate trees) and the internal `GetClades`/`CollectClades` (clade key extraction and non-trivial filtering).

### 5.2 Current Behavior

- The reference tree is built once from the original data; its non-trivial clades define the keys of the result. Replicate trees only contribute counts; clades unique to replicates are never reported.
- Column resampling uses a single `Random(seed)` instance; the seed is an explicit parameter (default 42) so callers can reproduce or vary results.
- **Search reuse:** the repository suffix tree was evaluated and is **not used** — this unit performs RNG resampling, distance-matrix tree building, and leaf-set clade matching; there is no substring/pattern search, so the suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Resampling alignment columns with replacement to a pseudo-alignment of the same length L, keeping all taxa [1][3].
- Support = proportion of replicate trees containing a clade with the identical leaf-set, computed over the reference tree's non-trivial clades [2][3].
- A clade present in every replicate yields support 1.0 [1].

**Intentionally simplified:**

- Clades are compared as **rooted** subtree leaf-sets (UPGMA/NJ produce rooted trees here) rather than as unrooted bipartitions; **consequence:** identical to Biopython's clade-based `get_support`, but a user expecting unrooted-bipartition counting on an explicitly unrooted tree should note the rooting convention [3].
- Support is returned as a **proportion in [0,1]** rather than a percentage; **consequence:** multiply by 100 to obtain the value reported in Felsenstein (1985)/Lemoine (2018) [1][2].

**Not implemented:**

- Majority-rule consensus tree construction from the replicates; **users should rely on:** the per-clade support map plus the reference tree (no current in-repo consensus-tree builder) [1].
- Transfer Bootstrap Expectation (TBE) gradual support; **users should rely on:** external tools (booster) [2].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Rooted-clade scoring vs unrooted bipartitions | Assumption | Affects interpretation only; matches Biopython | accepted | ASM-02; Evidence Assumption 1 |
| 2 | Support as proportion [0,1] vs percentage | Assumption | Units/labeling only | accepted | ×100 = published percentage; Evidence Assumption 2 |
| 3 | `seed` parameter added to a previously seed-hardcoded method | Deviation | Enables deterministic tests; default 42 preserves prior behavior | fixed | Existing callers unaffected (seed defaulted) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `sequences == null` | `ArgumentNullException` | input contract |
| fewer than 2 sequences | `ArgumentException` | cannot build a tree |
| `replicates < 1` | `ArgumentException` | denominator must be ≥1 |
| all-identical sequences | every reported clade has support 1.0 | zero-distance matrix reproduces one topology each replicate (INV-05) |
| two well-separated groups | support 1.0 for each group | distances invariant under column resampling (Evidence Dataset 1) |

### 6.2 Limitations

- Distance-based bootstrap inherits UPGMA's molecular-clock assumption; support reflects the chosen distance/tree method, not maximum-likelihood or parsimony.
- Support is sampling-noise sensitive for small B; the value is only reproducible for a fixed seed.
- No consensus tree is produced; only per-clade support of the reference tree.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var sequences = new Dictionary<string, string>
{
    ["A"] = "AAAAAAAAAA",
    ["B"] = "AAAAAAAAAA",
    ["C"] = "GGGGGGGGGG",
    ["D"] = "GGGGGGGGGG",
};
var support = PhylogeneticAnalyzer.Bootstrap(
    sequences, replicates: 100, seed: 42);
// support["A|B"] == 1.0 and support["C|D"] == 1.0
// (every resampled column set leaves d(A,B)=d(C,D)=0 and A/B vs C/G saturated,
//  so every replicate recovers the same {A,B},{C,D} topology — Felsenstein 1985)
```

### 7.2 Performance Baseline

This is an O(B·n³) algorithm (B replicates, each an O(n³) tree build). Measured baseline (informal, `dotnet test --no-build`, .NET runtime, Apple-silicon dev machine, 2026-06-13): the PHYLO-BOOT-001 fixture of 10 tests — which includes runs of 100, 50, 40, 30, 25, 20 and 10 replicates on 3–4 taxon alignments — completes in ≈10 ms total. The property-based invariant `Bootstrap_KnownReplicateCount_SupportsAreQuantizedToCountOverReplicates` (INV-2) is the O(n²)+ property test required by the Definition of Done.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PhylogeneticAnalyzer_Bootstrap_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PhylogeneticAnalyzer_Bootstrap_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [PHYLO-BOOT-001-Evidence.md](../../../docs/Evidence/PHYLO-BOOT-001-Evidence.md)
- Related algorithms: [Tree_Construction](Tree_Construction.md), [Tree_Comparison](Tree_Comparison.md)

## 8. References

1. Felsenstein, J. 1985. Confidence Limits on Phylogenies: An Approach Using the Bootstrap. Evolution 39(4):783–791. https://doi.org/10.1111/j.1558-5646.1985.tb00420.x (text retrieved via https://www.osti.gov/biblio/6044842).
2. Lemoine, F., et al. 2018. Renewing Felsenstein's phylogenetic bootstrap in the era of big data. Nature 556:452–456. https://pmc.ncbi.nlm.nih.gov/articles/PMC6030568/
3. Biopython contributors. Bio.Phylo.Consensus (`bootstrap`, `bootstrap_trees`, `get_support`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/Phylo/Consensus.py (accessed 2026-06-13).
