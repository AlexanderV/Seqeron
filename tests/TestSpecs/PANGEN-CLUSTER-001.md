# Test Specification: PANGEN-CLUSTER-001

**Test Unit ID:** PANGEN-CLUSTER-001
**Area:** PanGenome
**Algorithm:** Gene Clustering (homolog grouping by global sequence identity, CD-HIT greedy incremental model)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Li W, Godzik A (2006) Cd-hit. Bioinformatics 22(13):1658 | 3 | https://doi.org/10.1093/bioinformatics/btl158 | 2026-06-13 |
| 2 | CD-HIT User's Guide (`-c`/`-G` definitions) | 3 | https://vcru.wisc.edu/simonlab/bioinformatics/programs/cd-hit/cdhit-user-guide.pdf | 2026-06-13 |
| 3 | CD-HIT Algorithm wiki (greedy incremental) | 3 | https://github.com/weizhongli/cdhit/wiki/1.-Algorithm | 2026-06-13 |
| 4 | Page AJ et al. (2015) Roary. Bioinformatics 31(22):3691 | 1 | https://doi.org/10.1093/bioinformatics/btv421 | 2026-06-13 |
| 5 | EMBOSS needle manual (percent-identity convention) | 3 | https://galaxy-iuc.github.io/emboss-5.0-docs/needle.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Global sequence identity = number of identical residues in alignment / full length of the **shorter** sequence (CD-HIT default `-G 1`) — Source 2.
2. Identity cutoff `-c` default is 0.9; a sequence is grouped when identity meets the cutoff — Source 2.
3. Greedy incremental clustering: sort long→short; longest becomes first representative; each query compared to existing representatives and grouped into the first one it is similar to, else becomes a new representative — Source 3.
4. Each cluster has one representative sequence; member identity is measured against the representative (not all-pairs) — Source 2.
5. Pan-genome ortholog grouping is identity-based (Roary BLASTP default 95%), confirming sequence identity — not k-mer similarity — is the clustering basis — Source 4.

### 1.3 Documented Corner Cases

- Representative is the longest member (long→short sort) — Source 3.
- Greedy first-match assignment; deterministic for fixed input — Source 3.
- Threshold inclusive (`>=`); `idThreshold = 1.0` clusters only exact-identity-over-shorter-length sequences — Source 2.

### 1.4 Known Failure Modes / Pitfalls

1. Using k-mer/Jaccard similarity in place of true sequence identity misclassifies near-identical genes (the prior implementation's defect) — corrected per Source 2.
2. Internal indels: ungapped identity underestimates CD-HIT's gapped identity — documented ASSUMPTION.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClusterGenes(genomes, identityThreshold)` | PanGenomeAnalyzer | Canonical | Greedy CD-HIT clustering by global identity |
| `CalculateSequenceIdentity(seq1, seq2)` | PanGenomeAnalyzer | Internal | Tested indirectly via `ClusterGenes` (private; identity values asserted through AverageIdentity and cluster membership) |
| `CreatePresenceAbsenceMatrix(genomes, clusters)` | PanGenomeAnalyzer | Delegate | Matrix derived from clusters; smoke verification only |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every input gene appears in exactly one cluster (partition) | Yes | Source 3 (each sequence is representative or redundant) |
| INV-2 | Sum of cluster sizes = total input genes | Yes | Source 3 |
| INV-3 | 0 ≤ global identity ≤ 1; identical sequences = 1.0, disjoint = 0.0 | Yes | Source 2 |
| INV-4 | Sequences with pairwise identity ≥ threshold to a representative are grouped together | Yes | Source 3 |
| INV-5 | GenomeCount = number of distinct genomes contributing a member | Yes | implementation contract / Source 4 |
| INV-6 | Singleton cluster AverageIdentity = 1.0 | Yes | self-identity (Source 2) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Identical sequences cluster | Same sequence in 3 genomes, threshold 0.9 | 1 cluster, GenomeCount 3, AverageIdentity 1.0 | Source 2,3 |
| M2 | Disjoint sequences separate | 3 fully different sequences | 3 singleton clusters | Source 2,3 |
| M3 | Substitution below threshold separates | 2 seqs at 0.875 identity, threshold 0.9 | 2 clusters | Source 2 |
| M4 | Substitution above threshold groups | 2 seqs at 0.9 identity, threshold 0.9 | 1 cluster (inclusive cutoff) | Source 2 |
| M5 | Length-difference identity over shorter length | seq vs same-prefix-longer seq, threshold 1.0 | identity 1.0 → grouped; representative is the longer | Source 2,3 |
| M6 | Greedy multi-cluster | hand-derived dataset (R,Q1,Q2,Q3), threshold 0.8 | exactly 2 clusters: {Q2,R,Q1},{Q3} | Source 3 (Evidence dataset) |
| M7 | Threshold lowered merges | same two near-identical seqs at 0.7 vs 0.9 threshold | 0.9→2 clusters, 0.7→1 cluster | Source 2 |
| M8 | AverageIdentity within cluster | cluster of two seqs at 0.875 (threshold 0.8) | AverageIdentity = 0.875 | Source 2 |
| M9 | Partition invariant | INV-1/INV-2 over mixed input | sum of GeneIds counts = total input genes | Source 3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty genomes | no genomes | empty cluster sequence | input validation |
| S2 | Null genomes | null argument | empty (no throw) | matches ConstructPanGenome null contract |
| S3 | Genome with no genes | empty inner list | no clusters | validation |
| S4 | Singleton AverageIdentity | one gene total | 1 cluster, AverageIdentity 1.0 | INV-6 |
| S5 | CreatePresenceAbsenceMatrix delegation | clusters → matrix | one row per genome; PresentGenes correct | delegate smoke |
| S6 | Empty-sequence identity branches | empty+empty vs empty+non-empty | empties cluster (1.0); non-empty separate (0.0) | identity edge contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | call twice | identical cluster count and membership | greedy long→short stable order |
| C2 | Property (INV-3) | varied sequence pairs via clusters | all AverageIdentity in [0,1] | O(g²·s) property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior dedicated test file for `ClusterGenes`. `ClusterGenes` was exercised only indirectly via `PanGenomeAnalyzer_ConstructPanGenome_Tests.cs` (PANGEN-CORE-001), which uses only identical-or-fully-distinct sequences and never asserts identity values, cluster representative choice, or threshold boundaries.
- New canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ClusterGenes_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new file |
| M2 | ❌ Missing | new file |
| M3 | ❌ Missing | new file |
| M4 | ❌ Missing | new file |
| M5 | ❌ Missing | new file |
| M6 | ❌ Missing | new file |
| M7 | ❌ Missing | new file |
| M8 | ❌ Missing | new file |
| M9 | ❌ Missing | new file |
| S1 | ❌ Missing | new file |
| S2 | ❌ Missing | new file |
| S3 | ❌ Missing | new file |
| S4 | ❌ Missing | new file |
| S5 | ❌ Missing | new file |
| S6 | ❌ Missing | new file |
| C1 | ❌ Missing | new file |
| C2 | ❌ Missing | new file |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ClusterGenes_Tests.cs` — all ClusterGenes + identity + matrix-delegate tests for this unit.
- **Remove:** nothing; PANGEN-CORE-001 tests stay (different unit, ConstructPanGenome).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| PanGenomeAnalyzer_ClusterGenes_Tests.cs | Canonical (this unit) | 17 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented | ✅ Done |
| 13 | S4 | ❌ Missing | Implemented | ✅ Done |
| 14 | S5 | ❌ Missing | Implemented | ✅ Done |
| 15 | S6 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |
| 17 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | ClusterGenes_IdenticalSequences_SingleClusterAllGenomes |
| M2 | ✅ Covered | ClusterGenes_DisjointSequences_SeparateClusters |
| M3 | ✅ Covered | ClusterGenes_SubstitutionBelowThreshold_SeparateClusters |
| M4 | ✅ Covered | ClusterGenes_SubstitutionAtThreshold_SingleCluster |
| M5 | ✅ Covered | ClusterGenes_LengthDifferenceFullIdentity_GroupsAndPicksLongestRepresentative |
| M6 | ✅ Covered | ClusterGenes_GreedyHandDerivedDataset_TwoClusters |
| M7 | ✅ Covered | ClusterGenes_LoweredThreshold_MergesNearIdentical |
| M8 | ✅ Covered | ClusterGenes_TwoMembers_AverageIdentityMatchesGlobalIdentity |
| M9 | ✅ Covered | ClusterGenes_MixedInput_ClustersPartitionAllGenes |
| S1 | ✅ Covered | ClusterGenes_EmptyGenomes_ReturnsEmpty |
| S2 | ✅ Covered | ClusterGenes_NullGenomes_ReturnsEmpty |
| S3 | ✅ Covered | ClusterGenes_GenomeWithNoGenes_ReturnsEmpty |
| S4 | ✅ Covered | ClusterGenes_SingleGene_SingletonAverageIdentityOne |
| S5 | ✅ Covered | CreatePresenceAbsenceMatrix_FromClusters_RowPerGenome |
| S6 | ✅ Covered | ClusterGenes_EmptyAndNonEmptySequences_IdentityBranchesHonoured |
| C1 | ✅ Covered | ClusterGenes_CalledTwice_ProducesIdenticalClustering |
| C2 | ✅ Covered | ClusterGenes_VariedPairs_AverageIdentityInUnitInterval |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Ungapped alignment (identity over shared prefix / shorter length; no internal indels) | CalculateSequenceIdentity; M5, M8 |
| 2 | Homolog clustering only (no paralog/synteny splitting) | ClusterGenes scope |

---

## 7. Open Questions / Decisions

1. **Decision:** Canonical conformance target is CD-HIT global identity (shorter-length denominator), not EMBOSS alignment-length identity, because the canonical method is a greedy identity clusterer matching the CD-HIT model. Documented in algorithm doc §2.2 / §5.3.
2. **Decision:** The checklist names `GeneratePresenceAbsenceMatrix`; the existing, sibling-consistent API is `CreatePresenceAbsenceMatrix`. Kept the existing name (API-shape only, non-correctness-affecting) and noted the discrepancy here and in the registry.
