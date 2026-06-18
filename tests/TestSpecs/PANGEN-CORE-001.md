# Test Specification: PANGEN-CORE-001

**Test Unit ID:** PANGEN-CORE-001
**Area:** PanGenome
**Algorithm:** Core / Accessory / Unique genome construction, genome fluidity, open/closed classification
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Tettelin et al. (2005), PNAS 102:13950 | 1 | https://doi.org/10.1073/pnas.0506758102 | 2026-06-13 |
| 2 | Tettelin et al. (2008), Curr Opin Microbiol 11:472 | 1 | https://doi.org/10.1016/j.mib.2008.09.006 | 2026-06-13 |
| 3 | Kislyuk et al. (2011), BMC Genomics 12:32 | 1 | https://doi.org/10.1186/1471-2164-12-32 | 2026-06-13 |
| 4 | Page et al. (2015), Roary, Bioinformatics 31:3691 | 3 | https://doi.org/10.1093/bioinformatics/btv421 | 2026-06-13 |
| 5 | micropan `heaps()`/`fluidity()` (CRAN) | 3 | https://rdrr.io/cran/micropan/man/heaps.html | 2026-06-13 |
| 6 | Wikipedia, Pan-genome (citing primaries) | 4 | https://en.wikipedia.org/wiki/Pan-genome | 2026-06-13 |

### 1.2 Key Evidence Points

1. Core genome = gene families present in all genomes; accessory/dispensable = present in some but not all; unique/strain-specific = in exactly one genome — Tettelin (2005, 2008).
2. Operational core threshold is a fraction of genomes (Roary default 99%); membership = present in ≥ coreFraction of genomes, i.e. occupancy / N ≥ coreFraction (a fractional/percentage test, NOT floor(coreFraction · N)) — Page et al. (2015): "a gene being in at least 99% of samples".
3. Genome fluidity `φ = [2/(N(N−1))]·Σ_{k<l}(U_k+U_l)/(M_k+M_l)`, range 0..1; 0 = identical gene content, 1 = disjoint — Kislyuk (2011).
4. Open pan-genome ⟺ Heaps'-law decay exponent alpha < 1 (new genes keep accumulating); closed ⟺ alpha > 1 — Tettelin (2008), micropan.

### 1.3 Documented Corner Cases

- Empty input → empty pan-genome.
- Single genome (N=1): no pairs ⇒ fluidity 0; every cluster has occupancy 1 (unique).
- Fluidity pair with `M_k+M_l = 0` contributes 0 (undefined term, neutral element).
- Openness fit needs ≥ 3 genomes; below that openness is not determinable from the new-gene curve (defaults to Closed).

### 1.4 Known Failure Modes / Pitfalls

1. Using an ad-hoc "unique fraction > 0.1" heuristic for open/closed instead of the Heaps'-law decay exponent — not source-backed; corrected in this unit. — Tettelin (2008), micropan.
2. Confusing the per-pair symmetric-difference fluidity with a global symmetric difference; must be averaged over pairs with the `2/(N(N−1))` factor. — Kislyuk (2011).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ConstructPanGenome(genomes, identityThreshold, coreFraction)` | PanGenomeAnalyzer | Canonical | Partitions clusters into core/accessory/unique; computes fluidity + openness |
| `ClusterGenes(genomes, identityThreshold)` | PanGenomeAnalyzer | Internal | Occupancy source; exercised via ConstructPanGenome and directly for occupancy invariants |
| `GetCoreGeneClusters(clusters, totalGenomes, threshold)` | PanGenomeAnalyzer | Canonical | Core-gene identification (`IdentifyCoreGenes` referent in Registry) |
| `CalculateGenomeFluidity` (private, via ConstructPanGenome.Statistics.GenomeFluidity) | PanGenomeAnalyzer | Internal | Kislyuk formula |
| `DeterminePanGenomeType` (private, via ConstructPanGenome.Statistics.Type) | PanGenomeAnalyzer | Internal | Heaps decay-exponent openness |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | core ∪ accessory ∪ unique = all clusters; the three sets are disjoint | Yes | Tettelin (2005, 2008) |
| INV-02 | CoreGeneCount + AccessoryGeneCount + UniqueGeneCount = TotalGenes | Yes | Tettelin (2008) |
| INV-03 | A cluster is core ⟺ occupancy / N ≥ coreFraction (present in ≥ coreFraction of genomes); unique ⟺ occupancy = 1; accessory otherwise | Yes | Page et al. (2015) |
| INV-04 | 0 ≤ GenomeFluidity ≤ 1 | Yes | Kislyuk (2011) |
| INV-05 | Identical gene content across all genomes ⇒ fluidity = 0; pairwise-disjoint ⇒ fluidity = 1 | Yes | Kislyuk (2011) |
| INV-06 | CoreFraction = CoreGeneCount / TotalGenes (0 when TotalGenes = 0) | Yes | definitional |
| INV-07 | Open ⟺ Heaps decay exponent alpha < 1 (requires N ≥ 3; else Closed) | Yes | Tettelin (2008), micropan |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Partition core/accessory/unique | 3 genomes: c1 all, c2 two, c3/c4/c5 one each; coreFraction 1.0 | Core={c1}, Accessory has the 2-occupancy cluster, Unique count = 3 | Tettelin (2005,2008); Page (2015) |
| M2 | INV-02 counts sum to TotalGenes | same input | Core+Accessory+Unique = TotalGenes | Tettelin (2008) |
| M3 | Fluidity exact value | hand-derived 3-genome example (Evidence dataset) | GenomeFluidity = 0.5555555555… (10/18) within 1e-10 | Kislyuk (2011) |
| M4 | Fluidity = 0 for identical gene content | all genomes share the same single cluster | GenomeFluidity = 0 | Kislyuk (2011) |
| M5 | Fluidity = 1 for pairwise-disjoint content | each genome has a distinct cluster only | GenomeFluidity = 1 | Kislyuk (2011) |
| M6 | INV-04 bounds | M3 input | 0 ≤ fluidity ≤ 1 | Kislyuk (2011) |
| M7 | Open classification | ≥3 genomes, each adding new unique genes (decay alpha < 1) | Type = Open | Tettelin (2008); micropan |
| M8 | Closed classification | ≥3 genomes with shared core, few/no new genes after first (decay alpha > 1) | Type = Closed | Tettelin (2008); micropan |
| M9 | Core-gene identification threshold | GetCoreGeneClusters with threshold 1.0 over occupancy {3,2,1}/3 | only the occupancy-3 cluster returned | Page (2015) |
| M10 | INV-06 core fraction | M1 input | CoreFraction = CoreGeneCount / TotalGenes | definitional |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty input | empty genomes dict | empty result, TotalGenomes=0, fluidity 0, Type Closed | corner case |
| S2 | Null input | null genomes | empty result (no throw) | matches existing contract |
| S3 | Single genome | N=1, two clusters | no pairs ⇒ fluidity 0; both clusters unique (occupancy 1) | corner case |
| S4 | Core threshold boundary | 3 genomes, coreFraction 0.99 → core iff occupancy/3 ≥ 0.99 → only occupancy 3; occupancy 2 (66.7%) is accessory | only the 3/3 cluster counted core; 2/3 accessory | Page (2015) fractional "≥ 99% of samples" |
| S4b | Core threshold float boundary | N=100, coreFraction 0.99 → 99/100 (99%) core, 98/100 (98%) not | 99-occupancy cluster core, 98 not | Page (2015) exact 99% boundary |
| S5 | GetCoreGeneClusters empty | empty cluster list | empty | trivial |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | ConstructPanGenome called twice on same input | identical Statistics | reproducibility |
| C2 | Property: fluidity in [0,1] on random-ish structured inputs | several structured inputs | always in [0,1] | INV-04 property (O(N²) algo) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzerTests.cs` — pre-existing fixture covering ConstructPanGenome, ClusterGenes, presence/absence, Heaps, core/accessory/markers. Assertions are predominantly permissive (`GreaterThan`, `GreaterThanOrEqualTo`, `Or.EqualTo`), and the open/closed test asserts "Open OR Closed" (no exact value). No exact fluidity value, no exact partition counts, no source-backed openness check.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 partition exact | ❌ Missing | existing only checks `CoreGeneCount > 0` |
| M2 counts sum | ❌ Missing | not present |
| M3 fluidity exact | ❌ Missing | existing only checks `>= 0` (S2-style weak) |
| M4 fluidity 0 | ❌ Missing | not present |
| M5 fluidity 1 | ❌ Missing | not present |
| M6 fluidity bounds | ⚠ Weak | `ConstructPanGenome_CalculatesGenomeFluidity` uses `GreaterThanOrEqualTo(0)` only |
| M7 open | ⚠ Weak | `ConstructPanGenome_DeterminesPanGenomeType` allows Open OR Closed |
| M8 closed | ❌ Missing | not present |
| M9 core-gene identification | ⚠ Weak | `GetCoreGeneClusters_FiltersByThreshold` exists with exact value — but in old fixture; re-author in canonical file |
| M10 core fraction | ⚠ Weak | `ConstructPanGenome_CalculatesCoreFraction` only checks 0..1 range |
| S1 empty | ✅ Covered | `ConstructPanGenome_EmptyInput_ReturnsEmptyResult` (kept logic, re-authored exact) |
| S2 null | ❌ Missing | not present |
| S3 single genome | ⚠ Weak | `ConstructPanGenome_SingleGenome_AllUnique` only checks `TotalGenes > 0` |
| S4 core boundary | ❌ Missing | not present |
| S5 GetCoreGeneClusters empty | ❌ Missing | not present |
| C1 determinism | ❌ Missing | not present |
| C2 fluidity property | ❌ Missing | not present |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_ConstructPanGenome_Tests.cs` — all PANGEN-CORE-001 cases (partition, fluidity, openness, core-gene identification) with exact evidence-based values.
- **Remove:** from the legacy `PanGenomeAnalyzerTests.cs`, the ConstructPanGenome/GetCoreGeneClusters/fluidity/type tests that this unit now owns canonically (the weak/permissive ones), to avoid duplicate ownership. Tests for out-of-scope methods (Heaps fit, presence/absence matrix, accessory analysis, markers, core alignment) remain in the legacy file (other units).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PanGenomeAnalyzer_ConstructPanGenome_Tests.cs` | Canonical PANGEN-CORE-001 | 17 |
| `PanGenomeAnalyzerTests.cs` | Legacy, out-of-scope methods only | reduced |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact partition test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented sum invariant | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented exact fluidity 10/18 | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented fluidity=0 | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented fluidity=1 | ✅ Done |
| 6 | M6 | ⚠ Weak | Re-authored exact bounds | ✅ Done |
| 7 | M7 | ⚠ Weak | Re-authored exact Open | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented exact Closed | ✅ Done |
| 9 | M9 | ⚠ Weak | Re-authored in canonical file | ✅ Done |
| 10 | M10 | ⚠ Weak | Re-authored exact fraction | ✅ Done |
| 11 | S1 | ✅ Covered | Re-authored exact | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented null | ✅ Done |
| 13 | S3 | ⚠ Weak | Re-authored exact | ✅ Done |
| 14 | S4 | ❌ Missing | Implemented boundary | ✅ Done |
| 15 | S5 | ❌ Missing | Implemented empty | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented determinism | ✅ Done |
| 17 | C2 | ❌ Missing | Implemented property test | ✅ Done |

**Total items:** 17
**✅ Done:** 17 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact partition asserted |
| M2 | ✅ | sum invariant asserted |
| M3 | ✅ | fluidity = 10/18 within 1e-10 |
| M4 | ✅ | fluidity = 0 |
| M5 | ✅ | fluidity = 1 |
| M6 | ✅ | 0 ≤ fluidity ≤ 1 |
| M7 | ✅ | Type = Open exact |
| M8 | ✅ | Type = Closed exact |
| M9 | ✅ | only occupancy-3 cluster |
| M10 | ✅ | CoreFraction exact |
| S1 | ✅ | empty result exact |
| S2 | ✅ | null → empty |
| S3 | ✅ | single genome exact |
| S4 | ✅ | floor boundary |
| S5 | ✅ | empty clusters |
| C1 | ✅ | determinism |
| C2 | ✅ | property bounds |

**Total in-scope cases:** 17 | **✅:** 17

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Clustering identity metric (k-mer Jaccard) is upstream of partition logic; tests use unambiguous identical/disjoint sequences | M1, M3–M8 inputs |
| 2 | Empty-pair (`M_k+M_l=0`) fluidity term = 0 | fluidity edge cases |

---

## 7. Open Questions / Decisions

1. `DeterminePanGenomeType` previously used an unsourced `uniqueFraction > 0.1` heuristic. Decision: replaced with the source-backed Heaps'-law decay-exponent criterion (open ⟺ alpha < 1, requires N ≥ 3 else Closed). This is a correctness-affecting fix within `ConstructPanGenome`'s scope.
