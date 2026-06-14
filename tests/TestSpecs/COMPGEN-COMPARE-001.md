# Test Specification: COMPGEN-COMPARE-001

**Test Unit ID:** COMPGEN-COMPARE-001
**Area:** Comparative
**Algorithm:** Comprehensive two-genome comparison — core/dispensable (conserved vs genome-specific) gene partition and syntenic-gene fraction
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Tettelin et al. (2005), PNAS 102(39):13950–13955 — pan-genome / core / dispensable | 1 | https://doi.org/10.1073/pnas.0506758102 (https://pmc.ncbi.nlm.nih.gov/articles/PMC1216834/) | 2026-06-14 |
| 2 | Moreno-Hagelsieb & Latimer (2008), Bioinformatics 24(3):319–324 — RBH ortholog criterion | 1 | https://doi.org/10.1093/bioinformatics/btm585 | 2026-06-14 |
| 3 | Wang et al. (2012), NAR 40(7):e49 — MCScanX collinearity (syntenic blocks) | 1 | https://doi.org/10.1093/nar/gkr1293 | 2026-06-14 |
| 4 | Synteny overview (fraction of syntenic genes metric) | 4 | https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/synteny ; https://en.wikipedia.org/wiki/Synteny | 2026-06-14 |

### 1.2 Key Evidence Points

1. Core genome = "genes present in all strains"; a gene shared by both genomes is conserved (core) — Tettelin (2005).
2. Dispensable genome = "genes absent from one or more strains and genes that are unique to each strain"; a gene present in only one genome is genome-specific — Tettelin (2005).
3. Shared (conserved) genes are operationalised as reciprocal best hits between the two genomes — Moreno-Hagelsieb & Latimer (2008); a one-directional best hit is not shared.
4. Overall synteny = fraction of genes inside syntenic blocks — "the fraction of syntenic genes is a metric used to measure synteny conservation" (Synteny overview); blocks come from MCScanX chains scoring ≥250 (≥5 collinear anchors) — Wang et al. (2012).

### 1.3 Documented Corner Cases

- All genes shared → both genome-specific counts = 0 (all core). Tettelin (2005).
- No genes shared → conserved = 0; every gene genome-specific. Tettelin (2005).
- Empty genome → no ortholog pairs → conserved = 0. (RBH corner case.)
- Fewer than 5 collinear orthologs → no syntenic block → OverallSynteny = 0 even with conserved orthologs. MCScanX threshold (Wang et al. 2012).

### 1.4 Known Failure Modes / Pitfalls

1. Counting a one-directional best hit as a shared/core gene (would understate genome-specific counts) — Moreno-Hagelsieb & Latimer (2008).
2. Reporting nonzero synteny for too-few collinear anchors (block threshold ignored) — Wang et al. (2012).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CompareGenomes(genome1Genes, genome2Genes, minOrthologIdentity, minSyntenicBlockSize)` | ComparativeGenomics | **Canonical** | Aggregator: produces the conserved/genome-specific partition + OverallSynteny; orthologs/blocks/rearrangements delegated to validated sub-units. |
| `FindReciprocalBestHits` | ComparativeGenomics | **Delegate** | Validated by COMPGEN-RBH-001; here only as the source of conserved pairs. |
| `FindSyntenicBlocks` | ComparativeGenomics | **Delegate** | Validated by COMPGEN-SYNTENY-001; source of OverallSynteny numerator. |
| `DetectRearrangements` | ComparativeGenomics | **Delegate** | Validated by COMPGEN-REARR-001. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `ConservedGenes == Orthologs.Count` (the conserved/core set equals the reciprocal-best-hit pairs). | Yes | Source 1 (core = shared) + Source 2 (shared = RBH) |
| INV-2 | `ConservedGenes + GenomeSpecificGenes1 == |genome1|` and `ConservedGenes + GenomeSpecificGenes2 == |genome2|` (every gene is core or that genome's dispensable gene; an RBH matching maps each gene at most once). | Yes | Source 1 (core ∪ dispensable = all genes) |
| INV-3 | `0 ≤ OverallSynteny ≤ 1` (it is a fraction, clamped to 1). | Yes | Source 4 (fraction metric) |
| INV-4 | Swapping genome1 and genome2 leaves `ConservedGenes` unchanged and swaps `GenomeSpecificGenes1` ↔ `GenomeSpecificGenes2` (the RBH matching is symmetric). | Yes | Source 2 (reciprocal/symmetric) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | One shared + one unique each | g1={a1=S, b1=U1}, g2={c2=S, d2=U2}; S identical in both, U1≠U2 unique | Conserved=1, Specific1=1, Specific2=1, Orthologs.Count=1 | Source 1 (core/dispensable), Source 2 (RBH) |
| M2 | Disjoint content | g1, g2 share no sequence (2 genes each) | Conserved=0, Specific1=2, Specific2=2, Orthologs empty | Source 1 ("unique to each strain") |
| M3 | Identical content, 5 collinear + 1 unique each | g1, g2 each = 5 distinct shared seqs S₀…S₄ in same order + 1 unique gene | Conserved=5, Specific1=1, Specific2=1, OverallSynteny=5/6=0.8333…, 1 syntenic block of 5 genes, 0 rearrangements | Source 1 (core) + Source 4 (fraction of syntenic genes) + Source 3 (block) |
| M4 | Empty genomes | both gene lists empty | Conserved=0, Specific1=0, Specific2=0, OverallSynteny=0, all collections empty | Corner case (no ortholog pairs) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Conserved but sub-block-size collinear | 3 shared collinear orthologs (< 5) + unique genes | Conserved=3, OverallSynteny=0 (no block reported) | Documents Assumption 2 / MCScanX threshold |
| S2 | Symmetry | M3 inputs with genome1/genome2 swapped | Conserved unchanged; Specific1 ↔ Specific2 | INV-4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | All genes shared | g1, g2 = same 2 sequences (no uniques) | Specific1=0, Specific2=0, Conserved=2 | Tettelin: all core, dispensable empty |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomicsTests.cs` contains a `#region CompareGenomes Tests` with three tests: `CompareGenomes_SimilarGenomes_ReturnsComprehensiveResult` (permissive: `GreaterThan`/`GreaterThanOrEqualTo`), `CompareGenomes_EmptyGenomes_HandlesGracefully` (exact zeros), `CompareGenomes_CompletelyDifferent_ReturnsNoConservation` (exact). These are migrated/replaced by the canonical unit file.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (one shared/one unique) | ❌ Missing | Not present (this is the partition case the prior attempt got wrong). |
| M2 (disjoint) | ⚠ Weak | `CompareGenomes_CompletelyDifferent` exists with exact values but no assertion messages and only 2 genes; rewrite as evidence-based. |
| M3 (identical, 5 collinear, synteny) | ❌ Missing | No synteny-fraction assertion anywhere. |
| M4 (empty) | ⚠ Weak | `CompareGenomes_EmptyGenomes` exists but no messages, no OverallSynteny/collection checks; rewrite. |
| S1 (sub-block collinear) | ❌ Missing | — |
| S2 (symmetry) | ❌ Missing | — |
| C1 (all shared) | ❌ Missing | — |
| `CompareGenomes_SimilarGenomes_ReturnsComprehensiveResult` | 🔁 Duplicate | Permissive; superseded by M1/M3. Remove. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CompareGenomes_Tests.cs` — all COMPGEN-COMPARE-001 tests (M1–M4, S1–S2, C1).
- **Remove:** the entire `#region CompareGenomes Tests` from `ComparativeGenomicsTests.cs` (3 tests) — replaced by the canonical file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_CompareGenomes_Tests.cs` | Canonical COMPGEN-COMPARE-001 | 7 |
| `ComparativeGenomicsTests.cs` (`CompareGenomes` region) | Removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented evidence-based test | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewritten from scratch in canonical file | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented (partition + synteny fraction + block + rearr) | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewritten with full collection + synteny checks | ✅ Done |
| 5 | S1 | ❌ Missing | Implemented | ✅ Done |
| 6 | S2 | ❌ Missing | Implemented (symmetry) | ✅ Done |
| 7 | C1 | ❌ Missing | Implemented | ✅ Done |
| 8 | Duplicate `SimilarGenomes` | 🔁 Duplicate | Removed from `ComparativeGenomicsTests.cs` | ✅ Done |

**Total items:** 8
**✅ Done:** 8 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `CompareGenomes_OneSharedOneUnique_PartitionsCoreAndSpecific` |
| M2 | ✅ | `CompareGenomes_DisjointContent_NoCoreAllSpecific` |
| M3 | ✅ | `CompareGenomes_IdenticalCollinearContent_ReportsCoreAndSyntenyFraction` |
| M4 | ✅ | `CompareGenomes_EmptyGenomes_ReturnsEmptyPartition` |
| S1 | ✅ | `CompareGenomes_FewCollinearOrthologs_ConservedButZeroSynteny` |
| S2 | ✅ | `CompareGenomes_SwappedGenomes_SwapsGenomeSpecificCounts` |
| C1 | ✅ | `CompareGenomes_AllGenesShared_NoGenomeSpecificGenes` |
| Duplicate | ✅ | Removed from `ComparativeGenomicsTests.cs` |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Alignment-free 5-mer Jaccard similarity (id ≥0.3, cov ≥0.5) maps Tettelin's 50%/50% gate (inherited from COMPGEN-RBH-001). Does not affect the partition tested here (identical→pass, disjoint→fail). | M1, M2, M3, S1, S2, C1 |
| 2 | OverallSynteny uses the MCScanX block threshold (≥5 collinear anchors / score ≥250). | M3, S1 |

---

## 7. Open Questions / Decisions

1. None. `CompareGenomes` is an aggregator over already-validated sub-units; this unit verifies the core/dispensable partition and the syntenic-gene fraction, all source-derived.
