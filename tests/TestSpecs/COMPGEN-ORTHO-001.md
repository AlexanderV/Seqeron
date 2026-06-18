# Test Specification: COMPGEN-ORTHO-001

**Test Unit ID:** COMPGEN-ORTHO-001
**Area:** Comparative
**Algorithm:** Ortholog identification (Reciprocal Best Hits) and paralog (in-paralog) identification
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Fitch WM 1970, *Syst. Zool.* 19:99–106 (orthology/paralogy definitions) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3178060/ | 2026-06-13 |
| 2 | Tatusov, Koonin, Lipman 1997, *Science* 278:631–637 (COG, symmetrical best hits) | 1 | https://doi.org/10.1126/science.278.5338.631 | 2026-06-13 |
| 3 | Moreno-Hagelsieb & Latimer 2008, *Bioinformatics* 24:319–324 (RBH definition, ≥50% coverage) | 1 | https://doi.org/10.1093/bioinformatics/btm585 | 2026-06-13 |
| 4 | Remm, Storm, Sonnhammer 2001, *J. Mol. Biol.* 314:1041–1052 (RBH seed, in-paralog rule) | 1 | https://doi.org/10.1006/jmbi.2000.5197 | 2026-06-13 |
| 5 | Ondov et al. 2016, *Genome Biol.* 17:132 (alignment-free k-mer similarity) | 1 | https://doi.org/10.1186/s13059-016-0997-x | 2026-06-13 |

### 1.2 Key Evidence Points

1. RBH = two genes in two genomes are orthologs iff each is the other's best hit — Moreno-Hagelsieb & Latimer (2008).
2. Best hit = candidate with maximum similarity score; ties broken deterministically — Moreno-Hagelsieb & Latimer (2008).
3. A one-directional best hit is NOT an ortholog (symmetry required) — Tatusov et al. (1997).
4. Coverage ≥ 50% of a sequence required — Moreno-Hagelsieb & Latimer (2008).
5. Paralogs split by gene duplication within a genome; in-paralogs are within-genome mutual best hits — Fitch (1970), Remm et al. (2001).

### 1.3 Documented Corner Cases

- Tie in best hit broken deterministically (Moreno-Hagelsieb 2008).
- Coverage filter rejects short high-scoring matches (Moreno-Hagelsieb 2008).
- Non-reciprocal best hit excluded (Tatusov 1997).
- Empty/single-gene genome → no pair (definitions require two genes).

### 1.4 Known Failure Modes / Pitfalls

1. Accepting a one-directional best hit as an ortholog (the pre-existing `FindOrthologs` defect) — Tatusov et al. (1997).
2. Spurious orthology from a short shared region without coverage gate — Moreno-Hagelsieb & Latimer (2008).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindOrthologs(genome1Genes, genome2Genes, minIdentity, minCoverage)` | ComparativeGenomics | Canonical | RBH; rewritten to require reciprocity |
| `FindParalogs(genes, minIdentity, minCoverage)` | ComparativeGenomics | Canonical | within-genome mutual best hits (in-paralogs); new method |
| `FindReciprocalBestHits(...)` | ComparativeGenomics | Internal | shares the RBH core; covered indirectly |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned ortholog pair is reciprocal: gene1's best qualifying hit is gene2 AND gene2's best qualifying hit is gene1 | Yes | Moreno-Hagelsieb & Latimer (2008); Tatusov et al. (1997) |
| INV-2 | RBH yields a matching: no genome-1 gene appears in two pairs and no genome-2 gene appears in two pairs | Yes | RBH definition (best hit is unique) — Moreno-Hagelsieb (2008) |
| INV-3 | Every paralog pair is an unordered within-genome pair of distinct genes that are mutual best hits | Yes | Fitch (1970); Remm et al. (2001) |
| INV-4 | Result is deterministic and order-independent for fixed input | Yes | Deterministic best-hit ranking with tie-break (Moreno-Hagelsieb 2008) |
| INV-5 | Pairs below the identity or coverage threshold are excluded | Yes | Moreno-Hagelsieb & Latimer (2008) ≥50% coverage |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | RBH mutual best hits | Two genomes, two identical gene pairs | 2 ortholog pairs {a1↔b1, a2↔b2} | Moreno-Hagelsieb (2008) RBH def |
| M2 | Non-reciprocal excluded | a1 best=b1; b2 also best=a1 but a1's best=b1 | 1 ortholog pair {a1↔b1}; b2 excluded | Tatusov (1997) symmetry |
| M3 | Below threshold excluded | Dissimilar sequences below minIdentity/coverage | 0 ortholog pairs | Moreno-Hagelsieb (2008) ≥50% coverage |
| M4 | Empty genome | genome2 empty | 0 ortholog pairs | def requires a pair |
| M5 | Null inputs (orthologs) | null genome1 / null genome2 | `ArgumentNullException` | repository contract |
| M6 | Paralog mutual best hit | G1 has duplicate p1,p2 + unrelated q1 | 1 paralog pair {p1↔p2} | Fitch (1970); Remm (2001) |
| M7 | Single-gene genome (paralogs) | one gene | 0 paralog pairs | needs two genes |
| M8 | Null inputs (paralogs) | null genes | `ArgumentNullException` | repository contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Matching property (INV-2) | Several genes, multiple RBH | No gene id appears twice on either side | INV-2 |
| S2 | Symmetry of pairs (INV-1) | Each returned pair verified reciprocal by re-checking best hits | All pairs reciprocal | INV-1 |
| S3 | Genes without sequence skipped | Gene with empty Sequence | Not paired | similarity undefined |
| S4 | Paralog identity not self-paired | gene not paired with itself | p1≠p2 in pair | INV-3 |
| S5 | FindReciprocalBestHits == FindOrthologs | public RBH entry point directly | same matching {a1↔b1, a2↔b2} | delegation; Moreno-Hagelsieb (2008) |
| S6 | Null inputs (FindReciprocalBestHits) | null genome1 / null genome2 | `ArgumentNullException` | repository contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | Run twice, compare pair sets | Identical | INV-4 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No `ComparativeGenomics_FindOrthologs_*` or `*_FindParalogs_*` test file exists.
- `FindOrthologs` existed in `ComparativeGenomics.cs` but returned **one-directional** best hits (non-conforming to RBH; violated INV-1). `FindParalogs` did not exist.
- Sibling: `ComparativeGenomics_FindSyntenicBlocks_Tests.cs` (COMPGEN-SYNTENY-001).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new — guards the RBH reciprocity defect |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new (FindParalogs did not exist) |
| M7 | ❌ Missing | new |
| M8 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| C1 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindOrthologs_Tests.cs` — all FindOrthologs and FindParalogs tests for this unit, `#region` per method.
- **Remove:** nothing (no prior tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_FindOrthologs_Tests.cs` | Canonical for COMPGEN-ORTHO-001 | 13 |

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
| 9 | S1 | ❌ Missing | Implemented | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `FindOrthologs_MutualBestHits_ReturnsBothPairs` |
| M2 | ✅ | `FindOrthologs_OneDirectionalBestHit_NotReturnedAsOrtholog` |
| M3 | ✅ | `FindOrthologs_BelowThreshold_ReturnsNoPairs` |
| M4 | ✅ | `FindOrthologs_EmptyGenome_ReturnsEmpty` |
| M5 | ✅ | `FindOrthologs_NullInput_Throws` |
| M6 | ✅ | `FindParalogs_DuplicateGene_ReturnsParalogPair` |
| M7 | ✅ | `FindParalogs_SingleGene_ReturnsEmpty` |
| M8 | ✅ | `FindParalogs_NullInput_Throws` |
| S1 | ✅ | `FindOrthologs_MultipleGenes_YieldsMatchingNoGeneTwice` |
| S2 | ✅ | `FindOrthologs_AllPairs_AreReciprocal` |
| S3 | ✅ | `FindOrthologs_GeneWithoutSequence_IsSkipped` |
| S4 | ✅ | `FindParalogs_PairGenesAreDistinct` |
| S5 | ✅ | `FindReciprocalBestHits_SameAsFindOrthologs` (added 2026-06-15) |
| S6 | ✅ | `FindReciprocalBestHits_NullInput_Throws` (added 2026-06-15) |
| C1 | ✅ | `FindOrthologs_RunTwice_IsDeterministic` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Best-hit ranking uses alignment-free 5-mer Jaccard similarity (no bit-score available in project); reciprocity, coverage gate, and thresholds are source-backed. Order-preserving for the test datasets. | Score metric in FindOrthologs / FindParalogs |

---

## 7. Open Questions / Decisions

1. **Decision:** `FindOrthologs` was rewritten from a one-directional best hit to RBH to conform to Tatusov (1997) and Moreno-Hagelsieb (2008). The pre-existing behavior violated INV-1 and was a defect.
2. **Decision:** `FindParalogs(genes, minIdentity, minCoverage)` added (Registry-listed method that did not exist) implementing within-genome mutual best hits per Fitch (1970) / Remm et al. (2001).
3. **Assumption 1** (score metric) is documented and order-preserving for the evidence datasets; the correctness-critical reciprocity logic is source-backed, so it does not block completion.
