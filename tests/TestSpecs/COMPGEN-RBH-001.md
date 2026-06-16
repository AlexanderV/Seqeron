# Test Specification: COMPGEN-RBH-001

**Test Unit ID:** COMPGEN-RBH-001
**Area:** Comparative
**Algorithm:** Reciprocal Best Hits (RBH / bidirectional best hits) for ortholog identification
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Moreno-Hagelsieb G, Latimer K (2008). Choosing BLAST options for better detection of orthologs as reciprocal best hits. *Bioinformatics* 24(3):319–324 | 1 | https://doi.org/10.1093/bioinformatics/btm585 | 2026-06-14 |
| 2 | Tatusov RL, Koonin EV, Lipman DJ (1997). A genomic perspective on protein families. *Science* 278:631–637 (method via NCBI Handbook NBK21090) | 1 | https://doi.org/10.1126/science.278.5338.631 ; https://www.ncbi.nlm.nih.gov/books/NBK21090/ | 2026-06-14 |
| 3 | Ondov BD et al. (2016). Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol.* 17:132 | 1 | https://doi.org/10.1186/s13059-016-0997-x | 2026-06-14 |

### 1.2 Key Evidence Points

1. RBH definition: "two genes residing in two different genomes are deemed orthologs if their protein products find each other as the best hit in the opposite genome." — Moreno-Hagelsieb & Latimer (2008).
2. Best hit = maximum bit-score, ties broken by smallest E-value ("highest to lowest bit-score, then, if the bit-scores were identical, from smallest to highest E-values"). — Moreno-Hagelsieb & Latimer (2008).
3. Qualifying gate: "coverage of at least 50% of any of the protein sequences" and "maximum E-value threshold of 1×10−6". — Moreno-Hagelsieb & Latimer (2008).
4. COGs are built from mutually consistent genome-specific best hits ("Detect triangles of mutually consistent, genome-specific best hits"); the pairwise case is the reciprocal best hit. — Tatusov et al. (1997) via NCBI Handbook.

### 1.3 Documented Corner Cases

- Tie in best hit must be resolved deterministically (Moreno-Hagelsieb 2008).
- A short high-scoring local match is rejected unless ≥ 50% coverage (Moreno-Hagelsieb 2008).
- Non-reciprocity: A→B best but B→C≠A best ⇒ A–B is not an ortholog (Tatusov 1997; Moreno-Hagelsieb 2008).
- A gene with no qualifying hit yields no pair; an empty genome yields no pairs.

### 1.4 Known Failure Modes / Pitfalls

1. Reporting a one-directional best hit as an ortholog (missing reciprocity check). — Tatusov (1997).
2. Reporting placeholder hit metrics (hardcoded coverage 1.0, alignment length 0) instead of the actual hit's coverage/identity/length. — defect class corrected in this unit; metrics must reflect the real hit per Moreno-Hagelsieb (2008).
3. Ignoring the coverage gate so spurious short matches become orthologs. — Moreno-Hagelsieb (2008).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindReciprocalBestHits(genome1Genes, genome2Genes, minIdentity, minCoverage)` | ComparativeGenomics | **Canonical** | RBH ortholog identification (dedicated entry point) |
| `FindOrthologs(...)` | ComparativeGenomics | **Delegate** | Delegates to `FindReciprocalBestHits`; deep tests live under COMPGEN-ORTHO-001 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned pair is reciprocal: g1's best hit is g2 AND g2's best hit is g1. | Yes | Moreno-Hagelsieb & Latimer (2008); Tatusov et al. (1997) |
| INV-2 | The result is a matching: no genome-1 gene and no genome-2 gene appears in two pairs. | Yes | Best-hit reciprocity (Moreno-Hagelsieb 2008) |
| INV-3 | Each pair carries the actual hit identity and coverage (identical sequences ⇒ identity = coverage = 1.0), not placeholders. | Yes | Moreno-Hagelsieb (2008) best-hit metrics |
| INV-4 | Deterministic and order-independent for fixed input (deterministic tie-break). | Yes | Moreno-Hagelsieb (2008) tie-break by E-value |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Mutual best hits returned | Two independent identical pairs (a1≡b1, a2≡b2) | 2 pairs, identity 1.0 each, a1↔b1, a2↔b2 | Moreno-Hagelsieb & Latimer (2008) RBH definition |
| M2 | Non-reciprocity excluded | a1≡b1, b2 is a1's superstring (Jaccard 0.667) | 1 pair (a1↔b1); b2 excluded | Tatusov (1997); Moreno-Hagelsieb (2008) |
| M3 | Actual hit metrics | Identical pair returns real coverage/identity/length | Identity = 1.0, Coverage = 1.0, AlignmentLength = 14 | Moreno-Hagelsieb (2008) best-hit metrics |
| M4 | Below-threshold rejected | minIdentity = 1.0 on a Jaccard-0.667 superstring pair | empty (no qualifying hit) | Moreno-Hagelsieb (2008) ≥50% coverage / significance gate |
| M5 | Empty genome | genome2 empty | empty | RBH requires a between-genome pair |
| M6 | Null inputs throw | null genome1 / null genome2 | `ArgumentNullException` | Repository contract (sibling methods) |
| M7 | Coverage gate (independent of identity gate) | id 0.375 (≥0.3) but coverage 6/11=0.5455; minCoverage=0.6 | kept by default; empty at minCoverage 0.6 | Moreno-Hagelsieb & Latimer (2008) ≥50% coverage gate |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Matching property | Three independent identical pairs | each gene id distinct across pairs | INV-2 |
| S2 | Reciprocity property | Pairs from forward call are reciprocal under reversed call | every pair symmetric | INV-1 |
| S3 | Gene without sequence skipped | a2 has null sequence | a2 never paired | similarity undefined without sequence |
| S4 | Sequence shorter than k=5 | identical 4-nt sequences | empty (similarity 0) | empty k-mer set ⇒ never qualifies |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | run twice / reversed input order | identical pair set | INV-4 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomicsTests.cs` contained a `#region FindReciprocalBestHits Tests` with two weak tests (`...MutualBestMatches_ReturnsRBH`, `...NoMutualBest_ReturnsEmpty`).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| MutualBestMatches_ReturnsRBH | ⚠ Weak | `Has.Count`/`Any` with no assertion messages; no identity/coverage/length check |
| NoMutualBest_ReturnsEmpty | ⚠ Weak | shape-only; no message; not tied to a documented gate value |
| M1 mutual best hits | ❌ Missing | not covered with exact values + messages |
| M2 non-reciprocity | ❌ Missing | the defect class was untested |
| M3 actual hit metrics | ❌ Missing | the placeholder-coverage defect was untested |
| M4 below-threshold | ❌ Missing | gate not exercised by an exact threshold |
| M5 empty genome | ❌ Missing | |
| M6 null throws | ❌ Missing | old impl did not even null-check |
| S1 matching | ❌ Missing | |
| S2 reciprocity | ❌ Missing | |
| S3 no-sequence skipped | ❌ Missing | |
| C1 determinism | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindReciprocalBestHits_Tests.cs` — all RBH evidence-based tests for this unit.
- **Remove:** the `#region FindReciprocalBestHits Tests` block from `ComparativeGenomicsTests.cs` (two weak tests); replace with a pointer comment.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_FindReciprocalBestHits_Tests.cs` | Canonical (this unit) | 11 |
| `ComparativeGenomicsTests.cs` | weak RBH region removed | 0 (for RBH) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | MutualBestMatches_ReturnsRBH | ⚠ Weak | removed; rewritten as M1 with exact values + messages | ✅ Done |
| 2 | NoMutualBest_ReturnsEmpty | ⚠ Weak | removed; subsumed by M4/M5 with exact gate | ✅ Done |
| 3 | M1 | ❌ Missing | implemented | ✅ Done |
| 4 | M2 | ❌ Missing | implemented | ✅ Done |
| 5 | M3 | ❌ Missing | implemented | ✅ Done |
| 6 | M4 | ❌ Missing | implemented | ✅ Done |
| 7 | M5 | ❌ Missing | implemented | ✅ Done |
| 8 | M6 | ❌ Missing | implemented | ✅ Done |
| 9 | S1 | ❌ Missing | implemented | ✅ Done |
| 10 | S2 | ❌ Missing | implemented | ✅ Done |
| 11 | S3 | ❌ Missing | implemented | ✅ Done |
| 12 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 mutual best hits | ✅ Covered | `...MutualBestHits_ReturnsBothPairs` |
| M2 non-reciprocity | ✅ Covered | `...OneDirectionalBestHit_NotReturned` |
| M3 actual hit metrics | ✅ Covered | `...IdenticalPair_ReportsActualCoverageAndLength` |
| M4 below-threshold | ✅ Covered | `...AboveMinIdentity_RejectsNonIdenticalPair` |
| M5 empty genome | ✅ Covered | `...EmptyGenome_ReturnsEmpty` |
| M6 null throws | ✅ Covered | `...NullInput_Throws` |
| S1 matching | ✅ Covered | `...MultipleGenes_YieldsMatching` |
| S2 reciprocity | ✅ Covered | `...AllPairs_AreReciprocal` |
| S3 no-sequence skipped | ✅ Covered | `...GeneWithoutSequence_IsSkipped` |
| C1 determinism | ✅ Covered | `...ReversedInputOrder_IsDeterministic` (+ run-twice) |

Count of ✅ = 10 distinct in-scope cases (M1–M6, S1–S3, C1); all covered.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Alignment-free 5-mer Jaccard similarity is the best-hit ranking metric (no alignment bit-score available in this project). Reciprocity rule, deterministic tie-break, coverage gate, and minimum-similarity gate are source-backed; metric only affects near-identical tie ordering. | M1–M4, ranking |

---

## 7. Open Questions / Decisions

1. **Decision:** `FindReciprocalBestHits` and `FindOrthologs` implement the identical RBH criterion. To prevent divergence, `FindOrthologs` now delegates to `FindReciprocalBestHits`. Deep ortholog/paralog tests stay under COMPGEN-ORTHO-001; this unit owns the dedicated RBH entry point.
2. **Decision:** the prior `FindReciprocalBestHits` was nonconforming (no coverage gate, no deterministic tie-break, placeholder `Coverage=1.0`/`AlignmentLength=0`, `Identity` set to the score product). Corrected to the canonical RBH (Phase 5).
