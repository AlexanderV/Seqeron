# Test Specification: META-FUNC-001

**Test Unit ID:** META-FUNC-001
**Area:** Metagenomics
**Algorithm:** Functional Prediction (homology-based annotation transfer + pathway over-representation)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Altschul et al., NCBI BLAST tutorial — Statistics of Sequence Similarity Scores | 2 | https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html | 2026-06-13 |
| 2 | NCBI C BLAST Toolkit `blast_stat.c` (`blosum62_values` ungapped λ, K, H) | 3 | https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/algo/blast/core/blast_stat.c | 2026-06-13 |
| 3 | NCBI BLOSUM62 substitution matrix file | 2 | https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62 | 2026-06-13 |
| 4 | PNNL Proteomics Data Analysis §8.2 Over-Representation Analysis (hypergeometric) | 4 | https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. E-value of an HSP: `E = K·m·n·e^(−λS)` — source 1.
2. Bit score: `S' = (λS − ln K)/ln 2`; and equivalently `E = m·n·2^(−S')` — source 1.
3. Ungapped BLOSUM62 parameters: λ = 0.3176, K = 0.134, H = 0.4012 — source 2.
4. BLOSUM62 diagonal (self-match) scores (e.g. W = 11, C = 9, A = 4) — source 3.
5. Hypergeometric over-representation p-value: `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)` with N,M,n,x as background/annotated/query/overlap counts — source 4.
6. ORA worked example: N=8000, M=400, n=100, x=20 → P = 7.88 × 10⁻⁸ — source 4.

### 1.3 Documented Corner Cases

- ORA: x = 0 (no overlap) ⇒ p-value = 1 (empty sum). Degenerate margins M = 0 or n = 0 ⇒ p-value = 1 — source 4.
- BLAST: Karlin-Altschul EVD requires a positive raw score; an exact BLOSUM62 self-match of a non-empty segment always scores S > 0 — source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Treating the lower-tail probability as enrichment — ORA reports only the upper tail P(X ≥ x) — source 4.
2. Using gapped λ/K with an ungapped score (or vice-versa) — λ/K are scoring-system specific; the ungapped row of `blast_stat.c` is used here — source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictFunctions(genes, database)` | MetagenomicsAnalyzer | Canonical | Homology transfer; computes BLOSUM62 raw score → bit score → E-value; best hit per gene. |
| `FindPathwayEnrichment(functions, pathwayDb)` | MetagenomicsAnalyzer | Canonical | Hypergeometric over-representation p-value per pathway. |
| `FunctionalBitScore(rawScore, m, n)` (exposed helper) | MetagenomicsAnalyzer | Internal | Bit-score/E-value computation; tested indirectly + directly. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Bit score `S' = (λS − ln K)/ln 2` is strictly increasing in raw score S | Yes | Source 1 (formula is linear in S, λ > 0) |
| INV-2 | `K·m·n·e^(−λS) = m·n·2^(−S')` for the cited λ, K | Yes | Source 1 (algebraic identity) |
| INV-3 | E-value is strictly decreasing in raw score S | Yes | Source 1 (E ∝ e^(−λS), λ > 0) |
| INV-4 | ORA p-value ∈ [0, 1] for all inputs | Yes | Source 4 (probability) |
| INV-5 | ORA p-value = 1 when overlap x = 0, or M = 0, or n = 0 | Yes | Source 4 (empty sum / degenerate margins) |
| INV-6 | The annotation transferred to a gene is the matching database entry with the lowest E-value (best hit) | Yes | Source 1 (E-value ranks significance) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Bit score of exact match | Gene contains signature "WWW" (S = 3×11 = 33) | BitScore = 18.0202932787533 (Within 1e-9) | Source 1 + 2 + 3 |
| M2 | E-value of exact match | Same "WWW" hit, m = n = 3 | EValue = 3.3852730346546e-5 (Within 1e-15) | Source 1 |
| M3 | Best-hit selection | Gene matches two entries (signatures "AAAA" S=16 and "WW" S=22); transfer entry with lower E-value | Annotation = the "WW" entry (lower E-value / higher bit score) | Source 1 (INV-6) |
| M4 | ORA worked example | N=8000, M=400, n=100, x=20 | p-value = 7.88e-8 (Within 1e-9) | Source 4 |
| M5 | ORA no-overlap | Pathway with x = 0 interesting genes | p-value = 1.0 | Source 4 (INV-5) |
| M6 | Annotation fields | Hit transfers Function/Pathway/KO of the matched database entry | Returned FunctionalAnnotation carries the DB entry's Function, Pathway, KoNumber | Source 1 (homology transfer) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty protein sequence | Gene with empty/whitespace sequence | No annotation emitted | Robustness |
| S2 | No matching signature | Gene matches no DB entry | No annotation emitted | Homology transfer requires a hit |
| S3 | Empty / null database | Empty or null function database | No annotations, no throw | Robustness |
| S4 | Null gene list | Null `genes` argument | ArgumentNullException | Mirrors sibling validation |
| S5 | ORA degenerate margin M=0 | Pathway with 0 annotated genes | p-value = 1.0 | INV-5 |
| S6 | ORA empty query set n=0 | No interesting genes | p-value = 1.0 (each pathway) | INV-5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Bit-score monotonicity | Higher raw score ⇒ higher bit score ⇒ lower E-value | Strictly increasing bit / decreasing E | INV-1, INV-3 |
| C2 | ORA p-value bounds | Random valid inputs | p ∈ [0,1] | INV-4 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/` contains MetagenomicsAnalyzer tests for AlphaDiversity, BetaDiversity, GenomeBinning, TaxonomicClassification, TaxonomicProfile, plus a legacy `MetagenomicsAnalyzerTests.cs`. None target `PredictFunctions` or `FindPathwayEnrichment`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | No bit-score test |
| M2 | ❌ Missing | No E-value test |
| M3 | ❌ Missing | No best-hit test |
| M4 | ❌ Missing | No ORA test (method did not exist) |
| M5 | ❌ Missing | — |
| M6 | ❌ Missing | — |
| S1 | ❌ Missing | — |
| S2 | ❌ Missing | — |
| S3 | ❌ Missing | — |
| S4 | ❌ Missing | — |
| S5 | ❌ Missing | — |
| S6 | ❌ Missing | — |
| C1 | ❌ Missing | — |
| C2 | ❌ Missing | — |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_PredictFunctions_Tests.cs` — all META-FUNC-001 cases.
- **Remove:** none (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| MetagenomicsAnalyzer_PredictFunctions_Tests.cs | Canonical META-FUNC-001 fixture | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented bit-score test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented E-value test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented best-hit test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented ORA worked-example test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented ORA no-overlap test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented annotation-fields test | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented empty-sequence test | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented no-match test | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented empty/null DB test | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented null-gene-list test | ✅ Done |
| 11 | S5 | ❌ Missing | Implemented M=0 margin test | ✅ Done |
| 12 | S6 | ❌ Missing | Implemented n=0 query test | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented monotonicity property test | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented p-value bounds property test | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Bit score 18.0202932787533 verified |
| M2 | ✅ Covered | E-value 3.3852730346546e-5 verified |
| M3 | ✅ Covered | Best hit = lower E-value entry |
| M4 | ✅ Covered | ORA = 7.88e-8 verified |
| M5 | ✅ Covered | x=0 ⇒ p=1 |
| M6 | ✅ Covered | Function/Pathway/KO transferred |
| S1 | ✅ Covered | Empty sequence ⇒ no annotation |
| S2 | ✅ Covered | No match ⇒ no annotation |
| S3 | ✅ Covered | Empty/null DB ⇒ no annotation |
| S4 | ✅ Covered | Null genes ⇒ ArgumentNullException |
| S5 | ✅ Covered | M=0 ⇒ p=1 |
| S6 | ✅ Covered | n=0 ⇒ p=1 |
| C1 | ✅ Covered | Monotonicity holds |
| C2 | ✅ Covered | p ∈ [0,1] |

Total in-scope cases: 14. ✅ count: 14.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Ungapped exact-match scoring model (signature occurrence → BLOSUM62 self-match score) | Implementation of `PredictFunctions`; affects which hits found, not the cited bit-score/E-value formulas |

---

## 7. Open Questions / Decisions

1. The Registry signature `PredictFunctions(genes, database)` is realized as `PredictFunctions(IEnumerable<(GeneId, ProteinSequence)>, IReadOnlyDictionary<signature,(Function,Pathway,Ko)>)` (the pre-existing API shape), now made evidence-based by computing source-backed bit score / E-value instead of the previous invented constants (`EValue: 1e-10`, `BitScore: motif.Length*2.0`).
2. `FindPathwayEnrichment(functions, pathwayDb)` did not previously exist; it is added per the Registry and implements the cited hypergeometric ORA test.
