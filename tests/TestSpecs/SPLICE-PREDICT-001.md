# Test Specification: SPLICE-PREDICT-001

**Test Unit ID:** SPLICE-PREDICT-001
**Area:** Splicing
**Algorithm:** Gene Structure Prediction
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-12

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| S1 | Gilbert (1978) | 1 | doi:10.1038/271501a0 | 2026-02-12 |
| S2 | Breathnach & Chambon (1981) | 1 | doi:10.1146/annurev.bi.50.070181.002025 | 2026-02-12 |
| S3 | Shapiro & Senapathy (1987) | 1 | doi:10.1093/nar/15.17.7155 | 2026-02-12 |
| S4 | Burge, Tuschl & Sharp (1999) | 1 | ISBN 0-87969-380-6 | 2026-02-12 |
| S5 | Wikipedia — Intron | 4 | https://en.wikipedia.org/wiki/Intron | 2026-02-12 |
| S6 | Wikipedia — Exon | 4 | https://en.wikipedia.org/wiki/Exon | 2026-02-12 |
| S7 | Sakharkar et al. (2002) | 1 | doi:10.1093/nar/30.1.191 | 2026-02-12 |
| S8 | Alberts et al. (2002) | 1 | Molecular Biology of the Cell, 4th ed. | 2026-02-12 |

### 1.2 Key Evidence Points

1. Exons are expressed regions retained in mature mRNA; introns are excised — S1
2. >99% of spliceosomal introns follow the GT-AG rule — S2
3. Donor consensus: MAG|GURAGU; Acceptor consensus: (Y)nNCAG|G — S2
4. Intron types: U2 (GT-AG, ~99%), U12 (AT-AC, ~0.5%), GC-AG (~0.5%) — S4
5. Exon types: Initial, Internal, Terminal, Single — S1, S6
6. Average 5.48 exons per protein-coding gene — S7
7. Exon phase = cumulative length of preceding exons mod 3 — S8

### 1.3 Documented Corner Cases

1. **Single-exon genes:** Genes with no introns produce a single exon — S6
2. **Empty input:** No sequence → no structure (trivial)
3. **Overlapping intron candidates:** Multiple donor-acceptor pairings may overlap;
   non-overlapping selection required — implementation heuristic
4. **Short sequences:** Sequences too short to contain both splice sites produce
   no introns

### 1.4 Known Failure Modes / Pitfalls

1. **Score threshold too high:** May miss real introns, defaulting to single-exon — S3
2. **Greedy selection suboptimal:** May select a high-scoring short intron that
   blocks a biologically correct longer one — implementation-specific

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictGeneStructure` | SpliceSitePredictor | Canonical | Primary orchestrator |
| `PredictIntrons` | SpliceSitePredictor | Canonical | Intron prediction by donor-acceptor pairing |
| `SelectNonOverlappingIntrons` | SpliceSitePredictor | Internal | Tested via PredictGeneStructure |
| `DeriveExons` | SpliceSitePredictor | Internal | Tested via PredictGeneStructure |
| `GenerateSplicedSequence` | SpliceSitePredictor | Internal | Tested via PredictGeneStructure |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Introns start at donor position, end at acceptor position | Yes | S2 |
| INV-2 | Selected introns do not overlap | Yes | Algorithm design |
| INV-3 | Exon positions + intron positions cover entire sequence | Yes | Definition |
| INV-4 | Spliced sequence = concatenation of exon sequences | Yes | S1 |
| INV-5 | Spliced sequence length = total sequence length - total intron length | Yes | S1 |
| INV-6 | Exon phase = (sum of preceding exon lengths) mod 3 | Yes | S8 |
| INV-7 | Overall score = mean of selected intron scores (0 if none) | Yes | Implementation |
| INV-8 | All intron lengths ≥ minIntronLength | Yes | Parameter contract |
| INV-9 | All intron lengths ≤ maxIntronLength | Yes | Parameter contract |
| INV-10 | Score values in [0, 1] | Yes | Normalization |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | EmptyNull_ReturnsEmpty | Empty/null → empty GeneStructure | 0 exons, 0 introns, "" spliced, 0 score | Trivial |
| M2 | SingleExon_NoIntrons | No splice sites found → single exon | 1 exon of type Single, 0 introns | S1, S6 |
| M3 | TwoExon_OneIntron | Clear GT-AG intron → 2 exons, 1 intron | Initial + Terminal exons, 1 intron | S2 |
| M4 | SplicedSequence_ExcludesIntron | Spliced seq = exons only | Length = total - intron length | S1 |
| M5 | IntronMinLength_Respected | Short introns filtered | No intron shorter than minIntronLength | S2, INV-8 |
| M6 | IntronMaxLength_Respected | Long introns filtered | No intron longer than maxIntronLength | INV-9 |
| M7 | IntronType_U2_GTAG | GT-AG intron classified as U2 | Intron.Type == U2 | S4 |
| M8 | ExonTypes_Correct | Multi-exon assigns Initial/Terminal | First=Initial, Last=Terminal | S1, S6 |
| M9 | ExonPhase_Tracks | Phase = cumulative length mod 3 | Phase values computed correctly | S8, INV-6 |
| M10 | ScoreRange | Intron and overall scores in [0,1] | All scores between 0 and 1 | INV-10 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | NonOverlapping_Selection | Overlapping introns → only non-overlapping selected | No position used twice | INV-2 |
| S2 | DNA_T_Equivalence | DNA input (T) same as RNA (U) | Same structure result | Implementation converts T→U |
| S3 | IntronSequence_Correct | Intron.Sequence matches substring | Sequence equals input[start..end] | Definition |
| S4 | Threshold_Filtering | Higher minScore → fewer introns | Subset relationship | Parameter semantics |
| S5 | CaseInsensitive | Lowercase input produces same result | Same structure | Implementation uses ToUpperInvariant |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | OverallScore_MeanOfIntrons | Overall score = mean of intron scores | Math verification | INV-7 |
| C2 | SplicedSequence_NoIntrons_EqualsInput | No introns → spliced = original | Identity | Trivial |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Tests found in `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictorTests.cs`
- 5 intron prediction tests, 5 gene structure tests in dedicated regions

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| PredictIntrons_ValidDonorAcceptor_FindsIntron | ⚠ Weak | Asserts Count >= 0 (always true) |
| PredictIntrons_RespectsMinLength | ⚠ Weak | Uses `‖ Count == 0` fallback |
| PredictIntrons_RespectsMaxLength | ⚠ Weak | Uses `‖ Count == 0` fallback |
| PredictIntrons_IncludesIntronSequence | ⚠ Weak | Guarded by `if`, asserts Not.Empty |
| PredictIntrons_ClassifiesIntronType | ⚠ Weak | Guarded by `if`, allows Unknown |
| PredictGeneStructure_SingleExon_NoIntrons | ⚠ Weak | Asserts Count == 0 OR > 0 (tautology) |
| PredictGeneStructure_ReturnsExons | ⚠ Weak | Asserts Is.Not.Null only |
| PredictGeneStructure_GeneratesSplicedSequence | ⚠ Weak | Guarded by `if` |
| PredictGeneStructure_EmptySequence_ReturnsEmpty | ✅ Covered | Good empty check |
| PredictGeneStructure_ExonPhase_Calculated | ⚠ Weak | Not.Null.Or.EqualTo(null) (tautology) |
| M3 TwoExon OneIntron | ❌ Missing | |
| M4 SplicedSequence excludes intron | ❌ Missing | |
| M7 IntronType U2 | ❌ Missing | |
| M8 ExonTypes assignment | ❌ Missing | |
| M9 ExonPhase calculation | ❌ Missing | |
| M10 Score range | ❌ Missing | |
| S1 NonOverlapping | ❌ Missing | |
| S2 DNA T equivalence | ❌ Missing | |
| S3 IntronSequence | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_GeneStructure_Tests.cs`
  — all SPLICE-PREDICT-001 tests consolidated here
- **Remove from shared file:** 5 intron tests + 5 gene structure tests from
  `SpliceSitePredictorTests.cs` (all weak)
- **Retain in shared file:** Alternative splicing, branch point, MaxEntScore,
  IsWithinCodingRegion, input handling, integration tests (other Test Units)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| SpliceSitePredictor_GeneStructure_Tests.cs | Canonical for SPLICE-PREDICT-001 | 17 |
| SpliceSitePredictorTests.cs | Shared (other TUs) | ~19 |
| SpliceSitePredictor_DonorSite_Tests.cs | SPLICE-DONOR-001 | 17 |
| SpliceSitePredictor_AcceptorSite_Tests.cs | SPLICE-ACCEPTOR-001 | 17 |

---

## 6. Assumption Register

**Total assumptions:** 4

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Greedy non-overlapping selection is acceptable heuristic | S1, Algorithm doc |
| A2 | Exon phase starts at 0 (no UTR offset) | M9 |
| A3 | Overall score = arithmetic mean of intron scores | C1 |
| A4 | Default branch point score of 0.3 when none found | Scoring tests |

---

## 7. Open Questions / Decisions

None — all behaviors confirmed by evidence or documented as assumptions.
