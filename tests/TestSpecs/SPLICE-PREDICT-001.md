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

- Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_GeneStructure_Tests.cs`
- All SPLICE-PREDICT-001 tests consolidated in canonical file
- 21 tests covering all MUST, SHOULD, and COULD cases (0 missing, 0 weak, 0 duplicate)

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 EmptyNull_ReturnsEmpty | ✅ Covered | 2 tests: asserts all fields of empty structure (empty + null) |
| M2 SingleExon_NoIntrons | ✅ Covered | Exact count, type=Single, start=0, end=len−1 |
| M3 TwoExon OneIntron | ✅ Covered | Exact counts (1 intron, 2 exons), exact intron position/length, GU…AG sequence, INV-3 coverage |
| M4 SplicedSequence excludes intron | ✅ Covered | INV-4 (concat of exon sequences), INV-5 (length check), exact Exon1+Exon2 |
| M5 IntronMinLength | ✅ Covered | Non-vacuous: shows introns found at 60, filtered at 200 |
| M6 IntronMaxLength | ✅ Covered | Non-vacuous: shows introns found at default, filtered at 70 |
| M7 IntronType U2 | ✅ Covered | Asserts U2 for GU-donor introns per Burge et al. (1999) |
| M8 ExonTypes assignment | ✅ Covered | Exact count prerequisite, Initial + Terminal |
| M9 ExonPhase calculation | ✅ Covered | Exact count prerequisite, Phase[0]=0, Phase[1]=Exon1.Length%3=2, general formula |
| M10 Score range | ✅ Covered | 2 tests: non-empty prerequisite + all scores in [0,1] |
| S1 NonOverlapping | ✅ Covered | No position overlap in selected introns |
| S2 DNA T equivalence | ✅ Covered | Same intron/exon counts, spliced sequence, and overall score |
| S3 IntronSequence | ✅ Covered | Sequence matches input substring |
| S4 Threshold_Filtering | ✅ Covered | Higher minScore → fewer introns |
| S5 CaseInsensitive | ✅ Covered | Same intron/exon counts, spliced sequence, and overall score |
| C1 OverallScore_MeanOfIntrons | ✅ Covered | Non-empty prerequisite + math verification (ε=1e-10) |
| C2 SplicedSequence_NoIntrons_EqualsInput | ✅ Covered | Explicit 0-intron assertion + identity check |

### 5.3 Final State

| File | Role | Test Count |
|------|------|------------|
| SpliceSitePredictor_GeneStructure_Tests.cs | Canonical for SPLICE-PREDICT-001 | 21 |
| SpliceSitePredictorTests.cs | Shared (other TUs) | ~19 |
| SpliceSitePredictor_DonorSite_Tests.cs | SPLICE-DONOR-001 | 17 |
| SpliceSitePredictor_AcceptorSite_Tests.cs | SPLICE-ACCEPTOR-001 | 17 |

---

## 6. Assumption Register

**Total assumptions:** 0

All former assumptions have been resolved:

| # | Former Assumption | Resolution |
|---|-------------------|------------|
| ~A1~ | Greedy non-overlapping selection | **Design Decision** — external sources prescribe the biological constraint (introns don't overlap — S2, S4) but not the computational algorithm. Greedy-by-score is a valid simplified approach; no literature standard exists for overlap resolution in computational predictors. |
| ~A2~ | Exon phase starts at 0 | **Mathematical Fact** — phase = (sum of preceding exon lengths) mod 3; with 0 preceding exons the result is trivially 0. Standard convention per Alberts et al. (2002). Removed. |
| ~A3~ | Overall score = arithmetic mean | **Definition** — no biological standard exists for a gene structure quality metric. The arithmetic mean of intron scores is a defined internal metric, not a biological assumption. |
| ~A4~ | Default branch point score 0.3 | **Fixed** — arbitrary magic constant removed from code. When no branch point is found, combined score is now (donor + acceptor) / 2 instead of (donor + acceptor + 0.3) / 3. |

---

## 7. Open Questions / Decisions

None — all behaviors confirmed by evidence. No residual assumptions.
