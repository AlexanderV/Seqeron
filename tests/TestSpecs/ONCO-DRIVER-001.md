# Test Specification: ONCO-DRIVER-001

**Test Unit ID:** ONCO-DRIVER-001
**Area:** Oncology
**Algorithm:** Driver Mutation Detection (20/20 rule)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Vogelstein et al. 2013, Cancer Genome Landscapes, Science | 1 | https://doi.org/10.1126/science.1235122 | 2026-06-14 |
| 2 | Tokheim & Karchin 2020, 20/20+, Bioinformatics | 1 | https://doi.org/10.1093/bioinformatics/btz759 | 2026-06-14 |
| 3 | Schroeder et al. 2014, OncodriveROLE, Bioinformatics | 1 | https://doi.org/10.1093/bioinformatics/btu467 | 2026-06-14 |
| 4 | Miller et al. 2017, mutational hotspots, Oncotarget | 1 | https://doi.org/10.18632/oncotarget.15514 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Oncogene if >20% of a gene's mutations are missense at recurrent positions — Vogelstein 2013; Tokheim 2020 ("OGs have >20% mutations causing missense changes at recurrent positions").
2. Tumor suppressor if >20% of a gene's mutations are inactivating/truncating — Vogelstein 2013; Tokheim 2020 ("TSGs have >20% mutations causing inactivating changes").
3. Truncating types = nonsense (stop gain/loss), frameshift indel, splice donor/acceptor — Schroeder 2014; Miller 2017.
4. A recurrent position = same protein position observed ≥ 2 times — Miller 2017.
5. IDH1 = oncogene archetype (all mutations at codon 132) — Vogelstein 2013 / Miller 2017.

### 1.3 Documented Corner Cases

- Low-recurrence genes may satisfy neither criterion and stay unclassified (OncodriveROLE 2014).
- Passenger truncations can appear in oncogenes; rule is a heuristic, not a test (Tokheim 2020).
- Single-codon recurrence (IDH1) → recurrent-missense fraction ≈ 1.0 (Vogelstein 2013).

### 1.4 Known Failure Modes / Pitfalls

1. Lowly recurrent drivers missed — Schroeder 2014.
2. Drift-induced truncations in OGs mislead — Tokheim 2020.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifyGene(mutations)` | OncologyAnalyzer | Canonical | Applies the 20/20 rule to one gene's mutations → DriverGeneRole |
| `IdentifyDriverMutations(variants, hotspots)` | OncologyAnalyzer | Canonical | Per-gene 20/20 classification + hotspot match; returns driver subset |
| `ScoreDriverPotential(mutations)` | OncologyAnalyzer | Canonical | 20/20 driver-signal fraction in [0,1] (max of the two criterion fractions) |
| `MatchCancerHotspots(variant, hotspots)` | OncologyAnalyzer | Canonical | Caller-supplied hotspot-set membership lookup |

<!-- Type values: -->
<!-- **Canonical** — deep evidence-based testing -->
<!-- **Delegate** — smoke verification only (1–2 tests proving delegation) -->
<!-- **Internal** — tested indirectly via canonical methods -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | driver_mutations ⊆ input variants | Yes | Registry invariant; Vogelstein 2013 |
| INV-2 | Oncogene ⟺ recurrent-missense fraction > 0.20 (and ≥ TSG fraction) | Yes | Vogelstein 2013; Tokheim 2020 |
| INV-3 | TumorSuppressor ⟺ truncating fraction > 0.20 (and > OG fraction) | Yes | Vogelstein 2013; Tokheim 2020 |
| INV-4 | 0 ≤ ScoreDriverPotential ≤ 1 | Yes | fraction of mutations, bounded |
| INV-5 | A recurrent missense position requires the same position ≥ 2 times | Yes | Miller 2017 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | IDH1 all-recurrent missense | 10 missense all at codon 132 | role=Oncogene; recurrentMissenseFraction=1.00; truncatingFraction=0.00 | Vogelstein 2013; Miller 2017 |
| M2 | Dispersed truncating TSG | 5 nonsense + 2 frameshift + 1 splice (distinct pos) + 2 missense (distinct) | role=TumorSuppressor; truncatingFraction=0.80 | Vogelstein 2013; Schroeder 2014 |
| M3 | TSG threshold boundary (=0.20) | 2 truncating + 8 distinct missense | role≠TumorSuppressor (strict >0.20) | Vogelstein 2013 ">20%" |
| M4 | TSG just above boundary (0.30) | 3 truncating + 7 distinct missense | role=TumorSuppressor; truncatingFraction=0.30 | Vogelstein 2013 |
| M5 | MatchCancerHotspots hit | variant at a (gene,pos) in hotspot set | returns true | Miller 2017 recurrent-position |
| M6 | MatchCancerHotspots miss | variant not in hotspot set | returns false | Miller 2017 |
| M7 | IdentifyDriverMutations subset | mixed driver + non-driver variants | returns only variants in driver genes; ⊆ input | INV-1; Vogelstein 2013 |
| M8 | ScoreDriverPotential = recurrent-missense fraction | IDH1 set | 1.00 within 1e-10 | Vogelstein 2013 |
| M9 | ScoreDriverPotential = truncating fraction | dispersed TSG set | 0.80 within 1e-10 | Vogelstein 2013 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Low-recurrence → Ambiguous | 5 muts: 1 nonsense (frac 0.20) + 4 missense at distinct pos (recurrent frac 0.00) | role=Ambiguous (neither fraction >0.20) | OncodriveROLE failure mode |
| S2 | Singleton missense not recurrent | 3 missense at 3 distinct positions | recurrentMissenseFraction=0.00 → not Oncogene | Miller 2017 INV-5 |
| S3 | Recurrence needs ≥2 at same pos | 2 missense same pos + 8 distinct | recurrentMissenseFraction=0.20 → not >0.20 → not Oncogene | Miller 2017; strict > |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty gene mutations | ClassifyGene([]) | role=Ambiguous, both fractions 0 | no evidence |
| C2 | Null input | null variants / mutations | ArgumentNullException | contract |
| C3 | Property: driver ⊆ input over random-ish panel | mixed panel | every returned variant present in input | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- ONCO-DRIVER-001 is a new unit; no prior tests existed for `IdentifyDriverMutations`/`ClassifyGene`/`ScoreDriverPotential`/`MatchCancerHotspots`. All planned cases start ❌ Missing.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M9 | ❌ Missing | new unit |
| S1–S3 | ❌ Missing | new unit |
| C1–C3 | ❌ Missing | new unit |

<!-- Status values: ✅ Covered, ⚠ Weak, ❌ Missing, 🔁 Duplicate -->

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_IdentifyDriverMutations_Tests.cs` — all cases for this unit.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_IdentifyDriverMutations_Tests.cs` | Canonical fixture for ONCO-DRIVER-001 | 20 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | ClassifyGene_AllRecurrentMissense_ClassifiesOncogene | ✅ Done |
| 2 | M2 | ❌ Missing | ClassifyGene_DispersedTruncating_ClassifiesTumorSuppressor | ✅ Done |
| 3 | M3 | ❌ Missing | ClassifyGene_TruncatingFractionAtThreshold_IsNotTumorSuppressor | ✅ Done |
| 4 | M4 | ❌ Missing | ClassifyGene_TruncatingFractionAboveThreshold_ClassifiesTumorSuppressor | ✅ Done |
| 5 | M5 | ❌ Missing | MatchCancerHotspots_KnownHotspot_ReturnsTrue (+ integration) | ✅ Done |
| 6 | M6 | ❌ Missing | MatchCancerHotspots_NotInSet_ReturnsFalse | ✅ Done |
| 7 | M7 | ❌ Missing | IdentifyDriverMutations_MixedPanel_ReturnsOnlyDriverGeneMutations | ✅ Done |
| 8 | M8 | ❌ Missing | ScoreDriverPotential_RecurrentMissenseGene_ReturnsRecurrentFraction | ✅ Done |
| 9 | M9 | ❌ Missing | ScoreDriverPotential_TruncatingGene_ReturnsTruncatingFraction | ✅ Done |
| 10 | S1 | ❌ Missing | ClassifyGene_LowRecurrence_ClassifiesAmbiguous | ✅ Done |
| 11 | S2 | ❌ Missing | ClassifyGene_AllMissenseDistinctPositions_NotRecurrentNotOncogene | ✅ Done |
| 12 | S3 | ❌ Missing | ClassifyGene_RecurrentMissenseFractionAtThreshold_IsNotOncogene | ✅ Done |
| 13 | C1 | ❌ Missing | ClassifyGene_EmptyMutations_ReturnsAmbiguousZeroFractions | ✅ Done |
| 14 | C2 | ❌ Missing | *_NullInput_Throws (ClassifyGene/Score/Identify/MatchHotspots) | ✅ Done |
| 15 | C3 | ❌ Missing | IdentifyDriverMutations_Property_ResultIsSubsetInInputOrder | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Oncogene call, f_OG=1.00 |
| M2 | ✅ | TSG call, f_TSG=0.80 |
| M3 | ✅ | exactly 0.20 not TSG |
| M4 | ✅ | 0.30 → TSG |
| M5 | ✅ | hotspot hit (+ rescue integration) |
| M6 | ✅ | hotspot miss (gene & position) |
| M7 | ✅ | driver subset, INV-01 |
| M8 | ✅ | score = f_OG |
| M9 | ✅ | score = f_TSG |
| S1 | ✅ | low recurrence → Ambiguous |
| S2 | ✅ | distinct missense not recurrent |
| S3 | ✅ | f_OG exactly 0.20 not Oncogene |
| C1 | ✅ | empty → Ambiguous, fractions 0 |
| C2 | ✅ | null guards on all four methods |
| C3 | ✅ | property: subset + order |

**Total in-scope cases: 15 — ✅ 15.**

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Both-criteria tie-break by dominant fraction; exact tie → Ambiguous | ClassifyGene |
| 2 | Strict `>` 0.20 threshold (Vogelstein/Tokheim wording) | M3, S3, ClassifyGene |
| 3 | ScoreDriverPotential = 20/20 driver-signal fraction (no CADD/SIFT/PolyPhen model) | ScoreDriverPotential |

---

## 7. Open Questions / Decisions

1. CADD/SIFT/PolyPhen pathogenicity scores are external trained models; not reproducible from retrievable formulas. Decision: ScoreDriverPotential returns the transparent 20/20 driver-signal fraction; external scores are caller-supplied / not implemented (documented in §5.3 of the algorithm doc).
