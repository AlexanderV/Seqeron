# Test Specification: ONCO-NEO-001

**Test Unit ID:** ONCO-NEO-001
**Area:** Oncology
**Algorithm:** Neoantigen Candidate Peptide Window Generation (somatic missense mutation)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Hundal et al. (2020), pVACtools, Cancer Immunol Res 8(3):409–420 | 1 | https://doi.org/10.1158/2326-6066.CIR-19-0401 | 2026-06-14 |
| 2 | Li et al. (2020), ProGeo-neo, BMC Med Genomics 13:52 | 1 | https://doi.org/10.1186/s12920-020-0683-4 | 2026-06-14 |
| 3 | Jurtz et al. (2017), NetMHCpan-4.0, J Immunol 199(9):3360–3368 | 1 | https://doi.org/10.4049/jimmunol.1700893 | 2026-06-14 |
| 4 | NetMHCpan-4.1 web service (DTU) | 3 | https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/ | 2026-06-14 |
| 5 | Wells et al. (2020), TESLA, Cell 183(3):818–834 | 1 | https://doi.org/10.1016/j.cell.2020.09.015 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Class I candidate peptides are 8–11-mers — "8–11-mer for Class I MHC" — Hundal et al. (2020); class I benchmark "lengths 8-11" — pVACtools docs.
2. Candidate peptides are built by centring the substitution with flanks "on each side" of the mutation (21-mer, 10 flanks each side) and extracting the windows that overlap the mutated residue — Li et al. (2020).
3. Each mutant peptide is paired with the wild-type peptide at the same coordinates (the agretope); the two differ only at the substituted residue — Wells et al. (2020) DAI; Hundal et al. (2020).
4. Length 9 is the dominant class I ligand length but enumeration spans 8–11 because length preference varies by allele — NetMHCpan-4.1 service.

### 1.3 Documented Corner Cases

- Mutation near a terminus: fewer than k windows of length k exist; only the windows that fit while spanning the mutation are produced (pVACtools centres / ProGeo-neo builds 21-mer "if possible").
- Variant producing no amino-acid change: no candidate peptide (only protein-altering variants — pVACtools).
- Length outside 8–11: outside the class I canonical band (NetMHCpan).

### 1.4 Known Failure Modes / Pitfalls

1. Off-by-one in window enumeration (missing the window ending exactly at the mutation, or one past the end) — derived from the spanning definition.
2. Treating a non-substitution (mutant == wild-type) as a neoantigen — would emit identical mutant/WT peptides; rejected as invalid input.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GenerateNeoantigenPeptides(wildTypeProtein, mutantResidue, mutationPosition, minLength, maxLength)` | OncologyAnalyzer | **Canonical** | Window enumeration + agretope pairing |

> Binding-affinity scoring (`ScoreNeoantigens` / IC50) is intentionally NOT implemented in this unit
> (caller-supplied / ONCO-MHC-001) — see Evidence Assumption 2.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned peptide has `Length` in `[minLength, maxLength]` and 8–11 by default | Yes | Hundal (2020) |
| INV-2 | Every peptide spans the mutation: `0 ≤ MutationOffset < Length` and `StartPosition + MutationOffset == mutationPosition` | Yes | Li (2020) |
| INV-3 | `MutantPeptide` and `WildTypePeptide` have equal length and differ at exactly one index — the `MutationOffset` | Yes | Wells (2020); Hundal (2020) |
| INV-4 | `MutantPeptide[MutationOffset]` == mutant residue; `WildTypePeptide[MutationOffset]` == original wild-type residue | Yes | Wells (2020) |
| INV-5 | For length k with the mutation ≥ k−1 residues from both ends, exactly k windows of length k are produced | Yes | windowing definition (Li 2020) |
| INV-6 | Peptides are ordered by length ascending, then start position ascending; count never exceeds Σ k over the range | Yes | implementation contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

Protein `MKTAYIAKQRSTVWLNDEFGH` (L=21), missense `Y5C` unless stated.

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Interior k=9 window count | One length, interior mutation | exactly 9 windows? No — p=5 is < k−1 from start, so 5 windows; start positions 1..5; each spans position 5 | Li (2020) windowing |
| M2 | Default range total | k=8..11, Y5C | 20 peptides total (5 per length) | Hundal (2020) 8–11; Li (2020) |
| M3 | Mutant/WT pairing | First 8-mer window | mutant `MKTACIAK`, WT `MKTAYIAK`, offset 4, differ only at index 4 (C vs Y) | Wells (2020); Hundal (2020) |
| M4 | All windows span mutation | every peptide in default result | `StartPosition + MutationOffset == 5`; `MutantPeptide[offset]=='C'`, `WildTypePeptide[offset]=='Y'` | Li (2020); Wells (2020) |
| M5 | Fully interior k=9 count | mutation `V10A` (p=10) so ≥ k−1 from both ends, k=9 | exactly 9 windows of length 9 | Li (2020); INV-5 |
| M6 | Terminal mutation truncation | `M1V` (p=1), k=9 | exactly 1 window: mutant `VKTAYIAKQ`, WT `MKTAYIAKQ`, start 1, offset 0 | ProGeo-neo "if possible" |
| M7 | C-terminal mutation | `H21R` (p=21), k=8 | exactly 1 window ending at 21: start 14, mutant ends `…NDEFGR`, offset 7 | ProGeo-neo "if possible" |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null protein | `wildTypeProtein` null | `ArgumentNullException` | input validation |
| S2 | Empty protein | `""` | `ArgumentException` | input validation |
| S3 | Non-substitution | mutant residue == wild-type at position | `ArgumentException` | failure mode 2 |
| S4 | Position out of range | position 0 or 22 (L=21) | `ArgumentOutOfRangeException` | bounds |
| S5 | Invalid length range | minLength 0, or maxLength < minLength | `ArgumentException` | bounds |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Length exceeds protein | protein length 5, request k=8..11 | empty result (no window fits) | bounds handling |
| C2 | Single length subset | minLength=maxLength=10, Y5C | only 10-mers returned, all length 10 | range respected |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior tests for `GenerateNeoantigenPeptides` (new method). Sibling fixtures live in
  `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_*_Tests.cs`. New canonical file created:
  `OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs` — all cases above.
- **Remove:** none (no pre-existing tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_GenerateNeoantigenPeptides_Tests.cs | Canonical | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (property over result set) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | window count + start positions asserted |
| M2 | ✅ | total = 20 asserted |
| M3 | ✅ | exact mutant/WT strings + offset |
| M4 | ✅ | property over all peptides |
| M5 | ✅ | interior k=9 → 9 windows |
| M6 | ✅ | N-terminal single window |
| M7 | ✅ | C-terminal single window |
| S1 | ✅ | ArgumentNullException |
| S2 | ✅ | ArgumentException |
| S3 | ✅ | ArgumentException |
| S4 | ✅ | ArgumentOutOfRangeException (×2) |
| S5 | ✅ | ArgumentException (×2) |
| C1 | ✅ | empty result |
| C2 | ✅ | single length subset |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Single-residue missense substitution only (no indel/frameshift/fusion neopeptides) | scope of method |
| 2 | Binding affinity / IC50 out of scope (caller-supplied, ONCO-MHC-001) | scope of method |

---

## 7. Open Questions / Decisions

1. Decision: class name is `OncologyAnalyzer` (per task instruction and existing project layout), not the
   `NeoantigenPredictor` placeholder named in the checklist Registry — recorded as a conflict; external/task
   requirement wins. Method placed alongside sibling oncology methods.
2. Decision: binding-affinity scoring deferred to ONCO-MHC-001; this unit produces only the well-defined
   candidate-peptide windows + agretope pairs (no fabricated MHC model).
