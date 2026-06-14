# Test Specification: COMPGEN-ANI-001

**Test Unit ID:** COMPGEN-ANI-001
**Area:** Comparative
**Algorithm:** Average Nucleotide Identity (ANI), ANIb (Goris et al. 2007)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Goris et al. (2007) DNA-DNA hybridization values and whole-genome similarities. IJSEM 57:81-91 | 1 | https://doi.org/10.1099/ijs.0.64483-0 | 2026-06-14 |
| 2 | Konstantinidis & Tiedje (2005) Genomic insights that advance the species definition. PNAS 102:2567-2572 | 1 | https://doi.org/10.1073/pnas.0409727102 | 2026-06-14 |
| 3 | Lee et al. (2016) OrthoANI. IJSEM 66:1100-1103 | 1 | https://doi.org/10.1099/ijsem.0.000760 | 2026-06-14 |
| 4 | pyani ANIb reference implementation (Pritchard et al.) | 3 | https://pyani.readthedocs.io/en/latest/api/pyani.anib.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. Query genome is cut into consecutive (non-overlapping) 1020 nt fragments. — Goris 2007 (Methods)
2. Each fragment's best match against the reference is found via BLASTN; identity is **recalculated over the whole fragment length**. — Goris 2007 / pyani
3. ANI = mean identity of matches with **>30 % overall identity over an alignable region of ≥70 % of their length**. — Goris 2007 (verbatim)
4. ANI ≈ 95 % (Goris) / ≈ 94 % (Konstantinidis & Tiedje) corresponds to the 70 % DDH species boundary. — sources 1, 2

### 1.3 Documented Corner Cases

- Fragments below the identity or alignable-region cut-off are discarded and do not enter the mean (Goris 2007 / pyani).
- Per-fragment identity is over the fragment length, not over only the aligned sub-region (Goris 2007 "recalculated to an identity along the entire sequence").
- ANI is direction-dependent (query is the fragmented genome); pyani notes non-symmetrical matrices.

### 1.4 Known Failure Modes / Pitfalls

1. Computing identity as longest-common-substring length over fragment length (NOT nucleotide identity) — the pre-existing implementation defect this unit corrects. — Goris 2007 definition.
2. Counting all fragments instead of only conserved (qualifying) fragments overstates divergence. — Goris 2007.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction)` | ComparativeGenomics | **Canonical** | Mean per-fragment nucleotide identity under the 30 %/70 % cut-offs |
| `BestUngappedFragmentMatch(...)` | ComparativeGenomics (private) | **Internal** | Tested indirectly via CalculateANI |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ANI is a fraction in [0, 1] | Yes | Identity = matches/fragmentLength ∈ [0,1]; mean preserves bound (Goris 2007) |
| INV-2 | Identical sequences → ANI = 1.0 | Yes | Every fragment is a perfect substring, identity 1.0 (Goris 2007) |
| INV-3 | Only fragments with identity > minIdentity AND alignable fraction ≥ minAlignableFraction contribute | Yes | Goris 2007 cut-off clause |
| INV-4 | Fragmentation is consecutive and non-overlapping; trailing partial fragment (< fragmentLength) is ignored | Yes | Goris 2007 "consecutive 1020 nt fragments" |
| INV-5 | No qualifying fragment / empty / null input → 0 | Yes | Definition (mean over empty set undefined → 0) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Identical genomes | query = reference, fragLen 4 | ANI = 1.0 | INV-2, Goris 2007 |
| M2 | One substituted base | query `AAAACCCCGGGGTTTA` vs R, fragLen 4; last frag `TTTA`=3/4 | ANI = 0.9375 | Goris 2007 recalculated identity |
| M3 | Half-identity fragment | query `AAAACCCCGGGGAATT` vs R; last frag `AATT`=2/4 | ANI = 0.875 | Goris 2007 mean identity |
| M4 | Identity cut-off excludes fragment | query `AAAACGTC`, ref `AAAAAAAA`, fragLen 4; frag2=0.0 not >0.30 | ANI = 1.0 (only frag1) | Goris 2007 ">30 % identity" |
| M5 | Alignable cut-off excludes fragment | query `AAAA`, ref `AA` (ref < frag) | ANI = 0 | Goris 2007 "≥70 % alignable" |
| M6 | Consecutive non-overlapping fragmentation | query of length 10, fragLen 4 → 2 fragments, trailing 2 nt ignored | uses exactly 2 fragments | INV-4, Goris 2007 |
| M7 | Null / empty inputs | null or "" for either genome | ANI = 0 | Validation contract |
| M8 | Non-positive fragmentLength | fragmentLength = 0 | ArgumentOutOfRangeException | Validation contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Range invariant | ANI of arbitrary divergent pair | 0 ≤ ANI ≤ 1 | INV-1 |
| S2 | Query shorter than fragment | query length 3, fragLen 4 | ANI = 0 | No fragment fits |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Custom minIdentity keeps low-identity frag | lower minIdentity below a frag's identity | frag now contributes | Parameter exposure |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for `CalculateANI`. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `CalculateANI` / `ANI` — none found. Pre-existing implementation in `ComparativeGenomics.cs` used an LCS-length proxy (defect), with no tests.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8 | ❌ Missing | No tests existed for ANI |
| S1, S2 | ❌ Missing | No tests existed |
| C1 | ❌ Missing | No tests existed |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs` — all cases for this unit.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| ComparativeGenomics_CalculateANI_Tests.cs | Canonical | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented identical-genomes test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented one-mismatch (0.9375) test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented half-identity (0.875) test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented identity-cutoff exclusion test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented alignable-cutoff exclusion test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented non-overlapping fragmentation test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented null/empty tests | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented invalid fragmentLength test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented range invariant property test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented query-shorter-than-fragment test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented custom minIdentity test | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact value 1.0 |
| M2 | ✅ Covered | Exact value 0.9375 |
| M3 | ✅ Covered | Exact value 0.875 |
| M4 | ✅ Covered | Exact value 1.0 with frag2 excluded |
| M5 | ✅ Covered | Exact value 0.0 |
| M6 | ✅ Covered | 2 fragments used |
| M7 | ✅ Covered | 0 for null/empty |
| M8 | ✅ Covered | Throws ArgumentOutOfRangeException |
| S1 | ✅ Covered | 0 ≤ ANI ≤ 1 |
| S2 | ✅ Covered | 0.0 |
| C1 | ✅ Covered | Custom minIdentity |

**In-scope cases:** 11 — **✅:** 11

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Ungapped local alignment for the per-fragment best match (vs gapped BLASTN); per-fragment identity = matching bases / fragmentLength. The mean-identity ANI math under the 30 %/70 % cut-offs is implemented and tested exactly. | Implementation §5.3; M1–M5 |

---

## 7. Open Questions / Decisions

1. The reciprocal/averaged ANI (mean of both directions) and gapped-alignment refinement are out of scope for this unit (single-direction, ungapped); documented as a simplification in the algorithm doc §5.3. The well-defined mean-identity math is the contract under test.
