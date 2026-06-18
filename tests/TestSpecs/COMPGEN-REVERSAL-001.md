# Test Specification: COMPGEN-REVERSAL-001

**Test Unit ID:** COMPGEN-REVERSAL-001
**Area:** Comparative
**Algorithm:** Reversal Distance (breakpoint-based lower bound)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Bafna V, Pevzner PA (1998). Sorting by Transpositions. SIAM J. Discrete Math. 11(2):224–240. | 1 | https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf | 2026-06-14 |
| 2 | Hunter College CompBio Lecture 16 — sorting by reversals. | 2 | https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf | 2026-06-14 |
| 3 | Hübotter J (2020). On Sorting by Reversals. | 4 | https://jonhue.github.io/min-sbr/paper.pdf | 2026-06-14 |
| 4 | Bergeron, Mixtacki, Stoye (2009). The Inversion Distance Problem. | 1 | https://gi.cebitec.uni-bielefeld.de/_media/teaching/2018winter/cg/inversionbergeron.pdf | 2026-06-14 |

### 1.2 Key Evidence Points

1. A pair (π_i, π_{i+1}) of the extended permutation is a breakpoint iff π_{i+1} ≠ π_i + 1; the identity is the only permutation with 0 breakpoints — Source 1, §2.
2. Extended permutation = (0, π_1, …, π_n, n+1) — Source 2.
3. A reversal removes at most two breakpoints: b(α) − b(αρ) ≤ 2, hence b(α) ≤ 2t and d(α) ≥ b(α)/2 — Source 2.
4. Unsigned breakpoint: |π_{i+1} − π_i| ≠ 1; same bound b(π)/2 ≤ d_r(π) — Source 3.
5. Reversal distance is symmetric: d_β(α) = d_α(β); target may be taken as identity WLOG — Source 2.
6. Worked example (signed) α=(−2,−3,+1,+6,−5,−4), b(α)=6 — Source 2.

### 1.3 Documented Corner Cases

- Identity permutation ⇒ 0 breakpoints ⇒ distance 0 (Source 1).
- The breakpoint bound is a *lower bound*, not exact distance (Source 1: "a permutation with few breakpoints may be more distant…"; Source 2: "This lower bound is not very tight.").

### 1.4 Known Failure Modes / Pitfalls

1. Treating the breakpoint bound as the exact reversal distance — it is only a lower bound (Sources 1, 2).
2. Reversal that removes 2 breakpoints may not reduce true distance ("not necessarily one that makes progress") — Source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateReversalDistance(IReadOnlyList<int> permutation1, IReadOnlyList<int> permutation2)` | ComparativeGenomics | **Canonical** | Unsigned breakpoint-based lower bound ⌈b/2⌉. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | d(π, π) = 0 for any permutation (identity ⇒ 0 breakpoints). | Yes | Source 1 §2 |
| INV-2 | Result ≥ 0 (a count of reversals). | Yes | Source 1/2 (d ≥ b/2 ≥ 0) |
| INV-3 | Symmetry: d(α, β) = d(β, α). | Yes | Source 2 ("d_β(α) = d_α(β)") |
| INV-4 | Result equals ⌈b/2⌉ where b is the unsigned breakpoint count of the extended relative permutation. | Yes | Sources 1, 2, 3 |
| INV-5 | Result is a lower bound: it never exceeds the number of reversals actually applied to build the input. | Yes | Source 2 (b(α) ≤ 2t ⇒ d ≥ b/2) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Identity equal orders | perm1 = perm2 = [1,2,3,4,5] | 0 | Source 1 §2 (identity ⇒ 0 breakpoints) |
| M2 | Hunter worked example (unsigned) | perm1=[2,3,1,6,5,4], perm2=[1,2,3,4,5,6] (b=4) | 2 | Source 2 worked example + unsigned bp def (b=4 ⇒ ⌈4/2⌉) |
| M3 | Fully reversed | perm1=[4,3,2,1], perm2=[1,2,3,4] (b=2) | 1 | Source 1 §2 bp def (b=2 ⇒ ⌈2/2⌉=1) |
| M4 | Single adjacent swap | perm1=[1,2,4,3], perm2=[1,2,3,4] (b=2, derivation in §7.2) | 1 | Source 1/3 unsigned bp def (b=2 ⇒ ⌈2/2⌉=1) |
| M5 | Lower-bound property | distance ≤ number of reversals applied to derive perm1 | result ≤ applied reversals | Source 2 (d ≥ b/2 ⇒ bound never exceeds true t) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty inputs | perm1 = perm2 = [] | 0 | n ≤ 1 ⇒ no breakpoint |
| S2 | Single element | perm1 = perm2 = [7] | 0 | n ≤ 1 ⇒ no breakpoint |
| S3 | Unequal lengths throw | perm1=[1,2], perm2=[1] | `ArgumentException` | Distance undefined across different marker sets |
| S4 | Symmetry | d([2,3,1,6,5,4],[1..6]) == d([1..6],[2,3,1,6,5,4]) | equal | INV-3 / Source 2 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Non-zero target labels | perm1=[30,10,20], perm2=[10,20,30] | depends on relative order | Confirms relative-permutation remapping works for arbitrary labels |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Existing `CalculateReversalDistance` lives in `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs`.
- No test file existed for this method (other `ComparativeGenomics_*` test files cover orthologs, synteny, rearrangements, RBH, CompareGenomes — none reference `CalculateReversalDistance`).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | No prior test |
| M2 | ❌ Missing | No prior test |
| M3 | ❌ Missing | No prior test |
| M4 | ❌ Missing | No prior test |
| M5 | ❌ Missing | No prior test |
| S1 | ❌ Missing | No prior test |
| S2 | ❌ Missing | No prior test |
| S3 | ❌ Missing | No prior test |
| S4 | ❌ Missing | No prior test |
| C1 | ❌ Missing | No prior test |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateReversalDistance_Tests.cs` — all cases for this unit.
- **Remove:** nothing (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_CalculateReversalDistance_Tests.cs` | Canonical | 10 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented (property) | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented | ✅ Done |
| 9 | S4 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Identity ⇒ 0 |
| M2 | ✅ Covered | Hunter unsigned example ⇒ 2 |
| M3 | ✅ Covered | Fully reversed ⇒ 1 |
| M4 | ✅ Covered | Single end swap ⇒ 1 |
| M5 | ✅ Covered | Property: bound ≤ applied reversals |
| S1 | ✅ Covered | Empty ⇒ 0 |
| S2 | ✅ Covered | Single ⇒ 0 |
| S3 | ✅ Covered | Unequal ⇒ ArgumentException |
| S4 | ✅ Covered | Symmetry |
| C1 | ✅ Covered | Arbitrary labels |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Integer rounding ⌈b/2⌉ = `(b+1)/2` is the integer form of d ≥ b/2 (tightest integer lower bound). | Implementation return; M2, M3, M4 |
| 2 | Unequal-length inputs throw `ArgumentException` (distance undefined across different marker sets; not separately specified by sources). | S3 |

---

## 7. Open Questions / Decisions

1. **Decision:** The method computes the **unsigned** breakpoint lower bound, not the exact signed Hannenhalli–Pevzner distance. The Hunter signed worked example (b=6) is therefore not used directly; the unsigned specialization (b=4) is. Implementation Status documented as **Simplified**. (M4 expected value derived below.)
2. **M4 derivation:** perm1=[1,2,4,3], target=[1,2,3,4] ⇒ relative perm = [0,1,3,2]; extended [−1,0,1,3,2,4]: boundaries −1→0 adj, 0→1 adj, 1→3 (BP), 3→2 adj, 2→4 (BP) ⇒ b=2 ⇒ ⌈2/2⌉ = **1**.
