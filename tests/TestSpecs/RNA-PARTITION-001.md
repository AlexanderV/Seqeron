# Test Specification: RNA-PARTITION-001

**Test Unit ID:** RNA-PARTITION-001
**Area:** RnaStructure
**Algorithm:** RNA Partition Function (McCaskill) and Boltzmann Structure Probability
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | McCaskill JS (1990) *Biopolymers* 29:1105-1119 | 1 | https://doi.org/10.1002/bip.360290621 (PMID 1695107, via https://pubmed.ncbi.nlm.nih.gov/1695107/) | 2026-06-14 |
| 2 | Will S, MIT 18.417 — McCaskill inside recursion slides | 1 | https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf | 2026-06-14 |
| 3 | Will S, MIT 18.417 — McCaskill base-pair-probability slides | 1 | https://math.mit.edu/classes/18.417/Slides/mccaskill2.pdf | 2026-06-14 |
| 4 | Freiburg RNA Tools — McCaskill teaching tool (simplified model) | 3 | https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=McCaskill | 2026-06-14 |
| 5 | ViennaRNA — Partition Function reference (pf_fold) | 3 | https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/pf_fold.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. `Z = Σ_S exp(−E(S)/RT)` over all pseudoknot-free structures; `Z = Q_{1n}` — Source 1, 2, 4.
2. Inside recursion `Q_ij = Q_{i,j-1} + Σ_{i≤k<j−m} Q_{i,k-1}·Q^b_{kj}`, base case `Q_ij = 1` for `i ≥ j−m` — Source 2.
3. `Q^b_ij` is the partition function of `[i..j]` restricted to structures with `(i,j)` paired; in the simplified model each pair contributes a fixed `E_bp` so `Q^b_ij = exp(−β E_bp)·Q_{i+1,j-1}` — Source 2, 4.
4. Base-pair probability (external decomposition) `P[i,j] = Q_{1,i-1}·Q^b_{ij}·Q_{j+1,n} / Q_{1n}` — Source 3.
5. Structure probability `Pr[P|S] = Z⁻¹ exp(−βE(P))`, i.e. `p = exp(−βE)/Z`, with `β = 1/RT`, `R ≈ 1.987e-3 kcal/(mol·K)` — Source 1, 5.
6. Complexity O(n³) time, O(n²) space — Source 1, 2.
7. Default folding temperature 37 °C = 310.15 K — Source 5.

### 1.3 Documented Corner Cases

- Pairs with `j − i ≤ m` (m = minimum hairpin loop) are forbidden → `Q^b = 0` (Source 2).
- A sub-sequence too short to contain a pair has `Q = 1` (one empty structure) (Source 2).
- Only Watson-Crick (A-U, G-C) and GU pairs are admissible (Source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Ambiguous decomposition would over-count structures; the `Q`/`Q^b` split must be disjoint (Source 2).
2. Using the full Turner loop model vs. the simplified fixed-per-pair `E_bp` changes the numeric `Z` but not the invariants — documented as ASSUMPTION (Source 4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculatePartitionFunction(string, double, double)` | RnaSecondaryStructure | **Canonical** | McCaskill O(n³) DP; returns Z and base-pair probabilities |
| `CalculateStructureProbability(double, double, double)` | RnaSecondaryStructure | **Canonical** | Boltzmann probability `p = exp(−βE)/Z` |
| `GenerateRandomRna(int, double)` | RnaSecondaryStructure | **Delegate** | Forwards to seeded overload with `new Random()` |
| `GenerateRandomRna(int, Random, double)` | RnaSecondaryStructure | **Canonical** | Deterministic random generation (seeded) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `Z ≥ 1` for any sequence (the empty structure always contributes weight 1) | Yes | Source 2 (base case `Q=1`) |
| INV-2 | Every base-pair probability `P[i,j] ∈ [0,1]` | Yes | Source 1, 3 (it is a probability) |
| INV-3 | When `E_bp = 0`, `Z` equals the count of admissible pseudoknot-free structures | Yes | Source 2, 4 (`exp(0)=1`) |
| INV-4 | `Z` is strictly monotonically increasing as `E_bp` decreases (more favourable pairing) | Yes | Sum of increasing exponential weights (Source 1) |
| INV-5 | `CalculateStructureProbability(E, E) = 1` and `= exp(−ΔE/RT)` for `E_ensemble−E_struct=−ΔE` | Yes | Source 1, 5 (`p=exp(−βE)/Z`) |
| INV-6 | `GenerateRandomRna` is deterministic for a fixed seed; output length and alphabet (A/C/G/U) correct | Yes | Reproducibility convention |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Z no-pair | `Z("AAAA")` — no canonical pair | `Z = 1` | INV-1; Source 2 |
| M2 | Z short-span | `Z("GC")` — only pair has span ≤ m | `Z = 1` | Source 2 (min-loop) |
| M3 | Z count GGGGCCCC | `Z("GGGGCCCC", E_bp=0)` | `Z = 16` (# structures) | Evidence dataset (recurrence + brute force) |
| M4 | Z count GGGAAACCC | `Z("GGGAAACCC", E_bp=0)` | `Z = 20` (# structures) | Evidence dataset (recurrence + brute force) |
| M5 | Single-pair probability | `P[0,5]("GAAAAC", E_bp=0)` | `Z=2`, `P[0,5]=0.5` | Source 3 outside recursion |
| M6 | Probability spectrum (all 9 pairs incl. nested) | `P("GGGAAACCC", E_bp=0)` | `P[0,8]=6/20`, `P[1,7]=4/20`, `P[2,6]=6/20`, `P[1,6]=P[2,7]=P[1,8]=3/20`, `P[0,6]=P[2,8]=1/20`, `P[0,7]=3/20` | Brute-force enumeration; outside recursion |
| M6b | Probability spectrum (all 10 pairs incl. nested) | `P("GGGGCCCC", E_bp=0)` | `P[1,5]=P[2,6]=3/16`, `P[0,7]=4/16`, … (all 10) | Brute-force enumeration; outside recursion |
| M6c | Weighted probabilities (nested) | `P("GGGGCCCC", E_bp=-1)` | `Z=180.0183…`, `P[1,5]=P[2,6]=0.31334323`, … | Weighted brute-force enumeration |
| M6d | Per-base pairing ≤ 1 | `Σ_{pairs at p} P ≤ 1` for several seqs | every position ≤ 1 | McCaskill ensemble property |
| M7 | Probabilities in [0,1] | All `P[i,j]` for a folding sequence | every value in [0,1] | INV-2 |
| M8 | Boltzmann identity | `CalculateStructureProbability(−5,−5)` | `= 1.0` | INV-5; Source 1 |
| M9 | Boltzmann value | `CalculateStructureProbability(−5,−6)` | `= exp(−1/RT) = 0.197370910785…` | INV-5; Source 1, 5 |
| M10 | Null input | `CalculatePartitionFunction(null)` | throws `ArgumentNullException` | contract |
| M11 | Empty input | `CalculatePartitionFunction("")` | `Z = 1`, no pairs | INV-1; Source 2 |
| M12 | Seeded determinism | `GenerateRandomRna(50, new Random(42))` twice | identical strings, length 50, alphabet ⊆ {A,C,G,U} | INV-6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Monotonicity in E_bp | `Z(E_bp=−2) > Z(E_bp=−1) > Z(E_bp=0)` for `GGGAAACCC` | strictly increasing | INV-4 |
| S2 | Temperature validation | `CalculatePartitionFunction("GGGAAACCC", temperature:0)` | throws `ArgumentOutOfRangeException` | contract |
| S3 | GenerateRandomRna GC content | length & alphabet correct for length 0 and 100 | empty string / 100 RNA bases | INV-6 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property-based invariants | random sequences (fixed seed) up to len 30 | `Z ≥ 1` and all `P[i,j] ∈ [0,1]` | INV-1, INV-2; O(n³) property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for RNA-PARTITION-001. `CalculateStructureProbability(double,double,double)` and `GenerateRandomRna(int,double)` pre-existed in `RnaSecondaryStructure.cs` but had no dedicated tests. `CalculatePartitionFunction` and the seeded `GenerateRandomRna` overload are new in this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12, S1–S3, C1 | ❌ Missing | brand-new unit, no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PartitionFunction_Tests.cs` — all cases for this unit.
- **Remove:** none (no pre-existing tests for these methods).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|-----------|
| `RnaSecondaryStructure_PartitionFunction_Tests.cs` | canonical | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | M11 | ❌ Missing | implemented | ✅ Done |
| 12 | M12 | ❌ Missing | implemented | ✅ Done |
| 13 | S1 | ❌ Missing | implemented | ✅ Done |
| 14 | S2 | ❌ Missing | implemented | ✅ Done |
| 15 | S3 | ❌ Missing | implemented | ✅ Done |
| 16 | C1 | ❌ Missing | implemented (property-based) | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact `Z=1` |
| M2 | ✅ Covered | exact `Z=1` |
| M3 | ✅ Covered | exact `Z=16` |
| M4 | ✅ Covered | exact `Z=20` |
| M5 | ✅ Covered | exact `P=0.5` |
| M6 | ✅ Covered | exact probability spectrum |
| M7 | ✅ Covered | bounds checked |
| M8 | ✅ Covered | exact `1.0` |
| M9 | ✅ Covered | exact `exp(−1/RT)` |
| M10 | ✅ Covered | ArgumentNullException |
| M11 | ✅ Covered | `Z=1`, empty pairs |
| M12 | ✅ Covered | determinism + alphabet |
| S1 | ✅ Covered | strict monotonicity |
| S2 | ✅ Covered | ArgumentOutOfRangeException |
| S3 | ✅ Covered | length/alphabet |
| C1 | ✅ Covered | property-based invariants |

Total in-scope cases: 16; ✅ = 16.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Simplified fixed-per-pair energy model (`E_bp`) rather than full Turner 2004 loop energies; the partition-function recurrence, probability formula, and invariants remain fully conformant to McCaskill 1990 | §1.4, all Z/P expected values |

---

## 7. Open Questions / Decisions

1. Decision: expected `Z` values use `E_bp = 0` so `Z` reduces to a pure structure count, derivable independently of the implementation (recurrence + exhaustive enumeration). This makes the MUST tests evidence-based rather than implementation-echoing.
2. Decision: the energy-model simplification is documented (ASM-01 in the algorithm doc); a full Turner-parameter partition function is out of scope for this unit.
