# Test Specification: KMER-GENERATE-001

**Test Unit ID:** KMER-GENERATE-001
**Area:** K-mer
**Algorithm:** K-mer Generation (enumerate all possible k-mers over an alphabet)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — K-mer | 4 | https://en.wikipedia.org/wiki/K-mer | 2026-06-14 |
| 2 | BioInfoLogics — k-mer counting, part I (Clavijo, 2018) | 3 | https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/ | 2026-06-14 |
| 3 | Python Standard Library — itertools.product | 3 | https://docs.python.org/3/library/itertools.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. A k-mer is a length-k string over the sequence alphabet ("substrings of length k contained within a biological sequence"). — Source 1.
2. The number of *all possible* k-mers over an alphabet of n monomers is **n^k**; for DNA {A,C,G,T} this is **4^k**. "there exist n^k total possible k-mers, where n is number of possible monomers" — Source 1; "the possible combinations of k positions are computed as 4^k" — Source 2.
3. Generating all k-mers is the k-fold Cartesian product of the alphabet: `product(A, repeat=k)`. — Source 3.
4. If the alphabet is supplied sorted, the products are emitted in **lexicographic order** with the rightmost position advancing fastest ("if the input's iterables are sorted, the product tuples are emitted in sorted order"; "product(range(2), repeat=3) → 000 001 010 011 100 101 110 111"). — Source 3. {A,C,G,T} is already sorted ⇒ DNA k-mers come out AAA, AAC, …, TTT.

### 1.3 Documented Corner Cases

- **k ≤ 0:** no valid k-mer length (definition: substrings *of length k*). — Source 1.
- **Single-letter alphabet:** exactly one k-mer (the homopolymer), since 1^k = 1. — Source 3.
- **Unsorted alphabet:** still yields all n^k strings, but in the alphabet's positional order (lexicographic only when alphabet sorted). — Source 3.

### 1.4 Known Failure Modes / Pitfalls

1. Emitting fewer/more than n^k strings, or duplicates — violates the Cartesian-product/n^k count. — Sources 1–3.
2. Wrong position varying fastest — would break the documented lexicographic order on a sorted alphabet. — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GenerateAllKmers(int k, string alphabet = "ACGT")` | KmerAnalyzer | **Canonical** | Returns all n^k k-mers; lexicographic when alphabet sorted |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Output count = alphabet.Length^k (4^k for default DNA alphabet) | Yes | Source 1 (n^k), Source 2 (4^k) |
| INV-2 | All emitted k-mers are distinct (the output is a set of size n^k) | Yes | Source 3 (each k-mer is a unique length-k tuple) |
| INV-3 | Every emitted k-mer has length exactly k and uses only alphabet characters | Yes | Source 1 (length-k string over alphabet) |
| INV-4 | For a sorted alphabet, emission order is lexicographic (rightmost position fastest) | Yes | Source 3 (odometer ordering) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | k=1 default DNA | All monomers | Exactly [A, C, G, T] in order (count 4^1=4) | Source 1 (DNA monomers, n^k) |
| M2 | k=2 default DNA | All 2-mers, lexicographic | [AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT] (count 16) | Source 1 + Source 3 |
| M3 | k=3 default DNA | Count and boundaries | 64 k-mers; first AAA, second AAC, last TTT | Source 1 (4^3); Source 3 (ordering) |
| M4 | count = 4^k for k=1..6 | Cardinality of universe | counts 4,16,64,256,1024,4096 | Source 1/2 (4^k) |
| M5 | custom 20-letter protein alphabet, k=2 | n^k generalisation | count = 400 | Source 1 (n^k generalised) |
| M6 | no duplicates / is the Cartesian set | Distinctness | distinct count = 4^k; set equals product set | Source 3 (unique tuples) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | single-letter alphabet "A", k=4 | 1^k degenerate | Exactly [AAAA] (count 1) | Source 3 |
| S2 | every k-mer length = k, alphabet chars only | INV-3 | all length-2 over {A,C,G,T} | Source 1 |
| S3 | k ≤ 0 (k=0 and k=-1) | invalid length | throws ArgumentOutOfRangeException | definition |
| S4 | null / empty alphabet | no alphabet | throws ArgumentException | derived from contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | unsorted alphabet "TGCA", k=1 | ordering follows alphabet | [T,G,C,A] (not sorted) | Source 3 (order depends on alphabet) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file targets `GenerateAllKmers`. `KmerAnalyzerTests.cs` and the per-method `KmerAnalyzer_*_Tests.cs` files cover other methods (CountKmers, Frequency, Find, Distance, Unique/MinCount, Async). `GenerateAllKmers` is implemented in `KmerAnalyzer.cs` (recursive Cartesian product, prefix built left-to-right) but untested.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_GenerateAllKmers_Tests.cs` — all GenerateAllKmers tests.
- **Remove:** none (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `KmerAnalyzer_GenerateAllKmers_Tests.cs` | Canonical for KMER-GENERATE-001 | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact ordered-list assertion | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact 16-element ordered list | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented count + first/second/last | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented 4^k for k=1..6 | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented 20^2 = 400 | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented distinctness/set-equality | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented 1^4 = {AAAA} | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented length/alphabet membership | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented k=0 and k=-1 throw | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented null/empty alphabet throw | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented unsorted-alphabet ordering | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact ordered list [A,C,G,T] |
| M2 | ✅ | Exact 16-element lexicographic list |
| M3 | ✅ | Count 64; first AAA, second AAC, last TTT |
| M4 | ✅ | 4^k for k=1..6 |
| M5 | ✅ | 20-letter alphabet → 400 |
| M6 | ✅ | distinct count = 4^k; set equality |
| S1 | ✅ | "A", k=4 → {AAAA} |
| S2 | ✅ | all length-k, alphabet-only |
| S3 | ✅ | k=0, k=-1 throw ArgumentOutOfRangeException |
| S4 | ✅ | null/empty alphabet throw ArgumentException |
| C1 | ✅ | "TGCA", k=1 → [T,G,C,A] |

Total in-scope cases: 11. ✅ count: 11.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default alphabet "ACGT" (upper case); lexicographic order guaranteed only when supplied alphabet is sorted (itertools.product property) | Default-parameter behaviour; M1–M4, C1 |

---

## 7. Open Questions / Decisions

1. None. Behaviour (n^k universe, Cartesian-product enumeration, sorted-alphabet lexicographic order) is fully determined by the retrieved sources.
