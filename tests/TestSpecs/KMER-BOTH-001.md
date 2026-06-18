# Test Specification: KMER-BOTH-001

**Test Unit ID:** KMER-BOTH-001
**Area:** K-mer
**Algorithm:** K-mer counting over both strands (forward + reverse-complement) of double-stranded DNA
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Anvar et al. (2014), Genome Biology 15:555 (kPAL) | 1 | https://doi.org/10.1186/s13059-014-0555-3 | 2026-06-14 |
| 2 | kPAL documentation — Methodology | 3 | https://kpal.readthedocs.io/en/latest/method.html | 2026-06-14 |
| 3 | Shporer et al. (2016), Inversion symmetry of DNA k-mer counts | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5006273/ | 2026-06-14 |
| 4 | Marçais & Kingsford (2011), Bioinformatics 27(6):764–770 (Jellyfish) | 1 | https://doi.org/10.1093/bioinformatics/btr011 | 2026-06-14 |
| 5 | Mash issue #45 — canonical k-mer (contrast) | 3 | https://github.com/marbl/Mash/issues/45 | 2026-06-14 |

### 1.2 Key Evidence Points

1. Both-strand count is additive: "adding the values of each k-mer to its reverse complement" — kPAL Methodology (source 2); "a sum of k-mers and their reverse complements" — Anvar et al. 2014 (source 1).
2. Counting w on the reverse-complement strand (read 5'→3') equals counting RC(w) on the forward strand ⇒ count[w] = forward[w] + forward[RC(w)] — Shporer et al. 2016 (source 3).
3. A length-L sequence has L − k + 1 overlapping k-mers per strand ⇒ both-strand grand total = 2·(L − k + 1) — Marçais & Kingsford 2011 (source 4).
4. Reverse-complement palindromes (RC(w)=w, e.g. AT, GC, ACGT) receive both contributions on one key ⇒ count = 2·forward[w] — derived from source 3.
5. This is NOT canonical collapsing (Jellyfish -C / Mash, sources 4–5), which keeps only the lexicographically smaller representative; here every observed k-mer keeps a key and w, RC(w) carry equal counts.

### 1.3 Documented Corner Cases

- k > L (L − k + 1 ≤ 0) ⇒ no windows ⇒ empty result (source 4 occurrence definition).
- Palindromic k-mer ⇒ count doubled on a single key (source 3).
- k = L ⇒ exactly one window per strand (source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing additive both-strand counting with canonical collapsing — they give different keys/counts (sources 1–5).
2. Forgetting that palindromic k-mers double on one key rather than splitting (source 3).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CountKmersBothStrands(string sequence, int k)` | KmerAnalyzer | Canonical | Additive forward + reverse-complement count per k-mer |
| `CountKmersBothStrands(DnaSequence dna, int k)` | KmerAnalyzer | Delegate | Thin wrapper over the string overload |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | count[w] = forward[w] + forward[RC(w)] for every key w | Yes | Shporer et al. 2016 (source 3) |
| INV-2 | Σ counts = 2·(L − k + 1) when L ≥ k, else 0 | Yes | Marçais & Kingsford 2011 (source 4) |
| INV-3 | count[w] = count[RC(w)] (profile is strand-symmetric) | Yes | Anvar et al. 2014 (source 1) / source 3 |
| INV-4 | Palindromic w (RC(w)=w) ⇒ count[w] = 2·forward[w] | Yes | derived from source 3 |
| INV-5 | All counts ≥ 1 and integer | Yes | k-mer occurrence definition (source 4) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Worked example ATGGC, k=2 | Mixed palindromic + non-palindromic | {AT:2,TG:1,GG:1,GC:2,CC:1,CA:1} | sources 1,3 |
| M2 | Palindromic ACGT, k=2 | Every 2-mer is its own RC | {AC:2,CG:2,GT:2} | source 3 (RC(w)=w) |
| M3 | Non-palindromic AAA, k=2 | RC adds complementary key | {AA:2,TT:2} | source 1 (sum w+RC) |
| M4 | Grand total invariant | ATGGC, k=2 | Σ = 2·(5−2+1) = 8 | source 4 (INV-2) |
| M5 | Identity vs forward+RC | For ATGGC k=2, each key = forward[w]+forward[RC(w)] | all keys match | source 3 (INV-1) |
| M6 | Strand symmetry | count[w] = count[RC(w)] for ATGGC k=2 | holds for every key | source 1/3 (INV-3) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity | lowercase "atggc" == uppercase | equal dictionaries | sibling CountKmers upper-cases |
| S2 | k = L | "ATGC", k=4 | {ATGC:1, GCAT:1} (RC(ATGC)=GCAT) | one window per strand |
| S3 | DnaSequence overload delegates | DnaSequence("ATGGC") == string overload | equal dictionaries | delegate smoke |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty sequence | "" , k=2 | empty dictionary | window count ≤ 0 |
| C2 | Null sequence | null, k=2 | empty dictionary | null-safe |
| C3 | k > L | "AC", k=5 | empty dictionary | L − k + 1 ≤ 0 |
| C4 | k ≤ 0 | "ACGT", k=0 | throws ArgumentOutOfRangeException | API contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmers_Tests.cs` contained a `#region Both Strands (CountKmersBothStrands)` (3 tests under KMER-COUNT-001) covering this method with weak comments and partial assertions.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (ATGGC exact dict) | ❌ Missing | no full-dictionary exact-match test existed |
| M2 (ACGT palindrome) | ⚠ Weak | old test asserted only 3 keys, no completeness/total |
| M3 (AAA non-palindrome) | ⚠ Weak | old test asserted only AA/TT, no completeness |
| M4 (grand-total invariant) | ⚠ Weak | old test asserted total only, no per-key values |
| M5 (forward+RC identity) | ❌ Missing | not present |
| M6 (strand symmetry) | ❌ Missing | not present |
| S1 (case-insensitivity) | ❌ Missing | not present |
| S2 (k=L) | ❌ Missing | not present |
| S3 (DnaSequence delegate) | ❌ Missing | not present |
| C1–C4 (edge/error) | ❌ Missing | not present |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_CountKmersBothStrands_Tests.cs` — all KMER-BOTH-001 tests live here.
- **Remove:** the `#region Both Strands (CountKmersBothStrands)` (3 tests) from `KmerAnalyzer_CountKmers_Tests.cs` to eliminate duplicate coverage of this method (it belongs to KMER-BOTH-001, not KMER-COUNT-001).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `KmerAnalyzer_CountKmersBothStrands_Tests.cs` | Canonical KMER-BOTH-001 | 13 |
| `KmerAnalyzer_CountKmers_Tests.cs` | KMER-COUNT-001 (both-strand region removed) | (unchanged minus 3) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented exact-dictionary test | ✅ Done |
| 2 | M2 | ⚠ Weak | rewrote with full dict + total | ✅ Done |
| 3 | M3 | ⚠ Weak | rewrote with full dict + total | ✅ Done |
| 4 | M4 | ⚠ Weak | rewrote per-key + total invariant | ✅ Done |
| 5 | M5 | ❌ Missing | implemented forward+RC identity | ✅ Done |
| 6 | M6 | ❌ Missing | implemented strand-symmetry | ✅ Done |
| 7 | S1 | ❌ Missing | implemented case-insensitivity | ✅ Done |
| 8 | S2 | ❌ Missing | implemented k=L | ✅ Done |
| 9 | S3 | ❌ Missing | implemented DnaSequence delegate smoke | ✅ Done |
| 10 | C1 | ❌ Missing | implemented empty | ✅ Done |
| 11 | C2 | ❌ Missing | implemented null | ✅ Done |
| 12 | C3 | ❌ Missing | implemented k>L | ✅ Done |
| 13 | C4 | ❌ Missing | implemented k≤0 throws | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact dict {AT:2,TG:1,GG:1,GC:2,CC:1,CA:1} |
| M2 | ✅ | {AC:2,CG:2,GT:2} with completeness |
| M3 | ✅ | {AA:2,TT:2} with completeness |
| M4 | ✅ | Σ = 8 = 2·(5−2+1) |
| M5 | ✅ | per-key forward+RC identity |
| M6 | ✅ | strand symmetry count[w]=count[RC(w)] |
| S1 | ✅ | lowercase == uppercase |
| S2 | ✅ | k=L → {ATGC:1, GCAT:1} |
| S3 | ✅ | DnaSequence overload == string overload |
| C1 | ✅ | empty → empty |
| C2 | ✅ | null → empty |
| C3 | ✅ | k>L → empty |
| C4 | ✅ | k≤0 → throws |

Total in-scope cases: 13. ✅ count: 13.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty/short input (incl. k > L) ⇒ empty dictionary | C1, C2, C3 |
| 2 | k ≤ 0 ⇒ ArgumentOutOfRangeException (API-shape, matches sibling methods) | C4 |

Both are boundary/API-shape only; neither changes output for valid in-range input (not correctness-affecting on the algorithm's defined domain).

---

## 7. Open Questions / Decisions

1. Decision: implement the additive both-strand (kPAL "balance") semantics, NOT canonical collapsing — the method name and registry ("Forward + reverse complement"), and the existing implementation, both denote the additive view. Canonical collapsing is a separate, non-implemented variant noted in the algorithm doc §5.3.
2. Decision: added a `string` overload to match the registry signature `CountKmersBothStrands(sequence, k)`; the pre-existing `DnaSequence` overload now delegates to it.
