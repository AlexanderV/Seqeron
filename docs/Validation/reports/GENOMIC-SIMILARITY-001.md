# Validation Report: GENOMIC-SIMILARITY-001 — Sequence Similarity (k-mer Jaccard index)

- **Validated:** 2026-06-16   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int kmerSize = 5)` (private helper `GenomicAnalyzer.GetKmers(string, int)`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Jaccard index — Wikipedia** (citing Jaccard 1901 primary). WebFetched https://en.wikipedia.org/wiki/Jaccard_index this session. Confirms verbatim:
   - Formula: `J(A,B) = |A∩B| / |A∪B| = |A∩B| / (|A|+|B|−|A∩B|)`.
   - Range: "By definition, 0 ≤ J(A,B) ≤ 1".
   - Disjoint: "If the sets A and B have no elements in common … |A∩B|=0 and therefore J(A,B)=0".
   - Identical: "the two sets are equal. In that case A∩B=A∪B=A=B, so then J(A,B)=1".
   - Undefined case: "The definition is not well-defined when μ(A∪B)=0 or μ(A∪B)=∞" (set form does not assign a value to empty union).
   - Jaccard distance: `d_J(A,B) = 1 − J(A,B) = (|A∪B| − |A∩B|) / |A∪B|`.
2. **Ondov et al. 2016, Mash (Genome Biology 17:132)** — WebFetched https://pmc.ncbi.nlm.nih.gov/articles/PMC4915045/ this session. Confirms verbatim the k-mer-set application: "The Jaccard index is simply the fraction of shared hashes … out of all distinct hashes in A and B," and the MinHash sketch estimate `j = x/s'`. Establishes that k-mers are compared as **distinct sets** (hashes), repeats collapse.

### Formula check

The description (`Sequence_Similarity.md` §2.2, Evidence, TestSpec §1.2) states `J(A,B) = |A∩B| / |A∪B| = |A∩B| / (|A|+|B|−|A∩B|)`, scaled ×100. This is verbatim identical to the Wikipedia/Jaccard-1901 formula and the Ondov et al. k-mer-set definition. Range [0,1]→[0,100], identical→100, disjoint→0 all match the sourced extremes.

### Edge-case semantics check

- **Empty union** (both sequences empty, or both shorter than k): the sources declare J undefined for empty union (Wikipedia: "not well-defined when μ(A∪B)=0"). The unit therefore documents the `0.0` return as an explicit **implementation convention (ASM-1)**, not a sourced value. This is the correct, honest treatment — no source is misquoted as mandating 0.
- **One side empty:** `|A∩B|=0` over a non-empty union → J=0 → 0.0. This is a genuine sourced value (disjoint case), not the undefined case.
- **Distinct-set semantics:** repeats within a sequence counted once — directly from Ondov et al. ("distinct hashes").
- **×100 scaling** — flagged (ASM-2) as a presentation convention; does not alter ordering.

### Independent cross-check (numbers)

Independent Python `set`-based reimplementation built from the sourced formula (NOT from the C# code) reproduced every expected value this session:

| seq1 | seq2 | k | \|A∩B\| | \|A∪B\| | J×100 |
|------|------|---|---------|---------|-------|
| ACGTACGT | ACGTACGA | 3 | 4 | 5 | 80.0 |
| ACGTACGT | ACGTACGT | 3 | 4 | 4 | 100.0 |
| AAAAA | CCCCC | 3 | 0 | 2 | 0.0 |
| ACGT | ACGA | 3 | 1 | 3 | 33.33333333333333 (100/3) |
| AAAAAA | AAAA | 3 | 1 | 1 | 100.0 (set semantics; 4 vs 2 k-mer counts collapse) |
| ACGTAC | "" | 3 | 0 | 4 | 0.0 |
| "" | "" | 3 | 0 | 0 | 0.0 (empty union convention) |
| AC | GT | 3 | 0 | 0 | 0.0 (both shorter than k) |

All match the formula exactly. `100/3` confirmed = 33.333333333333336 (double).

### Findings / divergences

None substantive. Two flagged assumptions (empty-union→0; ×100 scaling) are correctly documented as conventions rather than sourced literature values. Stage A **PASS**.

## Stage B — Implementation

### Code path reviewed

- `GenomicAnalyzer.CalculateSimilarity` — `GenomicAnalyzer.cs:358–384`.
- `GenomicAnalyzer.GetKmers` — `GenomicAnalyzer.cs:502–510`.
- `DnaSequence` constructor — `DnaSequence.cs:22–33` (empty allowed → empty string; non-empty normalized via `ToUpperInvariant` + alphabet validation).

### Formula realised correctly? (evidence)

- `GetKmers` builds a `HashSet<string>` over the sliding window `for (i = 0; i <= sequence.Length - k; i++)` → distinct-set semantics; when `length < k` the loop body never runs → empty set (matches the empty-k-mer-set edge cases). Correct.
- `intersection = kmers1.Count(kmers2.Contains)` = `|A∩B|`; `union = |A| + |B| − |A∩B|` (inclusion–exclusion) = `|A∪B|`. Correct, single pass.
- `return union == 0 ? 0.0 : (double)intersection / union * 100.0` — exactly `J×100`, with the documented empty-union convention. `(double)` cast prevents integer division; no overflow on stated ranges. Correct.

### Cross-verification table recomputed vs code

The full unfiltered test suite (which exercises all eight rows above through the live code) passes; the asserted values are identical to the independent Python table. Code output ≡ sourced values.

### Variant/delegate consistency

Single public overload; no `*Fast` variant or instance delegate. Nothing to cross-check.

### Test quality audit (HARD gate)

File: `GenomicAnalyzer_CalculateSimilarity_Tests.cs` (13 tests).

- **Sourced, not echoes:** M1=80.0, M2=100.0, M3=0.0, M4=100/3, M5=100.0 are taken from the Jaccard formula / Ondov set semantics, verified by independent computation. M4 is a non-round irrational-in-decimal value (a rounding or wrong-k impl fails it). M5 discriminates **set vs multiset** (a bag impl would not return 100 for AAAAAA vs AAAA). M1's own comment notes an identity-only impl would return 100 — i.e. the test would fail a deliberately-wrong implementation.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` everywhere a value is known; the only `Is.InRange(0,100)` (C1) is a legitimate *invariant* over varied inputs (INV-1), not a substitute for a known exact value. No skipped/ignored/commented tests; `1e-10` tolerance is tight and justified (needed for 100/3; harmless for representable 80.0/100.0).
- **Coverage:** the single public method's happy path + all Stage-A branches and edges — partial overlap (M1), identical (M2), disjoint (M3), non-integer fraction (M4), set semantics (M5), symmetry (M6), both empty (S1), both shorter than k (S2), null seq1/seq2 (S3/S4), kmerSize=0 (S5), range invariant (C1), one side empty (C2). Both validation branches and the empty-union branch are exercised.
- **Honest green:** full unfiltered suite **6619 passed, 0 failed, 0 skipped** (the one `MFE_Benchmark_AllScenarios` skip is unrelated/pre-existing and outside this unit); `dotnet build` 0 warnings, 0 errors.

Test-quality gate: **PASS**.

### Findings / defects

None. No code or test change was required.

## Verdict & follow-ups

- **Stage A: PASS.** Formula, range, extremes, and distinct-set semantics independently confirmed against Wikipedia/Jaccard-1901 and Ondov et al. 2016 (both fetched this session). Empty-union→0 and ×100 are correctly documented as conventions.
- **Stage B: PASS.** Code realises `J×100` over distinct k-mer `HashSet`s exactly; tests assert sourced exact values and cover every branch/edge.
- **End-state: ✅ CLEAN** — no defect found; algorithm fully functional. No files changed except validation docs.
- **Test-quality gate: PASS.**
- **Full suite:** 6619 passed / 0 failed; build 0 errors, 0 warnings.
