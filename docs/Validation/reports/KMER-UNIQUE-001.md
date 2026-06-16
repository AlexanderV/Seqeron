# Validation Report: KMER-UNIQUE-001 — Unique K-mers / K-mers with Minimum Count

- **Validated:** 2026-06-16   **Area:** K-mer
- **Canonical method(s):** `KmerAnalyzer.FindUniqueKmers(string, int)`, `KmerAnalyzer.FindKmersWithMinCount(string, int, int)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened (retrieved this session)

1. **Wikipedia — K-mer** (`https://en.wikipedia.org/wiki/K-mer`, WebFetch 2026-06-16). Confirms verbatim:
   - "In bioinformatics, k-mers are substrings of length k contained within a biological sequence."
   - "a sequence of length L will have L − k + 1 k-mers."
   - "there exist n^k total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)."
   - AGAT worked example: k=1 → A,G,A,T; k=2 → AG,GA,AT (3); k=3 → AGA,GAT (2); k=4 → AGAT (1). Overlapping, step-1 window.
2. **BioInfoLogics — k-mer counting, part I** (`https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/`, WebFetch 2026-06-16). Confirms verbatim:
   - "Distinct k-mers are counted only once, even if they appear more times."
   - "Unique k-mers are those that appear only once." (frequency exactly 1)
   - Worked table for **ATCGATCAC, k=3**: ATC total=2 (unique=0); TCG, CGA, GAT, TCA, CAC each total=1 (unique=1). **Sum: 7 total, 6 distinct, 5 unique.**
3. **Compeau & Pevzner — Bioinformatics Algorithms**, via Rosalind BA1A (`https://rosalind.info/problems/ba1a/`, WebSearch 2026-06-16). Confirms `Count(Text, Pattern)` = number of times `Pattern` appears as a substring of `Text`, **overlapping occurrences counted**: Count(CGATATATCCATAG, ATA) = 3 (not 2). "Most frequent k-mer" maximizes Count → basis for count-descending ranking and the `Count ≥ t` min-count filter.
4. **Independent terminology cross-check** (KAT/Jellyfish/KMC literature, WebSearch 2026-06-16): distinct = different sequences present; unique = singletons = appear exactly once. phi-X 9-mer spectrum: 4972 singletons (unique) out of 5170 distinct (4972+189+8+1). Independently confirms unique ⊊ distinct.

### Formula check
- Total overlapping k-mers = **L − k + 1** (Source 1). ✓
- Unique = `{P : Count(Text,P) = 1}` (Source 2). ✓
- Min-count = `{(P, Count) : Count ≥ t}`, recurrent k-mers (Source 3). ✓
- Possible k-mer universe = n^k (Source 1; relevant to sibling generator, not these two methods).

### Edge-case semantics check (all sourced)
- k > L ⇒ L − k + 1 ≤ 0 ⇒ zero k-mers (Source 1). ✓
- Empty sequence ⇒ no k-mers. ✓
- k ≤ 0 ⇒ invalid (k-mer length must be positive, "substrings of length k"). ✓
- Homopolymer (single distinct k-mer with count > 1) ⇒ zero unique k-mers (Source 2). ✓
- minCount ≤ 1 ⇒ predicate `Count ≥ t` satisfied by every observed k-mer ⇒ all distinct k-mers (consistent extension; Assumption, not contradicted by any source). ✓

### Independent cross-check (hand computation, this session)
- **ATCGATCAC, k=3** (7 windows): ATC,TCG,CGA,GAT,ATC,TCA,CAC → counts ATC=2, others=1 → 6 distinct, 5 unique {TCG,CGA,GAT,TCA,CAC}. Matches BioInfoLogics table exactly.
- **AGAT, k=2**: AG,GA,AT all count 1 → unique = {AG,GA,AT}. Matches Wikipedia.
- **ACGTACGT, k=4** (5 windows): ACGT,CGTA,GTAC,TACG,ACGT → ACGT=2, others=1. So minCount=2 → {(ACGT,2)}; minCount=1 → 4 distinct, ACGT first. Verified by hand + Python Counter.
- **AAAACGTAAA, k=2** (S2 ordering case): AA=5, AC/CG/GT/TA=1 — confirms a genuine multi-count ordering boundary (5→1).

### Findings / divergences
None. All definitions, the L−k+1 formula, the distinct-vs-unique distinction, overlapping `Count`, and every worked value in the TestSpec/Evidence trace to the external sources retrieved this session. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`:
- `CountKmers(string,int)` (L20–42): null/empty → empty dict; `k ≤ 0` → `ArgumentOutOfRangeException`; `k > Length` → empty dict; `ToUpperInvariant`; overlapping step-1 window `i = 0 .. Length-k`, `Substring(i,k)`, `TryAdd`/increment. Exactly L−k+1 windows. ✓
- `FindUniqueKmers` (L253–257): `CountKmers` then `Where(value == UniqueKmerCount /*=1*/).Select(key)`. Realises `{P : count = 1}`. ✓
- `FindKmersWithMinCount` (L274–282): `CountKmers` then `Where(value >= minCount).Select((key,value)).OrderByDescending(value)`. Realises `{(P,c) : c ≥ t}` ordered count-descending. ✓

### Formula realised correctly?
Yes. Both methods delegate to the single overlapping-count pass and apply the exact validated predicates. Counts are exact (Dictionary increments), not approximate. Invalid-k throws via `CountKmers`; empty/k>L return empty (the `Where` over an empty dict yields empty). Ordering uses `OrderByDescending` (INV-4).

### Cross-verification table recomputed vs code (via tests, exact)
| Input | Method | Expected (sourced) | Code result |
|-------|--------|--------------------|-------------|
| ATCGATCAC k=3 | FindUniqueKmers | {TCG,CGA,GAT,TCA,CAC}, ATC excluded | ✓ |
| AGAT k=2 | FindUniqueKmers | {AG,GA,AT} | ✓ |
| AAAAA k=3 | FindUniqueKmers | ∅ | ✓ |
| AGAT k=1 | FindUniqueKmers | {G,T} (A count 2) | ✓ |
| ACGTACGT k=4 min=2 | FindKmersWithMinCount | {(ACGT,2)} | ✓ |
| ACGTACGT k=4 min=1 | FindKmersWithMinCount | 4 distinct, ACGT(2) first, CGTA/GTAC/TACG(1) | ✓ |
| ACGTACGT k=4 min=3 | FindKmersWithMinCount | ∅ | ✓ |
| empty / k>L | both | ∅ | ✓ |
| k ≤ 0 | both | ArgumentOutOfRangeException | ✓ |

### Variant/delegate consistency
Both methods share the same `CountKmers` source of truth; `FindKmersWithMinCount(seq,k,1)` keys == `CountKmers` keys (M12/INV-5) and unique set == count-1 keys of the map (M11/INV-1), cross-checked against `CountKmers` independently in the tests.

### Test quality audit (canonical file `KmerAnalyzer_FindUniqueAndMinCount_Tests.cs`, 15 tests)
- **Sourced exact values, not code echoes:** M1/M2/S1 use `Is.EquivalentTo` exact sets; M4/M5 assert exact members and exact counts (2 and 1); C1 asserts exact `{G,T}`. These would FAIL a deliberately-wrong implementation (e.g. one returning distinct instead of unique, or non-overlapping counts).
- **No green-washing:** exact equality/membership and exact counts used wherever a value is known; `ArgumentOutOfRangeException` asserted by exact type (M9/M10); no skips, no widened tolerances, no comment-outs.
- **Coverage:** both public methods; all Stage-A branches — happy path, repeated-k-mer exclusion, homopolymer, empty, k>L, k≤0, case-normalization (S1), count-descending ordering on a real multi-count case (S2, AA=5 then 1s — non-vacuous), threshold-above-max (S3), monomers k=1 (C1), and two independent cross-checks vs `CountKmers` (M11/INV-1, M12/INV-5).
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6607`; `dotnet build` 0 errors. The 4 build warnings (NUnit2007) are in an unrelated file (`ApproximateMatcher_EditDistance_Tests.cs`), not in this unit's files.

### Findings / defects
- **No defect.** The legacy weak `FindKmersWithMinCount` tests were correctly removed from `KmerAnalyzerTests.cs` (only a redirect comment remains).
- **Observation (out of scope, not a defect for this unit):** the sibling KMER-FIND-001 file `KmerAnalyzer_Find_Tests.cs` also tests `FindUniqueKmers` and uses `Does.Contain`/`Does.Not.Contain` (subset assertions) for the ACGTACGT case. That is a property of the KMER-FIND-001 unit; KMER-UNIQUE-001's own canonical coverage of the same scenario uses exact `Is.EquivalentTo`, so this unit's assertions are not weakened. Noted for the KMER-FIND-001 owner; no action taken here.

## Verdict & follow-ups
- **Stage A: PASS** — description matches authoritative sources retrieved this session; every non-trivial value independently reproduced by hand.
- **Stage B: PASS** — code faithfully realises the validated formulas and edge semantics; tests assert exact sourced values and cover all branches.
- **Test-quality gate: PASS** — sourced (not code-echoed) exact expectations, no green-washing, full public surface + all Stage-A edge/error cases covered, honest green (full suite 6607 passed / 0 failed).
- **End-state: ✅ CLEAN** — no defect found; no code or test changes required.
