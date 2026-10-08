# KMER-COUNT-001: K-mer Counting Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | KMER-COUNT-001 |
| **Area** | K-mer Analysis |
| **Canonical Method** | `KmerAnalyzer.CountKmers(string, int)` |
| **Complexity** | O(n) |
| **Invariant** | Sum of counts = n − k + 1 |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `CountKmers(string, int)` | KmerAnalyzer | Canonical | Deep |
| `CountKmersSpan(ReadOnlySpan<char>, int)` | SequenceExtensions | Span variant | Deep |
| `CountKmersBothStrands(DnaSequence, int)` | KmerAnalyzer | Both strands | Deep |
| `CountKmers(DnaSequence, int)` | KmerAnalyzer | Wrapper | Smoke |
| `CountKmers(string, int, CancellationToken, IProgress<double>)` | KmerAnalyzer | Async delegate | Smoke |
| `CountKmers(string, int, KmerCountingOptions, CancellationToken, IProgress<double>)` | KmerAnalyzer | ACGT-only / canonical (Jellyfish) | Deep |
| `DistinctKmers(string, int[, KmerCountingOptions])` | KmerAnalyzer | Distinct set | Smoke |
| `GetKmerSpectrum` / `AnalyzeKmers` with `KmerCountingOptions` | KmerAnalyzer | Jellyfish histo/stats on -C | Deep |

## Evidence Sources

1. **Wikipedia — K-mer:** Definition, pseudocode, L − k + 1 formula, k-mer tables, 4^k bound
   - URL: https://en.wikipedia.org/wiki/K-mer
2. **Rosalind — K-mer Composition (KMER):** 4-mer composition sample dataset (415 bp, 209 unique 4-mers, sum = 412)
   - URL: https://rosalind.info/problems/kmer/
3. **Rosalind — Clump Finding (BA1E):** K-mer clump definition
   - URL: https://rosalind.info/problems/ba1e/
4. **Jellyfish source + executed binary 2.3.1:** `mer_iterator.hpp` (non-ACGT window reset; canonical `m_ < rcm_`), `mer_dna.hpp` (`codes[256]`), `count_main_cmdline.yaggo` (`-C`)
   - URL: https://github.com/gmarcais/Jellyfish

## Test Categories

### MUST Tests (Evidence-Backed)

| ID | Test | Evidence |
|----|------|----------|
| M1 | Empty sequence returns empty dictionary | Wikipedia pseudocode: loop 0 to L−k+1 yields nothing when L=0 |
| M2 | k > sequence length returns empty dictionary | Wikipedia: L − k + 1 becomes ≤ 0 |
| M3 | k ≤ 0 throws ArgumentOutOfRangeException (all APIs) | k must be positive for valid k-mer definition |
| M4 | null sequence returns empty dictionary | Defensive programming |
| M5 | Total count invariant: sum(counts) = L − k + 1 | Wikipedia: "a sequence of length L will have L − k + 1 k-mers" |
| M6 | Homopolymer: single k-mer, count = L − k + 1 | Wikipedia: all same bases example |
| M7 | Case-insensitive counting (all APIs) | Wikipedia/algorithm norm: k-mers case-insensitive |
| M8 | Distinct k-mers counted correctly | Rosalind KMER problem |
| M9 | Overlapping k-mers counted correctly | Wikipedia sliding window pseudocode |
| M10 | CountKmersSpan produces same results as CountKmers (including case) | API consistency |
| M11 | CountKmersBothStrands combines forward + reverse complement | DNA double-strand property |
| M12 | Wikipedia example: "ATGG" → ATG, TGG | Wikipedia lead diagram caption |
| M13 | Wikipedia table: "GTAGAGCTGT" exact k-mers for k=2 and k=4 | Wikipedia Introduction table |
| M14 | Rosalind KMER sample: specific 4-mer counts + exact unique count (209) | Rosalind sample output |
| M15 | CancellationToken overload normalizes case | API consistency with canonical method |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Mixed case input normalized (exact counts) | Robustness for real-world data |
| S2 | Non-DNA characters handled (IUPAC N, exact counts) | Genomic data often contains N |
| S3 | k = 1 counts individual nucleotides | Edge case at minimum valid k |
| S4 | k = sequence length yields single k-mer | Boundary condition |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Cancellation token stops operation | Async API contract |
| C2 | Progress reporting works | Async API contract |

## Test File

| File | Tests | Coverage |
|------|-------|----------|
| KmerAnalyzer_CountKmers_Tests.cs | 36 tests | All MUST/SHOULD tests covered |

## Coverage Classification

All tests use exact values derived from theory (Wikipedia/Rosalind), not from implementation output.

| Test Method | Classification | Notes |
|-------------|---------------|-------|
| CountKmers_EmptySequence_ReturnsEmptyDictionary | ✅ Covered | M1 |
| CountKmers_NullSequence_ReturnsEmptyDictionary | ✅ Covered | M4 |
| CountKmers_KLargerThanSequence_ReturnsEmptyDictionary | ✅ Covered | M2 |
| CountKmers_InvalidK_ThrowsArgumentOutOfRangeException | ✅ Covered | M3 (k=0, -1, -10) |
| CountKmers_KEqualSequenceLength_ReturnsSingleKmer | ✅ Covered | S4 |
| CountKmers_TotalCountInvariant_SumEqualsLMinusKPlusOne | ✅ Covered | M5 (4 cases) |
| CountKmers_TotalCountInvariant_HoldsForAllValidK | ✅ Covered | M5 (k=1..12) |
| CountKmers_SimpleSequence_CountsDistinctKmersCorrectly | ✅ Covered | M8 |
| CountKmers_Homopolymer_SingleKmerWithCorrectCount | ✅ Covered | M6 (4 cases) |
| CountKmers_OverlappingKmers_AllCounted | ✅ Covered | M9 |
| CountKmers_KEqualsOne_CountsNucleotides | ✅ Covered | S3 |
| CountKmers_LowercaseSequence_NormalizedToUppercase | ✅ Covered | M7 (exact counts) |
| CountKmers_MixedCase_TreatedAsSameKmer | ✅ Covered | S1 (all k-mers verified) |
| CountKmers_WithAmbiguousBase_CountedAsIs | ✅ Covered | S2 (exact counts) |
| CountKmersSpan_ProducesSameResultAsCountKmers | ✅ Covered | M10 |
| CountKmersSpan_EmptySpan_ReturnsEmptyDictionary | ✅ Covered | M10 edge |
| CountKmersSpan_KLargerThanSpan_ReturnsEmptyDictionary | ✅ Covered | M10 edge |
| CountKmersSpan_InvalidK_ThrowsArgumentOutOfRangeException | ✅ Covered | M3 (Span API, k=0, -1) |
| CountKmersSpan_KEqualToLength_ReturnsSingleKmer | ✅ Covered | M10 boundary |
| CountKmersBothStrands_CombinesForwardAndReverseComplement | ✅ Covered | M11 palindromic |
| CountKmersBothStrands_NonPalindromicSequence_AddsNewKmers | ✅ Covered | M11 non-palindromic |
| CountKmersBothStrands_TotalCountInvariant | ✅ Covered | M11 invariant |
| CountKmers_DnaSequence_DelegatesToStringVersion | ✅ Covered | Wrapper smoke |
| CountKmers_WikipediaExample_ATGG_TwoThreeMers | ✅ Covered | M12 |
| CountKmers_WikipediaTable_GTAGAGCTGT_TwoMers | ✅ Covered | M13 k=2 (all 7 k-mers) |
| CountKmers_WikipediaTable_GTAGAGCTGT_TotalKmersPerK | ✅ Covered | M13 k=1..10 |
| CountKmers_WikipediaTable_GTAGAGCTGT_FourMers | ✅ Covered | M13 k=4 (all 7 k-mers) |
| CountKmers_RosalindKmerSample_SpecificCounts | ✅ Covered | M14 (7 spot checks + total) |
| CountKmers_RosalindKmerSample_ExactUniqueCount | ✅ Covered | M14 (exact = 209) |
| CountKmersSpan_LowercaseInput_NormalizesToUppercase | ✅ Covered | M7/M10 regression |
| CountKmersSpan_MixedCase_MatchesCountKmers | ✅ Covered | M10 case regression |
| CountKmers_CancellationOverload_NormalizesCase | ✅ Covered | M15 |
| CountKmers_RosalindKmerSample_FullCompositionArray_AllOverloads | ✅ Covered | M14/M10/M7: full 256-value Rosalind array locked for string, cancellation, async, DnaSequence, Span and lower-case input (review 2026-09) |

## Deviations and Assumptions

- Non-ACGT symbols (e.g. `N`) are counted literally by the option-less overloads (generic string k-mer definition, Wikipedia). The Jellyfish convention, which drops every window containing a non-ACGT base, is `KmerCountingOptions.AcgtOnly`, and Jellyfish `-C` is `KmerCountingOptions.Canonical` (which implies ACGT-only).
- `CountKmers(string,…)` returns empty for null/empty input *before* validating `k`; `CountKmersSpan` validates `k` first (empty span with `k ≤ 0` throws). Documented overload asymmetry.
- `CountKmers(DnaSequence null, k)` throws `NullReferenceException` (no explicit null guard).

## Review 2026-09 (batch B06)

Stage A PASS-with-notes, Stage B PASS-with-notes. Rosalind KMER sample output (256 values) reproduced exactly by Python `collections.Counter` and Biopython `Seq.count_overlap` (sum 412, 209 non-zero, max CAGT = 8) and by every C# overload. The synchronous `CountKmers(string,int)` now delegates to the cancellation-aware overload (single counting loop).

## Audit round 1 (B06 WP1) — ACGT-only and canonical counting

Test file: `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_CountingOptions_Tests.cs` (36 cases). All values below come from the executed Jellyfish 2.3.1 binary and are not derived from the code.

| ID | Test | Expected (Jellyfish 2.3.1) |
|----|------|----------------------------|
| O1 | `-C` full `dump -c` tables, 8 inputs (incl. lower case, N, R, U, homopolymer) | e.g. ATGATG k=3 → ATC:1 ATG:2 TCA:1; `acgtNNacgtacgRtTTGCAnA` k=3 → AAA:1 ACG:5 CAA:1 GCA:2 GTA:2 |
| O2 | ACGT-only `dump -c` tables, 4 inputs | `ACGTNACGT` k=4 → ACGT:2; `AAUUAA` k=2 → AA:2 |
| O3 | `stats` + `histo` rows on plain/`-C` counts (9 rows), plus the `DistinctKmers` size | Rosalind k=4 `-C`: 23/130/412/10, histo 1:23 … 10:2 |
| O4 | Mean and entropy on `-C` tables | scipy `entropy(base=2)`: 6.779144227048732, 4.1343361131944505, 2.0403733936884962 |
| O5 | `-C` with `-L 2` | BA1B: Distinct 4, Total 11, Max 4 |
| O6 | Invariants: default options equal the literal overloads; `Canonical` implies ACGT-only; keys = ordinal min(w, RC(w)); canonical table invariant under reverse complement of the input; sum = number of all-ACGT windows | — |
| O7 | Contracts: null/empty/k>L/all-N give empty results; k ≤ 0 throws (ParamName "k"); cancellation and progress (0.0, 1.0); `DistinctKmers` returns a caller-owned ordinal set | — |

## Audit round 1 (B06 WP4) — single span-lookup loop, CountKmersSpan contract, null DnaSequence, parallel count

Test file: `Unit/Analysis/KmerAnalyzer_ParallelAndBackgroundD2_Tests.cs`.

| ID | Test | Evidence |
|----|------|----------|
| P1 | `CountKmers` (literal and ACGT-only) equals an independent naive substring counter on 300 random inputs (alphabet ACGTacgtNR, L ≤ 3000, k ≤ 8) | k-mer definition; Jellyfish window rule (WP1) |
| P2 | `KmerAnalyzer.CountKmersSpan` equals `CountKmers` for every k incl. k > L on mixed-case/IUPAC input | same contract |
| P3 | Empty span with k ≤ 0 → empty (as `CountKmers`); non-empty with k ≤ 0 → `ArgumentOutOfRangeException("k")` | `ValidateKmerLength` contract (F2) |
| P4 | `CountKmers(DnaSequence)`, `CountKmers(DnaSequence, k, ct)`, `CountKmersBothStrands(DnaSequence)` with null → `ArgumentNullException("dna")` | .NET Framework Design Guidelines, argument validation |
| P5 | `CountKmersParallel` equals the serial count on 8 random inputs × 3 option sets × degrees 2/3/4/8/−1 (140k–420k windows) | counting is a sum over windows (any partition) |

