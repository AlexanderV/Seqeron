# Test Specification: PAT-IUPAC-001

**Test Unit ID:** PAT-IUPAC-001  
**Area:** Pattern Matching  
**Algorithm:** IUPAC Degenerate Motif Matching  
**Status:** ☑ Complete  
**Owner:** Algorithm QA Architect  
**Last Updated:** 2026-01-22  

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Nucleic acid notation | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-01-22 |
| Bioinformatics.org: IUPAC codes | https://www.bioinformatics.org/sms/iupac.html | 2026-01-22 |
| IUPAC-IUB Commission (1970) | Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents. Biochemistry 9(20):4022–4027 | Reference |
| NC-IUB (1984) | Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences. NAR 13(9):3021–3030 | Reference |

### 1.2 Algorithm Description

#### IUPAC Nucleotide Codes (Wikipedia, Bioinformatics.org)

The IUPAC notation includes eleven "ambiguity" or "degenerate" characters representing combinations of the four DNA bases. These were designed to encode positional variations for reporting DNA sequencing errors, consensus sequences, or single-nucleotide polymorphisms.

**Standard Codes (4 bases):**
- A = Adenine
- C = Cytosine
- G = Guanine
- T = Thymine (U = Uracil for RNA)

**Two-base Ambiguity Codes (6 codes):**
| Code | Mnemonic | Represents | Source |
|------|----------|------------|--------|
| R | puRine | A or G | Wikipedia, Bioinformatics.org |
| Y | pYrimidine | C or T | Wikipedia, Bioinformatics.org |
| S | Strong | G or C | Wikipedia, Bioinformatics.org |
| W | Weak | A or T | Wikipedia, Bioinformatics.org |
| K | Keto | G or T | Wikipedia, Bioinformatics.org |
| M | aMino | A or C | Wikipedia, Bioinformatics.org |

**Three-base Ambiguity Codes (4 codes):**
| Code | Mnemonic | Represents | Source |
|------|----------|------------|--------|
| B | not A | C or G or T | Wikipedia, Bioinformatics.org |
| D | not C | A or G or T | Wikipedia, Bioinformatics.org |
| H | not G | A or C or T | Wikipedia, Bioinformatics.org |
| V | not T | A or C or G | Wikipedia, Bioinformatics.org |

**Four-base Ambiguity Code (1 code):**
| Code | Mnemonic | Represents | Source |
|------|----------|------------|--------|
| N | aNy | A or C or G or T | Wikipedia, Bioinformatics.org |

#### Degenerate Pattern Matching

Degenerate pattern matching extends exact pattern matching to handle IUPAC ambiguity codes in the pattern. At each position, the sequence nucleotide must be one of the bases represented by the pattern's IUPAC code.

**Algorithm (Brute Force):**
1. For each position i in sequence where pattern fits:
2. For each position j in pattern:
   - If seq[i+j] matches IUPAC code pattern[j], continue
   - Else, break (no match at position i)
3. If all positions match, report position i

**Complexity:** O(n × m) where n = sequence length, m = pattern length

### 1.3 Reference Examples from Evidence

#### Wikipedia IUPAC Table Examples

| Nucleotide | Matches IUPAC Code | Source |
|------------|-------------------|--------|
| A | A, R, W, M, D, H, V, N | Wikipedia table |
| C | C, Y, S, M, B, H, V, N | Wikipedia table |
| G | G, R, S, K, B, D, V, N | Wikipedia table |
| T | T, Y, W, K, B, D, H, N | Wikipedia table |

#### Bioinformatics.org Examples

From the IUPAC codes reference:
- R = A or G (puRine)
- Y = C or T (pYrimidine)
- N = any base

### 1.4 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Standard base (A,C,G,T) in pattern | Exact match only | IUPAC 1970 |
| N in pattern | Matches any A,C,G,T | IUPAC 1970 |
| All ambiguity codes | Match according to table | IUPAC 1970 |
| Mixed standard + IUPAC | Each position matched independently | Definition |
| Empty pattern | Return empty/no matches | Implementation |
| Empty sequence | Return empty | Standard |
| Pattern longer than sequence | Return empty | Standard |
| Null input | ArgumentNullException | Implementation contract |
| Lowercase input | Case-insensitive matching | Implementation (ASSUMPTION) |
| Unknown character in pattern | Exact match or no match | Implementation (ASSUMPTION) |

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDegenerateMotif(DnaSequence, string)` | MotifFinder | **Canonical** | Pattern matching |
| `FindDegenerateMotif(DnaSequence, string, CancellationToken)` | MotifFinder | Variant | Cancellable |
| `FindDegenerateMotif(string, string, CancellationToken)` | MotifFinder | Variant | String API |
| `MatchesIupac(char, char)` | IupacHelper | **Canonical** | IUPAC code matching |

---

## 3. Invariants

| ID | Invariant | Verifiable |
|----|-----------|------------|
| INV-1 | Standard base codes match only themselves | Yes |
| INV-2 | N matches all four standard bases (A, C, G, T) | Yes |
| INV-3 | Each ambiguity code matches exactly the specified bases (IUPAC table) | Yes |
| INV-4 | Each ambiguity code does NOT match excluded bases | Yes |
| INV-5 | Result positions are in range [0, seq.Length - pattern.Length] | Yes |
| INV-6 | For all results r: all positions in r.MatchedSequence satisfy IUPAC pattern | Yes |
| INV-7 | Case-insensitive matching (both sequence and pattern normalized) | Yes |
| INV-8 | MatchesIupac(n, 'N') = true for all n ∈ {A, C, G, T} | Yes |
| INV-9 | MatchesIupac is symmetric for standard bases | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

#### MatchesIupac Tests (IupacHelper)

| ID | Test Case | Input (nucleotide, code) | Expected | Evidence |
|----|-----------|--------------------------|----------|----------|
| M1 | Standard A matches A | 'A', 'A' | true | IUPAC |
| M2 | Standard C matches C | 'C', 'C' | true | IUPAC |
| M3 | Standard G matches G | 'G', 'G' | true | IUPAC |
| M4 | Standard T matches T | 'T', 'T' | true | IUPAC |
| M5 | A does not match T | 'A', 'T' | false | IUPAC |
| M6 | N matches A | 'A', 'N' | true | IUPAC |
| M7 | N matches C | 'C', 'N' | true | IUPAC |
| M8 | N matches G | 'G', 'N' | true | IUPAC |
| M9 | N matches T | 'T', 'N' | true | IUPAC |
| M10 | R (purine) matches A | 'A', 'R' | true | Wikipedia |
| M11 | R (purine) matches G | 'G', 'R' | true | Wikipedia |
| M12 | R (purine) does NOT match C | 'C', 'R' | false | Wikipedia |
| M13 | R (purine) does NOT match T | 'T', 'R' | false | Wikipedia |
| M14 | Y (pyrimidine) matches C | 'C', 'Y' | true | Wikipedia |
| M15 | Y (pyrimidine) matches T | 'T', 'Y' | true | Wikipedia |
| M16 | Y (pyrimidine) does NOT match A | 'A', 'Y' | false | Wikipedia |
| M17 | Y (pyrimidine) does NOT match G | 'G', 'Y' | false | Wikipedia |
| M18 | S (strong) matches G | 'G', 'S' | true | Wikipedia |
| M19 | S (strong) matches C | 'C', 'S' | true | Wikipedia |
| M20 | S (strong) does NOT match A | 'A', 'S' | false | Wikipedia |
| M21 | S (strong) does NOT match T | 'T', 'S' | false | Wikipedia |
| M22 | W (weak) matches A | 'A', 'W' | true | Wikipedia |
| M23 | W (weak) matches T | 'T', 'W' | true | Wikipedia |
| M24 | W (weak) does NOT match G | 'G', 'W' | false | Wikipedia |
| M25 | W (weak) does NOT match C | 'C', 'W' | false | Wikipedia |
| M26 | K (keto) matches G | 'G', 'K' | true | Wikipedia |
| M27 | K (keto) matches T | 'T', 'K' | true | Wikipedia |
| M28 | K (keto) does NOT match A | 'A', 'K' | false | Wikipedia |
| M29 | K (keto) does NOT match C | 'C', 'K' | false | Wikipedia |
| M30 | M (amino) matches A | 'A', 'M' | true | Wikipedia |
| M31 | M (amino) matches C | 'C', 'M' | true | Wikipedia |
| M32 | M (amino) does NOT match G | 'G', 'M' | false | Wikipedia |
| M33 | M (amino) does NOT match T | 'T', 'M' | false | Wikipedia |
| M34 | B (not A) matches C | 'C', 'B' | true | Wikipedia |
| M35 | B (not A) matches G | 'G', 'B' | true | Wikipedia |
| M36 | B (not A) matches T | 'T', 'B' | true | Wikipedia |
| M37 | B (not A) does NOT match A | 'A', 'B' | false | Wikipedia |
| M38 | D (not C) matches A | 'A', 'D' | true | Wikipedia |
| M39 | D (not C) matches G | 'G', 'D' | true | Wikipedia |
| M40 | D (not C) matches T | 'T', 'D' | true | Wikipedia |
| M41 | D (not C) does NOT match C | 'C', 'D' | false | Wikipedia |
| M42 | H (not G) matches A | 'A', 'H' | true | Wikipedia |
| M43 | H (not G) matches C | 'C', 'H' | true | Wikipedia |
| M44 | H (not G) matches T | 'T', 'H' | true | Wikipedia |
| M45 | H (not G) does NOT match G | 'G', 'H' | false | Wikipedia |
| M46 | V (not T) matches A | 'A', 'V' | true | Wikipedia |
| M47 | V (not T) matches C | 'C', 'V' | true | Wikipedia |
| M48 | V (not T) matches G | 'G', 'V' | true | Wikipedia |
| M49 | V (not T) does NOT match T | 'T', 'V' | false | Wikipedia |

#### FindDegenerateMotif Tests

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M50 | Purine R matches A and G | seq="ATGC", motif="R" | positions [0, 2] | Wikipedia |
| M51 | Pyrimidine Y matches C and T | seq="ATGC", motif="Y" | positions [1, 3] | Wikipedia |
| M52 | Any N matches all | seq="ACGT", motif="N" | positions [0, 1, 2, 3] | Wikipedia |
| M53 | Mixed pattern RTG | seq="ATGCATGC", motif="RTG" | positions [0, 4] | Wikipedia |
| M54 | E-box CANNTG | seq="CAGCTG", motif="CANNTG" | position [0] | Biology |
| M55 | No match returns empty | seq="AAAA", motif="GGG" | [] | Standard |
| M56 | Empty pattern returns empty | seq="ACGT", motif="" | [] | Standard |
| M57 | Empty sequence returns empty | seq="", motif="ATG" | [] | Standard |
| M58 | Pattern longer than sequence | seq="AC", motif="ACGT" | [] | Standard |
| M59 | Null sequence throws | null, "ATG" | ArgumentNullException | Contract |
| M60 | Case insensitive matching | seq="ATGC", motif="atgc" | position [0] | Implementation |
| M61 | Result contains matched sequence | seq="CAGCTG", motif="CANNTG" | MatchedSequence="CAGCTG" | API |
| M62 | Result positions are valid | any valid | All in [0, len-patLen] | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | W (weak) pattern finds AT alternation | seq="ATATAT", motif="WWW" | 4 positions | Wikipedia |
| S2 | S (strong) pattern finds GC regions | seq="GCGCGC", motif="SSS" | 4 positions | Wikipedia |
| S3 | K (keto) pattern finds GT alternation | seq="GTGTGT", motif="KK" | 5 positions | Wikipedia |
| S4 | M (amino) pattern finds AC alternation | seq="ACACAC", motif="MM" | 5 positions | Wikipedia |
| S5 | Pattern at end of sequence | seq="ATGCATG", motif="ATG" | [0, 4] | Standard |
| S6 | Overlapping matches | seq="ARAR", motif="ARA" | position [0] | Standard |
| S7 | Unknown char in pattern matches exactly | seq="XYZX", motif="X" | positions [0, 3] | ASSUMPTION |
| S8 | All B positions (not A) | seq="CGT", motif="BBB" | position [0] | Wikipedia |
| S9 | All D positions (not C) | seq="AGT", motif="DDD" | position [0] | Wikipedia |
| S10 | All H positions (not G) | seq="ACT", motif="HHH" | position [0] | Wikipedia |
| S11 | All V positions (not T) | seq="ACG", motif="VVV" | position [0] | Wikipedia |
| S12 | Cancellation token works | long sequence | Completes or cancels | API |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| C1 | Large sequence performance | 10000+ chars | Completes in time | Performance |
| C2 | Restriction site pattern | GAATTC, degenerate | Correct positions | Bioinformatics |

---

## 5. ASSUMPTIONS

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Case-insensitive comparison | DNA sequences conventionally use case-insensitive matching |
| A2 | Unknown pattern characters match exactly | Fallback behavior if IUPAC code not recognized |
| A3 | Empty pattern returns empty collection | Implementation choice; safe default |

---

## 6. Test Consolidation Plan

### Canonical Test File

- **Location:** `SuffixTree.Genomics.Tests/IupacMotifMatchingTests.cs`
- **Contains:** All MatchesIupac and FindDegenerateMotif MUST/SHOULD tests for PAT-IUPAC-001
- **Status:** New file for comprehensive IUPAC testing

### Existing Tests to Migrate/Consolidate

From `MotifFinderTests.cs`:
- `FindDegenerateMotif_PurineR_MatchesAG` → Keep as smoke reference (duplicate logic in canonical)
- `FindDegenerateMotif_PyrimidineY_MatchesCT` → Keep as smoke reference
- `FindDegenerateMotif_AnyN_MatchesAll` → Keep as smoke reference
- `FindDegenerateMotif_WeakW_MatchesAT` → Keep as smoke reference
- `FindDegenerateMotif_ReturnsMatchedSequence` → Keep as smoke reference
- `FindDegenerateMotif_NullSequence_ThrowsException` → Move to canonical file

From `PerformanceExtensionsTests.cs`:
- `FindDegenerateMotif_WithCancellation_CompletesNormally` → Retain as smoke test

### Actions

1. Create new canonical file `IupacMotifMatchingTests.cs`
2. Implement comprehensive IUPAC code tests (all 15 codes)
3. Implement comprehensive FindDegenerateMotif tests
4. Move null check test from MotifFinderTests.cs
5. Keep minimal smoke tests in MotifFinderTests.cs (mark as delegation)
6. Keep cancellation smoke test in PerformanceExtensionsTests.cs

---

## 7. Audit Results

### Coverage Analysis of Existing Tests

| Category | Status | Notes |
|----------|--------|-------|
| Standard bases (A,C,G,T) matching | ❌ Missing | No IupacHelper tests |
| N matches all | ⚠️ Partial | FindDegenerateMotif test only |
| R (purine) | ⚠️ Partial | FindDegenerateMotif test only |
| Y (pyrimidine) | ⚠️ Partial | FindDegenerateMotif test only |
| S (strong) | ❌ Missing | No tests |
| W (weak) | ⚠️ Partial | FindDegenerateMotif test only |
| K (keto) | ❌ Missing | No tests |
| M (amino) | ❌ Missing | No tests |
| B (not A) | ❌ Missing | No tests |
| D (not C) | ❌ Missing | No tests |
| H (not G) | ❌ Missing | No tests |
| V (not T) | ❌ Missing | No tests |
| Negative tests (code does NOT match) | ❌ Missing | Critical gap |
| Empty pattern | ❌ Missing | Edge case |
| Empty sequence | ❌ Missing | Edge case |
| Pattern too long | ❌ Missing | Edge case |
| Case insensitivity | ❌ Missing | Need explicit test |
| Cancellation | ✅ Covered | PerformanceExtensionsTests.cs |
| Null input | ✅ Covered | MotifFinderTests.cs |

### Weak/Redundant Tests

- Existing tests cover only positive cases (code matches)
- No tests verify that codes do NOT match excluded bases
- No direct tests for IupacHelper.MatchesIupac

### Duplicates Found

- None across files

---

## 8. Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Test IupacHelper directly? | Yes | Canonical helper, needs comprehensive testing |
| Separate positive and negative tests? | Yes | Clear verification of both matching and non-matching |
| Include biological motif examples? | Yes (E-box) | Real-world validation |

---
