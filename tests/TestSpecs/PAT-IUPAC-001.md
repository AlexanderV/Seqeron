# Test Specification: PAT-IUPAC-001

**Test Unit ID:** PAT-IUPAC-001
**Area:** Pattern Matching
**Algorithm:** IUPAC Degenerate Motif Matching
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-02

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Nucleic acid notation | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-01-22 (verified 2026-03-02) |
| Bioinformatics.org: IUPAC codes | https://www.bioinformatics.org/sms/iupac.html | 2026-01-22 (verified 2026-03-02) |
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
| Lowercase input | Case-insensitive matching | Guaranteed: DnaSequence normalizes to uppercase; FindDegenerateMotif normalizes pattern |
| Invalid IUPAC code in pattern | ArgumentException | IUPAC-IUB 1970: only 15 standard codes defined |

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
| M53 | Mixed pattern RTG | seq="ATGCGTGC", motif="RTG" | positions [0, 4], matched ["ATG","GTG"] | Wikipedia |
| M54 | E-box CANNTG | seq="CAGCTG", motif="CANNTG" | position [0] | Biology |
| M55 | No match returns empty | seq="AAAA", motif="GGG" | [] | Standard |
| M56 | Empty pattern returns empty | seq="ACGT", motif="" | [] | Standard |
| M57 | Empty sequence returns empty | seq="", motif="ATG" | [] | Standard |
| M58 | Pattern longer than sequence | seq="AC", motif="ACGT" | [] | Standard |
| M59 | Null sequence throws | null, "ATG" | ArgumentNullException | Contract |
| M60 | Case insensitive matching | seq="ATGC", motif="atgc" | position [0] | Guaranteed |
| M62 | Result positions are valid | seq="ACGTACGTACGT", motif="NNN" | count=10, positions [0..9] | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | W (weak) pattern finds AT alternation | seq="ATATAT", motif="WWW" | positions [0, 1, 2, 3] | Wikipedia |
| S2 | S (strong) pattern finds GC regions | seq="GCGCGC", motif="SSS" | positions [0, 1, 2, 3] | Wikipedia |
| S3 | K (keto) pattern finds GT alternation | seq="GTGTGT", motif="KK" | positions [0, 1, 2, 3, 4] | Wikipedia |
| S4 | M (amino) pattern finds AC alternation | seq="ACACAC", motif="MM" | positions [0, 1, 2, 3, 4] | Wikipedia |
| S5 | Pattern at end of sequence | seq="ATGCATG", motif="ATG" | positions [0, 4] | Standard |
| S6 | Overlapping IUPAC matches | seq="AAGAG", motif="RRG" | positions [0, 2], matched ["AAG","GAG"] | Standard |
| S7 | Invalid IUPAC code in pattern throws | seq="ACGT", motif="AXG" | ArgumentException | IUPAC-IUB 1970 |
| S8 | All B positions (not A) | seq="CGT", motif="BBB" | position [0] | Wikipedia |
| S9 | All D positions (not C) | seq="AGT", motif="DDD" | position [0] | Wikipedia |
| S10 | All H positions (not G) | seq="ACT", motif="HHH" | position [0] | Wikipedia |
| S11 | All V positions (not T) | seq="ACG", motif="VVV" | position [0] | Wikipedia |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| C1 | Large sequence performance | 10000+ chars | Completes in time | Performance |
| C2 | Restriction site pattern | GAATTC, degenerate | Correct positions | Bioinformatics |

---

## 5. Coverage Classification

**Total: 74 tests in canonical file (was 85 before classification — 13 duplicates removed, 2 missing tests added)**

### 5.1 Summary

| Category | Count |
|----------|-------|
| ❌ Missing → Implemented | 2 (M57, S6) |
| ⚠ Weak → Strengthened | 6 (S1, S2, S3, S4, S5, M62) |
| 🔁 Duplicate → Removed | 17 (11 INV invariants, M61, S12 in canonical; PurineR, WeakW, MatchedSequence, NullSequence in MotifFinderTests.cs) |
| ✅ Covered | All remaining |

### 5.2 Classification Detail

#### ❌ Missing (Implemented)

| ID | Test | Action |
|----|------|--------|
| M57 | Empty sequence returns empty | Implemented `FindDegenerateMotif_EmptySequence_ReturnsEmpty` |
| S6 | Overlapping IUPAC matches | Implemented `FindDegenerateMotif_OverlappingMatches_AllReported` (fixed spec: "ARAR" invalid DNA → "AAGAG" with motif "RRG") |

#### ⚠ Weak (Strengthened)

| ID | Test | Before | After |
|----|------|--------|-------|
| S1 | W finds AT alternation | `count == 4` | `count == 4` + `positions == [0,1,2,3]` |
| S2 | S finds GC regions | `count == 4` | `count == 4` + `positions == [0,1,2,3]` |
| S3 | K finds GT alternation | `count == 5` | `count == 5` + `positions == [0,1,2,3,4]` |
| S4 | M finds AC alternation | `count == 5` | `count == 5` + `positions == [0,1,2,3,4]` |
| S5 | Pattern at end of sequence | `Does.Contain(3)` on seq="GCGATG" | `positions == [0,4]` on seq="ATGCATG" (aligned with spec) |
| M62 | Result positions valid | `matches.All(...) Is.True` | `count == 10` + `positions == [0..9]` |

#### ⚠ Weak (Strengthened — external files)

| File | Test | Before | After |
|------|------|--------|-------|
| PatternMatchingProperties.cs | IupacMatch_N_MatchesAnyBase | `Is.GreaterThan(0)` | `Is.EqualTo(5)` |

#### 🔁 Duplicate (Removed from IupacMotifMatchingTests.cs)

| Test | Reason |
|------|--------|
| `MatchesIupac_N_MatchesAllBases_Invariant` (INV-8) | Exact duplicate of M6+M7+M8+M9 |
| `MatchesIupac_R_MatchesExactlyAG_Invariant` (INV-3) | Exact duplicate of M10+M11+M12+M13 |
| `MatchesIupac_Y_MatchesExactlyCT_Invariant` | Exact duplicate of M14+M15+M16+M17 |
| `MatchesIupac_S_MatchesExactlyGC_Invariant` | Exact duplicate of M18+M19+M20+M21 |
| `MatchesIupac_W_MatchesExactlyAT_Invariant` | Exact duplicate of M22+M23+M24+M25 |
| `MatchesIupac_K_MatchesExactlyGT_Invariant` | Exact duplicate of M26+M27+M28+M29 |
| `MatchesIupac_M_MatchesExactlyAC_Invariant` | Exact duplicate of M30+M31+M32+M33 |
| `MatchesIupac_B_MatchesExactlyCGT_Invariant` | Exact duplicate of M34+M35+M36+M37 |
| `MatchesIupac_D_MatchesExactlyAGT_Invariant` | Exact duplicate of M38+M39+M40+M41 |
| `MatchesIupac_H_MatchesExactlyACT_Invariant` | Exact duplicate of M42+M43+M44+M45 |
| `MatchesIupac_V_MatchesExactlyACG_Invariant` | Exact duplicate of M46+M47+M48+M49 |
| `FindDegenerateMotif_ResultContainsMatchedSequence` (M61) | Exact duplicate of M54 (same input, same assertion) |
| `FindDegenerateMotif_WithCancellationToken_CompletesNormally` (S12) | Duplicate of `PerformanceExtensionsTests.FindDegenerateMotif_WithCancellation_CompletesNormally` |

#### 🔁 Duplicate (Removed from MotifFinderTests.cs)

| Test | Reason |
|------|--------|
| `FindDegenerateMotif_PurineR_MatchesAG` | Weaker duplicate of M53 (same pattern, weaker assertions) |
| `FindDegenerateMotif_WeakW_MatchesAT` | Exact duplicate of S1 (same input, same assertion) |
| `FindDegenerateMotif_ReturnsMatchedSequence` | Exact duplicate of M54 (same input, same assertion) |
| `FindDegenerateMotif_NullSequence_ThrowsException` | Exact duplicate of M59 |

### 5.3 Test File Map

| File | Tests | Role |
|------|-------|------|
| `IupacMotifMatchingTests.cs` | 74 | Canonical — all MUST/SHOULD tests |
| `MotifFinderTests.cs` | 2 smoke | YAT, NNG — unique scenarios not in canonical |
| `PatternMatchingProperties.cs` | 3 property | N/R/bounds property-based tests |
| `PatternMatchingSnapshotTests.cs` | 1 snapshot | RCGT snapshot |
| `PerformanceExtensionsTests.cs` | 1 | Cancellation token test |

---

## 6. Source Verification Audit (2026-03-02)

### 6.1 Sources Fetched

Both online sources were fetched and their content cross-referenced:

1. **Wikipedia: Nucleic acid notation** — Table "Single nucleobase and nucleoside" lists all 15 codes with Description, Symbol, Bases represented (columns A, C, G, T), and Complementary bases. Also states: "eleven 'ambiguity' or 'degenerate' characters associated with every possible combination of the four DNA bases."
2. **Bioinformatics.org: IUPAC codes** — Clean lookup table: 4 standard + 6 two-base + 4 three-base + 1 four-base + gap symbols.
3. **IUPAC-IUB Commission (1970)** and **NC-IUB (1984)** — Referenced via Wikipedia footnotes [1] and [2]. Wikipedia explicitly cites Biochemistry 9(20):4022–4027 (1970) and NAR 13(9):3021–3030 (1984).

### 6.2 Cross-Reference Results

All 15 IUPAC codes verified across 5 artifacts:

| Code | Wikipedia | Bioinformatics.org | IupacHelper.cs | MotifFinder IupacCodes | MotifFinder Core Switch |
|------|-----------|-------------------|----------------|----------------------|------------------------|
| A=A | ✅ | ✅ | ✅ | ✅ | ✅ |
| C=C | ✅ | ✅ | ✅ | ✅ | ✅ |
| G=G | ✅ | ✅ | ✅ | ✅ | ✅ |
| T=T | ✅ | ✅ | ✅ | ✅ | ✅ |
| R=A,G | ✅ purine | ✅ A or G | ✅ | ✅ "AG" | ✅ |
| Y=C,T | ✅ pyrimidine | ✅ C or T | ✅ | ✅ "CT" | ✅ |
| S=G,C | ✅ strong | ✅ G or C | ✅ | ✅ "GC" | ✅ |
| W=A,T | ✅ weak | ✅ A or T | ✅ | ✅ "AT" | ✅ |
| K=G,T | ✅ keto | ✅ G or T | ✅ | ✅ "GT" | ✅ |
| M=A,C | ✅ amino | ✅ A or C | ✅ | ✅ "AC" | ✅ |
| B=C,G,T | ✅ not A | ✅ C or G or T | ✅ | ✅ "CGT" | ✅ |
| D=A,G,T | ✅ not C | ✅ A or G or T | ✅ | ✅ "AGT" | ✅ |
| H=A,C,T | ✅ not G | ✅ A or C or T | ✅ | ✅ "ACT" | ✅ |
| V=A,C,G | ✅ not T | ✅ A or C or G | ✅ | ✅ "ACG" | ✅ |
| N=A,C,G,T | ✅ any | ✅ any base | ✅ | ✅ "ACGT" | ✅ |

Wikipedia full match matrix verified (each nucleotide matches exactly 8 IUPAC codes):
- A → {A, R, W, M, D, H, V, N} — tested by M1, M10, M22, M30, M38, M42, M46, M6
- C → {C, Y, S, M, B, H, V, N} — tested by M2, M14, M19, M31, M34, M43, M47, M7
- G → {G, R, S, K, B, D, V, N} — tested by M3, M11, M18, M26, M35, M39, M48, M8
- T → {T, Y, W, K, B, D, H, N} — tested by M4, M15, M23, M27, M36, M40, M44, M9

### 6.3 Source-Only Items Not Implemented (Justified)

| Source Item | Decision | Reason |
|-------------|----------|--------|
| U (Uracil) | Not supported | DNA-only system; DnaSequence validates A,C,G,T |
| Gap symbols (. or -) | Not supported | Not applicable for DNA sequence matching |
| Complement pairs (R↔Y, S↔S, etc.) | Not tested | Not part of matching algorithm; belongs to complement function |

### 6.4 Discrepancies Found

**None.** Zero discrepancies between authoritative sources and implementation/tests/spec.

---
