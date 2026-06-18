# Test Specification: MOTIF-CONS-001

**Test Unit ID:** MOTIF-CONS-001
**Area:** Matching
**Algorithm:** Consensus Sequence from a Multiple Alignment (most-frequent residue)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — Consensus sequence (citing Schneider & Stephens 1990) | 4 | https://en.wikipedia.org/wiki/Consensus_sequence | 2026-06-13 |
| 2 | Rosalind — Consensus and Profile (CONS) | 5 | https://rosalind.info/problems/cons/ | 2026-06-13 |
| 3 | EMBOSS `cons` documentation (Rice et al. 2000) | 3 | https://www.bioinformatics.nl/cgi-bin/emboss/help/cons | 2026-06-13 |
| 4 | Los Alamos HIV DB — Advanced Consensus explanation | 5 | https://hfv.lanl.gov/content/sequence/CONSENSUS/AdvConExplain.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Consensus = "the calculated sequence of most frequent residues … found at each position in a sequence alignment." — Source 1.
2. The jth consensus symbol "corresponds to the symbol having the maximum value in the j-th column of the profile matrix." — Source 2.
3. Worked example (7 strings of length 8) yields consensus `ATGCAACT` with profile A=`5 1 0 0 5 5 0 0`, C=`0 0 1 4 2 0 6 1`, G=`1 1 6 3 0 1 0 0`, T=`1 5 0 0 0 1 1 6`. — Source 2.
4. On a tie, "the residue letter occurring earlier in the alphabet was chosen" (alphabetical tie-break). — Source 4 (and Geneious manual cited therein).
5. The consensus is defined over equal-length strings. — Source 2.

### 1.3 Documented Corner Cases

- Equal-length precondition (Source 2): consensus defined only for equal-length aligned strings.
- Tie columns (Source 2): more than one valid consensus; a deterministic implementation must fix a rule (alphabetical, per Source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Non-deterministic output if tie-breaking is unspecified — resolved by alphabetical tie-break. — Source 4.
2. Misapplying a plurality threshold (EMBOSS) to the plain most-frequent definition — out of scope; this method takes no threshold (see §6). — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CreateConsensusFromAlignment(IEnumerable<string>)` | MotifFinder | Canonical | Most-frequent residue per column; alphabetical tie-break |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Output length equals the common input string length | Yes | Source 2 (one symbol per column) |
| INV-2 | Each output symbol is a symbol present in its column with maximum count | Yes | Source 2 |
| INV-3 | On a tie, the alphabetically-earliest tied symbol is chosen (A<C<G<T) | Yes | Source 4 |
| INV-4 | Identical input sequences ⇒ output equals that sequence | Yes | Source 1, 2 |
| INV-5 | Deterministic: same input ⇒ same output | Yes | INV-3 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Rosalind CONS sample | 7 strings of length 8 | `ATGCAACT` | Source 2 sample |
| M2 | Identical sequences | 3× `ACGT` | `ACGT` | Source 1, 2 (INV-4) |
| M3 | Alphabetical tie-break | `AT`,`GT` (col1 A vs G tie) | `AT` | Source 4 (INV-3) |
| M4 | Single sequence | `["GATTACA"]` | `GATTACA` | Source 2 (most-frequent of one) |
| M5 | Case-insensitive | `acgt`,`ACGT` | `ACGT` | Sibling convention (`ToUpperInvariant`) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null input | `null` collection | `ArgumentNullException` | Library convention |
| S2 | Empty collection | `[]` | `""` | Mirrors `GenerateConsensus` |
| S3 | Unequal lengths | `AC`,`ACG` | `ArgumentException` | Source 2 equal-length precondition |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Invalid character | `AX` | `ArgumentException` | Alphabet validation as in `CreatePwm` |
| C2 | Three-way majority over tie | `A,A,C` column | `A` (count 2 > 1) | Pure majority, no tie |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `CreateConsensusFromAlignment`: no existing tests. `MotifFinder_*` test files cover PWM, IUPAC, regulatory motifs only.
- No existing implementation of `CreateConsensusFromAlignment` in `MotifFinder.cs` (a separate IUPAC-degenerate `GenerateConsensus` exists; it is not this unit).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_CreateConsensusFromAlignment_Tests.cs` — all cases for this unit.
- **Remove:** none (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_CreateConsensusFromAlignment_Tests.cs` | Canonical | 10 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented Rosalind sample test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented identical-sequences test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented alphabetical tie-break test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented single-sequence test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented case-insensitive test | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented null test | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented empty-collection test | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented unequal-length test | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented invalid-character test | ✅ Done |
| 10 | C2 | ❌ Missing | Implemented majority-over-tie test | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Rosalind sample → `ATGCAACT` |
| M2 | ✅ Covered | identical → unchanged |
| M3 | ✅ Covered | tie → A |
| M4 | ✅ Covered | single → unchanged |
| M5 | ✅ Covered | lowercase normalised |
| S1 | ✅ Covered | ArgumentNullException |
| S2 | ✅ Covered | "" |
| S3 | ✅ Covered | ArgumentException |
| C1 | ✅ Covered | ArgumentException |
| C2 | ✅ Covered | majority A |

In-scope cases: 10. ✅ Covered: 10.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Alphabetical tie-break (A<C<G<T) for determinism (Geneious/LANL rule) | INV-3, M3 |
| 2 | Pure most-frequent consensus, no plurality threshold (matches Registry signature) | §6 scope |

---

## 7. Open Questions / Decisions

1. Decision: tie-breaking fixed to alphabetical order to guarantee determinism (Source 4). Rosalind permits any tied symbol, so the rank-5 worked example (no decisive ties) remains conformant.
