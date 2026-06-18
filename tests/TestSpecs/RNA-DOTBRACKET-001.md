# Test Specification: RNA-DOTBRACKET-001

**Test Unit ID:** RNA-DOTBRACKET-001
**Area:** RnaStructure
**Algorithm:** Dot-Bracket (extended WUSS) Notation — parsing and validation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | ViennaRNA — RNA Structure Notations | 3 | https://viennarna.readthedocs.io/en/latest/io/rna_structures.html | 2026-06-14 |
| 2 | ViennaRNA — Dot-Bracket Notation | 3 | https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/dotbracket.html | 2026-06-14 |
| 3 | ViennaRNA — WUSS notation | 3 | https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/wuss.html | 2026-06-14 |
| 4 | Nawrocki & Eddy (2013), Infernal 1.1 (WUSS) | 1 | https://doi.org/10.1093/bioinformatics/btt509 | 2026-06-14 |
| 5 | Rfam Documentation — Glossary (WUSS) | 5 | https://docs.rfam.org/en/latest/glossary.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. `()` = base pair, `.` = unpaired; structure must be balanced/nested — Source 1, 2.
2. Extended families `<>`, `{}`, `[]` and uppercase/lowercase letter pairs are independent pairing systems (need not nest with each other) → enable pseudoknots — Source 1, 2, 3.
3. Equivalent crossing examples `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`, `AAAA{{{{....aaaa}}}}`; uppercase letter = 5' opener, lowercase = 3' closer — Source 1, 2.
4. A closing symbol must match an opening of the SAME family ("partners must match up") — Source 3, 4.
5. Non-bracket WUSS symbols `-`, `,`, `:`, `.` are all single-stranded (unpaired) — Source 5.

### 1.3 Documented Corner Cases

- Crossing families (`([)]`) must be matched on separate stacks (Source 1/3).
- Mismatched families (`(]`) are not well-formed even when one open + one close are present (Source 3/4).
- Non-bracket WUSS symbols and dots are unpaired, not errors (Source 5).

### 1.4 Known Failure Modes / Pitfalls

1. Single shared stack mis-pairs crossing families (e.g. `([)]` → wrong pairs) — Source 1/3.
2. Counting all bracket types in one balance counter accepts mismatched families like `(]` — Source 3/4.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ParseDotBracket(string)` | RnaSecondaryStructure | **Canonical** | returns (Position1, Position2) 0-based pairs |
| `ValidateDotBracket(string)` | RnaSecondaryStructure | **Canonical** | well-formedness per-family |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Each returned pair has Position1 < Position2 (opener precedes closer) | Yes | Source 1/2 (opening before closing) |
| INV-2 | A returned pair's two indices carry the same bracket family / letter case-pair | Yes | Source 3/4 ("partners must match up") |
| INV-3 | For a string `ValidateDotBracket` accepts, `ParseDotBracket` returns exactly (#opening symbols) pairs | Yes | Source 1 (balanced) |
| INV-4 | `ValidateDotBracket` is true iff every closer matches an earlier opener of the same family and no opener is left unclosed | Yes | Source 3/4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Parse simple hairpin | `((((....))))` | pairs {(0,11),(1,10),(2,9),(3,8)} | Source 1/2 example |
| M2 | Parse crossing families | `([)]` | pairs {(0,2),(1,3)} | Source 1/3 independent stacks |
| M3 | Parse bracket-vs-letter equivalence | `<<<<[[[[....>>>>]]]]` and `((((AAAA....))))aaaa` yield identical pair sets (families nest LIFO) | both → {(0,15),(1,14),(2,13),(3,12),(4,19),(5,18),(6,17),(7,16)} | Source 1/2 equivalent examples |
| M4 | Validate balanced/nested | `(((...)))`, `(([[]]))`, `([)]`, `....`, `""` | all true | Source 1 (balanced) + Source 3 (crossing families) |
| M5 | Validate malformed | `(((...)`, `...)`, `)(`, `(]` | all false | Source 1 (balanced) + Source 3/4 (matching partners) |
| M6 | Parse count == openers (INV-3) | `(([[]]))` → 4 pairs | exactly 4 pairs | Source 1 balanced |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Letter pair direction | `AAAA....aaaa` | pairs {(0,11),(1,10),(2,9),(3,8)} (uppercase opens) | Source 1/2 `AAAA...aaaa` |
| S2 | Non-bracket WUSS symbols unpaired | `<<<-->>>` and `((,,))` | `<<<-->>>`→3 pairs (0,7),(1,6),(2,5); `((,,))`→2 pairs | Source 5 single-stranded symbols |
| S3 | Best-effort parse of stray closer | `())` | yields only (0,1) | Evidence Assumption (contract) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty / all-dots / null | `""`, `.....`, null | parse empty; validate("")/validate(null) true | Evidence Assumption |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructureTests.cs` contained 5 pre-existing dot-bracket tests (`ParseDotBracket_SimpleStructure`, `_EmptyStructure`, `_MultipleBrackets`, `ValidateDotBracket_Balanced`, `_Unbalanced`).
- No `RnaSecondaryStructure_ParseDotBracket_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (simple hairpin positions) | ⚠ Weak | old test asserts 3 pairs (no exact full set / nesting), pre-template |
| M2 (crossing families) | ❌ Missing | old code/tests would mis-pair `([)]` |
| M3 (bracket/letter equivalence) | ❌ Missing | letters not handled by old code |
| M4 (validate balanced incl. crossing) | ⚠ Weak | old test lacks crossing-family case |
| M5 (validate malformed incl. `(]`) | ⚠ Weak | old test never checks mismatched families |
| M6 (count == openers) | ⚠ Weak | old `_MultipleBrackets` checks count only, no positions |
| S1 (letter direction) | ❌ Missing | — |
| S2 (WUSS unpaired symbols) | ❌ Missing | — |
| S3 (best-effort stray closer) | ❌ Missing | — |
| C1 (empty/all-dots/null) | ⚠ Weak | old empty test exists but no null |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_ParseDotBracket_Tests.cs` — all M/S/C cases with exact evidence values.
- **Remove:** the 5 pre-existing dot-bracket tests in `RnaSecondaryStructureTests.cs` (superseded; avoid duplicate coverage).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| RnaSecondaryStructure_ParseDotBracket_Tests.cs | canonical | 11 |
| RnaSecondaryStructureTests.cs | other RNA units; dot-bracket tests removed | n/a |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with exact full pair set | ✅ Done |
| 2 | M2 | ❌ Missing | implemented crossing-family parse | ✅ Done |
| 3 | M3 | ❌ Missing | implemented bracket/letter equivalence | ✅ Done |
| 4 | M4 | ⚠ Weak | rewrote incl. `([)]` crossing | ✅ Done |
| 5 | M5 | ⚠ Weak | rewrote incl. `(]` mismatch | ✅ Done |
| 6 | M6 | ⚠ Weak | rewrote count + positions | ✅ Done |
| 7 | S1 | ❌ Missing | implemented letter direction | ✅ Done |
| 8 | S2 | ❌ Missing | implemented WUSS unpaired symbols | ✅ Done |
| 9 | S3 | ❌ Missing | implemented best-effort parse | ✅ Done |
| 10 | C1 | ⚠ Weak | rewrote incl. null | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact pair set asserted |
| M2 | ✅ | (0,2),(1,3) asserted |
| M3 | ✅ | bracket & letter sets equal |
| M4 | ✅ | balanced + crossing accepted |
| M5 | ✅ | unclosed/unopened/reversed/mismatched rejected |
| M6 | ✅ | 4 pairs + positions |
| S1 | ✅ | uppercase opener verified |
| S2 | ✅ | WUSS symbols unpaired |
| S3 | ✅ | stray closer dropped |
| C1 | ✅ | empty/all-dots/null |

In-scope cases: 10. ✅ = 10.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Best-effort parse of malformed input (drop unmatched closers, no throw) | S3 |
| 2 | Empty/null = valid pair-free structure | M4, C1 |

---

## 7. Open Questions / Decisions

1. Decision: each bracket family and each letter case-pair is matched on its own stack (independent pairing systems) per ViennaRNA/WUSS — resolved from evidence, no open questions.
