# Test Specification: ANNOT-REPEAT-001

**Test Unit ID:** ANNOT-REPEAT-001
**Area:** Annotation
**Algorithm:** Repetitive Element Detection and Classification
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wikipedia — Tandem repeat (cites Duitama 2014; MeSH) | 4 | https://en.wikipedia.org/wiki/Tandem_repeat | 2026-06-13 |
| 2 | Wikipedia — Inverted repeat (cites Ussery 2008; Ye 2014) | 4 | https://en.wikipedia.org/wiki/Inverted_repeat | 2026-06-13 |
| 3 | Hampson et al. (2021) IUPACpal, BMC Bioinformatics | 1 | https://doi.org/10.1186/s12859-021-03983-2 | 2026-06-13 |
| 4 | RepeatMasker documentation (Smit/Hubley/Green; Repbase) | 3 | https://www.repeatmasker.org/webrepeatmaskerhelp.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. A tandem repeat is "a pattern of one or more nucleotides … repeated and the repetitions are directly adjacent to each other" (head-to-tail), minimum two copies — Source 1.
2. STR/microsatellite motif size is 1–6 bp; minisatellite 10–60 bp — Source 1.
3. An inverted repeat is the form **WGW̄ᴿ**: left arm W, gap G (|G|≥0), right arm = reverse complement of W; gap 0 ⇒ palindrome — Sources 2, 3.
4. RepeatMasker classes: SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity, RNA, Unknown; class assigned by best homology match to a library entry — Source 4.

### 1.3 Documented Corner Cases

- Single motif occurrence (copies = 1) is not a tandem repeat (Source 1).
- Non-primitive units (e.g. "AA") should collapse to the primitive period "A" (Source 1).
- Inverted repeat with gap 0 is a reverse-complement palindrome (Sources 2, 3).
- No library match above threshold ⇒ query is Unclassified/Unknown (Source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Double-counting tandem repeats by reporting non-primitive units — Source 1.
2. Reporting an inverted repeat whose right arm is not the exact reverse complement of the left arm — Sources 2, 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindRepetitiveElements(string, int minRepeatLength, int minCopies)` | GenomeAnnotator | Canonical | Detects tandem + inverted repeats |
| `ClassifyRepeat(string sequence, IReadOnlyDictionary<string,string> repeatDb)` | GenomeAnnotator | Canonical | Assigns RepeatMasker-style class by library match, motif-size fallback |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported tandem repeat's `sequence` equals `dnaSequence[start..end]` (0-based, end-exclusive) and is exactly `copies × unitLen` bp | Yes | Source 1 (head-to-tail) |
| INV-2 | Every reported tandem repeat has ≥ `minCopies` (≥2) adjacent identical primitive units | Yes | Source 1 |
| INV-3 | For every reported inverted repeat, right arm = reverse complement of left arm | Yes | Sources 2, 3 (WGW̄ᴿ) |
| INV-4 | `ClassifyRepeat` returns a value from the source class vocabulary (library class or Simple_repeat/Unknown) | Yes | Source 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Tandem repeat detected | `ATTCGATTCGATTCG`, minLen 5, minCopies 2 | one tandem_repeat: unit `ATTCG`, start 0, end 15, sequence = full string | Source 1 worked example |
| M2 | Single copy not a repeat | `ATTCGAAAAA` motif `ATTCG` once | no tandem_repeat for `ATTCG` | Source 1 ("two or more") |
| M3 | Inverted repeat gap 0 | `GAATTC`, minArm 3 | inverted_repeat left `GAA`, right `TTC` = revcomp(`GAA`), start 0 end 6 | Sources 2,3 |
| M4 | Gapped inverted repeat | `TTACGAAAAAACGTAA`, minArm 5 | inverted_repeat left `TTACG`, right `CGTAA` = revcomp, gap 6, span [0,16) | Sources 2,3 |
| M5 | Classify by library match | seq contains library entry `GGCCGGGCGCGGTGGCTCAC`→`SINE/Alu` | returns `SINE/Alu` | Source 4 best-match |
| M6 | Classify fallback | seq `CACACACA` not in library | returns `Simple_repeat` (dinucleotide) | Source 4 + Source 1 STR 2bp |
| M7 | Classify no match non-simple | random 30bp not in library, not simple | returns `Unknown` | Source 4 Unclassified |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Primitive unit preference | `AAAAAA`, minLen 4, minCopies 2 | mononucleotide unit `A`, not `AA`/`AAA` | Source 1 non-primitive |
| S2 | Null sequence | `FindRepetitiveElements(null!)` | throws `ArgumentNullException` | Sibling analyzer contract |
| S3 | Empty sequence | `FindRepetitiveElements("")` | empty result (no throw) | Sibling contract |
| S4 | ClassifyRepeat null args | null sequence / null db | throws `ArgumentNullException` | Sibling contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | INV-1/INV-2 property | random + tandem inputs | every tandem repeat span = copies×unit, copies≥minCopies, sequence==slice | property-based, O(n²) invariant |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file targets `GenomeAnnotator.FindRepetitiveElements` or `ClassifyRepeat`. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` (GenomeAnnotator_*.cs, GenomeAnnotatorTests.cs). Related repeat tests exist for `RepeatFinder`/`GenomicAnalyzer` but cover different methods/classes.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_FindRepetitiveElements_Tests.cs` — all M/S/C cases for both methods, `#region` per method.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| GenomeAnnotator_FindRepetitiveElements_Tests.cs | canonical, this unit | 12 |

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
| 8 | S1 | ❌ Missing | implemented | ✅ Done |
| 9 | S2 | ❌ Missing | implemented | ✅ Done |
| 10 | S3 | ❌ Missing | implemented | ✅ Done |
| 11 | S4 | ❌ Missing | implemented | ✅ Done |
| 12 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | tandem repeat worked example |
| M2 | ✅ Covered | single copy excluded |
| M3 | ✅ Covered | gap-0 inverted repeat |
| M4 | ✅ Covered | gapped inverted repeat |
| M5 | ✅ Covered | library best-match class |
| M6 | ✅ Covered | simple-repeat fallback |
| M7 | ✅ Covered | unknown fallback |
| S1 | ✅ Covered | primitive unit |
| S2 | ✅ Covered | null throws |
| S3 | ✅ Covered | empty empty |
| S4 | ✅ Covered | null args throw |
| C1 | ✅ Covered | INV-1/INV-2 property |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `ClassifyRepeat` matches the library by exact substring containment of a library element within the query (one-directional, element ⊆ query) rather than Smith-Waterman homology; the class vocabulary is source-backed | M5, M6, M7 |

---

## 7. Open Questions / Decisions

1. **Decided:** suffix tree evaluated, not used for repeat *discovery* (structural unit discovery / reverse-complement arm pairing is not a fixed-pattern occurrence query). Recorded in algorithm doc §5.2.
2. **Decided:** `ClassifyRepeat` uses exact-substring library matching (Assumption 1) — documented as Simplified/Framework limitation, not an invented constant.
