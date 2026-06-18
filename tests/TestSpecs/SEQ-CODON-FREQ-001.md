# Test Specification: SEQ-CODON-FREQ-001

**Test Unit ID:** SEQ-CODON-FREQ-001
**Area:** Statistics
**Algorithm:** Codon Frequencies (non-overlapping in-frame triplet usage)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Nakamura, Gojobori, Ikemura (2000), Nucleic Acids Res 28(1):292 | 1 | https://doi.org/10.1093/nar/28.1.292 | 2026-06-14 |
| 2 | Kazusa Codon Usage Database (CUTG) README | 5 | https://www.kazusa.or.jp/codon/readme_codon.html | 2026-06-14 |
| 3 | EMBOSS `cusp` documentation | 3 | https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html | 2026-06-14 |
| 4 | Wikipedia, Codon usage bias (cited primaries) | 4 | https://en.wikipedia.org/wiki/Codon_usage_bias | 2026-06-14 |

### 1.2 Key Evidence Points

1. Codon frequency = (count of a codon / total counted codons); Kazusa reports this scaled per thousand, "calculated by summing up the numbers of codons used." — Source 2.
2. Codons are non-overlapping in-frame triplets read from CDS. — Sources 2, 3.
3. Codons containing ambiguous (non-ACGT) bases are excluded from the count. — Source 2.
4. EMBOSS cusp sample output: Σ Number = 386 codons; CGC=22 → 56.995 per thousand; 22/386×1000 = 56.995, confirming count/total = per-thousand ÷ 1000. — Source 3.
5. cusp "Fraction" (per-amino-acid proportion) is a distinct metric and is NOT the value computed here. — Source 3.

### 1.3 Documented Corner Cases

- Ambiguous / non-ACGT codon excluded from count (Source 2).
- Trailing 1–2 leftover bases form no codon and are ignored (count/total definition, Source 2).
- Reading frame offset changes the codon multiset (non-overlapping triplets from frame start, Source 3).
- All-ambiguous or shorter-than-3 input → total = 0 → empty table (only well-defined result; no authoritative source defines a non-empty result here — see §6).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing per-amino-acid Fraction with count/total frequency — EMBOSS cusp documents both; only count/total is in scope (Source 3).
2. Counting overlapping triplets instead of non-overlapping in-frame triplets (Sources 2, 3).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateCodonFrequencies(string dnaSequence, int readingFrame = 0)` | SequenceStatistics | Canonical | count/total over non-overlapping in-frame triplets |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Each frequency is in (0, 1]; only observed codons appear as keys | Yes | count/total definition, Source 2 |
| INV-02 | Frequencies over all counted codons sum to 1.0 (when ≥1 valid codon) | Yes | count/total normalization, Source 2 |
| INV-03 | Codons with any non-ACGT base never appear and never affect the total | Yes | Source 2 ("ambiguous excluded") |
| INV-04 | Result is independent of input letter case | Yes | codons are case-independent; impl upper-cases |
| INV-05 | count/total fraction = Kazusa per-thousand frequency ÷ 1000 | Yes | Source 3 (cusp 22/386×1000 = 56.995) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Frame-0 exact frequencies | `ATGATGAAA` frame 0 → ATG, ATG, AAA | ATG=2/3, AAA=1/3, keys={ATG,AAA} | Source 2 count/total |
| M2 | Reading-frame offset | `ATGATGAAA` frame 1 → TGA, TGA | TGA=1.0, keys={TGA} | Source 3 non-overlap from frame |
| M3 | Sum to one (INV-02) | `ATGATGAAA` frame 0 | Σ freq = 1.0 | Source 2 normalization |
| M4 | Non-ACGT excluded (INV-03) | `ATGNNNAAA` frame 0 → ATG, AAA | ATG=1/2, AAA=1/2, no NNN key | Source 2 |
| M5 | cusp cross-check (INV-05) | Multiset of 386 codons matching cusp: build CGC×22 etc.; minimal proxy `CGCCGCCGC`+filler reproducing fraction | a codon with count k over total n gives k/n = per-thousand/1000 | Source 3 (22/386=0.056995) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Trailing bases ignored | `ATGAA` frame 0 | ATG=1.0, keys={ATG} | remainder rule |
| S2 | Case-insensitive (INV-04) | `atgaaa` equals `ATGAAA` | ATG=1/2, AAA=1/2 | upper-cased |
| S3 | All-ambiguous → empty | `NNNNNN` frame 0 | empty dictionary (total=0) | zero-codon corner |
| S4 | Values in (0,1] (INV-01) | `ATGATGAAA` frame 0 | all values > 0 and ≤ 1 | bound |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | null/empty/short guard | null, "", "AT" | empty dictionary | shared guard contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatisticsTests.cs` — legacy `#region Codon Frequency Tests` with 3 weak tests (`Is.Not.Empty`, `ContainsKey`, `Is.Not.EquivalentTo`) — no exact values.
- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateDinucleotide_Tests.cs` (SEQ-DINUC-001) — `#region CalculateCodonFrequencies` with 6 evidence-based tests for this method, bundled into the dinucleotide unit.
- Implementation already present: `SequenceStatistics.CalculateCodonFrequencies` (SequenceStatistics.cs).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 frame-0 exact | 🔁 Duplicate | exists in SEQ-DINUC-001 file (M4); belongs to this unit — relocate to canonical file |
| M2 frame offset | 🔁 Duplicate | exists in SEQ-DINUC-001 file (M5); relocate |
| M3 sum to one | 🔁 Duplicate | exists in SEQ-DINUC-001 file (M6); relocate |
| M4 non-ACGT excluded | 🔁 Duplicate | exists in SEQ-DINUC-001 file (S5); relocate |
| M5 cusp cross-check | ❌ Missing | not present anywhere |
| S1 trailing ignored | 🔁 Duplicate | exists in SEQ-DINUC-001 file (S6); relocate |
| S2 case-insensitive | ❌ Missing | not present |
| S3 all-ambiguous empty | ❌ Missing | not present |
| S4 values bound | ❌ Missing | not present |
| C1 guards | 🔁 Duplicate | exists in SEQ-DINUC-001 file (S3) + weak legacy; relocate, drop legacy |
| legacy `_ReturnsFrequencies` | ⚠ Weak | `Is.Not.Empty`/`ContainsKey` — remove |
| legacy `_DifferentReadingFrame` | ⚠ Weak | `Is.Not.EquivalentTo` permissive — remove |
| legacy `_ShortSequence_ReturnsEmpty` | 🔁 Duplicate | covered by C1 — remove |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateCodonFrequencies_Tests.cs` — owns all `CalculateCodonFrequencies` tests for SEQ-CODON-FREQ-001 (M1–M5, S1–S4, C1).
- **Remove:** the entire `#region CalculateCodonFrequencies` from `SequenceStatistics_CalculateDinucleotide_Tests.cs` (those tests belong to this unit, not SEQ-DINUC-001); update that TestSpec note.
- **Remove:** the legacy `#region Codon Frequency Tests` (3 weak tests) from `SequenceStatisticsTests.cs`.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateCodonFrequencies_Tests.cs` | Canonical (this unit) | 10 |
| `SequenceStatistics_CalculateDinucleotide_Tests.cs` | SEQ-DINUC-001 (codon region removed) | unchanged minus 6 |
| `SequenceStatisticsTests.cs` | legacy smoke (codon region removed) | unchanged minus 3 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | 🔁 Duplicate | moved into canonical file | ✅ Done |
| 2 | M2 | 🔁 Duplicate | moved into canonical file | ✅ Done |
| 3 | M3 | 🔁 Duplicate | moved into canonical file | ✅ Done |
| 4 | M4 | 🔁 Duplicate | moved into canonical file | ✅ Done |
| 5 | M5 | ❌ Missing | implemented (cusp cross-check) | ✅ Done |
| 6 | S1 | 🔁 Duplicate | moved into canonical file | ✅ Done |
| 7 | S2 | ❌ Missing | implemented (case-insensitive) | ✅ Done |
| 8 | S3 | ❌ Missing | implemented (all-ambiguous empty) | ✅ Done |
| 9 | S4 | ❌ Missing | implemented (values bound) | ✅ Done |
| 10 | C1 | 🔁 Duplicate | moved into canonical file; legacy removed | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | canonical file |
| M2 | ✅ Covered | canonical file |
| M3 | ✅ Covered | canonical file |
| M4 | ✅ Covered | canonical file |
| M5 | ✅ Covered | canonical file |
| S1 | ✅ Covered | canonical file |
| S2 | ✅ Covered | canonical file |
| S3 | ✅ Covered | canonical file |
| S4 | ✅ Covered | canonical file |
| C1 | ✅ Covered | canonical file (legacy + dinuc duplicates removed) |

In-scope cases: 10. ✅ count: 10. No ❌ / ⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Zero valid codons (total=0) → empty table. Non-correctness-affecting for any input with ≥1 valid codon; only well-defined choice consistent with count/total. | S3 |

---

## 7. Open Questions / Decisions

1. Decision: `CalculateCodonFrequencies` was previously implemented and tested inside SEQ-DINUC-001's file. It is the dedicated method of THIS unit, so its tests are relocated to the canonical file for SEQ-CODON-FREQ-001 and weak legacy duplicates are removed (one canonical file per unit). No behavioral change to the implementation.
2. The method is not a search/matching operation; the repository suffix tree is N/A (linear single pass).
