# Test Specification: MOTIF-REGULATORY-001

**Test Unit ID:** MOTIF-REGULATORY-001
**Area:** Matching
**Algorithm:** Regulatory Elements (scan a DNA sequence for known regulatory consensus motifs)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-29

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Bucher (1990) J Mol Biol — eukaryotic Pol II promoter weight matrices (TATA, CCAAT) | 1 | https://doi.org/10.1016/0022-2836(90)90223-9 | 2026-06-14 |
| 2 | Harley & Reynolds (1987) Nucleic Acids Res — E. coli -10/-35 hexamers | 1 | https://doi.org/10.1093/nar/15.5.2343 | 2026-06-14 |
| 3 | Dynan & Tjian (1983) Cell 35:79; Gidoni, Dynan & Tjian (1984) Nature 312:409 — Sp1 GC box GGGCGG (replaces Lundin et al. 1994, which is about yeast MIG1) | 1 | WebSearch record (2026-09-29) | 2026-09-29 |
| 4 | Kozak (1987) Nucleic Acids Res — Kozak gccRccATGG | 1 | https://doi.org/10.1093/nar/15.20.8125 | 2026-06-14 |
| 5 | Proudfoot & Brownlee (1976) Nature — poly(A) AATAAA | 1 | https://doi.org/10.1038/263211a0 | 2026-06-14 |
| 6 | Massari & Murre (2000) Mol Cell Biol — E-box CANNTG | 1 | https://doi.org/10.1128/MCB.20.2.429-440.2000 | 2026-06-14 |
| 7 | Lee, Mitchell & Tjian (1987) Cell — AP-1 TGACTCA | 1 | https://pubmed.ncbi.nlm.nih.gov/3034433/ | 2026-06-14 |
| 8 | Sen & Baltimore (1986) Cell — NF-κB κB site | 1 | https://doi.org/10.1016/0092-8674(86)90346-6 | 2026-06-14 |
| 9 | Montminy et al. (1986) PNAS — CREB CRE TGACGTCA | 1 | https://doi.org/10.1073/pnas.83.18.6682 | 2026-06-14 |
| 10 | Shine & Dalgarno (1974) PNAS 71:1342 — 16S rRNA 3' ACCUCCUUA; SD = AGGAGG | 1 | https://doi.org/10.1073/pnas.71.4.1342 (WebSearch record) | 2026-09-29 |
| 11 | Angel et al. (1987) Cell 49:729 — TRE TGAGTCA (collagenase), consensus TGA(C/G)TCA | 1 | https://doi.org/10.1016/0092-8674(87)90611-8 (WebSearch record) | 2026-09-29 |
| 12 | Gilmore (2006) Oncogene 25:6680 — κB consensus GGGRNWYYCC | 2 | https://doi.org/10.1038/sj.onc.1209954 (WebSearch record) | 2026-09-29 |
| 13 | Biopython 1.88 `Bio.SeqUtils.nt_search` — reference scan (IUPAC in pattern) | ref impl | installed package | 2026-09-29 |

### 1.2 Key Evidence Points

1. TATA box consensus = `TATAAA` — Bucher (1990) / Wikipedia TATA box.
2. -10 (Pribnow) box = `TATAAT`, -35 box = `TTGACA` — Harley & Reynolds (1987).
3. CCAAT box pentanucleotide = `CCAAT` (~30% of promoters) — Bucher (1990).
4. GC box = `GGGCGG` — Dynan & Tjian (1983); Gidoni, Dynan & Tjian (1984).
5. Kozak consensus = `GCCGCC(A/G)CCATGG` = `GCCGCCRCCATGG`; -3 purine, +4 G — Kozak (1987) abstract.
6. Shine-Dalgarno = `AGGAGG`, complementary to 3' 16S rRNA — Shine & Dalgarno (1974).
7. Poly(A) signal = `AATAAA` — Proudfoot & Brownlee (1976).
8. E-box = `CANNTG` (IUPAC N), canonical palindrome `CACGTG` — Massari & Murre (2000).
9. AP-1 (TRE) consensus = `TGA(C/G)TCA` = `TGASTCA` — Lee, Mitchell & Tjian (1987); the collagenase TRE is `TGAGTCA` (Angel et al. 1987), the reverse complement of `TGACTCA`.
10. NF-κB κB consensus = `GGGRNWYYCC` (Gilmore 2006); includes the Ig κ site `GGGACTTTCC` (Sen & Baltimore 1986).
11. CREB CRE palindrome = `TGACGTCA` — Montminy et al. (1986).

### 1.3 Documented Corner Cases

- E-box `CANNTG`, Kozak `GCCGCCRCCATGG`, AP-1 `TGASTCA` and κB `GGGRNWYYCC` are degenerate; matching uses the canonical IUPAC scan (`IupacHelper.MatchesIupac`), identical to Biopython `nt_search`.
- Only the given strand is scanned; AP-1/E-box/CREB patterns are self-reverse-complementary.
- Multiple occurrences of one element and occurrences of several distinct elements may co-exist; each is reported with a 0-based start position.
- Empty sequence → no occurrences. Null sequence → `ArgumentNullException`.

### 1.4 Known Failure Modes / Pitfalls

1. Scanning AP-1 as only one of `TGACTCA`/`TGAGTCA` misses half of the TREs: the consensus is `TGA(C/G)TCA` and the two strings are reverse complements. (An earlier revision declared `TGAGTCA` wrong; corrected 2026-09.)
2. Confusing eukaryotic TATA (`TATAAA`) with the prokaryotic -10 hexamer (`TATAAT`) — distinct elements, differ only in the last base — Harley & Reynolds (1987).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindRegulatoryElements(DnaSequence)` | MotifFinder | **Canonical** | Scans for the 12 known regulatory consensus motifs. |
| `KnownMotifs` (consensus constants) | MotifFinder.KnownMotifs | **Internal** | Consensus strings verified equal to cited sources. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported element's `Sequence` (the matched substring) has the same length as its `Pattern`, and occurs at `Position` (0-based) in the input. | Yes | Scanning-window definition |
| INV-2 | Each reported element's `Sequence` IUPAC-matches its `Pattern`. | Yes | IUPAC matching (FindDegenerateMotif) |
| INV-3 | A consensus string and its element name are exactly those in the cited primary sources (no fabricated constants). | Yes | Sources §1.1 |
| INV-4 | Scanning is exhaustive: all 0-based start positions `0 <= i <= n - m` that match are reported (no missed/extra occurrences). | Yes | Window definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | TATA box detected | `GGGTATAAAGGG` | one "TATA Box" at pos 3, Pattern `TATAAA`, Sequence `TATAAA` | Bucher (1990) |
| M2 | -10 box detected | `CCTATAATCC` | one "-10 Box" at pos 2, Pattern `TATAAT` | Harley & Reynolds (1987) |
| M3 | -35 box detected | `AATTGACAGG` | one "-35 Box" at pos 2, Pattern `TTGACA` | Harley & Reynolds (1987) |
| M4 | CAAT box detected | `GGCCAATGG` | one "CAAT Box" at pos 2, Pattern `CCAAT` | Bucher (1990) |
| M5 | GC box detected | `AAGGGCGGTT` | one "GC Box" at pos 2, Pattern `GGGCGG` | Dynan & Tjian (1983) |
| M6 | Kozak detected | `TTGCCGCCACCATGGAA` | one "Kozak" at pos 2, Pattern `GCCGCCRCCATGG`, Sequence `GCCGCCACCATGG` | Kozak (1987) |
| M7 | Shine-Dalgarno detected | `TTAGGAGGTTT` | one "Shine-Dalgarno" at pos 2, Pattern `AGGAGG` | Shine & Dalgarno (1974) |
| M8 | Poly(A) signal detected | `CCAATAAACC` | one "Poly(A) Signal" at pos 2, Pattern `AATAAA` | Proudfoot & Brownlee (1976) |
| M9 | E-box degenerate match | `GGCACGTGGG` | one "E-box" at pos 2, Pattern `CANNTG`, Sequence `CACGTG` | Massari & Murre (2000) |
| M10 | AP-1 detected | `AATGACTCAGG` | one "AP-1" at pos 2, Pattern `TGASTCA`, Sequence `TGACTCA` | Lee, Mitchell & Tjian (1987) |
| M11 | Collagenase TRE detected | `AATGAGTCAGG` | one "AP-1" at pos 2, Sequence `TGAGTCA` (nt_search → [2]); `AATGATTCAGG` → no AP-1 | Angel et al. (1987) |
| M12 | NF-κB detected | `AAGGGACTTTCCAA` | one "NF-κB" at pos 2, Pattern `GGGRNWYYCC`, Sequence `GGGACTTTCC` | Gilmore (2006); Sen & Baltimore (1986) |
| M13 | CREB detected | `CCTGACGTCAGG` | one "CREB" at pos 2, Pattern `TGACGTCA` | Montminy et al. (1986) |
| M14 | Null sequence | `null` | `ArgumentNullException` | Contract |
| M15 | Empty sequence | `""` | empty result | Window definition |
| M16 | Constants equal source consensus | `KnownMotifs.*` | each constant equals its cited consensus string | §1.1 sources / INV-3 |
| M17 | Kozak -3 purine | `TTGCCGCCGCCATGGAA` / `TTGCCGCCTCCATGGAA` | Kozak at 2 (`GCCGCCGCCATGG`) / none | Kozak (1987); nt_search |
| M18 | κB variant sites | `AAGGGAAATTCCAA` / `AAGGGACGTTCCAA` | NF-κB at 2 (`GGGAAATTCC`) / none | Gilmore (2006); nt_search |
| M19 | Full output = Biopython nt_search | CCAATAAACC, GTATAATATAAA, CACATGTG, TATAAAAGGAGG, AATAAACGAATAAA, GCCGCCACCATG | exact (Name:Pos:Seq) lists in library order | Biopython 1.88 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Multiple occurrences of one element | `AATAAACGAATAAA` | two "Poly(A) Signal" hits at pos 0 and 8 | INV-4 exhaustiveness |
| S2 | Multiple distinct elements | `TATAAA` + `AGGAGG` in one sequence | both "TATA Box" and "Shine-Dalgarno" reported | corner case |
| S3 | No regulatory element present | `GGGGCCCCGGGG` (no consensus) | empty result | corner case |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Matched substring length invariant | any detected element | `Sequence.Length == Pattern.Length` | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinderTests.cs` contained a legacy "Regulatory Element Tests" region (`FindRegulatoryElements_FindsTataBox/FindsPolyASignal/FindsEBox/ReturnsDescription`, `KnownMotifs_ContainsExpectedPatterns`) plus `FindRegulatoryElements_NullSequence_ThrowsException` in the Edge Cases region.
- No canonical `MotifFinder_FindRegulatoryElements_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 TATA box | ⚠ Weak | legacy used permissive `.Any(e => e.Name=="TATA Box")`, no position/pattern check |
| M2 -10 box | ❌ Missing | element did not exist before this unit |
| M3 -35 box | ❌ Missing | element did not exist before this unit |
| M4 CAAT box | ❌ Missing | not tested |
| M5 GC box | ❌ Missing | not tested |
| M6 Kozak | ❌ Missing | not tested |
| M7 Shine-Dalgarno | ❌ Missing | only asserted via constant |
| M8 Poly(A) signal | ⚠ Weak | legacy `.Any()` only |
| M9 E-box degenerate | ⚠ Weak | legacy `.Any()` only |
| M10 AP-1 corrected | ❌ Missing | not tested; pattern was wrong |
| M11 AP-1 wrong-pattern regression | ❌ Missing | not tested |
| M12 NF-κB | ❌ Missing | not tested |
| M13 CREB | ❌ Missing | not tested |
| M14 null | ✅ Covered | legacy null test (moved to canonical file) |
| M15 empty | ❌ Missing | not tested |
| M16 constants | ⚠ Weak | legacy checked only 3 constants |
| S1 multiple occurrences | ❌ Missing | not tested |
| S2 multiple distinct elements | ❌ Missing | not tested |
| S3 no element | ❌ Missing | not tested |
| C1 length invariant | ❌ Missing | not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_FindRegulatoryElements_Tests.cs` — all MUST/SHOULD/COULD cases.
- **Remove:** the "Regulatory Element Tests" region and the `FindRegulatoryElements_NullSequence_ThrowsException` test from `MotifFinderTests.cs` (duplicates/weak), replaced by a NOTE pointer.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MotifFinder_FindRegulatoryElements_Tests.cs` | Canonical MOTIF-REGULATORY-001 | 20 |
| `MotifFinderTests.cs` | Other MotifFinder methods (regulatory region removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with exact pos/pattern/sequence | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ⚠ Weak | rewrote with exact pos | ✅ Done |
| 9 | M9 | ⚠ Weak | rewrote with exact pos/sequence | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | M11 | ❌ Missing | implemented (regression) | ✅ Done |
| 12 | M12 | ❌ Missing | implemented | ✅ Done |
| 13 | M13 | ❌ Missing | implemented | ✅ Done |
| 14 | M14 | ✅ Covered | moved to canonical file | ✅ Done |
| 15 | M15 | ❌ Missing | implemented | ✅ Done |
| 16 | M16 | ⚠ Weak | rewrote — all 12 constants vs sources | ✅ Done |
| 17 | S1 | ❌ Missing | implemented | ✅ Done |
| 18 | S2 | ❌ Missing | implemented | ✅ Done |
| 19 | S3 | ❌ Missing | implemented | ✅ Done |
| 20 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M16 | ✅ | canonical file, exact evidence-based values |
| S1–S3 | ✅ | canonical file |
| C1 | ✅ | canonical file (invariant) |

All 20 in-scope cases ✅.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | TATA/CCAAT/GC box scanned as core consensus strings rather than Bucher (1990) weight matrices; single-strand scan (declared in the algorithm doc §5.3). | M1, M4, M5 |

(2026-09: former assumptions "NF-κB as GGGACTTTCC" and "Kozak as GCCGCCACCATGG" removed — the published IUPAC consensus is now scanned.)

## 7. Open Questions / Decisions

1. **Decision (revised 2026-09):** AP-1 = `TGASTCA` (TGA(C/G)TCA). The 2026-06 change `TGAGTCA` → `TGACTCA` only swapped one strand's representation of the same site for the other; both are now reported.
2. **Decision:** Added prokaryotic -10 (`TATAAT`) and -35 (`TTGACA`) hexamers (Harley & Reynolds 1987), which the task note lists as expected regulatory elements and which were absent before.
3. **Decision:** Used the per-position IUPAC scan (`FindDegenerateMotif`); the suffix tree was not used because the IUPAC E-box pattern is degenerate (not an exact substring), so a single linear scan over a short pattern set is the correct algorithm.

## 8. Review 2026-09 (B05 follow-up) — strands and Bucher matrices

Tests: `Unit/Analysis/MotifFinder_RegulatoryStrands_Tests.cs`, `Unit/Analysis/MotifFinder_PwmStrandsAndThresholds_Tests.cs` (Bucher part), `Metamorphic/MotifPwmMetamorphicTests.cs` (REG-BOTH).

| ID | Test | Locked values / invariant |
|----|------|---------------------------|
| R1 | `FindRegulatoryElements_BothStrands_EqualsBiopythonNtSearch` (2 cases) | Biopython `nt_search` on both strands (Evidence, 2026-09 follow-up) |
| R2 | `FindRegulatoryElements_SingleStrand_EqualsLegacyOverload` | `bothStrands=false` ≡ legacy overload, all '+' |
| R3 | `FindRegulatoryElements_BothStrands_NoPalindromeDuplicates_StrandSpecificPlusOnly` | AP-1/E-box/CREB once; TATA/poly(A) reverse orientation not reported |
| R4 | `OrientationIndependentRegulatoryElements_AreTheSourcedSet` | CAAT, GC, E-box, AP-1, NF-κB, CREB |
| R5 | `BucherMatrix_PssmAndThresholds_EqualBiopython` (4 matrices) | JASPAR pseudocounts, pssm max/min/consensus/mean, `threshold_fpr(1e-3/1e-4)`, `threshold_patser`, `threshold_balanced` |
| R6 | `FromCounts_JasparTataBox_EqualsBiopythonJasparPssm` | all 60 POL012.1 PSSM cells |
| R7 | `FindPromoterElementsByMatrix_EqualsBiopythonSearch` | Biopython `pssm.search` at `threshold_fpr(1e-3)` (both strands for CCAAT/GC) |
| R8 | `RegulatoryBothStrands_ReverseComplementMirror` (4 seeds) | hits(s) ↔ hits(revcomp s) mirrored |


## 9. Review 2026-09 (B05 audit group D, part 2) — σ70 promoter pairing and Promoter Calculator

Tests: `Unit/Analysis/MotifFinder_Sigma70Promoters_Tests.cs`; MCP `AlphabetPwmAndSigma70Tests` (find_sigma70_promoters, predict_sigma70_promoters).
References: Harley & Reynolds 1987 (TTGACA / TATAAT, spacer 15–21, 17 ± 1); independent Python brute force (300 random cases, identical); Promoter Calculator v1.0 reference Python (hsalis/SalisLabCode; 60 random sequences, 4634 per-TSS predictions, bit-identical).

| ID | Test | Locked values / invariant |
|----|------|---------------------------|
| S1 | `FindSigma70_LacUv5_PairsTttacaTataatWithSpacer18` | −35 TTTACA at 68 (1 mm), −10 TATAAT at 92, spacer 18, deviation 1 |
| S2 | `FindSigma70_ConsensusPromoter_SpacerBounds` | spacer 17 found; 15–16 / 22 excluded unless allowed |
| S3 | `FindSigma70_MinusStrand_ForwardCoordinatesAndOrder` | minus-strand boxes, forward coordinates |
| S4 | `FindSigma70_MismatchTolerance_EqualsBruteForceEnumeration`, `FindSigma70_Guards` | every pair within limits; contracts |
| S5 | `PromoterCalculator_LacUv5_BestForwardStateEqualsReference` | 65 + 65 TSSs; best TSS 108 TTTACA·18·TATAAT, ΔG −3.044998584365458 (all 9 terms), rate 6123.861140216731 |
| S6 | `PromoterCalculator_LacUv5_ReverseStrandAndInVitroEqualReference` | reverse best key 63 CACACA/TAAAGT ΔG −1.9519008907397968; in-vitro rate 504.40616830647446 |
| S7 | `PromoterCalculator_TooShortSequence_NoPredictions` | 77 nt → none; 78 nt → TSS 58 (+), 20 (−) |
