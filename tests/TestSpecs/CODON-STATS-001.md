# Test Specification: CODON-STATS-001

**Test Unit ID:** CODON-STATS-001
**Area:** Codon
**Algorithm:** Codon Usage Statistics (GetStatistics, CalculateCai, EColiOptimalCodons, HumanOptimalCodons)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-28

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sharp & Li (1987), CAI, NAR 15:1281–1295 | 1 | https://doi.org/10.1093/nar/15.3.1281 | 2026-06-13 |
| 2 | Peden (1999) CodonW thesis (GC3s definition §1.8.2.1.3) | 1–3 | https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf | 2026-06-13 |
| 3 | CodonW codon-usage indices (CAI exclusions) | 3 | https://codonw.sourceforge.net/Indices.html | 2026-06-13 |
| 4 | seqinr `cai` documentation | 3 | https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html | 2026-06-13 |
| 5 | EMBOSS `cusp` (GC1/GC2/GC3) | 3 | https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp | 2026-06-13 |
| 6 | Biopython v1.79 `SharpEcoliIndex` | 3 | https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py | 2026-06-13 |
| 7 | Kazusa H. sapiens [gbpri] codon usage | 5 | https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606 | 2026-06-13 |
| 8 | CodonW 1.4.4 source (`codon_us.c` `gc_out`, `how_synon`) + compiled binary | 3 | https://codonw.sourceforge.net | 2026-09-28 |
| 9 | EMBOSS 6.6 `ajcod.c` (`ajCodSetTripletsS`, `ajCodWrite`) | 3 | EMBOSS source | 2026-09-28 |
| 10 | CodonW 1.4.4 `codon_us.c` `fop_out` / `cbi_out`, `codonW.h` `fop[]` (8 built-in optimal-codon sets), `commline.c` (`-fop -cbi -f_type -fop_file -cbi_file -code`) + compiled binary | 3 | Ubuntu `codonw_1.4.4.orig.tar.gz` | 2026-10-14 |
| 11 | Ikemura (1981) J. Mol. Biol. 151:389–409 (Fop) | 1 | https://doi.org/10.1016/0022-2836(81)90003-6 | 2026-10-14 |
| 12 | Bennetzen & Hall (1982) J. Biol. Chem. 257:3026–3031 (CBI) | 1 | https://doi.org/10.1016/S0021-9258(19)81068-2 | 2026-10-14 |

### 1.2 Key Evidence Points

1. `w_i = f_i / max(f_j)` over the synonymous family; `CAI = exp[(1/L) Σ ln w_i]` — Sharp & Li 1987 / Wikipedia.
2. CAI excludes non-synonymous (single-codon Met, Trp) and termination codons — seqinr; CodonW Indices.
3. GC3s = "frequency of G or C at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)" — Peden §1.8.2.1.3.
4. GC1/GC2/GC3 = GC content at codon positions 1/2/3 ("1st/2nd/3rd letter GC") — EMBOSS cusp.
5. E. coli reference w values (CTG=1, GCC=0.122, CGT=1, AGG=0.002 …) — Biopython SharpEcoliIndex.
6. Human RSCU derived from Kazusa per-thousand frequencies (CTG≈2.3713, GCC≈1.5988, Met/Trp=1.0).
7. CodonW Fop = opt/(opt+common+rare) and CBI = (opt − Σ naa/da)/(tot − Σ naa/da) (Σ over optimal codons) over amino acids with > 1 codon in the code, not stop, with ≥ 1 optimal codon (`has_opt_info`); denominator 0 → prints 0.000; modified Fop (opt − rare)/total offered interactively, not clamped — CodonW 1.4.4 source [10].

### 1.3 Documented Corner Cases

- A sequence of only Met/Trp/stop codons has no scorable codon → CAI 0 and GC3s 0 (empty denominator). Source: seqinr / Peden.
- Empty / null input → zeroed statistics or `ArgumentNullException` (implementation contract; not specified by the literature).

### 1.4 Known Failure Modes / Pitfalls

1. Treating GC3s as GC at *all* third positions (including Met/Trp/stop) — contradicts Peden §1.8.2.1.3.
2. Including Met/Trp/stop codons in CAI — contradicts Sharp & Li 1987 / CodonW.
3. `ln(0)` when a codon's w is 0 — w < 0.0001 is scored as 0.01 (CodonW `cai_out`, seqinr `cai`, Bulmer 1988). *(Review 2026-09, CODON-CAI-001 F13: zero-w codons used to be skipped, which raised CAI.)*

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetStatistics(string / DnaSequence [, GeneticCode])` | CodonUsageAnalyzer | Canonical | Aggregates counts, RSCU, ENC, GC1/2/3, GC3s, total codons; code default NCBI 1 |
| `GetStatistics(string / DnaSequence, GeneticCode, CodonGcConvention)` | CodonUsageAnalyzer | Canonical (GC convention) | `EmbossCusp` (= the overloads above) or `CodonW` (GC1/2/3 over sense codons, stops of the code excluded — CodonW `gc_out`); review 2026-09 F37 |
| `CodonUsageStatistics.IsEncCalculable` | record | Derived | `Enc > 0`; false where CodonW prints Nc `*****` (F36) |
| `CalculateCai(string, refRscu)` / `CalculateCai(DnaSequence, refRscu)` | CodonUsageAnalyzer | Canonical | Sharp & Li 1987 CAI |
| `EColiOptimalCodons` (property) | CodonUsageAnalyzer | Reference | Sharp & Li 1987 w table |
| `HumanOptimalCodons` (property) | CodonUsageAnalyzer | Reference | Kazusa-derived human RSCU |
| `CodonUsageStatistics.OverallGc` | record | Internal | (GC1+GC2+GC3)/3 |
| `CalculateFop` / `TryCalculateFop(string / DnaSequence, …, optimalCodons, code, rareCodons)` | CodonUsageAnalyzer | Canonical (A3-13) | CodonW 1.4.4 `-fop` (null set = E. coli; `rareCodons` = modified Fop) |
| `CalculateCbi` / `TryCalculateCbi(string / DnaSequence, …, optimalCodons, code)` | CodonUsageAnalyzer | Canonical (A3-13) | CodonW 1.4.4 `-cbi` |
| `GetStatistics(string, GeneticCode, CodonGcConvention, optimalCodons)` → `Fop` / `Cbi` | CodonUsageAnalyzer | Opt-in (A3-13) | null when not calculable / not requested |
| `CodonWOptimalCodonSets` (`All`, 8 sets, `TryGetByName`) | static class | Reference (A3-13) | CodonW `fop[]`, `-f_type 0..7` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `0 ≤ CalculateCai ≤ 1` for any sequence and reference | Yes | Sharp & Li 1987 (geometric mean of w∈[0,1]) |
| INV-2 | CAI of an all-optimal sequence (each codon = family's w-max) = 1.0 | Yes | Sharp & Li 1987 |
| INV-3 | GC3s ignores stop codons and single-codon amino acids of the genetic code (Met ATG, Trp TGG in table 1); denominator counts only synonymous codons | Yes | Peden §1.8.2.1.3; CodonW `gc_out` |
| INV-4 | `0 ≤ GC1, GC2, GC3, GC3s ≤ 100` | Yes | counts/positions ratio ×100 |
| INV-5 | `TotalCodons` = number of valid (ACGT-only) codons in frame | Yes | EMBOSS cusp "Number" |
| INV-6 | `OverallGc` = (GC1+GC2+GC3)/3 | Yes | record definition |
| INV-7 | `CodonW` convention: GC1/2/3 denominators = sense codons (CodonW `Len_aa`); `OverallGc` = CodonW `GC`; equal to `EmbossCusp` for stop-free genes | Yes | CodonW 1.4.4 `gc_out` (`ca == 11` skipped) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CAI all-optimal E. coli | `CTGATCGTTGCTCGTAAA` vs `EColiOptimalCodons` | 1.0 (Within 1e-10) | Sharp & Li 1987 |
| M2 | CAI geometric-mean derivation | `GCTGCC` vs E. coli (GCT w=1, GCC w=0.122) | 0.34928498393146 | Sharp & Li 1987 √(1×0.122) |
| M3 | CAI suboptimal derivation | `CTAATAGTC` vs E. coli (0.007,0.003,0.066) | 0.01114947479545 | Sharp & Li 1987 ∛(·) |
| M4 | CAI excludes Met/Trp/stop | `ATGTGGTAA` vs E. coli | 0.0 (no scorable codon) | seqinr / CodonW |
| M5 | GC3s excludes Met/Trp/stop | `ATGGCA`: GC3s vs GC3 | GC3s=0.0, GC3=50.0 | Peden §1.8.2.1.3 |
| M6 | GC3s synonymous fraction | `GCCGCA` (2 Ala; 3rd C,A) | GC3s=50.0 | Peden §1.8.2.1.3 |
| M7 | GC1/GC2/GC3 per position | `CTGGTTAAA` | 66.6667 / 0.0 / 33.3333 | EMBOSS cusp |
| M8 | TotalCodons & counts | `CTGCTGGTTAAA` → 4 codons, CTG=2 | TotalCodons=4, counts["CTG"]=2 | EMBOSS cusp |
| M9 | E. coli table values | `EColiOptimalCodons` | CTG=1.0, GCC=0.122, CGT=1.0, AGG=0.002, TTT=0.296 | Biopython SharpEcoliIndex |
| M10 | Human table values | `HumanOptimalCodons` | CTG=2.3713, GCC=1.5988, ATG=1.0, GTG=1.8517 | Kazusa-derived RSCU |
| M11 | RSCU populated in stats | `GetStatistics("CTGCTG")` Rscu["CTG"] | =6.0 (only Leu codon present, 6-fold) | Sharp et al. 1986 RSCU |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty string stats | `GetStatistics("")` | all-zero `CodonUsageStatistics` | input contract |
| S2 | Empty string CAI | `CalculateCai("", EColiOptimalCodons)` | 0.0 | input contract |
| S3 | Null DnaSequence stats | `GetStatistics((DnaSequence)null)` | `ArgumentNullException` | input contract |
| S4 | Null DnaSequence CAI | `CalculateCai((DnaSequence)null, ref)` | `ArgumentNullException` | input contract |
| S5 | Null reference CAI | `CalculateCai(seq, null)` | `ArgumentNullException` | input contract |
| S6 | CAI in [0,1] | arbitrary CDS vs E. coli | 0 ≤ CAI ≤ 1 | INV-1 |
| S7 | Trailing partial codon ignored | `GCCGCAG` (7 nt) | only 2 codons counted | in-frame parsing |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | OverallGc average | `GetStatistics("CTGGTTAAA").OverallGc` | (66.6667+0+33.3333)/3 = 33.3333 | INV-6 |
| C2 | Lowercase normalized | `GetStatistics("ctg")` | same as "CTG" | case-insensitive |
| R1 | GC3s genetic-code aware (review 2026-09 F18) | `ATAATG`, `TGATGGGCC`, `AGAAGGCTG`, `augcccuuaugg` × NCBI 1/2 | CodonW -gc3s ×100: 0/50, 100/66.67, 66.67/100, 50/75 | CodonW 1.4.4 binary |
| R2 | RSCU/ENC in stats use the same code | NCBI 2 | = `CalculateRscu/CalculateEnc(…, code)` | consistency |
| R3 | Stops + ambiguous triplet | `ATGGCAGCCTAAtggNNNgcg` | Total 6, GC1 50, GC2/GC3 66.67 (cusp, stops included), GC3s 66.67 (CodonW 0.667) | EMBOSS ajcod.c; CodonW |
| R4 | DnaSequence overload / null code | NCBI 2 | same as string; `ArgumentNullException` | input contract |
| R11 | CodonW convention, terminal stop (F37) | `ATGAAACTGTTCGCCAAATAA`, NCBI 1 | GC1/2/3 = 2/6, 1/6, 4/6 (×100), OverallGc 7/18 (CodonW GC 0.389), GC3s 60; cusp OverallGc 7/21, TotalCodons 7 both | CodonW 1.4.4 binary `-base -gc -gc3s` |
| R12 | CodonW convention, random gene r53 (N triplet, in-frame stops, trailing base) | NCBI 1 / 2 | Len_aa 48/46; GC1/2/3 = 27,25,23 / 27,23,23 over Len_aa; CodonW GC 0.521/0.529, GC3s 0.479/0.500 | CodonW `-code 0/1` |
| R13 | CodonW convention, internal stops + NNN | `ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA`, NCBI 1 / 2 | Len_aa 10/11; GC 0.600/0.576; GC1 0.600/0.545; GC2 0.500/0.545; GC3 0.700/0.636; GC3s 0.625/0.636; TotalCodons 14 | CodonW `-code 0/1` |
| R14 | Default = EMBOSS cusp (re-verified on the binary) | `ATGGCTTAAGCCGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA` | Total 14; Coding GC 47.62, GC1 42.86, GC2 42.86, GC3 57.14 | EMBOSS 6.6.0.0 `cusp` |
| R15 | Conventions agree on a stop-free gene; stops-only gene → CodonW-mode GC 0 (CodonW writes no value) | M3 gene; `TAATAGTGA` | equal; 0 | CodonW `gc_out` "too short" |
| R16 | Contracts: DnaSequence overload; undefined enum → `ArgumentOutOfRangeException`; null code/sequence; MCP `gc_convention` (default EmbossCusp), unknown value → `ArgumentException` | — | — | input contract |
| R17 | Fop/CBI vs CodonW binary (A3-13) | 70 genes × 8 built-in `-f_type` sets + 3 user `-fop_file` sets × `-code 0..7` (NCBI 1,2,3,4,5,6,10,9), incl. modified Fop (B. subtilis + user1) — `TestData/Codon/codonw_fop_cbi.tsv` | 7 280 Fop + 6 160 CBI values within ±0.0005 (CodonW %5.3f); 0 mismatches | CodonW 1.4.4 binary [10] |
| R18 | Hand-derived Fop/CBI | `ATGCTGCTGCTGCTGAAAAAGTTTTTCGGCGGTGGATAA`, default E. coli | Fop 8/11, CBI 23/41 (N_ran 25/6); CodonW 0.727 / 0.561 | [10][11][12] |
| R19 | Families follow the code; other set | same gene, NCBI 2, B. subtilis | Fop 3/11, CBI −0.055 (CodonW `-f_type 1 -code 1`) | CodonW binary |
| R20 | Not calculable | Met/Trp/stop only, stops only, `ATGNNN`, empty | Try false, value 0 (CodonW 0.000) | `fop_out`/`cbi_out` |
| R21 | CBI with only-optimal families | Phe {TTT,TTC} both optimal, `TTTTTCTTT` | CBI Try false (0/0, CodonW 0.000), Fop 1 | `cbi_out` |
| R22 | Modified Fop negative (not clamped) | oracle B. subtilis rows | value < 0 equals CodonW; original Fop ≥ 0 | CodonW binary |
| R23 | Built-in sets + lookup + validation | 8 sets (24/19/15/19/22/20/20/22 optimal; B. subtilis 10 rare); name lookup; RNA/lower-case codons; invalid codon or optimal∩rare → `ArgumentException` | — | codonW.h `fop[]` |
| R24 | Opt-in in statistics + MCP | `GetStatistics(…, optimalCodons)` Fop/Cbi; other overloads null; MCP `include_fop_cbi`, `fop_cbi_organism`, `optimal_codons` | — | additive API |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzerTests.cs` — legacy fixture containing weak/permissive tests for `GetStatistics`, `CalculateCai`, `EColiOptimalCodons`, `HumanOptimalCodons` (regions "CAI", "GetStatistics", "Optimal Codon Tables", "Edge Cases"). Also contains tests for `CountCodons`, `CalculateRscu`, `CalculateEnc` belonging to **other** units (CODON-RSCU-001, CODON-ENC-001).
- `tests/Seqeron/Seqeron.Genomics.Tests/Properties/CodonProperties.cs` — property `INV-CAI3` calls `CalculateCai` (bounds in [0,1]); kept (other unit's property file).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| CAI all-optimal (M1) | ⚠ Weak | legacy `CalculateCai_OptimalCodons_ReturnsHigh` used `Within(0.01)` and a non-derived sequence |
| CAI geo-mean (M2) | ❌ Missing | no exact-derivation test existed |
| CAI suboptimal (M3) | ⚠ Weak | legacy used `Is.LessThan(0.5)` (permissive) |
| CAI excludes Met/Trp/stop (M4) | ❌ Missing | not covered (was a defect) |
| GC3s excl. Met/Trp/stop (M5) | ❌ Missing | not covered (was a defect: GC3s==GC3) |
| GC3s synonymous fraction (M6) | ❌ Missing | not covered |
| GC1/GC2/GC3 (M7) | ⚠ Weak | legacy used `Is.GreaterThan(0)` |
| TotalCodons/counts (M8) | ⚠ Weak | legacy checked `Not.Empty`, `EqualTo(5)` on undocumented seq |
| E. coli table (M9) | ⚠ Weak | legacy checked `Count==64` and ordering only |
| Human table (M10) | ⚠ Weak | legacy checked `Count==64` only |
| RSCU in stats (M11) | ❌ Missing | not covered |
| Empty stats (S1) | ⚠ Weak | legacy checked TotalCodons/Enc only |
| Empty CAI (S2) | ✅→re-do | legacy `EqualTo(0)`; re-implemented in canonical file |
| Null stats/CAI/ref (S3–S5) | ✅→re-do | legacy present; re-implemented in canonical file |
| CAI in [0,1] (S6) | ⚠ Weak | legacy range test |
| trailing partial codon (S7) | ❌ Missing | not covered |
| OverallGc (C1) | ⚠ Weak | legacy `Within(0.01)` ok but moved to canonical |
| lowercase (C2) | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_GetStatistics_Tests.cs` — all CODON-STATS-001 cases (GetStatistics, CalculateCai, EColiOptimalCodons, HumanOptimalCodons).
- **Remove:** the `GetStatistics`/`CAI`/`Optimal Codon Tables`/relevant Edge-Case tests from legacy `CodonUsageAnalyzerTests.cs` (superseded). Keep the legacy `CountCodons`/`CalculateRscu`/`CalculateEnc` tests (other units).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `CodonUsageAnalyzer_GetStatistics_Tests.cs` | canonical CODON-STATS-001 | 21 |
| `CodonUsageAnalyzerTests.cs` | legacy; other-unit tests only | (reduced) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with derived seq + Within(1e-10) | ✅ Done |
| 2 | M2 | ❌ Missing | implemented √0.122 derivation | ✅ Done |
| 3 | M3 | ⚠ Weak | rewrote with ∛(·) exact value | ✅ Done |
| 4 | M4 | ❌ Missing | implemented exclusion test | ✅ Done |
| 5 | M5 | ❌ Missing | implemented GC3s≠GC3 test | ✅ Done |
| 6 | M6 | ❌ Missing | implemented GC3s synonymous fraction | ✅ Done |
| 7 | M7 | ⚠ Weak | rewrote with exact per-position values | ✅ Done |
| 8 | M8 | ⚠ Weak | rewrote with derived counts | ✅ Done |
| 9 | M9 | ⚠ Weak | rewrote with exact Sharp&Li w values | ✅ Done |
| 10 | M10 | ⚠ Weak | rewrote with exact Kazusa RSCU values | ✅ Done |
| 11 | M11 | ❌ Missing | implemented RSCU-in-stats test | ✅ Done |
| 12 | S1 | ⚠ Weak | rewrote full zeroed-struct check | ✅ Done |
| 13 | S2 | ✅ | re-implemented in canonical file | ✅ Done |
| 14 | S3 | ✅ | re-implemented in canonical file | ✅ Done |
| 15 | S4 | ✅ | re-implemented in canonical file | ✅ Done |
| 16 | S5 | ✅ | re-implemented in canonical file | ✅ Done |
| 17 | S6 | ⚠ Weak | kept as bounds property test | ✅ Done |
| 18 | S7 | ❌ Missing | implemented partial-codon test | ✅ Done |
| 19 | C1 | ⚠ Weak | rewrote with exact value | ✅ Done |
| 20 | C2 | ❌ Missing | implemented lowercase test | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 1.0 Within(1e-10) |
| M2 | ✅ | √0.122 exact |
| M3 | ✅ | ∛(·) exact |
| M4 | ✅ | exclusion → 0.0 |
| M5 | ✅ | GC3s=0, GC3=50 |
| M6 | ✅ | GC3s=50 |
| M7 | ✅ | exact per-position |
| M8 | ✅ | counts/total exact |
| M9 | ✅ | Sharp&Li w exact |
| M10 | ✅ | Kazusa RSCU exact |
| M11 | ✅ | RSCU["CTG"]=6.0 |
| S1 | ✅ | full zeroed struct |
| S2 | ✅ | 0.0 |
| S3 | ✅ | throws |
| S4 | ✅ | throws |
| S5 | ✅ | throws |
| S6 | ✅ | bounds |
| S7 | ✅ | partial codon ignored |
| C1 | ✅ | exact average |
| C2 | ✅ | lowercase normalized |

All in-scope cases ✅. Count of ✅ = 20 = total in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | GC3s reported as a percentage (×100) for consistency with GC1/GC2/GC3 (CodonW uses a fraction). | M5, M6 expected values |
| 2 | ~~Zero-w codons are skipped~~ — resolved 2026-09 (CODON-CAI-001 F13): w < 0.0001 → 0.01 as CodonW/seqinr. | — |

---

## 7. Open Questions / Decisions

1. **Decision:** `EColiOptimalCodons` replaced with Sharp & Li 1987 w-values (Biopython SharpEcoliIndex) and `HumanOptimalCodons` with Kazusa-derived RSCU; prior values were untraceable (defect). CAI rescales by family max, so passing w (max 1.0) reproduces w.
2. **Decision:** Fixed two implementation defects — GC3s now restricted to synonymous third positions (Peden), and CAI now excludes Met/Trp/stop (seqinr/CodonW).
3. **Decision (review 2026-09, F37, audit A1-14):** the default GC1/GC2/GC3/OverallGc stay EMBOSS cusp (stops included); an additive `CodonGcConvention.CodonW` option gives CodonW's numbers (stop codons of the genetic code excluded, CodonW `gc_out`). Differential vs the CodonW 1.4.4 binary: 304 genes (in-frame stops, N/R/Y, lower-case, RNA, trailing bases) × CodonW codes 0–7 — GC, GC1, GC2, GC3, GC3s all equal at 3 dp on every gene with output (2220 of the 2432 gene×code cases differ from the cusp values).
4. **Note (2026-09, F37):** the EMBOSS 6.6 `cusp` *binary* counts an ambiguous triplet as a codon (`NNN` appears as `AAA`; `ajCodBase` maps N to its lowest constituent) and so includes it in the GC denominators (R3 gene: GC1 42.86 % on 7 codons); this artefact contradicts the documented `ajCodSetTripletsS` contract and is not reproduced (same decision as CODON-USAGE-001 F9). On ACGT genes the default equals the binary (R14).
5. **Decision (review 2026-10, A3-13):** Fop and CBI are implemented as opt-in APIs exactly as CodonW 1.4.4 (Ikemura 1981 / Bennetzen & Hall 1982 as implemented by Peden); CodonW's 8 built-in optimal-codon sets are exposed as `CodonWOptimalCodonSets`. CodonW computes in single precision and prints 3 dp; values are locked to ±0.0005. CodonW prints 0.000 for a zero denominator; the Try APIs report it as not calculable. The modified Fop is reproduced unclamped (CodonW's comment says negatives are adjusted to 0, the code and binary do not).
6. None open.
