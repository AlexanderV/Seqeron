# Test Specification: COMPGEN-ANI-001

**Test Unit ID:** COMPGEN-ANI-001
**Area:** Comparative
**Algorithm:** Average Nucleotide Identity (ANI), ANIb (Goris et al. 2007)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-23

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Goris et al. (2007) DNA-DNA hybridization values and whole-genome similarities. IJSEM 57:81-91 | 1 | https://doi.org/10.1099/ijs.0.64483-0 | 2026-06-14 |
| 2 | Konstantinidis & Tiedje (2005) Genomic insights that advance the species definition. PNAS 102:2567-2572 | 1 | https://doi.org/10.1073/pnas.0409727102 | 2026-06-14 |
| 3 | Lee et al. (2016) OrthoANI. IJSEM 66:1100-1103 | 1 | https://doi.org/10.1099/ijsem.0.000760 | 2026-06-14 |
| 4 | pyani ANIb reference implementation (Pritchard et al.) | 3 | https://pyani.readthedocs.io/en/latest/api/pyani.anib.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. Query genome is cut into consecutive (non-overlapping) 1020 nt fragments. — Goris 2007 (Methods)
2. Each fragment's best match against the reference is found via BLASTN; identity is **recalculated over the whole fragment length**. — Goris 2007 / pyani
3. ANI = mean identity of matches with **>30 % overall identity over an alignable region of ≥70 % of their length**. — Goris 2007 (verbatim)
4. ANI ≈ 95 % (Goris) / ≈ 94 % (Konstantinidis & Tiedje) corresponds to the 70 % DDH species boundary. — sources 1, 2
5. The best BLASTN match is **gapped** (pyani `-xdrop_gap_final 150`; `ani_alnlen = blast_alnlen - blast_gaps`); identity AND coverage are recalculated over the query-fragment length (`ani_pid = ani_alnids/qlen`, `ani_coverage = ani_alnlen/qlen`). — pyani anib source (2026-06-23)
6. The search is performed in **both directions**: "reverse searching ... was also performed to provide reciprocal values"; the symmetric ANIb value is the mean of the two directions. — Goris 2007 / pyani (2026-06-23)

### 1.3 Documented Corner Cases

- Fragments below the identity or alignable-region cut-off are discarded and do not enter the mean (Goris 2007 / pyani).
- Per-fragment identity is over the fragment length, not over only the aligned sub-region (Goris 2007 "recalculated to an identity along the entire sequence").
- ANI is direction-dependent (query is the fragmented genome); pyani notes non-symmetrical matrices.

### 1.4 Known Failure Modes / Pitfalls

1. Computing identity as longest-common-substring length over fragment length (NOT nucleotide identity) — the pre-existing implementation defect this unit corrects. — Goris 2007 definition.
2. Counting all fragments instead of only conserved (qualifying) fragments overstates divergence. — Goris 2007.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction, gapped)` | ComparativeGenomics | **Canonical** | Mean per-fragment nucleotide identity under the 30 %/70 % cut-offs; ungapped (default) or gapped placement |
| `CalculateReciprocalAni(genomeA, genomeB, fragmentLength, minIdentity, minAlignableFraction, gapped)` | ComparativeGenomics | **Canonical** | Reciprocal (symmetric) ANI = mean of both directions |
| `BestUngappedFragmentMatch(...)` | ComparativeGenomics (private) | **Internal** | Tested indirectly via CalculateANI (ungapped) |
| `BestGappedFragmentMatch(...)` | ComparativeGenomics (private) | **Internal** | Smith-Waterman placement (via `SequenceAligner`); tested indirectly via CalculateANI (gapped) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | ANI is a fraction in [0, 1] | Yes | Identity = matches/fragmentLength ∈ [0,1]; mean preserves bound (Goris 2007) |
| INV-2 | Identical sequences → ANI = 1.0 | Yes | Every fragment is a perfect substring, identity 1.0 (Goris 2007) |
| INV-3 | Only fragments with identity > minIdentity AND alignable fraction ≥ minAlignableFraction contribute | Yes | Goris 2007 cut-off clause |
| INV-4 | Fragmentation is consecutive and non-overlapping; trailing partial fragment (< fragmentLength) is ignored | Yes | Goris 2007 "consecutive 1020 nt fragments" |
| INV-5 | No qualifying fragment / empty / null input → 0 | Yes | Definition (mean over empty set undefined → 0) |
| INV-6 | Gapped placement identity ≥ ungapped identity for the same fragment (gaps recover indels) | Yes | Goris gapped BLASTN; pyani `ani_alnlen = blast_alnlen - blast_gaps` |
| INV-7 | Reciprocal ANI is symmetric: ANI(A,B) = ANI(B,A) | Yes | Mean of both directions is order-independent (Goris reverse searching) |
| INV-8 | Reciprocal ANI = (ANI(A→B) + ANI(B→A)) / 2 | Yes | Goris "reverse searching ... reciprocal values" |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Identical genomes | query = reference, fragLen 4 | ANI = 1.0 | INV-2, Goris 2007 |
| M2 | One substituted base | query `AAAACCCCGGGGTTTA` vs R, fragLen 4; last frag `TTTA`=3/4 | ANI = 0.9375 | Goris 2007 recalculated identity |
| M3 | Half-identity fragment | query `AAAACCCCGGGGAATT` vs R; last frag `AATT`=2/4 | ANI = 0.875 | Goris 2007 mean identity |
| M4 | Identity cut-off excludes fragment | query `AAAACGTC`, ref `AAAAAAAA`, fragLen 4; frag2=0.0 not >0.30 | ANI = 1.0 (only frag1) | Goris 2007 ">30 % identity" |
| M5 | Alignable cut-off excludes fragment | query `AAAA`, ref `AA` (ref < frag) | ANI = 0 | Goris 2007 "≥70 % alignable" |
| M6 | Consecutive non-overlapping fragmentation | query of length 10, fragLen 4 → 2 fragments, trailing 2 nt ignored | uses exactly 2 fragments | INV-4, Goris 2007 |
| M7 | Null / empty inputs | null or "" for either genome | ANI = 0 | Validation contract |
| M8 | Non-positive fragmentLength | fragmentLength = 0 | ArgumentOutOfRangeException | Validation contract |
| G1 | Gapped identical genomes | query = reference, fragLen 4, gapped | ANI = 1.0 | INV-02, Goris 2007 |
| G2 | Gapped recovers indel | query `AAAACCCC`, ref `AAAATCCCC`, fragLen 8 | gapped = 1.0 > ungapped = 0.875 | INV-06; gapped dataset, Goris/pyani |
| G3 | Gapped alignable cut-off | query `AAAACCCC`, ref `AA`, fragLen 8, gapped | ANI = 0 (coverage 0.25 < 0.70) | Goris ≥70 % alignable |
| G4 | Gapped identity cut-off | query `AAAACGTC`, ref `AAAAAAAA`, fragLen 4, gapped | ANI = 1.0 (CGTC excluded) | Goris >30 % identity |
| R1 | Reciprocal identical genomes | A = B, fragLen 4 | reciprocal = 1.0 | INV-08, Goris reverse searching |
| R2 | Reciprocal symmetry | A `AAAACCCCGGGGTTTT`, B `...TTTA` | ANI(A,B) = ANI(B,A) | INV-07 |
| R3 | Reciprocal = mean of directions | A `AAAACGTC`, B `AAAAAAAA`, fragLen 4 | reciprocal = (1.0+1.0)/2 = 1.0 | INV-08, reciprocal dataset |
| R4 | Reciprocal null / empty | null or "" | reciprocal = 0 | Validation contract |
| R5 | Reciprocal non-positive fragmentLength | fragmentLength = 0 | ArgumentOutOfRangeException | Validation contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Range invariant | ANI of arbitrary divergent pair | 0 ≤ ANI ≤ 1 | INV-1 |
| S2 | Query shorter than fragment | query length 3, fragLen 4 | ANI = 0 | No fragment fits |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Custom minIdentity keeps low-identity frag | lower minIdentity below a frag's identity | frag now contributes | Parameter exposure |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Original unit (2026-06-14): no existing test file; created `ComparativeGenomics_CalculateANI_Tests.cs` (M1–M8, S1, S2, C1).
- Limitation fix (2026-06-23): added gapped placement (`gapped` param) and `CalculateReciprocalAni`. New cases G1–G4 (gapped) and R1–R5 (reciprocal) added to the same canonical file.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M8 | ✅ Covered | Ungapped ANI (pre-existing) |
| S1, S2, C1 | ✅ Covered | Pre-existing edge/parameter cases |
| G1–G4 | ❌ Missing | Gapped placement cases (this fix) |
| R1–R5 | ❌ Missing | Reciprocal ANI cases (this fix) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs` — all cases for this unit.
- **Remove:** nothing (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| ComparativeGenomics_CalculateANI_Tests.cs | Canonical | 20 (11 ungapped/edge + 4 gapped + 5 reciprocal) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented identical-genomes test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented one-mismatch (0.9375) test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented half-identity (0.875) test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented identity-cutoff exclusion test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented alignable-cutoff exclusion test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented non-overlapping fragmentation test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented null/empty tests | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented invalid fragmentLength test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented range invariant property test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented query-shorter-than-fragment test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented custom minIdentity test | ✅ Done |
| 12 | G1 | ❌ Missing | Gapped identical genomes → 1.0 | ✅ Done |
| 13 | G2 | ❌ Missing | Gapped recovers indel (1.0 > ungapped 0.875) | ✅ Done |
| 14 | G3 | ❌ Missing | Gapped alignable cut-off → 0 | ✅ Done |
| 15 | G4 | ❌ Missing | Gapped identity cut-off excludes fragment | ✅ Done |
| 16 | R1 | ❌ Missing | Reciprocal identical → 1.0 | ✅ Done |
| 17 | R2 | ❌ Missing | Reciprocal symmetry | ✅ Done |
| 18 | R3 | ❌ Missing | Reciprocal = mean of directions | ✅ Done |
| 19 | R4 | ❌ Missing | Reciprocal null/empty → 0 | ✅ Done |
| 20 | R5 | ❌ Missing | Reciprocal non-positive fragmentLength throws | ✅ Done |

**Total items:** 20 (9 new in this fix)
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact value 1.0 |
| M2 | ✅ Covered | Exact value 0.9375 |
| M3 | ✅ Covered | Exact value 0.875 |
| M4 | ✅ Covered | Exact value 1.0 with frag2 excluded |
| M5 | ✅ Covered | Exact value 0.0 |
| M6 | ✅ Covered | 2 fragments used |
| M7 | ✅ Covered | 0 for null/empty |
| M8 | ✅ Covered | Throws ArgumentOutOfRangeException |
| S1 | ✅ Covered | 0 ≤ ANI ≤ 1 |
| S2 | ✅ Covered | 0.0 |
| C1 | ✅ Covered | Custom minIdentity |
| G1 | ✅ Covered | Gapped 1.0 |
| G2 | ✅ Covered | Gapped 1.0 > ungapped 0.875 (indel) |
| G3 | ✅ Covered | Gapped coverage cut-off → 0 |
| G4 | ✅ Covered | Gapped identity cut-off → 1.0 |
| R1 | ✅ Covered | Reciprocal 1.0 |
| R2 | ✅ Covered | Reciprocal symmetric |
| R3 | ✅ Covered | Reciprocal = mean of directions |
| R4 | ✅ Covered | Reciprocal 0 for null/empty |
| R5 | ✅ Covered | Reciprocal throws on fragmentLength 0 |

**In-scope cases:** 20 — **✅:** 20

---

## 6. Assumption Register

**Total assumptions:** 0 unresolved (the prior ungapped assumption is RESOLVED).

| # | Item | Type | Used In |
|---|------|------|---------|
| 1 | Gapped path uses the library Smith-Waterman aligner (`SequenceAligner.LocalAlign`, BLAST DNA scoring), not the NCBI BLASTN engine. Full DP (more sensitive than BLAST heuristic); identity/coverage recalculated over the fragment length per pyani. Not correctness-affecting for the formal definition; real-genome numerics may differ slightly by engine. | Decision (resolved) | Implementation §5.3; G1–G4 |

---

## 7. Open Questions / Decisions

1. **RESOLVED (this fix):** gapped per-fragment placement (`gapped: true`) and reciprocal/averaged ANI (`CalculateReciprocalAni`) are now implemented per Goris 2007 (gapped BLASTN; reverse searching for reciprocal values) and pyani (identity/coverage over query length). The earlier single-direction/ungapped limitation no longer applies.
2. The OrthoANI orthologous-best-pair variant (both genomes fragmented, intersected) remains out of scope; users needing it should use OrthoANI/FastANI (algorithm doc §5.3 "Not implemented").
