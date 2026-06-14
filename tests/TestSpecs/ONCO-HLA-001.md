# Test Specification: ONCO-HLA-001

**Test Unit ID:** ONCO-HLA-001
**Area:** Oncology
**Algorithm:** HLA allele nomenclature parsing/validation + allele-specific HLA LOH (LOHHLA) classification
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-15

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | WHO Nomenclature Committee — Naming Alleles (IPD-IMGT/HLA) | 2 | https://hla.alleles.org/pages/nomenclature/naming_alleles/ | 2026-06-15 |
| 2 | Marsh SGE et al. (2010), Nomenclature for factors of the HLA system, Tissue Antigens 75(4):291–455 | 1 | https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1399-0039.2010.01466.x | 2026-06-15 |
| 3 | McGranahan N et al. (2017), Allele-Specific HLA Loss…, Cell 171(6):1259–1271 (LOHHLA) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/ | 2026-06-15 |
| 4 | mskcc/lohhla — LOHHLAscript.R (reference impl) | 3 | https://raw.githubusercontent.com/mskcc/lohhla/master/LOHHLAscript.R | 2026-06-15 |

### 1.2 Key Evidence Points

1. HLA allele name = `HLA-<Gene>*F1:F2[:F3[:F4]][suffix]`; gene separated by `*`, fields by colons — Source 1, 2.
2. F1 = type/allele group; F2 = specific HLA protein/subtype; F3 = synonymous coding substitutions; F4 = non-coding differences — Source 1.
3. Minimum is two fields ("all alleles receive at least a four digit name … the first two sets of digits"); maximum four — Source 1, 2.
4. Optional trailing expression suffix ∈ {N, L, S, C, A, Q} — Source 1.
5. HLA LOH: an HLA allele with copy number **< 0.5** is "subject to loss, … indicative of LOH" — Source 3 (verbatim).
6. To avoid over-calling LOH, allelic imbalance must be significant: paired Student's t-test **p < 0.01** — Source 3 (verbatim), confirmed by paired `t.test(..., paired=TRUE)` in Source 4.

### 1.3 Documented Corner Cases

- Single-field name (`HLA-A*02`) is incomplete → invalid (two-field minimum) — Source 1, 2.
- Five+ fields → invalid (four-field maximum) — Source 1, 2.
- Trailing letter outside {N,L,S,C,A,Q} → invalid — Source 1.
- Raw allele CN < 0.5 but imbalance p ≥ 0.01 → NOT LOH (over-calling guard) — Source 3.
- Both alleles CN ≥ 0.5 → heterozygous retained, no LOH — Source 3.

### 1.4 Known Failure Modes / Pitfalls

1. Over-calling LOH from copy-number noise without the imbalance significance guard — Source 3.
2. Treating `HLA-A*02` (one field) as a complete allele — Source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ParseHlaAllele(string)` | OncologyAnalyzer | Canonical | Parse + validate WHO HLA nomenclature into gene/fields/suffix. |
| `TryParseHlaAllele(string, out HlaAllele)` | OncologyAnalyzer | Delegate | Non-throwing wrapper over `ParseHlaAllele`. |
| `DetectHlaLoh(HlaAlleleCopyNumber)` | OncologyAnalyzer | Canonical | LOHHLA classification from caller-supplied per-allele CN + imbalance p. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A parsed allele has 2 ≤ field count ≤ 4 | Yes | Source 1, 2 |
| INV-2 | Parsed gene + reconstructed `HLA-...` name round-trips the normalized input | Yes | Source 1 |
| INV-3 | HLA LOH is called ⇔ (exactly one allele CN < 0.5) ∧ (imbalance p < 0.01) | Yes | Source 3 |
| INV-4 | When LOH is called, the lost allele is the one with CN < 0.5 | Yes | Source 3 |
| INV-5 | Copy-number boundary is strict: CN = 0.5 is retained, p = 0.01 is not significant | Yes | Source 3 (verbatim `< 0.5`, `p < 0.01`) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Parse two-field | `HLA-A*02:01` | Gene A, F1=02, F2=01, no F3/F4/suffix | Source 1, 2 |
| M2 | Parse three-field | `HLA-B*07:02:01` | Gene B, F1=07, F2=02, F3=01 | Source 1 |
| M3 | Parse four-field | `HLA-C*07:02:01:03` | Gene C, fields 07/02/01/03 | Source 1 |
| M4 | Parse with suffix | `HLA-A*24:02:01:02L` | fields 24/02/01/02, suffix L (Low) | Source 1 |
| M5 | Reject missing prefix | `A*02:01` | throws FormatException | Source 1 |
| M6 | Reject single field | `HLA-A*02` | throws FormatException (two-field min) | Source 1, 2 |
| M7 | Reject five fields | `HLA-A*02:01:01:01:01` | throws FormatException (four-field max) | Source 1, 2 |
| M8 | Reject bad suffix | `HLA-A*02:01X` | throws FormatException (X∉{N,L,S,C,A,Q}) | Source 1 |
| M9 | LOH — lose allele 2 | CN(1.8, 0.30), p=0.001 | LOH=true, lost=Allele2 | Source 3 |
| M10 | LOH — lose allele 1 | CN(0.10, 1.50), p=0.0005 | LOH=true, lost=Allele1 | Source 3 |
| M11 | No LOH — both retained | CN(1.10, 0.90), p=0.30 | LOH=false | Source 3 |
| M12 | No LOH — imbalance guard | CN(1.60, 0.40), p=0.05 | LOH=false (p ≥ 0.01) | Source 3 |
| M13 | Boundary CN = 0.5 | CN(1.50, 0.50), p=0.001 | LOH=false (0.5 not < 0.5) | Source 3 |
| M14 | Boundary p = 0.01 | CN(1.70, 0.40), p=0.01 | LOH=false (0.01 not < 0.01) | Source 3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | TryParse valid | `HLA-A*02:01` | returns true, populated allele | Delegate smoke |
| S2 | TryParse invalid | `bad` | returns false, default allele | Delegate smoke |
| S3 | Null/empty parse | null / "" / "   " | throws ArgumentException-family | Validation |
| S4 | Negative CN | CN(-1, 1), p=0.001 | throws ArgumentException | Validation (CN ≥ 0) |
| S5 | p out of [0,1] | CN(1,0.3), p=1.5 | throws ArgumentException | Validation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Both alleles < 0.5 | CN(0.20, 0.30), p=0.001 | HomozygousLoss (assumption) | ASSUMPTION row |
| C2 | Lowercase gene normalized | `hla-a*02:01` | Gene A (uppercased) | Normalization |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for HLA nomenclature parsing or HLA allele-specific LOH. The only LOH tests in the suite (`OncologyAnalyzer_DetectLOH_Tests.cs`) cover genome-wide HRD-LOH (ONCO-LOH-001), a different algorithm (segment-based scarHRD), not allele-specific HLA LOH. No `HlaTyper` class exists; the checklist's `TypeHLA(bamFile)` caller is out of scope (no retrievable, well-defined behavior without a trained HLA genotyping model / reference DB).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14 | ❌ Missing | New unit; no prior tests. |
| S1–S5 | ❌ Missing | New unit. |
| C1–C2 | ❌ Missing | New unit. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_HlaAnalysis_Tests.cs` — all parsing + LOH tests for this unit.
- **Remove:** none (no pre-existing HLA tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_HlaAnalysis_Tests.cs` | Canonical fixture for ONCO-HLA-001 | 21 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | S3 | ❌ Missing | Implemented | ✅ Done |
| 18 | S4 | ❌ Missing | Implemented | ✅ Done |
| 19 | S5 | ❌ Missing | Implemented | ✅ Done |
| 20 | C1 | ❌ Missing | Implemented | ✅ Done |
| 21 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 21
**✅ Done:** 21 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1–M14 | ✅ Covered | Evidence-based exact assertions. |
| S1–S5 | ✅ Covered | Delegate smoke + validation. |
| C1–C2 | ✅ Covered | Assumption boundary + normalization. |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Both alleles CN < 0.5 with significant imbalance → reported as `HomozygousLoss`, not allele-specific LOH (source defines lost-allele as one CN < 0.5 only). | C1 |

---

## 7. Open Questions / Decisions

1. The checklist by-area definition lists `HlaTyper.TypeHLA(bamFile)` (a full HLA genotyping caller). That caller requires a trained model / IPD-IMGT/HLA reference allele database whose behavior is not retrievable as an exact specification; fabricating it would violate the evidence-first policy. Per the "external evidence wins" rule, this unit instead implements the two retrievable, formally-specified pieces (nomenclature parsing per WHO; allele-specific HLA LOH per LOHHLA). Conflict recorded; checklist canonical-method note updated accordingly.
