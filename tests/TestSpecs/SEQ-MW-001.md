# Test Specification: SEQ-MW-001

**Test Unit ID:** SEQ-MW-001
**Area:** Statistics
**Algorithm:** Molecular Weight Calculation (protein and nucleotide)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Expasy Compute pI/Mw documentation | 2 | https://web.expasy.org/compute_pi/pi_tool-doc.html | 2026-06-13 |
| 2 | Expasy FindMod — average residue masses | 2 | https://web.expasy.org/findmod/findmod_masses.html | 2026-06-13 |
| 3 | Expasy ProtParam documentation | 2 | https://web.expasy.org/protparam/protparam-doc.html | 2026-06-13 |
| 4 | Biopython `Bio/SeqUtils/__init__.py` (`molecular_weight`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-13 |
| 5 | Biopython `Bio/Data/IUPACData.py` (weight tables) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py | 2026-06-13 |

### 1.2 Key Evidence Points

1. Protein Mw = sum of average amino-acid masses + average mass of one water molecule. — Expasy Compute pI/Mw (source 1).
2. Reference formula (all sequence types): `weight = sum(table[x]) − (len−1)·water`, water_avg = 18.0153 Da. — Biopython `molecular_weight` (source 4).
3. Average protein residue/free-amino-acid masses and DNA/RNA monophosphate masses are fixed tables. — Biopython IUPACData (source 5); Expasy FindMod (source 2).
4. Worked examples: AGC protein = 249.2874 (≈249.29), DNA = 949.6095 (≈949.61), RNA = 997.6077 (≈997.61). — Biopython docstring (source 4).
5. Nucleotide masses are monophosphates with an assumed 5' phosphate; one water removed per phosphodiester bond. — Biopython (sources 4, 5).

### 1.3 Documented Corner Cases

- Single monomer ⇒ zero bonds ⇒ result is the free monomer mass (Expasy formula; Biopython tables).
- Biopython allows "only unambiguous letters"; unknown letters are an error there. Seqeron resolves unknowns by skipping them (no invented mass) — recorded as ASSUMPTION-02.
- Protein has no double-stranded mode (Biopython raises).

### 1.4 Known Failure Modes / Pitfalls

1. Forgetting the `(len−1)·water` correction → over-counts nucleotide mass by one water per bond. — Biopython formula (source 4).
2. Using free-amino-acid masses without subtracting peptide-bond water → over-counts protein mass. — Expasy formula (source 1).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateMolecularWeight(string proteinSequence)` | SequenceStatistics | **Canonical** | Average-mass protein Mw |
| `CalculateNucleotideMolecularWeight(string sequence, bool isDna = true)` | SequenceStatistics | **Canonical** | Average-mass DNA/RNA Mw |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | MW(empty)=0 and MW(null)=0 | Yes | Implementation contract; sources define ≥1 monomer |
| INV-2 | MW > 0 for any non-empty sequence of recognized monomers | Yes | All table masses > water (sources 2,5) |
| INV-3 | For protein, MW(2 residues) = m1 + m2 − water (exactly one peptide-bond water removed) | Yes | Biopython formula (source 4) |
| INV-4 | For nucleotides, MW(2 monomers) = m1 + m2 − water (exactly one phosphodiester-bond water removed) | Yes | Biopython formula (source 4) |
| INV-5 | Case-insensitive: MW(seq) = MW(uppercase(seq)) | Yes | Sibling-method convention (`ToUpperInvariant`) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Protein AGC | Average-mass protein Mw of "AGC" | 249.2874 Da (`.Within(1e-3)`; docstring 249.29) | Sources 1,4,5 |
| M2 | DNA AGC | Average-mass DNA Mw of "AGC" | 949.6095 Da (docstring 949.61) | Sources 4,5 |
| M3 | RNA AGC | Average-mass RNA Mw of "AGC" | 997.6077 Da (docstring 997.61) | Sources 4,5 |
| M4 | Single Gly | Protein Mw of "G" (zero bonds) | 75.0666 Da | Sources 1,5 |
| M5 | Single DNA A | DNA Mw of "A" (zero bonds) | 331.2218 Da | Source 5 |
| M6 | Single RNA A | RNA Mw of "A" (zero bonds) | 347.2212 Da | Source 5 |
| M7 | Empty protein | `CalculateMolecularWeight("")` | 0 | INV-1 |
| M8 | Null protein | `CalculateMolecularWeight(null)` | 0 | INV-1 |
| M9 | Empty nucleotide | `CalculateNucleotideMolecularWeight("")` | 0 | INV-1 |
| M10 | Null nucleotide | `CalculateNucleotideMolecularWeight(null)` | 0 | INV-1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive protein | MW("agc") == MW("AGC") | 249.2874 Da | INV-5 |
| S2 | Case-insensitive DNA | MW("agc") == MW("AGC") | 949.6095 Da | INV-5 |
| S3 | Unknown protein symbol skipped | MW("AG*C") == MW("AGC") | 249.2874 Da | ASSUMPTION-02 |
| S4 | Unknown nucleotide symbol skipped | MW("AG*C", DNA) == MW("AGC", DNA) | 949.6095 Da | ASSUMPTION-02 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Bond-count invariant (protein) | MW("AG") = m_A + m_G − water | 89.0932+75.0666−18.0153 = 146.1445 Da | INV-3 |
| C2 | Bond-count invariant (DNA) | MW("AG", DNA) = m_A + m_G − water | 331.2218+347.2212−18.0153 = 660.4277 Da | INV-4 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatisticsTests.cs` — legacy fixture; grepped for MW coverage.
- No `{Class}_{Method}` canonical file exists for these two methods.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 protein AGC | ❌ Missing | New canonical file |
| M2 DNA AGC | ❌ Missing | New |
| M3 RNA AGC | ❌ Missing | New |
| M4 single Gly | ❌ Missing | New |
| M5 single DNA A | ❌ Missing | New |
| M6 single RNA A | ❌ Missing | New |
| M7 empty protein | ❌ Missing | New |
| M8 null protein | ❌ Missing | New |
| M9 empty nucleotide | ❌ Missing | New |
| M10 null nucleotide | ❌ Missing | New |
| S1 case protein | ❌ Missing | New |
| S2 case DNA | ❌ Missing | New |
| S3 unknown protein | ❌ Missing | New |
| S4 unknown nucleotide | ❌ Missing | New |
| C1 bond invariant protein | ❌ Missing | New |
| C2 bond invariant DNA | ❌ Missing | New |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateMolecularWeight_Tests.cs` — all SEQ-MW-001 cases.
- **Remove:** none (no pre-existing MW-specific tests found in the legacy fixture).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateMolecularWeight_Tests.cs` | Canonical SEQ-MW-001 fixture | 16 |

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
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | S4 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Implemented, exact value |
| M2 | ✅ | Implemented, exact value |
| M3 | ✅ | Implemented, exact value |
| M4 | ✅ | Implemented, exact value |
| M5 | ✅ | Implemented, exact value |
| M6 | ✅ | Implemented, exact value |
| M7 | ✅ | Implemented |
| M8 | ✅ | Implemented |
| M9 | ✅ | Implemented |
| M10 | ✅ | Implemented |
| S1 | ✅ | Implemented |
| S2 | ✅ | Implemented |
| S3 | ✅ | Implemented |
| S4 | ✅ | Implemented |
| C1 | ✅ | Implemented |
| C2 | ✅ | Implemented |

Total in-scope cases: 16. ✅ count: 16.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| ASSUMPTION-01 | Input is upper-cased before lookup (case folding; does not change cited values) | S1, S2 |
| ASSUMPTION-02 | Unknown symbols are skipped (no mass, no bond), deviating from Biopython's reject behavior | S3, S4 |

---

## 7. Open Questions / Decisions

1. Decision: nucleotide MW corrected to the Biopython average-mass tables and `(len−1)·water` phosphodiester correction (previous code summed raw monophosphate masses with no water removal — a defect vs. the reference). Confirmed by re-deriving the published AGC examples.
2. Decision: protein MW already followed the Expasy definition; constants upgraded to the 4-decimal Biopython/Expasy values for exactness.
