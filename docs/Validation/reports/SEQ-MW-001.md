# Validation Report: SEQ-MW-001 — Molecular Weight Calculation (protein & nucleotide)

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateMolecularWeight(string)` (average-mass protein Mw);
  `SequenceStatistics.CalculateNucleotideMolecularWeight(string, bool isDna = true)` (average-mass DNA/RNA Mw)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent retrieval)

1. **Biopython `Bio/Data/IUPACData.py`** (raw master) — fetched the verbatim weight tables.
   - `protein_weights` (average): A 89.0932, C 121.1582, D 133.1027, E 147.1293, F 165.1891,
     G 75.0666, H 155.1546, I 131.1729, K 146.1876, L 131.1729, M 149.2113, N 132.1179,
     O 255.3134, P 115.1305, Q 146.1445, R 174.201, S 105.0926, T 119.1192, U 168.0532,
     V 117.1463, W 204.2252, Y 181.1885.
   - `unambiguous_dna_weights`: A 331.2218, C 307.1971, G 347.2212, T 322.2085.
   - `unambiguous_rna_weights`: A 347.2212, C 323.1965, G 363.2206, U 324.1813.
   - These match the 20-AA table and DNA/RNA tables in the implementation **exactly**.
2. **Biopython `Bio/SeqUtils/__init__.py`** (`molecular_weight`) — confirmed:
   - water = `18.0153` (average), `18.010565` (monoisotopic).
   - single-strand formula `weight = sum(weight_table[x] for x in seq) − (len(seq) − 1) * water`
     (same formula for protein, DNA, RNA).
   - docstring worked examples for "AGC": DNA 949.61, RNA 997.61, protein 249.29.
   - unknown letters → `ValueError` (reject); `double_stranded` for protein → `ValueError`;
     `circular` subtracts one extra water.
3. **Expasy FindMod — Average masses of amino acid residues** (web.expasy.org). Confirmed:
   - H₂O average mass = **18.01524** (Biopython rounds to 18.0153).
   - Residue masses (in-chain): Ala 71.0788, Gly 57.0519, Cys 103.1388, Phe 147.1766; these are
     **residue** masses → free-amino-acid mass = residue + water (Ala 71.0788+18.0153 ≈ 89.094,
     matching Biopython 89.0932 to PubChem rounding).
4. **Expasy Compute pI/Mw doc** (cited in Evidence) — "Protein Mw = sum of average isotopic masses
   of amino acids + average isotopic mass of one water molecule." Algebraically identical to
   `Σ free-aa − (n−1)·water` (since `Σresidue + water = Σ(free-aa − water) + water = Σfree-aa − (n−1)·water`).

### Formula check

The implemented formula `Σ table[x] − (n−1)·W` with `W = 18.0153` is the exact Biopython
single-strand formula and the algebraic equivalent of the Expasy protein definition. PASS.

### Edge-case semantics

- Single monomer ⇒ (n−1)=0 ⇒ free monomer mass. Sourced (Expasy/Biopython). ✓
- null/empty ⇒ 0 — an API contract choice (sources define n≥1); documented. ✓
- Unknown symbols **skipped** (no mass, no bond) — a *deviation* from Biopython's reject-on-unknown,
  explicitly documented as ASSUMPTION-02. Defensible (no invented mass), recorded. ✓
- Average-only mass set; monoisotopic / double-stranded / circular out of scope — documented. ✓

### Independent cross-check (numbers)

Hand-computed from the externally-fetched Biopython tables (water 18.0153):

| Input | type | Recomputed | Biopython docstring |
|-------|------|-----------|---------------------|
| AGC | protein | 89.0932+75.0666+121.1582−2·18.0153 = **249.2874** | 249.29 ✓ |
| AGC | DNA | 331.2218+347.2212+307.1971−2·18.0153 = **949.6095** | 949.61 ✓ |
| AGC | RNA | 347.2212+363.2206+323.1965−2·18.0153 = **997.6077** | 997.61 ✓ |

Second, fully independent peptide cross-check (not from the docstring):

- **Met-enkephalin** (Tyr-Gly-Gly-Phe-Met, YGGFM): formula gives **573.66 Da**; Wikipedia
  (Met-enkephalin) lists average MW **573.66 Da**. ✓
- **Insulin A-chain linear** (GIVEQCCTSICSLYQLENYCN): formula gives **2383.7 Da**, consistent with
  the literature reduced-A-chain average mass (~2383.8). ✓

### Findings / divergences (Stage A)

- **NOTE-A (typo, FIXED).** The Evidence doc (line 104) and TestSpec (§1.2 pt 4 and M3 row) stated
  the RNA "AGC" value as **997.6177**. The correct value from the cited Biopython RNA table is
  **997.6077** (347.2212+363.2206+323.1965−2·18.0153). The off-by-0.01 came from a transcription
  slip in the docs; the **code and the M3 test already use the correct 997.6077**. Corrected the two
  docs this session. No code/biology error — hence Stage A = PASS-WITH-NOTES, not FAIL.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`:
- constants L149 (`AverageWaterMass = 18.0153`), L154–161 (AA table), L165–168 (DNA table),
  L172–175 (RNA table) — all byte-for-byte equal to the externally-fetched Biopython tables.
- `CalculateMolecularWeight` L188–210 and `CalculateNucleotideMolecularWeight` L223–247:
  null/empty→0; `ToUpperInvariant`; accumulate only table-recognized monomers + count; if count==0
  return 0; return `weight − (count−1)·W`.

### Formula realised correctly?

Yes — the code computes `Σ recognized-monomer-mass − (recognized_count − 1)·18.0153`, exactly the
validated formula. The `(count−1)` uses the *recognized* count, which is the correct behaviour given
the skip-unknown assumption (an unknown symbol forms no bond).

### Cross-verification table recomputed vs code (test run)

All 18 tests in the canonical fixture pass (`Failed: 0`). Exact values asserted:
protein AGC 249.2874; DNA AGC 949.6095; RNA AGC 997.6077; single G 75.0666; single dA 331.2218;
single rA 347.2212; AG 146.1445; DNA AG 660.4277; DNA ACGT 1253.8027; RNA ACGU 1303.7737 — each
traced to the Biopython tables above.

### Variant/delegate consistency

`AnalyzeSequence` (L116/L134) calls `CalculateMolecularWeight` for its `MolecularWeight` field —
same canonical method, no divergent reimplementation. No `*Fast` variant exists.

### Test quality audit (HARD gate)

Pre-existing fixture had 16 tests, all asserting **exact sourced values** with a tight `1e-4`
tolerance (not Greater/AtLeast/ranges) — good. Gaps found and fixed this session:

- The **"all-unknown ⇒ 0"** branch (code L205–206 / L242–243), a distinct path from the null/empty
  short-circuit and a documented edge case (spec §3.3), was untested. Added
  `CalculateMolecularWeight_AllUnknownSymbols_ReturnsZero` ("***"→0) and the nucleotide twin.
- The RNA **`U`** and DNA **`T`** table entries were never exercised (AGC omits both). Added exact-value
  tests `…_DnaACGT_…` (1253.8027) and `…_RnaACGU_…` (1303.7737), both traced to the Biopython tables.

After fixes the fixture has **18 tests**. No assertion was weakened, no tolerance widened, no test
skipped, no expected value bent to match code. Every expected number traces to an external source
retrieved this session (Biopython tables / Expasy / Wikipedia), not to the implementation's output.

**Gate result: PASS** — exact sourced expectations; all public methods + overloads and all Stage-A
branches (null, empty, single, dipeptide bond-count, lowercase, unknown-skip, all-unknown→0, full
alphabets, DNA & RNA) covered; full unfiltered suite green; changed file builds warning-free.

### Findings / defects (Stage B)

- No algorithm defect. Two **test-coverage gaps** (all-unknown branch; U/T table entries) — fixed.
- The 4 build warnings (NUnit2007) are pre-existing in an unrelated file
  (`ApproximateMatcher_EditDistance_Tests.cs`); the SEQ-MW-001 file builds 0-warning.

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** (one doc-only typo 997.6177→997.6077, fixed; biology/formula correct).
- **Stage B: PASS** (code faithful to the validated formula; tests strengthened to lock sourced values).
- **End-state: ✅ CLEAN** — no algorithm defect; the doc typo and the test-coverage gaps were fully
  fixed this session. `dotnet build` 0 errors (SEQ-MW-001 file 0 warnings); full unfiltered
  `dotnet test` = **6516 passed, 0 failed** (1 pre-existing skipped benchmark).

### Build/test evidence

```
dotnet build … Seqeron.Genomics.Tests.csproj -c Debug  → 0 Error(s)
dotnet test  … (full, unfiltered, --no-build)          → Failed: 0, Passed: 6516, Skipped: 1
```
