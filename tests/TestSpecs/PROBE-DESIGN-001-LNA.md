# Test Specification: PROBE-DESIGN-001 (LNA-adjusted NN Tm)

**Test Unit ID:** PROBE-DESIGN-001
**Area:** MolTools
**Algorithm:** LNA (locked nucleic acid)-adjusted nearest-neighbour melting temperature; citable MGB design rules
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-24

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | McTigue PM, Peterson RJ, Kahn JD (2004). Sequence-dependent thermodynamic parameters for LNA-DNA duplex formation. *Biochemistry* 43:5388–5405. | 1 | https://doi.org/10.1021/bi035976d | 2026-06-24 |
| 2 | MELTING 5 `McTigue2004lockedmn.xml` + `McTigue04LockedAcid.java` (reference impl., verbatim McTigue table) | 3 | github mohakjain/TmCalculator/MELTING5.2.0 ; rmelting inst/extdata | 2026-06-24 |
| 3 | rmelting tutorial worked example (`mct04`) | 3 | https://aravind-j.github.io/rmelting/articles/Tutorial.html | 2026-06-24 |
| 4 | Kutyavin IV et al. (2000). 3'-MGB-DNA probes increase sequence specificity. *Nucleic Acids Res* 28(2):655–661. | 1 | https://academic.oup.com/nar/article/28/2/655/1039630 | 2026-06-24 |
| 5 | SantaLucia J (1998). *PNAS* 95(4):1460–65 (base DNA NN model) | 1 | https://www.pnas.org/doi/10.1073/pnas.95.4.1460 | 2026-06-24 |

### 1.2 Key Evidence Points

1. McTigue 2004 reports ΔΔH°/ΔΔS° increments for all 32 LNA+DNA:DNA nearest neighbours; an internal LNA substitution **raises** the duplex Tm (substantial stabilization, sequence-dependent) — source 1, 2.
2. The LNA value is an **additive increment** to the underlying DNA NN ΔH°/ΔS° for the step containing the LNA base (MELTING `computeThermodynamics`: DNA NN sum, then `enthalpy += lockedAcidValue`) — source 2.
3. XML enthalpy/entropy are in **cal/mol** and **cal/(mol·K)** (MELTING) — source 2; so `992.0` → +0.992 kcal/mol.
4. **Terminal** LNA (duplex position 0 or last) is **not parameterised** → reject — source 2 (`isApplicable`).
5. MELTING worked example: `CCATTLGCTACC` (DNA `CCATTGCTACC`, LNA index 4), C=1e-4, Na=1 → Tm = 63.61426 °C — source 3.
6. **MGB:** attached at the **3' end**; enables **12–20mer** (shorter) probes; a 12mer MGB ≈ Tm of a 27mer unmodified; the quantitative MGB ΔTm is **empirical with no published formula** — source 4. ⇒ MGB *design rules* citable; quantitative MGB ΔTm is an honest residual.

### 1.3 Documented Corner Cases

- Terminal LNA position has no McTigue parameter (reject).
- Non-ACGT base → underlying DNA NN lookup fails → not-computable.
- LNA index out of range / empty / < 2 nt → not-computable.

### 1.4 Known Failure Modes / Pitfalls

1. Treating a terminal LNA as internal (wrong, no parameter) — source 2.
2. Confusing the two LNA orientations of a step (5'-locked `XLY` vs 3'-locked `XYL`) — both must map to the correct verbatim key — source 2.
3. Using a single average ΔΔ instead of the actual sequence-dependent value — source 1.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateNearestNeighborThermodynamicsLna(string sequence, IReadOnlyCollection<int> lnaPositions)` | `PrimerDesigner` | Canonical | LNA-adjusted ΔH°/ΔS° (base DNA NN + McTigue increments). |
| `CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection<int>, …)` | `PrimerDesigner` | Canonical | LNA-adjusted NN Tm (reuses the bimolecular Tm equation + salt corrections). |
| `EvaluateMgbProbeDesign(string sequence)` | `ProbeDesigner` | Canonical | Citable MGB design-rule check (3'-MGB; 12–20mer length window). |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Adding an internal LNA monomer raises the NN Tm vs the all-DNA duplex (stabilization). | Yes | McTigue 2004 (source 1) |
| INV-2 | LNA increment is additive to the base DNA NN: with no LNA positions, the LNA Tm equals the plain `CalculateMeltingTemperatureNN`. | Yes | MELTING (source 2) |
| INV-3 | A terminal LNA (index 0 or length−1) is not computable (null/NaN). | Yes | MELTING `isApplicable` (source 2) |
| INV-4 | The applied increment equals the verbatim McTigue value (cal/mol → kcal/mol) for the correct (NN step, locked-position) key. | Yes | source 2 (XML) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | LNA ΔH°/ΔS° of `CCATTGCTACC`, LNA index 4 | base DNA NN + `TTL/AA`(+2326,+8.1) + `TLG/AC`(−1540,−3.0) | ΔH° = −80.014 kcal/mol; ΔS° = −216.6 cal/(mol·K) | Evidence worked example (sources 2,5) |
| M2 | LNA Tm of that duplex (C=1e-4, Na=1, salt mode None) | bimolecular Tm equation, x=4 | Tm = 63.52759 °C (Within 1e-4); within 0.1 °C of MELTING 63.61426 | source 3 |
| M3 | LNA raises Tm vs all-DNA | same duplex, with vs without LNA index 4 | Tm(LNA) = 63.52759 > Tm(DNA) = 59.69230; ΔTm ≈ +3.84 | source 1 |
| M4 | No LNA positions ⇒ equals plain NN Tm | `CalculateMeltingTemperatureNNLna(seq, {})` vs `CalculateMeltingTemperatureNN(seq)` | exactly equal | source 2 (INV-2) |
| M5 | Terminal LNA rejected | LNA at index 0 and at last index | thermodynamics null; Tm NaN | source 2 (INV-3) |
| M6 | Negative-ΔΔ key applied with correct sign | `GGGCC...`-context step `GLG/CC` = (−2844, −6.7) | ΔH°/ΔS° lowered by exactly that increment vs base DNA | source 2 |
| M7 | MGB design rules | `EvaluateMgbProbeDesign` on a 15-mer and a 25-mer | 15-mer length-OK; 25-mer flagged "outside 12–20"; both report 3'-MGB placement guidance | source 4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | null / empty / 1-nt sequence | LNA thermo | null / NaN | guard |
| S2 | LNA index out of range / duplicate / unsorted | thermo with {7} on len-5, {2,2}, {4,1} | out-of-range → null; duplicates & order tolerated | guard |
| S3 | non-ACGT base | `CCANTGCTACC` LNA idx 4 | null / NaN | base NN fails |
| S4 | all 32 increment keys present | reflective/explicit count | 32 entries | completeness |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | two internal LNA positions | sum of both increments | both applied additively | linearity |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Existing probe tests: `ProbeDesigner_ProbeDesign_Tests.cs`, `ProbeDesigner_TaqMan_Tests.cs`, `ProbeDesigner_ProbeValidation_Tests.cs`, `ProbeDesigner_MutationKillers_Tests.cs`, `ProbeDesignerTests.cs`. NN Tm tests: `PrimerDesigner_NearestNeighborTm_Tests.cs`. None cover LNA-adjusted Tm or MGB rules.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M7, S1–S4, C1 | ❌ Missing | new LNA/MGB feature; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProbeDesigner_LnaTm_Tests.cs` — all LNA-Tm + MGB-rule tests.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProbeDesigner_LnaTm_Tests.cs` | LNA-adjusted NN Tm + MGB design rules | 12 |

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
| M1 | ✅ | exact ΔH°/ΔS° asserted |
| M2 | ✅ | Tm 63.52759 °C asserted; MELTING within 0.1 °C |
| M3 | ✅ | LNA > DNA Tm asserted |
| M4 | ✅ | equals plain NN Tm |
| M5 | ✅ | terminal LNA → null/NaN |
| M6 | ✅ | negative increment sign |
| M7 | ✅ | MGB length window + 3' placement |
| S1 | ✅ | null/empty/short guard |
| S2 | ✅ | index range / order |
| S3 | ✅ | non-ACGT |
| S4 | ✅ | 32-key completeness |
| C1 | ✅ | two LNA positions |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Base DNA NN = SantaLucia 1998 unified (library's existing model); the ~0.09 °C offset vs MELTING `mct04` is from this base-model choice, not the increments. | M2 (tolerance vs MELTING) |

---

## 7. Open Questions / Decisions

1. Quantitative MGB ΔTm has no published formula (Kutyavin 2000 describes it as empirical) → left as an honest residual; only the citable MGB *design rules* (3'-MGB, 12–20mer) are implemented.
