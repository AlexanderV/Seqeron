# Test Specification: PRIMER-TM-001-SPECIAL-LOOP

**Test Unit ID:** PRIMER-TM-001 (special tri/tetraloop hairpin bonus tables — bundled & applied)
**Area:** MolTools
**Algorithm:** Full Primer3 `ntthal` intramolecular-hairpin DP with bundled sequence-specific
special triloop / tetraloop stability bonuses
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | primer3 `triloop.dh`/`.ds` | 3 | https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_config/triloop.dh | 2026-06-25 |
| 2 | primer3 `tetraloop.dh`/`.ds` | 3 | …/primer3_config/tetraloop.dh | 2026-06-25 |
| 3 | primer3 `thal.c` `calc_hairpin` | 3 | …/libprimer3/thal.c | 2026-06-25 |
| 4 | SantaLucia & Hicks (2004) special hairpin loops | 1 | https://doi.org/10.1146/annurev.biophys.32.110601.141800 | 2026-06-25 |
| 5 | primer3-py 2.3.0 `calc_hairpin` | 3 | https://pypi.org/project/primer3-py/2.3.0/ | 2026-06-25 |

### 1.2 Key Evidence Points

1. The special triloop (16 rows) and tetraloop (76 rows) bonus tables are keyed on the **full loop
   string including the closing base pair** (5 chars triloop, 6 chars tetraloop) — thal.c
   `readTLoop` + bsearch key `numSeq1 + i`.
2. The bonus ΔH (cal/mol) and ΔS (cal/(K·mol)) are **added** to the loop ΔH°/ΔS° for `loopSize == 3`
   (triloop) / `loopSize == 4` (tetraloop) — thal.c `calc_hairpin` lines 2106-2127.
3. primer3-py `calc_hairpin` ground-truth ΔH/ΔS/ΔG/Tm captured this session (see Evidence).

### 1.3 Documented Corner Cases

- Only 3-nt and 4-nt loops have special tables; 5-nt+ loops never get a bonus.
- The key includes the closing pair — a matching 3/4-nt loop body with a different closing pair gets no bonus.
- Homopolymer / too-short input → ntthal `no_structure`.

### 1.4 Known Failure Modes / Pitfalls

1. Encoding the loop key by 3-bit packing collides on a 4-bit mask — must pack/extract consistently
   (the implementation uses 4-bit nibbles) — pitfall found and fixed this session.
2. The salt base-pair count N is counted over `bp[0..len−2]` (excludes the last index) — off-by-one
   shifts ΔS by one `saltCorrection` unit — thal.c calcHairpin.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateHairpinThermodynamicsNtthal(string, double)` | `PrimerDesigner` | **Canonical** | Full ntthal hairpin DP + bundled special-loop bonuses |
| `NtthalHairpin.Run(string, double)` | `NtthalHairpin` | **Internal** | Tested indirectly via the canonical method |
| `FindMostStableHairpin` / `CalculateHairpinMeltingTemperature` | `PrimerDesigner` | **Delegate (regression)** | Legacy Table-4 model — must be UNCHANGED |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A recognised tri/tetraloop hairpin's ΔH/ΔS/ΔG/Tm equals primer3 `calc_hairpin`. | Yes | Sources 1-5 |
| INV-2 | A non-special-loop hairpin's result is identical with or without the bundled tables. | Yes | Source 3 (bsearch no-match) + Source 5 |
| INV-3 | No hairpin (homopolymer / too short) → `null`. | Yes | Source 5 (structure_found=False) |
| INV-4 | The legacy Table-4 `FindMostStableHairpin` + `loopBonusDeltaG37` path is unchanged. | Yes | PRIMER-TM-001-HAIRPIN evidence |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | tetraloop CGAAAG | GGGGCGAAAGCCCC | dH=−40900, dS=−114.1872884299936, dG=−5484.812493437487, Tm=85.03347700825856 | Sources 1-5 |
| M2 | tetraloop GGGGAC | GGGGGGGACCCCC | dH=−34000, dS=−94.1872884299836, Tm=87.8328944728006 | Sources 2,5 |
| M3 | triloop CGAAG | GGGCGAAGCCC | dH=−27800, dS=−77.68485895331574, dG=−3706.040995629125, Tm=84.7060915802943 | Sources 1,3,5 |
| M4 | triloop GGAAC | GGGGGAACCCC | dH=−26000, dS=−73.18485895331571, Tm=82.11474153055735 | Sources 1,5 |
| S1 | non-special 4-nt | GGGCTTTTGCCC | dH=−32400, dS=−94.58485895332572, Tm=69.39954078842845 (unchanged) | Sources 3,5 |
| S2 | non-special 5-nt | GGGCAAAAAGCCC | dH=−30100, dS=−87.74485895332572, Tm=69.89004085311882 (unchanged) | Source 5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S3 | homopolymer | AAAAAAAAAAAA | null | ntthal no_structure |
| S4 | too short | GCGC | null | cannot close a ≥3-nt loop |
| S5 | invalid input | null / "" / "GGGCNAAGCCC" | null | guards |
| C1 | legacy unchanged | FindMostStableHairpin(GGGCTTTTGCCC) | Table-4 ΔH/ΔS/ΔG unchanged | regression |
| C2 | manual bonus path | FindMostStableHairpin(…, loopBonusDeltaG37=1.0) | ΔG°37 shifts by exactly +1.0 | regression |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| (none) | | | | |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Hairpin tests existed in `PrimerDesigner_HairpinTm_Tests.cs` (legacy Table-4 model). No tests
  existed for bundled special-loop bonuses or the ntthal hairpin DP. New canonical file:
  `tests/Seqeron/Seqeron.Genomics.Tests/PrimerDesigner_HairpinSpecialLoop_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 tetraloop CGAAAG | ❌ Missing | new method |
| M2 tetraloop GGGGAC | ❌ Missing | new method |
| M3 triloop CGAAG | ❌ Missing | new method |
| M4 triloop GGAAC | ❌ Missing | new method |
| S1 non-special 4-nt regression | ❌ Missing | new |
| S2 non-special 5-nt regression | ❌ Missing | new |
| S3 homopolymer | ❌ Missing | new |
| S4 too short | ❌ Missing | new |
| S5 invalid input | ❌ Missing | new |
| C1 legacy unchanged | ❌ Missing | new |
| C2 manual bonus path | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `PrimerDesigner_HairpinSpecialLoop_Tests.cs` — all M/S/C cases above.
- **Remove:** nothing (the legacy hairpin tests remain valid for the Table-4 model).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `PrimerDesigner_HairpinSpecialLoop_Tests.cs` | special-loop + ntthal parity + regression | 11 |
| `PrimerDesigner_HairpinTm_Tests.cs` | legacy Table-4 hairpin (unchanged) | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | Implemented | ✅ Done |
| 2 | M2 | ❌ | Implemented | ✅ Done |
| 3 | M3 | ❌ | Implemented | ✅ Done |
| 4 | M4 | ❌ | Implemented | ✅ Done |
| 5 | S1 | ❌ | Implemented | ✅ Done |
| 6 | S2 | ❌ | Implemented | ✅ Done |
| 7 | S3 | ❌ | Implemented | ✅ Done |
| 8 | S4 | ❌ | Implemented | ✅ Done |
| 9 | S5 | ❌ | Implemented | ✅ Done |
| 10 | C1 | ❌ | Implemented | ✅ Done |
| 11 | C2 | ❌ | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | parity to machine precision |
| M2 | ✅ | parity |
| M3 | ✅ | parity |
| M4 | ✅ | parity |
| S1 | ✅ | regression held |
| S2 | ✅ | regression held |
| S3 | ✅ | null |
| S4 | ✅ | null |
| S5 | ✅ | null guards |
| C1 | ✅ | legacy unchanged |
| C2 | ✅ | manual bonus shifts ΔG by +1.0 |

---

## 6. Assumption Register

**Total assumptions:** 0

| # | Assumption | Used In |
|---|-----------|---------|
| (none) | All tables/keying/convention taken verbatim from libprimer3 + cross-verified vs primer3-py. | — |

---

## 7. Open Questions / Decisions

1. None. The bundled special-loop bonus tables are applied via the full ntthal hairpin DP and
   reproduce primer3-py `calc_hairpin` to machine precision; the legacy Table-4 model is unchanged.
