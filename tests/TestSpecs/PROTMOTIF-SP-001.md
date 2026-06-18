# Test Specification: PROTMOTIF-SP-001

**Test Unit ID:** PROTMOTIF-SP-001
**Area:** ProteinMotif
**Algorithm:** Signal Peptide Cleavage-Site Prediction (von Heijne 1986 weight matrix)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | von Heijne G. (1986) A new method for predicting signal sequence cleavage sites. Nucleic Acids Res. 14(11):4683–4690 | 1 | https://doi.org/10.1093/nar/14.11.4683 | 2026-06-14 |
| 2 | EMBOSS 6.6.0 `sigcleave` application doc (embeds eukaryotic matrix, formula, worked example) | 3 | https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/sigcleave.html | 2026-06-14 |
| 3 | EMBOSS 6.6.0 `emboss/sigcleave.c` source (scoring loop + matrix transform) | 3 | https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/sigcleave.c | 2026-06-14 |
| 4 | EMBOSS 6.6.0 `data/Esig.euk`, `data/Esig.pro` (count matrices, verbatim) | 3 | https://raw.githubusercontent.com/lauringlab/CodonShuffle/master/lib/EMBOSS-6.6.0/emboss/data/Esig.euk | 2026-06-14 |
| 5 | UniProt P17644 (ACH2_DROME) FASTA | 5 | https://rest.uniprot.org/uniprotkb/P17644.fasta | 2026-06-14 |

### 1.2 Key Evidence Points

1. Score at a candidate site = sum of log-odds weights `ln(count/expect)` over positions −13..+2 (15 columns), natural log — Source 3 (`sigcleave_readSig`, scoring loop).
2. Zero counts are replaced by `1.0e-10` at columns −3 and −1 (a strong penalty) and by `1.0` elsewhere, before the log — Source 3.
3. The single prediction is the argmax of the weight over all positions; mature protein starts at that index — Source 3 (`maxweight`/`maxsite`).
4. Cleavage is between positions −1 and +1; +1 is the first residue of the mature protein — Source 2, 4.
5. Acceptance threshold `-minweight` default 3.5 (≈95% sensitivity/specificity) — Source 2.
6. Worked example ACH2_DROME: maximum score 13.739, mature-protein start residue 42 (`LLVLLLLCETVQA` = positions −13..−1) — Source 2; reproduced exactly with the implemented matrix.

### 1.3 Documented Corner Cases

- No intrinsic cutoff: a best site is always reported for in-window sequences; `-minweight` only flags likelihood (Source 2).
- Accuracy ≈75–80% for the cleavage site — heuristic, not exact (Source 2).
- Eukaryotic matrix is the default; prokaryotic matrix is selected by an option (Source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Using log base 10 instead of natural log changes every score — Source 3 uses C `log()` (natural). Verified by reproducing 13.739.
2. Forgetting the −3/−1 zero-count penalty (1.0e-10) inflates scores at sites with absent conserved residues — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `PredictSignalPeptide(string, bool, double)` | ProteinMotifFinder | Canonical | von Heijne weight-matrix score + argmax site |
| `BuildWeightMatrix(int[][], double[])` | ProteinMotifFinder | Internal | log-odds transform; tested via the public method's exact scores |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Returned `Score` equals the maximum over all candidate sites of `Σ ln(count/expect)` for positions −13..+2 | Yes | Source 3 |
| INV-2 | `CleavagePosition` (1-based mature start) ∈ [1, len]; cleavage is between `CleavagePosition−1` and `CleavagePosition` | Yes | Source 2,3 |
| INV-3 | `IsLikelySignalPeptide` ⇔ `Score ≥ minWeight` (default 3.5) | Yes | Source 2 |
| INV-4 | Result is independent of input letter case | Yes | Implementation upper-cases input |
| INV-5 | Inputs shorter than one full 15-residue window return `null` | Yes | Window requirement (ASSUMPTION-1) |
| INV-6 | `SignalSequence` = residues 1..`CleavagePosition−1`; `WindowSequence` length ≤ 15 | Yes | Source 3 window definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | ACH2_DROME cleavage position | Full UniProt P17644 sequence | `CleavagePosition == 42` | Source 2 worked example |
| M2 | ACH2_DROME score | Same input | `Score == 13.7390400704164` (±1e-3) | Source 2 + re-derivation |
| M3 | ACH2_DROME signal sequence | Same input | `SignalSequence` = first 41 residues; `WindowSequence == "LLVLLLLCETVQA"` ends at −1=A | Source 2 |
| M4 | Argmax over alternatives | Score of the runner-up site (mature start 39) is strictly lower | site-39 score `12.135…` < 13.739 | Source 2 ordering |
| M5 | IsLikely true | ACH2_DROME | `IsLikelySignalPeptide == true` (13.739 ≥ 3.5) | Source 2 threshold |
| M6 | Null input | `null` sequence | returns `null` | window/precondition |
| M7 | Empty input | `""` | returns `null` | window/precondition |
| M8 | Short input | 14-residue sequence | returns `null` | INV-5 |
| M9 | Minimum-length input | exactly 15 residues | returns non-null (one window scored) | INV-5 |
| M10 | Log-odds formula | Hand-built 15-residue window; expected = Σ of specific matrix cells | `Score` equals the hand-summed log-odds (±1e-9) | Source 3,4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity | lower vs upper ACH2_DROME | identical `CleavagePosition` and `Score` | INV-4 |
| S2 | minWeight raised above score | ACH2_DROME with `minWeight = 14.0` | `IsLikelySignalPeptide == false`, `Score` unchanged | INV-3 |
| S3 | Prokaryotic matrix selected | `prokaryote: true` on a sequence | scores against `Esig.pro`; result may differ from eukaryotic | Source 2,4 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Non-standard residue | Window containing `X` | no exception; `X` contributes 0 | Source 3 residue handling |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `ProteinMotifFinder_DomainPrediction_Tests.cs` previously contained signal-peptide tests (M7–M15, S5, S6, C2) built on the prior **fabricated** tripartite scoring (constants 0.95, 0.825; NRegion="MKRLL"; score range [0,1]). These encode invented values and were removed.
- `Properties/ProteinMotifProperties.cs` had `PredictSignalPeptide_CleavagePosition_WithinBounds` asserting `Score ∈ [0,1]` (no longer valid) — removed (replaced by C-bound coverage in the new file).
- `Snapshots/ProteinMotifSnapshotTests.cs` had a `PredictSignalPeptide_Snapshot` referencing the removed `Probability` field — updated to the new record and re-verified.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M10 | ❌ Missing | new canonical file |
| S1–S3 | ❌ Missing | new canonical file |
| C1 | ❌ Missing | new canonical file |
| (old) DomainPrediction signal tests | 🔁 Duplicate / ⚠ Weak | fabricated values — removed |
| (old) Properties score-range test | ⚠ Weak | invalid [0,1] bound — removed |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_PredictSignalPeptide_Tests.cs` — all SP-001 cases.
- **Remove:** signal-peptide tests from `ProteinMotifFinder_DomainPrediction_Tests.cs` (DOMAIN-001 keeps only domain tests); the score-range property test from `ProteinMotifProperties.cs`. Update the snapshot test/verified file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_PredictSignalPeptide_Tests.cs` | Canonical SP-001 | 14 |
| `ProteinMotifFinder_DomainPrediction_Tests.cs` | DOMAIN-001 (domains only) | unchanged domain tests |
| `ProteinMotifSnapshotTests.cs` | Snapshot (updated fields) | 1 (signal-peptide snapshot) |

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
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | old signal tests | 🔁/⚠ | Removed from DOMAIN-001 / Properties / snapshot updated | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact value |
| M2 | ✅ | Exact value (±1e-3) |
| M3 | ✅ | Exact strings |
| M4 | ✅ | Strict inequality vs known runner-up |
| M5 | ✅ | Boolean |
| M6 | ✅ | null |
| M7 | ✅ | null |
| M8 | ✅ | null |
| M9 | ✅ | non-null |
| M10 | ✅ | Hand-summed log-odds |
| S1 | ✅ | Case invariance |
| S2 | ✅ | Threshold boolean |
| S3 | ✅ | Prokaryotic distinct |
| C1 | ✅ | No-throw, X contributes 0 |

Total in-scope cases: 14. ✅ count: 14.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Minimum input length = one full 15-residue window; shorter inputs return `null` | INV-5, M8, M9 |

---

## 7. Open Questions / Decisions

1. **Decision:** The record was redesigned to the honest weight-matrix model (`CleavagePosition`, `Score`, `SignalSequence`, `WindowSequence`, `IsLikelySignalPeptide`), replacing the prior fabricated n/h/c tripartite fields. The MCP `SignalPeptideResult` and tool signature were updated accordingly (in-scope conformance correction).
2. **Decision:** Prior block resolved — the von Heijne (1986) matrix and the deterministic argmax selection model were retrieved in text from the EMBOSS reference implementation (data files + source) and the worked example reproduced exactly (13.739).
