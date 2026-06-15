# Test Specification: MIRNA-PAIR-001

**Test Unit ID:** MIRNA-PAIR-001
**Area:** MiRNA
**Algorithm:** MiRNA-Target Pairing Analysis (miRNA–mRNA duplex base pairing)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Agarwal et al. (2015) eLife — Predicting effective miRNA target sites | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/ | 2026-06-13 |
| 2 | Bartel (2009) Cell — MicroRNAs: Target Recognition | 1 | https://www.cell.com/fulltext/S0092-8674(09)00008-7 | 2026-06-13 |
| 3 | Transcriptome-Wide Non-Canonical Interactions (PMC4870184) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4870184/ | 2026-06-13 |
| 4 | Crick (1966) wobble — via Wikipedia "Wobble base pair" | 4 | https://en.wikipedia.org/wiki/Wobble_base_pair | 2026-06-13 |
| 5 | Lewis et al. (2005) Cell — Conserved seed pairing | 1 | https://pubmed.ncbi.nlm.nih.gov/15652477/ | 2026-06-13 |
| 6 | NNDB Turner 2004 parameters | 1/3 | https://rna.urmc.rochester.edu/NNDB/turner04/index.html | 2026-06-13 |

### 1.2 Key Evidence Points

1. Canonical pairing is Watson-Crick: A pairs with U, C pairs with G — Agarwal et al. (2015) "‘A' pairs with ‘U', and ‘C' pairs with ‘G'".
2. The only standard RNA wobble pair is G-U; it is distinguished from Watson-Crick — Crick (1966) via Wikipedia; PMC4870184 "solid lines indicate Watson-Crick base pairing" (vs G:U wobble).
3. The miRNA seed is nucleotides 2–8 and pairs antiparallel/Watson-Crick complementary to the mRNA 3′UTR — Bartel (2009); PMC4870184 "(position 2-8)".
4. Target sites are the reverse complement of the miRNA seed — Lewis et al. (2005) "conserved complementarity to the seed (nucleotides 2-7)".
5. RNA duplex free energy = sum of nearest-neighbor stacking energies (all negative at 37 °C for paired stacks) + initiation/terminal penalties — NNDB Turner 2004 / Xia et al. (1998).

### 1.3 Documented Corner Cases

- G:U is a wobble, NOT a Watson-Crick match — must be counted separately (PMC4870184).
- A opposite miRNA position 1 is Argonaute recognition, not a base pair (Agarwal et al. 2015) — not a "match" in a pure aligner.
- A-C, A-A, etc. do not pair → mismatches (Agarwal et al.: A pairs only with U).

### 1.4 Known Failure Modes / Pitfalls

1. Counting G:U as a Watson-Crick match would over-count complementarity — Crick (1966); PMC4870184.
2. Aligning miRNA to target in the same (parallel) orientation rather than antiparallel reverse-complement — Lewis et al. (2005).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AlignMiRnaToTarget(miRna, target)` | MiRnaAnalyzer | Canonical | Antiparallel duplex base-pairing alignment; counts matches/mismatches/wobbles, builds alignment string, estimates ΔG. |
| `GetReverseComplement(sequence)` | MiRnaAnalyzer | Canonical | RNA reverse complement (T→U), used to derive seed-match targets. |
| `CanPair(base1, base2)` | MiRnaAnalyzer | Canonical | Watson-Crick + G:U wobble pairing predicate. |
| `IsWobblePair(base1, base2)` | MiRnaAnalyzer | Canonical | G:U wobble predicate. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `CanPair` true ⇔ {A-U, U-A, G-C, C-G, G-U, U-G}; false otherwise | Yes | Agarwal et al. (A-U,C-G); Crick (1966) (G-U) |
| INV-2 | `IsWobblePair` true ⇔ {G-U, U-G}; Watson-Crick pairs are NOT wobble | Yes | PMC4870184; Crick (1966) |
| INV-3 | Every G:U/U:G counted by `IsWobblePair` is also pairable by `CanPair` (wobbles ⊆ pairs) | Yes | Crick (1966); CanPair definition |
| INV-4 | `GetReverseComplement` reverses and complements (A↔U, G↔C, T→A); length preserved; idempotent under double application for RNA | Yes | A-U/G-C rules (Agarwal); Lewis (2005) |
| INV-5 | In `AlignMiRnaToTarget`, Matches + Mismatches + GUWobbles = min(len(miRNA), len(target)); no gaps | Yes | Implementation contract (ungapped duplex over overlap) |
| INV-6 | A fully Watson-Crick-paired duplex has FreeEnergy ≤ 0 (stabilising); an all-mismatch duplex has FreeEnergy not stabilising (≥ 0 from stacking sum) | Yes | NNDB Turner 2004 (negative stacking energies) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CanPair Watson-Crick | A-U,U-A,G-C,C-G | all true | Agarwal et al. (2015) |
| M2 | CanPair wobble | G-U,U-G | true | Crick (1966) |
| M3 | CanPair non-pairs | A-A,A-C,A-G,C-U,C-C,G-G | all false | Agarwal et al. (A pairs only U; C only G) |
| M4 | CanPair case-insensitive + DNA T | a-u, A-T → true | true | normalisation contract; T↔U |
| M5 | IsWobblePair G:U | G-U, U-G | true | PMC4870184 |
| M6 | IsWobblePair Watson-Crick | A-U,G-C | false | PMC4870184 (wobble ≠ WC) |
| M7 | GetReverseComplement simple | `GAGGUAG` (let-7a seed) | `CUACCUC` | Lewis (2005) reverse-complement seed; A-U/G-C |
| M8 | GetReverseComplement DNA | `ACGT` | `ACGU` (RC of ACGT, T→U) | A-U/G-C, T→U |
| M9 | GetReverseComplement empty | `""` | `""` | defensive contract |
| M10 | AlignMiRnaToTarget perfect WC | `AAAA`/`UUUU` | Matches=4, Mismatches=0, GUWobbles=0, alignment `\|\|\|\|` | Watson-Crick A-U |
| M11 | AlignMiRnaToTarget wobble | `GGGG`/`UUUU` | Matches=0, GUWobbles=4, Mismatches=0, alignment `::::` | Crick (1966); G:U ≠ WC |
| M12 | AlignMiRnaToTarget all-mismatch | `AAAA`/`AAAA` | Matches=0, GUWobbles=0, Mismatches=4 | A pairs only U |
| M13 | AlignMiRnaToTarget empty miRNA | `""`/`AAAA` | empty duplex (all zero) | defensive contract |
| M14 | AlignMiRnaToTarget empty target | `AAAA`/`""` | empty duplex (all zero) | defensive contract |
| M15 | AlignMiRnaToTarget count invariant | mixed duplex | Matches+Mismatches+GUWobbles = min(len) | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Free energy sign (paired) | fully WC duplex `GCGCGCGC`/RC | FreeEnergy ≤ 0 | INV-6; Turner negative stacking |
| S2 | Free energy sign (mismatch) | `AAAA`/`AAAA` | FreeEnergy ≥ 0 (not stabilising) | INV-6 |
| S3 | DNA input normalised | `AAAA`/`TTTT` | Matches=4 (T→U) | T↔U normalisation |
| S4 | Unequal lengths | `AAAAAA`/`UUUU` (overlap 4) | alignment length = 4 | INV-5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | RC double application | RC(RC(seq)) == seq for RNA | equal | INV-4 |
| C2 | wobbles ⊆ pairs property | for all base pairs, IsWobblePair⇒CanPair | holds | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzerTests.cs` — has weak/partial tests for `GetReverseComplement`, `CanPair`, `IsWobblePair` (no assertion messages, partial cases).
- `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs` (MIRNA-TARGET-001) — touches `AlignMiRnaToTarget` (M-011..M-015) as smoke for target-site scoring; that is a different unit and stays.
- No canonical file for MIRNA-PAIR-001 exists.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 CanPair WC | ⚠ Weak | MiRnaAnalyzerTests TestCase has no message, partial pair set |
| M2 CanPair wobble | ⚠ Weak | covered by parametrised case w/o message |
| M3 CanPair non-pairs | ⚠ Weak | partial |
| M4 CanPair case/T | ❌ Missing | |
| M5 IsWobblePair G:U | ⚠ Weak | MiRnaAnalyzerTests, no message |
| M6 IsWobblePair WC | ⚠ Weak | MiRnaAnalyzerTests, no message |
| M7 RC let-7a seed | ❌ Missing | only `AUGC` checked elsewhere |
| M8 RC DNA | ❌ Missing | |
| M9 RC empty | ✅ Covered | exists but re-asserted here canonically |
| M10 Align perfect WC | ⚠ Weak | M-011 lacks exact alignment-string check |
| M11 Align wobble | ⚠ Weak | M-012 weak |
| M12 Align all-mismatch | ⚠ Weak | M-013 weak |
| M13 Align empty miRNA | ⚠ Weak | M-014 |
| M14 Align empty target | ⚠ Weak | M-015 |
| M15 count invariant | ❌ Missing | |
| S1 FE sign paired | ⚠ Weak | M-016 exists in TargetPrediction (different unit) |
| S2 FE sign mismatch | ❌ Missing | |
| S3 DNA normalise | ❌ Missing | |
| S4 unequal length | ❌ Missing | |
| C1 RC double | ❌ Missing | |
| C2 wobble⊆pair | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs` — all MIRNA-PAIR-001 cases (CanPair, IsWobblePair, GetReverseComplement, AlignMiRnaToTarget).
- **Remove:** none. The weak duplicates in `MiRnaAnalyzerTests.cs` belong to the legacy aggregate fixture and the `AlignMiRnaToTarget` smoke in `MiRnaAnalyzer_TargetPrediction_Tests.cs` belongs to MIRNA-TARGET-001; both stay (different unit / pre-template). The canonical evidence-based coverage for this unit lives in the new file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MiRnaAnalyzer_AlignMiRnaToTarget_Tests.cs` | Canonical MIRNA-PAIR-001 | 21 |
| `MiRnaAnalyzerTests.cs` | Legacy aggregate (untouched) | unchanged |
| `MiRnaAnalyzer_TargetPrediction_Tests.cs` | MIRNA-TARGET-001 (untouched) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote with exact pairs + messages | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote, full non-pair set | ✅ Done |
| 4 | M4 | ❌ Missing | Added case/T | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewrote | ✅ Done |
| 7 | M7 | ❌ Missing | Added let-7a seed RC | ✅ Done |
| 8 | M8 | ❌ Missing | Added DNA RC | ✅ Done |
| 9 | M9 | ✅ Covered | Re-asserted canonically | ✅ Done |
| 10 | M10 | ⚠ Weak | Rewrote with exact alignment string | ✅ Done |
| 11 | M11 | ⚠ Weak | Rewrote | ✅ Done |
| 12 | M12 | ⚠ Weak | Rewrote | ✅ Done |
| 13 | M13 | ⚠ Weak | Rewrote | ✅ Done |
| 14 | M14 | ⚠ Weak | Rewrote | ✅ Done |
| 15 | M15 | ❌ Missing | Added count invariant | ✅ Done |
| 16 | S1 | ⚠ Weak | Rewrote (sign only) | ✅ Done |
| 17 | S2 | ❌ Missing | Added | ✅ Done |
| 18 | S3 | ❌ Missing | Added | ✅ Done |
| 19 | S4 | ❌ Missing | Added | ✅ Done |
| 20 | C1 | ❌ Missing | Added | ✅ Done |
| 21 | C2 | ❌ Missing | Added property test | ✅ Done |

**Total items:** 21
**✅ Done:** 21 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | CanPair_WatsonCrickPairs_ReturnsTrue |
| M2 | ✅ | CanPair_WobblePairs_ReturnsTrue |
| M3 | ✅ | CanPair_NonPairs_ReturnsFalse |
| M4 | ✅ | CanPair_LowercaseAndDnaT_ReturnsTrue (now asserts CanPair('A','T')==true per contract §3.1/§6.1) |
| M4b | ✅ | IsWobblePair_DnaT_TreatedAsU (G-T/T-G wobble; A-T not wobble) |
| M5 | ✅ | IsWobblePair_GU_ReturnsTrue |
| M6 | ✅ | IsWobblePair_WatsonCrick_ReturnsFalse |
| M7 | ✅ | GetReverseComplement_Let7aSeed_ReturnsExpected |
| M8 | ✅ | GetReverseComplement_DnaInput_ReturnsRnaReverseComplement |
| M9 | ✅ | GetReverseComplement_Empty_ReturnsEmpty |
| M10 | ✅ | AlignMiRnaToTarget_PerfectComplement_AllWatsonCrickMatches |
| M11 | ✅ | AlignMiRnaToTarget_GUWobbleDuplex_CountedAsWobbles |
| M12 | ✅ | AlignMiRnaToTarget_NonPairingDuplex_AllMismatches |
| M13 | ✅ | AlignMiRnaToTarget_EmptyMiRna_ReturnsEmptyDuplex |
| M14 | ✅ | AlignMiRnaToTarget_EmptyTarget_ReturnsEmptyDuplex |
| M15 | ✅ | AlignMiRnaToTarget_CountInvariant_SumEqualsOverlapLength |
| S1 | ✅ | AlignMiRnaToTarget_FullyPairedDuplex_FreeEnergyNonPositive |
| S2 | ✅ | AlignMiRnaToTarget_AllMismatchDuplex_FreeEnergyNonNegative |
| S3 | ✅ | AlignMiRnaToTarget_DnaInput_NormalisedAndMatched |
| S4 | ✅ | AlignMiRnaToTarget_UnequalLengths_AlignsOverOverlap |
| C1 | ✅ | GetReverseComplement_DoubleApplication_RestoresRna |
| C2 | ✅ | WobblePairs_AreSubsetOfPairs_Property |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Turner 2004 stacking numeric values not re-retrieved this session (NNDB 404 to fetcher); free-energy tests assert sign only, not magnitude | S1, S2, INV-6 |

---

## 7. Open Questions / Decisions

1. The `FreeEnergy` magnitude is "Intentionally simplified" — refactored to sum the already-cited NNDB Turner 2004 nearest-neighbor stacking table (removing the previously invented −2.0/−1.0/−0.5/+0.5 per-position constants). Only its sign is validated, because exact per-stack kcal/mol values could not be independently re-opened this session.
