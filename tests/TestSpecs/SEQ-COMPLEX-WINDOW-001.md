# Test Specification: SEQ-COMPLEX-WINDOW-001

**Test Unit ID:** SEQ-COMPLEX-WINDOW-001
**Area:** Complexity
**Algorithm:** Windowed Sequence Complexity (sliding-window complexity profile)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-30

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Shannon (1948), A Mathematical Theory of Communication | 1 | https://doi.org/10.1002/j.1538-7305.1948.tb01338.x (via https://en.wikipedia.org/wiki/Entropy_(information_theory)) | 2026-06-14 |
| 2 | Troyanskaya et al. (2002), Sequence complexity profiles … fast algorithm | 1 | https://doi.org/10.1093/bioinformatics/18.5.679 (https://pubmed.ncbi.nlm.nih.gov/12050064/) | 2026-06-14 |
| 3 | Gabrielian & Bolshoy (1999), Sequence complexity and DNA curvature | 1 | https://doi.org/10.1016/S0097-8485(99)00007-8 (via Wikipedia LC) | 2026-06-14 |
| 4 | Wikipedia, Linguistic sequence complexity (cites Trifonov 1990, Troyanskaya 2002) | 4 | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | 2026-06-14 |
| 5 | BBTools BBMap 40.02 — `jgi/BBDuk.java` `maskLowEntropy`, `tracker/EntropyTracker.java`, `dna/AminoAcid.java` (source in the release `bbtools.jar`; release binary run as the oracle) | 2 (reference implementation) | https://sourceforge.net/projects/bbmap/ ; https://raw.githubusercontent.com/BioInfoTools/BBMap/master/current/structures/EntropyTracker.java | 2026-09-30 |
| 6 | Gabrielian & Bolshoy (1999) abstract (search snippet): LC word length "not in the range of 2 to N−1 but only up to W"; universalmotif `R/sequence_complexity.R` (`trifonov.max.word.size = 7`) | 1 / 3 | https://www.sciencedirect.com/science/article/abs/pii/S0097848599000078 ; https://raw.githubusercontent.com/bjmt/universalmotif/master/R/sequence_complexity.R | 2026-09-30 |

### 1.2 Key Evidence Points

1. Windowed complexity profiling = computing complexity as a function of position along a sequence using a sliding window — Troyanskaya et al. (2002).
2. Per-window Shannon entropy `H = -Σ p log₂ p` (bits); uniform DNA window ⇒ 2.0, homopolymer ⇒ 0 — Shannon (1948) via Wikipedia Entropy.
3. Per-window linguistic complexity (summation form) `LC = (Σ Vᵢ)/(Σ Vmax,i)`, `Vmax,i = min(4^i, N-i+1)`, range (0,1) — Wikipedia LC / Troyanskaya (2002); repo LC unit SEQ-COMPLEX-001.
4. Worked window `ACGTACGT` (maxWordLength=6): H=2.0, LC=23/29. Window `AAAAAAAA`: H=0.0, LC=6/29 — hand-derivation in Evidence §Test Datasets.
6. BBDuk low-entropy masking: window statistic = entropy of the w−k+1 overlapping k-mers normalised by ln(w−k+1) (defaults k=5, w=50), window fails when value < cutoff (float), only windows with no undefined base (`ns()<1`; A/C/G/T/U either case defined), masked set = union of failing windows, reads shorter than w untouched — BBDuk/EntropyTracker source (5).
7. Windowed LC word length is bounded by a free parameter W (Gabrielian & Bolshoy 1999); 6 is this library's default, now the `lcMaxWordLength` parameter — source (6).
5. Window enumeration emits only fully-contained windows (`i + w ≤ L`), advancing by `step`; count = floor((L-w)/s)+1 for L≥w — repository contract + profile definition.

### 1.3 Documented Corner Cases

- Uniform window ⇒ Shannon = log₂4 = 2.0 (Shannon 1948).
- Homopolymer window ⇒ Shannon = 0 (deterministic distribution, Shannon 1948).
- L < windowSize ⇒ no window emitted (empty profile); trailing partial fragment not emitted (repository contract / profile definition).

### 1.4 Known Failure Modes / Pitfalls

1. Off-by-one in window count or in WindowEnd (inclusive vs exclusive) — repository contract (0-based, end inclusive).
2. Per-window LC depends on the word-length cap `min(lcMaxWordLength, w)` (default 6) — Gabrielian & Bolshoy bounded-W form.
3. Treating the Shannon window scan as BBDuk's default rule: BBDuk scores normalised 5-mer entropy; the Shannon scan equals BBDuk only with `entropyk=1 entropy=t/log₂w` — source (5).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateWindowedComplexity(DnaSequence, int windowSize, int stepSize)` | SequenceComplexity | **Canonical** | Sliding-window driver returning `ComplexityPoint` per window |
| `CalculateWindowedComplexity(string, int, int, int lcMaxWordLength = 6)` | SequenceComplexity | **Canonical** | String overload; windows with a non-ACGTU symbol skipped (BBDuk `ns()<1`) |
| `FindLowComplexityRegions(DnaSequence \| string, int windowSize, double entropyThreshold)` | SequenceComplexity | **Canonical** | Region = union (maximal covered run) of step-1 windows with per-base Shannon entropy < threshold (BBDuk window-union reporting, 1-mer statistic); string overload never flags windows with non-ACGTU; threshold NaN/∞/<0 rejected; tests in `SequenceComplexityTests.cs` `#region Low Complexity Region Tests` + `SequenceComplexity_AuditWp5_Tests.cs` |
| `FindLowEntropyRegionsBbduk(string, double entropyCutoff, int windowSize = 50, int k = 5)` | SequenceComplexity | **Canonical** | Port of BBDuk `maskLowEntropy` (`entropymask=t`); tests in `SequenceComplexity_AuditWp5_Tests.cs` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Number of points = floor((L-w)/s)+1 for L≥w, else 0 | Yes | Sliding-window geometry (Troyanskaya 2002; repo contract) |
| INV-2 | For every point: WindowStart = i, WindowEnd = i+w-1, Position = i+w/2 (int division) | Yes | `ComplexityPoint` contract |
| INV-3 | 0 ≤ ShannonEntropy ≤ log₂4 = 2.0 per window | Yes | Shannon (1948): bounds 0..log₂(n) |
| INV-4 | 0 < LinguisticComplexity ≤ 1 for DNA windows | Yes | Wikipedia LC range (0,1) |
| INV-5 | Windows are non-overlapping iff stepSize ≥ windowSize; consecutive WindowStart differ by stepSize | Yes | repo contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Window count | L=24, w=8, s=8 | 3 points | INV-1; Evidence §Window enumeration geometry |
| M2 | Coordinates | L=24, w=8, s=8 | starts 0,8,16; ends 7,15,23; positions 4,12,20 | INV-2; Evidence §Window enumeration geometry |
| M3 | Uniform window entropy | window `ACGTACGT` | ShannonEntropy = 2.0 | Shannon (1948), uniform max = log₂4 |
| M4 | Homopolymer window entropy | window `AAAAAAAA` | ShannonEntropy = 0.0 | Shannon (1948), deterministic = 0 |
| M5 | Uniform window LC | window `ACGTACGT`, maxWord=6 | LinguisticComplexity = 23/29 = 0.7931034482758621 | Wikipedia LC; Evidence dataset |
| M6 | Homopolymer window LC | window `AAAAAAAA`, maxWord=6 | LinguisticComplexity = 6/29 = 0.20689655172413793 | Wikipedia LC; Evidence dataset |
| M7 | L < windowSize | L=5, w=8 | empty profile (0 points) | Corner case (no partial window); repo contract |
| M8 | Null DnaSequence | null input | ArgumentNullException | repository contract |
| M9 | windowSize < 1 | w=0 | ArgumentOutOfRangeException | repository contract |
| M10 | stepSize < 1 | s=0 | ArgumentOutOfRangeException | repository contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Overlapping step | L=24, w=8, s=4 | starts 0,4,8,12,16; count = floor(16/4)+1 = 5 | INV-1/INV-5 with s<w |
| S2 | Exact-fit single window | L=8, w=8, s=8 | exactly 1 point, start 0 end 7 | boundary L=w |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Bounds invariant | mixed sequence | every point: 0≤H≤2 and 0<LC≤1 | INV-3/INV-4 |

### 4.4 FindLowComplexityRegions (added 2026-09-29, review-2026-09 B04)

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| R1 | ATGC×20 + A×64 + ATGC×20, w=20, thr=0.5 | one region 79..145, length 67, MinEntropy 0 | BBDuk window union; Python reference |
| R2 | ATGCATGC + A×64 + ATGCATGC, w=32, thr=0.5 | one region 6..74, length 69 | idem |
| R3 | CAAAAACAAAAACAAACAAA, w=8, thr=0.6 | overlapping flagged windows merge: one region 1..15, MinEntropy 0.5435644431995964 | idem |
| R4 | A×10 + (ACGT)×3 + C×10, w=8, thr=1.0 | two disjoint regions 0..11 and 21..31 (trailing) | idem |
| R5 | null / windowSize 0 | ArgumentNullException / ArgumentOutOfRangeException at call time (eager) | repository contract |
| R6 | differential windowed profile (L=50, w=16, s=7) | per-window H and LC equal scipy/set reference to 1e-12 | Evidence §Python reference |

### 4.5 Completeness audit WP5 (2026-09-30, B04 F38/F39) — `SequenceComplexity_AuditWp5_Tests.cs`

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| W1 | `ACGTACGT`, w=8, lcMaxWordLength 8 / 1 / default | LC 13/16 / 1.0 / 23/29; H unchanged | brute force |
| W2 | string `acgtacgtaaNaaaaaacguacgt`, w=8, s=4 | points at 0, 12, 16 only (N windows skipped); U = T | source (5) `ns()<1` |
| W3 | string overload vs `DnaSequence` on ACGT (lower-cased input) | identical points / regions | contract |
| W4 | windowed `lcMaxWordLength 0`, string w/s 0, null | ArgumentOutOfRange / ArgumentNull (eager) | contract |
| W5 | `FindLowComplexityRegions` threshold NaN, ±∞, −0.1 (both overloads) | ArgumentOutOfRangeException | contract |
| W6 | ATGC×5 + A×20 + N + A×20 + ATGC×5, w=10, t=0.5 | regions 19..39 and 41..62 (N windows never flagged) | brute force |
| B1 | 30 bp + A×40 + 30 bp, defaults, cutoff 0.5 (also lower-case RNA spelling) | masked 11..88 | `bbduk.sh` 40.02 output |
| B2 | 40 bp + (AC)×40 + 30 bp, cutoff 0.7 / w=20 k=3 cutoff 0.6 | 17..143 / 34..125 | `bbduk.sh` output |
| B3 | 40 bp + (AAG)×10 + 30 bp, w=25 k=2 cutoff 0.55 | 32..74 | `bbduk.sh` output |
| B4 | MinEntropy of B1 / B2 | (float) 0.2674965368023891 / 0.1810425967800402 | −Σp ln p / ln 46 |
| B5 | N inside the 50-bp poly-A; 49-bp read with cutoff 1; empty | nothing masked | `bbduk.sh` output |
| B6 | k = 1, cutoff t/log₂w ≡ Shannon scan (w 8/16/20, t 0.5/1/1.3) | identical regions | §2.3 identity; 1 600 cases vs bbduk k=1 |
| B7 | null, cutoff −0.01 / 1.01 / NaN, k 0 / 16, w ≤ k | ArgumentNull / ArgumentOutOfRange | EntropyTracker assertions |
| X1 | differential vs compiled-release `bbduk.sh` 40.02 (harness, not a unit test) | 4 500 random cases + 20 long reads: 0 mismatches | Evidence 2026-09-30 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexityTests.cs` contains a `#region Windowed Complexity Tests` with three tests (`CalculateWindowedComplexity_ReturnsCorrectPointCount`, `_IncludesBothMetrics_ExactValues`, `_PositionsAreCorrect`) plus three failure-mode tests (`_NullSequence_ThrowsException`, `_ZeroWindowSize_ThrowsException`, `_ZeroStepSize_ThrowsException`).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (count) | ⚠ Weak | existing count test exists but no exact LC/coordinate triangulation; superseded by canonical file |
| M2 (coordinates) | ⚠ Weak | existing positions test only checks 2 positions, not start/end inclusivity |
| M3 (uniform entropy) | ⚠ Weak | existing test asserts entropy 2.0 but LC only `GreaterThan(0)` (permissive) |
| M4 (homopolymer entropy) | ❌ Missing | not covered |
| M5 (uniform LC exact) | ❌ Missing | existing uses `GreaterThan(0)` — permissive, not exact |
| M6 (homopolymer LC exact) | ❌ Missing | not covered |
| M7 (L<w empty) | ❌ Missing | not covered |
| M8 (null) | 🔁 Duplicate | exists in old file; moves to canonical file |
| M9 (windowSize<1) | 🔁 Duplicate | exists in old file; moves to canonical file |
| M10 (stepSize<1) | 🔁 Duplicate | exists in old file; moves to canonical file |
| S1 (overlap) | ❌ Missing | not covered |
| S2 (exact fit) | ❌ Missing | not covered |
| C1 (bounds) | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateWindowedComplexity_Tests.cs` — all M/S/C cases with exact evidence-derived values.
- **Remove:** the six windowed tests (3 in `#region Windowed Complexity Tests` + 3 failure-mode windowed tests) from `SequenceComplexityTests.cs` to avoid duplication; their region(s) are consolidated into the canonical file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceComplexity_CalculateWindowedComplexity_Tests.cs` | Canonical for this unit | 13 |
| `SequenceComplexityTests.cs` | Other complexity methods (windowed tests removed) | (unchanged minus 6) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewritten exact in canonical file | ✅ Done |
| 2 | M2 | ⚠ Weak | rewritten exact (start/end/position) | ✅ Done |
| 3 | M3 | ⚠ Weak | rewritten exact entropy 2.0 | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented exact 23/29 | ✅ Done |
| 6 | M6 | ❌ Missing | implemented exact 6/29 | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | 🔁 Duplicate | moved to canonical; removed from old file | ✅ Done |
| 9 | M9 | 🔁 Duplicate | moved to canonical; removed from old file | ✅ Done |
| 10 | M10 | 🔁 Duplicate | moved to canonical; removed from old file | ✅ Done |
| 11 | S1 | ❌ Missing | implemented | ✅ Done |
| 12 | S2 | ❌ Missing | implemented | ✅ Done |
| 13 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact count 3 |
| M2 | ✅ | starts/ends/positions exact |
| M3 | ✅ | entropy 2.0 exact |
| M4 | ✅ | entropy 0.0 exact |
| M5 | ✅ | LC 23/29 exact |
| M6 | ✅ | LC 6/29 exact |
| M7 | ✅ | empty profile |
| M8 | ✅ | ArgumentNullException |
| M9 | ✅ | ArgumentOutOfRangeException |
| M10 | ✅ | ArgumentOutOfRangeException |
| S1 | ✅ | overlap count 5 |
| S2 | ✅ | single window |
| C1 | ✅ | bounds invariant |

Total in-scope cases: 13; ✅ = 13.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Position = WindowStart + windowSize/2 (center, integer division) — non-correctness-affecting label | M2 |
| 2 | Default windowSize=64, stepSize=10 (caller-overridable; not value-affecting for an explicit window) | parameters |
| 3 | Per-window LC uses maxWordLength = min(6, windowSize) (repo LC efficiency cap, Gabrielian & Bolshoy) | M5, M6 |

---

## 7. Open Questions / Decisions

1. None. The method is a driver over the already-evidenced Shannon entropy (SEQ-COMPLEX-*) and linguistic complexity (SEQ-COMPLEX-001) units; per-window values are derived from those established formulas.
