# META-BIN-001: Genome Binning - Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | META-BIN-001 |
| **Algorithm** | Genome Binning (MAG assembly) |
| **Canonical Method** | `MetagenomicsAnalyzer.BinContigs(...)` |
| **Complexity** | O(n × k × i) |
| **Evidence Document** | [META-BIN-001-Evidence.md](../docs/Evidence/META-BIN-001-Evidence.md) |

## Method Under Test

```csharp
public static IEnumerable<GenomeBin> BinContigs(
    IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
    int numBins = 10,
    double minBinSize = 500000,
    double expectedGenomeSize = 4_000_000)
```

### Output Record

```csharp
public readonly record struct GenomeBin(
    string BinId,
    IReadOnlyList<string> ContigIds,
    double TotalLength,
    double GcContent,
    double Coverage,
    double Completeness,
    double Contamination,
    string PredictedTaxonomy);
```

## Test Categories

### Must Tests (Required for DoD)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | EmptyInput_ReturnsEmpty | Empty contigs must yield empty result | Standard edge case |
| M2 | SingleContig_BelowMinSize_ReturnsEmpty | Contigs below minBinSize excluded | Implementation behavior |
| M3 | ValidContigs_ReturnsNonEmptyBins | Non-trivial input must produce bins | Wikipedia |
| M4 | BinIds_AreUnique | Each bin must have unique ID | Invariant |
| M5 | Completeness_InValidRange | Must be [0, 100]; verified exact formula: min(totalLength/expectedGenomeSize×100, 100) | CheckM standard |
| M6 | Contamination_InValidRange | Must be [0, 100]; verified: uniform GC→0, max-variance GC→100 | CheckM standard |
| M7 | GcContent_InValidRange | Must be [0, 1]; verified exact values for pure high-GC (1.0) and low-GC (0.0) bins | Mathematical invariant |
| M8 | TotalLength_EqualsContigLengthSum | Length accounting correct | Mathematical invariant |
| M9 | MinBinSize_FiltersSmallBins | Bins below threshold excluded | Implementation spec |
| M10 | HighGcVsLowGc_SeparatesIntoDifferentBins | K-means separates extreme GC populations | Wikipedia, TETRA paper |
| M11 | ContigIds_PreservedInBins | All binned contigs traceable; bins disjoint | Data integrity |
| M12 | Coverage_PreservedInBins | Coverage = arithmetic mean of contig coverages; verified exact value | Implementation spec |

### Should Tests (Recommended)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | MultipleBins_FromDiverseInput | Diverse contigs create multiple bins | Wikipedia |
| S2 | NumBins_AffectsBinCount | k-means with k clusters produces at most k bins; verified both limits | Implementation |
| S3 | LargeDataset_CompletesWithinTimeout | Performance acceptable | Practical requirement |

### Could Tests (Optional)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C1 | AllSameGc_ClustersTogether | Uniform features converge to few bins |
| C2 | VeryShortSequences_HandledGracefully | Robustness |

## Invariants

1. **Empty Preservation**: `BinContigs([]) → []`
2. **Unique IDs**: `∀ bin1, bin2 ∈ result: bin1.BinId ≠ bin2.BinId`
3. **Metric Ranges**:
   - Completeness ∈ [0, 100]
   - Contamination ∈ [0, 100]
   - GcContent ∈ [0, 1]
4. **Size Filtering**: `∀ bin ∈ result: bin.TotalLength ≥ minBinSize`
5. **Length Accuracy**: `bin.TotalLength = Σ len(contig) for contig in bin`

## Test Data Design

### High GC Contigs (for separation test)
- Sequence: 70% G/C content (e.g., `GGGGGGGGCC` repeated)
- Expected: Group together in high-GC bin

### Low GC Contigs (for separation test)
- Sequence: 30% G/C content (e.g., `AAAATTTTGC` repeated)
- Expected: Group together in low-GC bin

### Mixed Contigs
- Combination of high and low GC
- Expected: Separate into different bins

## Existing Test Status

| Location | Status | Action |
|----------|--------|--------|
| MetagenomicsAnalyzer_GenomeBinning_Tests.cs | All Must/Should/Could tests implemented | Complete |
| MetagenomicsAnalyzerTests.cs | Smoke test removed (duplicate of M3) | Cleaned |
| MetagenomicsSnapshotTests.cs | Snapshot test for deterministic output | Independent |

## Coverage Classification Summary

| Classification | Count | Details |
|---------------|-------|----------|
| ✅ Covered | 11 | M1, M2, M4, M8, M9, M10, S1, S3, C1, C2, extra:ZeroCoverage |
| ⚠ Strengthened | 7 | M3 (Not.Empty), M5 (exact formula), M6 (exact formula), M7 (exact GC), M11 (disjoint), M12 (exact average), S2 (both limits) |
| 🔁 Removed | 1 | BinContigs_BasicFunctionality_ReturnsResults (subsumed by M3) |
| ❌ Missing | 0 | None |

## Audit Notes

- Current smoke tests adequate for basic functionality
- Need comprehensive invariant testing
- Need GC-based separation validation
- Need quality metric range validation

## Consolidation Plan

All tests implemented and verified in `MetagenomicsAnalyzer_GenomeBinning_Tests.cs`.
Duplicate smoke test in `MetagenomicsAnalyzerTests.cs` has been removed.
Snapshot test in `MetagenomicsSnapshotTests.cs` remains independent.

---

# Addendum (2026-06-24): TETRA z-score tetranucleotide signature (opt-in)

Covers the new opt-in methods `CalculateTetranucleotideZScores` and
`TetranucleotideZScoreCorrelation`. The default `BinContigs` raw-frequency path is unchanged.

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Teeling et al. 2004, TETRA (BMC Bioinformatics 5:163) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC529438/ | 2026-06-24 |
| 2 | Teeling et al. 2004, Environ Microbiol 6(9):938–947 | 1 | https://doi.org/10.1111/j.1462-2920.2004.00624.x | 2026-06-24 |
| 3 | Schbath 1997, J Comput Biol 4(2):189–192 (variance) | 1 | https://pubmed.ncbi.nlm.nih.gov/9228617/ | 2026-06-24 |
| 4 | Bohlin & Skjerve 2009 (PLOS ONE, expected-count form) | 1 | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008113 | 2026-06-24 |

### 1.2 Key Evidence Points

1. Expected count `E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4)/N(n2n3)` (maximal-order Markov) — source 1, 4.
2. Variance `var = E·[N(n2n3)−N(n1n2n3)][N(n2n3)−N(n2n3n4)]/N(n2n3)²`; z `=(N−E)/√var` — source 1, 3.
3. Sequence is extended by its reverse complement; signatures compared by Pearson correlation — source 1.

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateTetranucleotideZScores(string)` | MetagenomicsAnalyzer | Canonical | 256-component TETRA z-score signature |
| `TetranucleotideZScoreCorrelation(string,string)` | MetagenomicsAnalyzer | Canonical | Pearson r of two z-score vectors |

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-Z1 | Signature has exactly 256 components | Yes | source 1 (4^4 tetramers) |
| INV-Z2 | Self-correlation = 1.0 | Yes | Pearson identity |
| INV-Z3 | z=0 when N(n2n3)=0 or var≤0 | Yes | formula domain (source 1, 3) |

## 4. Test Cases

### 4.1 MUST Tests

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M-Z1 | Hand-derived z(ACGT) | `ACGTACGTGGCC`: E=3.2, var=0.128 | z = √5 = 2.2360679774997896 (±1e-10) | §1.2 pt 1–2 + Evidence dataset |
| M-Z2 | 256 components | any sequence | `Count==256`, contains AAAA/ACGT/TTTT | source 1 |
| M-Z3 | Absent middle dinucleotide | `AAAAAAAA`, tetramer ACGT (N(CG)=0) | z(ACGT)=0 | INV-Z3 |
| M-Z4 | Null/empty/single-base | null, "", "A" | all-zero 256-vector | corner cases |
| M-Z5 | Non-ACGT / case filtering | `acgtACGTGGCC` vs `ACGTACGTGGCC` | identical z(ACGT) | source 1 (ACGT words) |
| M-Z6 | Self-correlation | corr(s,s) | 1.0 (±1e-10) | INV-Z2 |
| M-Z7 | Similar > dissimilar | s1~s2 vs s1 vs AT-homopolymer | r_similar > r_dissimilar | source 1 (binning basis) |

### 4.2 SHOULD Tests

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| M-Z8 | Symmetry | corr(a,b)=corr(b,a) | equal (±1e-12) | Pearson symmetry |
| S-Z1 | Degenerate vs empty | corr(s,"") | 0, not NaN | zero-variance guard |

## 5. Audit / Coverage

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TetranucleotideZScore_Tests.cs`
- All MUST/SHOULD cases (M-Z1..M-Z8, S-Z1) implemented and ✅ Covered (9 tests). Remaining ❌/⚠ = 0.

## 6. Assumption Register

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Counting on the concatenated RC-extended strand reproduces the published E/var | M-Z1 |

## 7. Open Questions / Decisions

1. CheckM marker-gene completeness/contamination is now **implemented** as an opt-in addendum
   (see §8). The residual is narrowed to the full CheckM lineage-specific marker DB + reference
   tree (large gated data); a small CC0 universal ribosomal marker set + caller-supplied loader
   ship instead.

---

## 8. Addendum — CheckM-style marker-gene completeness/contamination (opt-in)

**Evidence:** [META-BIN-001-MarkerQC-Evidence.md](../../docs/Evidence/META-BIN-001-MarkerQC-Evidence.md)
**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`

### 8.1 Methods Under Test

| Method | Type | Notes |
|--------|------|-------|
| `EstimateBinQualityFromMarkerCounts(markerSets, markerCounts)` | Canonical | exact CheckM Eqs.1–2 |
| `DetectMarkers(proteins, markerHmms)` | Canonical | Plan7 Viterbi-bits ≥ GA1 ⇒ copy count |
| `EstimateBinQualityFromMarkers(proteins, markerSets, markerHmms)` | Canonical | detect + formula |
| `LoadBundledRibosomalMarkerHmms()` / `BundledRibosomalMarkerSets()` | Canonical | 9 CC0 markers |
| `LoadMarkerHmms(readers, thresholds?)` | Canonical | caller-supplied loader |

### 8.2 Invariants

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-M1 | Completeness = 100·(1/|M|)·Σ |s∩G_M|/|s| | Parks 2015 Eq.1 / `genomeCheck` |
| INV-M2 | Contamination = 100·(1/|M|)·Σ Σ(N_g−1)/|s| | Parks 2015 Eq.2 / `genomeCheck` |
| INV-M3 | A multi-copy marker counts once toward completeness | `genomeCheck` |
| INV-M4 | Empty marker sets / |M|=0 ⇒ no div-by-zero (0/0) | corner cases |

### 8.3 Test Cases (all MUST, all ✅ Covered)

| ID | Test Case | Expected | Evidence |
|----|-----------|----------|----------|
| M-FORMULA | synthetic bin {A,B},{C,D,E},{F}; A1,B0,C2,D1,E1,F1 | Comp 250/3 %, Cont 100/9 % | Eqs.1–2 hand-derived |
| M-ALLPRESENT | all single-copy present | 100% / 0% | Eq.1–2 |
| M-ALLABSENT | none found | 0% / 0% | Eq.1–2 |
| M-TRIPLICATE | {A}, A=3 | 100% / 200% | Eq.2 C_g=N−1 |
| M-EMPTYSET | one empty set | excluded from |M| | corner case 3 |
| M-NOMARKERSETS | |M|=0 | 0% / 0% | corner case 4 |
| H-TRUEPOSITIVE | bundled PF00410 vs E. coli uS8 | hit; 0 for all 8 others | UniProt P0A7W7 |
| H-COMPLETENESS | bin = {uS8}, 9 singleton sets | Comp 100/9 %, Cont 0 | integration |
| H-DUPLICATE | bin = {uS8,uS8} | Comp 100/9 %, Cont 100/9 % | Eq.2 |
| L-CALLERSUPPLIED | reader-loaded S8 keyed by ACC, detects uS8 | PF00410.25 = 1 | loader |
| (guards) | null args throw; loader counts | ArgumentNullException; 9 markers, GA1=24 | — |

**Coverage:** 14 tests, all ✅; Remaining ❌/⚠ = 0.

---

**Specification Version:** 4.0
**Created:** 2026-02-04
**Updated:** 2026-06-25 (CheckM marker-gene completeness/contamination implemented, opt-in;
residual narrowed to the full lineage DB + tree)
**Status:** Complete — TNF z-score + CheckM marker QC implemented, evidence-based, and verified
