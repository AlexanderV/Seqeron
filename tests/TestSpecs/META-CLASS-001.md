# META-CLASS-001: Taxonomic Classification Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-CLASS-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.ClassifyReads`, `MetagenomicsAnalyzer.BuildKmerDatabase` |
| **Complexity** | O(n × m) where n=reads, m=read length |
| **Invariant** | |output| = |input reads|; 0 ≤ Confidence ≤ 1; Confidence = C/Q per Kraken |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `ClassifyReads(reads, kmerDB, k)` | MetagenomicsAnalyzer | Canonical | Deep |
| `BuildKmerDatabase(refGenomes, k)` | MetagenomicsAnalyzer | DB construction | Deep |

## Evidence Sources

1. **Wikipedia — Metagenomics:** K-mer based binning, taxonomic classification principles
2. **Kraken 1 Manual (CCB JHU):** Canonical k-mer classification algorithm, k=31 default, ambiguous k-mer filtering
3. **Kraken 2 Manual (GitHub Wiki):** Confidence formula C/Q where C=clade k-mers, Q=non-ambiguous k-mers
4. **Wood & Salzberg (2014):** Kraken paper describing exact k-mer matching for classification

## Test Categories

### MUST Tests (Evidence-Backed)

#### BuildKmerDatabase

| ID | Test | Evidence |
|----|------|----------|
| M1 | Empty input returns empty database | Standard robustness |
| M2 | Sequence shorter than k produces no k-mers | Kraken manual: cannot extract k-mers |
| M3 | Valid reference produces k-mers in database | Core functionality (Kraken) |
| M4 | Uses canonical k-mers (lexicographically smaller) | Kraken paper: database size reduction |
| M5 | Non-ACGT characters excluded from k-mers | Kraken: only ACGT k-mers indexed |
| M6 | Database k-mer count = Σ(len - k + 1) for valid sequences | K-mer counting invariant |

#### ClassifyReads

| ID | Test | Evidence |
|----|------|----------|
| M7 | Empty sequence returns Unclassified with Confidence=0 | Kraken behavior |
| M8 | Sequence shorter than k returns Unclassified | Kraken: cannot extract k-mers |
| M9 | No database matches returns Unclassified | Kraken: no hits → unclassified |
| M10 | Matching k-mers classify to correct taxon | Core classification logic (Kraken) |
| M11 | Output count equals input read count | Output invariant |
| M12 | Confidence = C/Q (C=winning taxon k-mers, Q=non-ambiguous k-mers) | Kraken 1&2 confidence formula |
| M13 | TotalKmers = non-ambiguous k-mers (= len-k+1 for all-ACGT) | Kraken: Q excludes ambiguous k-mers |
| M14 | MatchedKmers ≤ TotalKmers; MatchedKmers = winning taxon count | Bound invariant + Kraken C/Q |
| M15 | Multiple reads all classified | Batch processing |
| M16 | Taxonomy string parsed correctly | Implementation: pipe/semicolon delimited |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Mixed case input handled (uppercased internally) | Robustness for real data |
| S2 | Multiple taxon matches resolves to highest count | LCA-like behavior |
| S3 | Canonical k-mer lookup in ClassifyReads (RC read matches canonical DB) | Kraken: canonicalize before lookup |
| S4 | Ambiguous nucleotides excluded from TotalKmers | Kraken: Q = non-ambiguous k-mers only |
| S5 | Multi-taxon confidence uses winning taxon count (C/Q) | Kraken 1&2: C = clade k-mers |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Performance with large database acceptable | Practical use case |
| C2 | Memory-efficient for large read sets | Streaming enumerable |

## Coverage Classification

All tests consolidated into `MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs` (27 tests).

### BuildKmerDatabase

| ID | Test | Status | Notes |
|----|------|--------|-------|
| M1 | EmptyInput_ReturnsEmptyDatabase | ✅ Covered | Exact assertion |
| M2 | SequenceShorterThanK_ReturnsEmpty | ✅ Covered | Uses k=31 |
| M3 | ValidReference_ProducesKmers | ✅ Covered | Exact count=6, all values verified (k=5, "ACGTACTGAC") |
| M4 | UsesCanonicalKmers | ✅ Covered | Forward A→canonical A, forward T→canonical A |
| M4b | CanonicalKmer_UsesReverseComplementWhenSmaller | ✅ Covered | Complementary to M4 |
| M5 | NonAcgtCharacters_Excluded | ✅ Covered | Regex validates all keys ACGT-only |
| M6 | KmerCount_FollowsFormula | ✅ Covered | Exact count=6=len-k+1 for non-repeating sequence |
| S1 | MixedCase_NormalizedToUppercase | ✅ Covered | Lowercase input, uppercase keys verified |

### ClassifyReads

| ID | Test | Status | Notes |
|----|------|--------|-------|
| M7 | EmptySequence_ReturnsUnclassified | ✅ Covered | All fields: Kingdom, Confidence, TotalKmers |
| M8 | SequenceShorterThanK_ReturnsUnclassified | ✅ Covered | Kingdom + TotalKmers=0 |
| M9 | NoMatch_ReturnsUnclassified | ✅ Covered | MatchedKmers=0, TotalKmers>0 |
| M10 | MatchingKmers_ClassifiesCorrectly | ✅ Covered | Exact: MatchedKmers=1, TotalKmers=8, Confidence=0.125 |
| M11 | OutputCountEqualsInputCount | ✅ Covered | 5 reads (matched, empty, short) |
| M12 | ConfidenceCalculation_IsCorrect | ✅ Covered | Exact: 1/2=0.5 (1 match out of 2 k-mers) |
| M13 | TotalKmers_MatchesFormula | ✅ Covered | 4 cases with exact values |
| M14 | MatchedKmers_BoundedByTotal | ✅ Covered | MatchedKmers ≤ TotalKmers |
| M15 | MultipleReads_AllClassified | ✅ Covered | Exact: Bacteria, Bacteria, Unclassified |
| M16 | TaxonomyParsing_PipeDelimited | ✅ Covered | All 7 ranks verified |
| M16b | TaxonomyParsing_SemicolonDelimited | ✅ Covered | 4 ranks verified |
| S1 | MixedCaseInput_Handled | ✅ Covered | Exact: Kingdom=Bacteria, MatchedKmers=2 |
| S2 | MultipleTaxonMatches_ResolvesToHighestCount | ✅ Covered | Phylum=Taxon2 (2 hits vs 1) |
| S3 | CanonicalKmerLookup_MatchesReverseComplement | ✅ Covered | RC lookup verified, MatchedKmers=2, Confidence=1.0 |
| S4 | AmbiguousNucleotides_ExcludedFromTotalKmers | ✅ Covered | N at pos 9: TotalKmers=1, Confidence=1.0 |
| S5 | MultiTaxon_ConfidenceUsesWinningTaxonCount | ✅ Covered | C=3, Q=4, Confidence=0.75 (Kraken formula) |

### Invariants

| Test | Status | Notes |
|------|--------|-------|
| AllOutputs_HaveValidConfidence | ✅ Covered | Range [0,1] for 4 diverse reads |
| UnclassifiedReads_HaveZeroMatchedKmers | ✅ Covered | Invariant: no hits → MatchedKmers=0 |
| ReadIdPreserved | ✅ Covered | 2 reads with distinct IDs |

### Test Structure

```
MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs (27 tests)
├── BuildKmerDatabase Tests (8 tests)
│   ├── M1:  EmptyInput_ReturnsEmptyDatabase
│   ├── M2:  SequenceShorterThanK_ReturnsEmpty
│   ├── M3:  ValidReference_ProducesKmers (exact count=6)
│   ├── M4:  UsesCanonicalKmers (A→A)
│   ├── M4b: CanonicalKmer_UsesReverseComplementWhenSmaller (T→A)
│   ├── M5:  NonAcgtCharacters_Excluded
│   ├── M6:  KmerCount_FollowsFormula (exact count=len-k+1)
│   └── S1:  MixedCase_NormalizedToUppercase
├── ClassifyReads Tests (16 tests)
│   ├── M7:  EmptySequence_ReturnsUnclassified
│   ├── M8:  SequenceShorterThanK_ReturnsUnclassified
│   ├── M9:  NoMatch_ReturnsUnclassified
│   ├── M10: MatchingKmers_ClassifiesCorrectly (exact: 1/8=0.125)
│   ├── M11: OutputCountEqualsInputCount
│   ├── M12: ConfidenceCalculation_IsCorrect (exact: 1/2=0.5)
│   ├── M13: TotalKmers_MatchesFormula (4 cases)
│   ├── M14: MatchedKmers_BoundedByTotal
│   ├── M15: MultipleReads_AllClassified (exact: Bacteria, Bacteria, Unclassified)
│   ├── M16: TaxonomyParsing_PipeDelimited (7 ranks)
│   ├── M16b:TaxonomyParsing_SemicolonDelimited
│   ├── S1:  MixedCaseInput_Handled (exact: Kingdom=Bacteria, Matched=2)
│   ├── S2:  MultipleTaxonMatches_ResolvesToHighestCount
│   ├── S3:  CanonicalKmerLookup_MatchesReverseComplement
│   ├── S4:  AmbiguousNucleotides_ExcludedFromTotalKmers
│   └── S5:  MultiTaxon_ConfidenceUsesWinningTaxonCount (C/Q=3/4)
└── Invariants Tests (3 tests)
    ├── AllOutputs_HaveValidConfidence [0,1]
    ├── UnclassifiedReads_HaveZeroMatchedKmers
    └── ReadIdPreserved
```

## Open Questions / Decisions

None. Algorithm behavior is well-documented.

## Test File

**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs`
