# META-CLASS-001: Taxonomic Classification Test Specification

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | META-CLASS-001 |
| **Area** | Metagenomics |
| **Canonical Methods** | `MetagenomicsAnalyzer.ClassifyReads`, `MetagenomicsAnalyzer.BuildKmerDatabase` |
| **Complexity** | O(n × m) where n=reads, m=read length |
| **Invariant** | |output| = |input reads|; 0 ≤ Confidence ≤ 1 |

## Methods Under Test

| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `ClassifyReads(reads, kmerDB, k)` | MetagenomicsAnalyzer | Canonical | Deep |
| `BuildKmerDatabase(refGenomes, k)` | MetagenomicsAnalyzer | DB construction | Deep |

## Evidence Sources

1. **Wikipedia — Metagenomics:** K-mer based binning, taxonomic classification principles
2. **Kraken Documentation (CCB JHU):** Canonical k-mer classification algorithm, k=31 default
3. **Wood & Salzberg (2014):** Kraken paper describing exact k-mer matching for classification

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
| M12 | Confidence = MatchedKmers / TotalKmers | Kraken confidence formula |
| M13 | TotalKmers = max(0, len - k + 1) | K-mer counting formula |
| M14 | MatchedKmers ≤ TotalKmers | Bound invariant |
| M15 | Multiple reads all classified | Batch processing |
| M16 | Taxonomy string parsed correctly | Implementation: pipe/semicolon delimited |

### SHOULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | Mixed case input handled (uppercased internally) | Robustness for real data |
| S2 | Multiple taxon matches resolves to highest count | LCA-like behavior |
| S3 | Null sequence treated as empty | Defensive programming |
| S4 | Large batch of reads processed correctly | Scalability |
| S5 | k=1 edge case works correctly | Minimum valid k |

### COULD Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | Performance with large database acceptable | Practical use case |
| C2 | Memory-efficient for large read sets | Streaming enumerable |

## Consolidation Plan

### Current Test Pool

| File | Tests | Status |
|------|-------|--------|
| MetagenomicsAnalyzerTests.cs | 9 tests for ClassifyReads/BuildKmerDatabase | Existing, needs consolidation |

### Existing Test Audit

| Existing Test | Coverage | Action |
|---------------|----------|--------|
| BuildKmerDatabase_CreatesDatabase | M3 (partial) | Strengthen assertions |
| BuildKmerDatabase_EmptyInput_ReturnsEmpty | M1 | Keep |
| BuildKmerDatabase_ShortSequence_IgnoresIt | M2 | Keep |
| BuildKmerDatabase_UsesCanonicalKmers | M4 | Strengthen |
| ClassifyReads_WithMatchingDatabase_ClassifiesCorrectly | M10 (partial) | Strengthen invariant checks |
| ClassifyReads_NoMatch_ReturnsUnclassified | M9 | Keep |
| ClassifyReads_EmptySequence_HandlesGracefully | M7 | Keep |
| ClassifyReads_ShortSequence_ReturnsUnclassified | M8 | Keep |
| ClassifyReads_MultipleReads_ClassifiesAll | M11, M15 (partial) | Strengthen |

### Consolidation Actions

1. **Create** `MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs` for META-CLASS-001
2. **Migrate** relevant tests from `MetagenomicsAnalyzerTests.cs`
3. ~~**Add** missing MUST tests (M5, M6, M12, M13, M14, M16)~~ ✅ Done
4. ~~**Strengthen** existing tests with invariant assertions~~ ✅ Done
5. ~~**Remove** tests from generic file once migrated~~ ✅ Done

### Test Structure

```
MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs
├── BuildKmerDatabase_Tests (region)
│   ├── EmptyInput_ReturnsEmpty
│   ├── ShortSequence_IgnoresIt
│   ├── ValidReference_ProducesKmers
│   ├── UsesCanonicalKmers
│   ├── ExcludesNonAcgtKmers
│   └── KmerCountMatchesFormula
├── ClassifyReads_Tests (region)
│   ├── EmptySequence_ReturnsUnclassified
│   ├── ShortSequence_ReturnsUnclassified
│   ├── NoMatch_ReturnsUnclassified
│   ├── MatchingKmers_ClassifiesCorrectly
│   ├── OutputCountEqualsInputCount
│   ├── ConfidenceCalculatedCorrectly
│   ├── TotalKmersMatchesFormula
│   ├── MatchedKmersBoundedByTotal
│   ├── MultipleReads_AllClassified
│   ├── TaxonomyParsedCorrectly
│   └── MixedCaseHandled
└── Invariants_Tests (region)
    └── PropertyBased_AllInvariantsHold
```

## Open Questions / Decisions

None. Algorithm behavior is well-documented.

## Test File

**Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_TaxonomicClassification_Tests.cs`
