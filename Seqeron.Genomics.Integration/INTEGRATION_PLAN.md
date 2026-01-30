# Seqeron.Genomics.Integration - Implementation Plan

## Overview

Integration library for Seqeron.Genomics, providing connections to external bioinformatics databases, web services, and structural format parsing.

**Goals:**
- Clean Architecture with clear layer separation
- Dependency Injection for all external dependencies
- Contract tests for external APIs
- Resilience patterns (retry, circuit breaker, timeout)
- Full async/await and CancellationToken support

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              Application Layer                               │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────────┐ │
│  │ BlastService    │  │ SequenceService │  │ StructureService            │ │
│  │ (orchestration) │  │ (orchestration) │  │ (orchestration)             │ │
│  └────────┬────────┘  └────────┬────────┘  └─────────────┬───────────────┘ │
└───────────┼────────────────────┼────────────────────────┼───────────────────┘
            │                    │                        │
┌───────────┼────────────────────┼────────────────────────┼───────────────────┐
│           │           Domain Layer                      │                    │
│  ┌────────▼────────┐  ┌────────▼────────┐  ┌───────────▼─────────────────┐ │
│  │ IBlastClient    │  │ ISequenceDb     │  │ IStructureRepository        │ │
│  │ IBlastResult    │  │ IEntrezClient   │  │ IStructureParser            │ │
│  └─────────────────┘  └─────────────────┘  └─────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
            │                    │                        │
┌───────────┼────────────────────┼────────────────────────┼───────────────────┐
│           │       Infrastructure Layer                  │                    │
│  ┌────────▼────────┐  ┌────────▼────────┐  ┌───────────▼─────────────────┐ │
│  │ NcbiBlastClient │  │ NcbiEntrezClient│  │ PdbFileParser               │ │
│  │ (HTTP impl)     │  │ UniProtClient   │  │ MmCifParser                 │ │
│  │                 │  │ EnsemblClient   │  │ RcsbPdbClient               │ │
│  └─────────────────┘  └─────────────────┘  └─────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
            │                    │                        │
┌───────────┴────────────────────┴────────────────────────┴───────────────────┐
│                           Cross-Cutting Concerns                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────────────────┐ │
│  │ Resilience      │  │ Caching         │  │ Logging                     │ │
│  │ (Polly)         │  │ (IMemoryCache)  │  │ (ILogger<T>)                │ │
│  └─────────────────┘  └─────────────────┘  └─────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Project Structure

```
Seqeron.Genomics.Integration/
├── Seqeron.Genomics.Integration.csproj
│
├── Abstractions/                          # Domain interfaces (no dependencies)
│   ├── Blast/
│   │   ├── IBlastClient.cs
│   │   ├── IBlastQuery.cs
│   │   └── IBlastResult.cs
│   ├── Sequences/
│   │   ├── ISequenceDatabase.cs
│   │   ├── IEntrezClient.cs
│   │   └── ISequenceFetcher.cs
│   ├── Structure/
│   │   ├── IStructureParser.cs
│   │   ├── IStructureRepository.cs
│   │   └── IStructureClient.cs
│   └── Common/
│       ├── IIntegrationClient.cs
│       ├── IRetryPolicy.cs
│       └── ICacheProvider.cs
│
├── Domain/                                # Domain models (POCOs, no logic)
│   ├── Blast/
│   │   ├── BlastProgram.cs               # enum: blastn, blastp, blastx, tblastn, tblastx
│   │   ├── BlastDatabase.cs              # enum: nr, nt, refseq_protein, etc.
│   │   ├── BlastQuery.cs                 # Query parameters
│   │   ├── BlastResult.cs                # Result container
│   │   ├── BlastHit.cs                   # Single hit
│   │   ├── BlastHsp.cs                   # High-scoring segment pair
│   │   └── BlastStatistics.cs            # E-value, bit score, etc.
│   ├── Sequences/
│   │   ├── EntrezRecord.cs               # GenBank/RefSeq record
│   │   ├── EntrezDatabase.cs             # enum: nucleotide, protein, gene, etc.
│   │   ├── UniProtEntry.cs               # UniProt protein entry
│   │   ├── EnsemblGene.cs                # Ensembl gene record
│   │   ├── EnsemblTranscript.cs          # Ensembl transcript
│   │   ├── SequenceFeature.cs            # Annotations (CDS, gene, etc.)
│   │   └── CrossReference.cs             # Database cross-references
│   ├── Structure/
│   │   ├── ProteinStructure.cs           # 3D structure container
│   │   ├── Chain.cs                      # Polypeptide chain
│   │   ├── Residue.cs                    # Amino acid residue
│   │   ├── Atom.cs                       # Atom with coordinates
│   │   ├── AtomType.cs                   # enum: CA, CB, N, C, O, etc.
│   │   ├── SecondaryStructureElement.cs  # Helix, sheet, coil
│   │   ├── SecondaryStructureType.cs     # enum: Helix, Sheet, Coil, Turn
│   │   ├── LigandInfo.cs                 # Bound ligands/cofactors
│   │   └── StructureMetadata.cs          # Resolution, method, date, etc.
│   └── Common/
│       ├── DatabaseIdentifier.cs         # Unified ID (DB + accession)
│       ├── Organism.cs                   # Taxonomy info
│       ├── Citation.cs                   # Publication reference
│       └── IntegrationError.cs           # Error details
│
├── Infrastructure/                        # External implementations
│   ├── Http/
│   │   ├── HttpClientFactory.cs          # Configured HttpClient creation
│   │   ├── RateLimiter.cs                # Per-service rate limiting
│   │   └── RetryHandler.cs               # Polly-based retry delegating handler
│   ├── Ncbi/
│   │   ├── NcbiBlastClient.cs            # BLAST REST API
│   │   ├── NcbiEntrezClient.cs           # Entrez E-utilities
│   │   ├── NcbiConfiguration.cs          # API keys, endpoints
│   │   └── Parsers/
│   │       ├── BlastXmlParser.cs         # BLAST XML output parser
│   │       ├── GenBankParser.cs          # GenBank flat file parser
│   │       └── FastaParser.cs            # FASTA response parser
│   ├── UniProt/
│   │   ├── UniProtClient.cs              # UniProt REST API
│   │   ├── UniProtConfiguration.cs
│   │   └── Parsers/
│   │       ├── UniProtXmlParser.cs
│   │       └── UniProtTxtParser.cs
│   ├── Ensembl/
│   │   ├── EnsemblClient.cs              # Ensembl REST API
│   │   ├── EnsemblConfiguration.cs
│   │   └── Parsers/
│   │       └── EnsemblJsonParser.cs
│   └── Structure/
│       ├── RcsbPdbClient.cs              # RCSB PDB REST API
│       ├── PdbConfiguration.cs
│       └── Parsers/
│           ├── PdbFileParser.cs          # Legacy PDB format
│           ├── MmCifParser.cs            # mmCIF format (modern)
│           └── StructureValidator.cs     # Validation utilities
│
├── Services/                              # Application services (orchestration)
│   ├── BlastService.cs                   # High-level BLAST operations
│   ├── SequenceFetchService.cs           # Fetch from multiple DBs
│   ├── StructureService.cs               # Structure retrieval & parsing
│   ├── CrossReferenceService.cs          # ID mapping between databases
│   └── BatchProcessingService.cs         # Bulk operations with throttling
│
├── Caching/                               # Caching infrastructure
│   ├── MemoryCacheProvider.cs
│   ├── FileCacheProvider.cs              # For large structures
│   ├── CacheKeyGenerator.cs
│   └── CacheOptions.cs
│
├── Configuration/                         # Configuration models
│   ├── IntegrationOptions.cs             # Root configuration
│   ├── NcbiOptions.cs
│   ├── UniProtOptions.cs
│   ├── EnsemblOptions.cs
│   ├── PdbOptions.cs
│   └── ResilienceOptions.cs              # Retry, timeout, circuit breaker
│
├── Extensions/                            # DI and fluent extensions
│   ├── ServiceCollectionExtensions.cs    # AddGenomicsIntegration()
│   ├── DnaSequenceExtensions.cs          # seq.BlastAsync(), seq.FetchHomologs()
│   ├── ProteinSequenceExtensions.cs      # prot.FetchStructure()
│   └── HttpClientBuilderExtensions.cs    # Resilience configuration
│
└── Exceptions/                            # Custom exceptions
    ├── IntegrationException.cs           # Base exception
    ├── RateLimitExceededException.cs
    ├── ServiceUnavailableException.cs
    ├── InvalidResponseException.cs
    └── ParseException.cs
```

---

## Test Project Structure

```
Seqeron.Genomics.Integration.Tests/
├── Seqeron.Genomics.Integration.Tests.csproj
│
├── Unit/                                  # Pure unit tests (no I/O)
│   ├── Domain/
│   │   ├── BlastResultTests.cs
│   │   ├── ProteinStructureTests.cs
│   │   └── EntrezRecordTests.cs
│   ├── Parsers/
│   │   ├── BlastXmlParserTests.cs
│   │   ├── GenBankParserTests.cs
│   │   ├── PdbFileParserTests.cs
│   │   ├── MmCifParserTests.cs
│   │   └── UniProtXmlParserTests.cs
│   ├── Services/
│   │   ├── BlastServiceTests.cs
│   │   ├── SequenceFetchServiceTests.cs
│   │   └── StructureServiceTests.cs
│   └── Caching/
│       └── CacheKeyGeneratorTests.cs
│
├── Contract/                              # API contract tests
│   ├── Ncbi/
│   │   ├── NcbiBlastContractTests.cs     # BLAST API response schema
│   │   ├── NcbiEntrezContractTests.cs    # Entrez response schema
│   │   └── TestData/
│   │       ├── blast_response_sample.xml
│   │       ├── genbank_sample.gb
│   │       └── entrez_esearch_sample.xml
│   ├── UniProt/
│   │   ├── UniProtContractTests.cs
│   │   └── TestData/
│   │       ├── uniprot_entry_sample.xml
│   │       └── uniprot_entry_sample.txt
│   ├── Ensembl/
│   │   ├── EnsemblContractTests.cs
│   │   └── TestData/
│   │       └── ensembl_gene_sample.json
│   └── Pdb/
│       ├── PdbContractTests.cs
│       └── TestData/
│           ├── 1crn.pdb                  # Crambin (small protein)
│           ├── 1crn.cif                  # Same in mmCIF
│           └── rcsb_search_sample.json
│
├── Integration/                           # Real API integration tests
│   ├── Ncbi/
│   │   ├── NcbiBlastIntegrationTests.cs  # [Category("Integration")]
│   │   └── NcbiEntrezIntegrationTests.cs
│   ├── UniProt/
│   │   └── UniProtIntegrationTests.cs
│   ├── Ensembl/
│   │   └── EnsemblIntegrationTests.cs
│   └── Pdb/
│       └── RcsbPdbIntegrationTests.cs
│
├── Resilience/                            # Fault tolerance tests
│   ├── RetryPolicyTests.cs
│   ├── CircuitBreakerTests.cs
│   ├── TimeoutTests.cs
│   └── RateLimiterTests.cs
│
├── Fixtures/                              # Shared test infrastructure
│   ├── MockHttpMessageHandler.cs
│   ├── TestDataLoader.cs
│   ├── IntegrationTestBase.cs
│   └── WireMockServerFixture.cs          # For contract tests
│
└── GlobalUsings.cs
```

---

## Phase 1: Foundation & Abstractions (Week 1)

### 1.1 Project Setup

**Tasks:**
- [ ] Create `Seqeron.Genomics.Integration.csproj`
- [ ] Create `Seqeron.Genomics.Integration.Tests.csproj`
- [ ] Configure NuGet dependencies
- [ ] Setup EditorConfig and analyzers

**Dependencies:**
```xml
<ItemGroup>
  <!-- HTTP & Resilience -->
  <PackageReference Include="Microsoft.Extensions.Http" Version="8.0.0" />
  <PackageReference Include="Microsoft.Extensions.Http.Polly" Version="8.0.0" />
  <PackageReference Include="Polly" Version="8.2.0" />

  <!-- Caching -->
  <PackageReference Include="Microsoft.Extensions.Caching.Memory" Version="8.0.0" />

  <!-- Configuration & DI -->
  <PackageReference Include="Microsoft.Extensions.Options" Version="8.0.0" />
  <PackageReference Include="Microsoft.Extensions.DependencyInjection.Abstractions" Version="8.0.0" />

  <!-- Logging -->
  <PackageReference Include="Microsoft.Extensions.Logging.Abstractions" Version="8.0.0" />

  <!-- Parsing -->
  <PackageReference Include="System.Text.Json" Version="8.0.0" />
</ItemGroup>
```

**Test Dependencies:**
```xml
<ItemGroup>
  <PackageReference Include="NUnit" Version="3.14.0" />
  <PackageReference Include="NUnit3TestAdapter" Version="4.5.0" />
  <PackageReference Include="Microsoft.NET.Test.Sdk" Version="17.8.0" />
  <PackageReference Include="Moq" Version="4.20.70" />
  <PackageReference Include="FluentAssertions" Version="6.12.0" />
  <PackageReference Include="WireMock.Net" Version="1.5.40" />
</ItemGroup>
```

**Files:**
```
Seqeron.Genomics.Integration.csproj
Seqeron.Genomics.Integration.Tests.csproj
```

**Tests (5):**
- [ ] ProjectCompilationTests.cs — project compiles without errors

---

### 1.2 Core Abstractions

**Tasks:**
- [ ] `IIntegrationClient` — base interface for all clients
- [ ] `Result<T>` — operation result (success/failure pattern)
- [ ] `IntegrationException` — exception hierarchy
- [ ] `DatabaseIdentifier` — unified identifier

**Files:**
```
Abstractions/Common/
├── IIntegrationClient.cs
├── Result.cs
├── IRetryPolicy.cs
└── ICacheProvider.cs

Domain/Common/
├── DatabaseIdentifier.cs
├── Organism.cs
├── Citation.cs
└── IntegrationError.cs

Exceptions/
├── IntegrationException.cs
├── RateLimitExceededException.cs
├── ServiceUnavailableException.cs
├── InvalidResponseException.cs
└── ParseException.cs
```

**Interfaces:**

```csharp
// IIntegrationClient.cs
public interface IIntegrationClient : IDisposable
{
    string ServiceName { get; }
    bool IsAvailable { get; }

    Task<bool> HealthCheckAsync(CancellationToken cancellationToken = default);
}

// Result.cs
public readonly struct Result<T>
{
    public bool IsSuccess { get; }
    public T? Value { get; }
    public IntegrationError? Error { get; }

    public static Result<T> Success(T value) => new(true, value, null);
    public static Result<T> Failure(IntegrationError error) => new(false, default, error);

    public TResult Match<TResult>(Func<T, TResult> onSuccess, Func<IntegrationError, TResult> onFailure);
}
```

**Tests (15):**
- [ ] ResultTests.cs (8 tests)
  - Success_ContainsValue
  - Failure_ContainsError
  - Match_OnSuccess_CallsSuccessFunc
  - Match_OnFailure_CallsFailureFunc
  - IsSuccess_WhenSuccess_ReturnsTrue
  - IsSuccess_WhenFailure_ReturnsFalse
  - Value_WhenFailure_ReturnsDefault
  - Error_WhenSuccess_ReturnsNull

- [ ] DatabaseIdentifierTests.cs (7 tests)
  - Constructor_ValidInput_CreatesIdentifier
  - Parse_ValidString_ParsesCorrectly
  - Parse_InvalidString_ThrowsArgumentException
  - Equals_SameIdentifier_ReturnsTrue
  - GetHashCode_SameIdentifier_ReturnsSameHash
  - ToString_ReturnsFormattedString
  - TryParse_InvalidInput_ReturnsFalse

---

## Phase 2: NCBI BLAST Integration (Week 2-3)

### 2.1 BLAST Domain Models

**Tasks:**
- [ ] `BlastProgram` — BLAST program enum
- [ ] `BlastDatabase` — database enum
- [ ] `BlastQuery` — query parameters
- [ ] `BlastResult` — result container
- [ ] `BlastHit` — single hit
- [ ] `BlastHsp` — High-scoring Segment Pair
- [ ] `BlastStatistics` — statistics

**Files:**
```
Domain/Blast/
├── BlastProgram.cs
├── BlastDatabase.cs
├── BlastQuery.cs
├── BlastResult.cs
├── BlastHit.cs
├── BlastHsp.cs
└── BlastStatistics.cs
```

**Models:**

```csharp
// BlastProgram.cs
public enum BlastProgram
{
    Blastn,      // nucleotide vs nucleotide
    Blastp,      // protein vs protein
    Blastx,      // translated nucleotide vs protein
    Tblastn,     // protein vs translated nucleotide
    Tblastx      // translated nucleotide vs translated nucleotide
}

// BlastQuery.cs
public sealed class BlastQuery
{
    public required string Sequence { get; init; }
    public required BlastProgram Program { get; init; }
    public required BlastDatabase Database { get; init; }
    public double ExpectThreshold { get; init; } = 10.0;
    public int MaxHits { get; init; } = 100;
    public int WordSize { get; init; } = 11;
    public string? EntrezQuery { get; init; }
    public GapCosts GapCosts { get; init; } = GapCosts.Default;
}

// BlastHsp.cs
public readonly record struct BlastHsp(
    int QueryStart,
    int QueryEnd,
    int HitStart,
    int HitEnd,
    string QuerySequence,
    string HitSequence,
    string Midline,
    double BitScore,
    double EValue,
    int Identity,
    int Positives,
    int Gaps,
    int AlignmentLength);
```

**Tests (20):**
- [ ] BlastQueryTests.cs (8 tests)
- [ ] BlastHitTests.cs (6 tests)
- [ ] BlastHspTests.cs (6 tests)

---

### 2.2 BLAST Client Interface & Implementation

**Tasks:**
- [ ] `IBlastClient` — interface
- [ ] `NcbiBlastClient` — implementation NCBI BLAST REST API
- [ ] `BlastXmlParser` — XML response parser
- [ ] Retry/timeout policies

**Files:**
```
Abstractions/Blast/
├── IBlastClient.cs
├── IBlastQuery.cs
└── IBlastResult.cs

Infrastructure/Ncbi/
├── NcbiBlastClient.cs
├── NcbiConfiguration.cs
└── Parsers/
    └── BlastXmlParser.cs

Configuration/
├── NcbiOptions.cs
└── ResilienceOptions.cs
```

**Interface:**

```csharp
// IBlastClient.cs
public interface IBlastClient : IIntegrationClient
{
    /// <summary>
    /// Submits a BLAST query and returns job ID for polling.
    /// </summary>
    Task<Result<string>> SubmitQueryAsync(
        BlastQuery query,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Checks the status of a submitted BLAST job.
    /// </summary>
    Task<Result<BlastJobStatus>> GetJobStatusAsync(
        string jobId,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Retrieves results when job is complete.
    /// </summary>
    Task<Result<BlastResult>> GetResultsAsync(
        string jobId,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Convenience method: submit, poll, and return results.
    /// </summary>
    Task<Result<BlastResult>> RunBlastAsync(
        BlastQuery query,
        TimeSpan? timeout = null,
        IProgress<BlastJobStatus>? progress = null,
        CancellationToken cancellationToken = default);
}

public enum BlastJobStatus
{
    Queued,
    Running,
    Completed,
    Failed,
    Unknown
}
```

**Tests (35):**

*Unit Tests:*
- [ ] BlastXmlParserTests.cs (15 tests)
  - Parse_ValidXml_ReturnsBlastResult
  - Parse_MultipleHits_ParsesAllHits
  - Parse_NoHits_ReturnsEmptyResult
  - Parse_WithHsps_ParsesAllHsps
  - Parse_Statistics_ParsesCorrectly
  - Parse_InvalidXml_ThrowsParseException
  - Parse_MissingRequiredField_ThrowsParseException
  - Parse_LargeResult_HandlesCorrectly
  - Parse_SpecialCharacters_HandlesCorrectly
  - Parse_Blastn_ParsesCorrectly
  - Parse_Blastp_ParsesCorrectly
  - Parse_Blastx_ParsesCorrectly
  - Parse_WithEntrezLinks_ParsesCorrectly
  - Parse_RealWorldSample_Insulin_ParsesCorrectly
  - Parse_RealWorldSample_P53_ParsesCorrectly

*Contract Tests:*
- [ ] NcbiBlastContractTests.cs (10 tests)
  - ResponseSchema_MatchesExpectedStructure
  - HitFields_AllPresent
  - HspFields_AllPresent
  - StatisticsFields_AllPresent
  - ErrorResponse_HasExpectedFormat
  - JobStatusResponse_HasExpectedFormat
  - MultipleIterations_ParseCorrectly
  - TaxonomyInfo_ParsesCorrectly
  - AccessionFormat_IsValid
  - EValueFormat_IsParseable

*Integration Tests:*
- [ ] NcbiBlastIntegrationTests.cs (10 tests) `[Category("Integration")]`
  - SubmitQuery_ValidSequence_ReturnsJobId
  - GetJobStatus_ValidJobId_ReturnsStatus
  - RunBlast_ShortSequence_ReturnsResults
  - RunBlast_InvalidSequence_ReturnsError
  - RunBlast_Timeout_ReturnsTimeoutError
  - RunBlast_RateLimited_RetriesCorrectly
  - RunBlast_WithEntrezQuery_FiltersResults
  - RunBlast_Blastn_ReturnsNucleotideHits
  - RunBlast_Blastp_ReturnsProteinHits
  - HealthCheck_ReturnsTrue

---

### 2.3 BLAST Service (Orchestration)

**Tasks:**
- [ ] `BlastService` — high-level service
- [ ] Integration with `DnaSequence` and `ProteinSequence`
- [ ] Extension methods for fluent API

**Files:**
```
Services/
└── BlastService.cs

Extensions/
├── DnaSequenceExtensions.cs
└── ProteinSequenceExtensions.cs
```

**Service:**

```csharp
// BlastService.cs
public class BlastService
{
    private readonly IBlastClient _blastClient;
    private readonly ILogger<BlastService> _logger;
    private readonly ICacheProvider _cache;

    public async Task<Result<BlastResult>> SearchHomologsAsync(
        ISequence sequence,
        BlastSearchOptions? options = null,
        CancellationToken cancellationToken = default);

    public async Task<Result<IReadOnlyList<BlastHit>>> FindSimilarSequencesAsync(
        ISequence sequence,
        double minIdentity = 0.8,
        double maxEValue = 1e-10,
        CancellationToken cancellationToken = default);
}

// DnaSequenceExtensions.cs
public static class DnaSequenceExtensions
{
    public static Task<Result<BlastResult>> BlastAsync(
        this DnaSequence sequence,
        IBlastClient client,
        BlastDatabase database = BlastDatabase.Nt,
        CancellationToken cancellationToken = default);
}
```

**Tests (15):**
- [ ] BlastServiceTests.cs (15 tests)
  - SearchHomologs_ValidDnaSequence_ReturnsResults
  - SearchHomologs_ValidProteinSequence_ReturnsResults
  - SearchHomologs_CachesResults
  - SearchHomologs_CancellationRequested_ThrowsOperationCanceled
  - FindSimilarSequences_FiltersCorrectly
  - FindSimilarSequences_NoMatches_ReturnsEmpty
  - SearchHomologs_ClientError_ReturnsFailure
  - SearchHomologs_WithProgress_ReportsProgress
  - Integration_DnaSequence_BlastAsync_Works
  - Integration_ProteinSequence_BlastAsync_Works

---

## Phase 3: NCBI Entrez Integration (Week 4-5)

### 3.1 Entrez Domain Models

**Tasks:**
- [ ] `EntrezDatabase` — database enum
- [ ] `EntrezRecord` — GenBank/RefSeq record
- [ ] `SequenceFeature` — annotations (CDS, gene, etc.)
- [ ] `CrossReference` — cross-references

**Files:**
```
Domain/Sequences/
├── EntrezRecord.cs
├── EntrezDatabase.cs
├── SequenceFeature.cs
├── FeatureType.cs
├── FeatureQualifier.cs
└── CrossReference.cs
```

**Models:**

```csharp
// EntrezDatabase.cs
public enum EntrezDatabase
{
    Nucleotide,   // nuccore
    Protein,      // protein
    Gene,         // gene
    PubMed,       // pubmed
    Taxonomy,     // taxonomy
    Structure,    // structure
    Snp,          // snp
    Biosample,    // biosample
    Sra           // sra
}

// EntrezRecord.cs
public sealed class EntrezRecord
{
    public required string Accession { get; init; }
    public required string Version { get; init; }
    public required EntrezDatabase Database { get; init; }
    public required string Definition { get; init; }
    public required string Sequence { get; init; }
    public required int Length { get; init; }
    public required Organism Organism { get; init; }
    public required DateOnly CreateDate { get; init; }
    public required DateOnly UpdateDate { get; init; }
    public IReadOnlyList<SequenceFeature> Features { get; init; } = [];
    public IReadOnlyList<CrossReference> CrossReferences { get; init; } = [];
    public IReadOnlyList<Citation> Citations { get; init; } = [];
    public string? Comment { get; init; }
    public string? Keywords { get; init; }
}

// SequenceFeature.cs
public sealed class SequenceFeature
{
    public required FeatureType Type { get; init; }
    public required int Start { get; init; }
    public required int End { get; init; }
    public required bool IsComplement { get; init; }
    public IReadOnlyDictionary<string, string> Qualifiers { get; init; } =
        ImmutableDictionary<string, string>.Empty;

    // Convenience properties
    public string? Gene => Qualifiers.GetValueOrDefault("gene");
    public string? Product => Qualifiers.GetValueOrDefault("product");
    public string? ProteinId => Qualifiers.GetValueOrDefault("protein_id");
    public string? Translation => Qualifiers.GetValueOrDefault("translation");
}
```

**Tests (20):**
- [ ] EntrezRecordTests.cs (10 tests)
- [ ] SequenceFeatureTests.cs (10 tests)

---

### 3.2 Entrez Client Interface & Implementation

**Tasks:**
- [ ] `IEntrezClient` — interface E-utilities
- [ ] `NcbiEntrezClient` — implementation
- [ ] `GenBankParser` — GenBank format parser
- [ ] Rate limiting (NCBI requires max 3 req/sec without API key, 10/sec with key)

**Files:**
```
Abstractions/Sequences/
├── IEntrezClient.cs
├── ISequenceDatabase.cs
└── ISequenceFetcher.cs

Infrastructure/Ncbi/
├── NcbiEntrezClient.cs
└── Parsers/
    ├── GenBankParser.cs
    └── FastaParser.cs
```

**Interface:**

```csharp
// IEntrezClient.cs
public interface IEntrezClient : IIntegrationClient
{
    /// <summary>
    /// Searches Entrez database and returns matching IDs.
    /// </summary>
    Task<Result<EntrezSearchResult>> SearchAsync(
        EntrezDatabase database,
        string query,
        int maxResults = 20,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Fetches full records by IDs.
    /// </summary>
    Task<Result<IReadOnlyList<EntrezRecord>>> FetchAsync(
        EntrezDatabase database,
        IEnumerable<string> ids,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Fetches a single record by accession.
    /// </summary>
    Task<Result<EntrezRecord>> FetchByAccessionAsync(
        string accession,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Fetches sequence only (FASTA format).
    /// </summary>
    Task<Result<string>> FetchSequenceAsync(
        string accession,
        CancellationToken cancellationToken = default);

    /// <summary>
    /// Gets linked records from another database.
    /// </summary>
    Task<Result<IReadOnlyList<string>>> GetLinksAsync(
        EntrezDatabase fromDb,
        EntrezDatabase toDb,
        IEnumerable<string> ids,
        CancellationToken cancellationToken = default);
}

public sealed class EntrezSearchResult
{
    public required int TotalCount { get; init; }
    public required IReadOnlyList<string> Ids { get; init; }
    public string? QueryTranslation { get; init; }
}
```

**Tests (40):**

*Unit Tests:*
- [ ] GenBankParserTests.cs (20 tests)
  - Parse_ValidGenBank_ReturnsRecord
  - Parse_WithFeatures_ParsesAllFeatures
  - Parse_CdsFeature_ExtractsTranslation
  - Parse_GeneFeature_ExtractsGeneName
  - Parse_ComplementFeature_SetsFlag
  - Parse_JoinLocation_HandlesCorrectly
  - Parse_MultipleRecords_ParsesAll
  - Parse_WithReferences_ParsesCitations
  - Parse_WithDbxref_ParsesCrossReferences
  - Parse_Organism_ParsesTaxonomy
  - Parse_InvalidFormat_ThrowsParseException
  - Parse_MissingLocus_ThrowsParseException
  - Parse_EmptySequence_ReturnsEmptyString
  - Parse_LongSequence_HandlesCorrectly
  - Parse_CircularMolecule_SetsFlag
  - Parse_RealSample_HumanInsulin_ParsesCorrectly
  - Parse_RealSample_EcoliGenome_ParsesCorrectly
  - Parse_RealSample_Plasmid_ParsesCorrectly
  - Parse_QualifiersWithQuotes_HandlesCorrectly
  - Parse_MultiLineQualifier_HandlesCorrectly

*Contract Tests:*
- [ ] NcbiEntrezContractTests.cs (10 tests)
  - ESearchResponse_MatchesSchema
  - EFetchGenBank_MatchesSchema
  - EFetchFasta_MatchesSchema
  - ELinkResponse_MatchesSchema
  - ErrorResponse_MatchesSchema
  - RateLimitResponse_MatchesSchema
  - AccessionFormats_AreValid
  - DateFormats_AreParseable
  - FeatureLocations_AreValid
  - CrossReferences_AreValid

*Integration Tests:*
- [ ] NcbiEntrezIntegrationTests.cs (10 tests) `[Category("Integration")]`
  - Search_ValidQuery_ReturnsIds
  - FetchByAccession_ValidAccession_ReturnsRecord
  - FetchSequence_ValidAccession_ReturnsFasta
  - GetLinks_NucleotideToProtein_ReturnsLinks
  - Fetch_MultipleIds_ReturnsAllRecords
  - Search_NoResults_ReturnsEmpty
  - Fetch_InvalidAccession_ReturnsError
  - RateLimit_ManyRequests_RespectsLimit
  - HealthCheck_ReturnsTrue
  - FetchWithApiKey_Works

---

### 3.3 Sequence Fetch Service

**Tasks:**
- [ ] `SequenceFetchService` — unified service
- [ ] Convert `EntrezRecord` → `DnaSequence`/`ProteinSequence`
- [ ] Batch operations with throttling

**Files:**
```
Services/
├── SequenceFetchService.cs
└── BatchProcessingService.cs
```

**Tests (15):**
- [ ] SequenceFetchServiceTests.cs (15 tests)

---

## Phase 4: UniProt Integration (Week 6)

### 4.1 UniProt Domain Models

**Tasks:**
- [ ] `UniProtEntry` — complete UniProt entry
- [ ] `UniProtFeature` — annotations (domains, sites, etc.)
- [ ] `GoAnnotation` — Gene Ontology terms
- [ ] `ProteinEvidence` — evidence level

**Files:**
```
Domain/Sequences/
├── UniProtEntry.cs
├── UniProtFeature.cs
├── UniProtFeatureType.cs
├── GoAnnotation.cs
├── GoCategory.cs
└── ProteinEvidence.cs
```

**Models:**

```csharp
// UniProtEntry.cs
public sealed class UniProtEntry
{
    public required string AccessionPrimary { get; init; }
    public IReadOnlyList<string> AccessionSecondary { get; init; } = [];
    public required string EntryName { get; init; }
    public required string ProteinName { get; init; }
    public required string GeneName { get; init; }
    public required Organism Organism { get; init; }
    public required string Sequence { get; init; }
    public required int Length { get; init; }
    public required ProteinEvidence Evidence { get; init; }
    public required DateOnly CreateDate { get; init; }
    public required DateOnly ModifyDate { get; init; }
    public IReadOnlyList<UniProtFeature> Features { get; init; } = [];
    public IReadOnlyList<GoAnnotation> GoAnnotations { get; init; } = [];
    public IReadOnlyList<CrossReference> CrossReferences { get; init; } = [];
    public IReadOnlyList<Citation> Citations { get; init; } = [];
    public string? Function { get; init; }
    public string? SubcellularLocation { get; init; }
    public IReadOnlyList<string> Keywords { get; init; } = [];
}

// UniProtFeature.cs
public sealed class UniProtFeature
{
    public required UniProtFeatureType Type { get; init; }
    public required int Start { get; init; }
    public required int End { get; init; }
    public string? Description { get; init; }
    public string? Evidence { get; init; }
    public string? Id { get; init; }
}

public enum UniProtFeatureType
{
    Chain,
    Domain,
    Repeat,
    Region,
    Motif,
    Site,
    ActiveSite,
    BindingSite,
    DnaBind,
    ZnFinger,
    ModifiedResidue,
    Lipidation,
    Glycosylation,
    DisulfideBond,
    CrossLink,
    Helix,
    Strand,
    Turn
}
```

**Tests (15):**
- [ ] UniProtEntryTests.cs (8 tests)
- [ ] UniProtFeatureTests.cs (7 tests)

---

### 4.2 UniProt Client

**Tasks:**
- [ ] `IUniProtClient` — interface
- [ ] `UniProtClient` — implementation REST API
- [ ] `UniProtXmlParser` — XML parser
- [ ] `UniProtTxtParser` — text format parser

**Files:**
```
Abstractions/Sequences/
└── IUniProtClient.cs

Infrastructure/UniProt/
├── UniProtClient.cs
├── UniProtConfiguration.cs
└── Parsers/
    ├── UniProtXmlParser.cs
    └── UniProtTxtParser.cs

Configuration/
└── UniProtOptions.cs
```

**Interface:**

```csharp
// IUniProtClient.cs
public interface IUniProtClient : IIntegrationClient
{
    Task<Result<UniProtEntry>> GetByAccessionAsync(
        string accession,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<UniProtEntry>>> SearchAsync(
        string query,
        int maxResults = 25,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<UniProtEntry>>> GetByGeneAsync(
        string geneName,
        string? organism = null,
        CancellationToken cancellationToken = default);

    Task<Result<string>> GetSequenceAsync(
        string accession,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<string>>> MapIdsAsync(
        IEnumerable<string> ids,
        string fromDb,
        string toDb,
        CancellationToken cancellationToken = default);
}
```

**Tests (30):**

*Unit Tests:*
- [ ] UniProtXmlParserTests.cs (12 tests)
- [ ] UniProtTxtParserTests.cs (8 tests)

*Contract Tests:*
- [ ] UniProtContractTests.cs (5 tests)

*Integration Tests:*
- [ ] UniProtIntegrationTests.cs (5 tests) `[Category("Integration")]`

---

## Phase 5: Ensembl Integration (Week 7)

### 5.1 Ensembl Domain Models

**Tasks:**
- [ ] `EnsemblGene` — gene
- [ ] `EnsemblTranscript` — transcript
- [ ] `EnsemblExon` — exon
- [ ] `EnsemblVariant` — variant

**Files:**
```
Domain/Sequences/
├── EnsemblGene.cs
├── EnsemblTranscript.cs
├── EnsemblExon.cs
├── EnsemblVariant.cs
├── Biotype.cs
└── Strand.cs
```

**Tests (15):**
- [ ] EnsemblGeneTests.cs (8 tests)
- [ ] EnsemblTranscriptTests.cs (7 tests)

---

### 5.2 Ensembl Client

**Tasks:**
- [ ] `IEnsemblClient` — interface
- [ ] `EnsemblClient` — implementation REST API
- [ ] `EnsemblJsonParser` — JSON parser

**Files:**
```
Abstractions/Sequences/
└── IEnsemblClient.cs

Infrastructure/Ensembl/
├── EnsemblClient.cs
├── EnsemblConfiguration.cs
└── Parsers/
    └── EnsemblJsonParser.cs

Configuration/
└── EnsemblOptions.cs
```

**Interface:**

```csharp
// IEnsemblClient.cs
public interface IEnsemblClient : IIntegrationClient
{
    Task<Result<EnsemblGene>> GetGeneByIdAsync(
        string ensemblId,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<EnsemblGene>>> SearchGenesAsync(
        string query,
        string species = "homo_sapiens",
        CancellationToken cancellationToken = default);

    Task<Result<EnsemblTranscript>> GetTranscriptByIdAsync(
        string ensemblId,
        CancellationToken cancellationToken = default);

    Task<Result<string>> GetSequenceAsync(
        string ensemblId,
        SequenceType type = SequenceType.Genomic,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<EnsemblVariant>>> GetVariantsAsync(
        string region,
        string species = "homo_sapiens",
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<CrossReference>>> GetCrossReferencesAsync(
        string ensemblId,
        CancellationToken cancellationToken = default);
}
```

**Tests (25):**

*Unit Tests:*
- [ ] EnsemblJsonParserTests.cs (10 tests)

*Contract Tests:*
- [ ] EnsemblContractTests.cs (8 tests)

*Integration Tests:*
- [ ] EnsemblIntegrationTests.cs (7 tests) `[Category("Integration")]`

---

## Phase 6: PDB Structure Integration (Week 8-9)

### 6.1 Structure Domain Models

**Tasks:**
- [ ] `ProteinStructure` — 3D structure
- [ ] `Chain` — chain
- [ ] `Residue` — residue
- [ ] `Atom` — atom with coordinates
- [ ] `SecondaryStructureElement` — secondary structure

**Files:**
```
Domain/Structure/
├── ProteinStructure.cs
├── Chain.cs
├── Residue.cs
├── Atom.cs
├── AtomType.cs
├── SecondaryStructureElement.cs
├── SecondaryStructureType.cs
├── LigandInfo.cs
├── StructureMetadata.cs
├── ExperimentalMethod.cs
└── Coordinates3D.cs
```

**Models:**

```csharp
// ProteinStructure.cs
public sealed class ProteinStructure
{
    public required string PdbId { get; init; }
    public required StructureMetadata Metadata { get; init; }
    public IReadOnlyList<Chain> Chains { get; init; } = [];
    public IReadOnlyList<LigandInfo> Ligands { get; init; } = [];
    public IReadOnlyList<SecondaryStructureElement> SecondaryStructure { get; init; } = [];

    // Computed properties
    public int TotalAtoms => Chains.Sum(c => c.Atoms.Count);
    public int TotalResidues => Chains.Sum(c => c.Residues.Count);
    public IEnumerable<Atom> AllAtoms => Chains.SelectMany(c => c.Atoms);
}

// Chain.cs
public sealed class Chain
{
    public required char ChainId { get; init; }
    public required string EntityType { get; init; }  // polymer, non-polymer, water
    public IReadOnlyList<Residue> Residues { get; init; } = [];
    public IReadOnlyList<Atom> Atoms { get; init; } = [];
    public string? Sequence { get; init; }
    public Organism? Organism { get; init; }
}

// Atom.cs
public readonly record struct Atom(
    int Serial,
    AtomType Type,
    string Name,
    char AltLoc,
    string ResidueName,
    char ChainId,
    int ResidueSeq,
    Coordinates3D Coordinates,
    double Occupancy,
    double BFactor,
    string Element);

// Coordinates3D.cs
public readonly record struct Coordinates3D(double X, double Y, double Z)
{
    public double DistanceTo(Coordinates3D other) =>
        Math.Sqrt(
            Math.Pow(X - other.X, 2) +
            Math.Pow(Y - other.Y, 2) +
            Math.Pow(Z - other.Z, 2));
}

// StructureMetadata.cs
public sealed class StructureMetadata
{
    public required string Title { get; init; }
    public required ExperimentalMethod Method { get; init; }
    public required double? Resolution { get; init; }  // Angstroms
    public required DateOnly ReleaseDate { get; init; }
    public IReadOnlyList<string> Authors { get; init; } = [];
    public IReadOnlyList<Citation> Citations { get; init; } = [];
    public string? Keywords { get; init; }
}

public enum ExperimentalMethod
{
    XRay,
    NMR,
    CryoEM,
    NeutronDiffraction,
    ElectronDiffraction,
    Other
}
```

**Tests (25):**
- [ ] ProteinStructureTests.cs (10 tests)
- [ ] ChainTests.cs (5 tests)
- [ ] AtomTests.cs (5 tests)
- [ ] Coordinates3DTests.cs (5 tests)

---

### 6.2 PDB File Parser

**Tasks:**
- [ ] `PdbFileParser` — legacy PDB format parser
- [ ] `MmCifParser` — mmCIF format parser
- [ ] `StructureValidator` — validation

**Files:**
```
Abstractions/Structure/
├── IStructureParser.cs
└── IStructureRepository.cs

Infrastructure/Structure/
└── Parsers/
    ├── PdbFileParser.cs
    ├── MmCifParser.cs
    └── StructureValidator.cs
```

**Parser Interface:**

```csharp
// IStructureParser.cs
public interface IStructureParser
{
    bool CanParse(string filePath);
    bool CanParse(Stream stream, string? hint = null);

    Result<ProteinStructure> Parse(string filePath);
    Result<ProteinStructure> Parse(Stream stream);
    Result<ProteinStructure> Parse(TextReader reader);

    Task<Result<ProteinStructure>> ParseAsync(
        Stream stream,
        CancellationToken cancellationToken = default);
}
```

**Tests (50):**

*Unit Tests:*
- [ ] PdbFileParserTests.cs (25 tests)
  - Parse_ValidPdb_ReturnsStructure
  - Parse_AtomRecords_ParsesCorrectly
  - Parse_HetatmRecords_ParsesCorrectly
  - Parse_MultipleChains_ParsesAll
  - Parse_SecondaryStructure_ParsesHelices
  - Parse_SecondaryStructure_ParsesSheets
  - Parse_Header_ParsesMetadata
  - Parse_Resolution_ParsesCorrectly
  - Parse_Model_ParsesFirstModel
  - Parse_AltLoc_HandlesCorrectly
  - Parse_MissingAtoms_HandlesGracefully
  - Parse_NegativeResidueNumbers_HandlesCorrectly
  - Parse_InsertionCodes_HandlesCorrectly
  - Parse_Anisou_Ignores
  - Parse_Conect_IgnoresForNow
  - Parse_RealStructure_1crn_ParsesCorrectly
  - Parse_RealStructure_1ubq_ParsesCorrectly
  - Parse_LargeStructure_Ribosome_HandlesCorrectly
  - Parse_NmrStructure_MultipleModels_ParsesFirst
  - Parse_InvalidFormat_ThrowsParseException
  - Parse_EmptyFile_ThrowsParseException
  - Parse_CorruptedFile_ThrowsParseException
  - Parse_Stream_WorksCorrectly
  - Parse_TextReader_WorksCorrectly
  - CanParse_PdbExtension_ReturnsTrue

- [ ] MmCifParserTests.cs (20 tests)
  - Parse_ValidMmCif_ReturnsStructure
  - Parse_AtomSite_ParsesCoordinates
  - Parse_EntityPoly_ParsesSequence
  - Parse_StructConf_ParsesHelices
  - Parse_StructSheet_ParsesSheets
  - Parse_Exptl_ParsesMethod
  - Parse_Refine_ParsesResolution
  - Parse_Citation_ParsesReferences
  - Parse_MultipleDataBlocks_ParsesFirst
  - Parse_CategoryLoop_ParsesCorrectly
  - Parse_QuotedStrings_HandlesCorrectly
  - Parse_Semicolons_HandlesMultiline
  - Parse_RealStructure_1crn_cif_ParsesCorrectly
  - Parse_RealStructure_Large_HandlesCorrectly
  - Parse_Stream_WorksCorrectly
  - Parse_InvalidFormat_ThrowsParseException
  - CanParse_CifExtension_ReturnsTrue
  - CanParse_MmcifExtension_ReturnsTrue
  - Parse_Async_WorksCorrectly
  - Parse_WithCancellation_Cancels

*Contract Tests:*
- [ ] PdbContractTests.cs (5 tests)
  - PdbFormat_AtomRecord_Has80Columns
  - PdbFormat_ChainId_Column22
  - PdbFormat_Coordinates_Columns31to54
  - MmCifFormat_ValidCategories
  - MmCifFormat_ValidDataTypes

---

### 6.3 RCSB PDB Client

**Tasks:**
- [ ] `IRcsbPdbClient` — interface
- [ ] `RcsbPdbClient` — implementation REST API
- [ ] Search structures by sequence

**Files:**
```
Abstractions/Structure/
└── IStructureClient.cs

Infrastructure/Structure/
├── RcsbPdbClient.cs
└── PdbConfiguration.cs

Configuration/
└── PdbOptions.cs
```

**Interface:**

```csharp
// IStructureClient.cs
public interface IStructureClient : IIntegrationClient
{
    Task<Result<ProteinStructure>> GetByIdAsync(
        string pdbId,
        CancellationToken cancellationToken = default);

    Task<Result<Stream>> DownloadAsync(
        string pdbId,
        StructureFormat format = StructureFormat.MmCif,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<PdbSearchHit>>> SearchBySequenceAsync(
        string sequence,
        double identityThreshold = 0.9,
        int maxResults = 10,
        CancellationToken cancellationToken = default);

    Task<Result<IReadOnlyList<PdbSearchHit>>> SearchByTextAsync(
        string query,
        int maxResults = 10,
        CancellationToken cancellationToken = default);

    Task<Result<StructureMetadata>> GetMetadataAsync(
        string pdbId,
        CancellationToken cancellationToken = default);
}

public enum StructureFormat
{
    Pdb,
    MmCif,
    Xml
}

public sealed class PdbSearchHit
{
    public required string PdbId { get; init; }
    public required double Score { get; init; }
    public string? Title { get; init; }
    public ExperimentalMethod? Method { get; init; }
    public double? Resolution { get; init; }
}
```

**Tests (20):**

*Unit Tests:*
- [ ] RcsbPdbClientTests.cs (5 tests)

*Contract Tests:*
- [ ] RcsbPdbContractTests.cs (8 tests)

*Integration Tests:*
- [ ] RcsbPdbIntegrationTests.cs (7 tests) `[Category("Integration")]`
  - GetById_ValidId_ReturnsStructure
  - Download_MmCif_ReturnsStream
  - SearchBySequence_ValidSequence_ReturnsHits
  - SearchByText_ValidQuery_ReturnsHits
  - GetMetadata_ValidId_ReturnsMetadata
  - GetById_InvalidId_ReturnsError
  - HealthCheck_ReturnsTrue

---

### 6.4 Structure Service

**Tasks:**
- [ ] `StructureService` — high-level service
- [ ] Integration with `ProteinSequence`
- [ ] Structure caching
- [ ] Extension methods

**Files:**
```
Services/
└── StructureService.cs

Extensions/
└── ProteinSequenceExtensions.cs (extend)
```

**Service:**

```csharp
// StructureService.cs
public class StructureService
{
    public async Task<Result<ProteinStructure>> GetStructureAsync(
        string pdbId,
        CancellationToken cancellationToken = default);

    public async Task<Result<IReadOnlyList<ProteinStructure>>> FindStructuresForSequenceAsync(
        ProteinSequence sequence,
        double minIdentity = 0.9,
        int maxResults = 5,
        CancellationToken cancellationToken = default);

    public async Task<Result<ProteinStructure>> ParseLocalFileAsync(
        string filePath,
        CancellationToken cancellationToken = default);

    public Result<ProteinStructure> ParseLocalFile(string filePath);

    // Analysis helpers
    public IEnumerable<Atom> GetBackboneAtoms(ProteinStructure structure);
    public double CalculateRmsd(ProteinStructure s1, ProteinStructure s2);
}

// Extension
public static class ProteinSequenceExtensions
{
    public static Task<Result<IReadOnlyList<ProteinStructure>>> FindStructuresAsync(
        this ProteinSequence sequence,
        IStructureClient client,
        double minIdentity = 0.9,
        CancellationToken cancellationToken = default);
}
```

**Tests (15):**
- [ ] StructureServiceTests.cs (15 tests)

---

## Phase 7: Cross-Cutting Concerns (Week 10)

### 7.1 Resilience Infrastructure

**Tasks:**
- [ ] `RetryHandler` — Polly retry policies
- [ ] `CircuitBreakerHandler` — circuit breaker
- [ ] `RateLimiter` — rate limiting per service
- [ ] `TimeoutHandler` — timeout policies

**Files:**
```
Infrastructure/Http/
├── HttpClientFactory.cs
├── RateLimiter.cs
├── RetryHandler.cs
├── CircuitBreakerHandler.cs
└── TimeoutHandler.cs

Configuration/
└── ResilienceOptions.cs
```

**Configuration:**

```csharp
// ResilienceOptions.cs
public sealed class ResilienceOptions
{
    public int MaxRetries { get; set; } = 3;
    public TimeSpan InitialRetryDelay { get; set; } = TimeSpan.FromSeconds(1);
    public TimeSpan MaxRetryDelay { get; set; } = TimeSpan.FromSeconds(30);
    public TimeSpan Timeout { get; set; } = TimeSpan.FromMinutes(2);
    public int CircuitBreakerThreshold { get; set; } = 5;
    public TimeSpan CircuitBreakerDuration { get; set; } = TimeSpan.FromMinutes(1);
}
```

**Tests (30):**
- [ ] RetryPolicyTests.cs (10 tests)
  - Retry_TransientError_RetriesCorrectly
  - Retry_MaxRetriesExceeded_ThrowsException
  - Retry_ExponentialBackoff_IncreaseDelay
  - Retry_NonTransientError_DoesNotRetry
  - Retry_CancellationRequested_Cancels

- [ ] CircuitBreakerTests.cs (10 tests)
  - CircuitBreaker_ThresholdExceeded_Opens
  - CircuitBreaker_Open_RejectsRequests
  - CircuitBreaker_HalfOpen_AllowsProbe
  - CircuitBreaker_Recovered_Closes
  - CircuitBreaker_Timeout_CountsAsFailure

- [ ] RateLimiterTests.cs (10 tests)
  - RateLimiter_UnderLimit_Allows
  - RateLimiter_OverLimit_Delays
  - RateLimiter_PerService_Isolated
  - RateLimiter_Concurrent_HandlesCorrectly
  - RateLimiter_Reset_ClearsState

---

### 7.2 Caching Infrastructure

**Tasks:**
- [ ] `MemoryCacheProvider` — in-memory cache
- [ ] `FileCacheProvider` — file cache for structures
- [ ] `CacheKeyGenerator` — key generation
- [ ] Cache invalidation policies

**Files:**
```
Caching/
├── ICacheProvider.cs
├── MemoryCacheProvider.cs
├── FileCacheProvider.cs
├── CacheKeyGenerator.cs
└── CacheOptions.cs
```

**Tests (15):**
- [ ] MemoryCacheProviderTests.cs (8 tests)
- [ ] FileCacheProviderTests.cs (7 tests)

---

### 7.3 Dependency Injection Setup

**Tasks:**
- [ ] `ServiceCollectionExtensions` — service registration
- [ ] Configuration binding
- [ ] Health checks registration

**Files:**
```
Extensions/
└── ServiceCollectionExtensions.cs
```

**Usage:**

```csharp
// ServiceCollectionExtensions.cs
public static class ServiceCollectionExtensions
{
    public static IServiceCollection AddGenomicsIntegration(
        this IServiceCollection services,
        Action<IntegrationOptions>? configure = null)
    {
        var options = new IntegrationOptions();
        configure?.Invoke(options);

        // Register options
        services.Configure<IntegrationOptions>(o => { /* copy */ });

        // Register HTTP clients with Polly
        services.AddHttpClient<IBlastClient, NcbiBlastClient>()
            .AddPolicyHandler(GetRetryPolicy())
            .AddPolicyHandler(GetCircuitBreakerPolicy());

        services.AddHttpClient<IEntrezClient, NcbiEntrezClient>()
            .AddPolicyHandler(GetRetryPolicy());

        // ... other clients

        // Register services
        services.AddScoped<BlastService>();
        services.AddScoped<SequenceFetchService>();
        services.AddScoped<StructureService>();

        // Register caching
        services.AddMemoryCache();
        services.AddSingleton<ICacheProvider, MemoryCacheProvider>();

        return services;
    }

    public static IServiceCollection AddGenomicsIntegrationWithApiKeys(
        this IServiceCollection services,
        string? ncbiApiKey = null,
        string? ensemblApiKey = null)
    {
        return services.AddGenomicsIntegration(options =>
        {
            options.Ncbi.ApiKey = ncbiApiKey;
            options.Ensembl.ApiKey = ensemblApiKey;
        });
    }
}
```

**Tests (10):**
- [ ] ServiceCollectionExtensionsTests.cs (10 tests)
  - AddGenomicsIntegration_RegistersAllServices
  - AddGenomicsIntegration_ConfiguresOptions
  - AddGenomicsIntegration_ConfiguresHttpClients
  - AddGenomicsIntegration_AddsCaching
  - ResolveBlastClient_ReturnsInstance
  - ResolveEntrezClient_ReturnsInstance
  - ResolveStructureService_ReturnsInstance
  - HttpClient_HasRetryPolicy
  - HttpClient_HasCircuitBreaker
  - HttpClient_HasTimeout

---

## Phase 8: Cross-Reference Service (Week 11)

### 8.1 ID Mapping Service

**Tasks:**
- [ ] `CrossReferenceService` — database ID mapping
- [ ] UniProt ID mapping API
- [ ] NCBI ID converter
- [ ] Ensembl xrefs

**Files:**
```
Services/
└── CrossReferenceService.cs
```

**Service:**

```csharp
// CrossReferenceService.cs
public class CrossReferenceService
{
    public async Task<Result<IReadOnlyDictionary<string, string>>> MapIdsAsync(
        IEnumerable<string> ids,
        DatabaseType fromDb,
        DatabaseType toDb,
        CancellationToken cancellationToken = default);

    public async Task<Result<IReadOnlyList<CrossReference>>> GetAllCrossReferencesAsync(
        string id,
        DatabaseType sourceDb,
        CancellationToken cancellationToken = default);

    public async Task<Result<string?>> ConvertIdAsync(
        string id,
        DatabaseType fromDb,
        DatabaseType toDb,
        CancellationToken cancellationToken = default);
}

public enum DatabaseType
{
    UniProtKB,
    RefSeq,
    GenBank,
    Ensembl,
    EnsemblGenome,
    Pdb,
    Pfam,
    InterPro,
    Go,
    Kegg,
    Ncbi_Gene,
    Ncbi_Taxonomy
}
```

**Tests (20):**
- [ ] CrossReferenceServiceTests.cs (20 tests)
  - MapIds_UniProtToRefSeq_Works
  - MapIds_RefSeqToUniProt_Works
  - MapIds_GenBankToEnsembl_Works
  - MapIds_EnsemblToPdb_Works
  - GetAllCrossReferences_ValidId_ReturnsAll
  - ConvertId_ValidConversion_ReturnsId
  - ConvertId_InvalidId_ReturnsNull
  - MapIds_BatchProcessing_Works
  - MapIds_MixedResults_HandlesPartialSuccess
  - MapIds_Caching_Works

---

## Summary

### Metrics

| Phase | Week | Components | Unit Tests | Contract Tests | Integration Tests | Total Tests |
|-------|------|------------|------------|----------------|-------------------|-------------|
| 1. Foundation | 1 | 15 | 20 | 0 | 0 | 20 |
| 2. BLAST | 2-3 | 12 | 15 | 10 | 10 | 35 |
| 3. Entrez | 4-5 | 15 | 35 | 10 | 10 | 55 |
| 4. UniProt | 6 | 10 | 20 | 5 | 5 | 30 |
| 5. Ensembl | 7 | 8 | 10 | 8 | 7 | 25 |
| 6. PDB | 8-9 | 20 | 50 | 5 | 7 | 62 |
| 7. Cross-cutting | 10 | 12 | 55 | 0 | 0 | 55 |
| 8. Cross-ref | 11 | 3 | 20 | 0 | 0 | 20 |
| **TOTAL** | **11 weeks** | **95 components** | **225** | **38** | **39** | **302 tests** |

### Test Categories

```
[Category("Unit")]        - 225 tests (fast, isolated)
[Category("Contract")]    - 38 tests (response schema validation)
[Category("Integration")] - 39 tests (real APIs, require network)
```

### CI/CD Configuration

```yaml
# Run unit tests always
- dotnet test --filter "Category=Unit"

# Run contract tests on PR
- dotnet test --filter "Category=Contract"

# Run integration tests nightly (with API keys)
- dotnet test --filter "Category=Integration"
```

---

### Definition of Done (per component)

- [ ] Implementation with XML documentation
- [ ] Nullable annotations
- [ ] Unit tests (min 5 per component)
- [ ] Contract tests for external APIs
- [ ] Integration tests (marked `[Category("Integration")]`)
- [ ] Async/await + CancellationToken support
- [ ] Resilience (retry, timeout, circuit breaker)
- [ ] Caching where applicable
- [ ] Logging via ILogger<T>
- [ ] Usage examples in XML docs

---

### Quick Start (after implementation)

```csharp
// Program.cs
var services = new ServiceCollection();
services.AddGenomicsIntegration(options =>
{
    options.Ncbi.ApiKey = "your-api-key";
    options.Ncbi.Email = "your@email.com";
    options.Resilience.MaxRetries = 3;
    options.Resilience.Timeout = TimeSpan.FromMinutes(5);
});

var provider = services.BuildServiceProvider();

// BLAST search
var blastService = provider.GetRequiredService<BlastService>();
var dna = new DnaSequence("ATGCGATCGATCGATCG...");
var result = await blastService.SearchHomologsAsync(dna);

// Fetch from GenBank
var entrez = provider.GetRequiredService<IEntrezClient>();
var record = await entrez.FetchByAccessionAsync("NM_000518");

// Get protein structure
var structureService = provider.GetRequiredService<StructureService>();
var structure = await structureService.GetStructureAsync("1CRN");

// Extension method style
var homologs = await dna.BlastAsync(blastClient);
var structures = await protein.FindStructuresAsync(structureClient);
```
