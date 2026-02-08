# TestSpec: PARSE-VCF-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **ID** | PARSE-VCF-001 |
| **Area** | FileIO |
| **Title** | VCF Parsing |
| **Canonical Class** | `VcfParser` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete |

---

## Methods Under Test

| Method | Type | Description |
|--------|------|-------------|
| `Parse(content)` | Canonical | Parse VCF content to records |
| `ParseFile(filePath)` | File I/O | Parse VCF from file |
| `ParseWithHeader(content)` | Extended | Parse with header extraction |
| `ClassifyVariant(record)` | Classification | Determine variant type |

### Supporting Methods

| Method | Category |
|--------|----------|
| `IsSNP(record)` | Type check |
| `IsIndel(record)` | Type check |
| `GetVariantLength(record)` | Measurement |
| `FilterByChrom(records, chrom)` | Filtering |
| `FilterByRegion(records, chrom, start, end)` | Filtering |
| `FilterByQuality(records, minQuality)` | Filtering |
| `FilterPassing(records)` | Filtering |
| `FilterSNPs(records)` | Filtering |
| `FilterIndels(records)` | Filtering |
| `FilterByInfo(records, key, predicate)` | Filtering |
| `GetGenotype(record, sampleIndex)` | Genotype |
| `IsHomRef(record, sampleIndex)` | Genotype |
| `IsHomAlt(record, sampleIndex)` | Genotype |
| `IsHet(record, sampleIndex)` | Genotype |
| `GetReadDepth(record, sampleIndex)` | Sample data |
| `GetAlleleDepth(record, sampleIndex)` | Sample data |
| `CalculateStatistics(records)` | Statistics |
| `CalculateTiTvRatio(records)` | Statistics |
| `WriteToStream(writer, records)` | Writing |
| `GetInfoValue(record, key)` | INFO access |
| `GetInfoInt(record, key)` | INFO access |
| `GetInfoDouble(record, key)` | INFO access |
| `HasInfoFlag(record, key)` | INFO access |

---

## Test Categories

### Must Tests (Evidence-Based)

| ID | Test | Source | Rationale |
|----|------|--------|-----------|
| M01 | Parse_EmptyContent_ReturnsEmpty | VCF Spec | Empty input handling |
| M02 | Parse_NullContent_ReturnsEmpty | VCF Spec | Null safety |
| M03 | Parse_ValidVcf_ReturnsRecords | Wikipedia | Basic parsing |
| M04 | Parse_SkipsMetadataLines | VCF Spec | ## lines are metadata |
| M05 | Parse_MultipleAlternates_ParsedCorrectly | VCF Spec | Multi-allelic sites |
| M06 | ParseWithHeader_ReturnsHeaderAndRecords | VCF Spec | Header extraction |
| M07 | ParseWithHeader_SampleNames_Parsed | VCF Spec | Sample column parsing |
| M08 | ParseWithHeader_InfoFields_Parsed | VCF Spec | INFO metadata parsing |
| M09 | ClassifyVariant_SNP_ReturnsSnp | Wikipedia | Single-base substitution |
| M10 | ClassifyVariant_Insertion_ReturnsInsertion | Wikipedia | Alt > Ref |
| M11 | ClassifyVariant_Deletion_ReturnsDeletion | Wikipedia | Ref > Alt |
| M12 | IsSNP_ForSnp_ReturnsTrue | Wikipedia | SNP detection |
| M13 | IsIndel_ForInsertion_ReturnsTrue | Wikipedia | Indel detection |
| M14 | GetVariantLength_Insertion_ReturnsCorrectLength | VCF Spec | Length calculation |
| M15 | GetVariantLength_Deletion_ReturnsCorrectLength | VCF Spec | Length calculation |
| M16 | FilterByChrom_ReturnsMatchingChromosome | VCF Spec | Chromosome filtering |
| M17 | FilterByRegion_ReturnsVariantsInRegion | VCF Spec | Position filtering |
| M18 | FilterByQuality_FiltersLowQuality | VCF Spec | Quality filtering |
| M19 | FilterPassing_ReturnsOnlyPassing | VCF Spec | PASS filter |
| M20 | FilterSNPs_ReturnsOnlySNPs | Wikipedia | Type filtering |
| M21 | FilterIndels_ReturnsOnlyIndels | Wikipedia | Type filtering |
| M22 | FilterByInfo_FiltersCorrectly | VCF Spec | INFO filtering |
| M23 | GetGenotype_ReturnsGenotype | VCF Spec | GT extraction |
| M24 | IsHomRef_ForHomRef_ReturnsTrue | Wikipedia | 0/0 detection |
| M25 | IsHomAlt_ForHomAlt_ReturnsTrue | Wikipedia | 1/1 detection |
| M26 | IsHet_ForHet_ReturnsTrue | Wikipedia | 0/1 detection |
| M27 | GetReadDepth_ReturnsDepth | VCF Spec | DP extraction |
| M28 | CalculateStatistics_ReturnsCorrectCounts | VCF Spec | Statistics |
| M29 | CalculateStatistics_PassingCount_Correct | VCF Spec | Pass count |
| M30 | CalculateStatistics_ChromosomeCounts_Correct | VCF Spec | Chrom distribution |
| M31 | CalculateTiTvRatio_ReturnsRatio | Danecek 2011 | Ti/Tv calculation |
| M32 | GetInfoValue_ReturnsValue | VCF Spec | INFO access |
| M33 | GetInfoValue_NonexistentKey_ReturnsNull | VCF Spec | Missing key handling |
| M34 | GetInfoInt_ReturnsInteger | VCF Spec | Integer parsing |
| M35 | GetInfoDouble_ReturnsDouble | VCF Spec | Float parsing |
| M36 | HasInfoFlag_ForExistingFlag_ReturnsTrue | VCF Spec | Flag detection |
| M37 | Parse_MissingQuality_QualIsNull | VCF Spec | "." quality handling |
| M38 | Parse_MissingId_IdIsDot | VCF Spec | "." ID handling |
| M39 | Parse_EmptyFilter_FilterIsEmpty | VCF Spec | "." filter handling |
| M40 | Parse_MultipleFilters_AllParsed | VCF Spec | Semicolon-separated |
| M41 | Parse_MalformedLine_Skips | VCF Spec | Error resilience |
| M42 | WriteAndRead_Roundtrip_PreservesData | Implementation | Data integrity |

### Should Tests (Extended Validation)

| ID | Test | Rationale |
|----|------|-----------|
| S01 | ParseFile_NonexistentFile_ReturnsEmpty | File error handling |
| S02 | ParseFile_ValidFile_ParsesRecords | File I/O integration |
| S03 | ParseFileWithHeader_ValidFile_ReturnsHeaderAndRecords | File I/O with header |
| S04 | WriteToStream_ProducesValidVcf | Output formatting |
| S05 | GetAlleleDepth_ReturnsDepths | AD field handling |
| S06 | ClassifyVariant_MNP_ReturnsMnp | MNP classification |
| S07 | ClassifyVariant_Symbolic_ReturnsSymbolic | <DEL> etc. handling |

### Could Tests (Nice-to-Have)

| ID | Test | Rationale |
|----|------|-----------|
| C01 | Parse_LargeFile_PerformsEfficiently | Performance baseline |
| C02 | CalculateAlleleFrequency_ReturnsFrequency | Population analysis |

---

## Test Data

### SimpleVcf
```vcf
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
chr1	100	rs123	A	G	99	PASS	DP=50	GT:DP	0/1:25	1/1:30
chr1	200	.	C	T	50	.	DP=30	GT:DP	0/0:15	0/1:20
chr2	300	rs456	G	A,C	80	PASS	DP=40	GT:DP	1/2:20	0/1:25
```

### VcfWithIndels
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	AT	99	PASS	.
chr1	200	.	ATG	A	99	PASS	.
chr1	300	.	A	G	99	PASS	.
```

### VcfWithFilters
```vcf
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	.
chr1	200	.	C	T	20	LowQual	.
chr1	300	.	G	A	10	LowQual;LowCov	.
```

---

## Invariants

1. **Position is 1-based positive integer**
2. **Quality is null or non-negative**
3. **FILTER empty or contains PASS indicates passing**
4. **Variant length = |len(ALT) - len(REF)|**
5. **Ti/Tv ratio is null when no transversions exist**

---

## Audit Notes

### Existing Test Analysis (VcfParserTests.cs)

**Current Coverage:** 46 tests covering:
- Basic parsing (4 tests)
- Header parsing (3 tests)
- Variant classification (7 tests)
- Filtering (7 tests)
- Genotype analysis (6 tests)
- Statistics (5 tests)
- INFO helpers (5 tests)
- Writing (2 tests)
- Edge cases (5 tests)
- File I/O (3 tests)

**Assessment:** Comprehensive coverage already exists. Tests are well-organized and follow NUnit conventions.

### MCP Wrapper Tests (VcfParseTests.cs, VcfUtilityTests.cs)

**Assessment:** These test the MCP layer wrappers. They provide smoke tests that delegate to VcfParser. Per the test pool policy, these remain in place with minimal verification tests.

### Consolidation Decision

- **Canonical tests:** VcfParserTests.cs (46 tests) - comprehensive, no changes needed
- **Wrapper tests:** Mcp.Parsers.Tests - smoke tests only, appropriate
- **Action:** Review existing tests for completeness against TestSpec; add any missing Must tests

---

## Missing Tests to Add (All Closed)

After audit, the following tests were added:

| Test | Category | Status |
|------|----------|--------|
| ClassifyVariant_MNP_ReturnsMnp | M | ✅ Covered |
| ClassifyVariant_Symbolic_ReturnsSymbolic | M | ✅ Covered |
| ClassifyVariant_Complex_ReturnsComplex | M | ✅ Covered |
| Parse_PhasedGenotypes_HandledCorrectly | M | ✅ Covered |
| IsHomAlt_MultiAllelic_DetectsCorrectly | M | ✅ Covered |

---

## Definition of Done

- [x] Evidence documented with sources
- [x] Algorithm documentation created
- [x] TestSpec with Must/Should/Could categories
- [x] Existing tests audited
- [x] Missing tests added
- [x] All tests pass
- [x] Zero warnings
- [x] Checklist updated
