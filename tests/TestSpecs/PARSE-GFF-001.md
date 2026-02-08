# TestSpec: PARSE-GFF-001 - GFF/GTF Parsing

## Test Unit Information
- **ID:** PARSE-GFF-001
- **Area:** FileIO
- **Class:** GffParser
- **Canonical Methods:** Parse, ParseFile, ToGff3
- **Evidence:** [PARSE-GFF-001-Evidence.md](../../docs/Evidence/PARSE-GFF-001-Evidence.md)

---

## Test Classification

### Must Tests (Evidence-Based)

| ID | Test Name | Source | Rationale |
|----|-----------|--------|-----------|
| M01 | Parse_GFF3_ReturnsCorrectRecords | GFF3 Spec | 9-column tab-delimited format |
| M02 | Parse_GTF_ReturnsCorrectRecords | UCSC/GTF Spec | GTF attribute format key "value" |
| M03 | Parse_EmptyContent_ReturnsEmpty | Standard behavior | Edge case handling |
| M04 | Parse_NullContent_ReturnsEmpty | Standard behavior | Null safety |
| M05 | Parse_SkipsComments | GFF3 Spec | Lines starting with # are comments |
| M06 | Parse_SkipsEmptyLines | GFF3 Spec | Blank lines ignored |
| M07 | Parse_GFF3Attributes_ParsedCorrectly | GFF3 Spec | key=value;key=value format |
| M08 | Parse_GTFAttributes_ParsedCorrectly | GTF Spec | key "value"; format |
| M09 | Parse_MalformedLine_Skips | GFF3 Spec | Lines with < 8 fields skipped |
| M10 | Parse_NoScore_ScoreIsNull | GFF3 Spec | "." means undefined |
| M11 | Parse_WithScore_ScoreIsParsed | GFF3 Spec | Floating point score |
| M12 | Parse_SpecialCharacters_Unescaped | GFF3 Spec RFC 3986 | URL percent-encoding |
| M13 | Parse_DetectsGFF3Version | GFF3 Spec | ##gff-version directive |
| M14 | Parse_1BasedCoordinates | GFF3 Spec | Start/end are 1-based |
| M15 | Parse_Phase_ParsedCorrectly | GFF3 Spec | CDS phase 0/1/2 |
| M16 | Parse_Strand_AllValidValues | GFF3 Spec | +, -, ., ? |
| M17 | FilterByType_ReturnsMatchingTypes | API Contract | Type filtering |
| M18 | FilterBySeqid_ReturnsMatchingChromosome | API Contract | Seqid filtering |
| M19 | FilterByRegion_ReturnsOverlappingFeatures | API Contract | Genomic region overlap |
| M20 | BuildGeneModels_CreatesHierarchy | GFF3 Spec | Parent-child relationships |
| M21 | CalculateStatistics_ReturnsCorrectCounts | API Contract | Feature counting |
| M22 | WriteToStream_GFF3Format_ValidOutput | GFF3 Spec | Correct output format |
| M23 | WriteAndRead_Roundtrip_PreservesData | Standard behavior | Data integrity |
| M24 | ParseFile_NonexistentFile_ReturnsEmpty | Error handling | Graceful failure |
| M25 | ParseFile_ValidFile_ParsesRecords | API Contract | File I/O |

### Should Tests (Quality/Robustness)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| S01 | FilterByType_MultipleTypes_ReturnsAll | Multiple type filtering |
| S02 | GetGenes_ReturnsOnlyGenes | Convenience method |
| S03 | GetExons_ReturnsOnlyExons | Convenience method |
| S04 | GetCDS_ReturnsOnlyCDS | Convenience method |
| S05 | BuildGeneModels_MultipleGenes_BuildsAll | Multiple gene support |
| S06 | CalculateStatistics_FeatureTypeCounts_Correct | Type breakdown |
| S07 | CalculateStatistics_SequenceIds_Listed | Seqid enumeration |
| S08 | GetAttribute_ExistingAttribute_ReturnsValue | Attribute access |
| S09 | GetAttribute_NonexistentAttribute_ReturnsNull | Missing attribute handling |
| S10 | GetGeneName_ReturnsGeneName | Name extraction (Name/gene_name/gene_id) |
| S11 | ExtractSequence_PlusStrand_ReturnsSequence | Sequence extraction |
| S12 | ExtractSequence_MinusStrand_ReturnsReverseComplement | Reverse complement |
| S13 | MergeOverlapping_MergesCorrectly | Feature merging |

### Could Tests (Extended Coverage)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C01 | Parse_MultipleAttributeValues | Comma-separated Parent values |
| C02 | Parse_CaseInsensitiveAttributeLookup | Attribute key case |
| C03 | Parse_LargeFile_Performance | Performance benchmark |
| C04 | Parse_CircularGenome_IsCircularAttribute | Circular genome support |
| C05 | FilterByRegion_NoOverlap_ReturnsEmpty | Region filtering edge case |

---

## Audit Results

### Existing Test File: GffParserTests.cs

| Test Name | Category | Status | Notes |
|-----------|----------|--------|-------|
| Parse_GFF3_ReturnsCorrectRecords | M01 | ✓ Covered | |
| Parse_GTF_ReturnsCorrectRecords | M02 | ✓ Covered | |
| Parse_EmptyContent_ReturnsEmpty | M03 | ✓ Covered | |
| Parse_NullContent_ReturnsEmpty | M04 | ✓ Covered | |
| Parse_SkipsComments | M05 | ✓ Covered | |
| Parse_SkipsEmptyLines | M06 | ✓ Covered | |
| Parse_GFF3Attributes_ParsedCorrectly | M07 | ✓ Covered | |
| Parse_GTFAttributes_ParsedCorrectly | M08 | ✓ Covered | |
| Parse_MalformedLine_Skips | M09 | ✓ Covered | |
| Parse_NoScore_ScoreIsNull | M10 | ✓ Covered | |
| Parse_WithScore_ScoreIsParsed | M11 | ✓ Covered | |
| Parse_SpecialCharacters_Unescaped | M12 | ✓ Covered | |
| FilterByType_ReturnsMatchingTypes | M17 | ✓ Covered | |
| FilterByType_MultipleTypes_ReturnsAll | S01 | ✓ Covered | |
| FilterBySeqid_ReturnsMatchingChromosome | M18 | ✓ Covered | |
| FilterByRegion_ReturnsOverlappingFeatures | M19 | ✓ Covered | |
| GetGenes_ReturnsOnlyGenes | S02 | ✓ Covered | |
| GetExons_ReturnsOnlyExons | S03 | ✓ Covered | |
| GetCDS_ReturnsOnlyCDS | S04 | ✓ Covered | |
| BuildGeneModels_CreatesHierarchy | M20 | ✓ Covered | |
| BuildGeneModels_MultipleGenes_BuildsAll | S05 | ✓ Covered | |
| CalculateStatistics_ReturnsCorrectCounts | M21 | ✓ Covered | |
| CalculateStatistics_FeatureTypeCounts_Correct | S06 | ✓ Covered | |
| CalculateStatistics_SequenceIds_Listed | S07 | ✓ Covered | |
| GetAttribute_ExistingAttribute_ReturnsValue | S08 | ✓ Covered | |
| GetAttribute_NonexistentAttribute_ReturnsNull | S09 | ✓ Covered | |
| GetGeneName_ReturnsGeneName | S10 | ✓ Covered | |
| WriteToStream_GFF3Format_ValidOutput | M22 | ✓ Covered | |
| WriteAndRead_Roundtrip_PreservesData | M23 | ✓ Covered | |
| ExtractSequence_PlusStrand_ReturnsSequence | S11 | ✓ Covered | |
| ExtractSequence_MinusStrand_ReturnsReverseComplement | S12 | ✓ Covered | |
| MergeOverlapping_MergesCorrectly | S13 | ✓ Covered | |
| ParseFile_NonexistentFile_ReturnsEmpty | M24 | ✓ Covered | |
| ParseFile_ValidFile_ParsesRecords | M25 | ✓ Covered | |

### Missing Tests (All Closed):

| Test ID | Description | Status |
|---------|-------------|--------|
| M13 | Parse_DetectsGFF3Version | ✅ Covered |
| M14 | Parse_1BasedCoordinates_Validated | ✅ Covered |
| M15 | Parse_Phase_ParsedCorrectly | ✅ Covered |
| M16 | Parse_Strand_AllValidValues | ✅ Covered |
| C01 | Parse_MultipleParentValues | ✅ Covered |
| C02 | Parse_AttributeCaseInsensitive | ✅ Covered |

### MCP Wrapper Tests: GffParseTests.cs

| Test Name | Status | Notes |
|-----------|--------|-------|
| GffParse_Schema_ValidatesCorrectly | ✓ Smoke test | Delegation verification |
| GffParse_Binding_ParsesRecords | ✓ Smoke test | |
| GffParse_Binding_ParsesAttributes | ✓ Smoke test | |
| GffParse_Binding_ParsesScore | ✓ Smoke test | |
| GffParse_Binding_CalculatesLength | ✓ Smoke test | |
| GffParse_Binding_RespectsFormat | ✓ Smoke test | |
| GffStatistics_* | ✓ Smoke tests | Statistics API |
| GffFilter_* | ✓ Smoke tests | Filtering API |

**Wrapper Assessment:** MCP tests provide appropriate smoke coverage; no deep logic tests needed there.

---

## Consolidation Plan

1. **Canonical test file:** `GffParserTests.cs` - ALL deep logic tests
2. **Wrapper test file:** `GffParseTests.cs` - Smoke tests only (current state is appropriate)
3. ~~**Action:** Add missing Must tests M13-M16 to canonical file~~ ✅ Done
4. ~~**Action:** Optionally add Could tests C01-C02~~ ✅ Done

---

## Test Data

### Standard GFF3 Test Data (from existing tests):
```gff3
##gff-version 3
chr1	ENSEMBL	gene	1000	5000	.	+	.	ID=gene1;Name=TestGene
chr1	ENSEMBL	mRNA	1000	5000	.	+	.	ID=transcript1;Parent=gene1
chr1	ENSEMBL	exon	1000	1500	.	+	.	ID=exon1;Parent=transcript1
chr1	ENSEMBL	CDS	1100	1500	.	+	0	ID=cds1;Parent=transcript1
```

### GTF Test Data:
```gtf
chr1	ENSEMBL	gene	1000	5000	.	+	.	gene_id "ENSG00001"; gene_name "TestGene";
chr1	ENSEMBL	transcript	1000	5000	.	+	.	gene_id "ENSG00001"; transcript_id "ENST00001";
```

---

## Open Questions

None - GFF3 specification is authoritative and well-documented.

---

## Sign-off

- **Author:** Algorithm QA Architect
- **Date:** 2026-02-05
- **Status:** Ready for implementation
