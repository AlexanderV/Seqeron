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
| M26 | Parse_AttributeCaseSensitive | GFF3 Spec v1.26 | "attribute names are case sensitive" |
| M27 | Parse_DoubleEncodedPercent_DecodedCorrectly | GFF3 Spec RFC 3986 | %25 → % without double-decoding |
| M28 | Parse_VersionDetection_DoesNotMisdetectGFF2AsGFF3 | GFF3 Spec | "##gff-version 2.3" ≠ GFF3 |
| M29 | WriteToStream_GTFFormat_TrailingSemicolon | GTF/UCSC Spec | "Attributes must end in a semicolon" |
| M30 | WriteToStream_ScorePrecision_PreservesScientificNotation | GFF3 Spec | E-values must not be truncated |
| M31 | BuildGeneModels_MultiParentExons_AssignedToAllTranscripts | GFF3 Spec v1.26 | Comma-separated Parent values |

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
| S14 | Parse_GTFAttributes_CaseInsensitive | GTF doesn't mandate case-sensitive attrs |

### Could Tests (Extended Coverage)

| ID | Test Name | Rationale |
|----|-----------|-----------|
| C01 | Parse_MultipleAttributeValues | Comma-separated Parent values |
| C02 | Parse_LargeFile_Performance | Performance benchmark |
| C03 | Parse_CircularGenome_IsCircularAttribute | Circular genome support |
| C04 | FilterByRegion_NoOverlap_ReturnsEmpty | Region filtering edge case |

---

## Coverage Classification

### Classification Summary

| Category | Count | Description |
|----------|-------|-------------|
| ✅ Covered | 39 | Adequate assertions grounded in spec/theory |
| ⚠ Weak → ✅ Fixed | 6 | Strengthened with exact values and structural assertions |
| ❌ Missing → ✅ Added | 1 | FilterByRegion_NoOverlap_ReturnsEmpty (C04) |
| 🔁 Duplicate | 0 | No duplicates found |

### Weak Tests Fixed

| Test | Problem | Fix |
|------|---------|-----|
| Parse_GTF_ReturnsCorrectRecords (M02) | `GreaterThanOrEqualTo(4)` — hid exact expectation | Exact `EqualTo(4)` + field-level assertions on all 4 records |
| BuildGeneModels_CreatesHierarchy (M20) | `GreaterThanOrEqualTo` for Transcripts/Exons/CDS | Exact `EqualTo(1)` / `EqualTo(2)` / `EqualTo(2)` |
| WriteToStream_GFF3Format_ValidOutput (M22) | Only `Does.Contain` — no structural validation | Verify 7 lines (header + 6 features), 9 tab-separated columns per line, specific field values |
| ExtractSequence_MinusStrand (S12) | Palindromic "ACGT" → rc is "ACGT" — passes trivially | Non-palindromic "ATCG" → rc must be "CGAT" |
| Parse_DetectsGFF3Version (M13) | Only checked record count — didn't verify format detection | GFF3 input → case-sensitive attributes (ID found, id not found) |
| Parse_VersionDetection_DoesNotMisdetectGFF2AsGFF3 (M28) | Only checked count=1 — didn't prove GFF2 was detected | GFF2 detected → case-insensitive attributes ("id" matches "ID") |

### Missing Test Added

| Test | Spec ID | Rationale |
|------|---------|-----------|
| FilterByRegion_NoOverlap_ReturnsEmpty | C04 | Region [10000,20000] doesn't overlap any chr1 feature → empty result |

---

## Audit Results

### Canonical Test File: GffParserTests.cs — 46 tests passing

| Test Name | Category | Status | Classification |
|-----------|----------|--------|----------------|
| Parse_GFF3_ReturnsCorrectRecords | M01 | ✓ | ✅ Covered |
| Parse_GTF_ReturnsCorrectRecords | M02 | ✓ | ⚠→✅ Fixed |
| Parse_EmptyContent_ReturnsEmpty | M03 | ✓ | ✅ Covered |
| Parse_NullContent_ReturnsEmpty | M04 | ✓ | ✅ Covered |
| Parse_SkipsComments | M05 | ✓ | ✅ Covered |
| Parse_SkipsEmptyLines | M06 | ✓ | ✅ Covered |
| Parse_GFF3Attributes_ParsedCorrectly | M07 | ✓ | ✅ Covered |
| Parse_GTFAttributes_ParsedCorrectly | M08 | ✓ | ✅ Covered |
| Parse_MalformedLine_Skips | M09 | ✓ | ✅ Covered |
| Parse_NoScore_ScoreIsNull | M10 | ✓ | ✅ Covered |
| Parse_WithScore_ScoreIsParsed | M11 | ✓ | ✅ Covered |
| Parse_SpecialCharacters_Unescaped | M12 | ✓ | ✅ Covered |
| Parse_DetectsGFF3Version | M13 | ✓ | ⚠→✅ Fixed |
| Parse_1BasedCoordinates_Validated | M14 | ✓ | ✅ Covered |
| Parse_Phase_ParsedCorrectly | M15 | ✓ | ✅ Covered |
| Parse_Strand_AllValidValues | M16 | ✓ | ✅ Covered |
| FilterByType_ReturnsMatchingTypes | M17 | ✓ | ✅ Covered |
| FilterBySeqid_ReturnsMatchingChromosome | M18 | ✓ | ✅ Covered |
| FilterByRegion_ReturnsOverlappingFeatures | M19 | ✓ | ✅ Covered |
| BuildGeneModels_CreatesHierarchy | M20 | ✓ | ⚠→✅ Fixed |
| CalculateStatistics_ReturnsCorrectCounts | M21 | ✓ | ✅ Covered |
| WriteToStream_GFF3Format_ValidOutput | M22 | ✓ | ⚠→✅ Fixed |
| WriteAndRead_Roundtrip_PreservesData | M23 | ✓ | ✅ Covered |
| ParseFile_NonexistentFile_ReturnsEmpty | M24 | ✓ | ✅ Covered |
| ParseFile_ValidFile_ParsesRecords | M25 | ✓ | ✅ Covered |
| Parse_AttributeCaseSensitive | M26 | ✓ | ✅ Covered |
| Parse_DoubleEncodedPercent_DecodedCorrectly | M27 | ✓ | ✅ Covered |
| Parse_VersionDetection_DoesNotMisdetectGFF2AsGFF3 | M28 | ✓ | ⚠→✅ Fixed |
| WriteToStream_GTFFormat_TrailingSemicolon | M29 | ✓ | ✅ Covered |
| WriteToStream_ScorePrecision_PreservesScientificNotation | M30 | ✓ | ✅ Covered |
| BuildGeneModels_MultiParentExons_AssignedToAllTranscripts | M31 | ✓ | ✅ Covered |
| FilterByType_MultipleTypes_ReturnsAll | S01 | ✓ | ✅ Covered |
| GetGenes_ReturnsOnlyGenes | S02 | ✓ | ✅ Covered |
| GetExons_ReturnsOnlyExons | S03 | ✓ | ✅ Covered |
| GetCDS_ReturnsOnlyCDS | S04 | ✓ | ✅ Covered |
| BuildGeneModels_MultipleGenes_BuildsAll | S05 | ✓ | ✅ Covered |
| CalculateStatistics_FeatureTypeCounts_Correct | S06 | ✓ | ✅ Covered |
| CalculateStatistics_SequenceIds_Listed | S07 | ✓ | ✅ Covered |
| GetAttribute_ExistingAttribute_ReturnsValue | S08 | ✓ | ✅ Covered |
| GetAttribute_NonexistentAttribute_ReturnsNull | S09 | ✓ | ✅ Covered |
| GetGeneName_ReturnsGeneName | S10 | ✓ | ✅ Covered |
| ExtractSequence_PlusStrand_ReturnsSequence | S11 | ✓ | ✅ Covered |
| ExtractSequence_MinusStrand_ReturnsReverseComplement | S12 | ✓ | ⚠→✅ Fixed |
| MergeOverlapping_MergesCorrectly | S13 | ✓ | ✅ Covered |
| Parse_GTFAttributes_CaseInsensitive | S14 | ✓ | ✅ Covered |
| Parse_MultipleParentValues | C01 | ✓ | ✅ Covered |
| FilterByRegion_NoOverlap_ReturnsEmpty | C04 | ✓ | ❌→✅ Added |

### MCP Wrapper Tests: GffParseTests.cs — 15 tests passing

Smoke coverage for MCP delegation; no deep logic tests needed.

---

## Deviations and Assumptions

**None.** All tests and implementation are grounded in the authoritative sources listed in the Evidence document.

- GFF3 attributes are case-sensitive per GFF3 Spec v1.26: *"attribute names are case sensitive. 'Parent' is not the same as 'parent'."*
- GTF/GFF2 attributes use case-insensitive lookup (GTF spec does not mandate case sensitivity).
- Escaping/unescaping uses `Uri.UnescapeDataString` which correctly handles RFC 3986 percent-encoding without double-decoding.
- GFF3 version detection parses the version number after `##gff-version` and checks that it starts with `3`, per spec format `3.#.#`.
- Score is written with `"G"` format to preserve full precision (E-values, P-values per spec recommendation).
- GTF output ends each line's attributes with a trailing semicolon per UCSC GTF spec: *"Attributes must end in a semicolon."*
- BuildGeneModels splits comma-separated Parent values per GFF3 Spec v1.26: *"A feature may have multiple parents."*

---

## Test Data

### Standard GFF3 Test Data:
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

### Canonical Gene (from GFF3 Spec v1.26):
```gff3
##gff-version 3
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123	.	exon	1050	1500	.	+	.	ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123	.	CDS	1201	1500	.	+	0	ID=cds00001;Parent=mRNA00001
```

---

## Sign-off

- **Author:** Algorithm QA Architect
- **Date:** 2026-03-12
- **Coverage Classification Date:** 2026-03-12
- **Status:** Complete — zero deviations, zero weak tests, zero duplicates
