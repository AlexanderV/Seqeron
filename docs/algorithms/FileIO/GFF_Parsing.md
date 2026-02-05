# GFF/GTF Parsing

## Overview

The General Feature Format (GFF) is a tab-delimited file format used for describing genes and other features of DNA, RNA, and protein sequences. This implementation supports GFF3 (the current standard), GTF (Gene Transfer Format), and legacy GFF2.

## Format Specification

### Source
- **GFF3 Specification v1.26** (Sequence Ontology Project, Lincoln Stein, 2020)
- **UCSC Genome Browser Format Guide**

### File Structure

GFF/GTF files are nine-column, tab-delimited plain text files.

```
seqid  source  type  start  end  score  strand  phase  attributes
```

### Column Definitions

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | seqid | string | Chromosome/scaffold identifier |
| 2 | source | string | Program/database that generated the feature |
| 3 | type | string | Feature type (gene, mRNA, exon, CDS, etc.) |
| 4 | start | int | 1-based start position |
| 5 | end | int | 1-based end position (inclusive) |
| 6 | score | float? | Confidence score; "." if undefined |
| 7 | strand | char | +, -, ., or ? |
| 8 | phase | int? | 0, 1, 2 for CDS; "." otherwise |
| 9 | attributes | dict | Key-value pairs |

### Coordinate System

- **1-based, fully closed** [start, end]
- Start is always ≤ end
- Zero-length features: start equals end

### Attribute Formats

**GFF3:** `key=value;key=value`
```
ID=gene00001;Name=EDEN;Parent=transcript001
```

**GTF:** `key "value"; key "value";`
```
gene_id "ENSG00001"; gene_name "TestGene";
```

### URL Escaping (RFC 3986)

Reserved characters must be percent-encoded:

| Character | Encoding |
|-----------|----------|
| Tab | %09 |
| Newline | %0A |
| Carriage Return | %0D |
| % | %25 |
| ; | %3B |
| = | %3D |
| & | %26 |
| , | %2C |

## Implementation

### Class: `GffParser`

**Namespace:** `Seqeron.Genomics.IO`

### Data Structures

```csharp
public readonly record struct GffRecord(
    string Seqid,
    string Source,
    string Type,
    int Start,
    int End,
    double? Score,
    char Strand,
    int? Phase,
    IReadOnlyDictionary<string, string> Attributes);

public enum GffFormat { GFF3, GTF, GFF2, Auto }
```

### Core Methods

| Method | Description |
|--------|-------------|
| `Parse(string content, GffFormat)` | Parse GFF records from text content |
| `Parse(TextReader reader, GffFormat)` | Parse GFF records from a stream |
| `ParseFile(string filePath, GffFormat)` | Parse GFF records from a file |

### Filtering Methods

| Method | Description |
|--------|-------------|
| `FilterByType(records, types[])` | Filter by feature type (case-insensitive) |
| `FilterBySeqid(records, seqid)` | Filter by sequence ID |
| `FilterByRegion(records, seqid, start, end)` | Filter by genomic region overlap |
| `GetGenes(records)` | Get all gene features |
| `GetExons(records)` | Get all exon features |
| `GetCDS(records)` | Get all CDS features |

### Gene Model Building

```csharp
public readonly record struct GeneModel(
    GffRecord Gene,
    IReadOnlyList<GffRecord> Transcripts,
    IReadOnlyList<GffRecord> Exons,
    IReadOnlyList<GffRecord> CDS,
    IReadOnlyList<GffRecord> UTRs);
```

| Method | Description |
|--------|-------------|
| `BuildGeneModels(records)` | Build hierarchical gene→transcript→exon/CDS models |
| `GetGeneName(record)` | Extract gene name from attributes |
| `GetAttribute(record, name)` | Get specific attribute value |

### Statistics

```csharp
public readonly record struct GffStatistics(
    int TotalFeatures,
    IReadOnlyDictionary<string, int> FeatureTypeCounts,
    IReadOnlyList<string> SequenceIds,
    IReadOnlyList<string> Sources,
    int GeneCount,
    int ExonCount);
```

### Writing Methods

| Method | Description |
|--------|-------------|
| `WriteToFile(filePath, records, format)` | Write records to file |
| `WriteToStream(writer, records, format)` | Write records to stream |

### Utility Methods

| Method | Description |
|--------|-------------|
| `ExtractSequence(record, reference)` | Extract sequence for a feature |
| `MergeOverlapping(records)` | Merge overlapping features |

## Algorithm Details

### Format Auto-Detection

1. Scan for `##gff-version` directive
2. If contains "3" → GFF3
3. Otherwise → GFF2
4. Attribute format (key=value vs key "value") determines GFF3 vs GTF

### Line Parsing

```
1. Skip empty lines
2. Skip comment lines (# prefix)
3. Process directives (## prefix)
4. Split line by tabs
5. Validate minimum 8 fields
6. Parse columns 1-8
7. Parse column 9 attributes based on format
8. Return GffRecord
```

### Attribute Parsing

**GFF3 Algorithm:**
```
1. Split by ';'
2. For each part:
   a. Find first '='
   b. key = part[0..=]
   c. value = part[=+1..end]
   d. URL-unescape both
```

**GTF Algorithm:**
```
1. Split by ';'
2. For each part:
   a. Trim whitespace
   b. Find first space
   c. key = part[0..space]
   d. value = part[space+1..end].Trim('"')
```

### Gene Model Building

1. Index all records by ID and Parent attributes
2. Identify top-level genes (type="gene")
3. For each gene:
   - Find direct children by Parent reference
   - Classify children as transcript, exon, CDS, UTR
   - For transcripts, recursively find their children
4. Return GeneModel with hierarchical structure

## Complexity

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| Parse | O(n) | O(n) |
| FilterByType | O(n) | O(m) where m = matched |
| FilterByRegion | O(n) | O(m) |
| BuildGeneModels | O(n) | O(n) |
| MergeOverlapping | O(n log n) | O(n) |

## References

1. **GFF3 Specification v1.26** - Sequence Ontology Project
   https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

2. **UCSC Genome Browser FAQ - GFF Format**
   https://genome.ucsc.edu/FAQ/FAQformat.html

3. **Wikipedia - General Feature Format**
   https://en.wikipedia.org/wiki/General_feature_format

4. **GTF2.2 Specification** - Washington University
   http://mblab.wustl.edu/GTF22.html
