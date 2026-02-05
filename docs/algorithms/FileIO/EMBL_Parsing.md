# EMBL Parsing

## Algorithm Overview

**Category**: FileIO/Parsing  
**Test Unit**: PARSE-EMBL-001  
**Implementation**: `Seqeron.Genomics.EmblParser`

## Description

The EMBL parser reads and interprets EMBL flat file format records, which is the format used by the European Nucleotide Archive (ENA) for storing nucleotide sequence data. The parser extracts all metadata, features, and sequence data from EMBL records.

## EMBL Format Specification

EMBL format is a line-based format where each line begins with a two-character line type code followed by three spaces. Records are terminated by `//`.

### Line Types

| Code | Description | Occurrence |
|------|-------------|------------|
| ID | Identification (accession, topology, length) | 1 per entry |
| AC | Accession number(s) | ≥1 per entry |
| SV | Sequence version (legacy, now in ID) | 0-1 per entry |
| DT | Date (created, updated) | 2 per entry |
| DE | Description | ≥1 per entry |
| KW | Keywords | ≥1 per entry |
| OS | Organism species | ≥1 per entry |
| OC | Organism classification | ≥1 per entry |
| OG | Organelle | 0-1 per entry |
| RN | Reference number | ≥1 per entry |
| RA | Reference authors | ≥0 per entry |
| RT | Reference title | ≥1 per entry |
| RL | Reference location | ≥1 per entry |
| RX | Reference cross-reference | ≥0 per entry |
| FH | Feature table header | 2 per entry |
| FT | Feature table data | ≥2 per entry |
| XX | Spacer line | Many per entry |
| SQ | Sequence header | 1 per entry |
| bb | Sequence data (blank prefix) | ≥1 per entry |
| // | Terminator | 1 per entry |

### ID Line Format

```
ID   <accession>; SV <version>; <topology>; <molecule>; <class>; <division>; <length> BP.
```

Example:
```
ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
```

### Feature Table Format

```
FT   <key>           <location>
FT                   /qualifier="value"
```

## API Reference

### Core Methods

```csharp
// Parse from string content
public static IEnumerable<EmblRecord> Parse(string content);

// Parse from file
public static IEnumerable<EmblRecord> ParseFile(string filePath);

// Parse feature locations
public static Location ParseLocation(string locationStr);

// Convert to GenBank format
public static GenBankParser.GenBankRecord ToGenBank(EmblRecord embl);

// Filter features by type
public static IEnumerable<Feature> GetFeatures(EmblRecord record, string featureKey);
public static IEnumerable<Feature> GetCDS(EmblRecord record);
public static IEnumerable<Feature> GetGenes(EmblRecord record);

// Extract subsequence
public static string ExtractSequence(EmblRecord record, Location location);
```

### Record Structure

```csharp
public readonly record struct EmblRecord(
    string Accession,
    string SequenceVersion,
    string DataClass,
    string MoleculeType,
    string Topology,
    string TaxonomicDivision,
    int SequenceLength,
    string Description,
    IReadOnlyList<string> Keywords,
    string Organism,
    IReadOnlyList<string> OrganismClassification,
    IReadOnlyList<Reference> References,
    IReadOnlyList<Feature> Features,
    string Sequence,
    IReadOnlyDictionary<string, string> AdditionalFields);
```

## Usage Examples

### Parse EMBL Content

```csharp
var content = File.ReadAllText("sequence.embl");
var records = EmblParser.Parse(content);

foreach (var record in records)
{
    Console.WriteLine($"Accession: {record.Accession}");
    Console.WriteLine($"Organism: {record.Organism}");
    Console.WriteLine($"Length: {record.SequenceLength} bp");
    Console.WriteLine($"Features: {record.Features.Count}");
}
```

### Extract CDS Features

```csharp
var record = EmblParser.Parse(content).First();
var cdsFeatures = EmblParser.GetCDS(record);

foreach (var cds in cdsFeatures)
{
    var product = cds.Qualifiers.GetValueOrDefault("product", "unknown");
    Console.WriteLine($"CDS: {cds.Location.RawLocation} - {product}");
}
```

### Parse Location Strings

```csharp
// Simple range
var loc1 = EmblParser.ParseLocation("100..200");
// loc1.Start = 100, loc1.End = 200

// Complement
var loc2 = EmblParser.ParseLocation("complement(100..200)");
// loc2.IsComplement = true

// Join
var loc3 = EmblParser.ParseLocation("join(1..50,60..100)");
// loc3.IsJoin = true, loc3.Parts.Count = 2
```

### Convert to GenBank

```csharp
var emblRecord = EmblParser.Parse(content).First();
var genBankRecord = EmblParser.ToGenBank(emblRecord);
```

## Location Syntax

| Pattern | Description | Example |
|---------|-------------|---------|
| `n..m` | Range from n to m | `100..200` |
| `n` | Single position | `467` |
| `<n..m` | Partial start (5') | `<1..200` |
| `n..>m` | Partial end (3') | `100..>500` |
| `complement(loc)` | Reverse complement | `complement(100..200)` |
| `join(loc1,loc2,...)` | Joined regions | `join(1..50,60..100)` |

## Data Classes

| Class | Description |
|-------|-------------|
| STD | Standard |
| CON | Constructed |
| PAT | Patent |
| EST | Expressed Sequence Tag |
| GSS | Genome Survey Sequence |
| HTG | High Throughput Genome |
| WGS | Whole Genome Shotgun |
| TSA | Transcriptome Shotgun Assembly |

## Taxonomic Divisions

| Code | Division |
|------|----------|
| PHG | Bacteriophage |
| FUN | Fungal |
| HUM | Human |
| INV | Invertebrate |
| PLN | Plant |
| PRO | Prokaryote |
| VRL | Viral |
| SYN | Synthetic |

## Test Coverage

- **Test File**: `tests/Seqeron/Seqeron.Genomics.Tests/EmblParserTests.cs`
- **Test Spec**: `tests/TestSpecs/PARSE-EMBL-001.md`
- **Evidence**: `docs/Evidence/PARSE-EMBL-001-Evidence.md`

### Test Categories

1. **Basic Parsing**: Valid record, empty content, null content
2. **ID Line**: Accession, version, topology, molecule type, data class, division, length
3. **Metadata**: Description, keywords, organism, classification
4. **References**: Number, authors, title, journal
5. **Features**: Extraction, qualifiers, filtering
6. **Sequence**: Extraction, normalization, case handling
7. **Locations**: Simple range, complement, join, partial
8. **Conversion**: ToGenBank preservation
9. **File Operations**: ParseFile, invalid paths

## References

- [EBI EMBL User Manual](https://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt)
- [INSDC Feature Table Definition](https://www.insdc.org/files/feature_table.html)
- [European Nucleotide Archive](https://www.ebi.ac.uk/ena/)

## Changelog

| Date | Version | Changes |
|------|---------|---------|
| 2025-01-28 | 1.0 | Initial documentation |
