# GenBank Parsing

## Overview

GenBank is the NIH genetic sequence database, maintained by the National Center for Biotechnology Information (NCBI). The GenBank flat file format is a text-based format for storing nucleotide sequences with associated metadata and annotations.

**Sources:**
- NCBI GenBank Sample Record: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
- Wikipedia GenBank: https://en.wikipedia.org/wiki/GenBank

---

## Format Structure

### Record Layout

A GenBank record consists of sections identified by keywords starting in column 1:

```
LOCUS       name         length bp    type    topology division date
DEFINITION  description
ACCESSION   accession_number
VERSION     accession.version
KEYWORDS    keywords or "."
SOURCE      organism_name
  ORGANISM  scientific_name
            taxonomy; lineage; here.
REFERENCE   n  (bases x to y)
  AUTHORS   author_list
  TITLE     title
  JOURNAL   journal
  PUBMED    pmid
FEATURES             Location/Qualifiers
     feature_key     location
                     /qualifier="value"
ORIGIN      
        1 sequence_data_in_60_char_lines
//
```

### LOCUS Line Fields

| Field | Position | Description |
|-------|----------|-------------|
| LOCUS | 1-5 | Keyword |
| Name | 13-28 | Locus name (unique identifier) |
| Length | 30-40 | Sequence length + "bp" or "aa" |
| Type | 45-47 | DNA, RNA, mRNA, etc. |
| Topology | 56-63 | linear, circular (optional) |
| Division | 65-67 | 3-letter GenBank division code |
| Date | 69-79 | DD-MMM-YYYY |

### GenBank Divisions

| Code | Description |
|------|-------------|
| PRI | Primate |
| ROD | Rodent |
| MAM | Other mammalian |
| VRT | Other vertebrate |
| INV | Invertebrate |
| PLN | Plant, fungal, algal |
| BCT | Bacterial |
| VRL | Viral |
| PHG | Bacteriophage |
| SYN | Synthetic |
| UNA | Unannotated |
| EST | Expressed Sequence Tag |
| PAT | Patent |
| STS | Sequence Tagged Site |
| GSS | Genome Survey Sequence |
| HTG | High-Throughput Genomic |
| HTC | High-Throughput cDNA |
| ENV | Environmental |

---

## Feature Locations

Feature locations use INSDC syntax:

| Syntax | Example | Description |
|--------|---------|-------------|
| n..m | 100..200 | Range from n to m |
| n | 42 | Single position |
| <n..m | <1..206 | Partial at 5' end |
| n..>m | 500..>600 | Partial at 3' end |
| complement() | complement(100..200) | Minus strand |
| join() | join(1..50,60..100) | Discontinuous regions |
| complement(join()) | complement(join(a..b,c..d)) | Minus strand discontinuous |

---

## Implementation: GenBankParser

### Namespace
`Seqeron.Genomics.IO`

### Public API

#### Records

```csharp
public readonly record struct GenBankRecord(
    string Locus,
    int SequenceLength,
    string MoleculeType,
    string Topology,
    string Division,
    DateTime? Date,
    string Definition,
    string Accession,
    string Version,
    IReadOnlyList<string> Keywords,
    string Organism,
    string Taxonomy,
    IReadOnlyList<Reference> References,
    IReadOnlyList<Feature> Features,
    string Sequence,
    IReadOnlyDictionary<string, string> AdditionalFields);

public readonly record struct Reference(
    int Number,
    string Authors,
    string Title,
    string Journal,
    string PubMed,
    int? BaseFrom,
    int? BaseTo);

public readonly record struct Feature(
    string Key,
    Location Location,
    IReadOnlyDictionary<string, string> Qualifiers);

public readonly record struct Location(
    int Start,
    int End,
    bool IsComplement,
    bool IsJoin,
    IReadOnlyList<(int Start, int End)> Parts,
    string RawLocation);
```

#### Core Methods

| Method | Description |
|--------|-------------|
| `Parse(string content)` | Parses GenBank text into records |
| `ParseFile(string filePath)` | Reads and parses GenBank file |
| `ParseLocation(string locationStr)` | Parses feature location string |

#### Utility Methods

| Method | Description |
|--------|-------------|
| `GetFeatures(record, featureKey)` | Extracts features by type |
| `GetCDS(record)` | Gets CDS features |
| `GetGenes(record)` | Gets gene features |
| `ExtractSequence(record, location)` | Extracts subsequence for location |
| `GetQualifier(feature, qualifierName)` | Gets qualifier value |
| `TranslateCDS(record, cds)` | Translates CDS to protein |

---

## Complexity

| Operation | Time Complexity |
|-----------|-----------------|
| Parse | O(n) where n = content length |
| ParseFile | O(n) |
| ParseLocation | O(k) where k = location string length |
| ExtractSequence | O(m) where m = extracted length |

---

## Edge Cases

### Handled by Implementation

1. **Empty/null content** → Returns empty enumerable
2. **Missing file** → Returns empty enumerable
3. **Multiple records** → Split by "//" delimiter
4. **Empty KEYWORDS** → Single "." returns empty list
5. **Multi-line fields** → Continuation lines merged
6. **Complement locations** → Sets IsComplement flag
7. **Join locations** → Populates Parts list
8. **Sequence normalization** → Uppercase, removes numbers/spaces

### Known Limitations

1. Does not validate sequence against declared length
2. Does not support order() location specifier
3. Does not preserve original sequence case

---

## References

- Benson DA et al. (2008). "GenBank". Nucleic Acids Research 36: D25-D30.
- NCBI GenBank: https://www.ncbi.nlm.nih.gov/genbank/
- INSDC Feature Table: https://www.insdc.org/documents/feature_table.html
