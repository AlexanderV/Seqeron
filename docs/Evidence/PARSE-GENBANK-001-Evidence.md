# PARSE-GENBANK-001 Evidence

## Test Unit
**ID:** PARSE-GENBANK-001  
**Title:** GenBank Flat File Parsing  
**Area:** FileIO  
**Status:** Research Complete

---

## Authoritative Sources

### Primary Sources

1. **NCBI GenBank Sample Record Documentation**
   - URL: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
   - Accessed: 2026-02-05
   - Provides official field-by-field explanation of GenBank flat file format
   - Contains canonical example record (accession U49845)

2. **Wikipedia: GenBank**
   - URL: https://en.wikipedia.org/wiki/GenBank
   - Accessed: 2026-02-05
   - Describes GenBank as NCBI's open access nucleotide sequence database
   - Notes: Database started 1982, doubles every ~18 months, 34+ trillion base pairs (Oct 2024)

3. **INSDC Feature Table Definition**
   - URL: https://www.insdc.org/documents/feature_table.html
   - Referenced by NCBI documentation
   - Defines all valid feature keys and qualifiers

---

## GenBank Flat File Format Specification

### Record Structure (from NCBI)
A GenBank record consists of multiple sections, each identified by a keyword:

| Section | Description | Required |
|---------|-------------|----------|
| LOCUS | Name, length, molecule type, topology, division, date | Yes |
| DEFINITION | Brief description of sequence | Yes |
| ACCESSION | Unique stable identifier | Yes |
| VERSION | Accession.version format (e.g., U49845.1) | Yes |
| KEYWORDS | Descriptive words (often empty: ".") | No |
| SOURCE | Organism common/abbreviated name | Yes |
| ORGANISM | Scientific name + taxonomy lineage | Yes (within SOURCE) |
| REFERENCE | Literature citations | No |
| FEATURES | Biological annotations with locations | No |
| ORIGIN | Sequence data section marker | Yes |
| // | Record terminator | Yes |

### LOCUS Line Format (from NCBI)
```
LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
```
Fields (whitespace-delimited):
1. "LOCUS" keyword
2. Locus name (max 16 chars, unique identifier)
3. Sequence length + "bp"/"aa"
4. Molecule type: DNA, RNA, mRNA, etc.
5. Topology (optional): linear, circular
6. GenBank division: PRI, ROD, MAM, VRT, INV, PLN, BCT, VRL, PHG, SYN, etc.
7. Modification date: DD-MMM-YYYY

### GenBank Divisions (18 total, from NCBI)
| Code | Description |
|------|-------------|
| PRI | Primate sequences |
| ROD | Rodent sequences |
| MAM | Other mammalian sequences |
| VRT | Other vertebrate sequences |
| INV | Invertebrate sequences |
| PLN | Plant, fungal, algal sequences |
| BCT | Bacterial sequences |
| VRL | Viral sequences |
| PHG | Bacteriophage sequences |
| SYN | Synthetic sequences |
| UNA | Unannotated sequences |
| EST | Expressed Sequence Tags |
| PAT | Patent sequences |
| STS | Sequence Tagged Sites |
| GSS | Genome Survey Sequences |
| HTG | High-Throughput Genomic |
| HTC | High-Throughput cDNA |
| ENV | Environmental sampling |

### Feature Location Syntax (from INSDC/NCBI)
| Format | Example | Meaning |
|--------|---------|---------|
| n..m | 100..200 | Range from position n to m |
| n | 42 | Single position |
| <n..m | <1..206 | Partial on 5' end |
| n..>m | 4821..>5028 | Partial on 3' end |
| complement(n..m) | complement(3300..4037) | Minus strand |
| join(ranges) | join(1..50,60..100) | Discontinuous feature |
| complement(join(...)) | complement(join(1..50,60..100)) | Combined |

### Sequence Origin Section (from NCBI)
- Starts with "ORIGIN" keyword
- Each line: position number + 60 nucleotides (6 groups of 10)
- Sequence is lowercase in standard format
- Numbers, spaces, and newlines should be stripped during parsing
- Uppercase normalization is common practice

---

## Test Datasets from Sources

### NCBI Sample Record U49845
```
LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds, and Axl2p
            (AXL2) and Rev7p (REV7) genes, complete cds.
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;
            Saccharomycetales; Saccharomycetaceae; Saccharomyces.
```

### Documented Edge Cases

1. **Empty KEYWORDS** - Single period "." indicates no keywords
2. **Multi-line fields** - DEFINITION, ORGANISM can span multiple lines with indentation
3. **Multiple records** - Files can contain multiple records separated by "//"
4. **Complement features** - Located on minus strand
5. **Join locations** - Features spanning non-contiguous regions
6. **Partial features** - Using < and > for incomplete sequences
7. **Circular topology** - For plasmids, chromosomes, etc.

---

## Known Corner Cases / Failure Modes

### From NCBI Documentation
1. Locus names are no longer guaranteed to follow historical naming conventions
2. Dates may be in different formats (DD-MMM-YYYY or DD-MMM-YY)
3. Keywords field often empty (just ".")
4. Taxonomy lineage may be abbreviated for long lineages
5. Some older records may lack certain fields

### From Implementation Analysis
1. Empty content → return empty enumerable
2. Null content → return empty enumerable (defensive)
3. Missing LOCUS line → skip record
4. Malformed location strings → should handle gracefully
5. Missing ORIGIN → empty sequence
6. Non-alphabetic characters in sequence → should be filtered

---

## Invariants

1. **Sequence length** - Declared length in LOCUS should match actual sequence length
2. **Record termination** - Every record ends with "//"
3. **Location bounds** - Start ≤ End (except for single positions where Start = End)
4. **Feature extraction** - Extracted sequence for feature should match location bounds
5. **Complement consistency** - Complement flag should be detected for all complement() locations

---

## Testing Methodology

### Recommended Test Categories (Evidence-Based)

1. **Header Parsing** - LOCUS line fields (name, length, type, topology, division, date)
2. **Metadata Extraction** - DEFINITION, ACCESSION, VERSION, KEYWORDS
3. **Organism/Taxonomy** - SOURCE and ORGANISM parsing
4. **Reference Parsing** - AUTHORS, TITLE, JOURNAL, PUBMED
5. **Feature Parsing** - Feature keys, locations, qualifiers
6. **Location Parsing** - Simple ranges, complement, join, partial
7. **Sequence Extraction** - ORIGIN parsing, normalization
8. **Multi-record Handling** - Multiple records in single file
9. **Edge Cases** - Empty/null content, missing fields, malformed data

---

## Sources Summary

| Source | Type | Reliability |
|--------|------|-------------|
| NCBI Sample Record | Official Documentation | Authoritative |
| Wikipedia GenBank | Encyclopedia | High (references NCBI) |
| INSDC Feature Table | Official Specification | Authoritative |

---

## Recommendations for Test Implementation

1. Use NCBI-documented field formats as ground truth
2. Test all 18 division codes
3. Test all location syntax variants (complement, join, partial)
4. Verify sequence length matches LOCUS declaration
5. Test multi-record parsing
6. Test defensive handling of malformed input
