# Evidence: PARSE-GFF-001 - GFF/GTF Parsing

## Test Unit
**ID:** PARSE-GFF-001  
**Area:** FileIO  
**Methods:** Parse(content), ParseFile(filePath), ToGff3(features)

---

## Authoritative Sources

### 1. Wikipedia - General Feature Format
**URL:** https://en.wikipedia.org/wiki/General_feature_format  
**Retrieved:** 2026-02-05

#### Key Specifications:
- GFF is a tab-delimited file format with **9 columns**
- All GFF formats (GFF2, GFF3, GTF) share the same structure for the first 8 fields
- GFF2/GTF has been deprecated in favor of GFF3

**Column Definitions:**
| Position | Field | Description |
|----------|-------|-------------|
| 1 | seqid | Name of the sequence where the feature is located |
| 2 | source | Algorithm or procedure that generated the feature |
| 3 | type | Feature type (e.g., "gene", "exon", "CDS") |
| 4 | start | Genomic start (1-based offset) |
| 5 | end | Genomic end (1-based, inclusive) |
| 6 | score | Numeric confidence value ("." for undefined) |
| 7 | strand | "+", "-", "." (unstranded), or "?" (unknown) |
| 8 | phase | 0, 1, or 2 for CDS features; "." for others |
| 9 | attributes | List of tag-value pairs |

---

### 2. UCSC Genome Browser - GFF/GTF Format
**URL:** https://genome.ucsc.edu/FAQ/FAQformat.html  
**Retrieved:** 2026-02-05

#### GFF2 Field Descriptions:
1. **seqname** - Chromosome or scaffold name
2. **source** - Program that generated this feature
3. **feature** - Feature type (e.g., "CDS", "start_codon", "stop_codon", "exon")
4. **start** - Starting position (1-based)
5. **end** - Ending position (inclusive)
6. **score** - Value between 0 and 1000 ("." if not available)
7. **strand** - "+", "-", or "." (don't know/don't care)
8. **frame** - 0-2 for coding exons; "." otherwise
9. **group** - All lines with same group are linked together

#### GTF Format Extensions:
- GTF is an extension of GFF2, backward compatible
- Attributes must end in semicolon, separated by exactly one space
- **Mandatory attributes:**
  - `gene_id value` - Globally unique identifier for genomic source
  - `transcript_id value` - Globally unique identifier for predicted transcript

---

### 3. Sequence Ontology - GFF3 Specification v1.26
**URL:** https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md  
**Author:** Lincoln Stein  
**Version:** 1.26 (18 August 2020)

#### GFF3 Format Overview:
- Nine-column, tab-delimited, plain text files
- UTF-8 encoding recommended for portability
- Lines beginning with `##` are directives (pragmas/meta-data)
- Lines beginning with `#` are comments
- Blank lines are ignored

#### Escaping Rules (RFC 3986):
Must escape:
- Tab: `%09`
- Newline: `%0A`
- Carriage return: `%0D`
- Percent: `%25`
- Control characters: `%00` through `%1F`, `%7F`

Column 9 reserved characters:
- Semicolon: `%3B`
- Equals: `%3D`
- Ampersand: `%26`
- Comma: `%2C`

#### Predefined Attribute Tags:
| Tag | Description |
|-----|-------------|
| ID | Feature identifier (required for features with children) |
| Name | Display name for the feature |
| Alias | Secondary name for the feature |
| Parent | Parent feature ID (establishes part-of relationship) |
| Target | Target of nucleotide-to-nucleotide alignment |
| Gap | Alignment format (CIGAR-like) |
| Derives_from | Temporal relationship between features |
| Note | Free text note |
| Dbxref | Database cross reference |
| Ontology_term | Cross reference to ontology term |
| Is_circular | Flag for circular features |

#### Phase Field (Column 8):
- Required for all CDS features
- Values: 0, 1, or 2
- Indicates number of bases to remove from 5' end to reach next codon
- Phase 0: codon begins at first nucleotide
- Phase 1: codon begins at second nucleotide
- Phase 2: codon begins at third nucleotide

#### Standard Directives:
- `##gff-version 3.x.x` - Must be first line
- `##sequence-region seqid start end` - Bounds for a sequence
- `##FASTA` - Indicates FASTA sequences follow
- `###` - All forward references resolved

---

## Test Datasets from Sources

### Canonical Gene Example (from GFF3 Spec):
```gff3
##gff-version 3
##sequence-region ctg123 1 1497228
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123 . exon            1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001
ctg123 . CDS             3000  3902  .  +  0  ID=cds00001;Parent=mRNA00001
```

### GTF Example (from Ensembl/UCSC):
```gtf
chr1	ENSEMBL	gene	1000	5000	.	+	.	gene_id "ENSG00001"; gene_name "TestGene";
chr1	ENSEMBL	transcript	1000	5000	.	+	.	gene_id "ENSG00001"; transcript_id "ENST00001";
chr1	ENSEMBL	exon	1000	1500	.	+	.	gene_id "ENSG00001"; transcript_id "ENST00001"; exon_number "1";
```

---

## Documented Corner Cases

### From GFF3 Specification:
1. **Zero-length features** - Start equals end, site is to the right of indicated base
2. **Circular genomes** - Use `Is_circular=true` attribute
3. **Discontinuous features** - Same ID on multiple lines
4. **Multiple parents** - Comma-separated Parent values
5. **Polycistronic transcripts** - Multiple genes parent single mRNA
6. **Trans-spliced transcripts** - mRNA spans disjunct genomic locations

### From UCSC:
1. **Fields must be tab-separated** - Spaces within fields are allowed
2. **Score interpretation** - Value 0-1000 determines gray shading

---

## Known Edge Cases

1. **Empty content** → Return empty collection
2. **Null content** → Return empty collection
3. **Malformed lines** (< 8 fields) → Skip line
4. **Missing score** (".") → Score is null
5. **Missing phase** (".") → Phase is null
6. **Unknown strand** ("?") → Valid strand character
7. **URL-encoded characters** → Must unescape
8. **Comment lines** (# prefix) → Skip
9. **Directive lines** (## prefix) → Process as metadata
10. **Empty lines** → Skip

---

## Testing Methodology

### Required Test Categories (per GFF3 Spec):
1. **Column parsing** - All 9 columns correctly extracted
2. **Attribute parsing** - GFF3 (key=value;) vs GTF (key "value";)
3. **Escaping/Unescaping** - RFC 3986 percent-encoding
4. **Hierarchical relationships** - Parent/child feature linking
5. **Format detection** - Auto-detection of GFF3 vs GTF vs GFF2
6. **Directives** - Version, sequence-region processing
7. **Round-trip** - Write and re-parse preserves data

---

## Summary

| Source | Key Contribution |
|--------|------------------|
| Wikipedia | General format overview, 9-column structure |
| UCSC | GFF2/GTF field descriptions, GTF mandatory attributes |
| Sequence Ontology | Authoritative GFF3 spec v1.26, escaping rules, attribute semantics |

**Confidence Level:** HIGH - All specifications are well-documented standards.
