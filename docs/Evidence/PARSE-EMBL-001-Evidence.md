# Evidence Document: PARSE-EMBL-001 (EMBL Parsing)

## Test Unit Information
- **Test Unit ID**: PARSE-EMBL-001
- **Algorithm Name**: EMBL Parsing
- **Category**: FileIO/Parsing
- **Implementation**: `EmblParser.Parse(string content)`, `EmblParser.ParseFile(string filePath)`
- **Location**: `src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs`

## Authoritative Sources

### Primary Source: EBI EMBL User Manual (Release 143, March 2020)
- **URL**: `https://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt`
- **Publisher**: European Bioinformatics Institute (EBI), European Nucleotide Archive (ENA)
- **Version**: Release 143, March 2020
- **Authority**: Official EMBL database format specification

### Secondary Source: INSDC Feature Table Definition v11.3
- **URL**: `https://www.insdc.org/files/feature_table.html`
- **Publisher**: International Nucleotide Sequence Database Collaboration (INSDC)
- **Version**: v11.3 (October 2024)
- **Authority**: Official DDBJ/ENA/GenBank feature table specification

## EMBL Format Specification

### Line Type Codes
Each line begins with a two-character line type code followed by three blanks. The actual information begins at character position 6.

| Code | Description | Occurrence |
|------|-------------|------------|
| ID | Identification | 1 per entry (first line) |
| AC | Accession number | ≥1 per entry |
| PR | Project identifier | 0 or 1 per entry |
| DT | Date | 2 per entry |
| DE | Description | ≥1 per entry |
| KW | Keyword | ≥1 per entry |
| OS | Organism species | ≥1 per entry |
| OC | Organism classification | ≥1 per entry |
| OG | Organelle | 0 or 1 per entry |
| RN | Reference number | ≥1 per entry |
| RC | Reference comment | ≥0 per entry |
| RP | Reference positions | ≥1 per entry |
| RX | Reference cross-reference | ≥0 per entry |
| RG | Reference group | ≥0 per entry |
| RA | Reference author(s) | ≥0 per entry |
| RT | Reference title | ≥1 per entry |
| RL | Reference location | ≥1 per entry |
| DR | Database cross-reference | ≥0 per entry |
| CC | Comments or notes | ≥0 per entry |
| AH | Assembly header | 0 or 1 per entry (TPA/TSA) |
| AS | Assembly information | 0 or ≥1 per entry (TPA/TSA) |
| FH | Feature table header | 2 per entry |
| FT | Feature table data | ≥2 per entry |
| CO | Contig/construct line | 0 or ≥1 per entry (CON) |
| XX | Spacer line | Many per entry |
| SQ | Sequence header | 1 per entry |
| bb | Sequence data | ≥1 per entry (blanks) |
| // | Termination line | 1 per entry (last line) |

### ID Line Format (Section 3.4.1)
```
ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
```

Token positions:
1. **Primary accession number**
2. **Sequence version number** (SV)
3. **Topology**: 'circular' or 'linear'
4. **Molecule type**: Same as mol_type qualifier value
5. **Data class**: CON, PAT, EST, GSS, HTC, HTG, WGS, TSA, STS, STD
6. **Taxonomic division**: PHG, ENV, FUN, HUM, INV, MAM, VRT, MUS, PLN, PRO, ROD, SYN, TGN, UNC, VRL
7. **Sequence length** in base pairs

**Example**:
```
ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
```

### AC Line Format (Section 3.4.2)
Semicolon-separated accession numbers. Primary accession first, secondary accessions follow.
```
AC   X56734; S46826;
AC   Y00001; X00001-X00005; X00008; Z00001-Z00005;
```

### DT Line Format (Section 3.4.4)
Two DT lines per entry:
```
DT   DD-MON-YYYY (Rel. #, Created)
DT   DD-MON-YYYY (Rel. #, Last updated, Version #)
```

### DE Line Format (Section 3.4.5)
Free-format description in ordinary English. May span multiple lines.
```
DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
```

### KW Line Format (Section 3.4.6)
Semicolon-separated keywords, terminated with period.
```
KW   beta-glucosidase.
KW   .   (empty keywords)
```

### OS Line Format (Section 3.4.7)
Organism species: Latin genus and species with common name in parentheses.
```
OS   Trifolium repens (white clover)
```

### OC Line Format (Section 3.4.8)
Taxonomic classification from general to specific, semicolon-separated, period-terminated.
```
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   Spermatophyta; Magnoliophyta; eudicotyledons; ...
```

### OG Line Format (Section 3.4.9)
Organelle for non-nuclear sequences.
```
OG   Plastid:Chloroplast
OG   Mitochondrion
OG   Plasmid pBR322
```

### Reference Block (Section 3.4.10)
Order: RN, RC, RP, RX, RG, RA, RT, RL
```
RN   [5]
RP   1-1859
RX   DOI; 10.1007/BF00039495.
RX   PUBMED; 1907511.
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
RT   "Nucleotide and derived amino acid sequence...";
RL   Plant Mol. Biol. 17(2):209-219(1991).
```

### FH/FT Lines (Section 3.4.15/3.4.16)
Feature header is fixed format for readability:
```
FH   Key             Location/Qualifiers
FH
```

Feature table data:
```
FT   source          1..1859
FT                   /organism="Trifolium repens"
FT                   /mol_type="mRNA"
FT   CDS             14..1495
FT                   /product="beta-glucosidase"
FT                   /protein_id="CAA40058.1"
```

### SQ Line Format (Section 3.4.17)
```
SQ   Sequence <length> BP; <A count> A; <C count> C; <G count> G; <T count> T; <other count> other;
```

**Example**:
```
SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
```

### Sequence Data Lines (Section 3.4.18)
- Line code is two blanks
- 60 bases per line in groups of 10 separated by spaces
- Sequence begins at position 6
- Direction: 5' to 3'
- Columns 73-80 contain right-justified base position numbers

**Example**:
```
     aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
```

### Terminator Line (Section 3.4.21)
```
//
```

## Feature Location Descriptors (INSDC Feature Table)

### Simple Locations
- `100..200` - Range from base 100 to 200
- `467` - Single base at position 467
- `123^124` - Site between bases 123 and 124

### Partial Locations
- `<1..200` - 5' partial (extends beyond known sequence)
- `100..>500` - 3' partial (extends beyond known sequence)
- `<1..>500` - Both ends partial

### Complement
- `complement(100..200)` - Complementary strand

### Join Operator
- `join(100..200,300..400)` - Discontinuous spans on same strand
- `join(complement(100..200),complement(300..400))` - Complementary discontinuous

### Order Operator (alternative to join)
- `order(100..200,300..400)` - Non-contiguous spans, order significant

### Remote References
- `J00194.1:100..200` - Location on another entry

## Data Classes (Section 3.1)

| Class | Definition |
|-------|------------|
| CON | Entry constructed from segment entry sequences |
| PAT | Patent |
| EST | Expressed Sequence Tag |
| GSS | Genome Survey Sequence |
| HTC | High Throughput CDNA sequencing |
| HTG | High Throughput Genome sequencing |
| WGS | Whole Genome Shotgun |
| TSA | Transcriptome Shotgun Assembly |
| STS | Sequence Tagged Site |
| STD | Standard (all entries not classified above) |

## Taxonomic Divisions (Section 3.2)

| Code | Division |
|------|----------|
| PHG | Bacteriophage |
| ENV | Environmental Sample |
| FUN | Fungal |
| HUM | Human |
| INV | Invertebrate |
| MAM | Other Mammal |
| VRT | Other Vertebrate |
| MUS | Mus musculus |
| PLN | Plant |
| PRO | Prokaryote |
| ROD | Other Rodent |
| SYN | Synthetic |
| TGN | Transgenic |
| UNC | Unclassified |
| VRL | Viral |

## Standard Base Codes (Appendix A)

| Code | Base | Description |
|------|------|-------------|
| G | Guanine | |
| A | Adenine | |
| T | Thymine | |
| C | Cytosine | |
| R | Purine | A or G |
| Y | Pyrimidine | C or T or U |
| M | Amino | A or C |
| K | Ketone | G or T |
| S | Strong | C or G |
| W | Weak | A or T |
| H | Not-G | A or C or T |
| B | Not-A | C or G or T |
| V | Not-T | A or C or G |
| D | Not-C | A or G or T |
| N | Any | A or C or G or T |

## Edge Cases and Special Handling

### Empty/Null Content
- Null content: Should throw or return empty result
- Empty string: Should return empty collection
- Whitespace only: Should return empty collection

### Minimal Valid Record
Minimum required lines:
- ID line (mandatory first line)
- AC line (≥1)
- DT lines (2)
- DE line (≥1)
- KW line (≥1, can be "KW   .")
- OS line (≥1)
- OC line (≥1)
- FT source feature (mandatory per INSDC)
- SQ line (1)
- Sequence data (≥1)
- // terminator (mandatory last line)

### Multiple Records
Multiple records separated by // terminators in single content string.

### Circular vs Linear Topology
ID line position 3: "circular" or "linear"

### Case Handling
- Sequence data: Stored as lowercase in EMBL, should normalize
- Line codes: Uppercase

### Line Continuation
- DE, KW, OS, OC, RA, RT, RL lines may span multiple lines
- FT qualifier values may span multiple lines (continuation starts at column 22)

## Test Categories

### Must Test (Critical)
1. Parse valid EMBL record - correct field extraction
2. Parse multiple records from single content
3. ID line parsing (accession, topology, molecule type, length)
4. Feature extraction with locations
5. Sequence extraction and normalization
6. Empty/null content handling
7. Record terminator (//) detection

### Should Test (Important)
1. Reference block parsing (RN, RA, RT, RL)
2. All metadata fields (DE, KW, OS, OC, OG)
3. Feature qualifier extraction
4. Location parsing (simple range, complement, join)
5. Partial locations (<, >)
6. Circular topology handling
7. ParseFile method

### Could Test (Nice to Have)
1. All data classes (CON, PAT, EST, etc.)
2. All taxonomic divisions
3. Remote location references
4. DR database cross-references
5. CC comments
6. ToGenBank conversion
7. Utility methods (GetCDS, GetGenes)

## Reference EMBL Record (from EBI Manual)

```
ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
XX
AC   X56734; S46826;
XX
DT   12-SEP-1991 (Rel. 29, Created)
DT   25-NOV-2005 (Rel. 85, Last updated, Version 11)
XX
DE   Trifolium repens mRNA for non-cyanogenic beta-glucosidase
XX
KW   beta-glucosidase.
XX
OS   Trifolium repens (white clover)
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons; rosids;
OC   fabids; Fabales; Fabaceae; Papilionoideae; Trifolieae; Trifolium.
XX
RN   [5]
RP   1-1859
RX   DOI; 10.1007/BF00039495.
RX   PUBMED; 1907511.
RA   Oxtoby E., Dunn M.A., Pancoro A., Hughes M.A.;
RT   "Nucleotide and derived amino acid sequence of the cyanogenic
RT   beta-glucosidase (linamarase) from white clover (Trifolium repens L.)";
RL   Plant Mol. Biol. 17(2):209-219(1991).
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..1859
FT                   /organism="Trifolium repens"
FT                   /mol_type="mRNA"
FT                   /clone_lib="lambda gt10"
FT                   /clone="TRE361"
FT                   /tissue_type="leaves"
FT                   /db_xref="taxon:3899"
FT   mRNA            1..1859
FT                   /experiment="experimental evidence, no additional details
FT                   recorded"
FT   CDS             14..1495
FT                   /product="beta-glucosidase"
FT                   /EC_number="3.2.1.21"
FT                   /note="non-cyanogenic"
FT                   /db_xref="GOA:P26204"
FT                   /db_xref="InterPro:IPR001360"
FT                   /db_xref="UniProtKB/Swiss-Prot:P26204"
FT                   /protein_id="CAA40058.1"
FT                   /translation="MDFIVAIFALFVISSFTITSTNAVEASTLLDIGNLSRSSFPRGFI
FT                   FGAGSSAYQFEGAVNEGGRGPSIWDTFTHKYPEKIRDGSNADITVDQYHRYKEDVGIMK..."
XX
SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
     aaacaaacca aatatggatt ttattgtagc catatttgct ctgtttgtta ttagctcatt        60
     cacaattact tccacaaatg cagttgaagc ttctactctt cttgacatag gtaacctgag       120
     ...
//
```

## Document History
- **Created**: 2025-01-28
- **Author**: Algorithm Testing Protocol
- **Sources Verified**: EBI User Manual (Release 143), INSDC Feature Table v11.3
