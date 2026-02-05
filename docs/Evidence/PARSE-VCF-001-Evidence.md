# Evidence: PARSE-VCF-001 - VCF Parsing

## Test Unit
**ID:** PARSE-VCF-001  
**Area:** FileIO  
**Algorithm:** VCF (Variant Call Format) Parsing

---

## Primary Sources

### 1. Wikipedia - Variant Call Format
**URL:** https://en.wikipedia.org/wiki/Variant_Call_Format  
**Accessed:** 2026-02-05

#### Key Information Extracted:

**VCF Format Overview:**
- VCF is a standard text file format used in bioinformatics for storing gene sequence variations
- Developed in 2010 for the 1000 Genomes Project
- Current version: VCFv4.5 (October 2024)
- Tab-delimited format with 8 mandatory columns

**Column Specification (1-based positions):**
| Column | Name | Description |
|--------|------|-------------|
| 1 | CHROM | Chromosome/sequence name |
| 2 | POS | 1-based position on the sequence |
| 3 | ID | Variant identifier (e.g., rsID) or "." if unknown |
| 4 | REF | Reference base(s) at this position |
| 5 | ALT | List of alternative alleles (comma-separated) |
| 6 | QUAL | Quality score (Phred-scaled) or "." if missing |
| 7 | FILTER | PASS if passed all filters, or list of failed filters (semicolon-separated) |
| 8 | INFO | Key-value pairs (semicolon-separated) |
| 9 | FORMAT | Sample field format specification (optional) |
| 10+ | SAMPLE | Sample genotype data (optional) |

**Header Lines:**
- Lines starting with `##` are metadata lines
- `##fileformat=VCFvX.X` specifies the VCF version
- `##INFO=<...>` defines INFO field meanings
- `##FORMAT=<...>` defines FORMAT field meanings
- `##FILTER=<...>` defines FILTER meanings
- `#CHROM` line is the column header line

**Common INFO Fields:**
| Field | Description |
|-------|-------------|
| DP | Total depth across samples |
| AF | Allele frequency |
| NS | Number of samples with data |
| DB | dbSNP membership (flag) |
| AA | Ancestral allele |

**Common FORMAT Fields:**
| Field | Description |
|-------|-------------|
| GT | Genotype |
| DP | Read depth |
| AD | Allele depth |
| GQ | Genotype quality |
| PL | Phred-scaled genotype likelihoods |

**Variant Classification Rules:**
- SNP: ref and alt are both single nucleotides with different values
- MNP: ref and alt are equal length (>1) with different values
- Insertion: alt longer than ref (commonly ref is 1 base)
- Deletion: ref longer than alt (commonly alt is 1 base)
- Symbolic: alt starts with `<` (e.g., `<DEL>`, `<INS>`, `<DUP>`)

**Genotype Notation:**
- `0/0`: Homozygous reference
- `0/1` or `1/0`: Heterozygous
- `1/1`: Homozygous alternate
- `|` indicates phased genotypes
- `/` indicates unphased genotypes
- `.` indicates missing data

---

### 2. Danecek et al. (2011) - The Variant Call Format and VCFtools
**Reference:** Danecek P, et al. "The variant call format and VCFtools." Bioinformatics. 2011;27(15):2156-2158.  
**DOI:** 10.1093/bioinformatics/btr330  
**PMID:** 21653522

#### Key Information:
- Authoritative paper defining the VCF format
- Establishes parsing and validation conventions
- Defines expected behavior for missing values (`.`)

---

### 3. SAMtools HTS-Specs VCF Specification
**URL:** https://samtools.github.io/hts-specs/VCFv4.3.pdf  
**Version:** VCF 4.3

#### Key Information:
- Official specification maintained by GA4GH
- Defines strict parsing rules
- Specifies error handling for malformed records

---

### 4. EMBL-EBI VCF Training
**URL:** https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/understanding-vcf-format/  
**Accessed:** 2026-02-05

#### Key Information:
- VCF is scalable and flexible
- Standard input for variant analysis tools like VEP
- Standard output from variant callers like GATK

---

## Edge Cases (Source-Derived)

### From VCF Specification (Wikipedia/Danecek 2011):

1. **Missing Values:**
   - Quality score "." → null/missing
   - ID "." → no identifier
   - Filter "." → no filters applied (equivalent to PASS in some contexts)
   - INFO "." → no info fields

2. **Multi-allelic Variants:**
   - ALT field can contain multiple alleles: `A,C,G`
   - Each represents a different alternate allele at the same position

3. **Complex Filters:**
   - Multiple failed filters: `LowQual;LowCov`
   - Single passing: `PASS`

4. **Genotype Edge Cases:**
   - Phased vs unphased: `0|1` vs `0/1`
   - Missing allele: `./0`, `0/.`, `./.`
   - Multi-allelic reference: `1/2` (both alleles are alternates)

5. **Symbolic Alleles:**
   - Structural variants: `<DEL>`, `<INS>`, `<DUP>`, `<INV>`, `<CNV>`
   - Breakend notation: `]chr:pos]` or `[chr:pos[`

6. **Malformed Records:**
   - Invalid position (non-integer) → skip record
   - Fewer than 8 columns → skip record

---

## Test Datasets

### Dataset 1: Basic VCF with Samples (Wikipedia Example)
```
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	100	rs123	A	G	99	PASS	DP=50	GT	0/1
chr1	200	.	C	T	50	.	DP=30	GT	0/0
chr2	300	rs456	G	A,C	80	PASS	DP=40	GT	1/2
```

### Dataset 2: Variant Types (Derived from Specification)
```
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	G	99	PASS	.        # SNP
chr1	200	.	AT	GC	99	PASS	.       # MNP
chr1	300	.	A	ATG	99	PASS	.      # Insertion
chr1	400	.	ATG	A	99	PASS	.      # Deletion
chr1	500	.	A	<DEL>	99	PASS	.    # Symbolic
```

### Dataset 3: Ti/Tv Transitions and Transversions
```
Transitions (purines A↔G, pyrimidines C↔T):
- A→G, G→A, C→T, T→C

Transversions (purine↔pyrimidine):
- A→C, A→T, G→C, G→T, C→A, C→G, T→A, T→G
```

---

## Testing Methodology

### Required Test Categories:

1. **Parsing Tests:**
   - Empty/null content → empty result
   - Valid VCF → correct records
   - Header extraction → correct metadata
   - Sample parsing → correct genotypes

2. **Variant Classification Tests:**
   - SNP detection (single-base substitution)
   - MNP detection (multi-nucleotide polymorphism)
   - Insertion detection (alt > ref)
   - Deletion detection (ref > alt)
   - Symbolic allele detection

3. **Filter Tests:**
   - Filter by chromosome
   - Filter by region (position range)
   - Filter by quality threshold
   - Filter PASS only
   - Filter by variant type

4. **Genotype Analysis Tests:**
   - Homozygous reference detection (0/0)
   - Homozygous alternate detection (1/1, 2/2)
   - Heterozygous detection (0/1, 1/2)
   - Phased vs unphased handling

5. **Statistics Tests:**
   - Variant counts by type
   - Chromosome distribution
   - Ti/Tv ratio calculation
   - Passing variant counts

6. **Edge Case Tests:**
   - Missing quality values
   - Missing IDs
   - Multi-allelic sites
   - Multiple filter values
   - Malformed lines (skip gracefully)

7. **I/O Tests:**
   - File parsing
   - Round-trip (write and re-read)

---

## Invariants (Source-Derived)

1. **Position Invariant:** POS is always 1-based positive integer
2. **Quality Invariant:** QUAL is non-negative when present, or null when "."
3. **Filter Invariant:** PASS indicates no failed filters
4. **Variant Length Invariant:** |len(ALT) - len(REF)| gives indel length
5. **Ti/Tv Ratio Invariant:** Typically 2.0-2.1 for whole genome, ~3.0 for exome in humans

---

## Test Coverage

### Reference Data Tests (Added 2026-02-05)
Tests validated against real-world VCF sources:

**1000 Genomes Project:**
- `Parse_1000GenomesFormat_CommonSnps_ParsedCorrectly` - validates common SNP format with real rsIDs
- `Parse_1000GenomesFormat_MultiAllelicSite` - validates multi-allelic variant parsing
- `Parse_1000GenomesFormat_InfoFields` - validates AF, DP, NS info field parsing

**ClinVar Database:**
- `Parse_ClinVarFormat_PathogenicVariants_ParsedCorrectly` - validates pathogenic variant format
- `Parse_ClinVarFormat_InfoFieldExtraction` - validates CLNSIG, CLNDN fields
- `Parse_ClinVarFormat_StructuralVariants` - validates symbolic allele parsing (<DEL>, <INS>)

---

## References Summary

| Source | Type | Key Contribution |
|--------|------|------------------|
| Wikipedia VCF | Encyclopedia | Format overview, column specs, common fields |
| Danecek et al. (2011) | Peer-reviewed | Authoritative format definition |
| SAMtools HTS-specs | Specification | Official format specification |
| EMBL-EBI Training | Educational | Usage context, tool integration |
| 1000 Genomes Project | Reference Data | Real-world SNP validation |
| ClinVar Database | Reference Data | Pathogenic variant validation |

---

## Notes

- VCF is the de facto standard for variant data exchange
- Implementation follows VCFv4.3 specification
- Parser gracefully handles malformed data by skipping invalid lines

---

## Change History
- **2026-02-05**: Added reference data tests from 1000 Genomes Project and ClinVar.
- **2026-02-05**: Initial documentation.
