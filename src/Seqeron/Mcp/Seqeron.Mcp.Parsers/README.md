# Seqeron.Mcp.Parsers

MCP server — **FASTA/FASTQ/GenBank/EMBL/GFF/BED/VCF parsing, conversion and quality utilities.**

Exposes **41 tools** wrapping the `Seqeron.Genomics` library. Every tool has an
explicit JSON input/output schema, a Schema+Binding test, and per-tool docs under
[`docs/mcp/tools/parsers/`](../../../../docs/mcp/tools/parsers) — see the
campaign ledger [`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Parsers
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Parsers"]`). See [`docs/mcp/README.md`](../../../../docs/mcp/README.md).

## Tools (41)

| Tool | Description |
|------|-------------|
| `bed_filter` | Filter BED records by chromosome, region, strand, length, or score. |
| `bed_intersect` | Find intersecting regions between two BED datasets. |
| `bed_merge` | Merge overlapping BED records into single intervals. |
| `bed_parse` | Parse BED format content into genomic region records. |
| `embl_features` | Extract features from EMBL records by feature type (gene, CDS, mRNA, etc.). |
| `embl_parse` | Parse EMBL flat file format into structured records. |
| `embl_statistics` | Calculate statistics for EMBL records. |
| `fasta_format` | Format sequence(s) to FASTA string. |
| `fasta_parse` | Parse FASTA format string into sequence entries. |
| `fasta_write` | Write sequence(s) to a FASTA file. |
| `fastq_detect_encoding` | Detect quality score encoding from FASTQ quality string. |
| `fastq_encode_quality` | Encode Phred quality scores to a quality string. |
| `fastq_error_to_phred` | Convert error probability to Phred quality score. |
| `fastq_filter` | Filter FASTQ reads by minimum average quality score. |
| `fastq_format` | Format a single FASTQ record to string format. |
| `fastq_parse` | Parse FASTQ format string into sequence entries with quality scores. |
| `fastq_phred_to_error` | Convert Phred quality score to error probability. |
| `fastq_statistics` | Calculate quality statistics for FASTQ data. |
| `fastq_trim_adapter` | Trim adapter sequences from FASTQ reads. |
| `fastq_trim_quality` | Trim low-quality bases from both ends of FASTQ reads. |
| `fastq_write` | Write FASTQ records to a file. |
| `genbank_extract_sequence` | Extract a subsequence from a GenBank record based on a feature location string. |
| `genbank_features` | Extract features from GenBank records by feature type (gene, CDS, mRNA, etc.). |
| `genbank_parse` | Parse GenBank flat file format into structured records. |
| `genbank_parse_location` | Parse a GenBank feature location string into its components. |
| `genbank_statistics` | Calculate statistics for GenBank records. |
| `gff_filter` | Filter GFF/GTF records by feature type, sequence ID, or genomic region. |
| `gff_parse` | Parse GFF3/GTF format content into feature records. |
| `gff_statistics` | Calculate statistics for GFF/GTF annotations. |
| `vcf_classify` | Classify variant type for a VCF record. |
| `vcf_filter` | Filter VCF variants by type, quality, chromosome, or PASS status. |
| `vcf_has_flag` | Check if a VCF INFO field flag is present. |
| `vcf_is_het` | Check if a genotype is heterozygous (different alleles). |
| `vcf_is_hom_alt` | Check if a genotype is homozygous alternate (e.g., 1/1, 2/2). |
| `vcf_is_hom_ref` | Check if a genotype is homozygous reference (0/0 or 0|0). |
| `vcf_is_indel` | Check if a variant is an indel (insertion or deletion). |
| `vcf_is_snp` | Check if a variant is a SNP (Single Nucleotide Polymorphism). |
| `vcf_parse` | Parse VCF (Variant Call Format) content into variant records. |
| `vcf_statistics` | Calculate statistics for VCF variants. |
| `vcf_variant_length` | Get the length difference of a variant (absolute difference between ref and alt lengths). |
| `vcf_write` | Write VCF records to a file. |
