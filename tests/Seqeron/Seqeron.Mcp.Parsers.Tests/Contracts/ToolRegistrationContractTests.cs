using Seqeron.Mcp.Parsers.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Parsers.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Parsers server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(ParsersTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "bed_filter",
        "bed_intersect",
        "bed_merge",
        "bed_parse",
        "embl_features",
        "embl_parse",
        "embl_statistics",
        "fasta_format",
        "fasta_parse",
        "fasta_write",
        "fastq_detect_encoding",
        "fastq_encode_quality",
        "fastq_error_to_phred",
        "fastq_filter",
        "fastq_format",
        "fastq_parse",
        "fastq_phred_to_error",
        "fastq_statistics",
        "fastq_trim_adapter",
        "fastq_trim_quality",
        "fastq_write",
        "genbank_extract_sequence",
        "genbank_features",
        "genbank_parse",
        "genbank_parse_location",
        "genbank_statistics",
        "gff_filter",
        "gff_parse",
        "gff_statistics",
        "vcf_classify",
        "vcf_filter",
        "vcf_has_flag",
        "vcf_is_het",
        "vcf_is_hom_alt",
        "vcf_is_hom_ref",
        "vcf_is_indel",
        "vcf_is_snp",
        "vcf_parse",
        "vcf_statistics",
        "vcf_variant_length",
        "vcf_write",
    };
}
