using Seqeron.Mcp.Sequence.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.Sequence server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(SequenceTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "amino_acid_composition",
        "complement_base",
        "complexity_compression_ratio",
        "complexity_dust_score",
        "complexity_kmer_entropy",
        "complexity_linguistic",
        "complexity_mask_low",
        "complexity_shannon",
        "dna_reverse_complement",
        "dna_validate",
        "gc_content",
        "hydrophobicity",
        "is_valid_dna",
        "is_valid_rna",
        "isoelectric_point",
        "iupac_code",
        "iupac_match",
        "iupac_matches",
        "kmer_analyze",
        "kmer_count",
        "kmer_distance",
        "kmer_entropy",
        "linguistic_complexity",
        "melting_temperature",
        "molecular_weight_nucleotide",
        "molecular_weight_protein",
        "nucleotide_composition",
        "protein_validate",
        "rna_from_dna",
        "rna_validate",
        "shannon_entropy",
        "summarize_sequence",
        "thermodynamics",
        "translate_dna",
        "translate_rna",
    };
}
