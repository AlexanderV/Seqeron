using Seqeron.Mcp.MolTools.Tools;
using Seqeron.Mcp.Tests.Contracts;
using NUnit.Framework;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Pins the Seqeron.Mcp.MolTools server's public MCP tool surface. Assertions live in
/// <see cref="ToolRegistrationContractTestsBase"/>; this fixture supplies the tool class and the
/// documented tool-name set.
/// </summary>
[TestFixture]
public sealed class ToolRegistrationContractTests : ToolRegistrationContractTestsBase
{
    protected override Type ToolsType => typeof(MolToolsTools);

    protected override string[] ExpectedToolNames { get; } =
    {
        "analyze_oligo",
        "blunt_cutters",
        "build_codon_table",
        "cai_from_organism_table",
        "codon_adaptation_index",
        "codon_usage_statistics",
        "compare_codon_usage",
        "compatible_enzymes",
        "count_codons",
        "crispr_specificity_score",
        "crispr_system_info",
        "design_antisense_probes",
        "design_guide_rnas",
        "design_molecular_beacon",
        "design_primers",
        "design_probes",
        "design_tiling_probes",
        "digest_summary",
        "effective_number_of_codons",
        "enzymes_by_cut_length",
        "enzymes_compatible",
        "evaluate_guide_rna",
        "evaluate_primer",
        "find_all_restriction_sites",
        "find_off_targets",
        "find_pam_sites",
        "find_rare_codons",
        "find_restriction_sites",
        "generate_primer_candidates",
        "get_enzyme",
        "hairpin_potential",
        "longest_dinucleotide_repeat",
        "longest_homopolymer",
        "oligo_concentration_from_absorbance",
        "oligo_extinction_coefficient",
        "optimize_codons",
        "primer_dimer",
        "primer_melting_temperature",
        "primer_melting_temperature_salt",
        "reduce_secondary_structure",
        "remove_restriction_sites",
        "restriction_digest",
        "restriction_map",
        "rscu",
        "sticky_cutters",
        "three_prime_stability",
        "validate_probe",
    };
}
