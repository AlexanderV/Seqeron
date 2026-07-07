using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Server contract / integration test: the Seqeron.Mcp.MolTools server must advertise EXACTLY its documented
/// tool surface, and every registered tool must meet the gold-standard binding shape
/// (explicit Name + Title + ReadOnly + Description). This pins the server's public MCP surface
/// so an accidental rename, removal, or un-annotated tool fails the build.
/// </summary>
[TestFixture]
public class ToolRegistrationContractTests
{
    private static readonly string[] ExpectedToolNames =
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

    [Test]
    public void MolToolsTools_AdvertisesExactlyTheDocumentedToolSurface()
    {
        var registered = GetRegisteredToolNames(typeof(MolToolsTools));

        Assert.That(registered, Is.EquivalentTo(ExpectedToolNames),
            "the server's registered tool names must match the documented surface exactly");
        Assert.That(registered.Length, Is.EqualTo(ExpectedToolNames.Length),
            "tool count must equal the documented surface");
        Assert.That(registered.Distinct().Count(), Is.EqualTo(registered.Length),
            "tool names must be unique");
    }

    [Test]
    public void AllRegisteredTools_AreGoldStandard_NameTitleReadOnlyDescription()
    {
        var methods = typeof(MolToolsTools).GetMethods(BindingFlags.Public | BindingFlags.Static);
        foreach (var method in methods)
        {
            var tool = method.GetCustomAttributes(inherit: false)
                .FirstOrDefault(a => string.Equals(a.GetType().Name, "McpServerToolAttribute", StringComparison.Ordinal));
            if (tool is null) continue;

            var name = tool.GetType().GetProperty("Name")?.GetValue(tool) as string;
            Assert.That(name, Is.Not.Null.And.Not.Empty, $"{method.Name} must have an explicit tool Name");

            var title = tool.GetType().GetProperty("Title")?.GetValue(tool) as string;
            Assert.That(title, Is.Not.Null.And.Not.Empty, $"{name} must have a non-empty Title");

            Assert.That(tool.GetType().GetProperty("ReadOnly")?.GetValue(tool), Is.Not.Null,
                $"{name} must set ReadOnly");

            var desc = method.GetCustomAttributes(inherit: false)
                .FirstOrDefault(a => string.Equals(a.GetType().Name, "DescriptionAttribute", StringComparison.Ordinal));
            var descText = desc?.GetType().GetProperty("Description")?.GetValue(desc) as string;
            Assert.That(descText, Is.Not.Null.And.Not.Empty, $"{name} must have a non-empty Description");
        }
    }

    private static string[] GetRegisteredToolNames(Type toolsType) =>
        toolsType.GetMethods(BindingFlags.Public | BindingFlags.Static)
            .Select(GetRegisteredToolName)
            .Where(n => n is not null)
            .Select(n => n!)
            .ToArray();

    private static string? GetRegisteredToolName(MethodInfo method)
    {
        var attribute = method.GetCustomAttributes(inherit: false)
            .FirstOrDefault(a => string.Equals(a.GetType().Name, "McpServerToolAttribute", StringComparison.Ordinal));
        if (attribute is null) return null;
        var explicitName = attribute.GetType().GetProperty("Name")?.GetValue(attribute) as string;
        return string.IsNullOrWhiteSpace(explicitName) ? method.Name : explicitName;
    }
}
