using System;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Server contract / integration test: the Seqeron.Mcp.Chromosome server must advertise EXACTLY its documented
/// tool surface, and every registered tool must meet the gold-standard binding shape
/// (explicit Name + Title + ReadOnly + Description). This pins the server's public MCP surface
/// so an accidental rename, removal, or un-annotated tool fails the build.
/// </summary>
[TestFixture]
public class ToolRegistrationContractTests
{
    private static readonly string[] ExpectedToolNames =
    {
        "analyze_centromere",
        "analyze_karyotype",
        "analyze_scaffolds",
        "analyze_telomeres",
        "arm_ratio",
        "assembly_statistics",
        "assess_completeness",
        "au_n",
        "classify_chromosome_by_arm_ratio",
        "compare_assemblies",
        "detect_aneuploidy",
        "detect_ploidy",
        "detect_rearrangements",
        "estimate_cell_divisions_from_telomere_length",
        "estimate_completeness_from_kmers",
        "estimate_telomere_length_from_ts_ratio",
        "extract_contigs",
        "find_gaps",
        "find_heterochromatin_regions",
        "find_repetitive_regions",
        "find_suspicious_regions",
        "find_syntenic_blocks_assemblies",
        "find_synteny_blocks",
        "find_tandem_repeats",
        "gap_distribution",
        "identify_whole_chromosome_aneuploidy",
        "length_distribution",
        "local_quality",
        "nx_curve",
        "nx_statistics",
        "predict_g_bands",
        "repeat_content",
    };

    [Test]
    public void ChromosomeTools_AdvertisesExactlyTheDocumentedToolSurface()
    {
        var registered = GetRegisteredToolNames(typeof(ChromosomeTools));

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
        var methods = typeof(ChromosomeTools).GetMethods(BindingFlags.Public | BindingFlags.Static);
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
