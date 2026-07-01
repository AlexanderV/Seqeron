using System;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// Server contract / integration test: the Seqeron.Mcp.Sequence server must advertise EXACTLY its documented
/// tool surface, and every registered tool must meet the gold-standard binding shape
/// (explicit Name + Title + ReadOnly + Description). This pins the server's public MCP surface
/// so an accidental rename, removal, or un-annotated tool fails the build.
/// </summary>
[TestFixture]
public class ToolRegistrationContractTests
{
    private static readonly string[] ExpectedToolNames =
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

    [Test]
    public void SequenceTools_AdvertisesExactlyTheDocumentedToolSurface()
    {
        var registered = GetRegisteredToolNames(typeof(SequenceTools));

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
        var methods = typeof(SequenceTools).GetMethods(BindingFlags.Public | BindingFlags.Static);
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
