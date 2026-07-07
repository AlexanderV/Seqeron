using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Server contract / integration test: the Seqeron.Mcp.Analysis server must advertise EXACTLY its documented
/// tool surface, and every registered tool must meet the gold-standard binding shape
/// (explicit Name + Title + ReadOnly + Description). This pins the server's public MCP surface
/// so an accidental rename, removal, or un-annotated tool fails the build.
/// </summary>
[TestFixture]
public class ToolRegistrationContractTests
{
    private static readonly string[] ExpectedToolNames =
    {
        "analyze_gc_content",
        "analyze_kmers",
        "at_skew",
        "base_pair_type",
        "bulge_loop_energy",
        "calculate_ani",
        "can_pair",
        "codon_frequencies",
        "compare_genomes",
        "compression_ratio",
        "count_kmers",
        "count_kmers_both_strands",
        "create_pwm",
        "cumulative_gc_skew",
        "dangling_end_energy",
        "detect_pseudoknots",
        "detect_rearrangements",
        "dinucleotide_frequencies",
        "dinucleotide_ratios",
        "discover_motifs",
        "disorder_propensity",
        "dust_score",
        "entropy_profile",
        "find_clumps",
        "find_common_regions",
        "find_conserved_clusters",
        "find_degenerate_motif",
        "find_direct_repeats",
        "find_exact_motif",
        "find_inverted_repeats",
        "find_known_motifs",
        "find_low_complexity_regions",
        "find_microsatellites",
        "find_motif",
        "find_motif_by_pattern",
        "find_motif_by_prosite",
        "find_open_reading_frames",
        "find_orthologs",
        "find_palindromes",
        "find_protein_domains",
        "find_protein_low_complexity_regions",
        "find_protein_motifs",
        "find_reciprocal_best_hits",
        "find_regulatory_elements",
        "find_repeats",
        "find_rna_inverted_repeats",
        "find_shared_motifs",
        "find_stem_loops",
        "find_syntenic_blocks",
        "find_tandem_repeats",
        "flush_coaxial_stacking",
        "gc_content_profile",
        "gc_skew",
        "generate_all_kmers",
        "generate_consensus",
        "generate_dot_plot",
        "hairpin_loop_energy",
        "hydrophobicity_profile",
        "internal_loop_energy",
        "is_disorder_promoting",
        "kmer_distance",
        "kmer_frequencies",
        "kmer_positions",
        "kmer_spectrum",
        "kmers_with_min_count",
        "mask_low_complexity",
        "minimum_free_energy",
        "mismatch_coaxial_stacking",
        "most_frequent_kmers",
        "multibranch_loop_energy",
        "parse_dot_bracket",
        "predict_chou_fasman",
        "predict_coiled_coils",
        "predict_disorder",
        "predict_low_complexity_seg",
        "predict_morfs",
        "predict_replication_origin",
        "predict_rna_structure",
        "predict_signal_peptide",
        "predict_transmembrane_helices",
        "prosite_to_regex",
        "reversal_distance",
        "rna_complement_base",
        "scan_with_pwm",
        "stem_energy",
        "tandem_repeat_summary",
        "terminal_mismatch_energy",
        "unique_kmers",
        "validate_dot_bracket",
        "windowed_complexity",
        "windowed_gc_skew",
    };

    [Test]
    public void AnalysisTools_AdvertisesExactlyTheDocumentedToolSurface()
    {
        var registered = GetRegisteredToolNames(typeof(AnalysisTools));

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
        var methods = typeof(AnalysisTools).GetMethods(BindingFlags.Public | BindingFlags.Static);
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
