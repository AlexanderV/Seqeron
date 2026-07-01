using System;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

/// <summary>
/// Server contract / integration test: the Seqeron.Mcp.Annotation server must advertise EXACTLY its documented
/// tool surface, and every registered tool must meet the gold-standard binding shape
/// (explicit Name + Title + ReadOnly + Description). This pins the server's public MCP surface
/// so an accidental rename, removal, or un-annotated tool fails the build.
/// </summary>
[TestFixture]
public class ToolRegistrationContractTests
{
    private static readonly string[] ExpectedToolNames =
    {
        "align_mirna_to_target",
        "analyze_target_context",
        "annotate_histone_modifications",
        "annotate_regulatory_elements",
        "annotate_svs",
        "annotate_variant_on_transcripts",
        "annotate_variants",
        "assemble_breakpoint_sequence",
        "build_coexpression_network",
        "calculate_conservation",
        "calculate_tpm",
        "call_variants",
        "call_variants_from_alignment",
        "can_pair",
        "classify_mutation",
        "classify_variant",
        "cluster_discordant_pairs",
        "cluster_genes_by_expression",
        "cluster_split_reads",
        "coding_potential",
        "codon_usage",
        "compare_seed_regions",
        "cpg_observed_expected",
        "create_mirna",
        "detect_alternative_splicing",
        "detect_differential_splicing",
        "detect_isoform_switching",
        "differential_expression",
        "enrichment_score",
        "epigenetic_age",
        "filter_svs",
        "find_acceptor_sites",
        "find_accessible_regions",
        "find_branch_points",
        "find_conserved_elements",
        "find_cpg_islands",
        "find_cpg_sites",
        "find_deletions",
        "find_discordant_pairs",
        "find_dmrs",
        "find_dominant_isoforms",
        "find_donor_sites",
        "find_indels",
        "find_insertions",
        "find_methylation_sites",
        "find_microhomology",
        "find_mirna_target_sites",
        "find_orfs",
        "find_pre_mirna_hairpins",
        "find_promoter_motifs",
        "find_repetitive_elements",
        "find_retained_intron_candidates",
        "find_ribosome_binding_sites",
        "find_similar_mirnas",
        "find_skipped_exon_events",
        "find_snps",
        "find_snps_direct",
        "find_split_reads",
        "format_vcf_info",
        "generate_seed_variants",
        "genotype_sv",
        "group_by_seed_family",
        "identify_cnvs",
        "impact_level",
        "is_within_coding_region",
        "is_wobble_pair",
        "log2_transform",
        "longest_orfs_per_frame",
        "maxent_score",
        "merge_overlapping_svs",
        "methylation_from_bisulfite",
        "methylation_profile",
        "mirna_seed_sequence",
        "normalize_variant",
        "over_representation_analysis",
        "parse_gff3",
        "parse_vcf_variant",
        "pearson_correlation",
        "perform_pca",
        "predict_chromatin_state",
        "predict_gene_structure",
        "predict_genes",
        "predict_imprinted_genes",
        "predict_introns",
        "predict_pathogenicity",
        "predict_tf_binding_change",
        "predict_variant_effect",
        "quantile_normalize",
        "rna_reverse_complement",
        "rnaseq_quality_metrics",
        "segment_copy_number",
        "simulate_bisulfite_conversion",
        "site_accessibility",
        "titv_ratio",
        "to_gff3",
        "variant_statistics",
        "variants_to_vcf",
    };

    [Test]
    public void AnnotationTools_AdvertisesExactlyTheDocumentedToolSurface()
    {
        var registered = GetRegisteredToolNames(typeof(AnnotationTools));

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
        var methods = typeof(AnnotationTools).GetMethods(BindingFlags.Public | BindingFlags.Static);
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
