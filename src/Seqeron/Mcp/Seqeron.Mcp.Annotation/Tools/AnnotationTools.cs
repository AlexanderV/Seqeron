using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Core;

namespace Seqeron.Mcp.Annotation.Tools;

/// <summary>
/// MCP tools for genome annotation (ORFs, gene prediction, GFF3, motifs, repeats, codon usage).
/// </summary>
[McpServerToolType]
public class AnnotationTools
{
    /// <summary>Find all open reading frames in a DNA sequence.</summary>
    [McpServerTool(Name = "find_orfs", Title = "Annotation — Find ORFs", ReadOnly = true)]
    [Description("Find all open reading frames (ORFs) in a DNA sequence across forward and reverse strands.")]
    public static FindOrfsResult FindOrfs(
        [Description("DNA sequence to search")] string dnaSequence,
        [Description("Minimum ORF length in amino acids")] int minLength = 100,
        [Description("Whether to also search the reverse complement")] bool searchBothStrands = true,
        [Description("Whether to require a start codon (ATG/GTG/TTG)")] bool requireStartCodon = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var orfs = GenomeAnnotator
            .FindOrfs(dnaSequence, minLength, searchBothStrands, requireStartCodon)
            .Select(ToDto)
            .ToList();
        return new FindOrfsResult(orfs);
    }

    /// <summary>Return the longest ORF in each reading frame.</summary>
    [McpServerTool(Name = "longest_orfs_per_frame", Title = "Annotation — Longest ORFs Per Frame", ReadOnly = true)]
    [Description("Return the longest ORF in each of the (up to) six reading frames.")]
    public static LongestOrfsPerFrameResult LongestOrfsPerFrame(
        [Description("DNA sequence to search")] string dnaSequence,
        [Description("Whether to also search the reverse complement (frames -1..-3)")] bool searchBothStrands = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var dict = GenomeAnnotator.FindLongestOrfsPerFrame(dnaSequence, searchBothStrands);
        var frames = dict
            .OrderBy(kv => kv.Key < 0 ? 100 - kv.Key : kv.Key)
            .Select(kv => new LongestOrfFrame(kv.Key, kv.Value.HasValue ? ToDto(kv.Value.Value) : null))
            .ToList();
        return new LongestOrfsPerFrameResult(frames);
    }

    /// <summary>Locate Shine–Dalgarno motifs upstream of ORFs.</summary>
    [McpServerTool(Name = "find_ribosome_binding_sites", Title = "Annotation — Find Ribosome Binding Sites", ReadOnly = true)]
    [Description("Locate Shine–Dalgarno (ribosome binding site) motifs upstream of ORFs.")]
    public static FindRibosomeBindingSitesResult FindRibosomeBindingSites(
        [Description("DNA sequence to scan")] string dnaSequence,
        [Description("Upstream window size (nt) to search before each ORF start")] int upstreamWindow = 20,
        [Description("Minimum distance from motif to start codon (nt)")] int minDistance = 4,
        [Description("Maximum distance from motif to start codon (nt)")] int maxDistance = 15)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var sites = GenomeAnnotator
            .FindRibosomeBindingSites(dnaSequence, upstreamWindow, minDistance, maxDistance)
            .Select(t => new RibosomeBindingSite(t.position, t.sequence, t.score))
            .ToList();
        return new FindRibosomeBindingSitesResult(sites);
    }

    /// <summary>Predict gene annotations from DNA using ORF-based heuristics.</summary>
    [McpServerTool(Name = "predict_genes", Title = "Annotation — Predict Genes", ReadOnly = true)]
    [Description("Predict gene annotations from a DNA sequence using ORF-based heuristics.")]
    public static PredictGenesResult PredictGenes(
        [Description("DNA sequence to annotate")] string dnaSequence,
        [Description("Minimum ORF length in amino acids")] int minOrfLength = 100,
        [Description("Gene identifier prefix")] string prefix = "gene")
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var genes = GenomeAnnotator
            .PredictGenes(dnaSequence, minOrfLength, prefix)
            .Select(ToDto)
            .ToList();
        return new PredictGenesResult(genes);
    }

    /// <summary>Parse a GFF3 annotation document into structured features.</summary>
    [McpServerTool(Name = "parse_gff3", Title = "Annotation — Parse GFF3", ReadOnly = true)]
    [Description("Parse a GFF3 annotation document into structured features.")]
    public static ParseGff3Result ParseGff3(
        [Description("GFF3 document text (lines separated by newlines)")] string gff3Text)
    {
        if (string.IsNullOrEmpty(gff3Text))
            throw new ArgumentException("GFF3 text cannot be null or empty", nameof(gff3Text));

        var lines = gff3Text.Split('\n');
        var features = GenomeAnnotator.ParseGff3(lines).Select(ToDto).ToList();
        return new ParseGff3Result(features);
    }

    /// <summary>Serialize gene annotations to GFF3 lines.</summary>
    [McpServerTool(Name = "to_gff3", Title = "Annotation — Export GFF3", ReadOnly = true)]
    [Description("Serialize gene annotations to GFF3 lines (with header).")]
    public static ToGff3Result ToGff3(
        [Description("Gene annotations to serialize")] IReadOnlyList<GeneAnnotationDto> annotations,
        [Description("Sequence identifier for column 1")] string seqId = "seq1")
    {
        if (annotations is null || annotations.Count == 0)
            throw new ArgumentException("Annotations cannot be null or empty", nameof(annotations));

        var source = annotations.Select(FromDto);
        var lines = GenomeAnnotator.ToGff3(source, seqId).ToList();
        return new ToGff3Result(lines);
    }

    /// <summary>Find -10 and -35 bacterial promoter motifs.</summary>
    [McpServerTool(Name = "find_promoter_motifs", Title = "Annotation — Find Promoter Motifs", ReadOnly = true)]
    [Description("Find -10 (Pribnow/TATA) and -35 box bacterial promoter motifs.")]
    public static FindPromoterMotifsResult FindPromoterMotifs(
        [Description("DNA sequence to scan")] string dnaSequence)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var motifs = GenomeAnnotator
            .FindPromoterMotifs(dnaSequence)
            .Select(t => new PromoterMotif(t.position, t.type, t.sequence, t.score))
            .ToList();
        return new FindPromoterMotifsResult(motifs);
    }

    /// <summary>Score the coding potential of a DNA sequence.</summary>
    [McpServerTool(Name = "coding_potential", Title = "Annotation — Coding Potential", ReadOnly = true)]
    [Description("Score the coding potential of a DNA sequence using stop-codon and GC3 heuristics.")]
    public static CodingPotentialResult CodingPotential(
        [Description("DNA sequence to score")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var score = GenomeAnnotator.CalculateCodingPotential(sequence);
        return new CodingPotentialResult(score);
    }

    /// <summary>Find tandem and inverted repeats.</summary>
    [McpServerTool(Name = "find_repetitive_elements", Title = "Annotation — Find Repetitive Elements", ReadOnly = true)]
    [Description("Find tandem and inverted repeats in a DNA sequence.")]
    public static FindRepetitiveElementsResult FindRepetitiveElements(
        [Description("DNA sequence to scan")] string dnaSequence,
        [Description("Minimum repeat length (nt)")] int minRepeatLength = 10,
        [Description("Minimum number of tandem copies")] int minCopies = 2)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var repeats = GenomeAnnotator
            .FindRepetitiveElements(dnaSequence, minRepeatLength, minCopies)
            .Select(t => new RepetitiveElement(t.start, t.end, t.type, t.sequence))
            .ToList();
        return new FindRepetitiveElementsResult(repeats);
    }

    /// <summary>Count codon occurrences in-frame.</summary>
    [McpServerTool(Name = "codon_usage", Title = "Annotation — Codon Usage", ReadOnly = true)]
    [Description("Count codon occurrences (in-frame, 5'→3') for a coding sequence.")]
    public static CodonUsageResult CodonUsage(
        [Description("Coding DNA sequence (length should be a multiple of 3)")] string dnaSequence)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var usage = GenomeAnnotator.GetCodonUsage(dnaSequence);
        return new CodonUsageResult(usage.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    // ================================
    // Mapping helpers
    // ================================

    private static OrfDto ToDto(GenomeAnnotator.OpenReadingFrame o) =>
        new(o.Start, o.End, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence);

    private static GeneAnnotationDto ToDto(GenomeAnnotator.GeneAnnotation g) =>
        new(g.GeneId, g.Start, g.End, g.Strand, g.Type, g.Product,
            g.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value));

    private static GenomicFeatureDto ToDto(GenomeAnnotator.GenomicFeature f) =>
        new(f.FeatureId, f.Type, f.Start, f.End, f.Strand, f.Score, f.Phase,
            f.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value));

    private static GenomeAnnotator.GeneAnnotation FromDto(GeneAnnotationDto g) =>
        new(g.GeneId, g.Start, g.End, g.Strand, g.Type, g.Product,
            g.Attributes ?? new Dictionary<string, string>());

    // ================================================================
    // VariantCaller tools (12)
    // ================================================================

    #region VariantCaller

    /// <summary>Detect SNPs and indels via global alignment.</summary>
    [McpServerTool(Name = "call_variants", Title = "Variant — Call Variants", ReadOnly = true)]
    [Description("Detect SNPs and indels between two DNA sequences using global alignment.")]
    public static CallVariantsResult CallVariants(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .CallVariants(new DnaSequence(reference), new DnaSequence(query))
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect variants from already-aligned sequences.</summary>
    [McpServerTool(Name = "call_variants_from_alignment", Title = "Variant — Call Variants From Alignment", ReadOnly = true)]
    [Description("Detect variants from already-aligned reference and query strings (gaps as '-').")]
    public static CallVariantsResult CallVariantsFromAlignment(
        [Description("Aligned reference sequence (gaps as '-')")] string alignedReference,
        [Description("Aligned query sequence (gaps as '-')")] string alignedQuery)
    {
        if (string.IsNullOrEmpty(alignedReference))
            throw new ArgumentException("Aligned reference cannot be null or empty", nameof(alignedReference));
        if (string.IsNullOrEmpty(alignedQuery))
            throw new ArgumentException("Aligned query cannot be null or empty", nameof(alignedQuery));

        var variants = VariantCaller
            .CallVariantsFromAlignment(alignedReference, alignedQuery)
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect only SNPs (alignment-based).</summary>
    [McpServerTool(Name = "find_snps", Title = "Variant — Find SNPs", ReadOnly = true)]
    [Description("Detect only SNPs between two DNA sequences using global alignment.")]
    public static CallVariantsResult FindSnps(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .FindSnps(new DnaSequence(reference), new DnaSequence(query))
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect SNPs by direct positional comparison.</summary>
    [McpServerTool(Name = "find_snps_direct", Title = "Variant — Find SNPs (Direct)", ReadOnly = true)]
    [Description("Detect SNPs by direct positional comparison without alignment.")]
    public static CallVariantsResult FindSnpsDirect(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .FindSnpsDirect(reference, query)
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect only insertions.</summary>
    [McpServerTool(Name = "find_insertions", Title = "Variant — Find Insertions", ReadOnly = true)]
    [Description("Detect only insertions between two DNA sequences.")]
    public static CallVariantsResult FindInsertions(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .FindInsertions(new DnaSequence(reference), new DnaSequence(query))
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect only deletions.</summary>
    [McpServerTool(Name = "find_deletions", Title = "Variant — Find Deletions", ReadOnly = true)]
    [Description("Detect only deletions between two DNA sequences.")]
    public static CallVariantsResult FindDeletions(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .FindDeletions(new DnaSequence(reference), new DnaSequence(query))
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Detect insertions and deletions (indels).</summary>
    [McpServerTool(Name = "find_indels", Title = "Variant — Find Indels", ReadOnly = true)]
    [Description("Detect insertions and deletions (indels) between two DNA sequences.")]
    public static CallVariantsResult FindIndels(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var variants = VariantCaller
            .FindIndels(new DnaSequence(reference), new DnaSequence(query))
            .Select(ToDto)
            .ToList();
        return new CallVariantsResult(variants);
    }

    /// <summary>Classify a SNP as transition or transversion.</summary>
    [McpServerTool(Name = "classify_mutation", Title = "Variant — Classify Mutation", ReadOnly = true)]
    [Description("Classify a SNP as transition, transversion, or other.")]
    public static ClassifyMutationResult ClassifyMutation(
        [Description("Variant to classify (must include type and ref/alt alleles)")] VariantDto variant)
    {
        if (variant is null)
            throw new ArgumentNullException(nameof(variant));
        if (string.IsNullOrEmpty(variant.ReferenceAllele) || string.IsNullOrEmpty(variant.AlternateAllele))
            throw new ArgumentException("Variant alleles cannot be empty", nameof(variant));

        var mutation = VariantCaller.ClassifyMutation(FromDto(variant));
        return new ClassifyMutationResult(mutation.ToString());
    }

    /// <summary>Compute Ti/Tv ratio.</summary>
    [McpServerTool(Name = "titv_ratio", Title = "Variant — Ti/Tv Ratio", ReadOnly = true)]
    [Description("Compute transition/transversion (Ti/Tv) ratio from a list of variants.")]
    public static TiTvRatioResult TiTvRatio(
        [Description("Variants (only SNP-typed entries are counted)")] IReadOnlyList<VariantDto> variants)
    {
        if (variants is null)
            throw new ArgumentNullException(nameof(variants));

        var ratio = VariantCaller.CalculateTiTvRatio(variants.Select(FromDto));
        return new TiTvRatioResult(ratio);
    }

    /// <summary>Compute summary variant statistics.</summary>
    [McpServerTool(Name = "variant_statistics", Title = "Variant — Statistics", ReadOnly = true)]
    [Description("Compute summary variant statistics between reference and query (totals, Ti/Tv, density).")]
    public static VariantStatisticsDto VariantStatistics(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var s = VariantCaller.CalculateStatistics(new DnaSequence(reference), new DnaSequence(query));
        return new VariantStatisticsDto(
            s.TotalVariants, s.Snps, s.Insertions, s.Deletions,
            s.TiTvRatio, s.VariantDensity, s.ReferenceLength, s.QueryLength);
    }

    /// <summary>Predict the protein-level effect of a variant.</summary>
    [McpServerTool(Name = "predict_variant_effect", Title = "Variant — Predict Effect", ReadOnly = true)]
    [Description("Predict the protein-level effect of a variant in a coding sequence (Synonymous/Missense/Nonsense/StopLoss/Frameshift/Unknown).")]
    public static PredictVariantEffectResult PredictVariantEffect(
        [Description("Variant to evaluate")] VariantDto variant,
        [Description("Coding DNA sequence (CDS)")] string codingSequence,
        [Description("Position of the variant within the coding sequence (0-based)")] int variantPosition)
    {
        if (variant is null)
            throw new ArgumentNullException(nameof(variant));
        if (string.IsNullOrEmpty(codingSequence))
            throw new ArgumentException("Coding sequence cannot be null or empty", nameof(codingSequence));

        var effect = VariantCaller.PredictEffect(
            FromDto(variant),
            new DnaSequence(codingSequence),
            variantPosition);
        return new PredictVariantEffectResult(effect.ToString());
    }

    /// <summary>Call and annotate all variants between reference and query.</summary>
    [McpServerTool(Name = "annotate_variants", Title = "Variant — Annotate Variants", ReadOnly = true)]
    [Description("Call and annotate all variants between reference and query (effect + Ti/Tv classification).")]
    public static AnnotateVariantsResult AnnotateVariants(
        [Description("Reference DNA sequence")] string reference,
        [Description("Query DNA sequence")] string query,
        [Description("Whether to interpret reference as a coding sequence and predict effects")] bool isCodingSequence = false)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));

        var annotated = VariantCaller
            .AnnotateVariants(new DnaSequence(reference), new DnaSequence(query), isCodingSequence)
            .Select(ToDto)
            .ToList();
        return new AnnotateVariantsResult(annotated);
    }

    /// <summary>Format variants as VCF v4.2 lines.</summary>
    [McpServerTool(Name = "variants_to_vcf", Title = "Variant — To VCF", ReadOnly = true)]
    [Description("Format variants as VCF v4.2 lines (header + records).")]
    public static VariantsToVcfResult VariantsToVcf(
        [Description("Variants to serialize")] IReadOnlyList<VariantDto> variants,
        [Description("Chromosome identifier for column 1 (CHROM)")] string chromosome = "chr1",
        [Description("Sample column name in the VCF header")] string sampleName = "SAMPLE")
    {
        if (variants is null)
            throw new ArgumentNullException(nameof(variants));

        var lines = VariantCaller
            .ToVcfLines(variants.Select(FromDto), chromosome, sampleName)
            .ToList();
        return new VariantsToVcfResult(lines);
    }

    // ================================
    // VariantCaller mapping helpers
    // ================================

    private static VariantDto ToDto(Variant v) =>
        new(v.Position, v.ReferenceAllele, v.AlternateAllele, v.Type.ToString(), v.QueryPosition);

    private static AnnotatedVariantDto ToDto(AnnotatedVariant a) =>
        new(ToDto(a.Variant), a.Effect.ToString(), a.MutationType.ToString());

    private static Variant FromDto(VariantDto v)
    {
        if (!Enum.TryParse<VariantType>(v.Type, ignoreCase: true, out var type))
            throw new ArgumentException(
                $"Unknown variant type '{v.Type}'. Expected one of: SNP, Insertion, Deletion, MNP, Complex.",
                nameof(v));
        return new Variant(v.Position, v.ReferenceAllele, v.AlternateAllele, type, v.QueryPosition);
    }

    #endregion

    // ================================================================
    // VariantAnnotator tools (11)
    // ================================================================

    #region VariantAnnotator

    /// <summary>Classify a (ref, alt) pair as SNV/MNV/Insertion/Deletion/Indel/Complex.</summary>
    [McpServerTool(Name = "classify_variant", Title = "Variant Annotator — Classify Variant", ReadOnly = true)]
    [Description("Classify a (ref, alt) allele pair as SNV / Insertion / Deletion / MNV / Indel / Complex.")]
    public static ClassifyVariantResult ClassifyVariant(
        [Description("Reference allele")] string reference,
        [Description("Alternate allele")] string alternate)
    {
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(alternate))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(alternate));

        var t = VariantAnnotator.ClassifyVariant(reference, alternate);
        return new ClassifyVariantResult(t.ToString());
    }

    /// <summary>Normalize a variant (left-trim prefix, right-trim suffix).</summary>
    [McpServerTool(Name = "normalize_variant", Title = "Variant Annotator — Normalize Variant", ReadOnly = true)]
    [Description("Normalize a variant by left-trimming common prefixes and right-trimming common suffixes; classifies type.")]
    public static AnnotatorVariantResult NormalizeVariant(
        [Description("Chromosome / contig identifier")] string chromosome,
        [Description("1-based genomic position of the reference allele")] int position,
        [Description("Reference allele")] string reference,
        [Description("Alternate allele")] string alternate,
        [Description("Optional reference sequence (currently unused by algorithm; reserved for left-alignment)")] string? referenceSequence = null)
    {
        if (string.IsNullOrEmpty(chromosome))
            throw new ArgumentException("Chromosome cannot be null or empty", nameof(chromosome));
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(alternate))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(alternate));

        var v = VariantAnnotator.NormalizeVariant(chromosome, position, reference, alternate, referenceSequence);
        return new AnnotatorVariantResult(ToDto(v));
    }

    /// <summary>VEP-like annotation of a variant against transcripts.</summary>
    [McpServerTool(Name = "annotate_variant_on_transcripts", Title = "Variant Annotator — Annotate On Transcripts", ReadOnly = true)]
    [Description("VEP-like annotation: predict consequence, impact, codon/AA change, SIFT/PolyPhen for one variant against transcripts.")]
    public static AnnotateVariantOnTranscriptsResult AnnotateVariantOnTranscripts(
        [Description("Variant to annotate")] AnnotatorVariantDto variant,
        [Description("Transcript models to annotate against")] IReadOnlyList<TranscriptDto> transcripts,
        [Description("Optional reference sequence (used for codon/AA change prediction)")] string? referenceSequence = null,
        [Description("Optional population frequencies map (e.g., AF, gnomAD_AF)")] Dictionary<string, double>? populationFrequencies = null)
    {
        if (variant is null)
            throw new ArgumentNullException(nameof(variant));
        if (string.IsNullOrEmpty(variant.Chromosome) || string.IsNullOrEmpty(variant.Reference) || string.IsNullOrEmpty(variant.Alternate))
            throw new ArgumentException("Variant chromosome / reference / alternate cannot be empty", nameof(variant));
        if (transcripts is null || transcripts.Count == 0)
            throw new ArgumentException("Transcripts cannot be null or empty", nameof(transcripts));

        var annotations = VariantAnnotator
            .AnnotateVariant(FromDto(variant), transcripts.Select(FromDto), referenceSequence, populationFrequencies)
            .Select(ToDto)
            .ToList();
        return new AnnotateVariantOnTranscriptsResult(annotations);
    }

    /// <summary>Map a consequence to its impact level (High/Moderate/Low/Modifier).</summary>
    [McpServerTool(Name = "impact_level", Title = "Variant Annotator — Impact Level", ReadOnly = true)]
    [Description("Map a ConsequenceType to an ImpactLevel (High/Moderate/Low/Modifier).")]
    public static ImpactLevelResult ImpactLevel(
        [Description("Consequence type name (e.g., MissenseVariant, StopGained, IntronVariant)")] string consequence)
    {
        if (string.IsNullOrEmpty(consequence))
            throw new ArgumentException("Consequence cannot be null or empty", nameof(consequence));
        if (!Enum.TryParse<VariantAnnotator.ConsequenceType>(consequence, ignoreCase: true, out var c))
            throw new ArgumentException($"Unknown consequence type '{consequence}'.", nameof(consequence));

        var impact = VariantAnnotator.GetImpactLevel(c);
        return new ImpactLevelResult(impact.ToString());
    }

    /// <summary>ACMG-like pathogenicity prediction from annotation + evidence.</summary>
    [McpServerTool(Name = "predict_pathogenicity", Title = "Variant Annotator — Predict Pathogenicity", ReadOnly = true)]
    [Description("ACMG-like pathogenicity prediction combining annotation, frequency, conservation, ClinVar, and functional evidence.")]
    public static PredictPathogenicityResult PredictPathogenicity(
        [Description("Variant annotation (one transcript)")] VariantAnnotationDto annotation,
        [Description("Population allele frequency (0..1)")] double? populationFrequency = null,
        [Description("Conservation score (e.g., PhyloP)")] double? conservationScore = null,
        [Description("Whether variant is present in ClinVar")] bool inClinvar = false,
        [Description("ClinVar clinical significance string")] string? clinvarSignificance = null,
        [Description("Optional functional evidence labels (e.g., 'LOF confirmed')")] IReadOnlyList<string>? functionalEvidence = null)
    {
        if (annotation is null)
            throw new ArgumentNullException(nameof(annotation));

        var pred = VariantAnnotator.PredictPathogenicity(
            FromDto(annotation),
            populationFrequency,
            conservationScore,
            inClinvar,
            clinvarSignificance,
            functionalEvidence);

        return new PredictPathogenicityResult(new PathogenicityPredictionDto(
            pred.Classification.ToString(),
            pred.ConfidenceScore,
            pred.EvidenceCriteria.ToList(),
            pred.ClinicalSignificance,
            pred.IsActionable));
    }

    /// <summary>PhyloP/PhastCons/GERP-like conservation scores per position.</summary>
    [McpServerTool(Name = "calculate_conservation", Title = "Variant Annotator — Calculate Conservation", ReadOnly = true)]
    [Description("Compute PhyloP-, PhastCons-, GERP-like conservation scores from per-position multi-species alleles.")]
    public static CalculateConservationResult CalculateConservation(
        [Description("Aligned positions with concatenated species alleles (one nucleotide per species).")]
        IReadOnlyList<ConservationPositionInputDto> positions)
    {
        if (positions is null || positions.Count == 0)
            throw new ArgumentException("Positions cannot be null or empty", nameof(positions));

        var input = positions.Select(p =>
            (p.Chromosome, p.Position, (IReadOnlyList<char>)(p.SpeciesAlleles ?? string.Empty).ToCharArray()));

        var scores = VariantAnnotator.CalculateConservation(input)
            .Select(s => new ConservationScoreDto(
                s.Chromosome, s.Position, s.PhyloP, s.PhastCons, s.Gerp, s.ConservedSpeciesCount))
            .ToList();
        return new CalculateConservationResult(scores);
    }

    /// <summary>Identify conserved elements from per-position scores.</summary>
    [McpServerTool(Name = "find_conserved_elements", Title = "Variant Annotator — Find Conserved Elements", ReadOnly = true)]
    [Description("Identify conserved genomic elements from per-position conservation scores (runs above threshold).")]
    public static FindConservedElementsResult FindConservedElements(
        [Description("Per-position conservation scores (PhastCons used for thresholding)")] IReadOnlyList<ConservationScoreDto> scores,
        [Description("PhastCons threshold for a conserved position")] double threshold = 0.8,
        [Description("Minimum element length (positions)")] int minLength = 20)
    {
        if (scores is null || scores.Count == 0)
            throw new ArgumentException("Scores cannot be null or empty", nameof(scores));

        var input = scores.Select(s => new VariantAnnotator.ConservationScore(
            s.Chromosome, s.Position, s.PhyloP, s.PhastCons, s.Gerp, s.ConservedSpeciesCount));

        var elements = VariantAnnotator.FindConservedElements(input, threshold, minLength)
            .Select(e => new ConservedElementDto(e.Chromosome, e.Start, e.End, e.Score))
            .ToList();
        return new FindConservedElementsResult(elements);
    }

    /// <summary>Find regulatory regions overlapping a variant.</summary>
    [McpServerTool(Name = "annotate_regulatory_elements", Title = "Variant Annotator — Annotate Regulatory Elements", ReadOnly = true)]
    [Description("Find regulatory regions overlapping a variant.")]
    public static AnnotateRegulatoryElementsResult AnnotateRegulatoryElements(
        [Description("Variant to annotate")] AnnotatorVariantDto variant,
        [Description("Regulatory regions (promoter/enhancer/silencer/etc.)")] IReadOnlyList<RegulatoryRegionInputDto> regulatoryRegions)
    {
        if (variant is null)
            throw new ArgumentNullException(nameof(variant));
        if (regulatoryRegions is null)
            throw new ArgumentNullException(nameof(regulatoryRegions));

        var input = regulatoryRegions.Select(r =>
            (r.Chromosome, r.Start, r.End, r.Type, r.CellType, r.Score, (IReadOnlyList<string>)(r.TranscriptionFactors ?? new List<string>())));

        var annotations = VariantAnnotator.AnnotateRegulatoryElements(FromDto(variant), input)
            .Select(a => new RegulatoryAnnotationDto(
                a.Chromosome, a.Start, a.End, a.FeatureType, a.CellType, a.Score, a.TranscriptionFactors.ToList()))
            .ToList();
        return new AnnotateRegulatoryElementsResult(annotations);
    }

    /// <summary>Predict TF binding score changes induced by a SNV.</summary>
    [McpServerTool(Name = "predict_tf_binding_change", Title = "Variant Annotator — Predict TF Binding Change", ReadOnly = true)]
    [Description("Predict transcription-factor binding score changes induced by a SNV against IUPAC motifs.")]
    public static PredictTfBindingChangeResult PredictTfBindingChange(
        [Description("Variant (must be type SNV)")] AnnotatorVariantDto variant,
        [Description("Transcription factor motifs (IUPAC, with score thresholds)")] IReadOnlyList<TfMotifInputDto> motifs,
        [Description("Reference sequence context centred near the variant")] string referenceContext,
        [Description("0-based offset of the variant within the reference context")] int contextOffset = 20)
    {
        if (variant is null)
            throw new ArgumentNullException(nameof(variant));
        if (string.IsNullOrEmpty(referenceContext))
            throw new ArgumentException("Reference context cannot be null or empty", nameof(referenceContext));
        if (motifs is null || motifs.Count == 0)
            throw new ArgumentException("Motifs cannot be null or empty", nameof(motifs));
        if (!string.Equals(variant.Type, nameof(VariantAnnotator.VariantType.SNV), StringComparison.OrdinalIgnoreCase))
            throw new ArgumentException("Variant must be of type SNV", nameof(variant));

        var motifInput = motifs.Select(m => (m.TfName, m.Motif, m.Threshold));
        var changes = VariantAnnotator.PredictTfBindingChange(FromDto(variant), motifInput, referenceContext, contextOffset)
            .Select(c => new TfBindingChangeDto(c.TfName, c.RefScore, c.AltScore, c.ScoreDifference))
            .ToList();
        return new PredictTfBindingChangeResult(changes);
    }

    /// <summary>Build a Variant record from VCF fields.</summary>
    [McpServerTool(Name = "parse_vcf_variant", Title = "Variant Annotator — Parse VCF Variant", ReadOnly = true)]
    [Description("Build a Variant record from VCF fields and classify its type.")]
    public static AnnotatorVariantResult ParseVcfVariant(
        [Description("Chromosome / contig (CHROM)")] string chromosome,
        [Description("1-based position (POS)")] int position,
        [Description("Variant ID (ID column; '.' if none)")] string id,
        [Description("Reference allele (REF)")] string reference,
        [Description("Alternate allele (ALT)")] string alternate,
        [Description("Quality score (QUAL); null if unknown")] double? quality = null)
    {
        if (string.IsNullOrEmpty(chromosome))
            throw new ArgumentException("Chromosome cannot be null or empty", nameof(chromosome));
        if (string.IsNullOrEmpty(reference))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(reference));
        if (string.IsNullOrEmpty(alternate))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(alternate));

        var v = VariantAnnotator.ParseVcfVariant(chromosome, position, id, reference, alternate, quality);
        return new AnnotatorVariantResult(ToDto(v));
    }

    /// <summary>Format an annotation as a VCF INFO field string.</summary>
    [McpServerTool(Name = "format_vcf_info", Title = "Variant Annotator — Format VCF INFO", ReadOnly = true)]
    [Description("Format a variant annotation as a VCF INFO field string (GENE/TRANSCRIPT/CONSEQUENCE/IMPACT/HGVSp/HGVSc/SIFT/POLYPHEN).")]
    public static FormatVcfInfoResult FormatVcfInfo(
        [Description("Variant annotation to serialise")] VariantAnnotationDto annotation)
    {
        if (annotation is null)
            throw new ArgumentNullException(nameof(annotation));

        var info = VariantAnnotator.FormatAsVcfInfo(FromDto(annotation));
        return new FormatVcfInfoResult(info);
    }

    // ================================
    // VariantAnnotator mapping helpers
    // ================================

    private static AnnotatorVariantDto ToDto(VariantAnnotator.Variant v) =>
        new(v.Chromosome, v.Position, v.Reference, v.Alternate, v.Type.ToString(), v.Quality, v.Id);

    private static VariantAnnotator.Variant FromDto(AnnotatorVariantDto v)
    {
        if (!Enum.TryParse<VariantAnnotator.VariantType>(v.Type, ignoreCase: true, out var type))
            throw new ArgumentException(
                $"Unknown variant type '{v.Type}'. Expected one of: SNV, Insertion, Deletion, MNV, Complex, Indel.",
                nameof(v));
        return new VariantAnnotator.Variant(v.Chromosome, v.Position, v.Reference, v.Alternate, type, v.Quality, v.Id);
    }

    private static VariantAnnotator.Transcript FromDto(TranscriptDto t)
    {
        var exons = (t.Exons ?? Array.Empty<GenomicIntervalDto>())
            .Select(e => (e.Start, e.End))
            .ToList();
        var codingExons = (t.CodingExons ?? Array.Empty<GenomicIntervalDto>())
            .Select(e => (e.Start, e.End))
            .ToList();
        return new VariantAnnotator.Transcript(
            t.TranscriptId, t.GeneId, t.GeneName, t.Chromosome,
            t.Start, t.End, t.Strand,
            exons, codingExons,
            t.CdsStart, t.CdsEnd);
    }

    private static VariantAnnotationDto ToDto(VariantAnnotator.VariantAnnotation a) =>
        new(
            ToDto(a.Variant),
            a.TranscriptId, a.GeneId, a.GeneName,
            a.Consequence.ToString(), a.Impact.ToString(),
            a.CodonChange, a.AminoAcidChange,
            a.ProteinPosition, a.CdsPosition,
            a.SiftScore, a.PolyphenScore, a.CaddScore,
            a.ExistingVariation,
            a.PopulationFrequencies?.ToDictionary(kv => kv.Key, kv => kv.Value));

    private static VariantAnnotator.VariantAnnotation FromDto(VariantAnnotationDto a)
    {
        if (!Enum.TryParse<VariantAnnotator.ConsequenceType>(a.Consequence, ignoreCase: true, out var consequence))
            throw new ArgumentException($"Unknown consequence type '{a.Consequence}'.", nameof(a));
        if (!Enum.TryParse<VariantAnnotator.ImpactLevel>(a.Impact, ignoreCase: true, out var impact))
            throw new ArgumentException($"Unknown impact level '{a.Impact}'.", nameof(a));

        return new VariantAnnotator.VariantAnnotation(
            FromDto(a.Variant),
            a.TranscriptId, a.GeneId, a.GeneName,
            consequence, impact,
            a.CodonChange, a.AminoAcidChange,
            a.ProteinPosition, a.CdsPosition,
            a.SiftScore, a.PolyphenScore, a.CaddScore,
            a.ExistingVariation,
            a.PopulationFrequencies);
    }

    #endregion

    // ================================================================
    // EpigeneticsAnalyzer tools (13)
    // ================================================================

    #region EpigeneticsAnalyzer

    /// <summary>Return all CpG dinucleotide start positions.</summary>
    [McpServerTool(Name = "find_cpg_sites", Title = "Epigenetics — Find CpG Sites", ReadOnly = true)]
    [Description("Return all CpG dinucleotide start positions in a sequence.")]
    public static FindCpGSitesResult FindCpGSites(
        [Description("Nucleotide sequence to scan")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var positions = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();
        return new FindCpGSitesResult(positions);
    }

    /// <summary>Find candidate cytosine methylation sites and classify context.</summary>
    [McpServerTool(Name = "find_methylation_sites", Title = "Epigenetics — Find Methylation Sites", ReadOnly = true)]
    [Description("Find candidate cytosine methylation sites and classify context (CpG/CHG/CHH).")]
    public static FindMethylationSitesResult FindMethylationSites(
        [Description("Nucleotide sequence to scan")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var sites = EpigeneticsAnalyzer.FindMethylationSites(sequence).Select(ToDto).ToList();
        return new FindMethylationSitesResult(sites);
    }

    /// <summary>Compute CpG observed/expected ratio.</summary>
    [McpServerTool(Name = "cpg_observed_expected", Title = "Epigenetics — CpG O/E Ratio", ReadOnly = true)]
    [Description("Compute CpG observed/expected ratio for a sequence.")]
    public static CpGObservedExpectedResult CpGObservedExpected(
        [Description("Nucleotide sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);
        return new CpGObservedExpectedResult(ratio);
    }

    /// <summary>Identify CpG islands using Gardiner-Garden &amp; Frommer criteria.</summary>
    [McpServerTool(Name = "find_cpg_islands", Title = "Epigenetics — Find CpG Islands", ReadOnly = true)]
    [Description("Identify CpG islands using Gardiner-Garden & Frommer criteria (length, GC content, CpG O/E).")]
    public static FindCpGIslandsResult FindCpGIslands(
        [Description("Nucleotide sequence")] string sequence,
        [Description("Minimum island length (nt)")] int minLength = 200,
        [Description("Minimum GC fraction (0..1)")] double minGc = 0.5,
        [Description("Minimum CpG observed/expected ratio")] double minCpGRatio = 0.6)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var islands = EpigeneticsAnalyzer
            .FindCpGIslands(sequence, minLength, minGc, minCpGRatio)
            .Select(t => new CpGIslandDto(t.Start, t.End, t.GcContent, t.CpGRatio))
            .ToList();
        return new FindCpGIslandsResult(islands);
    }

    /// <summary>Simulate bisulfite conversion of a DNA sequence.</summary>
    [McpServerTool(Name = "simulate_bisulfite_conversion", Title = "Epigenetics — Simulate Bisulfite Conversion", ReadOnly = true)]
    [Description("Simulate bisulfite conversion (unmethylated C → T) given a set of methylated positions to protect.")]
    public static SimulateBisulfiteConversionResult SimulateBisulfiteConversion(
        [Description("Nucleotide sequence")] string sequence,
        [Description("0-based positions of methylated cytosines (protected from conversion); null for fully unmethylated")] IReadOnlyList<int>? methylatedPositions = null)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        IReadOnlySet<int>? set = methylatedPositions is null ? null : new HashSet<int>(methylatedPositions);
        var converted = EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, set);
        return new SimulateBisulfiteConversionResult(converted);
    }

    /// <summary>Compute per-CpG methylation levels from bisulfite reads.</summary>
    [McpServerTool(Name = "methylation_from_bisulfite", Title = "Epigenetics — Methylation From Bisulfite", ReadOnly = true)]
    [Description("Compute per-CpG methylation levels from bisulfite sequencing reads aligned to a reference.")]
    public static FindMethylationSitesResult MethylationFromBisulfite(
        [Description("Reference DNA sequence")] string referenceSequence,
        [Description("Bisulfite reads with 0-based start positions on the reference")] IReadOnlyList<BisulfiteReadDto> bisulfiteReads)
    {
        if (string.IsNullOrEmpty(referenceSequence))
            throw new ArgumentException("Reference sequence cannot be null or empty", nameof(referenceSequence));
        if (bisulfiteReads is null || bisulfiteReads.Count == 0)
            throw new ArgumentException("Bisulfite reads cannot be null or empty", nameof(bisulfiteReads));

        var reads = bisulfiteReads.Select(r => (r.ReadSequence, r.StartPosition));
        var sites = EpigeneticsAnalyzer
            .CalculateMethylationFromBisulfite(referenceSequence, reads)
            .Select(ToDto)
            .ToList();
        return new FindMethylationSitesResult(sites);
    }

    /// <summary>Aggregate methylation sites into a global / CpG / CHG / CHH profile.</summary>
    [McpServerTool(Name = "methylation_profile", Title = "Epigenetics — Methylation Profile", ReadOnly = true)]
    [Description("Aggregate methylation sites into a global/CpG/CHG/CHH profile with per-position methylation.")]
    public static MethylationProfileDto MethylationProfile(
        [Description("Methylation sites to aggregate")] IReadOnlyList<MethylationSiteDto> sites)
    {
        if (sites is null)
            throw new ArgumentNullException(nameof(sites));

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites.Select(FromDto));
        return ToDto(profile);
    }

    /// <summary>Identify differentially methylated regions between two samples.</summary>
    [McpServerTool(Name = "find_dmrs", Title = "Epigenetics — Find DMRs", ReadOnly = true)]
    [Description("Identify differentially methylated regions (DMRs) between two methylation samples using a sliding window.")]
    public static FindDmrsResult FindDmrs(
        [Description("Methylation sites for sample 1")] IReadOnlyList<MethylationSiteDto> sample1,
        [Description("Methylation sites for sample 2")] IReadOnlyList<MethylationSiteDto> sample2,
        [Description("Window size (bp)")] int windowSize = 1000,
        [Description("Minimum mean methylation difference for a DMR (0..1)")] double minDifference = 0.25,
        [Description("Minimum CpG count per window")] int minCpGCount = 3)
    {
        if (sample1 is null)
            throw new ArgumentNullException(nameof(sample1));
        if (sample2 is null)
            throw new ArgumentNullException(nameof(sample2));
        if (sample1.Count == 0 && sample2.Count == 0)
            throw new ArgumentException("Both samples cannot be empty", nameof(sample1));

        var regions = EpigeneticsAnalyzer
            .FindDMRs(
                sample1.Select(FromDto),
                sample2.Select(FromDto),
                windowSize, minDifference, minCpGCount)
            .Select(r => new DmrDto(r.Start, r.End, r.MeanDifference, r.PValue, r.CpGCount, r.Annotation))
            .ToList();
        return new FindDmrsResult(regions);
    }

    /// <summary>Predict chromatin state from histone modification signal levels.</summary>
    [McpServerTool(Name = "predict_chromatin_state", Title = "Epigenetics — Predict Chromatin State", ReadOnly = true)]
    [Description("Predict chromatin state from H3K4me3 / H3K4me1 / H3K27ac / H3K36me3 / H3K27me3 / H3K9me3 signal levels (0..1).")]
    public static PredictChromatinStateResult PredictChromatinState(
        [Description("H3K4me3 signal (active promoter), 0..1")] double h3k4me3,
        [Description("H3K4me1 signal (enhancer), 0..1")] double h3k4me1,
        [Description("H3K27ac signal (active enhancer/promoter), 0..1")] double h3k27ac,
        [Description("H3K36me3 signal (transcription), 0..1")] double h3k36me3,
        [Description("H3K27me3 signal (Polycomb repression), 0..1")] double h3k27me3,
        [Description("H3K9me3 signal (heterochromatin), 0..1")] double h3k9me3)
    {
        ValidateSignal(h3k4me3, nameof(h3k4me3));
        ValidateSignal(h3k4me1, nameof(h3k4me1));
        ValidateSignal(h3k27ac, nameof(h3k27ac));
        ValidateSignal(h3k36me3, nameof(h3k36me3));
        ValidateSignal(h3k27me3, nameof(h3k27me3));
        ValidateSignal(h3k9me3, nameof(h3k9me3));

        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3, h3k4me1, h3k27ac, h3k36me3, h3k27me3, h3k9me3);
        return new PredictChromatinStateResult(state.ToString());
    }

    /// <summary>Annotate intervals with predicted chromatin state from a single histone mark.</summary>
    [McpServerTool(Name = "annotate_histone_modifications", Title = "Epigenetics — Annotate Histone Modifications", ReadOnly = true)]
    [Description("Annotate intervals with predicted chromatin state from a histone mark and signal level.")]
    public static AnnotateHistoneModificationsResult AnnotateHistoneModifications(
        [Description("Histone modification intervals (start, end, mark, signal)")] IReadOnlyList<HistoneModificationInputDto> modifications)
    {
        if (modifications is null)
            throw new ArgumentNullException(nameof(modifications));

        var input = modifications.Select(m => (m.Start, m.End, m.Mark, m.Signal));
        var annotations = EpigeneticsAnalyzer
            .AnnotateHistoneModifications(input)
            .Select(h => new HistoneModificationAnnotationDto(
                h.Start, h.End, h.Mark, h.Signal, h.PredictedState.ToString()))
            .ToList();
        return new AnnotateHistoneModificationsResult(annotations);
    }

    /// <summary>Identify accessible chromatin regions (ATAC-seq-like) from per-position signal.</summary>
    [McpServerTool(Name = "find_accessible_regions", Title = "Epigenetics — Find Accessible Regions", ReadOnly = true)]
    [Description("Identify accessible chromatin regions (ATAC-seq-like peaks) from per-position accessibility signal.")]
    public static FindAccessibleRegionsResult FindAccessibleRegions(
        [Description("Per-position accessibility signal")] IReadOnlyList<AccessibilitySignalDto> accessibilitySignal,
        [Description("Signal threshold for an accessible position")] double threshold = 0.5,
        [Description("Minimum region width (bp)")] int minWidth = 100,
        [Description("Maximum gap between accessible positions before splitting a region (bp)")] int maxGap = 50)
    {
        if (accessibilitySignal is null)
            throw new ArgumentNullException(nameof(accessibilitySignal));

        var input = accessibilitySignal.Select(s => (s.Position, s.Signal));
        var regions = EpigeneticsAnalyzer
            .FindAccessibleRegions(input, threshold, minWidth, maxGap)
            .Select(r => new AccessibilityRegionDto(
                r.Start, r.End, r.AccessibilityScore, r.PeakType, r.NearbyGenes.ToList()))
            .ToList();
        return new FindAccessibleRegionsResult(regions);
    }

    /// <summary>Predict imprinted genes from allele-specific methylation differences.</summary>
    [McpServerTool(Name = "predict_imprinted_genes", Title = "Epigenetics — Predict Imprinted Genes", ReadOnly = true)]
    [Description("Predict imprinted genes from allele-specific (maternal vs paternal) methylation differences.")]
    public static PredictImprintedGenesResult PredictImprintedGenes(
        [Description("Per-gene allele-specific methylation values")] IReadOnlyList<ImprintedGeneInputDto> genes,
        [Description("Minimum methylation difference between alleles")] double minDifference = 0.4)
    {
        if (genes is null)
            throw new ArgumentNullException(nameof(genes));

        var input = genes.Select(g =>
            (g.GeneId, g.Start, g.End, g.MaternalMethylation, g.PaternalMethylation));
        var imprinted = EpigeneticsAnalyzer
            .PredictImprintedGenes(input, minDifference)
            .Select(g => new ImprintedGeneDto(
                g.GeneId, g.Start, g.End, g.ImprintingScore, g.ParentalOrigin, g.HasDMR))
            .ToList();
        return new PredictImprintedGenesResult(imprinted);
    }

    /// <summary>Estimate epigenetic age (Horvath-clock-style) from CpG methylation.</summary>
    [McpServerTool(Name = "epigenetic_age", Title = "Epigenetics — Epigenetic Age", ReadOnly = true)]
    [Description("Estimate epigenetic age (Horvath-clock-style) from methylation values at clock CpGs.")]
    public static EpigeneticAgeResult EpigeneticAge(
        [Description("Methylation values at clock CpGs (CpG ID → methylation 0..1)")] Dictionary<string, double> methylationAtClockCpGs,
        [Description("Clock coefficients (CpG ID → coefficient); required — caller supplies the published clock table")] Dictionary<string, double> coefficients,
        [Description("Model intercept added to the weighted sum before the Horvath inverse transform")] double intercept = 0.0)
    {
        if (methylationAtClockCpGs is null || methylationAtClockCpGs.Count == 0)
            throw new ArgumentException("Methylation map cannot be null or empty", nameof(methylationAtClockCpGs));
        if (coefficients is null || coefficients.Count == 0)
            throw new ArgumentException("Clock coefficient table cannot be null or empty", nameof(coefficients));

        var age = EpigeneticsAnalyzer.CalculateEpigeneticAge(methylationAtClockCpGs, coefficients, intercept);
        return new EpigeneticAgeResult(age);
    }

    // ================================
    // EpigeneticsAnalyzer mapping helpers
    // ================================

    private static MethylationSiteDto ToDto(EpigeneticsAnalyzer.MethylationSite s) =>
        new(s.Position, s.Type.ToString(), s.Context, s.MethylationLevel, s.Coverage);

    private static EpigeneticsAnalyzer.MethylationSite FromDto(MethylationSiteDto s)
    {
        if (!Enum.TryParse<EpigeneticsAnalyzer.MethylationType>(s.Type, ignoreCase: true, out var type))
            throw new ArgumentException(
                $"Unknown methylation type '{s.Type}'. Expected one of: CpG, CHG, CHH, N6A, N4C.",
                nameof(s));
        return new EpigeneticsAnalyzer.MethylationSite(
            s.Position, type, s.Context, s.MethylationLevel, s.Coverage);
    }

    private static MethylationProfileDto ToDto(EpigeneticsAnalyzer.MethylationProfile p) =>
        new(
            p.GlobalMethylation,
            p.CpGMethylation,
            p.CHGMethylation,
            p.CHHMethylation,
            p.TotalCpGSites,
            p.MethylatedCpGSites,
            p.MethylationByPosition
                .Select(x => new MethylationByPositionDto(x.Position, x.Level))
                .ToList());

    private static void ValidateSignal(double value, string name)
    {
        if (value < 0.0 || value > 1.0)
            throw new ArgumentOutOfRangeException(name, value, "Histone signal must be in [0, 1].");
    }

    #endregion

    // ================================================================
    // MiRnaAnalyzer tools (14)
    // ================================================================

    #region MiRnaAnalyzer

    /// <summary>Extract the canonical seed region (positions 2–8) from a miRNA sequence.</summary>
    [McpServerTool(Name = "mirna_seed_sequence", Title = "miRNA — Seed Sequence", ReadOnly = true)]
    [Description("Extract the canonical seed region (positions 2-8) from a miRNA sequence.")]
    public static MiRnaSeedResult MiRnaSeedSequence(
        [Description("miRNA nucleotide sequence (must be at least 8 nt long)")] string miRnaSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(miRnaSequence));
        if (miRnaSequence.Length < 8)
            throw new ArgumentException("miRNA sequence must be at least 8 nt long to extract seed", nameof(miRnaSequence));

        return new MiRnaSeedResult(MiRnaAnalyzer.GetSeedSequence(miRnaSequence));
    }

    /// <summary>Build a MiRna record (T→U normalized) with seed metadata.</summary>
    [McpServerTool(Name = "create_mirna", Title = "miRNA — Create Record", ReadOnly = true)]
    [Description("Build a MiRna record (T→U normalized) with seed metadata.")]
    public static CreateMiRnaResult CreateMiRna(
        [Description("miRNA name / identifier")] string name,
        [Description("miRNA nucleotide sequence")] string sequence)
    {
        if (string.IsNullOrEmpty(name))
            throw new ArgumentException("Name cannot be null or empty", nameof(name));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        return new CreateMiRnaResult(ToDto(MiRnaAnalyzer.CreateMiRna(name, sequence)));
    }

    /// <summary>Compare seed regions of two miRNAs (matches, mismatches, family identity).</summary>
    [McpServerTool(Name = "compare_seed_regions", Title = "miRNA — Compare Seed Regions", ReadOnly = true)]
    [Description("Compare seed regions of two miRNAs and report Hamming-style match counts and seed-family identity.")]
    public static CompareSeedRegionsResult CompareSeedRegions(
        [Description("First miRNA")] MiRnaDto miRna1,
        [Description("Second miRNA")] MiRnaDto miRna2)
    {
        if (miRna1 is null) throw new ArgumentNullException(nameof(miRna1));
        if (miRna2 is null) throw new ArgumentNullException(nameof(miRna2));
        if (string.IsNullOrEmpty(miRna1.Sequence) || string.IsNullOrEmpty(miRna2.Sequence))
            throw new ArgumentException("Both miRNAs must have non-empty sequences");

        var cmp = MiRnaAnalyzer.CompareSeedRegions(FromDto(miRna1), FromDto(miRna2));
        return new CompareSeedRegionsResult(cmp.Matches, cmp.Mismatches, cmp.IsSameFamily);
    }

    /// <summary>Find miRNA target sites in an mRNA (8mer/7mer-m8/7mer-A1/6mer/offset-6mer).</summary>
    [McpServerTool(Name = "find_mirna_target_sites", Title = "miRNA — Find Target Sites", ReadOnly = true)]
    [Description("Find miRNA target sites in an mRNA (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer per Bartel 2009 and TargetScan).")]
    public static FindMiRnaTargetSitesResult FindMiRnaTargetSites(
        [Description("mRNA nucleotide sequence (DNA T→U is applied internally)")] string mRnaSequence,
        [Description("Query miRNA")] MiRnaDto miRna,
        [Description("Minimum site score to report")] double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (miRna is null || string.IsNullOrEmpty(miRna.Sequence))
            throw new ArgumentException("miRNA cannot be null or empty", nameof(miRna));

        var sites = MiRnaAnalyzer
            .FindTargetSites(mRnaSequence, FromDto(miRna), minScore)
            .Select(ToDto)
            .ToList();
        return new FindMiRnaTargetSitesResult(sites);
    }

    /// <summary>Reverse-complement of an RNA sequence (A↔U, G↔C, unknowns → N).</summary>
    [McpServerTool(Name = "rna_reverse_complement", Title = "miRNA — RNA Reverse Complement", ReadOnly = true)]
    [Description("Reverse-complement of an RNA sequence (A↔U, G↔C, unknown bases → N).")]
    public static RnaReverseComplementResult RnaReverseComplement(
        [Description("RNA (or DNA) nucleotide sequence")] string rnaSequence)
    {
        if (string.IsNullOrEmpty(rnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(rnaSequence));

        return new RnaReverseComplementResult(MiRnaAnalyzer.GetReverseComplement(rnaSequence));
    }

    /// <summary>Test whether two RNA bases can pair (Watson–Crick or G-U wobble).</summary>
    [McpServerTool(Name = "can_pair", Title = "miRNA — Can Pair", ReadOnly = true)]
    [Description("Test whether two RNA bases can pair (Watson-Crick A-U, G-C, or G-U wobble).")]
    public static CanPairResult CanPair(
        [Description("First base (single character A/C/G/U/T)")] string base1,
        [Description("Second base (single character A/C/G/U/T)")] string base2)
    {
        if (base1 is null || base1.Length != 1)
            throw new ArgumentException("base1 must be a single character", nameof(base1));
        if (base2 is null || base2.Length != 1)
            throw new ArgumentException("base2 must be a single character", nameof(base2));

        return new CanPairResult(MiRnaAnalyzer.CanPair(base1[0], base2[0]));
    }

    /// <summary>Test whether two bases form a G-U wobble pair.</summary>
    [McpServerTool(Name = "is_wobble_pair", Title = "miRNA — Is Wobble Pair", ReadOnly = true)]
    [Description("Test whether two bases form a G-U wobble pair.")]
    public static IsWobblePairResult IsWobblePair(
        [Description("First base (single character)")] string base1,
        [Description("Second base (single character)")] string base2)
    {
        if (base1 is null || base1.Length != 1)
            throw new ArgumentException("base1 must be a single character", nameof(base1));
        if (base2 is null || base2.Length != 1)
            throw new ArgumentException("base2 must be a single character", nameof(base2));

        return new IsWobblePairResult(MiRnaAnalyzer.IsWobblePair(base1[0], base2[0]));
    }

    /// <summary>Align a miRNA against a target sequence and compute duplex statistics.</summary>
    [McpServerTool(Name = "align_mirna_to_target", Title = "miRNA — Align To Target", ReadOnly = true)]
    [Description("Align a miRNA against a target sequence and compute duplex statistics (matches, mismatches, G-U wobbles, gaps, free energy).")]
    public static AlignMiRnaToTargetResult AlignMiRnaToTarget(
        [Description("miRNA nucleotide sequence")] string miRnaSequence,
        [Description("Target nucleotide sequence")] string targetSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence))
            throw new ArgumentException("miRNA sequence cannot be null or empty", nameof(miRnaSequence));
        if (string.IsNullOrEmpty(targetSequence))
            throw new ArgumentException("Target sequence cannot be null or empty", nameof(targetSequence));

        var d = MiRnaAnalyzer.AlignMiRnaToTarget(miRnaSequence, targetSequence);
        return new AlignMiRnaToTargetResult(
            d.MiRnaSequence,
            d.TargetSequence,
            d.AlignmentString,
            d.Matches,
            d.Mismatches,
            d.GUWobbles,
            d.Gaps,
            d.FreeEnergy);
    }

    /// <summary>Identify pre-miRNA hairpin candidates.</summary>
    [McpServerTool(Name = "find_pre_mirna_hairpins", Title = "miRNA — Find Pre-miRNA Hairpins", ReadOnly = true)]
    [Description("Identify pre-miRNA hairpin candidates with stem/loop layout and Turner 2004 nearest-neighbor free energy.")]
    public static FindPreMiRnaHairpinsResult FindPreMiRnaHairpins(
        [Description("Nucleotide sequence to scan")] string sequence,
        [Description("Minimum hairpin length (nt)")] int minHairpinLength = 55,
        [Description("Maximum hairpin length (nt)")] int maxHairpinLength = 120,
        [Description("Mature miRNA length (nt)")] int matureLength = 22)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var hairpins = MiRnaAnalyzer
            .FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength)
            .Select(ToDto)
            .ToList();
        return new FindPreMiRnaHairpinsResult(hairpins);
    }

    /// <summary>Compute AU content and positional context score around a target site.</summary>
    [McpServerTool(Name = "analyze_target_context", Title = "miRNA — Analyze Target Context", ReadOnly = true)]
    [Description("Compute AU content and positional context score around a miRNA target site within an mRNA.")]
    public static AnalyzeTargetContextResult AnalyzeTargetContext(
        [Description("mRNA nucleotide sequence")] string mRnaSequence,
        [Description("0-based inclusive start of the target site")] int targetStart,
        [Description("0-based inclusive end of the target site")] int targetEnd,
        [Description("Flanking context window (nt) on each side of the site")] int contextWindow = 30)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (targetStart < 0 || targetEnd >= mRnaSequence.Length || targetStart > targetEnd)
            throw new ArgumentOutOfRangeException(nameof(targetStart),
                $"Target indices [{targetStart}, {targetEnd}] are out of range for sequence of length {mRnaSequence.Length}");

        var ctx = MiRnaAnalyzer.AnalyzeTargetContext(mRnaSequence, targetStart, targetEnd, contextWindow);
        return new AnalyzeTargetContextResult(
            ctx.AuContent, ctx.NearStart, ctx.NearEnd, ctx.ContextScore);
    }

    /// <summary>Estimate site accessibility from local secondary-structure density.</summary>
    [McpServerTool(Name = "site_accessibility", Title = "miRNA — Site Accessibility", ReadOnly = true)]
    [Description("Estimate miRNA target site accessibility from local secondary-structure density (lower structure ⇒ higher accessibility).")]
    public static SiteAccessibilityResult SiteAccessibility(
        [Description("mRNA nucleotide sequence")] string mRnaSequence,
        [Description("0-based inclusive site start")] int siteStart,
        [Description("0-based inclusive site end")] int siteEnd)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (siteStart < 0 || siteEnd >= mRnaSequence.Length || siteStart > siteEnd)
            throw new ArgumentOutOfRangeException(nameof(siteStart),
                $"Site indices [{siteStart}, {siteEnd}] are out of range for sequence of length {mRnaSequence.Length}");

        var acc = MiRnaAnalyzer.CalculateSiteAccessibility(mRnaSequence, siteStart, siteEnd);
        return new SiteAccessibilityResult(acc);
    }

    /// <summary>Group miRNAs by identical seed sequence.</summary>
    [McpServerTool(Name = "group_by_seed_family", Title = "miRNA — Group By Seed Family", ReadOnly = true)]
    [Description("Group miRNAs by identical seed sequence (seed family).")]
    public static GroupBySeedFamilyResult GroupBySeedFamily(
        [Description("miRNAs to group")] IReadOnlyList<MiRnaDto> miRnas)
    {
        if (miRnas is null) throw new ArgumentNullException(nameof(miRnas));

        var families = MiRnaAnalyzer
            .GroupBySeedFamily(miRnas.Select(FromDto))
            .Select(g => new SeedFamilyDto(g.SeedFamily, g.Members.Select(ToDto).ToList()))
            .ToList();
        return new GroupBySeedFamilyResult(families);
    }

    /// <summary>Find miRNAs in a database with seed regions within `maxMismatches` of a query.</summary>
    [McpServerTool(Name = "find_similar_mirnas", Title = "miRNA — Find Similar miRNAs", ReadOnly = true)]
    [Description("Find miRNAs in a database whose seed region is within `maxMismatches` Hamming distance of a query miRNA.")]
    public static FindSimilarMiRnasResult FindSimilarMiRnas(
        [Description("Query miRNA")] MiRnaDto query,
        [Description("miRNA database to search")] IReadOnlyList<MiRnaDto> database,
        [Description("Maximum allowed seed mismatches")] int maxMismatches = 1)
    {
        if (query is null) throw new ArgumentNullException(nameof(query));
        if (database is null) throw new ArgumentNullException(nameof(database));

        var matches = MiRnaAnalyzer
            .FindSimilarMiRnas(FromDto(query), database.Select(FromDto), maxMismatches)
            .Select(ToDto)
            .ToList();
        return new FindSimilarMiRnasResult(matches);
    }

    /// <summary>Enumerate single-nucleotide variants of a seed sequence.</summary>
    [McpServerTool(Name = "generate_seed_variants", Title = "miRNA — Generate Seed Variants", ReadOnly = true)]
    [Description("Enumerate single-nucleotide variants of a seed sequence (the original plus all single-position substitutions over A/C/G/U).")]
    public static GenerateSeedVariantsResult GenerateSeedVariants(
        [Description("Seed nucleotide sequence")] string seedSequence,
        [Description("Reserved for future wobble-aware expansion (currently unused by the underlying algorithm)")] bool includeWobble = true)
    {
        if (string.IsNullOrEmpty(seedSequence))
            throw new ArgumentException("Seed sequence cannot be null or empty", nameof(seedSequence));

        var variants = MiRnaAnalyzer.GenerateSeedVariants(seedSequence, includeWobble).ToList();
        return new GenerateSeedVariantsResult(variants);
    }

    // ================================
    // MiRnaAnalyzer mapping helpers
    // ================================

    private static MiRnaDto ToDto(MiRnaAnalyzer.MiRna m) =>
        new(m.Name, m.Sequence, m.SeedSequence, m.SeedStart, m.SeedEnd);

    private static MiRnaAnalyzer.MiRna FromDto(MiRnaDto m) =>
        new(m.Name, m.Sequence, m.SeedSequence, m.SeedStart, m.SeedEnd);

    private static TargetSiteDto ToDto(MiRnaAnalyzer.TargetSite s) =>
        new(s.Start, s.End, s.TargetSequence, s.MiRnaName, s.Type.ToString(),
            s.SeedMatchLength, s.Score, s.FreeEnergy, s.Alignment);

    private static PreMiRnaHairpinDto ToDto(MiRnaAnalyzer.PreMiRna p) =>
        new(p.Start, p.End, p.Sequence, p.MatureSequence, p.StarSequence, p.Structure, p.FreeEnergy);

    #endregion

    #region SpliceSitePredictor

    /// <summary>Find candidate 5' (donor) splice sites.</summary>
    [McpServerTool(Name = "find_donor_sites", Title = "Splice — Find Donor Sites", ReadOnly = true)]
    [Description("Find candidate 5' (donor) splice sites (canonical GT/GU, optional GC and U12 AT/AU).")]
    public static FindDonorSitesResult FindDonorSites(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Minimum PWM score [0,1] for a site to be reported")] double minScore = 0.5,
        [Description("Whether to include non-canonical (GC, U12 AT/AU) donor motifs")] bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var sites = SpliceSitePredictor
            .FindDonorSites(sequence, minScore, includeNonCanonical)
            .Select(ToDto)
            .ToList();
        return new FindDonorSitesResult(sites);
    }

    /// <summary>Find candidate 3' (acceptor) splice sites.</summary>
    [McpServerTool(Name = "find_acceptor_sites", Title = "Splice — Find Acceptor Sites", ReadOnly = true)]
    [Description("Find candidate 3' (acceptor) splice sites (canonical AG, optional U12 AC).")]
    public static FindAcceptorSitesResult FindAcceptorSites(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Minimum PWM score [0,1] for a site to be reported")] double minScore = 0.5,
        [Description("Whether to include non-canonical (U12 AC) acceptor motifs")] bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var sites = SpliceSitePredictor
            .FindAcceptorSites(sequence, minScore, includeNonCanonical)
            .Select(ToDto)
            .ToList();
        return new FindAcceptorSitesResult(sites);
    }

    /// <summary>Find branch-point candidates using a YNYURAC-like PWM.</summary>
    [McpServerTool(Name = "find_branch_points", Title = "Splice — Find Branch Points", ReadOnly = true)]
    [Description("Find branch-point candidates using a YNYURAC-like position weight matrix.")]
    public static FindBranchPointsResult FindBranchPoints(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Inclusive start index of the search window")] int searchStart = 0,
        [Description("Exclusive-style end index of the search window; -1 scans to end of sequence")] int searchEnd = -1,
        [Description("Minimum PWM score [0,1] for a branch point to be reported")] double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var sites = SpliceSitePredictor
            .FindBranchPoints(sequence, searchStart, searchEnd, minScore)
            .Select(ToDto)
            .ToList();
        return new FindBranchPointsResult(sites);
    }

    /// <summary>Predict introns by pairing donor/acceptor sites.</summary>
    [McpServerTool(Name = "predict_introns", Title = "Splice — Predict Introns", ReadOnly = true)]
    [Description("Predict introns by pairing donor and acceptor splice sites and locating branch points.")]
    public static PredictIntronsResult PredictIntrons(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum intron length (nt)")] int minIntronLength = 60,
        [Description("Maximum intron length (nt)")] int maxIntronLength = 100000,
        [Description("Minimum combined score for intron to be reported")] double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var introns = SpliceSitePredictor
            .PredictIntrons(sequence, minIntronLength, maxIntronLength, minScore)
            .Select(ToDto)
            .ToList();
        return new PredictIntronsResult(introns);
    }

    /// <summary>Predict exon/intron gene structure.</summary>
    [McpServerTool(Name = "predict_gene_structure", Title = "Splice — Predict Gene Structure", ReadOnly = true)]
    [Description("Predict exon/intron gene structure with greedy non-overlapping intron selection.")]
    public static PredictGeneStructureResult PredictGeneStructure(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum exon length (nt)")] int minExonLength = 30,
        [Description("Minimum intron length (nt)")] int minIntronLength = 60,
        [Description("Minimum combined score for an intron to be considered")] double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var gs = SpliceSitePredictor.PredictGeneStructure(sequence, minExonLength, minIntronLength, minScore);
        var exons = gs.Exons.Select(ToDto).ToList();
        var introns = gs.Introns.Select(ToDto).ToList();
        return new PredictGeneStructureResult(exons, introns, gs.SplicedSequence, gs.OverallScore);
    }

    /// <summary>Detect candidate alternative splicing patterns.</summary>
    [McpServerTool(Name = "detect_alternative_splicing", Title = "Splice — Detect Alternative Splicing", ReadOnly = true)]
    [Description("Detect candidate alternative splicing patterns (exon skipping, alternative 5'/3' splice sites).")]
    public static DetectAlternativeSplicingResult DetectAlternativeSplicing(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum splice-site score for events to be considered")] double minScore = 0.4)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var events = SpliceSitePredictor
            .DetectAlternativeSplicing(sequence, minScore)
            .Select(t => new AlternativeSplicingEventDto(t.Type, t.Position, t.Description))
            .ToList();
        return new DetectAlternativeSplicingResult(events);
    }

    /// <summary>Find short, moderately-scored retained-intron candidates.</summary>
    [McpServerTool(Name = "find_retained_intron_candidates", Title = "Splice — Find Retained Intron Candidates", ReadOnly = true)]
    [Description("Find short, moderately-scored intron candidates likely to be retained in some transcripts.")]
    public static FindRetainedIntronCandidatesResult FindRetainedIntronCandidates(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum combined score for an intron candidate")] double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var introns = SpliceSitePredictor
            .FindRetainedIntronCandidates(sequence, minScore)
            .Select(ToDto)
            .ToList();
        return new FindRetainedIntronCandidatesResult(introns);
    }

    /// <summary>Compute a MaxEntScan-like log-likelihood score for a splice motif.</summary>
    [McpServerTool(Name = "maxent_score", Title = "Splice — MaxEnt Score", ReadOnly = true)]
    [Description("Compute a MaxEntScan-like log-likelihood score for a donor or acceptor motif.")]
    public static MaxEntScoreResult MaxentScore(
        [Description("Splice-site motif sequence (RNA or DNA)")] string motif,
        [Description("Splice-site type: \"Donor\" or \"Acceptor\"")] string type)
    {
        if (string.IsNullOrEmpty(motif))
            throw new ArgumentException("Motif cannot be null or empty", nameof(motif));
        if (string.IsNullOrEmpty(type))
            throw new ArgumentException("Type cannot be null or empty", nameof(type));

        if (!Enum.TryParse<SpliceSitePredictor.SpliceSiteType>(type, ignoreCase: true, out var siteType))
            throw new ArgumentException(
                $"Unknown splice-site type '{type}'. Expected 'Donor' or 'Acceptor'.", nameof(type));

        var score = SpliceSitePredictor.CalculateMaxEntScore(motif, siteType);
        return new MaxEntScoreResult(score);
    }

    /// <summary>Heuristic check whether a position lies in a coding region.</summary>
    [McpServerTool(Name = "is_within_coding_region", Title = "Splice — Is Within Coding Region", ReadOnly = true)]
    [Description("Heuristic check whether a sequence position lies in a coding region downstream of an in-frame start codon.")]
    public static IsWithinCodingRegionResult IsWithinCodingRegion(
        [Description("DNA/RNA sequence")] string sequence,
        [Description("Zero-based position to test")] int position,
        [Description("Reading frame (0, 1, or 2) relative to the upstream start codon")] int frame = 0)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (position < 0 || position >= sequence.Length)
            throw new ArgumentOutOfRangeException(nameof(position),
                $"Position {position} is out of range for sequence of length {sequence.Length}.");

        var isCoding = SpliceSitePredictor.IsWithinCodingRegion(sequence, position, frame);
        return new IsWithinCodingRegionResult(isCoding);
    }

    // ================================
    // SpliceSitePredictor mapping helpers
    // ================================

    private static SpliceSiteDto ToDto(SpliceSitePredictor.SpliceSite s) =>
        new(s.Position, s.Type.ToString(), s.Motif, s.Score, s.Confidence);

    private static IntronDto ToDto(SpliceSitePredictor.Intron i) =>
        new(
            i.Start,
            i.End,
            i.Length,
            ToDto(i.DonorSite),
            ToDto(i.AcceptorSite),
            i.BranchPoint.HasValue ? ToDto(i.BranchPoint.Value) : null,
            i.Sequence,
            i.Type.ToString(),
            i.Score);

    private static ExonDto ToDto(SpliceSitePredictor.Exon e) =>
        new(e.Start, e.End, e.Length, e.Type.ToString(), e.Phase, e.Sequence);

    #endregion

    // ================================================================
    // StructuralVariantAnalyzer tools (12)
    // ================================================================

    #region StructuralVariantAnalyzer

    /// <summary>Identify discordant read pairs (interchromosomal, abnormal insert, abnormal orientation).</summary>
    [McpServerTool(Name = "find_discordant_pairs", Title = "SV — Find Discordant Pairs", ReadOnly = true)]
    [Description("Identify discordant read pairs (interchromosomal, abnormal insert size, or abnormal orientation).")]
    public static FindDiscordantPairsResult FindDiscordantPairs(
        [Description("Read pairs to evaluate")] IReadOnlyList<ReadPairInputDto> readPairs,
        [Description("Expected insert size (nt) for the library")] int expectedInsertSize = 400,
        [Description("Insert-size standard deviation (nt)")] int insertSizeStdDev = 50,
        [Description("Hard maximum insert size (nt) above which the pair is always discordant")] int maxInsertSize = 10000)
    {
        if (readPairs is null) throw new ArgumentNullException(nameof(readPairs));

        var input = readPairs.Select(p =>
            (p.ReadId, p.Chr1, p.Pos1, p.Strand1, p.Chr2, p.Pos2, p.Strand2, p.InsertSize));
        var pairs = StructuralVariantAnalyzer
            .FindDiscordantPairs(input, expectedInsertSize, insertSizeStdDev, maxInsertSize)
            .Select(ToDto)
            .ToList();
        return new FindDiscordantPairsResult(pairs);
    }

    /// <summary>Cluster discordant read pairs into structural-variant candidates.</summary>
    [McpServerTool(Name = "cluster_discordant_pairs", Title = "SV — Cluster Discordant Pairs", ReadOnly = true)]
    [Description("Cluster discordant read pairs into structural-variant (SV) candidates.")]
    public static ClusterDiscordantPairsResult ClusterDiscordantPairs(
        [Description("Discordant read pairs to cluster")] IReadOnlyList<ReadPairSignatureDto> discordantPairs,
        [Description("Maximum distance (nt) between pair anchors in a cluster")] int clusterDistance = 500,
        [Description("Minimum supporting pairs per cluster")] int minSupport = 3)
    {
        if (discordantPairs is null) throw new ArgumentNullException(nameof(discordantPairs));

        var variants = StructuralVariantAnalyzer
            .ClusterDiscordantPairs(discordantPairs.Select(FromDto), clusterDistance, minSupport)
            .Select(ToDto)
            .ToList();
        return new ClusterDiscordantPairsResult(variants);
    }

    /// <summary>Find split reads from soft-clipped CIGAR alignments.</summary>
    [McpServerTool(Name = "find_split_reads", Title = "SV — Find Split Reads", ReadOnly = true)]
    [Description("Find split reads from soft-clipped CIGAR alignments.")]
    public static FindSplitReadsResult FindSplitReads(
        [Description("Aligned reads with CIGAR strings")] IReadOnlyList<AlignmentInputDto> alignments,
        [Description("Minimum soft-clip length (nt) to call a split read")] int minClipLength = 20)
    {
        if (alignments is null) throw new ArgumentNullException(nameof(alignments));

        var input = alignments.Select(a => (a.ReadId, a.Chromosome, a.Position, a.Cigar, a.Sequence));
        var reads = StructuralVariantAnalyzer
            .FindSplitReads(input, minClipLength)
            .Select(ToDto)
            .ToList();
        return new FindSplitReadsResult(reads);
    }

    /// <summary>Cluster split reads into breakpoint candidates.</summary>
    [McpServerTool(Name = "cluster_split_reads", Title = "SV — Cluster Split Reads", ReadOnly = true)]
    [Description("Cluster split reads into breakpoint candidates.")]
    public static ClusterSplitReadsResult ClusterSplitReads(
        [Description("Split reads to cluster")] IReadOnlyList<SplitReadDto> splitReads,
        [Description("Maximum distance (nt) between primary positions in a cluster")] int clusterDistance = 10,
        [Description("Minimum supporting reads per cluster")] int minSupport = 2)
    {
        if (splitReads is null) throw new ArgumentNullException(nameof(splitReads));

        var breakpoints = StructuralVariantAnalyzer
            .ClusterSplitReads(splitReads.Select(FromDto), clusterDistance, minSupport)
            .Select(ToDto)
            .ToList();
        return new ClusterSplitReadsResult(breakpoints);
    }

    /// <summary>Segment copy-number probe data into log-ratio plateaus.</summary>
    [McpServerTool(Name = "segment_copy_number", Title = "SV — Segment Copy Number", ReadOnly = true)]
    [Description("Segment copy-number probe data into log-ratio plateaus (CBS-like).")]
    public static SegmentCopyNumberResult SegmentCopyNumber(
        [Description("Copy-number probe observations")] IReadOnlyList<CopyNumberProbeDto> probes,
        [Description("Log-ratio change threshold for starting a new segment")] double changeThreshold = 0.3,
        [Description("Minimum probes per segment")] int minProbes = 5)
    {
        if (probes is null) throw new ArgumentNullException(nameof(probes));

        var input = probes.Select(p => (p.Chromosome, p.Position, p.LogRatio, p.Baf));
        var segments = StructuralVariantAnalyzer
            .SegmentCopyNumber(input, changeThreshold, minProbes)
            .Select(ToDto)
            .ToList();
        return new SegmentCopyNumberResult(segments);
    }

    /// <summary>Convert non-baseline copy-number segments into deletion / duplication SVs.</summary>
    [McpServerTool(Name = "identify_cnvs", Title = "SV — Identify CNVs", ReadOnly = true)]
    [Description("Convert non-baseline copy-number segments into deletion / duplication structural variants.")]
    public static IdentifyCnvsResult IdentifyCnvs(
        [Description("Copy-number segments")] IReadOnlyList<CopyNumberSegmentDto> segments,
        [Description("Diploid baseline copy number")] int normalCopyNumber = 2,
        [Description("Minimum segment length (nt) to emit as a CNV")] int minLength = 10000)
    {
        if (segments is null) throw new ArgumentNullException(nameof(segments));

        var variants = StructuralVariantAnalyzer
            .IdentifyCNVs(segments.Select(FromDto), normalCopyNumber, minLength)
            .Select(ToDto)
            .ToList();
        return new IdentifyCnvsResult(variants);
    }

    /// <summary>Merge structural variants whose reciprocal overlap exceeds a threshold.</summary>
    [McpServerTool(Name = "merge_overlapping_svs", Title = "SV — Merge Overlapping SVs", ReadOnly = true)]
    [Description("Merge structural variants of the same type whose reciprocal overlap fraction exceeds a threshold.")]
    public static MergeOverlappingSvsResult MergeOverlappingSvs(
        [Description("Structural variants to merge")] IReadOnlyList<StructuralVariantDto> variants,
        [Description("Minimum overlap fraction (0..1) to merge two SVs")] double overlapFraction = 0.5)
    {
        if (variants is null) throw new ArgumentNullException(nameof(variants));

        var merged = StructuralVariantAnalyzer
            .MergeOverlappingSVs(variants.Select(FromDto), overlapFraction)
            .Select(ToDto)
            .ToList();
        return new MergeOverlappingSvsResult(merged);
    }

    /// <summary>Filter SVs by quality, support, and length bounds.</summary>
    [McpServerTool(Name = "filter_svs", Title = "SV — Filter SVs", ReadOnly = true)]
    [Description("Filter structural variants by quality, supporting-read count, and length bounds.")]
    public static FilterSvsResult FilterSvs(
        [Description("Structural variants to filter")] IReadOnlyList<StructuralVariantDto> variants,
        [Description("Minimum quality score")] double minQuality = 20,
        [Description("Minimum supporting reads")] int minSupport = 2,
        [Description("Minimum SV length (nt)")] int minLength = 50,
        [Description("Maximum SV length (nt)")] int maxLength = 100_000_000)
    {
        if (variants is null) throw new ArgumentNullException(nameof(variants));

        var filtered = StructuralVariantAnalyzer
            .FilterSVs(variants.Select(FromDto), minQuality, minSupport, minLength, maxLength)
            .Select(ToDto)
            .ToList();
        return new FilterSvsResult(filtered);
    }

    /// <summary>Annotate SVs with overlapping genes and exons and a coarse impact level.</summary>
    [McpServerTool(Name = "annotate_svs", Title = "SV — Annotate SVs", ReadOnly = true)]
    [Description("Annotate structural variants with overlapping genes and exons and a coarse functional impact (HIGH / MODERATE / MODIFIER / LOW).")]
    public static AnnotateSvsResult AnnotateSvs(
        [Description("Structural variants to annotate")] IReadOnlyList<StructuralVariantDto> variants,
        [Description("Gene models with exon intervals")] IReadOnlyList<SvGeneInputDto> genes)
    {
        if (variants is null) throw new ArgumentNullException(nameof(variants));
        if (genes is null) throw new ArgumentNullException(nameof(genes));

        var geneInput = genes.Select(g =>
            (g.GeneId, g.Chromosome, g.Start, g.End,
             (IReadOnlyList<(int Start, int End)>)((g.Exons ?? Array.Empty<SvExonIntervalDto>())
                 .Select(e => (e.Start, e.End)).ToList())));

        var annotations = StructuralVariantAnalyzer
            .AnnotateSVs(variants.Select(FromDto), geneInput)
            .Select(ToDto)
            .ToList();
        return new AnnotateSvsResult(annotations);
    }

    /// <summary>Genotype an SV from reference / alternate read counts.</summary>
    [McpServerTool(Name = "genotype_sv", Title = "SV — Genotype SV", ReadOnly = true)]
    [Description("Genotype a structural variant (0/0, 0/1, 1/1, ./.) from reference and alternate supporting-read counts.")]
    public static GenotypeSvResult GenotypeSv(
        [Description("Structural variant to genotype")] StructuralVariantDto sv,
        [Description("Reference-supporting reads")] int refReads,
        [Description("Alternate-supporting reads")] int altReads,
        [Description("Total reads spanning the variant locus")] int totalReads)
    {
        if (sv is null) throw new ArgumentNullException(nameof(sv));
        if (refReads < 0) throw new ArgumentOutOfRangeException(nameof(refReads), refReads, "refReads must be non-negative");
        if (altReads < 0) throw new ArgumentOutOfRangeException(nameof(altReads), altReads, "altReads must be non-negative");
        if (totalReads < 0) throw new ArgumentOutOfRangeException(nameof(totalReads), totalReads, "totalReads must be non-negative");

        var (genotype, quality) = StructuralVariantAnalyzer.GenotypeSV(FromDto(sv), refReads, altReads, totalReads);
        return new GenotypeSvResult(genotype, quality);
    }

    /// <summary>Heuristically assemble a breakpoint junction from split-read clipped sequences.</summary>
    [McpServerTool(Name = "assemble_breakpoint_sequence", Title = "SV — Assemble Breakpoint Sequence", ReadOnly = true)]
    [Description("Heuristically assemble a breakpoint-junction sequence from split-read clipped fragments (returns null if no reads).")]
    public static AssembleBreakpointSequenceResult AssembleBreakpointSequence(
        [Description("Split reads supporting a breakpoint")] IReadOnlyList<SplitReadDto> splitReads,
        [Description("Minimum overlap (nt) between assembled fragments")] int minOverlap = 10)
    {
        if (splitReads is null) throw new ArgumentNullException(nameof(splitReads));

        var seq = StructuralVariantAnalyzer.AssembleBreakpointSequence(splitReads.Select(FromDto), minOverlap);
        return new AssembleBreakpointSequenceResult(seq);
    }

    /// <summary>Find shared microhomology at the junction of two flanking sequences.</summary>
    [McpServerTool(Name = "find_microhomology", Title = "SV — Find Microhomology", ReadOnly = true)]
    [Description("Find the longest shared microhomology between the 3' end of a left flank and the 5' end of a right flank (up to maxLength nt).")]
    public static FindMicrohomologyResult FindMicrohomology(
        [Description("Left flanking sequence (5' side of the junction)")] string leftFlank,
        [Description("Right flanking sequence (3' side of the junction)")] string rightFlank,
        [Description("Maximum microhomology length to search (nt)")] int maxLength = 20)
    {
        if (string.IsNullOrEmpty(leftFlank))
            throw new ArgumentException("Left flank cannot be null or empty", nameof(leftFlank));
        if (string.IsNullOrEmpty(rightFlank))
            throw new ArgumentException("Right flank cannot be null or empty", nameof(rightFlank));

        var (len, seq) = StructuralVariantAnalyzer.FindMicrohomology(leftFlank, rightFlank, maxLength);
        return new FindMicrohomologyResult(len, seq);
    }

    // ================================
    // StructuralVariantAnalyzer mapping helpers
    // ================================

    private static StructuralVariantDto ToDto(StructuralVariantAnalyzer.StructuralVariant v) =>
        new(v.Id, v.Chromosome, v.Start, v.End, v.Type.ToString(),
            v.Length, v.Quality, v.SupportingReads, v.InsertedSequence);

    private static StructuralVariantAnalyzer.StructuralVariant FromDto(StructuralVariantDto v)
    {
        if (!Enum.TryParse<StructuralVariantAnalyzer.SVType>(v.Type, ignoreCase: true, out var type))
            throw new ArgumentException(
                $"Unknown SV type '{v.Type}'. Expected one of: Deletion, Duplication, Inversion, Insertion, Translocation, ComplexRearrangement, CopyNumberVariation.",
                nameof(v));
        return new StructuralVariantAnalyzer.StructuralVariant(
            v.Id, v.Chromosome, v.Start, v.End, type,
            v.Length, v.Quality, v.SupportingReads, v.InsertedSequence);
    }

    private static ReadPairSignatureDto ToDto(StructuralVariantAnalyzer.ReadPairSignature p) =>
        new(p.ReadId, p.Chromosome1, p.Position1, p.Strand1,
            p.Chromosome2, p.Position2, p.Strand2, p.InsertSize, p.IsDiscordant);

    private static StructuralVariantAnalyzer.ReadPairSignature FromDto(ReadPairSignatureDto p) =>
        new(p.ReadId, p.Chromosome1, p.Position1, p.Strand1,
            p.Chromosome2, p.Position2, p.Strand2, p.InsertSize, p.IsDiscordant);

    private static SplitReadDto ToDto(StructuralVariantAnalyzer.SplitRead r) =>
        new(r.ReadId, r.Chromosome, r.PrimaryPosition, r.SupplementaryPosition, r.ClipLength, r.ClippedSequence);

    private static StructuralVariantAnalyzer.SplitRead FromDto(SplitReadDto r) =>
        new(r.ReadId, r.Chromosome, r.PrimaryPosition, r.SupplementaryPosition, r.ClipLength, r.ClippedSequence);

    private static BreakpointDto ToDto(StructuralVariantAnalyzer.Breakpoint b) =>
        new(b.Chromosome1, b.Position1, b.Strand1, b.Chromosome2, b.Position2, b.Strand2, b.SupportingReads, b.Quality);

    private static CopyNumberSegmentDto ToDto(StructuralVariantAnalyzer.CopyNumberSegment s) =>
        new(s.Chromosome, s.Start, s.End, s.LogRatio, s.CopyNumber, s.BAlleleFrequency, s.ProbeCount);

    private static StructuralVariantAnalyzer.CopyNumberSegment FromDto(CopyNumberSegmentDto s) =>
        new(s.Chromosome, s.Start, s.End, s.LogRatio, s.CopyNumber, s.BAlleleFrequency, s.ProbeCount);

    private static SvAnnotationDto ToDto(StructuralVariantAnalyzer.SVAnnotation a) =>
        new(a.SVId,
            a.AffectedGenes.ToList(),
            a.AffectedExons.ToList(),
            a.FunctionalImpact,
            a.PopulationFrequency,
            a.IsPathogenic);

    #endregion

    #region TranscriptomeAnalyzer

    /// <summary>Compute Transcripts-Per-Million (and FPKM) from raw counts and gene lengths.</summary>
    [McpServerTool(Name = "calculate_tpm", Title = "Transcriptome — Calculate TPM", ReadOnly = true)]
    [Description("Compute Transcripts-Per-Million (TPM) and FPKM expression values from per-gene raw counts and gene lengths.")]
    public static CalculateTpmResult CalculateTpm(
        [Description("Per-gene raw counts and lengths")] IReadOnlyList<GeneCountInputDto> geneCounts)
    {
        if (geneCounts is null) throw new ArgumentNullException(nameof(geneCounts));
        if (geneCounts.Count == 0)
            throw new ArgumentException("geneCounts cannot be empty.", nameof(geneCounts));

        var input = geneCounts.Select(g => (g.GeneId, g.RawCount, g.Length));
        var expressions = TranscriptomeAnalyzer
            .CalculateTPM(input)
            .Select(e => new GeneExpressionDto(e.GeneId, e.RawCount, e.TPM, e.FPKM, e.Length))
            .ToList();
        return new CalculateTpmResult(expressions);
    }

    /// <summary>Quantile-normalize multiple expression vectors of equal length.</summary>
    [McpServerTool(Name = "quantile_normalize", Title = "Transcriptome — Quantile Normalize", ReadOnly = true)]
    [Description("Quantile-normalize multiple expression vectors of equal length so that each sample shares the same value distribution.")]
    public static QuantileNormalizeResult QuantileNormalize(
        [Description("Expression vectors per sample (all of equal length)")] IReadOnlyList<IReadOnlyList<double>> samples)
    {
        if (samples is null) throw new ArgumentNullException(nameof(samples));
        if (samples.Count == 0)
            throw new ArgumentException("samples cannot be empty.", nameof(samples));

        int expectedLength = samples[0]?.Count ?? 0;
        for (int i = 0; i < samples.Count; i++)
        {
            if (samples[i] is null)
                throw new ArgumentException($"samples[{i}] is null.", nameof(samples));
            if (samples[i].Count != expectedLength)
                throw new ArgumentException(
                    $"All samples must have the same length; samples[0]={expectedLength}, samples[{i}]={samples[i].Count}.",
                    nameof(samples));
        }

        var normalized = TranscriptomeAnalyzer
            .QuantileNormalize(samples.Select(s => (IEnumerable<double>)s))
            .Select(v => (IReadOnlyList<double>)v.ToList())
            .ToList();
        return new QuantileNormalizeResult(normalized);
    }

    /// <summary>Log2-transform values with a pseudocount.</summary>
    [McpServerTool(Name = "log2_transform", Title = "Transcriptome — Log2 Transform", ReadOnly = true)]
    [Description("Apply log2(x + pseudocount) to each value (used to stabilize variance of count-like expression data).")]
    public static Log2TransformResult Log2Transform(
        [Description("Values to transform")] IReadOnlyList<double> values,
        [Description("Pseudocount added before taking log2 (avoids log(0))")] double pseudocount = 1.0)
    {
        if (values is null) throw new ArgumentNullException(nameof(values));

        var transformed = TranscriptomeAnalyzer
            .Log2Transform(values, pseudocount)
            .ToList();
        return new Log2TransformResult(transformed);
    }

    /// <summary>Differential expression analysis with log2 fold change, t-test, and BH correction.</summary>
    [McpServerTool(Name = "differential_expression", Title = "Transcriptome — Differential Expression", ReadOnly = true)]
    [Description("Run a simple two-group differential-expression analysis using log2 fold change, Welch-style t-test p-values, and Benjamini–Hochberg FDR correction.")]
    public static DifferentialExpressionResult DifferentialExpression(
        [Description("Per-gene replicates for the two groups")] IReadOnlyList<DiffExprInputDto> expressionData,
        [Description("Minimum |log2 fold change| required for significance")] double foldChangeThreshold = 1.0,
        [Description("Adjusted-p-value threshold for significance")] double pValueThreshold = 0.05)
    {
        if (expressionData is null) throw new ArgumentNullException(nameof(expressionData));
        if (expressionData.Count == 0)
            throw new ArgumentException("expressionData cannot be empty.", nameof(expressionData));

        var input = expressionData.Select(d => (d.GeneId, d.Group1, d.Group2));
        var results = TranscriptomeAnalyzer
            .AnalyzeDifferentialExpression(input, foldChangeThreshold, pValueThreshold)
            .Select(r => new DifferentialExpressionDto(
                r.GeneId, r.Log2FoldChange, r.PValue, r.AdjustedPValue, r.IsSignificant, r.Regulation))
            .ToList();
        return new DifferentialExpressionResult(results);
    }

    /// <summary>Pathway over-representation analysis (hypergeometric / Fisher approximation).</summary>
    [McpServerTool(Name = "over_representation_analysis", Title = "Transcriptome — Over-Representation Analysis", ReadOnly = true)]
    [Description("Pathway / gene-set over-representation analysis (ORA) using a hypergeometric / Fisher's exact-test approximation.")]
    public static OverRepresentationAnalysisResult OverRepresentationAnalysis(
        [Description("Differentially expressed gene IDs")] IReadOnlyList<string> differentiallyExpressedGenes,
        [Description("Pathway / gene-set definitions")] IReadOnlyList<PathwayInputDto> pathways,
        [Description("Total background gene count (universe size)")] int backgroundGeneCount)
    {
        if (differentiallyExpressedGenes is null) throw new ArgumentNullException(nameof(differentiallyExpressedGenes));
        if (pathways is null) throw new ArgumentNullException(nameof(pathways));
        if (differentiallyExpressedGenes.Count == 0)
            throw new ArgumentException("differentiallyExpressedGenes cannot be empty.", nameof(differentiallyExpressedGenes));
        if (backgroundGeneCount <= 0)
            throw new ArgumentException("backgroundGeneCount must be positive.", nameof(backgroundGeneCount));

        var deSet = new HashSet<string>(differentiallyExpressedGenes);
        var pathwayInput = pathways.Select(p =>
            (p.PathwayId, p.PathwayName, (IReadOnlySet<string>)new HashSet<string>(p.Genes)));
        var results = TranscriptomeAnalyzer
            .PerformOverRepresentationAnalysis(deSet, pathwayInput, backgroundGeneCount)
            .Select(r => new EnrichmentResultDto(
                r.PathwayId, r.PathwayName, r.GenesInPathway, r.OverlappingGenes,
                r.EnrichmentScore, r.PValue, r.Genes.ToList()))
            .ToList();
        return new OverRepresentationAnalysisResult(results);
    }

    /// <summary>Compute a GSEA-like running-sum enrichment score over a ranked gene list.</summary>
    [McpServerTool(Name = "enrichment_score", Title = "Transcriptome — Enrichment Score", ReadOnly = true)]
    [Description("Compute a GSEA-like running-sum enrichment score for a gene set against a ranked gene list.")]
    public static EnrichmentScoreResult EnrichmentScore(
        [Description("Gene IDs ordered by ranking metric (e.g., descending log2 fold change)")] IReadOnlyList<string> rankedGenes,
        [Description("Gene IDs in the gene set being tested")] IReadOnlyList<string> geneSet)
    {
        if (rankedGenes is null) throw new ArgumentNullException(nameof(rankedGenes));
        if (geneSet is null) throw new ArgumentNullException(nameof(geneSet));
        if (rankedGenes.Count == 0)
            throw new ArgumentException("rankedGenes cannot be empty.", nameof(rankedGenes));
        if (geneSet.Count == 0)
            throw new ArgumentException("geneSet cannot be empty.", nameof(geneSet));

        var score = TranscriptomeAnalyzer.CalculateEnrichmentScore(rankedGenes, new HashSet<string>(geneSet));
        return new EnrichmentScoreResult(score);
    }

    /// <summary>Compute Percent Spliced In (PSI) for skipped-exon candidates.</summary>
    [McpServerTool(Name = "find_skipped_exon_events", Title = "Transcriptome — Find Skipped Exon Events", ReadOnly = true)]
    [Description("Compute Percent Spliced In (PSI = inclusion / (inclusion + skipping)) for candidate skipped-exon events.")]
    public static FindSkippedExonEventsResult FindSkippedExonEvents(
        [Description("Per-exon inclusion and skipping read counts")] IReadOnlyList<SkippedExonInputDto> exonData)
    {
        if (exonData is null) throw new ArgumentNullException(nameof(exonData));

        var input = exonData.Select(e => (e.GeneId, e.ExonStart, e.ExonEnd, e.InclusionReads, e.SkippingReads));
        var events = TranscriptomeAnalyzer
            .FindSkippedExonEvents(input)
            .Select(ToDto)
            .ToList();
        return new FindSkippedExonEventsResult(events);
    }

    /// <summary>Detect splicing events with |delta-PSI| above a threshold across two conditions.</summary>
    [McpServerTool(Name = "detect_differential_splicing", Title = "Transcriptome — Detect Differential Splicing", ReadOnly = true)]
    [Description("Detect splicing events whose PSI differs between two conditions by more than the given delta-PSI threshold.")]
    public static DetectDifferentialSplicingResult DetectDifferentialSplicing(
        [Description("Per-event PSI values in the two conditions")] IReadOnlyList<SplicingComparisonInputDto> splicingData,
        [Description("Minimum |PSI(cond2) - PSI(cond1)| to report an event")] double deltaPsiThreshold = 0.1)
    {
        if (splicingData is null) throw new ArgumentNullException(nameof(splicingData));

        var input = splicingData.Select(s => (s.GeneId, s.Start, s.End, s.PsiCondition1, s.PsiCondition2));
        var events = TranscriptomeAnalyzer
            .DetectDifferentialSplicing(input, deltaPsiThreshold)
            .Select(ToDto)
            .ToList();
        return new DetectDifferentialSplicingResult(events);
    }

    /// <summary>Identify the dominant transcript isoform per gene with its dominance ratio.</summary>
    [McpServerTool(Name = "find_dominant_isoforms", Title = "Transcriptome — Find Dominant Isoforms", ReadOnly = true)]
    [Description("For each gene, identify the most-expressed transcript isoform and report its share of the total per-gene expression.")]
    public static FindDominantIsoformsResult FindDominantIsoforms(
        [Description("Transcript isoforms with expression values")] IReadOnlyList<TranscriptIsoformDto> isoforms)
    {
        if (isoforms is null) throw new ArgumentNullException(nameof(isoforms));

        var input = isoforms.Select(FromDto);
        var dominants = TranscriptomeAnalyzer
            .FindDominantIsoforms(input)
            .Select(t => new DominantIsoformDto(t.GeneId, ToDto(t.DominantIsoform), t.DominanceRatio))
            .ToList();
        return new FindDominantIsoformsResult(dominants);
    }

    /// <summary>Detect isoform-usage switching between two conditions.</summary>
    [McpServerTool(Name = "detect_isoform_switching", Title = "Transcriptome — Detect Isoform Switching", ReadOnly = true)]
    [Description("Detect genes where the dominant transcript isoform switches between two conditions based on changes in usage proportions.")]
    public static DetectIsoformSwitchingResult DetectIsoformSwitching(
        [Description("Per-isoform expression in two conditions")] IReadOnlyList<IsoformExpressionInputDto> isoformData,
        [Description("Minimum |delta-usage| (0..1) for an isoform to be considered switching")] double switchThreshold = 0.3)
    {
        if (isoformData is null) throw new ArgumentNullException(nameof(isoformData));

        var input = isoformData.Select(d => (FromDto(d.Isoform), d.Expression1, d.Expression2));
        var switches = TranscriptomeAnalyzer
            .DetectIsoformSwitching(input, switchThreshold)
            .Select(t => new IsoformSwitchDto(t.GeneId, t.TranscriptId1, t.TranscriptId2, t.SwitchScore))
            .ToList();
        return new DetectIsoformSwitchingResult(switches);
    }

    /// <summary>Compute Pearson correlation between two expression vectors.</summary>
    [McpServerTool(Name = "pearson_correlation", Title = "Transcriptome — Pearson Correlation", ReadOnly = true)]
    [Description("Compute the Pearson product-moment correlation coefficient between two expression vectors of equal length.")]
    public static PearsonCorrelationResult PearsonCorrelation(
        [Description("First expression vector")] IReadOnlyList<double> expression1,
        [Description("Second expression vector (same length as expression1)")] IReadOnlyList<double> expression2)
    {
        if (expression1 is null) throw new ArgumentNullException(nameof(expression1));
        if (expression2 is null) throw new ArgumentNullException(nameof(expression2));
        if (expression1.Count == 0 || expression2.Count == 0)
            throw new ArgumentException("Expression vectors cannot be empty.");
        if (expression1.Count != expression2.Count)
            throw new ArgumentException(
                $"Expression vectors must have equal length; got {expression1.Count} vs {expression2.Count}.");

        var corr = TranscriptomeAnalyzer.CalculatePearsonCorrelation(expression1, expression2);
        return new PearsonCorrelationResult(corr);
    }

    /// <summary>Build a co-expression network keeping pairs with |correlation| ≥ threshold.</summary>
    [McpServerTool(Name = "build_coexpression_network", Title = "Transcriptome — Build Co-Expression Network", ReadOnly = true)]
    [Description("Build a gene co-expression network by emitting all gene pairs whose Pearson correlation magnitude meets or exceeds the threshold.")]
    public static BuildCoexpressionNetworkResult BuildCoexpressionNetwork(
        [Description("Per-gene expression profiles across samples")] IReadOnlyList<GeneProfileInputDto> geneProfiles,
        [Description("Minimum |Pearson correlation| required for an edge")] double correlationThreshold = 0.7)
    {
        if (geneProfiles is null) throw new ArgumentNullException(nameof(geneProfiles));

        var input = geneProfiles.Select(g => (g.GeneId, g.Expression));
        var edges = TranscriptomeAnalyzer
            .BuildCoExpressionNetwork(input, correlationThreshold)
            .Select(e => new CoExpressionEdgeDto(e.Gene1, e.Gene2, e.Correlation))
            .ToList();
        return new BuildCoexpressionNetworkResult(edges);
    }

    /// <summary>Cluster genes by expression-profile correlation (k-means-like).</summary>
    [McpServerTool(Name = "cluster_genes_by_expression", Title = "Transcriptome — Cluster Genes By Expression", ReadOnly = true)]
    [Description("Cluster genes by similarity of their expression profiles using a k-means-like procedure on Pearson correlation.")]
    public static ClusterGenesByExpressionResult ClusterGenesByExpression(
        [Description("Per-gene expression profiles across samples")] IReadOnlyList<GeneProfileInputDto> geneProfiles,
        [Description("Number of clusters to form")] int numClusters = 5,
        [Description("Minimum mean within-cluster correlation (currently informational)")] double correlationThreshold = 0.5)
    {
        if (geneProfiles is null) throw new ArgumentNullException(nameof(geneProfiles));

        var input = geneProfiles.Select(g => (g.GeneId, g.Expression));
        var clusters = TranscriptomeAnalyzer
            .ClusterGenesByExpression(input, numClusters, correlationThreshold)
            .Select(c => new CoExpressionClusterDto(
                c.ClusterId,
                c.Genes.ToList(),
                c.MeanCorrelation,
                c.RepresentativeGene,
                c.EnrichedFunctions.ToList()))
            .ToList();
        return new ClusterGenesByExpressionResult(clusters);
    }

    /// <summary>Compute basic RNA-seq quality-control metrics.</summary>
    [McpServerTool(Name = "rnaseq_quality_metrics", Title = "Transcriptome — RNA-seq QC Metrics", ReadOnly = true)]
    [Description("Compute basic RNA-seq QC metrics: mapping rate, exonic rate, rRNA rate, and number of detected genes.")]
    public static RnaseqQualityMetricsResult RnaseqQualityMetrics(
        [Description("Total sequenced reads")] double totalReads,
        [Description("Reads mapped to the reference")] double mappedReads,
        [Description("Reads mapped to exonic regions")] double exonicReads,
        [Description("Reads mapped to rRNA loci")] double rRnaReads,
        [Description("Per-gene read counts (used to count detected genes with count > 0)")] IReadOnlyList<double> geneCounts)
    {
        if (geneCounts is null) throw new ArgumentNullException(nameof(geneCounts));
        if (totalReads < 0) throw new ArgumentException("totalReads cannot be negative.", nameof(totalReads));
        if (mappedReads < 0) throw new ArgumentException("mappedReads cannot be negative.", nameof(mappedReads));
        if (exonicReads < 0) throw new ArgumentException("exonicReads cannot be negative.", nameof(exonicReads));
        if (rRnaReads < 0) throw new ArgumentException("rRnaReads cannot be negative.", nameof(rRnaReads));

        var (mappingRate, exonicRate, rrnaRate, detectedGenes) =
            TranscriptomeAnalyzer.CalculateQualityMetrics(totalReads, mappedReads, exonicReads, rRnaReads, geneCounts);
        return new RnaseqQualityMetricsResult(mappingRate, exonicRate, rrnaRate, detectedGenes);
    }

    /// <summary>Project samples onto top-2 principal components from the most variable genes.</summary>
    [McpServerTool(Name = "perform_pca", Title = "Transcriptome — Perform PCA", ReadOnly = true)]
    [Description("Project samples onto the first two principal components computed from the most variable genes (approximate PCA).")]
    public static PerformPcaResult PerformPca(
        [Description("Per-sample expression vectors (all of equal length)")] IReadOnlyList<SampleExpressionInputDto> samples,
        [Description("Number of top-variable genes to use")] int topGenes = 500)
    {
        if (samples is null) throw new ArgumentNullException(nameof(samples));
        if (samples.Count == 0)
            throw new ArgumentException("samples cannot be empty.", nameof(samples));

        int expectedLength = samples[0]?.Expression?.Count ?? 0;
        for (int i = 0; i < samples.Count; i++)
        {
            if (samples[i] is null || samples[i].Expression is null)
                throw new ArgumentException($"samples[{i}] is null.", nameof(samples));
            if (samples[i].Expression.Count != expectedLength)
                throw new ArgumentException(
                    $"All sample expression vectors must have the same length; samples[0]={expectedLength}, samples[{i}]={samples[i].Expression.Count}.",
                    nameof(samples));
        }

        var input = samples.Select(s => (s.SampleId, s.Expression));
        var points = TranscriptomeAnalyzer
            .PerformPCA(input, topGenes)
            .Select(p => new PcaPointDto(p.SampleId, p.PC1, p.PC2))
            .ToList();
        return new PerformPcaResult(points);
    }

    // ================================
    // TranscriptomeAnalyzer mapping helpers
    // ================================

    private static SplicingEventDto ToDto(TranscriptomeAnalyzer.SplicingEvent e) =>
        new(e.GeneId, e.EventType, e.Start, e.End, e.InclusionLevel, e.DeltaPSI);

    private static TranscriptIsoformDto ToDto(TranscriptomeAnalyzer.TranscriptIsoform i) =>
        new(i.TranscriptId, i.GeneId, i.Length, i.ExonCount, i.Expression, i.IsProteinCoding,
            i.Exons.Select(e => new ExonRangeDto(e.Start, e.End)).ToList());

    private static TranscriptomeAnalyzer.TranscriptIsoform FromDto(TranscriptIsoformDto i) =>
        new(i.TranscriptId, i.GeneId, i.Length, i.ExonCount, i.Expression, i.IsProteinCoding,
            i.Exons.Select(e => (e.Start, e.End)).ToList());

    #endregion
}
