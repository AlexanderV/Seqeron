using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Annotation.Tools;

/// <summary>
/// MCP tools for genome annotation (ORFs, gene prediction, GFF3, motifs, repeats, codon usage).
/// </summary>
[McpServerToolType]
public class AnnotationTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private AnnotationTools() { }

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
    [Description("Locate Shine–Dalgarno (ribosome binding site) motifs upstream of ORFs. score = motif length / 6 (scoreMode MotifLength, default) or the Prodigal log-odds weight of the hit's motif/spacer bin trained on the input genome itself (scoreMode LogOdds; needs >= 20000 nt; Hyatt et al. 2010).")]
    public static FindRibosomeBindingSitesResult FindRibosomeBindingSites(
        [Description("DNA sequence to scan")] string dnaSequence,
        [Description("Upstream window size (nt) to search before each ORF start")] int upstreamWindow = 20,
        [Description("Minimum distance from motif to start codon (nt)")] int minDistance = 4,
        [Description("Maximum distance from motif to start codon (nt)")] int maxDistance = 15,
        [Description("Score: MotifLength (motif length / 6, default) or LogOdds (Prodigal bin log-odds trained on the input genome)")] string scoreMode = "MotifLength")
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));
        if (!Enum.TryParse<RbsScoreMode>(scoreMode, ignoreCase: true, out var mode) || !Enum.IsDefined(mode))
            throw new ArgumentException("scoreMode must be MotifLength or LogOdds", nameof(scoreMode));

        var sites = GenomeAnnotator
            .FindRibosomeBindingSites(dnaSequence, upstreamWindow, minDistance, maxDistance, mode)
            .Select(t => new RibosomeBindingSite(t.position, t.sequence, t.score))
            .ToList();
        return new FindRibosomeBindingSitesResult(sites);
    }

    /// <summary>Predict gene annotations from DNA using ORF-based heuristics.</summary>
    [McpServerTool(Name = "predict_genes", Title = "Annotation — Predict Genes", ReadOnly = true)]
    [Description("Predict gene annotations from a DNA sequence using ORF-based heuristics. startSelection LongestOrf (default): most upstream ATG/GTG/TTG per stop; ProdigalScore: per stop the start with the best Prodigal coding+RBS log-odds+start-type+upstream score; ProdigalDynamicProgramming: Prodigal's complete gene model. Prodigal modes train on the input genome (>= 20000 nt) and are numerically identical to Prodigal V2.6.3.")]
    public static PredictGenesResult PredictGenes(
        [Description("DNA sequence to annotate")] string dnaSequence,
        [Description("Minimum ORF length in amino acids")] int minOrfLength = 100,
        [Description("Gene identifier prefix")] string prefix = "gene",
        [Description("Start selection: LongestOrf (default), ProdigalScore or ProdigalDynamicProgramming")] string startSelection = "LongestOrf")
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));
        if (!Enum.TryParse<StartSelection>(startSelection, ignoreCase: true, out var selection) || !Enum.IsDefined(selection))
            throw new ArgumentException("startSelection must be LongestOrf, ProdigalScore or ProdigalDynamicProgramming", nameof(startSelection));

        var genes = GenomeAnnotator
            .PredictGenes(dnaSequence, minOrfLength, prefix, selection)
            .Select(ToDto)
            .ToList();
        return new PredictGenesResult(genes);
    }

    /// <summary>Parse a GFF3 annotation document into structured features.</summary>
    [McpServerTool(Name = "parse_gff3", Title = "Annotation — Parse GFF3", ReadOnly = true)]
    [Description("Parse a GFF3 (SO spec v1.26) document into features: 1-based inclusive coordinates, percent-decoded seqid/source/type/attributes, multi-valued attributes (e.g. Parent=a,b) as value lists in AttributeValues. Stops at ##FASTA (or a '>' line); spec-invalid feature lines (not 9 columns, start<1 or start>end, strand not +-.?, phase not 0/1/2, malformed or repeated column-9 tags) are skipped.")]
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
    [Description("Serialize gene annotations (0-based half-open) to GFF3 lines with the ##gff-version 3 header: start+1/end, source and score or '.', CDS phase accumulated per GeneId 5'->3', ID/product then attributes (AttributeValues written comma-separated), spec percent-encoding (seqid outside [a-zA-Z0-9.:^*$@!+_?-|], column-9 ;=&, escaped). Rejects strands other than +-.? and empty/inverted intervals.")]
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

    /// <summary>Find the -35 / -10 boxes of bacterial sigma70 promoters (delegates to GenomeAnnotator.FindPromoterMotifs → MotifFinder.FindSigma70Promoters).</summary>
    [McpServerTool(Name = "find_promoter_motifs", Title = "Annotation — Find Promoter Motifs", ReadOnly = true)]
    [Description("Find the -35 and -10 (Pribnow) boxes of bacterial sigma70 promoters: every forward-strand pair of a -35 hexamer within 2 mismatches of TTGACA and a -10 hexamer within 2 mismatches of TATAAT separated by a 15-21 bp spacer (Harley & Reynolds 1987); each box reported once, score = 1 - mismatches/6. For pairs/strand/spacer use find_sigma70_promoters; for a free-energy promoter score use predict_sigma70_promoters.")]
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
    [Description("Score the coding potential of a DNA sequence using the CPAT hexamer usage-bias log-likelihood (Wang et al. 2013). Positive = coding, negative = non-coding.")]
    public static CodingPotentialResult CodingPotential(
        [Description("DNA sequence to score")] string sequence,
        [Description("In-frame hexamer frequency table from a coding (CDS) training set: 6-mer -> count/proportion")] IReadOnlyDictionary<string, double> codingHexamerFrequencies,
        [Description("In-frame hexamer frequency table from a non-coding (background) training set: 6-mer -> count/proportion")] IReadOnlyDictionary<string, double> noncodingHexamerFrequencies)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var score = GenomeAnnotator.CalculateCodingPotential(
            sequence, codingHexamerFrequencies, noncodingHexamerFrequencies);
        return new CodingPotentialResult(score);
    }

    /// <summary>Fickett (1982) TESTCODE of a whole sequence.</summary>
    [McpServerTool(Name = "fickett_score", Title = "Annotation — Fickett TESTCODE", ReadOnly = true)]
    [Description("Fickett (1982) TESTCODE coding statistic of a whole DNA sequence (position asymmetry + base composition look-up tables). Call: < 0.74 non-coding, >= 0.95 coding, else no opinion. tableSet 'fickett1982' (published / EMBOSS tables, default) or 'cpat' (CPAT fickett_value, the feature of CPAT's logistic model).")]
    public static FickettScoreResult FickettScore(
        [Description("DNA sequence")] string sequence,
        [Description("'fickett1982' (default) or 'cpat'")] string tableSet = "fickett1982")
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var tables = tableSet?.Trim().ToLowerInvariant() switch
        {
            "fickett1982" => GenomeAnnotator.FickettTableSet.Fickett1982,
            "cpat" => GenomeAnnotator.FickettTableSet.Cpat,
            _ => throw new ArgumentException("tableSet must be 'fickett1982' or 'cpat'", nameof(tableSet)),
        };
        double score = GenomeAnnotator.CalculateFickettScore(sequence, tables);
        return new FickettScoreResult(score, GenomeAnnotator.ClassifyFickettScore(score).ToString());
    }

    /// <summary>Sliding-window Fickett TESTCODE (EMBOSS tcode).</summary>
    [McpServerTool(Name = "fickett_scan", Title = "Annotation — Fickett TESTCODE Scan", ReadOnly = true)]
    [Description("Slide a window along a DNA sequence and report the Fickett (1982) TESTCODE and coding call of each window, as EMBOSS tcode (window >= 200, default 200; step default 3). Coordinates are 0-based half-open.")]
    public static FickettScanResult FickettScan(
        [Description("DNA sequence")] string sequence,
        [Description("Window length (>= 200)")] int windowSize = 200,
        [Description("Window step (>= 1)")] int stepSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var windows = GenomeAnnotator.ScanFickettScore(sequence, windowSize, stepSize)
            .Select(w => new FickettWindowDto(w.Start, w.End, w.Score, w.Call.ToString()))
            .ToList();
        return new FickettScanResult(windows);
    }

    /// <summary>CPAT coding probability of a transcript.</summary>
    [McpServerTool(Name = "cpat_coding_probability", Title = "Annotation — CPAT Coding Probability", ReadOnly = true)]
    [Description("CPAT (Wang et al. 2013): find ORF candidates in a transcript, score each by transcript length, ORF length, Fickett TESTCODE and hexamer score, and return the logistic coding probability (coding when >= cutoff) plus the best candidate. Requires the hexamer tables and logistic coefficients of a CPAT model (e.g. CPAT's Human model: intercept -7.8284881840360567, mRNA -5.4831277888113823e-05, ORF 0.011187815405519514, Fickett 3.4471695641643225, Hexamer 6.068275924110611, cutoff 0.364).")]
    public static CpatCodingProbabilityResult CpatCodingProbability(
        [Description("Transcript DNA sequence (T, not U)")] string sequence,
        [Description("Coding hexamer table: 6-mer -> frequency")] IReadOnlyDictionary<string, double> codingHexamerFrequencies,
        [Description("Non-coding hexamer table: 6-mer -> frequency (same units)")] IReadOnlyDictionary<string, double> noncodingHexamerFrequencies,
        [Description("Logistic intercept")] double intercept,
        [Description("Coefficient of transcript length")] double mrnaLengthCoefficient,
        [Description("Coefficient of ORF length")] double orfLengthCoefficient,
        [Description("Coefficient of the Fickett score")] double fickettCoefficient,
        [Description("Coefficient of the hexamer score")] double hexamerCoefficient,
        [Description("Coding-probability cutoff (e.g. 0.364 for CPAT Human)")] double codingProbabilityCutoff,
        [Description("Also search the antisense strand (CPAT --antisense)")] bool searchAntisense = false,
        [Description("Number of longest ORF candidates scored (CPAT --top-orf, default 5)")] int topOrfCount = 5,
        [Description("ORFs must be longer than this many nt (CPAT --min-orf, default 75)")] int minOrfLength = 75,
        [Description("Start codons (default ATG); an empty list lets any non-stop codon open an ORF")] IReadOnlyList<string>? startCodons = null,
        [Description("Stop codons (default TAG, TAA, TGA)")] IReadOnlyList<string>? stopCodons = null,
        [Description("Pick the best candidate by length instead of probability (CPAT --best-orf l)")] bool bestByLength = false)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var model = new GenomeAnnotator.CpatModel(codingHexamerFrequencies, noncodingHexamerFrequencies,
            intercept, mrnaLengthCoefficient, orfLengthCoefficient, fickettCoefficient, hexamerCoefficient,
            codingProbabilityCutoff);
        var options = new GenomeAnnotator.CpatOptions { SearchAntisense = searchAntisense, TopOrfCount = topOrfCount, MinOrfLength = minOrfLength };
        if (startCodons is not null) options = options with { StartCodons = startCodons };
        if (stopCodons is not null) options = options with { StopCodons = stopCodons };

        var candidates = GenomeAnnotator.PredictCodingProbability(sequence, model, options);
        var best = GenomeAnnotator.PredictBestCodingOrf(sequence, model, options, bestByLength);
        var dtos = candidates.Select(p => new CpatOrfDto(p.Strand.ToString(), p.Frame, p.Start, p.End, p.OrfLength,
            p.OrfSequence, p.MrnaLength, p.FickettScore, p.HexamerScore, p.CodingProbability, p.IsCoding)).ToList();
        int bestIndex = best is null ? -1 : candidates.ToList().IndexOf(best);
        return new CpatCodingProbabilityResult(dtos, bestIndex);
    }

    /// <summary>Find tandem and inverted repeats.</summary>
    [McpServerTool(Name = "find_repetitive_elements", Title = "Annotation — Find Repetitive Elements", ReadOnly = true)]
    [Description("Find exact tandem repeats (maximal runs of a primitive ACGT unit of 1-60 bp, complete copies) and exact inverted repeats (maximal stems, EMBOSS palindrome model: arm >= max(minRepeatLength, 2) and <= 100, loop 0-100) in a DNA sequence.")]
    public static FindRepetitiveElementsResult FindRepetitiveElements(
        [Description("DNA sequence to scan")] string dnaSequence,
        [Description("Minimum tandem-repeat span and minimum inverted-repeat arm length (nt)")] int minRepeatLength = 10,
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
    [Description("Count codon occurrences (in-frame from position 0, non-overlapping, 5'→3') for a coding sequence. DNA or RNA, case-insensitive: U is read as T and keys are upper-case DNA codons; a triplet with any other symbol (N, IUPAC code, gap) is skipped without shifting the frame; a partial trailing codon is ignored (CodonW / EMBOSS cusp convention).")]
    public static CodonUsageResult CodonUsage(
        [Description("Coding DNA or RNA sequence (length should be a multiple of 3)")] string dnaSequence)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));

        var usage = GenomeAnnotator.GetCodonUsage(dnaSequence);
        return new CodonUsageResult(usage.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    /// <summary>Relative synonymous codon usage pooled over several coding sequences.</summary>
    [McpServerTool(Name = "pooled_rscu", Title = "Annotation — Pooled RSCU", ReadOnly = true)]
    [Description("Relative Synonymous Codon Usage (RSCU = n·x/Σx per synonymous family; Sharp, Tuohy & Mosurski 1986, CodonW) over a set of in-frame coding sequences. Codon counts are pooled across all sequences first (CodonW -totals); DNA or RNA, case-insensitive, U read as T; triplets with other symbols are skipped without shifting the frame; stop codons are excluded; an unobserved family gives RSCU 0. geneticCode is an NCBI translation table number (default 1).")]
    public static PooledRscuResult PooledRscu(
        [Description("In-frame coding DNA or RNA sequences (CDS); each is read from position 0 in steps of three")] string[] sequences,
        [Description("NCBI translation table number (default 1 = Standard)")] int geneticCode = 1)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("Sequences cannot be null or empty", nameof(sequences));

        var code = global::Seqeron.Genomics.Core.GeneticCode.GetByTableNumberOrThrow(geneticCode, nameof(geneticCode));
        var rscu = GenomeAnnotator.GetCodonUsage(sequences, code);
        return new PooledRscuResult(rscu.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    // ================================
    // Mapping helpers
    // ================================

    private static OrfDto ToDto(GenomeAnnotator.OpenReadingFrame o) =>
        new(o.Start, o.End, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence);

    private static GeneAnnotationDto ToDto(GenomeAnnotator.GeneAnnotation g) =>
        new(g.GeneId, g.Start, g.End, g.Strand, g.Type, g.Product,
            g.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value),
            g.Source, g.Score, ToListDictionary(g.AttributeValues));

    private static GenomicFeatureDto ToDto(GenomeAnnotator.GenomicFeature f) =>
        new(f.FeatureId, f.Type, f.Start, f.End, f.Strand, f.Score, f.Phase,
            f.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value),
            f.Seqid, f.Source, ToListDictionary(f.AttributeValues));

    private static GenomeAnnotator.GeneAnnotation FromDto(GeneAnnotationDto g) =>
        new(g.GeneId, g.Start, g.End, g.Strand, g.Type, g.Product,
            g.Attributes ?? new Dictionary<string, string>(),
            string.IsNullOrEmpty(g.Source) ? "." : g.Source,
            g.Score,
            g.AttributeValues?.ToDictionary(kv => kv.Key, kv => (IReadOnlyList<string>)kv.Value));

    private static Dictionary<string, List<string>>? ToListDictionary(
        IReadOnlyDictionary<string, IReadOnlyList<string>>? values) =>
        values?.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());

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
    [Description("Detect SNPs by direct positional comparison without alignment. Inputs must be pre-aligned and of equal length (unequal lengths are rejected); bases compare case-insensitively and gap columns are not SNPs.")]
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
        ArgumentNullException.ThrowIfNull(variant);
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
        ArgumentNullException.ThrowIfNull(variants);

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
        ArgumentNullException.ThrowIfNull(variant);
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
        ArgumentNullException.ThrowIfNull(variants);

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
        ArgumentNullException.ThrowIfNull(variant);
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
        ArgumentNullException.ThrowIfNull(annotation);

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
        ArgumentNullException.ThrowIfNull(variant);
        ArgumentNullException.ThrowIfNull(regulatoryRegions);

        var input = regulatoryRegions.Select(r =>
            (r.Chromosome, r.Start, r.End, r.Type, r.CellType, r.Score, r.TranscriptionFactors ?? new List<string>()));

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
        ArgumentNullException.ThrowIfNull(variant);
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
        ArgumentNullException.ThrowIfNull(annotation);

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
        ArgumentNullException.ThrowIfNull(sites);

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
        ArgumentNullException.ThrowIfNull(sample1);
        ArgumentNullException.ThrowIfNull(sample2);
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
        ArgumentNullException.ThrowIfNull(modifications);

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
        ArgumentNullException.ThrowIfNull(accessibilitySignal);

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
        ArgumentNullException.ThrowIfNull(genes);

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
        ArgumentNullException.ThrowIfNull(miRna1);
        ArgumentNullException.ThrowIfNull(miRna2);
        if (string.IsNullOrEmpty(miRna1.Sequence) || string.IsNullOrEmpty(miRna2.Sequence))
            throw new ArgumentException("Both miRNAs must have non-empty sequences");

        var cmp = MiRnaAnalyzer.CompareSeedRegions(FromDto(miRna1), FromDto(miRna2));
        return new CompareSeedRegionsResult(cmp.Matches, cmp.Mismatches, cmp.IsSameFamily);
    }

    /// <summary>Find miRNA target sites in an mRNA (8mer/7mer-m8/7mer-A1/6mer/offset-6mer).</summary>
    [McpServerTool(Name = "find_mirna_target_sites", Title = "miRNA — Find Target Sites", ReadOnly = true)]
    [Description("Find miRNA target sites in an mRNA treated as the 3'UTR (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer per Bartel 2009 and TargetScan). Score = max(0, -context score) of the fitted TargetScan 5.0 model (Grimson 2007: site type + 3' pairing + local AU + position); 6mer, offset-6mer and sites starting within the first 15 nt score 0.")]
    public static FindMiRnaTargetSitesResult FindMiRnaTargetSites(
        [Description("mRNA nucleotide sequence (DNA T→U is applied internally)")] string mRnaSequence,
        [Description("Query miRNA")] MiRnaDto miRna,
        [Description("Minimum site score to report (default 0 = every site)")] double minScore = 0.0)
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

    /// <summary>RNA reverse complement (Biopython reverse_complement_rna; IUPAC codes complemented).</summary>
    [McpServerTool(Name = "rna_reverse_complement", Title = "miRNA — RNA Reverse Complement", ReadOnly = true)]
    [Description("RNA reverse complement, as Biopython reverse_complement_rna: A↔U, G↔C, DNA T→A; IUPAC ambiguity codes are complemented (R↔Y, K↔M, B↔V, D↔H; S, W, N unchanged). Case-insensitive, upper-case output for recognized codes; other characters (gaps etc.) pass through.")]
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
    [Description("Align a miRNA antiparallel (ungapped) against a target sequence and compute duplex statistics: matches, mismatches, G-U wobbles, gaps (always 0) and the Turner 2004 nearest-neighbor duplex free energy (kcal/mol, 37 °C) of the reported duplex — intermolecular initiation, stacks, internal loops, terminal AU/GU, dangling ends / terminal mismatches.")]
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
    [Description("Identify pre-miRNA hairpin candidates. DEFAULT (useConsecutiveScreen=false): every window of length minHairpinLength..maxHairpinLength is folded with the canonical Zuker-Stiegler MFE engine (Turner 2004) and accepted when the MFE structure is a single hairpin with >= 16 bp (Ambros 2003), a 3-25 nt apical loop (Bartel 2004) and MFEI >= minMfei (Zhang 2006, default 0.85); bulges/internal loops are tolerated. Reports span, MFE dot-bracket, free energy (kcal/mol), mature (5p, matureLength nt) and star (Dicer 2-nt 3' overhang) sequences. COST: one MFE fold per window, ~20 ms per 55-120 nt window (parallelised over window starts) — for long inputs (more than a few hundred nt) set useConsecutiveScreen=true: the fast consecutive-pairing screen (ends must pair uninterrupted for >= 18 bp with a 3-25 nt loop, no bulges; blunt mature/star labels; minMfei ignored).")]
    public static FindPreMiRnaHairpinsResult FindPreMiRnaHairpins(
        [Description("Nucleotide sequence to scan")] string sequence,
        [Description("Minimum hairpin length (nt)")] int minHairpinLength = 55,
        [Description("Maximum hairpin length (nt)")] int maxHairpinLength = 120,
        [Description("Mature miRNA length (nt)")] int matureLength = 22,
        [Description("false (default) = MFE-fold based detection (accurate, ~20 ms per window); true = fast consecutive-pairing screen for long inputs (no bulges)")] bool useConsecutiveScreen = false,
        [Description("MFEI threshold of the default mode (Zhang 2006: 0.85); ignored when useConsecutiveScreen=true")] double minMfei = 0.85)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (!double.IsFinite(minMfei) || minMfei < 0)
            throw new ArgumentOutOfRangeException(nameof(minMfei), "minMfei must be a finite value >= 0");

        var hairpins = MiRnaAnalyzer
            .FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength, useConsecutiveScreen, minMfei)
            .Select(ToDto)
            .ToList();
        return new FindPreMiRnaHairpinsResult(hairpins);
    }

    /// <summary>randfold randomisation test (Bonnet et al. 2004) with a seeded dinucleotide shuffle.</summary>
    [McpServerTool(Name = "calculate_randfold_pvalue", Title = "miRNA — Randfold P-value", ReadOnly = true)]
    [Description("randfold randomisation test (Bonnet, Wuyts, Rouze & Van de Peer 2004, Bioinformatics 20:2911; dinucleotide mode randfold -d): the input is MFE-folded (canonical Turner 2004 engine), `iterations` Altschul-Erickson dinucleotide shuffles are folded, and P = (1 + #{shuffle dG <= input dG}) / (iterations + 1) (|diff| <= 1e-7 counts as equal, as randfold.c). Pre-miRNAs give small P. DETERMINISTIC: the shuffles are drawn from System.Random(seed), so the same sequence + iterations + seed always return the same P-value (default seed 42); the P-value is a Monte-Carlo estimate, vary the seed to see its spread. Cost: iterations + 1 MFE folds (iterations capped at 10000). found=false for an empty sequence.")]
    public static CalculateRandfoldPValueResult CalculateRandfoldPValue(
        [Description("Candidate nucleotide sequence (RNA or DNA; T read as U)")] string sequence,
        [Description("Number of dinucleotide shuffles N (randfold default 999; 0 gives P = 1; max 10000)")] int iterations = 999,
        [Description("Seed of the System.Random that drives the shuffles (default 42); fixed seed = reproducible P")] int seed = 42,
        [Description("Minimum hairpin loop size of the folder (NNDB minimum 3)")] int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (iterations < 0 || iterations > MaxRandfoldIterations)
            throw new ArgumentOutOfRangeException(nameof(iterations), $"iterations must be in [0, {MaxRandfoldIterations}]");
        if (minLoopSize < 0)
            throw new ArgumentOutOfRangeException(nameof(minLoopSize), "minLoopSize must be >= 0");

        var r = MiRnaAnalyzer.CalculateRandfoldPValue(sequence, iterations, new Random(seed), minLoopSize);
        return r is { } v
            ? new CalculateRandfoldPValueResult(true, v.FreeEnergy, v.PValue, v.Iterations, v.LowerOrEqualCount, seed)
            : new CalculateRandfoldPValueResult(false, 0, 1, iterations, 0, seed);
    }

    private const int MaxRandfoldIterations = 10000;

    /// <summary>In-repo logistic pre-miRNA classifier (natural vs shuffled background).</summary>
    [McpServerTool(Name = "classify_pre_mirna", Title = "miRNA — Classify Pre-miRNA (in-repo logistic model)", ReadOnly = true)]
    [Description("Classify a candidate hairpin as natural pre-miRNA vs dinucleotide-shuffled background with the bundled IN-REPO ridge-logistic model over MFE / AMFE / MFEI / GC / paired-fraction features (Bonnet 2004 shuffle convention, Zhang 2006 MFEI). This is an in-repo model (13 miRBase positives, 52 shuffles, 5 features), NOT a published classifier — use classify_pre_mirna_triplet_svm for the published Triplet-SVM (Xue 2005). Returns the extracted features, P(natural) and the call at `threshold`. found=false for an empty sequence.")]
    public static ClassifyPreMiRnaResult ClassifyPreMiRna(
        [Description("Candidate nucleotide sequence (RNA or DNA; T read as U)")] string sequence,
        [Description("Decision threshold on P(natural), default 0.5")] double threshold = 0.5,
        [Description("Minimum hairpin loop size of the folder (NNDB minimum 3)")] int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (!double.IsFinite(threshold) || threshold < 0 || threshold > 1)
            throw new ArgumentOutOfRangeException(nameof(threshold), "threshold must be in [0, 1]");
        if (minLoopSize < 0)
            throw new ArgumentOutOfRangeException(nameof(minLoopSize), "minLoopSize must be >= 0");

        var c = MiRnaAnalyzer.ClassifyPreMiRna(sequence, threshold, minLoopSize);
        if (c is not { } v)
            return new ClassifyPreMiRnaResult(false, null, 0, false);
        var f = v.Features;
        return new ClassifyPreMiRnaResult(true,
            new PreMiRnaFeaturesDto(f.FreeEnergy, f.Amfe, f.Mfei, f.GcContent, f.PairedFraction, f.StemBasePairs, f.LoopSize, f.Length, f.BasePairingPropensity),
            v.NaturalProbability, v.IsNatural);
    }

    /// <summary>Published Triplet-SVM pre-miRNA classifier (Xue et al. 2005).</summary>
    [McpServerTool(Name = "classify_pre_mirna_triplet_svm", Title = "miRNA — Classify Pre-miRNA (Triplet-SVM)", ReadOnly = true)]
    [Description("Published Triplet-SVM pre-miRNA classifier (Xue, Li, He, Liu, Li & Zhang 2005, BMC Bioinformatics 6:310) with the authors' own trained libsvm RBF model: retrieves the stem-loops with >= minBasePairs paired bases per arm (2_get_stemloop.pl), codes each single stem-loop as 32 triplet-element frequencies and returns the libsvm decision value (> 0 = real pre-miRNA). The structure is the canonical MFE fold of this library unless dotBracket ('.()' notation, same length) is supplied (the package used RNAfold). Returns one candidate per retrieved single-arm-pair stem-loop (1-based start/end); found=false (empty candidates) when none is retrieved — e.g. no stem-loop with >= minBasePairs pairs per arm (default 22), a non-ACGU query, or a multi-branch structure.")]
    public static ClassifyPreMiRnaTripletSvmResult ClassifyPreMiRnaTripletSvm(
        [Description("Query nucleotide sequence (RNA or DNA; T read as U)")] string sequence,
        [Description("Optional structure in '.()' notation (same length as sequence); folded with the canonical MFE engine when omitted")] string? dotBracket = null,
        [Description("Minimum paired bases per arm of a retrieved stem-loop (>= 1; default 22, the package example)")] int minBasePairs = MiRnaAnalyzer.TripletSvmDefaultMinBasePairs,
        [Description("Minimum hairpin loop size of the fold (>= 0; only used when folding; values below the NNDB minimum 3 are raised to 3)")] int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (minBasePairs < 1)
            throw new ArgumentOutOfRangeException(nameof(minBasePairs), "minBasePairs must be >= 1");
        if (minLoopSize < 0)
            throw new ArgumentOutOfRangeException(nameof(minLoopSize), "minLoopSize must be >= 0");

        var results = MiRnaAnalyzer.ClassifyPreMiRnaTripletSvm(sequence, dotBracket, minBasePairs, minLoopSize)
            .Select(r => new TripletSvmCandidateDto(r.Start, r.End, r.Sequence, r.DotBracket, r.Features, r.DecisionValue, r.IsPreMiRna))
            .ToList();
        return new ClassifyPreMiRnaTripletSvmResult(results.Count > 0, results);
    }

    private const int MaxMiPredShuffles = 10000;

    /// <summary>Published MiPred random-forest pre-miRNA classifier (Jiang et al. 2007).</summary>
    [McpServerTool(Name = "classify_pre_mirna_mipred", Title = "miRNA — Classify Pre-miRNA (MiPred)", ReadOnly = true)]
    [Description("Published MiPred pre-miRNA classifier (Jiang, Wu, Wang, Ma, Sun & Lu 2007, Nucleic Acids Res 35:W339-W344) with the authors' own trained 1000-tree random forest over 34 features (32 triplet-element frequencies, MFE, shuffle p-value). Steps: canonical MFE fold (Turner 2004; the package used RNAfold 1.x) -> package pre-judgement (MFE <= -15 kcal/mol, >= 18 base pairs, no multiloop, length 51-137 nt; when it fails the sequence is not classified: features = null, isPreMiRna = false) -> p-value feature = (1 + #{shuffle dG <= dG})/(shuffleCount + 1) from dinucleotide shuffles drawn from System.Random(seed) (statistical: fixed seed = reproducible, other seeds vary within sampling error; the forest was trained with shuffleCount = 999) -> forest votes; realProbability = realVotes/1000, real pre-miRNA when > 0.5. COST: one MFE fold plus shuffleCount shuffle-folds of a sequence of <= 137 nt, i.e. ~1000 folds for the default 999 (measured ~9 s in a Debug build for the 71 nt hsa-mir-21 hairpin; sequences failing the pre-judgement cost one fold). shuffleCount is capped at 10000. found=false when nothing remains of the sequence after the package clean-up (characters other than A/C/G/T/U are dropped).")]
    public static ClassifyPreMiRnaMiPredResult ClassifyPreMiRnaMiPred(
        [Description("Candidate nucleotide sequence (RNA or DNA, case-insensitive; characters other than A/C/G/T/U are dropped, T read as U)")] string sequence,
        [Description("Seed of the System.Random that drives the shuffles for the p-value feature (default 1 = library default); fixed seed = reproducible result")] int seed = MiRnaAnalyzer.MiPredDefaultSeed,
        [Description("Number of dinucleotide shuffles for the p-value feature (default 999 = package and training value, 1..10000; other values change the distribution the forest was trained on)")] int shuffleCount = MiRnaAnalyzer.MiPredDefaultShuffleCount,
        [Description("Minimum hairpin loop size of the folder (NNDB minimum 3)")] int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (shuffleCount < 1 || shuffleCount > MaxMiPredShuffles)
            throw new ArgumentOutOfRangeException(nameof(shuffleCount), $"shuffleCount must be in [1, {MaxMiPredShuffles}]");
        if (minLoopSize < 0)
            throw new ArgumentOutOfRangeException(nameof(minLoopSize), "minLoopSize must be >= 0");

        if (MiRnaAnalyzer.ClassifyPreMiRnaMiPred(sequence, new Random(seed), shuffleCount, minLoopSize) is not { } r)
            return new ClassifyPreMiRnaMiPredResult(false, "", "", 0, false, 0, false, 0, null, null, 0, 0, null, false, seed, shuffleCount);
        return new ClassifyPreMiRnaMiPredResult(true, r.Sequence, r.DotBracket, r.FreeEnergy,
            r.Check.IsHairpinLike, r.Check.BasePairs, r.Check.HasMultiloop, r.Sequence.Length,
            r.PValue, r.Features, r.RealVotes, r.TreeCount, r.RealProbability, r.IsPreMiRna, seed, shuffleCount);
    }

    /// <summary>The 34 MiPred features of a sequence + structure + MFE + p-value.</summary>
    [McpServerTool(Name = "compute_mipred_features", Title = "miRNA — MiPred Features", ReadOnly = true)]
    [Description("The 34 MiPred features V1..V34 (Jiang et al. 2007; package microRNAcheck_parallel.pl): V1..V32 = frequencies of the 32 triplet elements (centre nucleotide A,U,G,C x local structure '(((', '((.', '(.(', '(..', '.((', '.(.', '..(', '...'; closing brackets not counted, positions 2..n-1) divided by the sequence length, V33 = freeEnergy, V34 = pValue (both supplied by the caller). Pure counting, no folding or shuffling: instant. The structure must be '.()' notation of the same length as the sequence (A/C/G/T/U only).")]
    public static ComputeMiPredFeaturesResult ComputeMiPredFeatures(
        [Description("Nucleotide sequence (RNA or DNA, A/C/G/T/U only)")] string sequence,
        [Description("Secondary structure in '.()' notation, same length as the sequence")] string dotBracket,
        [Description("MFE (kcal/mol) -> V33")] double freeEnergy,
        [Description("Shuffle p-value in [0, 1] -> V34")] double pValue)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrEmpty(dotBracket))
            throw new ArgumentException("dotBracket cannot be null or empty", nameof(dotBracket));
        if (!double.IsFinite(freeEnergy))
            throw new ArgumentOutOfRangeException(nameof(freeEnergy), "freeEnergy must be finite");
        if (!double.IsFinite(pValue) || pValue < 0 || pValue > 1)
            throw new ArgumentOutOfRangeException(nameof(pValue), "pValue must be in [0, 1]");

        return new ComputeMiPredFeaturesResult(MiRnaAnalyzer.ComputeMiPredFeatures(sequence, dotBracket, freeEnergy, pValue));
    }

    /// <summary>Minimum-free-energy miRNA·target duplex (RNAduplex model, Turner 2004).</summary>
    [McpServerTool(Name = "fold_mirna_duplex", Title = "miRNA — Fold miRNA·Target Duplex", ReadOnly = true)]
    [Description("Minimum-free-energy miRNA·target duplex (RNAduplex / duplexfold model): one contiguous duplex with internal loops and bulges up to maxLoopSize per side, scored with Turner 2004 nearest-neighbor parameters (kcal/mol, 37 °C). Returns the dot-bracket in RNAduplex notation (miRNA '(' then '&', then target ')'), the 0-based miRNA/target spans, the base pairs (miRNA index, target index; 0-based, 5'→3') and the free energy. found=false when no base pair is possible. DNA T is read as U.")]
    public static FoldMiRnaDuplexResult FoldMiRnaDuplex(
        [Description("miRNA nucleotide sequence")] string miRnaSequence,
        [Description("Target nucleotide sequence")] string targetSequence,
        [Description("Maximum unpaired nucleotides per side of an internal loop / bulge (default 30, as RNAduplex)")] int maxLoopSize = 30)
    {
        if (string.IsNullOrEmpty(miRnaSequence))
            throw new ArgumentException("miRNA sequence cannot be null or empty", nameof(miRnaSequence));
        if (string.IsNullOrEmpty(targetSequence))
            throw new ArgumentException("Target sequence cannot be null or empty", nameof(targetSequence));
        if (maxLoopSize < 0)
            throw new ArgumentOutOfRangeException(nameof(maxLoopSize), "maxLoopSize must be >= 0");

        var d = MiRnaAnalyzer.FoldMiRnaDuplex(miRnaSequence, targetSequence, maxLoopSize);
        if (d is null)
            return new FoldMiRnaDuplexResult(false, null);

        var v = d.Value;
        return new FoldMiRnaDuplexResult(true, new MiRnaDuplexStructureDto(
            v.MiRnaSequence, v.TargetSequence, v.MiRnaStart, v.MiRnaEnd, v.TargetStart, v.TargetEnd, v.Structure,
            v.BasePairs.Select(p => new DuplexBasePairDto(p.MiRnaIndex, p.TargetIndex)).ToList(),
            v.FreeEnergy));
    }

    /// <summary>Predict Drosha and Dicer cleavage sites and the 5p/3p products of a pri-miRNA hairpin.</summary>
    [McpServerTool(Name = "predict_drosha_dicer_cleavage", Title = "miRNA — Predict Drosha/Dicer Cleavage", ReadOnly = true)]
    [Description("Predict Drosha (basal) and Dicer (apical) cleavage sites of a pri-miRNA hairpin by the published measuring rules applied to its stem (base pairs counted up the unbranched stem from basalJunction; coordinates 0-based, spans inclusive; 2-nt 3' overhang). droshaRule: BasalJunction (default; ~11 bp from the basal junction, Han 2006) or ApicalJunction (~22 bp from the loop-stem junction, Zeng 2005). dicerRule: FivePrimeCounting (default; 22-nt 5p product, Park 2011) or ThreePrimeCounting (22-nt 3p product). Also reports the Auyeung 2013 CNNC and basal UG motifs. dotBracket (optional) is the stem structure to measure on; default = canonical Turner 2004 MFE fold. predicted=false when the junction is not a 5'-arm stem base or the stem is too short.")]
    public static PredictDroshaDicerCleavageResult PredictDroshaDicerCleavage(
        [Description("pri-miRNA nucleotide sequence (RNA or DNA)")] string sequence,
        [Description("0-based index of the 5'-arm base of the basal junction pair (first base pair of the stem)")] int basalJunction,
        [Description("Optional dot-bracket structure (balanced, same length as the sequence); default = MFE fold")] string? dotBracket = null,
        [Description("Dicer counting rule: FivePrimeCounting (default) or ThreePrimeCounting")] string dicerRule = "FivePrimeCounting",
        [Description("Drosha counting rule: BasalJunction (default) or ApicalJunction")] string droshaRule = "BasalJunction")
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (basalJunction < 0 || basalJunction >= sequence.Length)
            throw new ArgumentOutOfRangeException(nameof(basalJunction),
                $"basalJunction {basalJunction} is out of range for sequence of length {sequence.Length}");
        if (!TryParseEnum<MiRnaAnalyzer.DicerCountingRule>(dicerRule, out var dicer))
            throw new ArgumentException(
                $"Unknown dicerRule '{dicerRule}'. Valid values: {string.Join(", ", Enum.GetNames<MiRnaAnalyzer.DicerCountingRule>())}", nameof(dicerRule));
        if (!TryParseEnum<MiRnaAnalyzer.DroshaCountingRule>(droshaRule, out var drosha))
            throw new ArgumentException(
                $"Unknown droshaRule '{droshaRule}'. Valid values: {string.Join(", ", Enum.GetNames<MiRnaAnalyzer.DroshaCountingRule>())}", nameof(droshaRule));

        var c = MiRnaAnalyzer.PredictDroshaDicerCleavage(sequence, basalJunction, dotBracket, dicer, drosha);
        if (c is null)
            return new PredictDroshaDicerCleavageResult(false, null);

        var v = c.Value;
        return new PredictDroshaDicerCleavageResult(true, new DroshaDicerCleavageDto(
            v.BasalJunction, v.ApicalJunction, v.DroshaCut5Prime, v.DroshaCut3Prime,
            v.MatureStart, v.MatureEnd, v.StarStart, v.StarEnd, v.MatureSequence, v.StarSequence,
            v.Sequence, v.Structure, v.ThreePrimeOverhang, v.HasCnncMotif, v.HasBasalUgMotif,
            v.DicerRule.ToString(), v.DroshaRule.ToString()));
    }

    /// <summary>Find centered miRNA target sites (Shin et al. 2010).</summary>
    [McpServerTool(Name = "find_centered_sites", Title = "miRNA — Find Centered Sites", ReadOnly = true)]
    [Description("Find centered miRNA target sites (Shin et al. 2010, Mol Cell 38:789): sites without perfect seed pairing (miRNA nt 2-7 not all Watson-Crick) that pair the miRNA centre with >= 11 contiguous Watson-Crick pairs (nt 4-14, 5-15 or 4-15; G:U wobble does not count). Start/End (0-based inclusive) span the paired run (11 or 12 nt = seedMatchLength). TargetScan has no coefficients for centered sites, so score is 0. The 3'-compensatory-pairing exclusion is not tested.")]
    public static FindCenteredSitesResult FindCenteredSites(
        [Description("mRNA nucleotide sequence (DNA T→U is applied internally)")] string mRnaSequence,
        [Description("Query miRNA (>= 15 nt)")] MiRnaDto miRna)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (miRna is null || string.IsNullOrEmpty(miRna.Sequence))
            throw new ArgumentException("miRNA cannot be null or empty", nameof(miRna));

        var sites = MiRnaAnalyzer.FindCenteredSites(mRnaSequence, FromDto(miRna)).Select(ToDto).ToList();
        return new FindCenteredSitesResult(sites);
    }

    /// <summary>TargetScan 5.0 / Grimson 2007 context score of a seed site.</summary>
    [McpServerTool(Name = "score_target_site_grimson", Title = "miRNA — Score Target Site (Grimson 2007)", ReadOnly = true)]
    [Description("TargetScan 5.0 / Grimson et al. (2007, Mol Cell 27:91) context score of an 8mer, 7mer-m8 or 7mer-A1 site in a 3'UTR (index 0 = first nt after the stop codon): site type + 3' pairing + local AU + position contributions in log2 fold-change units (negative = repressed). siteStart/siteEnd are 0-based inclusive (an 8mer includes its A1). contextScore is null when the 1-based site start is < 15 (TargetScan 'too_close'); scored=false for site types without coefficients (Seed6mer, Offset6mer, Centered).")]
    public static ScoreTargetSiteGrimsonResult ScoreTargetSiteGrimson(
        [Description("3'UTR / mRNA nucleotide sequence (DNA T→U is applied internally)")] string mRnaSequence,
        [Description("miRNA nucleotide sequence")] string miRnaSequence,
        [Description("0-based inclusive site start")] int siteStart,
        [Description("0-based inclusive site end")] int siteEnd,
        [Description("Site type: Seed8mer, Seed7merM8 or Seed7merA1 (others are not scored)")] string siteType)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (string.IsNullOrEmpty(miRnaSequence))
            throw new ArgumentException("miRNA sequence cannot be null or empty", nameof(miRnaSequence));
        if (siteStart < 0 || siteEnd >= mRnaSequence.Length || siteStart > siteEnd)
            throw new ArgumentOutOfRangeException(nameof(siteStart),
                $"Site indices [{siteStart}, {siteEnd}] are out of range for sequence of length {mRnaSequence.Length}");
        if (!TryParseEnum<MiRnaAnalyzer.TargetSiteType>(siteType, out var type))
            throw new ArgumentException(
                $"Unknown siteType '{siteType}'. Valid values: {string.Join(", ", Enum.GetNames<MiRnaAnalyzer.TargetSiteType>())}", nameof(siteType));

        var g = MiRnaAnalyzer.ScoreTargetSiteGrimson(mRnaSequence, miRnaSequence, siteStart, siteEnd, type);
        if (g is null)
            return new ScoreTargetSiteGrimsonResult(false, null);

        var v = g.Value;
        return new ScoreTargetSiteGrimsonResult(true, new GrimsonContextScoreDto(
            v.Type.ToString(), v.SiteTypeContribution, v.ThreePrimePairingContribution, v.LocalAuContribution,
            v.PositionContribution, v.ThreePrimePairingScore, v.LocalAuFraction, v.ContextScore));
    }

    /// <summary>TargetScan 7 context++ score (Agarwal et al. 2015) of a canonical seed site.</summary>
    [McpServerTool(Name = "score_target_site_context_plus_plus", Title = "miRNA — Score Target Site (context++)", ReadOnly = true)]
    [Description("TargetScan 7 context++ score (Agarwal et al. 2015, eLife 4:e05005; port of targetscan_70_context_scores.pl) of an 8mer, 7mer-m8, 7mer-A1 or 6mer site in a 3'UTR (siteStart/siteEnd 0-based inclusive; lengths 8/7/7/6). Returns every feature contribution (intercept, local AU, sRNA1/8, Site8, SA, 3P, Min_dist, Len_3UTR, Off6m, SPS, TA, Len_ORF, ORF8m, PCT), the unrounded partial sum and contextScore as TargetScan reports it (3-decimal contributions, per-site-type ceiling; null when the site starts within 14 nt of the CDS). SPS and TA default to the bundled TargetScan 7.0 seed-region table. All optional inputs are sparse; a feature without data contributes 0 and is listed in omittedFeatures. orfSequence derives Len_ORF and ORF8m as targetscan_count_8mers.pl does (orfLength/orf8mCount override it); plfoldUnpairedProbability is the exact RNAplfold -L 40 -W 80 -u 20 value (row utrStart+7, column 14) and replaces the library's single-window SA estimate; a score with any omitted feature is a guarded limitation (MIRNA-TARGET-001) and throws outside Permissive mode unless every optional feature is supplied; conservation (a Newick tree or a TargetScan bin tree 1-10, plus the species with the site, optional PCT sigmoid row) gives the Friedman 2009 branch length and PCT.")]
    public static ScoreTargetSiteContextPlusPlusResult ScoreTargetSiteContextPlusPlus(
        [Description("3'UTR nucleotide sequence (DNA T→U is applied internally)")] string mRnaSequence,
        [Description("miRNA nucleotide sequence (>= 8 nt)")] string miRnaSequence,
        [Description("0-based inclusive site start")] int siteStart,
        [Description("0-based inclusive site end (span length 8 for Seed8mer, 7 for Seed7merM8/Seed7merA1, 6 for Seed6mer)")] int siteEnd,
        [Description("Site type: Seed8mer, Seed7merM8, Seed7merA1 or Seed6mer")] string siteType,
        [Description("Optional miRNA name (label only)")] string? miRnaName = null,
        [Description("Optional seed-pairing stability SPS (kcal/mol); overrides the bundled TargetScan value")] double? sps = null,
        [Description("Optional target-site abundance TA (log10 site count; see compute_ta_3utr); overrides the bundled TargetScan value")] double? ta = null,
        [Description("Optional ORF (CDS) length in nt (> 0); overrides the length derived from orfSequence")] double? orfLength = null,
        [Description("Optional count of cognate 8mer sites in the ORF (>= 0); overrides the count derived from orfSequence")] int? orf8mCount = null,
        [Description("Optional ORF sequence (RNA/DNA, any case, '-'/'.' gaps accepted); Len_ORF and ORF8m are derived as targetscan_count_8mers.pl does")] string? orfSequence = null,
        [Description("Optional exact RNAplfold unpaired probability (0-1) of the site's 14-nt accessibility window; 0 or omitted ⇒ library estimate")] double? plfoldUnpairedProbability = null,
        [Description("Optional multi-species conservation input: newick OR targetScanBin (1-10), plus speciesWithSite, optional sigmoidParameters {b0,b1,b2,b3}")] ContextConservationInputDto? conservation = null)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            throw new ArgumentException("mRNA sequence cannot be null or empty", nameof(mRnaSequence));
        if (string.IsNullOrEmpty(miRnaSequence) || miRnaSequence.Length < 8)
            throw new ArgumentException("miRNA sequence must have at least 8 nt (seed region nt 2-8)", nameof(miRnaSequence));
        if (siteStart < 0 || siteEnd >= mRnaSequence.Length || siteStart > siteEnd)
            throw new ArgumentOutOfRangeException(nameof(siteStart),
                $"Site indices [{siteStart}, {siteEnd}] are out of range for sequence of length {mRnaSequence.Length}");
        if (!TryParseEnum<MiRnaAnalyzer.TargetSiteType>(siteType, out var type)
            || type is not (MiRnaAnalyzer.TargetSiteType.Seed8mer or MiRnaAnalyzer.TargetSiteType.Seed7merM8
                or MiRnaAnalyzer.TargetSiteType.Seed7merA1 or MiRnaAnalyzer.TargetSiteType.Seed6mer))
            throw new ArgumentException(
                $"Unknown or unsupported siteType '{siteType}'. Valid values: Seed8mer, Seed7merM8, Seed7merA1, Seed6mer", nameof(siteType));
        int length = type switch
        {
            MiRnaAnalyzer.TargetSiteType.Seed8mer => 8,
            MiRnaAnalyzer.TargetSiteType.Seed6mer => 6,
            _ => 7
        };
        if (siteEnd - siteStart + 1 != length)
            throw new ArgumentException($"A {type} site spans {length} nt; got [{siteStart}, {siteEnd}]", nameof(siteEnd));
        if (sps is { } s && !double.IsFinite(s))
            throw new ArgumentOutOfRangeException(nameof(sps), "sps must be finite");
        if (ta is { } a && !double.IsFinite(a))
            throw new ArgumentOutOfRangeException(nameof(ta), "ta must be finite");
        if (orfLength is { } ol && (!double.IsFinite(ol) || ol <= 0))
            throw new ArgumentOutOfRangeException(nameof(orfLength), "orfLength must be a finite value > 0");
        if (orf8mCount is < 0)
            throw new ArgumentOutOfRangeException(nameof(orf8mCount), "orf8mCount must be >= 0");
        if (plfoldUnpairedProbability is { } pf && (!double.IsFinite(pf) || pf < 0 || pf > 1))
            throw new ArgumentOutOfRangeException(nameof(plfoldUnpairedProbability), "plfoldUnpairedProbability must be in [0, 1]");

        MiRnaAnalyzer.PctConservation? pct = null;
        if (conservation is not null)
        {
            bool hasNewick = !string.IsNullOrWhiteSpace(conservation.Newick);
            if (hasNewick == conservation.TargetScanBin.HasValue)
                throw new ArgumentException("conservation requires exactly one of newick or targetScanBin", nameof(conservation));
            if (conservation.TargetScanBin is < 1 or > 10)
                throw new ArgumentOutOfRangeException(nameof(conservation), "targetScanBin must be 1-10");
            if (conservation.SpeciesWithSite is null || conservation.SpeciesWithSite.Count == 0)
                throw new ArgumentException("conservation.speciesWithSite cannot be null or empty", nameof(conservation));
            Genomics.Phylogenetics.PhylogeneticAnalyzer.PhyloNode tree;
            try
            {
                tree = hasNewick
                    ? Genomics.Phylogenetics.PhylogeneticAnalyzer.ParseNewick(conservation.Newick!)
                    : MiRnaAnalyzer.GetTargetScanBinTree(conservation.TargetScanBin!.Value);
            }
            catch (FormatException ex)
            {
                throw new ArgumentException($"Invalid Newick tree: {ex.Message}", nameof(conservation), ex);
            }
            MiRnaAnalyzer.PctSigmoidParameters? sig = conservation.SigmoidParameters is { } q
                ? new MiRnaAnalyzer.PctSigmoidParameters(q.B0, q.B1, q.B2, q.B3)
                : null;
            pct = new MiRnaAnalyzer.PctConservation(tree, conservation.SpeciesWithSite.ToList(), sig);
        }

        var miRna = MiRnaAnalyzer.CreateMiRna(miRnaName ?? "miRNA", miRnaSequence);
        var site = new MiRnaAnalyzer.TargetSite(siteStart, siteEnd,
            Seqeron.Genomics.Infrastructure.SequenceCase.ToUpperRna(mRnaSequence.Substring(siteStart, length)),
            miRna.Name, type, length, 0.0, 0.0, "");
        var inputs = new MiRnaAnalyzer.ContextPlusPlusInputs(
            Sps: sps, Ta: ta, OrfLength: orfLength, Orf8mCount: orf8mCount, Conservation: pct,
            PlfoldUnpairedProbability: plfoldUnpairedProbability, OrfSequence: orfSequence);

        var c = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mRnaSequence, miRna, site, inputs);
        return new ScoreTargetSiteContextPlusPlusResult(new ContextPlusPlusScoreDto(
            c.SiteType.ToString(), c.Intercept, c.LocalAuContribution, c.SRna1Contribution, c.SRna8Contribution,
            c.Site8Contribution, c.SaContribution, c.ThreePrimePairingContribution, c.MinDistContribution,
            c.Len3UtrContribution, c.Off6mContribution, c.SpsContribution, c.TaContribution, c.LenOrfContribution,
            c.Orf8mContribution, c.PctContribution, c.BranchLengthScore, c.Pct, c.PctIsNa, c.ContextScorePartial,
            c.ContextScore, c.OmittedFeatures.ToList()));
    }

    /// <summary>TA_3UTR (target-site abundance) of a miRNA seed over a set of 3'UTRs.</summary>
    [McpServerTool(Name = "compute_ta_3utr", Title = "miRNA — Compute TA_3UTR", ReadOnly = true)]
    [Description("TargetScan context++ feature TA_3UTR (Garcia et al. 2011; Agarwal et al. 2015): TA = log10(N), N = total number of non-overlapping canonical 8mer + 7mer-m8 + 7mer-A1 sites of the miRNA seed (nt 2-8) across the supplied 3'UTRs (bare 6mers and offset 6mers are not counted; scan is left-to-right greedy). Returns siteCount (N) and ta (0 when N = 0). Feed ta into score_target_site_context_plus_plus.")]
    public static ComputeTa3UtrResult ComputeTa3Utr(
        [Description("miRNA nucleotide sequence (>= 8 nt; seed = nt 2-8)")] string miRnaSequence,
        [Description("3'UTR sequences (RNA or DNA, any case) that define the transcriptome")] IReadOnlyList<string> threePrimeUtrs)
    {
        if (string.IsNullOrEmpty(miRnaSequence) || miRnaSequence.Length < 8)
            throw new ArgumentException("miRNA sequence must have at least 8 nt (seed region nt 2-8)", nameof(miRnaSequence));
        if (threePrimeUtrs is null)
            throw new ArgumentException("threePrimeUtrs cannot be null", nameof(threePrimeUtrs));

        var miRna = MiRnaAnalyzer.CreateMiRna("miRNA", miRnaSequence);
        long n = MiRnaAnalyzer.CountSeedSites3Utr(miRna, threePrimeUtrs);
        double ta = MiRnaAnalyzer.ComputeTa3Utr(miRna, threePrimeUtrs);
        return new ComputeTa3UtrResult(n, ta, threePrimeUtrs.Count);
    }

    /// <summary>Find and group seed-match sites across a multi-species 3'UTR alignment (TargetScan 7.0).</summary>
    [McpServerTool(Name = "find_targetscan_site_groups", Title = "miRNA — Find TargetScan Site Groups", ReadOnly = true)]
    [Description("Port of targetscan_70.pl: finds 8mer / 7mer-m8 / 7mer-A1 / 6mer seed-match sites of a miRNA family in every species of a multi-species 3'UTR alignment (alignedUtrs: species identifier → aligned UTR with '-' gaps; lower-case letters flag ORF overlap; the first 14 nt of each UTR are masked as TargetScan does) and groups sites that overlap across species. Coordinates are 1-based (msaStart/msaEnd alignment columns, utrStart/utrEnd ungapped UTR). Each group lists its sites, group type (e.g. '7mer-1a+8mer-1a'), species, and the species per site type that the branch-length/PCT step uses. speciesWithMiRna restricts the searched species (TargetScan 7: those annotated for the family); omit to search all.")]
    public static FindTargetScanSiteGroupsResult FindTargetScanSiteGroups(
        [Description("7-nt seed region (miRNA nt 2-8) or a mature miRNA of >= 8 nt")] string miRnaSeedOrMature,
        [Description("Species identifier → aligned 3'UTR (gaps '-'; lower case = ORF overlap)")] Dictionary<string, string> alignedUtrs,
        [Description("Optional species in which the miRNA family is annotated (only these are searched); default all")] IReadOnlyList<string>? speciesWithMiRna = null)
    {
        ValidateTargetScanInputs(miRnaSeedOrMature, alignedUtrs);
        var groups = MiRnaAnalyzer.FindTargetScanSiteGroups(miRnaSeedOrMature, alignedUtrs, speciesWithMiRna);
        return new FindTargetScanSiteGroupsResult(groups.Select(ToDto).ToList());
    }

    /// <summary>TargetScan 7.0 branch length, PCT and conservation call per site of a multi-species 3'UTR alignment.</summary>
    [McpServerTool(Name = "score_targetscan_conservation", Title = "miRNA — Score TargetScan Conservation", ReadOnly = true)]
    [Description("End-to-end TargetScan 7.0 conservation of a miRNA family's sites in a multi-species 3'UTR alignment: site groups (targetscan_70.pl) → UTR bin 1-10 from the median per-position branch length (targetscan_70_BL_bins.pl) → per-site branch length on the bin tree, PCT with the bundled TargetScan 7.0 sigmoid rows and the 'Conserved' call (targetscan_70_BL_PCT.pl). pct is null for TargetScan's 'NA' (seed not one of the 111 PCT families, or ORF-overlapping site). utrBin overrides the computed bin (medianBranchLength is then null); referenceSpecies defaults to 9606; thresholds is TargetScanHuman7 (default) or TargetScanMouse7.")]
    public static ScoreTargetScanConservationResult ScoreTargetScanConservation(
        [Description("7-nt seed region (miRNA nt 2-8) or a mature miRNA of >= 8 nt")] string miRnaSeedOrMature,
        [Description("Species identifier → aligned 3'UTR (gaps '-'; lower case = ORF overlap)")] Dictionary<string, string> alignedUtrs,
        [Description("Optional species in which the miRNA family is annotated (only these are searched); default all")] IReadOnlyList<string>? speciesWithMiRna = null,
        [Description("Reference species for the UTR bin (default 9606)")] string referenceSpecies = "9606",
        [Description("Optional TargetScan bin (1-10) to use instead of the computed one")] int? utrBin = null,
        [Description("Conservation thresholds: TargetScanHuman7 (default) or TargetScanMouse7")] string thresholds = "TargetScanHuman7")
    {
        ValidateTargetScanInputs(miRnaSeedOrMature, alignedUtrs);
        if (string.IsNullOrEmpty(referenceSpecies))
            throw new ArgumentException("referenceSpecies cannot be null or empty", nameof(referenceSpecies));
        if (utrBin is < 1 or > 10)
            throw new ArgumentOutOfRangeException(nameof(utrBin), "utrBin must be 1-10");
        MiRnaAnalyzer.BranchLengthConservationThresholds th;
        if (string.Equals(thresholds, "TargetScanHuman7", StringComparison.OrdinalIgnoreCase))
            th = MiRnaAnalyzer.BranchLengthConservationThresholds.TargetScanHuman7;
        else if (string.Equals(thresholds, "TargetScanMouse7", StringComparison.OrdinalIgnoreCase))
            th = MiRnaAnalyzer.BranchLengthConservationThresholds.TargetScanMouse7;
        else
            throw new ArgumentException($"Unknown thresholds '{thresholds}'. Valid values: TargetScanHuman7, TargetScanMouse7", nameof(thresholds));

        var r = MiRnaAnalyzer.ScoreTargetScanConservation(
            miRnaSeedOrMature, alignedUtrs, speciesWithMiRna, referenceSpecies, utrBin, th);
        return new ScoreTargetScanConservationResult(
            r.UtrBin.Bin,
            double.IsNaN(r.UtrBin.MedianBranchLength) ? null : r.UtrBin.MedianBranchLength,
            r.Groups.Select(ToDto).ToList(),
            r.Sites.Select(s => new TargetScanSiteConservationDto(
                ToDto(s.Site), s.GroupNumber, s.SpeciesWithSite.ToList(), s.BranchLength, s.Pct, s.Conserved)).ToList());
    }

    private static void ValidateTargetScanInputs(string miRnaSeedOrMature, Dictionary<string, string> alignedUtrs)
    {
        if (string.IsNullOrEmpty(miRnaSeedOrMature) || miRnaSeedOrMature.Length < 7)
            throw new ArgumentException("Expected the 7-nt seed region or a mature miRNA of >= 8 nt", nameof(miRnaSeedOrMature));
        if (alignedUtrs is null || alignedUtrs.Count == 0)
            throw new ArgumentException("alignedUtrs cannot be null or empty", nameof(alignedUtrs));
        foreach (var (species, utr) in alignedUtrs)
        {
            if (string.IsNullOrEmpty(species))
                throw new ArgumentException("A species key is empty", nameof(alignedUtrs));
            if (utr is null)
                throw new ArgumentException($"The aligned UTR of species '{species}' is null", nameof(alignedUtrs));
        }
    }

    private static TargetScanAlignedSiteDto ToDto(MiRnaAnalyzer.TargetScanAlignedSite s) =>
        new(s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd, s.SiteType.ToString(), s.OrfOverlap);

    private static TargetScanSiteGroupDto ToDto(MiRnaAnalyzer.TargetScanSiteGroup g) =>
        new(g.GroupNumber, g.GroupType, g.SiteTypes.Select(t => t.ToString()).ToList(), g.Species.ToList(),
            g.Sites.Select(ToDto).ToList(),
            Enum.GetValues<MiRnaAnalyzer.TargetSiteType>()
                .Select(t => new TargetScanSiteTypeSpeciesDto(t.ToString(), g.GetSpeciesWithSiteType(t).ToList()))
                .Where(x => x.Species.Count > 0).ToList());

    // Enum.TryParse accepts numeric strings and undefined values; accept defined names only (case-insensitive).
    private static bool TryParseEnum<T>(string? value, out T result) where T : struct, Enum =>
        Enum.TryParse(value, ignoreCase: true, out result) && Enum.IsDefined(result) && !int.TryParse(value, out _);

    /// <summary>Compute AU content and positional context score around a target site.</summary>
    [McpServerTool(Name = "analyze_target_context", Title = "miRNA — Analyze Target Context", ReadOnly = true)]
    [Description("Compute the local AU content and the UTR-context part of the TargetScan 5.0 / Grimson 2007 context score (local-AU + position contributions of a 7mer-m8 site, log2 units, negative = repressed; 0 when the site starts within the first 15 nt) around a miRNA target site within a 3'UTR.")]
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

    /// <summary>Site accessibility = McCaskill probability that the site is unpaired (Turner 2004).</summary>
    [McpServerTool(Name = "site_accessibility", Title = "miRNA — Site Accessibility", ReadOnly = true)]
    [Description("miRNA target site accessibility: McCaskill partition-function probability that the whole site is unpaired (Turner 2004 parameters), folding a window of max(80, site length) nt centred on the site (RNAplfold -W 80 context). Range 0–1; higher ⇒ more accessible.")]
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
        ArgumentNullException.ThrowIfNull(miRnas);

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
        ArgumentNullException.ThrowIfNull(query);
        ArgumentNullException.ThrowIfNull(database);

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
        [Description("Seed nucleotide sequence")] string seedSequence)
    {
        if (string.IsNullOrEmpty(seedSequence))
            throw new ArgumentException("Seed sequence cannot be null or empty", nameof(seedSequence));

        var variants = MiRnaAnalyzer.GenerateSeedVariants(seedSequence).ToList();
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
    [Description("Find candidate 5' (donor) splice sites: GU/GT (optional GC) scored with the Shapiro & Senapathy (1987) 5' ss score over positions -2..+6 (published 0-100 score / 100) or with MaxEntScan score5 (bits); optional U12-type /RTATCCTY donors. Position = index of the G of GT (first intron nt).")]
    public static FindDonorSitesResult FindDonorSites(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Minimum score for a site to be reported: [0,1] for shapiro-senapathy (and U12 match fraction), bits for maxent")] double minScore = 0.5,
        [Description("Also report GC donors and U12-type (/RTATCCTY: AT-AC or GT-AG) donors")] bool includeNonCanonical = false,
        [Description("Scoring model for GT/GC donors: 'shapiro-senapathy' (default) or 'maxent' (MaxEntScan score5, 9-nt window -3..+6)")] string model = "shapiro-senapathy",
        [Description("Shapiro & Senapathy (1987) Table 1 subtable: primate (default), total, rodent, othermammal, othervertebrate, invertebrate, plant, viral")] string category = "primate")
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var scoringModel = (model ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "shapiro-senapathy" or "shapirosenapathy" or "ss" => SpliceSitePredictor.DonorScoringModel.ShapiroSenapathy,
            "maxent" or "maxentscan" => SpliceSitePredictor.DonorScoringModel.MaxEntScan,
            _ => throw new ArgumentException(
                $"Unknown model '{model}'; expected 'shapiro-senapathy' or 'maxent'.", nameof(model))
        };

        var table = ParseShapiroSenapathyCategory(category);

        var sites = SpliceSitePredictor
            .FindDonorSites(sequence, minScore, includeNonCanonical, scoringModel, table)
            .Select(ToDto)
            .ToList();
        return new FindDonorSitesResult(sites);
    }

    private static SpliceSitePredictor.ShapiroSenapathyCategory ParseShapiroSenapathyCategory(string category)
    {
        if (string.IsNullOrWhiteSpace(category) || category.Trim().Any(char.IsDigit) ||
            !Enum.TryParse<SpliceSitePredictor.ShapiroSenapathyCategory>(category.Trim(), ignoreCase: true, out var table))
            throw new ArgumentException(
                $"Unknown category '{category}'; expected primate, total, rodent, othermammal, othervertebrate, invertebrate, plant or viral.",
                nameof(category));
        return table;
    }

    /// <summary>Find candidate 3' (acceptor) splice sites.</summary>
    [McpServerTool(Name = "find_acceptor_sites", Title = "Splice — Find Acceptor Sites", ReadOnly = true)]
    [Description("Find candidate 3' (acceptor) splice sites: AG scored with the Shapiro & Senapathy (1987) 3' ss score over positions -14..+1 (pyrimidine region best 8 of 10 + CAGG; published 0-100 score / 100) or with MaxEntScan score3 (bits, 23-nt window -20..+3); optional U12-type acceptors (AC or AG terminus scored by the U12 branch site UUCCUUAAC 8-21 nt upstream). Position = index of the G of AG (last intron nt).")]
    public static FindAcceptorSitesResult FindAcceptorSites(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Minimum score for a site to be reported: [0,1] for shapiro-senapathy (and U12 branch-site match fraction), bits for maxent")] double minScore = 0.5,
        [Description("Also report U12-type acceptors (AC or AG terminus with a UUCCUUAAC branch site 8-21 nt upstream)")] bool includeNonCanonical = false,
        [Description("Scoring model for AG acceptors: 'shapiro-senapathy' (default) or 'maxent' (MaxEntScan score3, 23-nt window -20..+3)")] string model = "shapiro-senapathy",
        [Description("Shapiro & Senapathy (1987) Table 1 subtable: primate (default), total, rodent, othermammal, othervertebrate, invertebrate, plant, viral")] string category = "primate")
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var scoringModel = (model ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "shapiro-senapathy" or "shapirosenapathy" or "ss" => SpliceSitePredictor.AcceptorScoringModel.ShapiroSenapathy,
            "maxent" or "maxentscan" => SpliceSitePredictor.AcceptorScoringModel.MaxEntScan,
            _ => throw new ArgumentException(
                $"Unknown model '{model}'; expected 'shapiro-senapathy' or 'maxent'.", nameof(model))
        };

        var table = ParseShapiroSenapathyCategory(category);

        var sites = SpliceSitePredictor
            .FindAcceptorSites(sequence, minScore, includeNonCanonical, scoringModel, table)
            .Select(ToDto)
            .ToList();
        return new FindAcceptorSitesResult(sites);
    }

    /// <summary>Find branch-point candidates (human yUnAy or yeast UACUAAC).</summary>
    [McpServerTool(Name = "find_branch_points", Title = "Splice — Find Branch Points", ReadOnly = true)]
    [Description("Find branch-point candidates: human yUnAy (Gao et al. 2008; conservation-weighted match of y@-3, U@-2, A@0, y@+1, score [0,1]) or yeast UACUAAC (match fraction). Position = index of the branch adenosine position; motifs with non-ACGTU characters are skipped.")]
    public static FindBranchPointsResult FindBranchPoints(
        [Description("RNA/DNA sequence to scan")] string sequence,
        [Description("Smallest motif start index scanned (inclusive)")] int searchStart = 0,
        [Description("Largest motif start index scanned (inclusive); -1 scans to the end of the sequence")] int searchEnd = -1,
        [Description("Minimum score [0,1] for a branch point to be reported")] double minScore = 0.5,
        [Description("Branch-point consensus: 'human' (yUnAy, default) or 'yeast' (UACUAAC)")] string model = "human")
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var bpModel = (model ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "human" or "yunay" => SpliceSitePredictor.BranchPointModel.HumanYUnAy,
            "yeast" or "uacuaac" or "tactaac" => SpliceSitePredictor.BranchPointModel.YeastUACUAAC,
            _ => throw new ArgumentException(
                $"Unknown model '{model}'; expected 'human' or 'yeast'.", nameof(model))
        };

        var sites = SpliceSitePredictor
            .FindBranchPoints(sequence, searchStart, searchEnd, minScore, bpModel)
            .Select(ToDto)
            .ToList();
        return new FindBranchPointsResult(sites);
    }

    /// <summary>Predict introns by pairing donor/acceptor sites.</summary>
    [McpServerTool(Name = "predict_introns", Title = "Splice — Predict Introns", ReadOnly = true)]
    [Description("Predict introns (0-based inclusive coordinates) by pairing donor and acceptor splice sites (GT-AG / GC-AG = U2, AT-AC or GT-AG with U12 donor + UUCCUUAAC branch site = U12) and locating Gao 2008 yUnAy branch points; score = mean(donor, acceptor, branch point) for U2.")]
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
    [Description("Predict exon/intron gene structure. model=greedy (default): splice-motif model, not coding-aware (greedy non-overlapping intron selection by score; exons >= minExonLength; exons + introns tile the sequence). model=splice_dp: exact dynamic program (GENSCAN/GeneSplicer style) maximising MaxEntScan 5'+3' splice scores (bits, minus intronPenalty per intron) subject to an open reading frame from the start codon (startPosition, default first AUG) to the first in-frame stop with no premature stop; returns the CDS and scores, errors when no open reading frame is reachable.")]
    public static PredictGeneStructureResult PredictGeneStructure(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum exon length (nt); for splice_dp the effective minimum is 2 (default 30 for greedy, 3 for splice_dp when omitted)")] int? minExonLength = null,
        [Description("Minimum intron length (nt)")] int minIntronLength = 60,
        [Description("greedy only: minimum combined score for an intron to be considered")] double minScore = 0.5,
        [Description("Model: greedy (default) or splice_dp")] string model = "greedy",
        [Description("splice_dp: 0-based first coding nucleotide (start codon); omitted = first AUG")] int? startPosition = null,
        [Description("splice_dp: maximum intron length (nt)")] int maxIntronLength = 100000,
        [Description("splice_dp: minimum MaxEntScan score5 (bits) of a donor")] double minDonorScore = 0.0,
        [Description("splice_dp: minimum MaxEntScan score3 (bits) of an acceptor")] double minAcceptorScore = 0.0,
        [Description("splice_dp: bits subtracted per intron (default 12)")] double intronPenalty = 12.0)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        switch (model?.Trim().ToLowerInvariant())
        {
            case "greedy":
            {
                var gs = SpliceSitePredictor.PredictGeneStructure(sequence, minExonLength ?? 30, minIntronLength, minScore);
                var exons = gs.Exons.Select(ToDto).ToList();
                var introns = gs.Introns.Select(ToDto).ToList();
                return new PredictGeneStructureResult(exons, introns, gs.SplicedSequence, gs.OverallScore);
            }
            case "splice_dp":
            {
                var dp = SpliceSitePredictor.PredictGeneStructureDp(sequence, startPosition, minExonLength ?? 3,
                    minIntronLength, maxIntronLength, minDonorScore, minAcceptorScore, null, intronPenalty)
                    ?? throw new ArgumentException(
                        "No open reading frame from the start codon to a stop codon is reachable with the splice candidates.",
                        nameof(sequence));
                var gs = dp.Structure;
                return new PredictGeneStructureResult(
                    gs.Exons.Select(ToDto).ToList(), gs.Introns.Select(ToDto).ToList(), gs.SplicedSequence, gs.OverallScore,
                    dp.CdsStart, dp.CodingSequence, dp.IntronScore);
            }
            default:
                throw new ArgumentException($"Unknown model '{model}'; use greedy or splice_dp.", nameof(model));
        }
    }

    /// <summary>Detect candidate alternative splicing patterns.</summary>
    [McpServerTool(Name = "detect_alternative_splicing", Title = "Splice — Detect Alternative Splicing", ReadOnly = true)]
    [Description("Detect candidate alternative splicing events relative to the predicted gene structure with SUPPA/rMATS definitions: ExonSkipping (SE), Alt5SS (A5SS), Alt3SS (A3SS), MutuallyExclusiveExons (MXE).")]
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

    /// <summary>Find retained-intron candidates by sequence-level features (open reading frame; optional length/weak-site/GC filters).</summary>
    [McpServerTool(Name = "find_retained_intron_candidates", Title = "Splice — Find Retained Intron Candidates", ReadOnly = true)]
    [Description("Predicted introns that keep the reading frame open (length multiple of 3, no in-frame stop in some phase; Jacob & Smith 2017). Optional filters for the Braunschweig 2014 features (shorter, weaker splice sites, higher GC) have no published cut-offs, so they are off unless set. Retention itself needs transcript evidence.")]
    public static FindRetainedIntronCandidatesResult FindRetainedIntronCandidates(
        [Description("RNA/DNA sequence to analyze")] string sequence,
        [Description("Minimum combined score for an intron candidate")] double minScore = 0.5,
        [Description("Optional maximum intron length (nt); omitted = no limit")] int? maxIntronLength = null,
        [Description("Optional: keep only introns with combined score below this value (weak sites)")] double? maxIntronScore = null,
        [Description("Optional minimum G+C fraction (0-1) of the intron")] double? minGcFraction = null,
        [Description("Require an open reading frame (no in-frame stop in some phase)")] bool requireOpenReadingFrame = true,
        [Description("Require intron length divisible by 3")] bool requireLengthMultipleOf3 = true)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var criteria = new SpliceSitePredictor.RetainedIntronCriteria(
            MaxIntronLength: maxIntronLength,
            MaxIntronScore: maxIntronScore,
            MinGcFraction: minGcFraction,
            RequireOpenReadingFrame: requireOpenReadingFrame,
            RequireLengthMultipleOf3: requireLengthMultipleOf3);
        var introns = SpliceSitePredictor
            .FindRetainedIntronCandidates(sequence, criteria, minScore)
            .Select(ToDto)
            .ToList();
        return new FindRetainedIntronCandidatesResult(introns);
    }

    /// <summary>Score a splice motif: Donor = MaxEntScan score5, Acceptor = MaxEntScan score3 (Yeo &amp; Burge 2004).</summary>
    [McpServerTool(Name = "maxent_score", Title = "Splice — MaxEnt Score", ReadOnly = true)]
    [Description("Score a splice motif. Donor: MaxEntScan score5 (Yeo & Burge 2004) in bits on the 9-nt window (3 exon + 6 intron nt, GT at positions 4-5), e.g. CAGGTAAGT = 10.86. Acceptor: MaxEntScan score3 in bits on the 23-nt window (20 intron + 3 exon nt, AG at positions 19-20), e.g. TTCCAAACGAACTTTTGTAGGGA = 2.89.")]
    public static MaxEntScoreResult MaxentScore(
        [Description("Splice-site motif sequence (RNA or DNA); Donor exactly 9 nt, Acceptor exactly 23 nt, of A/C/G/T(/U)")] string motif,
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
    [Description("Whether a position lies in an open reading frame opened by an upstream AUG with no in-frame stop before the position's codon (stop codon itself included).")]
    public static IsWithinCodingRegionResult IsWithinCodingRegion(
        [Description("DNA/RNA sequence")] string sequence,
        [Description("Zero-based position to test")] int position,
        [Description("Position of the queried nucleotide inside its codon (0, 1 or 2)")] int frame = 0)
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
        ArgumentNullException.ThrowIfNull(readPairs);

        var input = readPairs.Select(p =>
            (p.ReadId, p.Chr1, p.Pos1, p.Strand1, p.Chr2, p.Pos2, p.Strand2, p.InsertSize));
        var pairs = StructuralVariantAnalyzer
            .FindDiscordantPairs(input, expectedInsertSize, insertSizeStdDev, maxInsertSize: maxInsertSize)
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
        ArgumentNullException.ThrowIfNull(discordantPairs);

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
        ArgumentNullException.ThrowIfNull(alignments);

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
        ArgumentNullException.ThrowIfNull(splitReads);

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
        ArgumentNullException.ThrowIfNull(probes);

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
        ArgumentNullException.ThrowIfNull(segments);

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
        ArgumentNullException.ThrowIfNull(variants);

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
        ArgumentNullException.ThrowIfNull(variants);

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
        ArgumentNullException.ThrowIfNull(variants);
        ArgumentNullException.ThrowIfNull(genes);

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
        ArgumentNullException.ThrowIfNull(sv);
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
        ArgumentNullException.ThrowIfNull(splitReads);

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
        ArgumentNullException.ThrowIfNull(geneCounts);
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
        ArgumentNullException.ThrowIfNull(samples);
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
        ArgumentNullException.ThrowIfNull(values);

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
        ArgumentNullException.ThrowIfNull(expressionData);
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
        ArgumentNullException.ThrowIfNull(differentiallyExpressedGenes);
        ArgumentNullException.ThrowIfNull(pathways);
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
        ArgumentNullException.ThrowIfNull(rankedGenes);
        ArgumentNullException.ThrowIfNull(geneSet);
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
        ArgumentNullException.ThrowIfNull(exonData);

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
        ArgumentNullException.ThrowIfNull(splicingData);

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
        ArgumentNullException.ThrowIfNull(isoforms);

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
        ArgumentNullException.ThrowIfNull(isoformData);

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
        ArgumentNullException.ThrowIfNull(expression1);
        ArgumentNullException.ThrowIfNull(expression2);
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
        ArgumentNullException.ThrowIfNull(geneProfiles);

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
        ArgumentNullException.ThrowIfNull(geneProfiles);

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
        ArgumentNullException.ThrowIfNull(geneCounts);
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
        ArgumentNullException.ThrowIfNull(samples);
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
