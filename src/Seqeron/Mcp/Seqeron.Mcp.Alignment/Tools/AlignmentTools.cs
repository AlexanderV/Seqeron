using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

/// <summary>
/// MCP tools for sequence alignment, approximate pattern matching, and
/// genome assembly operations.
/// </summary>
[McpServerToolType]
public class AlignmentTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private AlignmentTools() { }

    #region SequenceAligner

    [McpServerTool(Name = "global_align", Title = "Pairwise — Global Alignment (Needleman–Wunsch)", ReadOnly = true)]
    [Description("Performs global pairwise alignment using the Needleman–Wunsch dynamic programming algorithm (aligns full length of both sequences end-to-end).")]
    public static AlignmentResultDto GlobalAlign(
        [Description("First sequence (DNA, uppercased internally).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open penalty.")] int gapOpen = -2,
        [Description("Gap-extend penalty.")] int gapExtend = -1)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, gapOpen, gapExtend);
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "local_align", Title = "Pairwise — Local Alignment (Smith–Waterman)", ReadOnly = true)]
    [Description("Performs local pairwise alignment using the Smith–Waterman algorithm (finds best-scoring substring alignment, with zero floor).")]
    public static AlignmentResultDto LocalAlign(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open penalty.")] int gapOpen = -2,
        [Description("Gap-extend penalty.")] int gapExtend = -1)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, gapOpen, gapExtend);
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.LocalAlign(sequence1, sequence2, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "semi_global_align", Title = "Pairwise — Semi-Global (Fitting) Alignment", ReadOnly = true)]
    [Description("Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2; useful for fitting a shorter query into a longer reference.")]
    public static AlignmentResultDto SemiGlobalAlign(
        [Description("Query sequence (typically shorter).")] string sequence1,
        [Description("Reference sequence (typically longer).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open penalty.")] int gapOpen = -2,
        [Description("Gap-extend penalty.")] int gapExtend = -1)
    {
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence1, out var d1))
            throw new ArgumentException("Invalid DNA in sequence1.", nameof(sequence1));
        if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequence2, out var d2))
            throw new ArgumentException("Invalid DNA in sequence2.", nameof(sequence2));

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, gapOpen, gapExtend);
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign(d1!, d2!, scoring);
        return ToAlignmentDto(r);
    }

    [McpServerTool(Name = "alignment_statistics", Title = "Alignment — Statistics", ReadOnly = true)]
    [Description("Computes match / mismatch / gap counts and percent identity, similarity, and gap percent for a previously produced pairwise alignment.")]
    public static AlignmentStatisticsDto AlignmentStatistics(
        [Description("Aligned representation of sequence1 (with '-' gaps).")] string alignedSequence1,
        [Description("Aligned representation of sequence2 (with '-' gaps).")] string alignedSequence2)
    {
        if (string.IsNullOrEmpty(alignedSequence1))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence1));
        if (string.IsNullOrEmpty(alignedSequence2))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence2));

        var a1 = alignedSequence1;
        var a2 = alignedSequence2;

        if (a1.Length != a2.Length)
            throw new ArgumentException("Aligned sequences must have equal length.", nameof(alignedSequence2));

        var ar = new global::Seqeron.Genomics.Infrastructure.AlignmentResult(
            a1, a2, 0,
            global::Seqeron.Genomics.Infrastructure.AlignmentType.Global,
            0, 0, 0, 0);

        var s = global::Seqeron.Genomics.Alignment.SequenceAligner.CalculateStatistics(ar);
        return new AlignmentStatisticsDto(
            s.Matches, s.Mismatches, s.Gaps, s.AlignmentLength,
            s.Identity, s.Similarity, s.GapPercent);
    }

    [McpServerTool(Name = "format_alignment", Title = "Alignment — Pretty Format", ReadOnly = true)]
    [Description("Renders a human-readable visual alignment (BLAST-style three-line block with '|' for matches, '.' for mismatches, ' ' for gaps), wrapped to a configurable line width.")]
    public static FormatAlignmentResult FormatAlignment(
        [Description("Aligned representation of sequence1.")] string alignedSequence1,
        [Description("Aligned representation of sequence2.")] string alignedSequence2,
        [Description("Maximum characters per line in the formatted output (>= 1).")] int lineWidth = 60)
    {
        if (string.IsNullOrEmpty(alignedSequence1))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence1));
        if (string.IsNullOrEmpty(alignedSequence2))
            throw new ArgumentException("Aligned sequence cannot be null or empty.", nameof(alignedSequence2));
        if (lineWidth < 1)
            throw new ArgumentOutOfRangeException(nameof(lineWidth), "lineWidth must be >= 1.");

        var ar = new global::Seqeron.Genomics.Infrastructure.AlignmentResult(
            alignedSequence1,
            alignedSequence2,
            0,
            global::Seqeron.Genomics.Infrastructure.AlignmentType.Global,
            0, 0, 0, 0);

        var formatted = global::Seqeron.Genomics.Alignment.SequenceAligner.FormatAlignment(ar, lineWidth);
        return new FormatAlignmentResult(formatted);
    }

    [McpServerTool(Name = "multiple_align", Title = "MSA — Anchor-Based Progressive Alignment", ReadOnly = true)]
    [Description("Anchor-based progressive multiple sequence alignment: picks a center sequence by 4-mer cosine similarity, builds a suffix tree on it, and reconciles per-sequence anchored alignments into a single MSA with consensus and sum-of-pairs score.")]
    public static MultipleAlignmentResultDto MultipleAlign(
        [Description("DNA sequences to align (at least one).")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap-open penalty.")] int gapOpen = -2,
        [Description("Gap-extend penalty.")] int gapExtend = -1)
    {
        if (sequences == null || sequences.Length == 0)
        {
            var e = global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult.Empty;
            return new MultipleAlignmentResultDto(e.AlignedSequences, e.Consensus, e.TotalScore);
        }

        var dnaList = new List<global::Seqeron.Genomics.Core.DnaSequence>(sequences.Length);
        for (int i = 0; i < sequences.Length; i++)
        {
            if (!global::Seqeron.Genomics.Core.DnaSequence.TryCreate(sequences[i], out var d))
                throw new ArgumentException($"Invalid DNA at sequences[{i}].", nameof(sequences));
            dnaList.Add(d!);
        }

        var scoring = new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, gapOpen, gapExtend);
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlign(dnaList, scoring);
        return new MultipleAlignmentResultDto(r.AlignedSequences, r.Consensus, r.TotalScore);
    }

    #endregion

    #region ApproximateMatcher

    [McpServerTool(Name = "find_with_mismatches", Title = "Approximate — Find with Mismatches (Hamming)", ReadOnly = true)]
    [Description("Finds all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). Reports each match position, matched substring, mismatch count, and 0-based mismatch positions.")]
    public static ApproximateMatchListResult FindWithMismatches(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern,
        [Description("Maximum number of allowed mismatches (>= 0).")] int maxMismatches)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));
        if (maxMismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches), "maxMismatches must be >= 0.");

        var items = new List<ApproximateMatchDto>();
        foreach (var m in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindWithMismatches(
            sequence, pattern, maxMismatches))
        {
            items.Add(ToApproximateMatchDto(m));
        }
        return new ApproximateMatchListResult(items.ToArray());
    }

    [McpServerTool(Name = "find_with_edits", Title = "Approximate — Find with Edits (Levenshtein)", ReadOnly = true)]
    [Description("Finds all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions) with variable-length sliding windows of length pattern.Length ± maxEdits.")]
    public static ApproximateMatchListResult FindWithEdits(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern,
        [Description("Maximum allowed edit distance (>= 0).")] int maxEdits)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));
        if (maxEdits < 0)
            throw new ArgumentOutOfRangeException(nameof(maxEdits), "maxEdits must be >= 0.");

        var items = new List<ApproximateMatchDto>();
        foreach (var m in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindWithEdits(
            sequence, pattern, maxEdits))
        {
            items.Add(ToApproximateMatchDto(m));
        }
        return new ApproximateMatchListResult(items.ToArray());
    }

    [McpServerTool(Name = "find_best_match", Title = "Approximate — Find Best (Minimum Hamming)", ReadOnly = true)]
    [Description("Returns the single best (minimum-Hamming-distance) fixed-length window of pattern inside sequence. Stops early on perfect (distance=0) match.")]
    public static FindBestMatchResult FindBestMatch(
        [Description("Sequence to search in.")] string sequence,
        [Description("Pattern to find.")] string pattern)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty.", nameof(pattern));

        var m = global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindBestMatch(
            sequence, pattern);
        if (m is null)
            return new FindBestMatchResult(null);
        return new FindBestMatchResult(ToApproximateMatchDto(m.Value));
    }

    [McpServerTool(Name = "frequent_kmers_with_mismatches", Title = "K-mers — Frequent with Mismatches", ReadOnly = true)]
    [Description("Finds the most-frequent k-mers within sequence allowing up to d mismatches (counts each k-mer plus all its DNA neighbors at Hamming distance ≤ d). Returns all k-mers tied for the maximum count.")]
    public static FrequentKmersResult FrequentKmersWithMismatches(
        [Description("Sequence to analyze.")] string sequence,
        [Description("K-mer length (> 0).")] int k,
        [Description("Maximum mismatches in neighborhood (>= 0).")] int d)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be > 0.");
        if (d < 0)
            throw new ArgumentOutOfRangeException(nameof(d), "d must be >= 0.");

        var items = new List<FrequentKmerItem>();
        foreach (var t in global::Seqeron.Genomics.Alignment.ApproximateMatcher.FindFrequentKmersWithMismatches(
            sequence, k, d))
        {
            items.Add(new FrequentKmerItem(t.Kmer, t.Count));
        }
        return new FrequentKmersResult(items.ToArray());
    }

    #endregion

    #region SequenceAssembler

    [McpServerTool(Name = "assemble_olc", Title = "Assembly — Overlap-Layout-Consensus", ReadOnly = true)]
    [Description("Assembles input reads into contigs using the overlap–layout–consensus approach: detects pairwise suffix–prefix overlaps above thresholds, builds the overlap graph, lays out contigs, and reports N50 / longest / total length.")]
    public static AssemblyResultDto AssembleOlc(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length to accept (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio for an overlap to be accepted (0.0-1.0).")] double minIdentity = 0.9,
        [Description("K-mer size (unused by OLC path; part of shared assembly parameters).")] int kmerSize = 31,
        [Description("Minimum contig length to keep (bp).")] int minContigLength = 100)
    {
        var p = new global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyParameters(
            minOverlap, minIdentity, kmerSize, minContigLength);
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.AssembleOLC(
            reads ?? Array.Empty<string>(), p);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "assemble_de_bruijn", Title = "Assembly — de Bruijn Graph", ReadOnly = true)]
    [Description("Assembles input reads using a de Bruijn graph: shreds reads into k-mers, builds the (k-1)-mer node graph, traces non-branching paths into contigs, then reports assembly statistics.")]
    public static AssemblyResultDto AssembleDeBruijn(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length (shared parameter; not the operative knob here).")] int minOverlap = 20,
        [Description("Minimum identity (shared parameter).")] double minIdentity = 0.9,
        [Description("K-mer size used to build the de Bruijn graph.")] int kmerSize = 31,
        [Description("Minimum contig length to keep (bp).")] int minContigLength = 100)
    {
        var p = new global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyParameters(
            minOverlap, minIdentity, kmerSize, minContigLength);
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.AssembleDeBruijn(
            reads ?? Array.Empty<string>(), p);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "find_all_overlaps", Title = "Assembly — All Pairwise Overlaps", ReadOnly = true)]
    [Description("Computes all suffix-of-i / prefix-of-j overlaps between read pairs above minOverlap length and minIdentity ratio.")]
    public static OverlapListResult FindAllOverlaps(
        [Description("Input sequence reads.")] string[] reads,
        [Description("Minimum overlap length (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio (0.0-1.0).")] double minIdentity = 0.9)
    {
        var src = global::Seqeron.Genomics.Alignment.SequenceAssembler.FindAllOverlaps(
            reads ?? Array.Empty<string>(), minOverlap, minIdentity);
        var arr = new OverlapDto[src.Count];
        for (int i = 0; i < src.Count; i++)
        {
            var o = src[i];
            arr[i] = new OverlapDto(o.ReadIndex1, o.ReadIndex2, o.OverlapLength, o.Position1, o.Position2);
        }
        return new OverlapListResult(arr);
    }

    [McpServerTool(Name = "find_overlap", Title = "Assembly — Pairwise Suffix–Prefix Overlap", ReadOnly = true)]
    [Description("Detects the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. Returns null if no qualifying overlap exists.")]
    public static FindOverlapResult FindOverlap(
        [Description("First sequence (suffix candidate).")] string sequence1,
        [Description("Second sequence (prefix candidate).")] string sequence2,
        [Description("Minimum overlap length (bp).")] int minOverlap = 20,
        [Description("Minimum identity ratio (0.0-1.0).")] double minIdentity = 0.9)
    {
        var t = global::Seqeron.Genomics.Alignment.SequenceAssembler.FindOverlap(
            sequence1 ?? string.Empty, sequence2 ?? string.Empty, minOverlap, minIdentity);
        if (t is null)
            return new FindOverlapResult(null);
        var v = t.Value;
        return new FindOverlapResult(new FindOverlapInfo(v.length, v.pos1, v.pos2));
    }

    [McpServerTool(Name = "sequence_identity", Title = "Sequence — Percent Identity (Equal Length)", ReadOnly = true)]
    [Description("Computes percent-identity (fraction in [0,1]) between two equal-length sequences via case-insensitive position-by-position comparison. Returns 0 for unequal-length inputs.")]
    public static IdentityResult SequenceIdentity(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence (must equal sequence1.Length).")] string sequence2)
    {
        var s1 = sequence1 ?? string.Empty;
        var s2 = sequence2 ?? string.Empty;
        var v = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateIdentity(s1, s2);
        return new IdentityResult(v);
    }

    [McpServerTool(Name = "assembly_stats", Title = "Assembly — N50 / Length Stats", ReadOnly = true)]
    [Description("Computes assembly statistics (N50, longest contig, total length, assembled-reads accounting) for a precomputed list of contigs.")]
    public static AssemblyResultDto AssemblyStats(
        [Description("Assembled contig sequences.")] string[] contigs,
        [Description("Total number of input reads (used to fill TotalReads).")] int totalReads)
    {
        var r = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateStats(
            contigs ?? Array.Empty<string>(), totalReads);
        return ToAssemblyDto(r);
    }

    [McpServerTool(Name = "merge_contigs", Title = "Assembly — Merge Contigs by Overlap", ReadOnly = true)]
    [Description("Concatenates two contigs, collapsing the specified suffix/prefix overlap of overlapLength bases. Invalid overlapLength falls back to plain concatenation.")]
    public static MergedContigResult MergeContigs(
        [Description("First contig (suffix donor).")] string contig1,
        [Description("Second contig (prefix donor).")] string contig2,
        [Description("Length of overlap to collapse.")] int overlapLength)
    {
        var s = global::Seqeron.Genomics.Alignment.SequenceAssembler.MergeContigs(
            contig1 ?? string.Empty, contig2 ?? string.Empty, overlapLength);
        return new MergedContigResult(s);
    }

    [McpServerTool(Name = "scaffold_contigs", Title = "Assembly — Scaffold by Paired Links", ReadOnly = true)]
    [Description("Builds scaffolds by joining contigs using paired-end link records, inserting gapCharacter between linked contigs to span the indicated gap size.")]
    public static ScaffoldsResult ScaffoldContigs(
        [Description("Contig sequences.")] string[] contigs,
        [Description("Paired-end links between contigs (by index).")] ContigLinkInput[] links,
        [Description("Single-character gap filler (e.g. \"N\").")] string gapCharacter = "N")
    {
        if (string.IsNullOrEmpty(gapCharacter) || gapCharacter.Length != 1)
            throw new ArgumentException("gapCharacter must be exactly one character.", nameof(gapCharacter));

        var srcLinks = new List<(int contig1, int contig2, int gapSize)>(links?.Length ?? 0);
        if (links != null)
        {
            foreach (var l in links)
                srcLinks.Add((l.Contig1, l.Contig2, l.GapSize));
        }

        var scaffolds = global::Seqeron.Genomics.Alignment.SequenceAssembler.Scaffold(
            contigs ?? Array.Empty<string>(), srcLinks, gapCharacter[0]);

        var arr = new string[scaffolds.Count];
        for (int i = 0; i < scaffolds.Count; i++)
            arr[i] = scaffolds[i];
        return new ScaffoldsResult(arr);
    }

    [McpServerTool(Name = "calculate_coverage", Title = "Assembly — Per-Base Coverage Depth", ReadOnly = true)]
    [Description("Maps each read to its best-matching position on the reference (≥ minOverlap matching bases) and returns per-base coverage depth as an integer array of length reference.Length.")]
    public static CoverageResult CalculateCoverage(
        [Description("Reference sequence.")] string reference,
        [Description("Reads to map onto reference.")] string[] reads,
        [Description("Minimum matching bases required to map a read.")] int minOverlap = 20)
    {
        var arr = global::Seqeron.Genomics.Alignment.SequenceAssembler.CalculateCoverage(
            reference ?? string.Empty, reads ?? Array.Empty<string>(), minOverlap);
        return new CoverageResult(arr);
    }

    [McpServerTool(Name = "compute_consensus", Title = "Assembly — Column-Majority Consensus", ReadOnly = true)]
    [Description("Builds a consensus sequence from a list of pre-aligned reads (same length, '-' / 'N' ignored) by majority vote per column. Columns with no informative bases yield 'N'.")]
    public static ConsensusResult ComputeConsensus(
        [Description("Pre-aligned reads of identical length.")] string[] alignedReads)
    {
        if (alignedReads != null && alignedReads.Length > 1)
        {
            int len = alignedReads[0]?.Length ?? 0;
            for (int i = 1; i < alignedReads.Length; i++)
            {
                if ((alignedReads[i]?.Length ?? 0) != len)
                    throw new ArgumentException(
                        $"All aligned reads must have the same length; alignedReads[{i}] differs.",
                        nameof(alignedReads));
            }
        }

        var s = global::Seqeron.Genomics.Alignment.SequenceAssembler.ComputeConsensus(
            alignedReads ?? Array.Empty<string>());
        return new ConsensusResult(s);
    }

    [McpServerTool(Name = "quality_trim_reads", Title = "Reads — Quality Trim (Phred+33)", ReadOnly = true)]
    [Description("Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. Reads shorter than minLength after trimming are filtered out.")]
    public static TrimmedReadsResult QualityTrimReads(
        [Description("Reads, each with sequence and matching Phred+33 quality string.")] QualityReadInput[] reads,
        [Description("Minimum acceptable quality score (Phred).")] int minQuality = 20,
        [Description("Minimum length after trimming for a read to be kept.")] int minLength = 50)
    {
        var input = new List<(string sequence, string quality)>(reads?.Length ?? 0);
        if (reads != null)
        {
            for (int i = 0; i < reads.Length; i++)
            {
                var r = reads[i];
                if (r is null)
                    throw new ArgumentException($"reads[{i}] is null.", nameof(reads));
                var seq = r.Sequence ?? string.Empty;
                var qual = r.Quality ?? string.Empty;
                if (seq.Length != qual.Length)
                    throw new ArgumentException(
                        $"reads[{i}]: sequence and quality must be the same length.", nameof(reads));
                input.Add((seq, qual));
            }
        }

        var trimmed = global::Seqeron.Genomics.Alignment.SequenceAssembler.QualityTrimReads(
            input, minQuality, minLength);

        var arr = new string[trimmed.Count];
        for (int i = 0; i < trimmed.Count; i++)
            arr[i] = trimmed[i];
        return new TrimmedReadsResult(arr);
    }

    [McpServerTool(Name = "error_correct_reads", Title = "Reads — k-mer Frequency Error Correction", ReadOnly = true)]
    [Description("Corrects single-base errors in reads using k-mer frequency: any k-mer occurring fewer than minKmerFrequency times is corrected by substituting its middle base if a higher-frequency k-mer can be obtained.")]
    public static CorrectedReadsResult ErrorCorrectReads(
        [Description("Input reads.")] string[] reads,
        [Description("K-mer size used to build the frequency table.")] int kmerSize = 21,
        [Description("Minimum k-mer frequency considered \"trusted\".")] int minKmerFrequency = 3)
    {
        var corrected = global::Seqeron.Genomics.Alignment.SequenceAssembler.ErrorCorrectReads(
            reads ?? Array.Empty<string>(), kmerSize, minKmerFrequency);

        var arr = new string[corrected.Count];
        for (int i = 0; i < corrected.Count; i++)
            arr[i] = corrected[i];
        return new CorrectedReadsResult(arr);
    }

    #endregion

    #region Mappers

    private static AlignmentResultDto ToAlignmentDto(global::Seqeron.Genomics.Infrastructure.AlignmentResult r) =>
        new(
            r.AlignedSequence1,
            r.AlignedSequence2,
            r.Score,
            r.AlignmentType.ToString(),
            r.StartPosition1,
            r.StartPosition2,
            r.EndPosition1,
            r.EndPosition2);

    private static AssemblyResultDto ToAssemblyDto(global::Seqeron.Genomics.Alignment.SequenceAssembler.AssemblyResult r)
    {
        var arr = new string[r.Contigs.Count];
        for (int i = 0; i < r.Contigs.Count; i++)
            arr[i] = r.Contigs[i];
        return new AssemblyResultDto(arr, r.TotalReads, r.AssembledReads, r.N50, r.LongestContig, r.TotalLength);
    }

    private static ApproximateMatchDto ToApproximateMatchDto(
        global::Seqeron.Genomics.Alignment.ApproximateMatchResult m)
    {
        var positions = new int[m.MismatchPositions.Count];
        for (int i = 0; i < m.MismatchPositions.Count; i++)
            positions[i] = m.MismatchPositions[i];
        return new ApproximateMatchDto(
            m.Position,
            m.MatchedSequence,
            m.Distance,
            positions,
            m.MismatchType.ToString());
    }

    #endregion
}
