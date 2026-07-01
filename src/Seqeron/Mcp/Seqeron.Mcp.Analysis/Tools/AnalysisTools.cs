using System.ComponentModel;
using System.Linq;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Mcp.Analysis.Tools;

/// <summary>
/// MCP tools wrapping <c>Seqeron.Genomics.Analysis</c>.
/// </summary>
[McpServerToolType]
public class AnalysisTools
{
    private static DnaSequence RequireDna(string sequence, string paramName)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", paramName);
        if (!DnaSequence.TryCreate(sequence, out var dna) || dna is null)
            throw new ArgumentException("Invalid DNA sequence", paramName);
        return dna;
    }

    #region KmerAnalyzer

    [McpServerTool(Name = "count_kmers", Title = "k-mers — Count All", ReadOnly = true)]
    [Description("Counts every k-mer (substring of length k) occurrence in a sequence. Returns k-mer → count.")]
    public static KmerCountsResult CountKmers(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length (>0).")] int k)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (k <= 0)
            throw new ArgumentException("k must be positive", nameof(k));

        var counts = KmerAnalyzer.CountKmers(sequence, k);
        return new KmerCountsResult(counts);
    }

    [McpServerTool(Name = "kmer_spectrum", Title = "k-mers — Spectrum", ReadOnly = true)]
    [Description("Frequency-of-frequencies: for each occurrence count, how many distinct k-mers reach that count.")]
    public static KmerSpectrumResult KmerSpectrum(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length.")] int k)
    {
        var spectrum = KmerAnalyzer.GetKmerSpectrum(sequence ?? string.Empty, k);
        return new KmerSpectrumResult(spectrum);
    }

    [McpServerTool(Name = "most_frequent_kmers", Title = "k-mers — Most Frequent", ReadOnly = true)]
    [Description("Returns all k-mers tied for the maximum occurrence count.")]
    public static KmerListResult MostFrequentKmers(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length.")] int k)
    {
        var kmers = KmerAnalyzer.FindMostFrequentKmers(sequence ?? string.Empty, k).ToArray();
        return new KmerListResult(kmers);
    }

    [McpServerTool(Name = "kmer_frequencies", Title = "k-mers — Normalized Frequencies", ReadOnly = true)]
    [Description("Normalized k-mer counts (each value in [0,1], summing to 1).")]
    public static KmerFrequenciesResult KmerFrequencies(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length.")] int k)
    {
        var freq = KmerAnalyzer.GetKmerFrequencies(sequence ?? string.Empty, k);
        return new KmerFrequenciesResult(freq);
    }

    [McpServerTool(Name = "kmer_distance", Title = "k-mers — Euclidean Distance", ReadOnly = true)]
    [Description("Euclidean distance between k-mer frequency vectors of two sequences. 0 means identical k-mer composition.")]
    public static KmerDistanceResult KmerDistance(
        [Description("First sequence.")] string seq1,
        [Description("Second sequence.")] string seq2,
        [Description("k-mer length.")] int k)
    {
        var d = KmerAnalyzer.KmerDistance(seq1 ?? string.Empty, seq2 ?? string.Empty, k);
        return new KmerDistanceResult(d);
    }

    [McpServerTool(Name = "unique_kmers", Title = "k-mers — Unique (Singletons)", ReadOnly = true)]
    [Description("k-mers that occur exactly once in the sequence.")]
    public static KmerListResult UniqueKmers(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length.")] int k)
    {
        var kmers = KmerAnalyzer.FindUniqueKmers(sequence ?? string.Empty, k).ToArray();
        return new KmerListResult(kmers);
    }

    [McpServerTool(Name = "kmers_with_min_count", Title = "k-mers — Min-Count Filter", ReadOnly = true)]
    [Description("k-mers occurring at least minCount times, sorted descending by count.")]
    public static KmersWithMinCountResult KmersWithMinCount(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length.")] int k,
        [Description("Minimum occurrence count.")] int minCount)
    {
        var items = KmerAnalyzer.FindKmersWithMinCount(sequence ?? string.Empty, k, minCount)
            .Select(t => new KmerCountItem(t.Kmer, t.Count))
            .ToArray();
        return new KmersWithMinCountResult(items);
    }

    [McpServerTool(Name = "generate_all_kmers", Title = "k-mers — Enumerate Alphabet Space", ReadOnly = true)]
    [Description("Enumerate the entire k-mer space for an alphabet (default \"ACGT\"). Result size = alphabet.Length^k.")]
    public static KmerListResult GenerateAllKmers(
        [Description("k-mer length (>0).")] int k,
        [Description("Alphabet (default \"ACGT\").")] string alphabet = "ACGT")
    {
        var kmers = KmerAnalyzer.GenerateAllKmers(k, alphabet).ToArray();
        return new KmerListResult(kmers);
    }

    [McpServerTool(Name = "find_clumps", Title = "k-mers — Find Clumps", ReadOnly = true)]
    [Description("Finds k-mers that occur at least minOccurrences times within any sliding window of size windowSize. Bioinformatics Algorithms Ch. 1.")]
    public static KmerListResult FindClumps(
        [Description("Sequence to scan.")] string sequence,
        [Description("k-mer length.")] int k,
        [Description("Sliding window size (>= k).")] int windowSize,
        [Description("Minimum occurrences within a window (>0).")] int minOccurrences)
    {
        var kmers = KmerAnalyzer.FindClumps(sequence ?? string.Empty, k, windowSize, minOccurrences).ToArray();
        return new KmerListResult(kmers);
    }

    [McpServerTool(Name = "kmer_positions", Title = "k-mers — Find Positions", ReadOnly = true)]
    [Description("Zero-based positions of all (overlapping) occurrences of a k-mer.")]
    public static KmerPositionsResult KmerPositions(
        [Description("Sequence to scan.")] string sequence,
        [Description("k-mer to locate.")] string kmer)
    {
        var positions = KmerAnalyzer.FindKmerPositions(sequence ?? string.Empty, kmer ?? string.Empty).ToArray();
        return new KmerPositionsResult(positions);
    }

    [McpServerTool(Name = "count_kmers_both_strands", Title = "k-mers — Count Both Strands", ReadOnly = true)]
    [Description("k-mer counts on the forward strand combined with counts on the reverse-complement strand.")]
    public static KmerCountsResult CountKmersBothStrands(
        [Description("DNA sequence.")] string sequence,
        [Description("k-mer length.")] int k)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var counts = KmerAnalyzer.CountKmersBothStrands(dna, k);
        return new KmerCountsResult(counts);
    }

    [McpServerTool(Name = "analyze_kmers", Title = "k-mers — Aggregate Statistics", ReadOnly = true)]
    [Description("Aggregate k-mer statistics: total, unique, min/max/avg count, and Shannon entropy.")]
    public static AnalyzeKmersResult AnalyzeKmers(
        [Description("Sequence to analyze.")] string sequence,
        [Description("k-mer length (>0).")] int k)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (k <= 0)
            throw new ArgumentException("k must be positive", nameof(k));

        var s = KmerAnalyzer.AnalyzeKmers(sequence, k);
        return new AnalyzeKmersResult(
            s.TotalKmers, s.UniqueKmers, s.MaxCount, s.MinCount, s.AverageCount, s.Entropy);
    }

    #endregion

    #region SequenceStatistics

    [McpServerTool(Name = "hydrophobicity_profile", Title = "Protein — Hydrophobicity Profile", ReadOnly = true)]
    [Description("Sliding-window Kyte-Doolittle hydropathy values for a protein sequence.")]
    public static DoubleProfileResult HydrophobicityProfile(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Sliding window size (default 9).")] int windowSize = 9)
    {
        var values = SequenceStatistics
            .CalculateHydrophobicityProfile(proteinSequence ?? string.Empty, windowSize)
            .ToArray();
        return new DoubleProfileResult(values);
    }

    [McpServerTool(Name = "dinucleotide_frequencies", Title = "Sequence — Dinucleotide Frequencies", ReadOnly = true)]
    [Description("Frequency of each adjacent dinucleotide over the alphabet A/T/G/C/U.")]
    public static DinucleotideFrequenciesResult DinucleotideFrequencies(
        [Description("Nucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(sequence);
        return new DinucleotideFrequenciesResult(new Dictionary<string, double>(freq));
    }

    [McpServerTool(Name = "dinucleotide_ratios", Title = "Sequence — Dinucleotide Obs/Exp Ratios", ReadOnly = true)]
    [Description("Observed/expected ratios for each dinucleotide (e.g., CpG ratio for CpG-island detection).")]
    public static DinucleotideRatiosResult DinucleotideRatios(
        [Description("Nucleotide sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var ratios = SequenceStatistics.CalculateDinucleotideRatios(sequence);
        return new DinucleotideRatiosResult(new Dictionary<string, double>(ratios));
    }

    [McpServerTool(Name = "codon_frequencies", Title = "DNA — Codon Frequencies", ReadOnly = true)]
    [Description("Codon usage frequencies in the specified reading frame (0..2).")]
    public static CodonFrequenciesResult CodonFrequencies(
        [Description("DNA sequence.")] string dnaSequence,
        [Description("Reading frame: 0, 1, or 2 (default 0).")] int readingFrame = 0)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(dnaSequence));
        if (readingFrame is < 0 or > 2)
            throw new ArgumentException("Reading frame must be 0, 1, or 2", nameof(readingFrame));

        var freq = SequenceStatistics.CalculateCodonFrequencies(dnaSequence, readingFrame);
        return new CodonFrequenciesResult(new Dictionary<string, double>(freq));
    }

    [McpServerTool(Name = "predict_chou_fasman", Title = "Protein — Chou-Fasman Propensities", ReadOnly = true)]
    [Description("Per-window helix/sheet/turn propensities for a protein sequence (Chou-Fasman parameters).")]
    public static PredictChouFasmanResult PredictChouFasman(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Sliding window size (default 7).")] int windowSize = 7)
    {
        var items = SequenceStatistics
            .PredictSecondaryStructure(proteinSequence ?? string.Empty, windowSize)
            .Select(t => new ChouFasmanItem(t.Helix, t.Sheet, t.Turn))
            .ToArray();
        return new PredictChouFasmanResult(items);
    }

    [McpServerTool(Name = "gc_content_profile", Title = "Sequence — GC Content Profile", ReadOnly = true)]
    [Description("GC content in sliding windows along the sequence.")]
    public static DoubleProfileResult GcContentProfile(
        [Description("Nucleotide sequence.")] string sequence,
        [Description("Window size (default 100).")] int windowSize = 100,
        [Description("Step size (default 1).")] int stepSize = 1)
    {
        var values = SequenceStatistics
            .CalculateGcContentProfile(sequence ?? string.Empty, windowSize, stepSize)
            .ToArray();
        return new DoubleProfileResult(values);
    }

    [McpServerTool(Name = "entropy_profile", Title = "Sequence — Shannon Entropy Profile", ReadOnly = true)]
    [Description("Shannon entropy in sliding windows along the sequence.")]
    public static DoubleProfileResult EntropyProfile(
        [Description("Sequence to analyze.")] string sequence,
        [Description("Window size (default 50).")] int windowSize = 50,
        [Description("Step size (default 1).")] int stepSize = 1)
    {
        var values = SequenceStatistics
            .CalculateEntropyProfile(sequence ?? string.Empty, windowSize, stepSize)
            .ToArray();
        return new DoubleProfileResult(values);
    }

    #endregion

    #region GenomicAnalyzer

    [McpServerTool(Name = "find_repeats", Title = "Genome — Repeated Substrings", ReadOnly = true)]
    [Description("All repeated substrings of length >= minLength in a DNA sequence, with their positions.")]
    public static FindRepeatsResult FindRepeats(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum repeat length.")] int minLength)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = GenomicAnalyzer.FindRepeats(dna, minLength)
            .Select(r => new RepeatItem(r.Sequence, r.Positions.ToArray(), r.Length, r.Count))
            .ToArray();
        return new FindRepeatsResult(items);
    }

    [McpServerTool(Name = "find_tandem_repeats", Title = "Genome — Tandem Repeats", ReadOnly = true)]
    [Description("Consecutive repeating units (e.g., ATGATGATG) of unit-length >= minUnitLength repeated >= minRepetitions times.")]
    public static FindTandemRepeatsResult FindTandemRepeats(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum repeat-unit length (default 2).")] int minUnitLength = 2,
        [Description("Minimum number of repetitions (default 2).")] int minRepetitions = 2)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength, minRepetitions)
            .Select(t => new TandemRepeatItem(t.Unit, t.Position, t.Repetitions, t.TotalLength))
            .ToArray();
        return new FindTandemRepeatsResult(items);
    }

    [McpServerTool(Name = "find_motif", Title = "Genome — Find Motif (Suffix Tree)", ReadOnly = true)]
    [Description("All exact occurrences of a motif (case-insensitive) in a DNA sequence via suffix tree.")]
    public static MotifPositionsResult FindMotif(
        [Description("DNA sequence.")] string sequence,
        [Description("Motif to locate.")] string motif)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var positions = GenomicAnalyzer.FindMotif(dna, motif ?? string.Empty);
        return new MotifPositionsResult(positions.ToArray());
    }

    [McpServerTool(Name = "find_known_motifs", Title = "Genome — Search Known Motif Set", ReadOnly = true)]
    [Description("Search for a user-provided set of motifs simultaneously. Only motifs with at least one hit are returned.")]
    public static FindKnownMotifsResult FindKnownMotifs(
        [Description("DNA sequence.")] string sequence,
        [Description("Set of motifs to search.")] string[] motifs)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var raw = GenomicAnalyzer.FindKnownMotifs(dna, motifs ?? Array.Empty<string>());
        var matches = new Dictionary<string, int[]>(raw.Count);
        foreach (var kvp in raw)
            matches[kvp.Key] = kvp.Value.ToArray();
        return new FindKnownMotifsResult(matches);
    }

    [McpServerTool(Name = "find_common_regions", Title = "Genome — Common Regions", ReadOnly = true)]
    [Description("All common regions between two DNA sequences with length >= minLength.")]
    public static FindCommonRegionsResult FindCommonRegions(
        [Description("First DNA sequence.")] string sequence1,
        [Description("Second DNA sequence.")] string sequence2,
        [Description("Minimum common-region length.")] int minLength)
    {
        var dna1 = RequireDna(sequence1, nameof(sequence1));
        var dna2 = RequireDna(sequence2, nameof(sequence2));
        var items = GenomicAnalyzer.FindCommonRegions(dna1, dna2, minLength)
            .Select(r => new CommonRegionItem(r.Sequence, r.PositionInFirst, r.PositionInSecond, r.Length))
            .ToArray();
        return new FindCommonRegionsResult(items);
    }

    [McpServerTool(Name = "find_open_reading_frames", Title = "Genome — Open Reading Frames", ReadOnly = true)]
    [Description("ORFs in all 6 frames (3 forward + 3 reverse-complement) starting ATG and ending TAA/TAG/TGA.")]
    public static FindOpenReadingFramesResult FindOpenReadingFrames(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum ORF length in nucleotides (default 100).")] int minLength = 100)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = GenomicAnalyzer.FindOpenReadingFrames(dna, minLength)
            .Select(o => new OrfItem(
                o.Sequence, o.Position, o.Frame, o.IsReverseComplement, o.Length, o.CodonCount))
            .ToArray();
        return new FindOpenReadingFramesResult(items);
    }

    #endregion

    #region RepeatFinder

    [McpServerTool(Name = "find_microsatellites", Title = "Repeats — Microsatellites (STR)", ReadOnly = true)]
    [Description("Short Tandem Repeats (STRs): 1-6 bp motif units repeated consecutively.")]
    public static FindMicrosatellitesResult FindMicrosatellites(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum unit length (default 1).")] int minUnitLength = 1,
        [Description("Maximum unit length (default 6).")] int maxUnitLength = 6,
        [Description("Minimum number of repeats (default 3).")] int minRepeats = 3)
    {
        var items = global::Seqeron.Genomics.Analysis.RepeatFinder
            .FindMicrosatellites(sequence ?? string.Empty, minUnitLength, maxUnitLength, minRepeats)
            .Select(m => new MicrosatelliteItem(
                m.Position, m.RepeatUnit, m.RepeatCount, m.TotalLength, m.RepeatType.ToString()))
            .ToArray();
        return new FindMicrosatellitesResult(items);
    }

    [McpServerTool(Name = "find_inverted_repeats", Title = "Repeats — Inverted Repeats / Hairpins", ReadOnly = true)]
    [Description("Sequences whose two arms are reverse-complement of each other (hairpin candidates).")]
    public static FindInvertedRepeatsResult FindInvertedRepeats(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum arm length (default 4).")] int minArmLength = 4,
        [Description("Maximum loop length (default 50).")] int maxLoopLength = 50,
        [Description("Minimum loop length (default 3).")] int minLoopLength = 3)
    {
        var items = global::Seqeron.Genomics.Analysis.RepeatFinder
            .FindInvertedRepeats(sequence ?? string.Empty, minArmLength, maxLoopLength, minLoopLength)
            .Select(r => new InvertedRepeatItem(
                r.LeftArmStart, r.RightArmStart, r.ArmLength, r.LoopLength,
                r.LeftArm, r.RightArm, r.Loop, r.CanFormHairpin, r.TotalLength))
            .ToArray();
        return new FindInvertedRepeatsResult(items);
    }

    [McpServerTool(Name = "find_direct_repeats", Title = "Repeats — Direct Repeats", ReadOnly = true)]
    [Description("Identical sequences appearing twice with a spacer between them.")]
    public static FindDirectRepeatsResult FindDirectRepeats(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum repeat length (default 5).")] int minLength = 5,
        [Description("Maximum repeat length (default 50).")] int maxLength = 50,
        [Description("Minimum spacing between the two copies (default 1).")] int minSpacing = 1)
    {
        var items = global::Seqeron.Genomics.Analysis.RepeatFinder
            .FindDirectRepeats(sequence ?? string.Empty, minLength, maxLength, minSpacing)
            .Select(r => new DirectRepeatItem(
                r.FirstPosition, r.SecondPosition, r.RepeatSequence, r.Length, r.Spacing))
            .ToArray();
        return new FindDirectRepeatsResult(items);
    }

    [McpServerTool(Name = "tandem_repeat_summary", Title = "Repeats — Tandem Repeat Summary", ReadOnly = true)]
    [Description("Aggregate statistics across all microsatellites in a DNA sequence.")]
    public static TandemRepeatSummaryResult TandemRepeatSummary(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum number of repeats (default 3).")] int minRepeats = 3)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var s = global::Seqeron.Genomics.Analysis.RepeatFinder.GetTandemRepeatSummary(dna, minRepeats);
        MicrosatelliteItem? longest = s.LongestRepeat is { } lr
            ? new MicrosatelliteItem(lr.Position, lr.RepeatUnit, lr.RepeatCount, lr.TotalLength, lr.RepeatType.ToString())
            : null;
        return new TandemRepeatSummaryResult(
            s.TotalRepeats,
            s.TotalRepeatBases,
            s.PercentageOfSequence,
            s.MononucleotideRepeats,
            s.DinucleotideRepeats,
            s.TrinucleotideRepeats,
            s.TetranucleotideRepeats,
            longest,
            s.MostFrequentUnit);
    }

    [McpServerTool(Name = "find_palindromes", Title = "Repeats — Palindromes (Restriction Sites)", ReadOnly = true)]
    [Description("Sequences identical to their reverse complement (restriction-site candidates). minLength must be even and >= 4.")]
    public static FindPalindromesResult FindPalindromes(
        [Description("DNA sequence.")] string sequence,
        [Description("Minimum palindrome length, even, >= 4 (default 4).")] int minLength = 4,
        [Description("Maximum palindrome length (default 12).")] int maxLength = 12)
    {
        var items = global::Seqeron.Genomics.Analysis.RepeatFinder
            .FindPalindromes(sequence ?? string.Empty, minLength, maxLength)
            .Select(p => new PalindromeItem(p.Position, p.Sequence, p.Length))
            .ToArray();
        return new FindPalindromesResult(items);
    }

    #endregion

    #region MotifFinder

    [McpServerTool(Name = "find_exact_motif", Title = "Motifs — Exact Match (DNA)", ReadOnly = true)]
    [Description("Exact-match motif positions in a DNA sequence via suffix tree.")]
    public static MotifPositionsResult FindExactMotif(
        [Description("DNA sequence to search.")] string sequence,
        [Description("Motif pattern.")] string motif)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var positions = global::Seqeron.Genomics.Analysis.MotifFinder
            .FindExactMotif(dna, motif ?? string.Empty)
            .ToArray();
        return new MotifPositionsResult(positions);
    }

    [McpServerTool(Name = "find_degenerate_motif", Title = "Motifs — Degenerate (IUPAC)", ReadOnly = true)]
    [Description("Motif search with IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V).")]
    public static FindDegenerateMotifResult FindDegenerateMotif(
        [Description("DNA sequence.")] string sequence,
        [Description("IUPAC motif pattern.")] string motif)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = global::Seqeron.Genomics.Analysis.MotifFinder
            .FindDegenerateMotif(dna, motif ?? string.Empty)
            .Select(m => new MotifMatchItem(m.Position, m.MatchedSequence, m.Pattern, m.Score))
            .ToArray();
        return new FindDegenerateMotifResult(items);
    }

    [McpServerTool(Name = "create_pwm", Title = "Motifs — Build PWM", ReadOnly = true)]
    [Description("Build a log-odds Position Weight Matrix (4×L; rows A,C,G,T) from aligned, equal-length DNA sequences.")]
    public static PwmResult CreatePwm(
        [Description("Aligned DNA sequences of equal length.")] string[] sequences,
        [Description("Pseudocount for smoothing (default 0.25).")] double pseudocount = 0.25)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(sequences));
        if (pseudocount < 0)
            throw new ArgumentOutOfRangeException(nameof(pseudocount), "Pseudocount must be non-negative.");

        var pwm = global::Seqeron.Genomics.Analysis.MotifFinder
            .CreatePwm(sequences, pseudocount);
        var jagged = MatrixToJagged(pwm.Matrix, 4, pwm.Length);
        return new PwmResult(jagged, pwm.Length, pwm.Consensus, pwm.MaxScore, pwm.MinScore);
    }

    [McpServerTool(Name = "scan_with_pwm", Title = "Motifs — Scan with PWM", ReadOnly = true)]
    [Description("Scan a DNA sequence with a 4×L Position Weight Matrix; rows = A,C,G,T. Returns matches scoring at or above threshold.")]
    public static ScanWithPwmResult ScanWithPwm(
        [Description("DNA sequence to scan.")] string sequence,
        [Description("Position Weight Matrix: matrix is jagged 4×L (rows A,C,G,T), length is L.")] PwmInput pwm,
        [Description("Minimum score threshold (default 0.0).")] double threshold = 0.0)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        if (pwm is null) throw new ArgumentNullException(nameof(pwm));
        var matrix2d = JaggedToMatrix(pwm.Matrix, 4, pwm.Length);
        var pwmObj = new global::Seqeron.Genomics.Analysis.PositionWeightMatrix(matrix2d, pwm.Length);
        var items = global::Seqeron.Genomics.Analysis.MotifFinder
            .ScanWithPwm(dna, pwmObj, threshold)
            .Select(m => new MotifMatchItem(m.Position, m.MatchedSequence, m.Pattern, m.Score))
            .ToArray();
        return new ScanWithPwmResult(items);
    }

    [McpServerTool(Name = "generate_consensus", Title = "Motifs — IUPAC Consensus", ReadOnly = true)]
    [Description("IUPAC consensus sequence from aligned equal-length DNA sequences (>25% per position threshold).")]
    public static ConsensusResult GenerateConsensus(
        [Description("Aligned DNA sequences of equal length.")] string[] sequences)
    {
        var consensus = global::Seqeron.Genomics.Analysis.MotifFinder
            .GenerateConsensus(sequences ?? Array.Empty<string>());
        return new ConsensusResult(consensus);
    }

    [McpServerTool(Name = "discover_motifs", Title = "Motifs — De Novo Discovery", ReadOnly = true)]
    [Description("Overrepresented k-mers (de novo motif discovery) in a DNA sequence.")]
    public static DiscoverMotifsResult DiscoverMotifs(
        [Description("DNA sequence.")] string sequence,
        [Description("k-mer length (default 6).")] int k = 6,
        [Description("Minimum occurrence count (default 2).")] int minCount = 2)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = global::Seqeron.Genomics.Analysis.MotifFinder
            .DiscoverMotifs(dna, k, minCount)
            .Select(m => new DiscoveredMotifItem(m.Sequence, m.Count, m.Positions.ToArray(), m.Enrichment))
            .ToArray();
        return new DiscoverMotifsResult(items);
    }

    [McpServerTool(Name = "find_shared_motifs", Title = "Motifs — Shared Across Sequences", ReadOnly = true)]
    [Description("k-mers present in at least minSequences of the input DNA sequences.")]
    public static FindSharedMotifsResult FindSharedMotifs(
        [Description("DNA sequences.")] string[] sequences,
        [Description("k-mer length (default 6).")] int k = 6,
        [Description("Minimum sequences containing the motif (default 2).")] int minSequences = 2)
    {
        var dnaList = (sequences ?? Array.Empty<string>())
            .Select(s => RequireDna(s, nameof(sequences)))
            .ToList();
        var items = global::Seqeron.Genomics.Analysis.MotifFinder
            .FindSharedMotifs(dnaList, k, minSequences)
            .Select(m => new SharedMotifItem(m.Sequence, m.SequenceIndices.ToArray(), m.Prevalence))
            .ToArray();
        return new FindSharedMotifsResult(items);
    }

    [McpServerTool(Name = "find_regulatory_elements", Title = "Motifs — Regulatory Elements", ReadOnly = true)]
    [Description("Scan for built-in regulatory motifs (TATA, CAAT, GC-box, Kozak, Shine-Dalgarno, poly(A), E-box, AP-1, NF-κB, CREB).")]
    public static FindRegulatoryElementsResult FindRegulatoryElements(
        [Description("DNA sequence.")] string sequence)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = global::Seqeron.Genomics.Analysis.MotifFinder
            .FindRegulatoryElements(dna)
            .Select(r => new RegulatoryElementItem(r.Name, r.Position, r.Sequence, r.Pattern, r.Description))
            .ToArray();
        return new FindRegulatoryElementsResult(items);
    }

    private static double[][] MatrixToJagged(double[,] matrix, int rows, int cols)
    {
        var result = new double[rows][];
        for (int r = 0; r < rows; r++)
        {
            var row = new double[cols];
            for (int c = 0; c < cols; c++) row[c] = matrix[r, c];
            result[r] = row;
        }
        return result;
    }

    private static double[,] JaggedToMatrix(double[][] jagged, int rows, int cols)
    {
        if (jagged is null || jagged.Length != rows)
            throw new ArgumentException($"Matrix must have {rows} rows.", nameof(jagged));
        var result = new double[rows, cols];
        for (int r = 0; r < rows; r++)
        {
            var row = jagged[r];
            if (row is null || row.Length != cols)
                throw new ArgumentException($"Row {r} must have length {cols}.", nameof(jagged));
            for (int c = 0; c < cols; c++) result[r, c] = row[c];
        }
        return result;
    }

    #endregion

    #region ProteinMotifFinder

    [McpServerTool(Name = "find_protein_motifs", Title = "Protein — Common Motifs (PROSITE)", ReadOnly = true)]
    [Description("Scan a protein sequence against the built-in PROSITE-style motif catalog (N-glycosylation, kinase phosphorylation sites, ATP/GTP P-loop, EF-hand, zinc finger, leucine zipper, NLS, NES, SIM, WW, SH3, etc.).")]
    public static FindProteinMotifsResult FindProteinMotifs(
        [Description("Protein sequence.")] string proteinSequence)
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .FindCommonMotifs(proteinSequence ?? string.Empty)
            .Select(m => new ProteinMotifMatchItem(m.Start, m.End, m.Sequence, m.MotifName, m.Pattern, m.Score, m.EValue))
            .ToArray();
        return new FindProteinMotifsResult(items);
    }

    [McpServerTool(Name = "find_motif_by_pattern", Title = "Protein — Find Motif by Regex", ReadOnly = true)]
    [Description("Find all (overlapping) regex pattern matches in a protein sequence.")]
    public static FindProteinMotifsResult FindMotifByPattern(
        [Description("Protein sequence.")] string proteinSequence,
        [Description(".NET regex pattern.")] string regexPattern,
        [Description("Motif display name (default \"Custom\").")] string motifName = "Custom",
        [Description("Pattern identifier (default empty).")] string patternId = "")
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .FindMotifByPattern(proteinSequence ?? string.Empty, regexPattern ?? string.Empty, motifName, patternId)
            .Select(m => new ProteinMotifMatchItem(m.Start, m.End, m.Sequence, m.MotifName, m.Pattern, m.Score, m.EValue))
            .ToArray();
        return new FindProteinMotifsResult(items);
    }

    [McpServerTool(Name = "find_motif_by_prosite", Title = "Protein — Find Motif by PROSITE Pattern", ReadOnly = true)]
    [Description("Convert a PROSITE pattern to regex internally, then scan a protein sequence.")]
    public static FindProteinMotifsResult FindMotifByProsite(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("PROSITE-format pattern.")] string prositePattern,
        [Description("Motif display name (default \"Custom\").")] string motifName = "Custom")
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .FindMotifByProsite(proteinSequence ?? string.Empty, prositePattern ?? string.Empty, motifName)
            .Select(m => new ProteinMotifMatchItem(m.Start, m.End, m.Sequence, m.MotifName, m.Pattern, m.Score, m.EValue))
            .ToArray();
        return new FindProteinMotifsResult(items);
    }

    [McpServerTool(Name = "prosite_to_regex", Title = "Protein — PROSITE → Regex", ReadOnly = true)]
    [Description("Translate a PROSITE pattern string to a .NET regex string.")]
    public static PrositeRegexResult PrositeToRegex(
        [Description("PROSITE-format pattern.")] string prositePattern)
    {
        var regex = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .ConvertPrositeToRegex(prositePattern ?? string.Empty);
        return new PrositeRegexResult(regex);
    }

    [McpServerTool(Name = "predict_signal_peptide", Title = "Protein — Signal Peptide (von Heijne)", ReadOnly = true)]
    [Description("von Heijne (1986) weight-matrix signal-peptide cleavage-site prediction (EMBOSS sigcleave).")]
    public static SignalPeptideResult PredictSignalPeptide(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Score against the prokaryotic matrix instead of the default eukaryotic one.")] bool prokaryote = false)
    {
        var sp = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .PredictSignalPeptide(proteinSequence ?? string.Empty, prokaryote);
        if (sp is null)
            return new SignalPeptideResult(false, 0, 0.0, "", "", false);
        var v = sp.Value;
        return new SignalPeptideResult(true, v.CleavagePosition, v.Score, v.SignalSequence, v.WindowSequence, v.IsLikelySignalPeptide);
    }

    [McpServerTool(Name = "predict_transmembrane_helices", Title = "Protein — Transmembrane Helices (Kyte-Doolittle)", ReadOnly = true)]
    [Description("Hydropathy-based transmembrane helix prediction (Kyte-Doolittle, ≥15 aa).")]
    public static PredictTransmembraneHelicesResult PredictTransmembraneHelices(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Sliding window size (default 19).")] int windowSize = 19,
        [Description("Hydropathy threshold (default 1.6).")] double threshold = 1.6)
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .PredictTransmembraneHelices(proteinSequence ?? string.Empty, windowSize, threshold)
            .Select(t => new RegionScoreItem(t.Start, t.End, t.Score))
            .ToArray();
        return new PredictTransmembraneHelicesResult(items);
    }

    [McpServerTool(Name = "predict_coiled_coils", Title = "Protein — Coiled-Coils", ReadOnly = true)]
    [Description("Heptad-repeat-based coiled-coil prediction.")]
    public static PredictCoiledCoilsResult PredictCoiledCoils(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Sliding window size (default 28).")] int windowSize = 28,
        [Description("Score threshold (default 0.5).")] double threshold = 0.5)
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .PredictCoiledCoils(proteinSequence ?? string.Empty, windowSize, threshold)
            .Select(t => new RegionScoreItem(t.Start, t.End, t.Score))
            .ToArray();
        return new PredictCoiledCoilsResult(items);
    }

    [McpServerTool(Name = "find_protein_low_complexity_regions", Title = "Protein — Low-Complexity Regions", ReadOnly = true)]
    [Description("Low-complexity regions in a protein via the SEG algorithm (Wootton & Federhen 1993): sliding-window Shannon entropy in bits/residue, two-pass trigger/extension.")]
    public static FindProteinLowComplexityRegionsResult FindProteinLowComplexityRegions(
        [Description("Protein sequence.")] string proteinSequence,
        [Description("Sliding window size W (default 12).")] int windowSize = 12,
        [Description("Trigger complexity K1 in bits/residue (default 2.2).")] double triggerComplexity = 2.2,
        [Description("Extension complexity K2 in bits/residue (default 2.5).")] double extensionComplexity = 2.5)
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .FindLowComplexityRegions(proteinSequence ?? string.Empty, windowSize, triggerComplexity, extensionComplexity)
            .Select(t => new ProteinLowComplexityItem(t.Start, t.End, t.Complexity))
            .ToArray();
        return new FindProteinLowComplexityRegionsResult(items);
    }

    [McpServerTool(Name = "find_protein_domains", Title = "Protein — Common Domains", ReadOnly = true)]
    [Description("Detect common protein domains with EXACT PROSITE patterns: zinc finger C2H2 (PS00028), WD-repeats (PS00678), kinase ATP-binding / Walker A P-loop (PS00017). Profile-only domains (SH3 PS50002, PDZ PS50106) are not detected (no deterministic pattern).")]
    public static FindProteinDomainsResult FindProteinDomains(
        [Description("Protein sequence.")] string proteinSequence)
    {
        var items = global::Seqeron.Genomics.Analysis.ProteinMotifFinder
            .FindDomains(proteinSequence ?? string.Empty)
            .Select(d => new ProteinDomainItem(d.Name, d.Accession, d.Start, d.End, d.Score, d.Description))
            .ToArray();
        return new FindProteinDomainsResult(items);
    }

    #endregion

    #region SequenceComplexity

    [McpServerTool(Name = "windowed_complexity", Title = "Complexity — Sliding Window (DNA)", ReadOnly = true)]
    [Description("Sliding-window Shannon entropy + linguistic complexity for a DNA sequence.")]
    public static WindowedComplexityResult WindowedComplexity(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size (default 64).")] int windowSize = 64,
        [Description("Step size (default 10).")] int stepSize = 10)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = global::Seqeron.Genomics.Analysis.SequenceComplexity
            .CalculateWindowedComplexity(dna, windowSize, stepSize)
            .Select(p => new ComplexityPointItem(p.Position, p.ShannonEntropy, p.LinguisticComplexity, p.WindowStart, p.WindowEnd))
            .ToArray();
        return new WindowedComplexityResult(items);
    }

    [McpServerTool(Name = "find_low_complexity_regions", Title = "Complexity — Low-Complexity Regions (DNA)", ReadOnly = true)]
    [Description("Entropy-thresholded contiguous low-complexity DNA regions.")]
    public static FindLowComplexityRegionsResult FindLowComplexityRegions(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size (default 64).")] int windowSize = 64,
        [Description("Entropy threshold (default 1.0).")] double entropyThreshold = 1.0)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var items = global::Seqeron.Genomics.Analysis.SequenceComplexity
            .FindLowComplexityRegions(dna, windowSize, entropyThreshold)
            .Select(r => new DnaLowComplexityItem(r.Start, r.End, r.Length, r.MinEntropy, r.Sequence))
            .ToArray();
        return new FindLowComplexityRegionsResult(items);
    }

    [McpServerTool(Name = "dust_score", Title = "Complexity — DUST Score", ReadOnly = true)]
    [Description("DUST low-complexity score (BLAST-style, triplet-based) for a DNA sequence.")]
    public static DustScoreResult DustScore(
        [Description("DNA sequence.")] string sequence,
        [Description("Word size (default 3).")] int wordSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (wordSize < 1)
            throw new ArgumentOutOfRangeException(nameof(wordSize), "Word size must be at least 1");

        var score = global::Seqeron.Genomics.Analysis.SequenceComplexity
            .CalculateDustScore(sequence, wordSize);
        return new DustScoreResult(score);
    }

    [McpServerTool(Name = "mask_low_complexity", Title = "Complexity — Mask Low-Complexity Windows", ReadOnly = true)]
    [Description("Mask low-complexity windows (DUST-driven) of a DNA sequence with a chosen character.")]
    public static MaskLowComplexityResult MaskLowComplexity(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size (default 64).")] int windowSize = 64,
        [Description("DUST threshold above which to mask (default 2.0).")] double threshold = 2.0,
        [Description("Mask character (default 'N').")] char maskChar = 'N')
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var masked = global::Seqeron.Genomics.Analysis.SequenceComplexity
            .MaskLowComplexity(dna, windowSize, threshold, maskChar);
        return new MaskLowComplexityResult(masked);
    }

    [McpServerTool(Name = "compression_ratio", Title = "Complexity — Compression Ratio", ReadOnly = true)]
    [Description("Estimate sequence repetitiveness as the normalized Lempel-Ziv complexity c/(n/log_b(n)). Lower values indicate more repetitive/less complex sequences.")]
    public static CompressionRatioResult CompressionRatio(
        [Description("Sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var ratio = global::Seqeron.Genomics.Analysis.SequenceComplexity
            .EstimateCompressionRatio(sequence);
        return new CompressionRatioResult(ratio);
    }

    #endregion

    #region ComparativeGenomics

    private static global::Seqeron.Genomics.Analysis.ComparativeGenomics.Gene ToGene(GeneInput g)
        => new(g.Id, g.GenomeId, g.Start, g.End, g.Strand, g.Sequence);

    private static IReadOnlyList<global::Seqeron.Genomics.Analysis.ComparativeGenomics.Gene> ToGenes(GeneInput[]? genes)
        => (genes ?? Array.Empty<GeneInput>()).Select(ToGene).ToList();

    private static SyntenicBlockItem MapBlock(global::Seqeron.Genomics.Analysis.ComparativeGenomics.SyntenicBlock b)
        => new(b.Genome1Id, b.Start1, b.End1, b.Genome2Id, b.Start2, b.End2, b.IsInverted, b.GeneCount, b.Identity);

    private static OrthologPairItem MapOrtholog(global::Seqeron.Genomics.Analysis.ComparativeGenomics.OrthologPair o)
        => new(o.Gene1Id, o.Gene2Id, o.Identity, o.Coverage, o.AlignmentLength);

    private static RearrangementEventItem MapRearrangement(global::Seqeron.Genomics.Analysis.ComparativeGenomics.RearrangementEvent e)
        => new(e.Type.ToString(), e.GenomeId, e.Position, e.Length, e.TargetPosition);

    [McpServerTool(Name = "find_syntenic_blocks", Title = "Comparative — Syntenic Blocks", ReadOnly = true)]
    [Description("Collinear runs of orthologous genes between two genomes.")]
    public static FindSyntenicBlocksResult FindSyntenicBlocks(
        [Description("Genes of genome 1.")] GeneInput[] genome1Genes,
        [Description("Genes of genome 2.")] GeneInput[] genome2Genes,
        [Description("Mapping from genome1 gene id to genome2 gene id.")] Dictionary<string, string> orthologMap,
        [Description("Minimum block size (default 3).")] int minBlockSize = 3,
        [Description("Maximum gap (default 5).")] int maxGap = 5)
    {
        var items = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .FindSyntenicBlocks(ToGenes(genome1Genes), ToGenes(genome2Genes),
                orthologMap ?? new Dictionary<string, string>(), minBlockSize, maxGap)
            .Select(MapBlock)
            .ToArray();
        return new FindSyntenicBlocksResult(items);
    }

    [McpServerTool(Name = "find_orthologs", Title = "Comparative — Orthologs (Best-Hit)", ReadOnly = true)]
    [Description("Best-hit ortholog pairs between two genomes by k-mer similarity (one-directional). Requires Gene.sequence populated.")]
    public static FindOrthologsResult FindOrthologs(
        [Description("Genes of genome 1 (with sequence).")] GeneInput[] genome1Genes,
        [Description("Genes of genome 2 (with sequence).")] GeneInput[] genome2Genes,
        [Description("Minimum identity (default 0.3).")] double minIdentity = 0.3,
        [Description("Minimum coverage (default 0.5).")] double minCoverage = 0.5)
    {
        var items = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .FindOrthologs(ToGenes(genome1Genes), ToGenes(genome2Genes), minIdentity, minCoverage)
            .Select(MapOrtholog)
            .ToArray();
        return new FindOrthologsResult(items);
    }

    [McpServerTool(Name = "find_reciprocal_best_hits", Title = "Comparative — Reciprocal Best Hits", ReadOnly = true)]
    [Description("Reciprocal best hits (RBH) for stricter ortholog identification. Requires Gene.sequence populated.")]
    public static FindReciprocalBestHitsResult FindReciprocalBestHits(
        [Description("Genes of genome 1 (with sequence).")] GeneInput[] genome1Genes,
        [Description("Genes of genome 2 (with sequence).")] GeneInput[] genome2Genes,
        [Description("Minimum identity (default 0.3).")] double minIdentity = 0.3)
    {
        var items = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .FindReciprocalBestHits(ToGenes(genome1Genes), ToGenes(genome2Genes), minIdentity)
            .Select(MapOrtholog)
            .ToArray();
        return new FindReciprocalBestHitsResult(items);
    }

    [McpServerTool(Name = "detect_rearrangements", Title = "Comparative — Rearrangements", ReadOnly = true)]
    [Description("Inversions, deletions, insertions inferred from gene-order disagreements between two genomes.")]
    public static DetectRearrangementsResult DetectRearrangements(
        [Description("Genes of genome 1.")] GeneInput[] genome1Genes,
        [Description("Genes of genome 2.")] GeneInput[] genome2Genes,
        [Description("Mapping from genome1 gene id to genome2 gene id.")] Dictionary<string, string> orthologMap)
    {
        if (genome1Genes is null || genome1Genes.Length == 0)
            throw new ArgumentException("Genome 1 must contain at least one gene", nameof(genome1Genes));
        if (genome2Genes is null || genome2Genes.Length == 0)
            throw new ArgumentException("Genome 2 must contain at least one gene", nameof(genome2Genes));
        if (orthologMap is null)
            throw new ArgumentException("Ortholog map cannot be null", nameof(orthologMap));

        var items = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .DetectRearrangements(ToGenes(genome1Genes), ToGenes(genome2Genes), orthologMap)
            .Select(MapRearrangement)
            .ToArray();
        return new DetectRearrangementsResult(items);
    }

    [McpServerTool(Name = "compare_genomes", Title = "Comparative — Full Comparison", ReadOnly = true)]
    [Description("End-to-end comparative pipeline: RBH orthologs + synteny + rearrangements + summary stats.")]
    public static CompareGenomesResult CompareGenomes(
        [Description("Genes of genome 1.")] GeneInput[] genome1Genes,
        [Description("Genes of genome 2.")] GeneInput[] genome2Genes,
        [Description("Minimum ortholog identity (default 0.3).")] double minOrthologIdentity = 0.3,
        [Description("Minimum syntenic block size (default 3).")] int minSyntenicBlockSize = 3)
    {
        if (genome1Genes is null || genome1Genes.Length == 0)
            throw new ArgumentException("Genome 1 must contain at least one gene", nameof(genome1Genes));
        if (genome2Genes is null || genome2Genes.Length == 0)
            throw new ArgumentException("Genome 2 must contain at least one gene", nameof(genome2Genes));

        var r = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .CompareGenomes(ToGenes(genome1Genes), ToGenes(genome2Genes), minOrthologIdentity, minSyntenicBlockSize);
        return new CompareGenomesResult(
            r.SyntenicBlocks.Select(MapBlock).ToArray(),
            r.Orthologs.Select(MapOrtholog).ToArray(),
            r.Rearrangements.Select(MapRearrangement).ToArray(),
            r.OverallSynteny,
            r.ConservedGenes,
            r.GenomeSpecificGenes1,
            r.GenomeSpecificGenes2);
    }

    [McpServerTool(Name = "reversal_distance", Title = "Comparative — Reversal Distance", ReadOnly = true)]
    [Description("Lower-bound reversal distance via breakpoint count for two equal-length permutations.")]
    public static ReversalDistanceResult ReversalDistance(
        [Description("Permutation 1.")] int[] permutation1,
        [Description("Permutation 2 (same length as permutation1).")] int[] permutation2)
    {
        var distance = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .CalculateReversalDistance(permutation1 ?? Array.Empty<int>(), permutation2 ?? Array.Empty<int>());
        return new ReversalDistanceResult(distance);
    }

    [McpServerTool(Name = "find_conserved_clusters", Title = "Comparative — Conserved Gene Clusters", ReadOnly = true)]
    [Description("Gene clusters preserved across multiple genomes. Returns lists of ortholog-group IDs per cluster.")]
    public static FindConservedClustersResult FindConservedClusters(
        [Description("Genes per genome (one array per genome).")] GeneInput[][] genomes,
        [Description("Mapping from gene id to ortholog-group id.")] Dictionary<string, string> orthologGroups,
        [Description("Minimum cluster size (default 3).")] int minClusterSize = 3,
        [Description("Maximum gap between cluster members (default 2).")] int maxGap = 2)
    {
        var genomeLists = (genomes ?? Array.Empty<GeneInput[]>())
            .Select(g => (IReadOnlyList<global::Seqeron.Genomics.Analysis.ComparativeGenomics.Gene>)ToGenes(g))
            .ToList();
        var clusters = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .FindConservedClusters(genomeLists, orthologGroups ?? new Dictionary<string, string>(), minClusterSize, maxGap)
            .Select(c => c.ToArray())
            .ToArray();
        return new FindConservedClustersResult(clusters);
    }

    [McpServerTool(Name = "calculate_ani", Title = "Comparative — Average Nucleotide Identity (ANI)", ReadOnly = true)]
    [Description("Average Nucleotide Identity (ANI) between two genome sequences (fragment-and-match).")]
    public static AniResult CalculateAni(
        [Description("Genome 1 sequence.")] string genome1Sequence,
        [Description("Genome 2 sequence.")] string genome2Sequence,
        [Description("Fragment size (default 1000).")] int fragmentSize = 1000,
        [Description("Minimum per-fragment identity (0-1) to keep a match (default 0.7).")] double minFragmentIdentity = 0.7)
    {
        if (string.IsNullOrEmpty(genome1Sequence))
            throw new ArgumentException("Genome 1 sequence cannot be null or empty", nameof(genome1Sequence));
        if (string.IsNullOrEmpty(genome2Sequence))
            throw new ArgumentException("Genome 2 sequence cannot be null or empty", nameof(genome2Sequence));
        if (fragmentSize <= 0)
            throw new ArgumentException("Fragment size must be positive", nameof(fragmentSize));

        var ani = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .CalculateANI(genome1Sequence, genome2Sequence, fragmentSize, minFragmentIdentity);
        return new AniResult(ani);
    }

    [McpServerTool(Name = "generate_dot_plot", Title = "Comparative — Dot Plot", ReadOnly = true)]
    [Description("Coordinates of matching k-mers between two sequences for dot-plot visualization.")]
    public static GenerateDotPlotResult GenerateDotPlot(
        [Description("Sequence 1.")] string sequence1,
        [Description("Sequence 2.")] string sequence2,
        [Description("Word size (default 10).")] int wordSize = 10,
        [Description("Step size (default 1).")] int stepSize = 1)
    {
        var points = global::Seqeron.Genomics.Analysis.ComparativeGenomics
            .GenerateDotPlot(sequence1 ?? string.Empty, sequence2 ?? string.Empty, wordSize, stepSize)
            .Select(p => new DotPlotPoint(p.x, p.y))
            .ToArray();
        return new GenerateDotPlotResult(points);
    }

    #endregion

    #region DisorderPredictor

    [McpServerTool(Name = "predict_disorder", Title = "Disorder — TOP-IDP Prediction", ReadOnly = true)]
    [Description("TOP-IDP disorder prediction (Campen 2008): per-residue scores plus contiguous IDR regions with confidence and subtype classification.")]
    public static PredictDisorderResult PredictDisorder(
        [Description("Protein sequence (single-letter amino acids).")] string sequence,
        [Description("Sliding window size (default 21).")] int windowSize = 21,
        [Description("Disorder threshold on TOP-IDP normalized score (default 0.542).")] double disorderThreshold = 0.542,
        [Description("Minimum length for a reported disordered region (default 5).")] int minRegionLength = 5)
    {
        var r = DisorderPredictor.PredictDisorder(sequence ?? string.Empty, windowSize, disorderThreshold, minRegionLength);
        var residues = r.ResiduePredictions
            .Select(p => new DisorderResiduePredictionItem(p.Position, p.Residue, p.DisorderScore, p.IsDisordered))
            .ToArray();
        var regions = r.DisorderedRegions
            .Select(g => new DisorderRegionItem(g.Start, g.End, g.MeanScore, g.Confidence, g.RegionType))
            .ToArray();
        return new PredictDisorderResult(r.Sequence, residues, regions, r.OverallDisorderContent, r.MeanDisorderScore);
    }

    [McpServerTool(Name = "predict_low_complexity_seg", Title = "Disorder — SEG Low-Complexity Regions", ReadOnly = true)]
    [Description("SEG algorithm (Wootton & Federhen 1993/1996) for low-complexity protein regions: trigger window K1 + extension K2.")]
    public static PredictLowComplexitySegResult PredictLowComplexitySeg(
        [Description("Protein sequence.")] string sequence,
        [Description("Trigger window length (default 12).")] int triggerWindow = 12,
        [Description("K1 trigger entropy threshold in bits (default 2.2).")] double triggerThreshold = 2.2,
        [Description("K2 extension entropy threshold in bits (default 2.5).")] double extensionThreshold = 2.5,
        [Description("Minimum reported region length (default 1).")] int minLength = 1)
    {
        var items = DisorderPredictor
            .PredictLowComplexityRegions(sequence ?? string.Empty, triggerWindow, triggerThreshold, extensionThreshold, minLength)
            .Select(t => new SegRegionItem(t.Start, t.End, t.Type))
            .ToArray();
        return new PredictLowComplexitySegResult(items);
    }

    [McpServerTool(Name = "predict_morfs", Title = "Disorder — MoRF Prediction", ReadOnly = true)]
    [Description("Predicts Molecular Recognition Features within IDRs by hydropathy enrichment (heuristic, Mohan 2006-inspired).")]
    public static PredictMorfsResult PredictMorfs(
        [Description("Protein sequence.")] string sequence,
        [Description("Minimum MoRF length (default 10).")] int minLength = 10,
        [Description("Maximum MoRF length (default 25).")] int maxLength = 25)
    {
        var items = DisorderPredictor
            .PredictMoRFs(sequence ?? string.Empty, minLength, maxLength)
            .Select(t => new MorfItem(t.Start, t.End, t.Score))
            .ToArray();
        return new PredictMorfsResult(items);
    }

    [McpServerTool(Name = "disorder_propensity", Title = "Disorder — TOP-IDP Propensity", ReadOnly = true)]
    [Description("Returns the TOP-IDP propensity value for a single amino acid (Campen 2008).")]
    public static DisorderPropensityResult DisorderPropensity(
        [Description("Single amino acid letter (length-1 string).")] string aminoAcid)
    {
        return new DisorderPropensityResult(DisorderPredictor.GetDisorderPropensity(RequireChar(aminoAcid, nameof(aminoAcid))));
    }

    [McpServerTool(Name = "is_disorder_promoting", Title = "Disorder — Is Promoting Residue", ReadOnly = true)]
    [Description("Whether an amino acid is in Dunker's disorder-promoting set {A, R, G, Q, S, P, E, K} (Dunker 2001).")]
    public static IsDisorderPromotingResult IsDisorderPromoting(
        [Description("Single amino acid letter (length-1 string).")] string aminoAcid)
    {
        return new IsDisorderPromotingResult(DisorderPredictor.IsDisorderPromoting(RequireChar(aminoAcid, nameof(aminoAcid))));
    }

    #endregion

    #region GcSkewCalculator

    [McpServerTool(Name = "gc_skew", Title = "GC Skew — Whole Sequence", ReadOnly = true)]
    [Description("Whole-sequence GC skew = (G - C) / (G + C). Range [-1, 1]; 0 if no G/C.")]
    public static GcSkewResult GcSkew(
        [Description("DNA sequence.")] string sequence)
    {
        return new GcSkewResult(GcSkewCalculator.CalculateGcSkew(sequence ?? string.Empty));
    }

    [McpServerTool(Name = "windowed_gc_skew", Title = "GC Skew — Sliding Window", ReadOnly = true)]
    [Description("Sliding-window GC skew along a sequence; each point reports center position and window bounds.")]
    public static WindowedGcSkewResult WindowedGcSkew(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size in bp (default 1000).")] int windowSize = 1000,
        [Description("Step size in bp (default 100).")] int stepSize = 100)
    {
        var items = GcSkewCalculator
            .CalculateWindowedGcSkew(sequence ?? string.Empty, windowSize, stepSize)
            .Select(p => new GcSkewPointItem(p.Position, p.GcSkew, p.WindowStart, p.WindowEnd))
            .ToArray();
        return new WindowedGcSkewResult(items);
    }

    [McpServerTool(Name = "cumulative_gc_skew", Title = "GC Skew — Cumulative", ReadOnly = true)]
    [Description("Cumulative GC skew along the sequence — minimum approximates origin and maximum approximates terminus of replication.")]
    public static CumulativeGcSkewResult CumulativeGcSkew(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size (default 1000).")] int windowSize = 1000)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), "Window size must be at least 1");

        var items = GcSkewCalculator
            .CalculateCumulativeGcSkew(sequence, windowSize)
            .Select(p => new CumulativeGcSkewPointItem(p.Position, p.GcSkew, p.CumulativeGcSkew))
            .ToArray();
        return new CumulativeGcSkewResult(items);
    }

    [McpServerTool(Name = "at_skew", Title = "AT Skew — Whole Sequence", ReadOnly = true)]
    [Description("Whole-sequence AT skew = (A - T) / (A + T). Range [-1, 1]; 0 if no A/T.")]
    public static AtSkewResult AtSkew(
        [Description("DNA sequence.")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        return new AtSkewResult(GcSkewCalculator.CalculateAtSkew(sequence));
    }

    [McpServerTool(Name = "predict_replication_origin", Title = "GC Skew — Predict Origin/Terminus", ReadOnly = true)]
    [Description("Predicts replication origin and terminus from cumulative GC skew extrema. Best on complete circular bacterial genomes.")]
    public static PredictReplicationOriginResult PredictReplicationOrigin(
        [Description("DNA sequence (ideally complete circular genome).")] string sequence)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var r = GcSkewCalculator.PredictReplicationOrigin(dna);
        return new PredictReplicationOriginResult(
            r.PredictedOrigin, r.PredictedTerminus, r.OriginSkew, r.TerminusSkew, r.IsSignificant);
    }

    [McpServerTool(Name = "analyze_gc_content", Title = "GC — Comprehensive Analysis", ReadOnly = true)]
    [Description("Comprehensive GC report: overall GC content, GC/AT skew, content/skew variances, and windowed GC profiles.")]
    public static AnalyzeGcContentResult AnalyzeGcContent(
        [Description("DNA sequence.")] string sequence,
        [Description("Window size (default 1000).")] int windowSize = 1000,
        [Description("Step size (default 100).")] int stepSize = 100)
    {
        var dna = RequireDna(sequence, nameof(sequence));
        var r = GcSkewCalculator.AnalyzeGcContent(dna, windowSize, stepSize);
        var skewPoints = r.WindowedGcSkew
            .Select(p => new GcSkewPointItem(p.Position, p.GcSkew, p.WindowStart, p.WindowEnd))
            .ToArray();
        var contentPoints = r.WindowedGcContent
            .Select(p => new GcContentPointItem(p.Position, p.GcContent, p.WindowStart, p.WindowEnd))
            .ToArray();
        return new AnalyzeGcContentResult(
            r.OverallGcContent, r.OverallGcSkew, r.OverallAtSkew,
            r.GcContentVariance, r.GcSkewVariance,
            skewPoints, contentPoints, r.SequenceLength);
    }

    #endregion

    #region RnaSecondaryStructure

    private static char RequireChar(string s, string param)
    {
        if (string.IsNullOrEmpty(s) || s.Length != 1)
            throw new ArgumentException("Expected a single character (length-1 string).", param);
        return s[0];
    }

    private static BasePairItem ToItem(RnaSecondaryStructure.BasePair bp) =>
        new(bp.Position1, bp.Position2, bp.Base1, bp.Base2, bp.Type.ToString());

    private static RnaSecondaryStructure.BasePair FromItem(BasePairItem b)
    {
        var t = b.Type switch
        {
            "WatsonCrick" => RnaSecondaryStructure.BasePairType.WatsonCrick,
            "Wobble" => RnaSecondaryStructure.BasePairType.Wobble,
            "NonCanonical" => RnaSecondaryStructure.BasePairType.NonCanonical,
            _ => throw new ArgumentException($"Unknown base pair type: {b.Type}", nameof(b))
        };
        return new RnaSecondaryStructure.BasePair(b.Position1, b.Position2, b.Base1, b.Base2, t);
    }

    private static StemItem ToItem(RnaSecondaryStructure.Stem s) =>
        new(s.Start5Prime, s.End5Prime, s.Start3Prime, s.End3Prime, s.Length,
            s.BasePairs.Select(ToItem).ToArray(), s.FreeEnergy);

    private static LoopItem ToItem(RnaSecondaryStructure.Loop l) =>
        new(l.Type.ToString(), l.Start, l.End, l.Size, l.Sequence);

    private static StemLoopItem ToItem(RnaSecondaryStructure.StemLoop sl) =>
        new(sl.Start, sl.End, ToItem(sl.Stem), ToItem(sl.Loop), sl.TotalFreeEnergy, sl.DotBracketNotation);

    private static PseudoknotItem ToItem(RnaSecondaryStructure.Pseudoknot p) =>
        new(p.Start1, p.End1, p.Start2, p.End2, p.CrossingPairs.Select(ToItem).ToArray());

    [McpServerTool(Name = "can_pair", Title = "RNA — Can Pair", ReadOnly = true)]
    [Description("Whether two RNA bases form a Watson-Crick (A-U, G-C) or wobble (G-U) pair.")]
    public static CanPairResult CanPair(
        [Description("First RNA base (length-1 string).")] string base1,
        [Description("Second RNA base (length-1 string).")] string base2)
    {
        return new CanPairResult(RnaSecondaryStructure.CanPair(RequireChar(base1, nameof(base1)), RequireChar(base2, nameof(base2))));
    }

    [McpServerTool(Name = "base_pair_type", Title = "RNA — Base-Pair Type", ReadOnly = true)]
    [Description("Classification of an RNA base-pair candidate: WatsonCrick, Wobble, or null if bases cannot pair.")]
    public static BasePairTypeResult BasePairType(
        [Description("First RNA base (length-1 string).")] string base1,
        [Description("Second RNA base (length-1 string).")] string base2)
    {
        var t = RnaSecondaryStructure.GetBasePairType(RequireChar(base1, nameof(base1)), RequireChar(base2, nameof(base2)));
        return new BasePairTypeResult(t?.ToString());
    }

    [McpServerTool(Name = "rna_complement_base", Title = "RNA — Complement Base", ReadOnly = true)]
    [Description("Returns the RNA complement (A↔U, G↔C) for a single base.")]
    public static RnaComplementBaseResult RnaComplementBase(
        [Description("RNA base (length-1 string).")] string @base)
    {
        return new RnaComplementBaseResult(RnaSecondaryStructure.GetComplement(RequireChar(@base, nameof(@base))));
    }

    [McpServerTool(Name = "find_stem_loops", Title = "RNA — Find Stem-Loops", ReadOnly = true)]
    [Description("Enumerates hairpin stem-loop candidates with stem, loop, and Turner 2004 free energy.")]
    public static FindStemLoopsResult FindStemLoops(
        [Description("RNA sequence.")] string rnaSequence,
        [Description("Minimum stem length (default 3).")] int minStemLength = 3,
        [Description("Minimum loop size (default 3).")] int minLoopSize = 3,
        [Description("Maximum loop size (default 10).")] int maxLoopSize = 10,
        [Description("Allow G-U wobble pairs in the stem (default true).")] bool allowWobble = true)
    {
        var items = RnaSecondaryStructure
            .FindStemLoops(rnaSequence ?? string.Empty, minStemLength, minLoopSize, maxLoopSize, allowWobble)
            .Select(ToItem)
            .ToArray();
        return new FindStemLoopsResult(items);
    }

    [McpServerTool(Name = "stem_energy", Title = "RNA — Stem Energy", ReadOnly = true)]
    [Description("Free energy of an RNA stem (Turner 2004 nearest-neighbor stacking + AU/GU terminal penalties).")]
    public static StemEnergyResult StemEnergy(
        [Description("RNA sequence containing the stem.")] string sequence,
        [Description("Base pairs forming the stem (5' → 3' order).")] BasePairItem[] basePairs)
    {
        var bps = (basePairs ?? Array.Empty<BasePairItem>()).Select(FromItem).ToArray();
        return new StemEnergyResult(RnaSecondaryStructure.CalculateStemEnergy(sequence ?? string.Empty, bps));
    }

    [McpServerTool(Name = "terminal_mismatch_energy", Title = "RNA — Terminal Mismatch Energy", ReadOnly = true)]
    [Description("Closing-pair × first-mismatch terminal stacking energy (Turner 2004).")]
    public static TerminalMismatchEnergyResult TerminalMismatchEnergy(
        [Description("Closing pair 5' base (length-1 string).")] string closingBase5,
        [Description("Closing pair 3' base (length-1 string).")] string closingBase3,
        [Description("5' mismatch base (length-1 string).")] string mismatch5,
        [Description("3' mismatch base (length-1 string).")] string mismatch3)
    {
        return new TerminalMismatchEnergyResult(RnaSecondaryStructure.GetTerminalMismatchEnergy(
            RequireChar(closingBase5, nameof(closingBase5)),
            RequireChar(closingBase3, nameof(closingBase3)),
            RequireChar(mismatch5, nameof(mismatch5)),
            RequireChar(mismatch3, nameof(mismatch3))));
    }

    [McpServerTool(Name = "dangling_end_energy", Title = "RNA — Dangling-End Energy", ReadOnly = true)]
    [Description("5' or 3' dangling-end stacking energy for an RNA helix end.")]
    public static DanglingEndEnergyResult DanglingEndEnergy(
        [Description("Closing pair 5' base (length-1 string).")] string closingBase5,
        [Description("Closing pair 3' base (length-1 string).")] string closingBase3,
        [Description("Dangling base (length-1 string).")] string danglingBase,
        [Description("True for 3' dangle, false for 5' dangle.")] bool is3Prime)
    {
        return new DanglingEndEnergyResult(RnaSecondaryStructure.GetDanglingEndEnergy(
            RequireChar(closingBase5, nameof(closingBase5)),
            RequireChar(closingBase3, nameof(closingBase3)),
            RequireChar(danglingBase, nameof(danglingBase)),
            is3Prime));
    }

    [McpServerTool(Name = "hairpin_loop_energy", Title = "RNA — Hairpin Loop Energy", ReadOnly = true)]
    [Description("Free energy of an RNA hairpin loop (Turner 2004 with special tri/tetra/hexaloops, terminal mismatch, all-C and special-GU adjustments).")]
    public static HairpinLoopEnergyResult HairpinLoopEnergy(
        [Description("Loop sequence (residues between the closing pair).")] string loopSequence,
        [Description("Closing pair 5' base (length-1 string).")] string closingBase5,
        [Description("Closing pair 3' base (length-1 string).")] string closingBase3,
        [Description("Apply special G-U closure bonus (default false).")] bool specialGUClosure = false)
    {
        return new HairpinLoopEnergyResult(RnaSecondaryStructure.CalculateHairpinLoopEnergy(
            loopSequence ?? string.Empty,
            RequireChar(closingBase5, nameof(closingBase5)),
            RequireChar(closingBase3, nameof(closingBase3)),
            specialGUClosure));
    }

    [McpServerTool(Name = "internal_loop_energy", Title = "RNA — Internal Loop Energy", ReadOnly = true)]
    [Description("Free energy of a generic RNA internal loop (Turner 2004; uses 1×1 lookup table when n1=n2=1).")]
    public static InternalLoopEnergyResult InternalLoopEnergy(
        [Description("Unpaired bases on the 5' side.")] int n1,
        [Description("Unpaired bases on the 3' side.")] int n2,
        [Description("Outer closing pair 5' base.")] string closingBase5_1,
        [Description("Outer closing pair 3' base.")] string closingBase3_1,
        [Description("Inner closing pair 5' base.")] string closingBase5_2,
        [Description("Inner closing pair 3' base.")] string closingBase3_2,
        [Description("First unpaired base adjacent to outer closing pair.")] string mismatch5_1,
        [Description("Last unpaired base adjacent to outer closing pair.")] string mismatch3_1,
        [Description("First unpaired base adjacent to inner closing pair.")] string mismatch5_2,
        [Description("Last unpaired base adjacent to inner closing pair.")] string mismatch3_2)
    {
        return new InternalLoopEnergyResult(RnaSecondaryStructure.CalculateInternalLoopEnergy(
            n1, n2,
            RequireChar(closingBase5_1, nameof(closingBase5_1)),
            RequireChar(closingBase3_1, nameof(closingBase3_1)),
            RequireChar(closingBase5_2, nameof(closingBase5_2)),
            RequireChar(closingBase3_2, nameof(closingBase3_2)),
            RequireChar(mismatch5_1, nameof(mismatch5_1)),
            RequireChar(mismatch3_1, nameof(mismatch3_1)),
            RequireChar(mismatch5_2, nameof(mismatch5_2)),
            RequireChar(mismatch3_2, nameof(mismatch3_2))));
    }

    [McpServerTool(Name = "bulge_loop_energy", Title = "RNA — Bulge Loop Energy", ReadOnly = true)]
    [Description("Free energy of an RNA bulge loop (special-C bonus, n=1 stacking, degeneracy entropy from numStates).")]
    public static BulgeLoopEnergyResult BulgeLoopEnergy(
        [Description("Number of unpaired bases in the bulge.")] int bulgeSize,
        [Description("Bulged base (length-1 string).")] string bulgedBase,
        [Description("5' base of the pair on the 5' side.")] string pair5_base1,
        [Description("3' base of the pair on the 5' side.")] string pair5_base2,
        [Description("5' base of the pair on the 3' side.")] string pair3_base1,
        [Description("3' base of the pair on the 3' side.")] string pair3_base2,
        [Description("Number of equivalent states for degeneracy entropy (default 1).")] int numStates = 1)
    {
        return new BulgeLoopEnergyResult(RnaSecondaryStructure.CalculateBulgeLoopEnergy(
            bulgeSize,
            RequireChar(bulgedBase, nameof(bulgedBase)),
            RequireChar(pair5_base1, nameof(pair5_base1)),
            RequireChar(pair5_base2, nameof(pair5_base2)),
            RequireChar(pair3_base1, nameof(pair3_base1)),
            RequireChar(pair3_base2, nameof(pair3_base2)),
            numStates));
    }

    [McpServerTool(Name = "multibranch_loop_energy", Title = "RNA — Multibranch Loop Energy", ReadOnly = true)]
    [Description("Free energy of an RNA multibranch loop (Turner 2004 affine model: offset + asymmetry + helix term + stacking + strain).")]
    public static MultibranchLoopEnergyResult MultibranchLoopEnergy(
        [Description("Number of helical branches.")] int numHelices,
        [Description("Total unpaired nucleotides in the loop.")] int numUnpaired,
        [Description("True if the junction has steric strain (default false).")] bool hasStrain = false,
        [Description("Pre-computed optimal stacking/dangling end energy (default 0).")] double stackingEnergy = 0.0)
    {
        return new MultibranchLoopEnergyResult(RnaSecondaryStructure.CalculateMultibranchLoopEnergy(
            numHelices, numUnpaired, hasStrain, stackingEnergy));
    }

    [McpServerTool(Name = "flush_coaxial_stacking", Title = "RNA — Flush Coaxial Stacking", ReadOnly = true)]
    [Description("Coaxial stacking energy for two RNA helices with no intervening unpaired bases.")]
    public static FlushCoaxialStackingResult FlushCoaxialStacking(
        [Description("5' base of first helix end.")] string base5_1,
        [Description("3' base of first helix end.")] string base3_1,
        [Description("5' base of second helix end.")] string base5_2,
        [Description("3' base of second helix end.")] string base3_2)
    {
        return new FlushCoaxialStackingResult(RnaSecondaryStructure.CalculateFlushCoaxialStacking(
            RequireChar(base5_1, nameof(base5_1)),
            RequireChar(base3_1, nameof(base3_1)),
            RequireChar(base5_2, nameof(base5_2)),
            RequireChar(base3_2, nameof(base3_2))));
    }

    [McpServerTool(Name = "mismatch_coaxial_stacking", Title = "RNA — Mismatch Coaxial Stacking", ReadOnly = true)]
    [Description("Mismatch-mediated coaxial stacking energy: terminal mismatch + base + WC/GU bonus.")]
    public static MismatchCoaxialStackingResult MismatchCoaxialStacking(
        [Description("Closing pair 5' base.")] string closingBase5,
        [Description("Closing pair 3' base.")] string closingBase3,
        [Description("5' mismatch base.")] string mismatch5,
        [Description("3' mismatch base.")] string mismatch3)
    {
        return new MismatchCoaxialStackingResult(RnaSecondaryStructure.CalculateMismatchCoaxialStacking(
            RequireChar(closingBase5, nameof(closingBase5)),
            RequireChar(closingBase3, nameof(closingBase3)),
            RequireChar(mismatch5, nameof(mismatch5)),
            RequireChar(mismatch3, nameof(mismatch3))));
    }

    [McpServerTool(Name = "minimum_free_energy", Title = "RNA — Minimum Free Energy (Zuker)", ReadOnly = true)]
    [Description("Zuker-style minimum free energy with Turner 2004 parameters (O(n³)). Returns kcal/mol.")]
    public static MinimumFreeEnergyResult MinimumFreeEnergy(
        [Description("RNA sequence.")] string rnaSequence,
        [Description("Minimum hairpin loop size (default 3; values <3 are clamped).")] int minLoopSize = 3)
    {
        return new MinimumFreeEnergyResult(RnaSecondaryStructure.CalculateMinimumFreeEnergy(rnaSequence ?? string.Empty, minLoopSize));
    }

    [McpServerTool(Name = "predict_rna_structure", Title = "RNA — Predict Secondary Structure", ReadOnly = true)]
    [Description("Greedy non-overlapping stem-loop selection — produces dot-bracket notation, base pairs, stems, pseudoknots, and total MFE.")]
    public static PredictRnaStructureResult PredictRnaStructure(
        [Description("RNA sequence.")] string rnaSequence,
        [Description("Minimum stem length (default 3).")] int minStemLength = 3,
        [Description("Minimum loop size (default 3).")] int minLoopSize = 3,
        [Description("Maximum loop size (default 10).")] int maxLoopSize = 10)
    {
        var s = RnaSecondaryStructure.PredictStructure(rnaSequence ?? string.Empty, minStemLength, minLoopSize, maxLoopSize);
        return new PredictRnaStructureResult(
            s.Sequence,
            s.DotBracket,
            s.BasePairs.Select(ToItem).ToArray(),
            s.StemLoops.Select(ToItem).ToArray(),
            s.Pseudoknots.Select(ToItem).ToArray(),
            s.MinimumFreeEnergy);
    }

    [McpServerTool(Name = "detect_pseudoknots", Title = "RNA — Detect Pseudoknots", ReadOnly = true)]
    [Description("Identifies pseudoknots as crossing base pairs: i < i' < j < j' for pairs (i,j) and (i',j').")]
    public static DetectPseudoknotsResult DetectPseudoknots(
        [Description("Base pairs to inspect for crossings.")] BasePairItem[] basePairs)
    {
        var bps = (basePairs ?? Array.Empty<BasePairItem>()).Select(FromItem).ToArray();
        var items = RnaSecondaryStructure.DetectPseudoknots(bps).Select(ToItem).ToArray();
        return new DetectPseudoknotsResult(items);
    }

    [McpServerTool(Name = "parse_dot_bracket", Title = "RNA — Parse Dot-Bracket", ReadOnly = true)]
    [Description("Parses dot-bracket notation into a list of base-pair coordinates. Supports ()/[]/{}/<> bracket families.")]
    public static ParseDotBracketResult ParseDotBracket(
        [Description("Dot-bracket secondary-structure string.")] string dotBracket)
    {
        var pairs = RnaSecondaryStructure
            .ParseDotBracket(dotBracket ?? string.Empty)
            .Select(t => new BasePairCoord(t.Position1, t.Position2))
            .ToArray();
        return new ParseDotBracketResult(pairs);
    }

    [McpServerTool(Name = "validate_dot_bracket", Title = "RNA — Validate Dot-Bracket", ReadOnly = true)]
    [Description("Validates that all bracket symbols in a dot-bracket string are balanced.")]
    public static ValidateDotBracketResult ValidateDotBracket(
        [Description("Dot-bracket secondary-structure string.")] string dotBracket)
    {
        return new ValidateDotBracketResult(RnaSecondaryStructure.ValidateDotBracket(dotBracket ?? string.Empty));
    }

    [McpServerTool(Name = "find_rna_inverted_repeats", Title = "RNA — Find Inverted Repeats", ReadOnly = true)]
    [Description("Finds antiparallel complementary regions (potential RNA hairpin stems). Distinct from DNA find_inverted_repeats.")]
    public static FindRnaInvertedRepeatsResult FindRnaInvertedRepeats(
        [Description("RNA sequence.")] string sequence,
        [Description("Minimum arm length (default 4).")] int minLength = 4,
        [Description("Minimum spacing between arms (default 3).")] int minSpacing = 3,
        [Description("Maximum spacing between arms (default 100).")] int maxSpacing = 100)
    {
        var items = RnaSecondaryStructure
            .FindInvertedRepeats(sequence ?? string.Empty, minLength, minSpacing, maxSpacing)
            .Select(t => new FindRnaInvertedRepeatsItem(t.Start1, t.End1, t.Start2, t.End2, t.Length))
            .ToArray();
        return new FindRnaInvertedRepeatsResult(items);
    }

    #endregion
}
