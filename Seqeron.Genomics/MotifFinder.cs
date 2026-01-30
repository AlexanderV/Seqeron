using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;

namespace Seqeron.Genomics;

/// <summary>
/// Finds conserved motifs and patterns in DNA sequences.
/// Supports exact and degenerate motif searching, position weight matrices, and consensus sequences.
/// </summary>
public static class MotifFinder
{
    #region Exact Motif Finding

    /// <summary>
    /// Finds all occurrences of an exact motif in a sequence.
    /// Uses SuffixTree for O(m+k) pattern matching where m=motif length, k=occurrences.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern to find.</param>
    /// <returns>Positions where the motif occurs.</returns>
    public static IEnumerable<int> FindExactMotif(DnaSequence sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (string.IsNullOrEmpty(motif)) yield break;

        string motifUpper = motif.ToUpperInvariant();

        // Use SuffixTree for efficient pattern matching
        var positions = sequence.SuffixTree.FindAllOccurrences(motifUpper);
        foreach (int pos in positions.OrderBy(p => p))
        {
            yield return pos;
        }
    }

    #endregion

    #region Degenerate Motif Finding

    /// <summary>
    /// IUPAC nucleotide codes for degenerate bases.
    /// </summary>
    private static readonly Dictionary<char, string> IupacCodes = new()
    {
        ['A'] = "A",
        ['T'] = "T",
        ['G'] = "G",
        ['C'] = "C",
        ['R'] = "AG",   // Purine
        ['Y'] = "CT",   // Pyrimidine
        ['S'] = "GC",   // Strong
        ['W'] = "AT",   // Weak
        ['K'] = "GT",   // Keto
        ['M'] = "AC",   // Amino
        ['B'] = "CGT",  // Not A
        ['D'] = "AGT",  // Not C
        ['H'] = "ACT",  // Not G
        ['V'] = "ACG",  // Not T
        ['N'] = "ACGT", // Any
    };

    /// <summary>
    /// Finds all occurrences of a degenerate motif using IUPAC codes.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern with IUPAC ambiguity codes.</param>
    /// <returns>Motif matches with positions.</returns>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(DnaSequence sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (string.IsNullOrEmpty(motif)) yield break;

        string seq = sequence.Sequence;
        string motifUpper = motif.ToUpperInvariant();

        for (int i = 0; i <= seq.Length - motifUpper.Length; i++)
        {
            bool matches = true;
            for (int j = 0; j < motifUpper.Length && matches; j++)
            {
                char motifChar = motifUpper[j];
                char seqChar = seq[i + j];

                if (IupacCodes.TryGetValue(motifChar, out string? allowed))
                {
                    matches = allowed.Contains(seqChar);
                }
                else
                {
                    matches = motifChar == seqChar;
                }
            }

            if (matches)
            {
                yield return new MotifMatch(
                    Position: i,
                    MatchedSequence: seq.Substring(i, motifUpper.Length),
                    Pattern: motifUpper,
                    Score: 1.0);
            }
        }
    }

    /// <summary>
    /// Finds all occurrences of a degenerate motif with cancellation support.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern with IUPAC ambiguity codes.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <returns>Motif matches with positions.</returns>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(
        DnaSequence sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CancellableOperations.FindDegenerateMotif(sequence.Sequence, motif, cancellationToken);
    }

    /// <summary>
    /// Finds degenerate motif in a raw string with cancellation support.
    /// </summary>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(
        string sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        return CancellableOperations.FindDegenerateMotif(sequence, motif, cancellationToken);
    }

    #endregion

    #region Position Weight Matrix

    /// <summary>
    /// Creates a Position Weight Matrix (PWM) from aligned sequences.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <param name="pseudocount">Pseudocount for smoothing (default: 0.25).</param>
    /// <returns>Position Weight Matrix.</returns>
    public static PositionWeightMatrix CreatePwm(IEnumerable<string> sequences, double pseudocount = 0.25)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.Select(s => s.ToUpperInvariant()).ToList();
        if (seqList.Count == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(sequences));

        int length = seqList[0].Length;
        if (!seqList.All(s => s.Length == length))
            throw new ArgumentException("All sequences must have the same length.", nameof(sequences));

        var matrix = new double[4, length]; // A, C, G, T
        int count = seqList.Count;

        // Count bases at each position
        foreach (var seq in seqList)
        {
            for (int i = 0; i < length; i++)
            {
                int baseIndex = seq[i] switch
                {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    _ => -1
                };

                if (baseIndex >= 0)
                    matrix[baseIndex, i]++;
            }
        }

        // Convert to log-odds with pseudocounts
        double background = 0.25; // Equal background frequencies
        for (int i = 0; i < length; i++)
        {
            for (int b = 0; b < 4; b++)
            {
                double freq = (matrix[b, i] + pseudocount) / (count + 4 * pseudocount);
                matrix[b, i] = Math.Log2(freq / background);
            }
        }

        return new PositionWeightMatrix(matrix, length);
    }

    /// <summary>
    /// Scans a sequence with a PWM and returns matches above threshold.
    /// </summary>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <param name="pwm">Position Weight Matrix.</param>
    /// <param name="threshold">Minimum score threshold.</param>
    /// <returns>Matches with scores.</returns>
    public static IEnumerable<MotifMatch> ScanWithPwm(
        DnaSequence sequence,
        PositionWeightMatrix pwm,
        double threshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);

        string seq = sequence.Sequence;
        int motifLen = pwm.Length;

        for (int i = 0; i <= seq.Length - motifLen; i++)
        {
            double score = 0;
            bool valid = true;

            for (int j = 0; j < motifLen && valid; j++)
            {
                int baseIndex = seq[i + j] switch
                {
                    'A' => 0,
                    'C' => 1,
                    'G' => 2,
                    'T' => 3,
                    _ => -1
                };

                if (baseIndex < 0)
                {
                    valid = false;
                }
                else
                {
                    score += pwm.Matrix[baseIndex, j];
                }
            }

            if (valid && score >= threshold)
            {
                yield return new MotifMatch(
                    Position: i,
                    MatchedSequence: seq.Substring(i, motifLen),
                    Pattern: pwm.Consensus,
                    Score: score);
            }
        }
    }

    #endregion

    #region Consensus Sequence

    /// <summary>
    /// Generates a consensus sequence from aligned sequences.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <returns>Consensus sequence using IUPAC codes.</returns>
    public static string GenerateConsensus(IEnumerable<string> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.Select(s => s.ToUpperInvariant()).ToList();
        if (seqList.Count == 0) return "";

        int length = seqList[0].Length;
        var consensus = new StringBuilder(length);

        for (int i = 0; i < length; i++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0 };

            foreach (var seq in seqList)
            {
                if (i < seq.Length && counts.ContainsKey(seq[i]))
                    counts[seq[i]]++;
            }

            consensus.Append(GetIupacCode(counts, seqList.Count));
        }

        return consensus.ToString();
    }

    private static char GetIupacCode(Dictionary<char, int> counts, int total)
    {
        double threshold = total * 0.25; // Base must be > 25% to be included

        var present = counts.Where(kv => kv.Value > threshold)
                           .Select(kv => kv.Key)
                           .OrderBy(c => c)
                           .ToList();

        if (present.Count == 0)
        {
            return counts.MaxBy(kv => kv.Value).Key;
        }

        string bases = string.Join("", present);

        return bases switch
        {
            "A" => 'A',
            "C" => 'C',
            "G" => 'G',
            "T" => 'T',
            "AG" => 'R',
            "CT" => 'Y',
            "CG" => 'S',
            "AT" => 'W',
            "GT" => 'K',
            "AC" => 'M',
            "CGT" => 'B',
            "AGT" => 'D',
            "ACT" => 'H',
            "ACG" => 'V',
            _ => 'N'
        };
    }

    #endregion

    #region Motif Discovery

    /// <summary>
    /// Discovers overrepresented k-mers that may represent motifs.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="k">K-mer length (default: 6).</param>
    /// <param name="minCount">Minimum occurrence count (default: 2).</param>
    /// <returns>Overrepresented k-mers with their counts and positions.</returns>
    public static IEnumerable<DiscoveredMotif> DiscoverMotifs(
        DnaSequence sequence,
        int k = 6,
        int minCount = 2)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

        string seq = sequence.Sequence;
        var kmerPositions = new Dictionary<string, List<int>>();

        // Count k-mers
        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);

            if (!kmerPositions.ContainsKey(kmer))
                kmerPositions[kmer] = new List<int>();

            kmerPositions[kmer].Add(i);
        }

        // Calculate expected frequency
        double expectedFreq = seq.Length - k + 1.0;
        double expectedCount = expectedFreq / Math.Pow(4, k);

        // Return overrepresented k-mers
        foreach (var (kmer, positions) in kmerPositions)
        {
            if (positions.Count >= minCount)
            {
                double enrichment = positions.Count / Math.Max(expectedCount, 0.1);

                yield return new DiscoveredMotif(
                    Sequence: kmer,
                    Count: positions.Count,
                    Positions: positions.AsReadOnly(),
                    Enrichment: enrichment);
            }
        }
    }

    /// <summary>
    /// Finds common motifs shared between multiple sequences.
    /// </summary>
    /// <param name="sequences">Collection of DNA sequences.</param>
    /// <param name="k">K-mer length (default: 6).</param>
    /// <param name="minSequences">Minimum sequences containing the motif (default: 2).</param>
    /// <returns>Shared motifs with occurrence information.</returns>
    public static IEnumerable<SharedMotif> FindSharedMotifs(
        IEnumerable<DnaSequence> sequences,
        int k = 6,
        int minSequences = 2)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k));

        var seqList = sequences.ToList();
        var kmerOccurrences = new Dictionary<string, List<int>>();

        // Find k-mers in each sequence
        for (int seqIdx = 0; seqIdx < seqList.Count; seqIdx++)
        {
            var seq = seqList[seqIdx].Sequence;
            var seenInSeq = new HashSet<string>();

            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);

                if (!seenInSeq.Contains(kmer))
                {
                    seenInSeq.Add(kmer);

                    if (!kmerOccurrences.ContainsKey(kmer))
                        kmerOccurrences[kmer] = new List<int>();

                    kmerOccurrences[kmer].Add(seqIdx);
                }
            }
        }

        // Return shared motifs
        foreach (var (kmer, seqIndices) in kmerOccurrences)
        {
            if (seqIndices.Count >= minSequences)
            {
                yield return new SharedMotif(
                    Sequence: kmer,
                    SequenceIndices: seqIndices.AsReadOnly(),
                    Prevalence: (double)seqIndices.Count / seqList.Count);
            }
        }
    }

    #endregion

    #region Regulatory Motif Patterns

    /// <summary>
    /// Known regulatory motif patterns.
    /// </summary>
    public static class KnownMotifs
    {
        /// <summary>TATA box consensus: TATAAA</summary>
        public const string TataBox = "TATAAA";

        /// <summary>CAAT box consensus: CCAAT</summary>
        public const string CaatBox = "CCAAT";

        /// <summary>GC box consensus: GGGCGG</summary>
        public const string GcBox = "GGGCGG";

        /// <summary>Kozak consensus: GCCGCCACC (around start codon)</summary>
        public const string Kozak = "GCCGCCACCATGG";

        /// <summary>Shine-Dalgarno (bacterial RBS): AGGAGG</summary>
        public const string ShineDalgarno = "AGGAGG";

        /// <summary>Poly(A) signal: AATAAA</summary>
        public const string PolyASignal = "AATAAA";

        /// <summary>E-box consensus: CANNTG (using IUPAC)</summary>
        public const string EBox = "CANNTG";

        /// <summary>AP-1 binding site: TGAGTCA</summary>
        public const string Ap1 = "TGAGTCA";

        /// <summary>NF-kB consensus: GGGACTTTCC</summary>
        public const string NfKb = "GGGACTTTCC";

        /// <summary>CREB binding site: TGACGTCA</summary>
        public const string Creb = "TGACGTCA";
    }

    /// <summary>
    /// Scans for known regulatory motifs.
    /// </summary>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <returns>Found regulatory elements.</returns>
    public static IEnumerable<RegulatoryElement> FindRegulatoryElements(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var patterns = new (string Name, string Pattern, string Description)[]
        {
            ("TATA Box", KnownMotifs.TataBox, "Core promoter element"),
            ("CAAT Box", KnownMotifs.CaatBox, "Promoter element"),
            ("GC Box", KnownMotifs.GcBox, "Sp1 binding site"),
            ("Kozak", KnownMotifs.Kozak, "Translation initiation"),
            ("Shine-Dalgarno", KnownMotifs.ShineDalgarno, "Bacterial ribosome binding"),
            ("Poly(A) Signal", KnownMotifs.PolyASignal, "Polyadenylation signal"),
            ("E-box", KnownMotifs.EBox, "bHLH transcription factor binding"),
            ("AP-1", KnownMotifs.Ap1, "AP-1 transcription factor binding"),
            ("NF-κB", KnownMotifs.NfKb, "NF-κB binding site"),
            ("CREB", KnownMotifs.Creb, "CREB transcription factor binding")
        };

        foreach (var (name, pattern, description) in patterns)
        {
            foreach (var match in FindDegenerateMotif(sequence, pattern))
            {
                yield return new RegulatoryElement(
                    Name: name,
                    Position: match.Position,
                    Sequence: match.MatchedSequence,
                    Pattern: pattern,
                    Description: description);
            }
        }
    }

    #endregion
}

/// <summary>
/// A motif match in a sequence.
/// </summary>
public readonly record struct MotifMatch(
    int Position,
    string MatchedSequence,
    string Pattern,
    double Score);

/// <summary>
/// A discovered motif from de novo discovery.
/// </summary>
public readonly record struct DiscoveredMotif(
    string Sequence,
    int Count,
    IReadOnlyList<int> Positions,
    double Enrichment);

/// <summary>
/// A motif shared between multiple sequences.
/// </summary>
public readonly record struct SharedMotif(
    string Sequence,
    IReadOnlyList<int> SequenceIndices,
    double Prevalence);

/// <summary>
/// A regulatory element found in a sequence.
/// </summary>
public readonly record struct RegulatoryElement(
    string Name,
    int Position,
    string Sequence,
    string Pattern,
    string Description);

/// <summary>
/// Position Weight Matrix for motif scoring.
/// </summary>
public sealed class PositionWeightMatrix
{
    public double[,] Matrix { get; }
    public int Length { get; }
    public string Consensus { get; }

    public PositionWeightMatrix(double[,] matrix, int length)
    {
        Matrix = matrix;
        Length = length;
        Consensus = GenerateConsensus();
    }

    private string GenerateConsensus()
    {
        var sb = new StringBuilder(Length);
        char[] bases = { 'A', 'C', 'G', 'T' };

        for (int i = 0; i < Length; i++)
        {
            int maxIdx = 0;
            double maxVal = Matrix[0, i];

            for (int b = 1; b < 4; b++)
            {
                if (Matrix[b, i] > maxVal)
                {
                    maxVal = Matrix[b, i];
                    maxIdx = b;
                }
            }

            sb.Append(bases[maxIdx]);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Gets the maximum possible score for this PWM.
    /// </summary>
    public double MaxScore
    {
        get
        {
            double max = 0;
            for (int i = 0; i < Length; i++)
            {
                double posMax = double.MinValue;
                for (int b = 0; b < 4; b++)
                    posMax = Math.Max(posMax, Matrix[b, i]);
                max += posMax;
            }
            return max;
        }
    }

    /// <summary>
    /// Gets the minimum possible score for this PWM.
    /// </summary>
    public double MinScore
    {
        get
        {
            double min = 0;
            for (int i = 0; i < Length; i++)
            {
                double posMin = double.MaxValue;
                for (int b = 0; b < 4; b++)
                    posMin = Math.Min(posMin, Matrix[b, i]);
                min += posMin;
            }
            return min;
        }
    }
}
