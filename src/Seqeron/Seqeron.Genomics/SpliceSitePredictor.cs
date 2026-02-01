using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Provides algorithms for RNA splice site prediction and intron/exon analysis.
/// </summary>
public static class SpliceSitePredictor
{
    #region Records and Types

    /// <summary>
    /// Represents a predicted splice site.
    /// </summary>
    public readonly record struct SpliceSite(
        int Position,
        SpliceSiteType Type,
        string Motif,
        double Score,
        double Confidence);

    /// <summary>
    /// Types of splice sites.
    /// </summary>
    public enum SpliceSiteType
    {
        Donor,      // 5' splice site (GT/GU)
        Acceptor,   // 3' splice site (AG)
        Branch,     // Branch point (A)
        U12Donor,   // Minor spliceosome donor (AT/AU)
        U12Acceptor // Minor spliceosome acceptor (AC)
    }

    /// <summary>
    /// Represents a predicted intron.
    /// </summary>
    public readonly record struct Intron(
        int Start,
        int End,
        int Length,
        SpliceSite DonorSite,
        SpliceSite AcceptorSite,
        SpliceSite? BranchPoint,
        string Sequence,
        IntronType Type,
        double Score);

    /// <summary>
    /// Types of introns.
    /// </summary>
    public enum IntronType
    {
        U2,     // Major spliceosome (GT-AG)
        U12,    // Minor spliceosome (AT-AC)
        GcAg,   // GC-AG variant
        Unknown
    }

    /// <summary>
    /// Represents an exon.
    /// </summary>
    public readonly record struct Exon(
        int Start,
        int End,
        int Length,
        ExonType Type,
        int? Phase,
        string Sequence);

    /// <summary>
    /// Types of exons.
    /// </summary>
    public enum ExonType
    {
        Initial,    // First exon (contains start codon)
        Internal,   // Internal exon
        Terminal,   // Last exon (contains stop codon)
        Single      // Single exon gene
    }

    /// <summary>
    /// Complete gene structure prediction.
    /// </summary>
    public readonly record struct GeneStructure(
        IReadOnlyList<Exon> Exons,
        IReadOnlyList<Intron> Introns,
        string SplicedSequence,
        double OverallScore);

    #endregion

    #region Position Weight Matrices

    // Consensus donor site: MAG|GURAGU (M=A/C, R=A/G)
    // Position -3 to +6 relative to splice site
    private static readonly Dictionary<int, Dictionary<char, double>> DonorPwm = new()
    {
        { -3, new Dictionary<char, double> { {'A', 0.35}, {'C', 0.35}, {'G', 0.15}, {'U', 0.15} } },
        { -2, new Dictionary<char, double> { {'A', 0.60}, {'C', 0.10}, {'G', 0.10}, {'U', 0.20} } },
        { -1, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.05}, {'G', 0.80}, {'U', 0.05} } },
        { 0,  new Dictionary<char, double> { {'A', 0.00}, {'C', 0.00}, {'G', 1.00}, {'U', 0.00} } }, // G
        { 1,  new Dictionary<char, double> { {'A', 0.00}, {'C', 0.00}, {'G', 0.00}, {'U', 1.00} } }, // U
        { 2,  new Dictionary<char, double> { {'A', 0.60}, {'C', 0.05}, {'G', 0.30}, {'U', 0.05} } },
        { 3,  new Dictionary<char, double> { {'A', 0.70}, {'C', 0.05}, {'G', 0.20}, {'U', 0.05} } },
        { 4,  new Dictionary<char, double> { {'A', 0.10}, {'C', 0.05}, {'G', 0.80}, {'U', 0.05} } },
        { 5,  new Dictionary<char, double> { {'A', 0.15}, {'C', 0.15}, {'G', 0.15}, {'U', 0.55} } }
    };

    // Consensus acceptor site: (Y)nNCAG|G
    // Position -15 to +1 relative to splice site
    private static readonly Dictionary<int, Dictionary<char, double>> AcceptorPwm = new()
    {
        { -15, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.30}, {'G', 0.10}, {'U', 0.50} } },
        { -10, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.30}, {'G', 0.10}, {'U', 0.50} } },
        { -5,  new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { -4,  new Dictionary<char, double> { {'A', 0.05}, {'C', 0.40}, {'G', 0.05}, {'U', 0.50} } },
        { -3,  new Dictionary<char, double> { {'A', 0.05}, {'C', 0.70}, {'G', 0.05}, {'U', 0.20} } },
        { -2,  new Dictionary<char, double> { {'A', 1.00}, {'C', 0.00}, {'G', 0.00}, {'U', 0.00} } }, // A
        { -1,  new Dictionary<char, double> { {'A', 0.00}, {'C', 0.00}, {'G', 1.00}, {'U', 0.00} } }, // G
        { 0,   new Dictionary<char, double> { {'A', 0.20}, {'C', 0.15}, {'G', 0.50}, {'U', 0.15} } }
    };

    // Branch point consensus: YNYURAC (Y=C/U, R=A/G, N=any)
    // Typically 18-40 nt upstream of 3' splice site
    private static readonly Dictionary<int, Dictionary<char, double>> BranchPointPwm = new()
    {
        { 0, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { 1, new Dictionary<char, double> { {'A', 0.25}, {'C', 0.25}, {'G', 0.25}, {'U', 0.25} } },
        { 2, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.40}, {'G', 0.10}, {'U', 0.40} } },
        { 3, new Dictionary<char, double> { {'A', 0.05}, {'C', 0.05}, {'G', 0.05}, {'U', 0.85} } },
        { 4, new Dictionary<char, double> { {'A', 0.60}, {'C', 0.05}, {'G', 0.30}, {'U', 0.05} } },
        { 5, new Dictionary<char, double> { {'A', 0.95}, {'C', 0.02}, {'G', 0.02}, {'U', 0.01} } }, // Branch A
        { 6, new Dictionary<char, double> { {'A', 0.10}, {'C', 0.60}, {'G', 0.10}, {'U', 0.20} } }
    };

    #endregion

    #region Splice Site Prediction

    /// <summary>
    /// Finds all potential donor (5') splice sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindDonorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 6)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 0; i <= upper.Length - 6; i++)
        {
            // Check for canonical GT/GU
            if (upper[i] == 'G' && upper[i + 1] == 'U')
            {
                double score = ScoreDonorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // Non-canonical GC
            else if (includeNonCanonical && upper[i] == 'G' && upper[i + 1] == 'C')
            {
                double score = ScoreDonorSite(upper, i) * 0.7; // Penalty for GC
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // U12 spliceosome AT/AU
            else if (includeNonCanonical && upper[i] == 'A' && (upper[i + 1] == 'U' || upper[i + 1] == 'T'))
            {
                double score = ScoreU12DonorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i,
                        Type: SpliceSiteType.U12Donor,
                        Motif: GetMotifContext(upper, i, 3, 6),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Finds all potential acceptor (3') splice sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindAcceptorSites(
        string sequence,
        double minScore = 0.5,
        bool includeNonCanonical = false)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 20)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 15; i < upper.Length - 1; i++)
        {
            // Check for canonical AG
            if (upper[i] == 'A' && upper[i + 1] == 'G')
            {
                double score = ScoreAcceptorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i + 1, // Position after the AG
                        Type: SpliceSiteType.Acceptor,
                        Motif: GetMotifContext(upper, i, 15, 2),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
            // U12 spliceosome AC
            else if (includeNonCanonical && upper[i] == 'A' && upper[i + 1] == 'C')
            {
                double score = ScoreU12AcceptorSite(upper, i);
                if (score >= minScore)
                {
                    yield return new SpliceSite(
                        Position: i + 1,
                        Type: SpliceSiteType.U12Acceptor,
                        Motif: GetMotifContext(upper, i, 15, 2),
                        Score: score,
                        Confidence: CalculateConfidence(score, 0.5, 1.0));
                }
            }
        }
    }

    /// <summary>
    /// Finds branch point sites in a sequence.
    /// </summary>
    public static IEnumerable<SpliceSite> FindBranchPoints(
        string sequence,
        int searchStart = 0,
        int searchEnd = -1,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 7)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        int end = searchEnd < 0 ? upper.Length - 7 : Math.Min(searchEnd, upper.Length - 7);

        for (int i = Math.Max(0, searchStart); i <= end; i++)
        {
            double score = ScoreBranchPoint(upper, i);
            if (score >= minScore)
            {
                yield return new SpliceSite(
                    Position: i + 5, // Position of branch A
                    Type: SpliceSiteType.Branch,
                    Motif: upper.Substring(i, 7),
                    Score: score,
                    Confidence: CalculateConfidence(score, 0.5, 1.0));
            }
        }
    }

    #endregion

    #region Scoring Functions

    private static double ScoreDonorSite(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        foreach (var (offset, weights) in DonorPwm)
        {
            int pos = position + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += Math.Log2(weight / 0.25 + 0.01); // Log-odds with pseudocount
                    count++;
                }
            }
        }

        // Normalize to 0-1 range
        double normalized = (score / count + 2) / 4; // Approximate normalization
        return Math.Max(0, Math.Min(1, normalized));
    }

    private static double ScoreAcceptorSite(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        // Score polypyrimidine tract
        int pptScore = 0;
        for (int i = position - 15; i < position - 3; i++)
        {
            if (i >= 0 && i < sequence.Length)
            {
                if (sequence[i] == 'C' || sequence[i] == 'U')
                    pptScore++;
            }
        }

        score += pptScore / 12.0 * 2; // PPT contribution

        // Score AG and context
        foreach (var (offset, weights) in AcceptorPwm)
        {
            int pos = position + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += Math.Log2(weight / 0.25 + 0.01);
                    count++;
                }
            }
        }

        double normalized = (score / (count + 1) + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    private static double ScoreBranchPoint(string sequence, int position)
    {
        double score = 0;
        int count = 0;

        foreach (var (offset, weights) in BranchPointPwm)
        {
            int pos = position + offset;
            if (pos >= 0 && pos < sequence.Length)
            {
                char b = sequence[pos];
                if (weights.TryGetValue(b, out double weight))
                {
                    score += Math.Log2(weight / 0.25 + 0.01);
                    count++;
                }
            }
        }

        double normalized = (score / count + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    private static double ScoreU12DonorSite(string sequence, int position)
    {
        // Simplified U12 scoring - /ATATCC/ consensus
        string motif = position + 6 <= sequence.Length
            ? sequence.Substring(position, 6)
            : "";

        int matches = 0;
        string consensus = "AUAUCC";
        for (int i = 0; i < Math.Min(motif.Length, consensus.Length); i++)
        {
            if (motif[i] == consensus[i])
                matches++;
        }

        return matches / 6.0;
    }

    private static double ScoreU12AcceptorSite(string sequence, int position)
    {
        // Simplified U12 acceptor - /AC/ dinucleotide
        if (position >= 0 && position + 1 < sequence.Length)
        {
            if (sequence[position] == 'A' && sequence[position + 1] == 'C')
                return 0.6;
        }
        return 0;
    }

    #endregion

    #region Intron Prediction

    /// <summary>
    /// Predicts introns by pairing donor and acceptor sites.
    /// </summary>
    public static IEnumerable<Intron> PredictIntrons(
        string sequence,
        int minIntronLength = 60,
        int maxIntronLength = 100000,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        var donors = FindDonorSites(upper, minScore * 0.8, true).ToList();
        var acceptors = FindAcceptorSites(upper, minScore * 0.8, true).ToList();

        foreach (var donor in donors)
        {
            foreach (var acceptor in acceptors.Where(a => a.Position > donor.Position + minIntronLength))
            {
                int intronLength = acceptor.Position - donor.Position;
                if (intronLength > maxIntronLength)
                    continue;

                // Determine intron type
                var type = DetermineIntronType(donor, acceptor);

                // Find branch point
                int searchStart = Math.Max(0, acceptor.Position - 50);
                int searchEnd = acceptor.Position - 18;
                var branchPoints = FindBranchPoints(upper, searchStart, searchEnd, 0.4).ToList();
                var bestBranch = branchPoints.OrderByDescending(b => b.Score).FirstOrDefault();

                double combinedScore = (donor.Score + acceptor.Score + (bestBranch.Score > 0 ? bestBranch.Score : 0.3)) / 3;

                if (combinedScore >= minScore)
                {
                    yield return new Intron(
                        Start: donor.Position,
                        End: acceptor.Position,
                        Length: intronLength,
                        DonorSite: donor,
                        AcceptorSite: acceptor,
                        BranchPoint: bestBranch.Score > 0 ? bestBranch : null,
                        Sequence: upper.Substring(donor.Position, intronLength),
                        Type: type,
                        Score: combinedScore);
                }
            }
        }
    }

    private static IntronType DetermineIntronType(SpliceSite donor, SpliceSite acceptor)
    {
        if (donor.Type == SpliceSiteType.U12Donor || acceptor.Type == SpliceSiteType.U12Acceptor)
            return IntronType.U12;

        string donorMotif = donor.Motif.ToUpperInvariant();
        if (donorMotif.Length >= 2)
        {
            string dinuc = donorMotif.Substring(0, 2);
            if (dinuc == "GC")
                return IntronType.GcAg;
            if (dinuc == "GU" || dinuc == "GT")
                return IntronType.U2;
        }

        return IntronType.Unknown;
    }

    #endregion

    #region Gene Structure Prediction

    /// <summary>
    /// Predicts exon/intron structure of a gene.
    /// </summary>
    public static GeneStructure PredictGeneStructure(
        string sequence,
        int minExonLength = 30,
        int minIntronLength = 60,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new GeneStructure(
                new List<Exon>(),
                new List<Intron>(),
                "",
                0);
        }

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Get all introns
        var introns = PredictIntrons(upper, minIntronLength, 100000, minScore)
            .OrderByDescending(i => i.Score)
            .ToList();

        // Select non-overlapping introns greedily
        var selectedIntrons = SelectNonOverlappingIntrons(introns);

        // Derive exons from intron positions
        var exons = DeriveExons(upper, selectedIntrons, minExonLength);

        // Generate spliced sequence
        string splicedSequence = GenerateSplicedSequence(upper, selectedIntrons);

        double overallScore = selectedIntrons.Count > 0
            ? selectedIntrons.Average(i => i.Score)
            : 0;

        return new GeneStructure(
            Exons: exons,
            Introns: selectedIntrons,
            SplicedSequence: splicedSequence,
            OverallScore: overallScore);
    }

    private static List<Intron> SelectNonOverlappingIntrons(List<Intron> introns)
    {
        var selected = new List<Intron>();
        var used = new HashSet<int>();

        foreach (var intron in introns)
        {
            bool overlaps = false;
            for (int pos = intron.Start; pos <= intron.End; pos++)
            {
                if (used.Contains(pos))
                {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps)
            {
                selected.Add(intron);
                for (int pos = intron.Start; pos <= intron.End; pos++)
                    used.Add(pos);
            }
        }

        return selected.OrderBy(i => i.Start).ToList();
    }

    private static List<Exon> DeriveExons(string sequence, List<Intron> introns, int minExonLength)
    {
        var exons = new List<Exon>();

        if (introns.Count == 0)
        {
            // Single exon gene
            exons.Add(new Exon(
                Start: 0,
                End: sequence.Length - 1,
                Length: sequence.Length,
                Type: ExonType.Single,
                Phase: 0,
                Sequence: sequence));
            return exons;
        }

        int currentPos = 0;

        for (int i = 0; i < introns.Count; i++)
        {
            var intron = introns[i];
            int exonEnd = intron.Start - 1;

            if (exonEnd - currentPos + 1 >= minExonLength)
            {
                var exonType = i == 0 ? ExonType.Initial : ExonType.Internal;
                exons.Add(new Exon(
                    Start: currentPos,
                    End: exonEnd,
                    Length: exonEnd - currentPos + 1,
                    Type: exonType,
                    Phase: CalculatePhase(exons),
                    Sequence: sequence.Substring(currentPos, exonEnd - currentPos + 1)));
            }

            currentPos = intron.End + 1;
        }

        // Terminal exon
        if (sequence.Length - currentPos >= minExonLength)
        {
            exons.Add(new Exon(
                Start: currentPos,
                End: sequence.Length - 1,
                Length: sequence.Length - currentPos,
                Type: ExonType.Terminal,
                Phase: CalculatePhase(exons),
                Sequence: sequence.Substring(currentPos)));
        }

        return exons;
    }

    private static int CalculatePhase(List<Exon> previousExons)
    {
        int totalLength = previousExons.Sum(e => e.Length);
        return totalLength % 3;
    }

    private static string GenerateSplicedSequence(string sequence, List<Intron> introns)
    {
        if (introns.Count == 0)
            return sequence;

        var sb = new StringBuilder();
        int currentPos = 0;

        foreach (var intron in introns.OrderBy(i => i.Start))
        {
            if (intron.Start > currentPos)
            {
                sb.Append(sequence.Substring(currentPos, intron.Start - currentPos));
            }
            currentPos = intron.End + 1;
        }

        if (currentPos < sequence.Length)
        {
            sb.Append(sequence.Substring(currentPos));
        }

        return sb.ToString();
    }

    #endregion

    #region Alternative Splicing

    /// <summary>
    /// Detects potential alternative splicing patterns.
    /// </summary>
    public static IEnumerable<(string Type, int Position, string Description)> DetectAlternativeSplicing(
        string sequence,
        double minScore = 0.4)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        var donors = FindDonorSites(sequence, minScore).ToList();
        var acceptors = FindAcceptorSites(sequence, minScore).ToList();

        // Exon skipping: multiple introns with shared boundaries
        for (int i = 0; i < donors.Count; i++)
        {
            int validAcceptors = acceptors.Count(a => a.Position > donors[i].Position + 60);
            if (validAcceptors > 1)
            {
                yield return ("ExonSkipping", donors[i].Position,
                    $"Potential exon skipping from donor at {donors[i].Position}");
            }
        }

        // Alternative 5' splice sites (multiple donors before same acceptor)
        var donorGroups = donors.GroupBy(d => d.Position / 50).Where(g => g.Count() > 1);
        foreach (var group in donorGroups)
        {
            yield return ("Alt5SS", group.First().Position,
                $"Alternative 5' splice sites at positions {string.Join(", ", group.Select(d => d.Position))}");
        }

        // Alternative 3' splice sites (multiple acceptors after same donor)
        var acceptorGroups = acceptors.GroupBy(a => a.Position / 50).Where(g => g.Count() > 1);
        foreach (var group in acceptorGroups)
        {
            yield return ("Alt3SS", group.First().Position,
                $"Alternative 3' splice sites at positions {string.Join(", ", group.Select(a => a.Position))}");
        }
    }

    /// <summary>
    /// Finds retained intron candidates (introns that might be retained in some transcripts).
    /// </summary>
    public static IEnumerable<Intron> FindRetainedIntronCandidates(
        string sequence,
        double minScore = 0.5)
    {
        var introns = PredictIntrons(sequence, 60, 500, minScore).ToList();

        // Short introns with moderate scores are candidates for retention
        foreach (var intron in introns.Where(i => i.Length < 500 && i.Score < 0.8))
        {
            yield return intron;
        }
    }

    #endregion

    #region Utility Methods

    private static string GetMotifContext(string sequence, int position, int upstream, int downstream)
    {
        int start = Math.Max(0, position - upstream);
        int end = Math.Min(sequence.Length, position + downstream);
        return sequence.Substring(start, end - start);
    }

    private static double CalculateConfidence(double score, double minExpected, double maxExpected)
    {
        return Math.Max(0, Math.Min(1, (score - minExpected) / (maxExpected - minExpected)));
    }

    /// <summary>
    /// Calculates MaxEntScan-like score for splice sites.
    /// </summary>
    public static double CalculateMaxEntScore(string motif, SpliceSiteType type)
    {
        if (string.IsNullOrEmpty(motif))
            return 0;

        string upper = motif.ToUpperInvariant().Replace('T', 'U');

        if (type == SpliceSiteType.Donor)
        {
            // Score 9-mer around GT
            double score = 0;
            for (int i = 0; i < upper.Length && i < 9; i++)
            {
                int offset = i - 3;
                if (DonorPwm.TryGetValue(offset, out var weights))
                {
                    if (weights.TryGetValue(upper[i], out double w))
                        score += Math.Log2(w + 0.01);
                }
            }
            return score;
        }
        else if (type == SpliceSiteType.Acceptor)
        {
            double score = 0;
            for (int i = 0; i < upper.Length; i++)
            {
                int offset = i - 15;
                if (AcceptorPwm.TryGetValue(offset, out var weights))
                {
                    if (weights.TryGetValue(upper[i], out double w))
                        score += Math.Log2(w + 0.01);
                }
            }
            return score;
        }

        return 0;
    }

    /// <summary>
    /// Checks if a sequence position is within a coding region (simple heuristic).
    /// </summary>
    public static bool IsWithinCodingRegion(string sequence, int position, int frame = 0)
    {
        if (position < 0 || position >= sequence.Length)
            return false;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Look for start codon upstream
        int searchStart = Math.Max(0, position - 300);
        for (int i = searchStart; i <= position - 3; i++)
        {
            if (upper.Substring(i, 3) == "AUG")
            {
                // Check if position is in frame
                return (position - i) % 3 == frame;
            }
        }

        return false;
    }

    #endregion
}
