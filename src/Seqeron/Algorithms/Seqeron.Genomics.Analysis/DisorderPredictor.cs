using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Predicts intrinsically disordered regions (IDRs) in proteins.
/// Scoring: TOP-IDP scale averaged over a sliding window — Campen et al. (2008) PMC2676888.
/// Classification: Dunker et al. (2001) PMID 11381529.
/// </summary>
public static class DisorderPredictor
{
    #region Constants

    // Kyte-Doolittle hydropathy scale
    private static readonly Dictionary<char, double> Hydropathy = new()
    {
        ['A'] = 1.8,
        ['R'] = -4.5,
        ['N'] = -3.5,
        ['D'] = -3.5,
        ['C'] = 2.5,
        ['Q'] = -3.5,
        ['E'] = -3.5,
        ['G'] = -0.4,
        ['H'] = -3.2,
        ['I'] = 4.5,
        ['L'] = 3.8,
        ['K'] = -3.9,
        ['M'] = 1.9,
        ['F'] = 2.8,
        ['P'] = -1.6,
        ['S'] = -0.8,
        ['T'] = -0.7,
        ['W'] = -0.9,
        ['Y'] = -1.3,
        ['V'] = 4.2
    };

    // TOP-IDP disorder propensity scale
    // Source: Campen et al. (2008) "TOP-IDP-Scale: A New Amino Acid Scale Measuring
    //   Propensity for Intrinsic Disorder" Protein Pept Lett 15(9):956-963.
    //   PMC2676888, PMID 18991772, Table 2.
    // Ranking (order→disorder): W,F,Y,I,M,L,V,N,C,T,A,G,R,D,H,Q,K,S,E,P
    private static readonly Dictionary<char, double> DisorderPropensity = new()
    {
        ['A'] = 0.060,
        ['R'] = 0.180,
        ['N'] = 0.007,
        ['D'] = 0.192,
        ['C'] = 0.020,
        ['Q'] = 0.318,
        ['E'] = 0.736,
        ['G'] = 0.166,
        ['H'] = 0.303,
        ['I'] = -0.486,
        ['L'] = -0.326,
        ['K'] = 0.586,
        ['M'] = -0.397,
        ['F'] = -0.697,
        ['P'] = 0.987,
        ['S'] = 0.341,
        ['T'] = 0.059,
        ['W'] = -0.884,
        ['Y'] = -0.510,
        ['V'] = -0.121
    };

    // Disorder-promoting amino acids — Dunker et al. (2001) J Mol Graph Model 19:26-59.
    // PMID 11381529. Confirmed by Campen et al. (2008) TOP-IDP ranking.
    private static readonly HashSet<char> DisorderPromotingSet =
        new() { 'A', 'R', 'G', 'Q', 'S', 'P', 'E', 'K' };

    // Order-promoting amino acids — Dunker et al. (2001).
    private static readonly HashSet<char> OrderPromotingSet =
        new() { 'W', 'C', 'F', 'I', 'Y', 'V', 'L', 'N' };

    // Ambiguous amino acids — Dunker et al. (2001).
    // Neither clearly disorder-promoting nor order-promoting.
    private static readonly HashSet<char> AmbiguousSet =
        new() { 'H', 'M', 'T', 'D' };

    // TOP-IDP normalization constants — Campen et al. (2008) Table 2.
    private const double TopIdpMin = -0.884;  // W (most order-promoting)
    private const double TopIdpMax = 0.987;   // P (most disorder-promoting)
    private const double TopIdpRange = TopIdpMax - TopIdpMin; // 1.871

    /// <summary>
    /// TOP-IDP prediction cutoff — Campen et al. (2008) PMC2676888.
    /// Based on maximum-likelihood methods.
    /// </summary>
    private const double TopIdpCutoff = 0.542;

    #endregion

    #region Records

    /// <summary>
    /// Disordered region prediction result.
    /// </summary>
    public readonly record struct DisorderedRegion(
        int Start,
        int End,
        double MeanScore,
        double Confidence,
        string RegionType);

    /// <summary>
    /// Per-residue disorder prediction.
    /// </summary>
    public readonly record struct ResiduePrediction(
        int Position,
        char Residue,
        double DisorderScore,
        bool IsDisordered);

    /// <summary>
    /// Full prediction result.
    /// </summary>
    public readonly record struct DisorderPredictionResult(
        string Sequence,
        IReadOnlyList<ResiduePrediction> ResiduePredictions,
        IReadOnlyList<DisorderedRegion> DisorderedRegions,
        double OverallDisorderContent,
        double MeanDisorderScore);

    #endregion

    #region Main Prediction Methods

    /// <summary>
    /// Predicts intrinsically disordered regions in a protein sequence.
    /// </summary>
    public static DisorderPredictionResult PredictDisorder(
        string sequence,
        int windowSize = 21,
        double disorderThreshold = TopIdpCutoff,
        int minRegionLength = 5)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new DisorderPredictionResult(
                "", new List<ResiduePrediction>(), new List<DisorderedRegion>(), 0, 0);
        }

        sequence = sequence.ToUpperInvariant();

        // Calculate per-residue disorder scores
        var residuePredictions = CalculatePerResidueScores(sequence, windowSize, disorderThreshold);

        // Identify disordered regions
        var regions = IdentifyDisorderedRegions(
            residuePredictions, disorderThreshold, minRegionLength).ToList();

        // Calculate overall statistics
        int disorderedCount = residuePredictions.Count(r => r.IsDisordered);
        double disorderContent = disorderedCount / (double)sequence.Length;
        double meanScore = residuePredictions.Average(r => r.DisorderScore);

        return new DisorderPredictionResult(
            sequence,
            residuePredictions,
            regions,
            disorderContent,
            meanScore);
    }

    /// <summary>
    /// Calculates per-residue disorder scores.
    /// </summary>
    private static List<ResiduePrediction> CalculatePerResidueScores(
        string sequence,
        int windowSize,
        double threshold)
    {
        var predictions = new List<ResiduePrediction>();
        int halfWindow = windowSize / 2;

        for (int i = 0; i < sequence.Length; i++)
        {
            int start = Math.Max(0, i - halfWindow);
            int end = Math.Min(sequence.Length, i + halfWindow + 1);
            string window = sequence[start..end];

            double score = CalculateDisorderScore(window);
            bool isDisordered = score >= threshold;

            predictions.Add(new ResiduePrediction(i, sequence[i], score, isDisordered));
        }

        return predictions;
    }

    /// <summary>
    /// Calculates disorder score for a window using the TOP-IDP prediction method.
    /// Score = average of normalized TOP-IDP values over the window.
    /// Source: Campen et al. (2008) PMC2676888, Section "TOP-IDP web-based prediction service".
    /// </summary>
    private static double CalculateDisorderScore(string window)
    {
        if (window.Length == 0)
            return 0;

        double sum = 0;
        int count = 0;
        foreach (char c in window)
        {
            if (DisorderPropensity.TryGetValue(c, out double prop))
            {
                sum += (prop - TopIdpMin) / TopIdpRange;
                count++;
            }
        }

        return count > 0 ? sum / count : 0;
    }

    /// <summary>
    /// Calculates sequence complexity (Shannon entropy, normalized to 0-1).
    /// Source: Shannon (1948).
    /// </summary>
    private static double CalculateSequenceComplexity(string sequence)
    {
        if (sequence.Length < 2)
            return 0;

        var counts = new Dictionary<char, int>();
        foreach (char c in sequence)
            counts[c] = counts.GetValueOrDefault(c) + 1;

        // Shannon entropy
        double entropy = 0;
        foreach (var count in counts.Values)
        {
            double p = count / (double)sequence.Length;
            if (p > 0)
                entropy -= p * Math.Log2(p);
        }

        // Normalize (max entropy for 20 amino acids)
        double maxEntropy = Math.Log2(Math.Min(20, sequence.Length));
        return maxEntropy > 0 ? entropy / maxEntropy : 0;
    }

    #endregion

    #region Region Identification

    /// <summary>
    /// Identifies contiguous disordered regions.
    /// </summary>
    private static IEnumerable<DisorderedRegion> IdentifyDisorderedRegions(
        List<ResiduePrediction> predictions,
        double threshold,
        int minLength)
    {
        int? regionStart = null;
        double scoreSum = 0;

        for (int i = 0; i < predictions.Count; i++)
        {
            var pred = predictions[i];

            if (pred.IsDisordered)
            {
                if (!regionStart.HasValue)
                {
                    regionStart = i;
                    scoreSum = 0;
                }
                scoreSum += pred.DisorderScore;
            }
            else if (regionStart.HasValue)
            {
                int length = i - regionStart.Value;
                if (length >= minLength)
                {
                    double meanScore = scoreSum / length;
                    string regionType = ClassifyDisorderedRegion(
                        predictions.Skip(regionStart.Value).Take(length).ToList());

                    yield return new DisorderedRegion(
                        regionStart.Value,
                        i - 1,
                        meanScore,
                        CalculateConfidence(meanScore, length),
                        regionType);
                }
                regionStart = null;
            }
        }

        // Handle region at end
        if (regionStart.HasValue)
        {
            int length = predictions.Count - regionStart.Value;
            if (length >= minLength)
            {
                double meanScore = scoreSum / length;
                string regionType = ClassifyDisorderedRegion(
                    predictions.Skip(regionStart.Value).Take(length).ToList());

                yield return new DisorderedRegion(
                    regionStart.Value,
                    predictions.Count - 1,
                    meanScore,
                    CalculateConfidence(meanScore, length),
                    regionType);
            }
        }
    }

    /// <summary>
    /// Classifies the type of disordered region.
    /// </summary>
    private static string ClassifyDisorderedRegion(List<ResiduePrediction> region)
    {
        var sequence = new string(region.Select(r => r.Residue).ToArray());

        // Check for specific patterns
        int proCount = sequence.Count(c => c == 'P');
        int gluAspCount = sequence.Count(c => c == 'E' || c == 'D');
        int serThrCount = sequence.Count(c => c == 'S' || c == 'T');
        int lysArgCount = sequence.Count(c => c == 'K' || c == 'R');

        double proFraction = proCount / (double)sequence.Length;
        double acidFraction = gluAspCount / (double)sequence.Length;
        double serThrFraction = serThrCount / (double)sequence.Length;
        double basicFraction = lysArgCount / (double)sequence.Length;

        if (proFraction > 0.25)
            return "Proline-rich";
        if (acidFraction > 0.25)
            return "Acidic";
        if (basicFraction > 0.25)
            return "Basic";
        if (serThrFraction > 0.25)
            return "Ser/Thr-rich";
        if (sequence.Length > 30)
            return "Long IDR";

        return "Standard IDR";
    }

    /// <summary>
    /// Calculates prediction confidence.
    /// </summary>
    private static double CalculateConfidence(double meanScore, int length)
    {
        // Distance from TOP-IDP cutoff, normalized to [0,1]
        double scoreConfidence = (meanScore - TopIdpCutoff) / (1.0 - TopIdpCutoff);
        double lengthConfidence = Math.Min(1.0, length / 20.0);

        return Math.Max(0, Math.Min(1, (scoreConfidence + lengthConfidence) / 2));
    }

    #endregion

    #region Specialized Predictions

    /// <summary>
    /// Predicts binding sites within disordered regions (MoRFs - Molecular Recognition Features).
    /// </summary>
    public static IEnumerable<(int Start, int End, double Score)> PredictMoRFs(
        string sequence,
        int minLength = 10,
        int maxLength = 30)
    {
        sequence = sequence.ToUpperInvariant();

        // MoRFs tend to have:
        // - Moderate disorder (transition regions)
        // - Hydrophobic residues that can form structure upon binding
        // - Located within or adjacent to disordered regions

        var prediction = PredictDisorder(sequence);

        foreach (var region in prediction.DisorderedRegions)
        {
            // Scan for potential MoRFs at region boundaries
            for (int len = minLength; len <= maxLength && len <= region.End - region.Start + 1; len++)
            {
                // Check start of region
                if (region.Start > 0)
                {
                    int start = region.Start;
                    string window = sequence.Substring(start, Math.Min(len, sequence.Length - start));
                    double score = CalculateMoRFScore(window);

                    if (score > 0.5)
                        yield return (start, start + window.Length - 1, score);
                }

                // Check end of region
                if (region.End < sequence.Length - 1)
                {
                    int end = region.End;
                    int start = Math.Max(0, end - len + 1);
                    string window = sequence.Substring(start, end - start + 1);
                    double score = CalculateMoRFScore(window);

                    if (score > 0.5)
                        yield return (start, end, score);
                }
            }
        }
    }

    /// <summary>
    /// Calculates MoRF score.
    /// </summary>
    private static double CalculateMoRFScore(string window)
    {
        // MoRFs have moderate hydrophobicity and disorder
        double hydro = 0;
        double disorder = 0;
        int count = 0;

        foreach (char c in window)
        {
            if (Hydropathy.TryGetValue(c, out double h))
            {
                hydro += h;
                count++;
            }
            if (DisorderPropensity.TryGetValue(c, out double d))
            {
                disorder += d;
            }
        }

        if (count == 0) return 0;

        double meanHydro = hydro / count;
        double meanDisorder = disorder / count;

        // MoRFs: moderate disorder, some hydrophobic content
        bool moderateDisorder = meanDisorder > -0.4 && meanDisorder < 0.5;
        bool hasHydrophobic = meanHydro > -1.0 && meanHydro < 1.0;

        if (moderateDisorder && hasHydrophobic)
            return 0.5 + 0.5 * (1 - Math.Abs(meanDisorder));

        return 0.3;
    }

    /// <summary>
    /// Predicts low complexity regions.
    /// </summary>
    public static IEnumerable<(int Start, int End, string Type)> PredictLowComplexityRegions(
        string sequence,
        int windowSize = 12,
        double complexityThreshold = 0.3,
        int minLength = 10)
    {
        sequence = sequence.ToUpperInvariant();

        int? regionStart = null;
        string? dominantAA = null;

        for (int i = 0; i <= sequence.Length - windowSize; i++)
        {
            string window = sequence.Substring(i, windowSize);
            double complexity = CalculateSequenceComplexity(window);

            if (complexity < complexityThreshold)
            {
                if (!regionStart.HasValue)
                {
                    regionStart = i;
                    // Find dominant amino acid
                    var counts = window.GroupBy(c => c)
                        .OrderByDescending(g => g.Count())
                        .First();
                    dominantAA = $"{counts.Key}-rich";
                }
            }
            else if (regionStart.HasValue)
            {
                int length = i - regionStart.Value + windowSize - 1;
                if (length >= minLength)
                {
                    yield return (regionStart.Value, regionStart.Value + length - 1, dominantAA!);
                }
                regionStart = null;
            }
        }

        if (regionStart.HasValue)
        {
            int length = sequence.Length - regionStart.Value;
            if (length >= minLength)
            {
                yield return (regionStart.Value, sequence.Length - 1, dominantAA!);
            }
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Gets disorder propensity for an amino acid.
    /// </summary>
    public static double GetDisorderPropensity(char aminoAcid)
    {
        return DisorderPropensity.GetValueOrDefault(char.ToUpperInvariant(aminoAcid), 0);
    }

    /// <summary>
    /// Checks if amino acid is disorder-promoting per Dunker et al. (2001).
    /// Returns true for {A, R, G, Q, S, P, E, K}.
    /// Returns false for order-promoting {W, C, F, I, Y, V, L, N}
    /// and ambiguous {H, M, T, D} residues.
    /// </summary>
    public static bool IsDisorderPromoting(char aminoAcid)
    {
        return DisorderPromotingSet.Contains(char.ToUpperInvariant(aminoAcid));
    }

    /// <summary>
    /// Gets list of disorder-promoting amino acids — Dunker et al. (2001).
    /// {A, E, G, K, P, Q, R, S}
    /// </summary>
    public static IReadOnlyList<char> DisorderPromotingAminoAcids =>
        DisorderPromotingSet.OrderBy(c => c).ToList();

    /// <summary>
    /// Gets list of order-promoting amino acids — Dunker et al. (2001).
    /// {C, F, I, L, N, V, W, Y}
    /// </summary>
    public static IReadOnlyList<char> OrderPromotingAminoAcids =>
        OrderPromotingSet.OrderBy(c => c).ToList();

    /// <summary>
    /// Gets list of ambiguous amino acids — Dunker et al. (2001).
    /// {D, H, M, T} — neither clearly disorder-promoting nor order-promoting.
    /// </summary>
    public static IReadOnlyList<char> AmbiguousAminoAcids =>
        AmbiguousSet.OrderBy(c => c).ToList();

    /// <summary>
    /// Calculates mean Kyte-Doolittle hydropathy for a sequence.
    /// Source: Kyte &amp; Doolittle (1982) J Mol Biol 157:105-132.
    /// </summary>
    public static double CalculateHydropathy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        sequence = sequence.ToUpperInvariant();
        double sum = 0;
        int count = 0;
        foreach (char c in sequence)
        {
            if (Hydropathy.TryGetValue(c, out double value))
            {
                sum += value;
                count++;
            }
        }
        return count > 0 ? sum / count : 0;
    }

    #endregion
}
