using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Predicts intrinsically disordered regions (IDRs) in proteins.
/// Uses hydropathy, charge distribution, and complexity-based methods.
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

    // Disorder propensity scale (based on DisProt statistics)
    private static readonly Dictionary<char, double> DisorderPropensity = new()
    {
        ['A'] = 0.06,
        ['R'] = 0.18,
        ['N'] = 0.01,
        ['D'] = 0.19,
        ['C'] = -0.20,
        ['Q'] = 0.32,
        ['E'] = 0.30,
        ['G'] = 0.16,
        ['H'] = -0.01,
        ['I'] = -0.49,
        ['L'] = -0.34,
        ['K'] = 0.27,
        ['M'] = -0.11,
        ['F'] = -0.42,
        ['P'] = 0.41,
        ['S'] = 0.15,
        ['T'] = 0.04,
        ['W'] = -0.49,
        ['Y'] = -0.31,
        ['V'] = -0.38
    };

    // Charge at pH 7
    private static readonly Dictionary<char, double> Charge = new()
    {
        ['A'] = 0,
        ['R'] = 1,
        ['N'] = 0,
        ['D'] = -1,
        ['C'] = 0,
        ['Q'] = 0,
        ['E'] = -1,
        ['G'] = 0,
        ['H'] = 0.1,
        ['I'] = 0,
        ['L'] = 0,
        ['K'] = 1,
        ['M'] = 0,
        ['F'] = 0,
        ['P'] = 0,
        ['S'] = 0,
        ['T'] = 0,
        ['W'] = 0,
        ['Y'] = 0,
        ['V'] = 0
    };

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
        double disorderThreshold = 0.5,
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
    /// Calculates combined disorder score for a window.
    /// </summary>
    private static double CalculateDisorderScore(string window)
    {
        if (window.Length == 0)
            return 0;

        // Component 1: Disorder propensity
        double propensityScore = 0;
        int validCount = 0;
        foreach (char c in window)
        {
            if (DisorderPropensity.TryGetValue(c, out double prop))
            {
                propensityScore += prop;
                validCount++;
            }
        }
        propensityScore = validCount > 0 ? propensityScore / validCount : 0;

        // Component 2: Low hydropathy (disordered regions are hydrophilic)
        double hydropathyScore = 0;
        validCount = 0;
        foreach (char c in window)
        {
            if (Hydropathy.TryGetValue(c, out double hydro))
            {
                hydropathyScore += hydro;
                validCount++;
            }
        }
        double meanHydropathy = validCount > 0 ? hydropathyScore / validCount : 0;
        // Convert to disorder score (negative hydropathy = more disorder)
        double hydroDisorder = (-meanHydropathy + 4.5) / 9.0; // Normalize to 0-1

        // Component 3: Net charge (high charge = more disorder)
        double chargeScore = CalculateNetCharge(window);
        double chargeDisorder = Math.Min(1.0, Math.Abs(chargeScore) / window.Length * 5);

        // Component 4: Low complexity (repetitive = more disorder)
        double complexity = CalculateSequenceComplexity(window);
        double complexityDisorder = 1.0 - complexity;

        // Component 5: Proline/Glycine content (disorder promoting)
        int proGly = window.Count(c => c == 'P' || c == 'G');
        double proGlyDisorder = proGly / (double)window.Length;

        // Combine scores with weights
        double combinedScore =
            0.35 * NormalizeScore(propensityScore, -0.5, 0.5) +
            0.25 * hydroDisorder +
            0.15 * chargeDisorder +
            0.15 * complexityDisorder +
            0.10 * proGlyDisorder;

        return Math.Max(0, Math.Min(1, combinedScore));
    }

    /// <summary>
    /// Normalizes score to 0-1 range.
    /// </summary>
    private static double NormalizeScore(double value, double min, double max)
    {
        return (value - min) / (max - min);
    }

    /// <summary>
    /// Calculates net charge of sequence.
    /// </summary>
    private static double CalculateNetCharge(string sequence)
    {
        double charge = 0;
        foreach (char c in sequence)
        {
            if (Charge.TryGetValue(c, out double q))
                charge += q;
        }
        return charge;
    }

    /// <summary>
    /// Calculates sequence complexity (0-1).
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
        // Higher score and longer regions = more confidence
        double scoreConfidence = (meanScore - 0.5) * 2; // 0.5-1.0 -> 0-1
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
        bool moderateDisorder = meanDisorder > -0.2 && meanDisorder < 0.3;
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
    /// Checks if amino acid is disorder-promoting.
    /// </summary>
    public static bool IsDisorderPromoting(char aminoAcid)
    {
        return GetDisorderPropensity(aminoAcid) > 0;
    }

    /// <summary>
    /// Gets list of disorder-promoting amino acids.
    /// </summary>
    public static IReadOnlyList<char> DisorderPromotingAminoAcids =>
        DisorderPropensity.Where(kv => kv.Value > 0).Select(kv => kv.Key).ToList();

    /// <summary>
    /// Gets list of order-promoting amino acids.
    /// </summary>
    public static IReadOnlyList<char> OrderPromotingAminoAcids =>
        DisorderPropensity.Where(kv => kv.Value < 0).Select(kv => kv.Key).ToList();

    #endregion
}
