using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Predicts intrinsically disordered regions (IDRs) in proteins.
/// Scoring: TOP-IDP scale averaged over a sliding window — Campen et al. (2008) PMC2676888.
/// Classification: Dunker et al. (2001) PMID 11381529.
/// </summary>
/// <remarks>
/// <para><b>Scope and limitations.</b>
/// This is a single-feature heuristic toolkit, not a competitive disorder predictor.
/// All components rely exclusively on amino-acid composition (TOP-IDP propensity,
/// Kyte-Doolittle hydropathy, Shannon entropy). No evolutionary information
/// (PSI-BLAST profiles), predicted secondary structure, or machine-learning models
/// are used.</para>
///
/// <para><b>Expected accuracy.</b>
/// Single-scale predictors such as TOP-IDP achieve AUC ≈ 0.65-0.72 on DisProt
/// benchmarks (Campen et al. 2008, PMC2676888). For comparison:
/// IUPred2A (Mészáros et al. 2018, NAR 46:W315-W320, PMID 29860432) — AUC ≈ 0.75-0.80;
/// flDPnn (Hu et al. 2021, NAR 49:e80, PMID 34023900) — AUC ≈ 0.85-0.90;
/// Critical Assessment of Intrinsic Disorder (CAID) benchmark
/// (Necci et al. 2021, Nat Methods 18:472-481, PMID 33875884) provides
/// community-wide accuracy comparisons.</para>
///
/// <para><b>Low complexity (SEG).</b>
/// The SEG implementation is an approximation (trigger + extend two-pass)
/// of the original recursive algorithm by Wootton &amp; Federhen (1993, 1996).
/// Results are similar for typical LC regions but may diverge at threshold boundaries.</para>
///
/// <para><b>MoRF prediction.</b>
/// The MoRF annotator is a hydropathy-enrichment heuristic inspired by
/// Mohan et al. (2006, J Mol Biol 362:1043-1059, PMID 16949612),
/// not a trained predictor. Dedicated tools such as MoRFchibi
/// (Malhis et al. 2016, Bioinformatics 32:1906-1913, PMID 27153720)
/// and ANCHOR2 (Mészáros et al. 2009, J Mol Biol 392:760-773, PMID 19654615)
/// achieve substantially higher accuracy.</para>
///
/// <para><b>Recommended use.</b>
/// Suitable for educational purposes, demonstration, and coarse first-pass screening.
/// For publication-grade or clinical annotation, use IUPred2A, MobiDB-lite
/// (Necci et al. 2017, Bioinformatics 33:1765-1767, PMID 28130230), or an
/// ML model trained and validated on DisProt (Hatos et al. 2020, NAR 48:D269-D276,
/// PMID 31713636).</para>
/// </remarks>
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

    // Pre-sorted cached lists for public API.
    private static readonly IReadOnlyList<char> CachedDisorderPromoting =
        new List<char>(DisorderPromotingSet.OrderBy(c => c));
    private static readonly IReadOnlyList<char> CachedOrderPromoting =
        new List<char>(OrderPromotingSet.OrderBy(c => c));
    private static readonly IReadOnlyList<char> CachedAmbiguous =
        new List<char>(AmbiguousSet.OrderBy(c => c));

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
    /// Shannon entropy in bits for a subsequence.
    /// H = −Σ p_i · log₂(p_i), summing over distinct residue types.
    /// SEG uses this formula directly (not normalized).
    /// Source: Shannon (1948); used by Wootton &amp; Federhen (1993, 1996).
    /// </summary>
    private static double CalculateShannonEntropy(string sequence, int start, int length)
    {
        Span<int> counts = stackalloc int[26];
        for (int i = start; i < start + length; i++)
        {
            int idx = sequence[i] - 'A';
            if ((uint)idx < 26)
                counts[idx]++;
        }

        double entropy = 0;
        for (int i = 0; i < 26; i++)
        {
            if (counts[i] > 0)
            {
                double p = counts[i] / (double)length;
                entropy -= p * Math.Log2(p);
            }
        }

        return entropy;
    }

    /// <summary>
    /// Mean Kyte-Doolittle hydropathy for a subsequence.
    /// Source: Kyte &amp; Doolittle (1982) J Mol Biol 157:105-132.
    /// </summary>
    private static double MeanHydropathy(string sequence, int start, int length)
    {
        double sum = 0;
        int count = 0;
        for (int i = start; i < start + length; i++)
        {
            if (Hydropathy.TryGetValue(sequence[i], out double h))
            {
                sum += h;
                count++;
            }
        }
        return count > 0 ? sum / count : 0;
    }

    /// <summary>
    /// Merges overlapping or adjacent (Start, End) segments.
    /// </summary>
    private static List<(int Start, int End)> MergeSegments(List<(int Start, int End)> segments)
    {
        if (segments.Count == 0)
            return segments;

        var sorted = segments.OrderBy(s => s.Start).ThenBy(s => s.End).ToList();
        var merged = new List<(int Start, int End)> { sorted[0] };

        for (int i = 1; i < sorted.Count; i++)
        {
            var last = merged[^1];
            var current = sorted[i];
            if (current.Start <= last.End + 1)
                merged[^1] = (last.Start, Math.Max(last.End, current.End));
            else
                merged.Add(current);
        }

        return merged;
    }

    /// <summary>
    /// Classifies a low complexity region by its dominant amino acid composition.
    /// </summary>
    private static string ClassifyLowComplexityType(string sequence, int start, int length)
    {
        var counts = new Dictionary<char, int>();
        for (int i = start; i < start + length; i++)
        {
            char c = sequence[i];
            counts[c] = counts.GetValueOrDefault(c) + 1;
        }

        var ranked = counts.OrderByDescending(kv => kv.Value).ToList();
        double topFraction = ranked[0].Value / (double)length;

        if (topFraction > 0.5)
            return $"{ranked[0].Key}-rich";

        if (ranked.Count >= 2)
            return $"{ranked[0].Key}/{ranked[1].Key}-rich";

        return $"{ranked[0].Key}-rich";
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
                        CalculateConfidence(meanScore),
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
                    CalculateConfidence(meanScore),
                    regionType);
            }
        }
    }

    /// <summary>
    /// Classifies the type of disordered region by compositional bias.
    /// Priority order (internal design decision, no single authoritative source):
    /// Proline → Acidic (E+D) → Basic (K+R) → Ser/Thr-rich (S+T) → Long IDR (&gt;30) → Standard IDR.
    /// Enrichment threshold: fraction &gt; 0.25 = 5× random single-AA frequency (1/20).
    /// AA groups follow standard biochemical classifications.
    /// Subtype names inspired by van der Lee et al. (2014) Chem Rev 114:6589–6631.
    /// Length threshold (&gt;30) from Ward et al. (2004) J Mol Biol 337:635–645.
    /// </summary>
    private static string ClassifyDisorderedRegion(List<ResiduePrediction> region)
    {
        var sequence = new string(region.Select(r => r.Residue).ToArray());

        int proCount = sequence.Count(c => c == 'P');
        int gluAspCount = sequence.Count(c => c == 'E' || c == 'D');
        int serThrCount = sequence.Count(c => c == 'S' || c == 'T');
        int lysArgCount = sequence.Count(c => c == 'K' || c == 'R');

        double proFraction = proCount / (double)sequence.Length;
        double acidFraction = gluAspCount / (double)sequence.Length;
        double serThrFraction = serThrCount / (double)sequence.Length;
        double basicFraction = lysArgCount / (double)sequence.Length;

        // Enrichment threshold: 0.25 = 5× random single-AA frequency (1/20 = 0.05).
        // Internal heuristic — no single published source defines this value.
        const double enrichmentThreshold = 0.25;

        if (proFraction > enrichmentThreshold)
            return "Proline-rich";
        if (acidFraction > enrichmentThreshold)
            return "Acidic";
        if (basicFraction > enrichmentThreshold)
            return "Basic";
        if (serThrFraction > enrichmentThreshold)
            return "Ser/Thr-rich";

        // Ward et al. (2004): >30 residues defines a "substantial" disordered segment.
        // Van der Lee et al. (2014): "44% of human protein-coding genes contain
        // substantial disordered segments of >30 amino acids in length."
        if (sequence.Length > 30)
            return "Long IDR";

        return "Standard IDR";
    }

    /// <summary>
    /// Calculates prediction confidence as the normalized distance from the
    /// TOP-IDP decision boundary (cutoff 0.542 from Campen et al. 2008).
    /// Formula: confidence = (meanScore − cutoff) / (1.0 − cutoff), clamped to [0, 1].
    /// The formula itself is an internal design decision — Campen defines only the cutoff.
    /// </summary>
    private static double CalculateConfidence(double meanScore)
    {
        // Normalized distance from the TOP-IDP cutoff (0.542).
        // Score 0.542 → confidence 0.0; score 1.0 → confidence 1.0.
        return Math.Max(0, Math.Min(1, (meanScore - TopIdpCutoff) / (1.0 - TopIdpCutoff)));
    }

    #endregion

    #region Specialized Predictions

    /// <summary>
    /// Predicts low complexity regions using the SEG algorithm.
    /// Two-pass approach: K1 trigger scan on small window, K2 extension.
    /// Source: Wootton &amp; Federhen (1993) Computers &amp; Chemistry 17(2):149-163.
    ///         Wootton &amp; Federhen (1996) Methods Enzymol 266:554-571.
    /// Default parameters (standard SEG): triggerWindow=12, K1=2.2 bits, K2=2.5 bits.
    /// </summary>
    public static IEnumerable<(int Start, int End, string Type)> PredictLowComplexityRegions(
        string sequence,
        int triggerWindow = 12,
        double triggerThreshold = 2.2,
        double extensionThreshold = 2.5,
        int minLength = 1)
    {
        sequence = sequence.ToUpperInvariant();
        if (sequence.Length < triggerWindow)
            yield break;

        // Phase 1: K1 trigger — slide window, mark positions with entropy ≤ K1
        var triggered = new bool[sequence.Length];
        for (int i = 0; i <= sequence.Length - triggerWindow; i++)
        {
            double entropy = CalculateShannonEntropy(sequence, i, triggerWindow);
            if (entropy <= triggerThreshold)
            {
                for (int j = i; j < i + triggerWindow; j++)
                    triggered[j] = true;
            }
        }

        // Phase 2: Collect contiguous triggered spans
        var spans = new List<(int Start, int End)>();
        int? spanStart = null;
        for (int i = 0; i < sequence.Length; i++)
        {
            if (triggered[i])
            {
                spanStart ??= i;
            }
            else if (spanStart.HasValue)
            {
                spans.Add((spanStart.Value, i - 1));
                spanStart = null;
            }
        }
        if (spanStart.HasValue)
            spans.Add((spanStart.Value, sequence.Length - 1));

        // Phase 3: K2 extension — extend each span while entropy ≤ K2
        var extended = new List<(int Start, int End)>();
        foreach (var (start, end) in spans)
        {
            int s = start;
            int e = end;

            while (s > 0 && CalculateShannonEntropy(sequence, s - 1, e - (s - 1) + 1) <= extensionThreshold)
                s--;

            while (e < sequence.Length - 1 && CalculateShannonEntropy(sequence, s, (e + 1) - s + 1) <= extensionThreshold)
                e++;

            extended.Add((s, e));
        }

        // Phase 4: Merge overlapping or adjacent segments
        var merged = MergeSegments(extended);

        // Phase 5: Filter by minimum length and classify
        foreach (var (start, end) in merged)
        {
            int length = end - start + 1;
            if (length >= minLength)
            {
                string type = ClassifyLowComplexityType(sequence, start, length);
                yield return (start, end, type);
            }
        }
    }

    /// <summary>
    /// Predicts Molecular Recognition Features (MoRFs) within disordered regions.
    /// MoRFs are short segments within IDRs that undergo disorder-to-order transition
    /// upon binding a partner protein. Identified here by elevated Kyte-Doolittle
    /// hydropathy relative to the surrounding IDR context.
    ///
    /// Definition: Mohan et al. (2006) J Mol Biol 362:1043-1059.
    /// Hydropathy scale: Kyte &amp; Doolittle (1982) J Mol Biol 157:105-132.
    /// Note: this is a heuristic annotation, not a machine-learning predictor.
    /// </summary>
    public static IEnumerable<(int Start, int End, double Score)> PredictMoRFs(
        string sequence,
        int minLength = 10,
        int maxLength = 25)
    {
        sequence = sequence.ToUpperInvariant();
        var prediction = PredictDisorder(sequence);

        var candidates = new List<(int Start, int End, double Score)>();

        foreach (var region in prediction.DisorderedRegions)
        {
            int regionLength = region.End - region.Start + 1;
            if (regionLength < minLength)
                continue;

            // Mean hydropathy of the entire IDR — context baseline
            double idrHydropathy = MeanHydropathy(sequence, region.Start, regionLength);

            // Scan all windows within the IDR
            int effectiveMaxLen = Math.Min(maxLength, regionLength);
            for (int len = minLength; len <= effectiveMaxLen; len++)
            {
                for (int s = region.Start; s + len - 1 <= region.End; s++)
                {
                    double windowHydropathy = MeanHydropathy(sequence, s, len);

                    // MoRF criterion: window more hydrophobic than IDR context
                    // Epsilon 0.01 KD units filters floating-point noise
                    if (windowHydropathy > idrHydropathy + 0.01)
                    {
                        double diff = windowHydropathy - idrHydropathy;
                        candidates.Add((s, s + len - 1, diff));
                    }
                }
            }
        }

        // Greedy non-overlapping selection: best score first
        var sorted = candidates.OrderByDescending(c => c.Score).ToList();
        var selected = new List<(int Start, int End, double Score)>();
        foreach (var cand in sorted)
        {
            bool overlaps = selected.Any(sel => cand.Start <= sel.End && cand.End >= sel.Start);
            if (!overlaps)
                selected.Add(cand);
        }

        // Normalize scores to [0, 1]: difference / 3.0, clamped
        // 3.0 KD units ≈ maximum plausible hydrophobic island vs IDR context
        selected = selected
            .Select(s => (s.Start, s.End, Math.Min(1.0, s.Score / 3.0)))
            .ToList();

        return selected.OrderBy(s => s.Start);
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
    public static IReadOnlyList<char> DisorderPromotingAminoAcids => CachedDisorderPromoting;

    /// <summary>
    /// Gets list of order-promoting amino acids — Dunker et al. (2001).
    /// {C, F, I, L, N, V, W, Y}
    /// </summary>
    public static IReadOnlyList<char> OrderPromotingAminoAcids => CachedOrderPromoting;

    /// <summary>
    /// Gets list of ambiguous amino acids — Dunker et al. (2001).
    /// {D, H, M, T} — neither clearly disorder-promoting nor order-promoting.
    /// </summary>
    public static IReadOnlyList<char> AmbiguousAminoAcids => CachedAmbiguous;

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
