namespace Seqeron.Genomics.Analysis;

// B09 A3-11: sourmash ANI estimates with confidence intervals (Hera, Pierce-Ward, Koslicki 2023, "Estimating the average nucleotide identity of genomes with FracMinHash";
// sourmash 4.9.4 sourmash/distance_utils.py, equations of github.com/KoslickiLab/mutation-rate-ci-calculator), ported statement by statement from the installed Python package.
public static partial class GenomicAnalyzer
{
    #region ANI with confidence intervals (sourmash distance_utils)

    /// <summary>
    /// sourmash <c>containment_to_distance</c>: the mutation-rate distance and (optionally) its confidence interval from a FracMinHash containment.
    /// </summary>
    /// <remarks>
    /// Point estimate d = 1 − C^(1/k) (0 for C = 1, 1 for C = 0). The interval (Hera et al. 2023): with q = 1 − (1 − d)^k, the number of mutated k-mers has mean L·q and variance
    /// given by <c>var_n_mutated</c>; the variance of the containment estimate is
    /// term₁·(L·E[N] − E[N²]) + Var(N)/L² with term₁ = (1 − 1/S)/((1/S)·L³·bias²), bias = 1 − (1 − 1/S)^L, and the interval ends are the roots d of
    /// (1 − d)^k ± z·√Var − C = 0 in [1e-7, 1 − 1e-7] (Brent's method, as <c>scipy.optimize.brentq</c>); z = Φ⁻¹(1 − α/2). When no sign change brackets the root (very small sketches)
    /// the interval is null (sourmash warns "Cannot estimate ANI confidence intervals"). <see cref="AniEstimate.ProbabilityNothingInCommon"/> = exp((L − E[N])·ln(1 − 1/S)).
    /// Verified against sourmash 4.9.4 on 1 120 oracle rows (docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md).
    /// </remarks>
    /// <param name="containment">Containment C in [0, 1].</param>
    /// <param name="ksize">k-mer size.</param>
    /// <param name="scaled">The scaled factor S ≥ 1 of the sketches.</param>
    /// <param name="nUniqueKmers">L, the number of distinct k-mers of the contained set (sourmash: sketch hashes × scaled).</param>
    /// <param name="confidence">Confidence level in (0, 1).</param>
    /// <param name="estimateConfidenceInterval">Compute the interval.</param>
    /// <param name="probabilityThreshold">Threshold of <see cref="AniEstimate.ProbabilityNothingInCommon"/>.</param>
    /// <exception cref="ArgumentOutOfRangeException">An argument is out of range.</exception>
    public static AniEstimate ContainmentToDistance(
        double containment, int ksize, double scaled, long nUniqueKmers, double confidence = 0.95, bool estimateConfidenceInterval = false, double probabilityThreshold = 1e-3)
    {
        ValidateAniArguments(containment, ksize, scaled, nUniqueKmers, "containment");
        if (!(confidence > 0.0 && confidence < 1.0))
            throw new ArgumentOutOfRangeException(nameof(confidence), confidence, "confidence must be in (0, 1).");
        double l = nUniqueKmers;
        double point;
        double? low = null, high = null;
        if (containment == 0.0)
        {
            point = 1.0;
            low = 1.0;
            high = 1.0;
        }
        else if (containment == 1.0)
        {
            point = 0.0;
            low = 0.0;
            high = 0.0;
        }
        else
        {
            point = 1.0 - Math.Pow(containment, 1.0 / ksize);
            if (estimateConfidenceInterval)
                (low, high) = ContainmentDistanceInterval(containment, ksize, scaled, l, confidence);
        }

        double pNothing = ProbabilityNothingInCommon(point, ksize, scaled, l);
        return new AniEstimate(point, low, high, pNothing, probabilityThreshold > 0.0 && pNothing > probabilityThreshold, false);
    }

    /// <summary>
    /// sourmash <c>jaccard_to_distance</c>: the point estimate d = 1 − (2J/(1 + J))^(1/k) and the lower bound of its approximation error
    /// (<see cref="AniEstimate.JaccardError"/> = L·Var(N)/(L + E[N])³; above <paramref name="errorThreshold"/> the ANI is withheld).
    /// </summary>
    /// <param name="jaccard">Jaccard index in [0, 1].</param>
    /// <param name="ksize">k-mer size.</param>
    /// <param name="scaled">The scaled factor S ≥ 1.</param>
    /// <param name="nUniqueKmers">L: distinct k-mers (sourmash: the rounded mean sketch size × scaled).</param>
    /// <param name="probabilityThreshold">Threshold of <see cref="AniEstimate.ProbabilityNothingInCommon"/>.</param>
    /// <param name="errorThreshold">Threshold of the error (sourmash default 1e-4).</param>
    /// <exception cref="ArgumentOutOfRangeException">An argument is out of range.</exception>
    public static AniEstimate JaccardToDistance(
        double jaccard, int ksize, double scaled, long nUniqueKmers, double probabilityThreshold = 1e-3, double errorThreshold = 1e-4)
    {
        ValidateAniArguments(jaccard, ksize, scaled, nUniqueKmers, "jaccard");
        double l = nUniqueKmers;
        double point, error;
        if (jaccard == 0.0)
        {
            point = 1.0;
            error = 0.0;
        }
        else if (jaccard == 1.0)
        {
            point = 0.0;
            error = 0.0;
        }
        else
        {
            point = 1.0 - Math.Pow(2.0 * jaccard / (1.0 + jaccard), 1.0 / ksize);
            double expN = ExpectedMutatedKmers(l, ksize, point);
            double varN = VarianceMutatedKmers(l, ksize, point);
            error = 1.0 * l * varN / Math.Pow(l + expN, 3);
        }

        double pNothing = ProbabilityNothingInCommon(point, ksize, scaled, l);
        return new AniEstimate(point, null, null, pNothing, probabilityThreshold > 0.0 && pNothing > probabilityThreshold, false, error, errorThreshold > 0.0 && error > errorThreshold);
    }

    /// <summary>
    /// sourmash <c>MinHash.size_is_accurate</c>: whether the estimate (hashes × scaled) of a set size is within <paramref name="relativeError"/> of the true size with probability at
    /// least <paramref name="confidence"/> (the sketch size is Binomial(set size, 1/scaled); <c>set_size_exact_prob</c>).
    /// </summary>
    /// <param name="sketchHashes">Number of hashes of the sketch.</param>
    /// <param name="scaled">The scaled factor S ≥ 1.</param>
    /// <param name="relativeError">Allowed relative error in [0, 1] (sourmash default 0.20).</param>
    /// <param name="confidence">Required probability in [0, 1] (default 0.95).</param>
    /// <exception cref="ArgumentOutOfRangeException">An argument is out of range.</exception>
    public static bool SketchSizeIsAccurate(int sketchHashes, double scaled, double relativeError = 0.20, double confidence = 0.95)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(sketchHashes);
        if (!(scaled >= 1.0)) throw new ArgumentOutOfRangeException(nameof(scaled), scaled, "scaled must be at least 1.");
        if (!(relativeError >= 0.0 && relativeError <= 1.0)) throw new ArgumentOutOfRangeException(nameof(relativeError), relativeError, "relativeError must be in [0, 1].");
        if (!(confidence >= 0.0 && confidence <= 1.0)) throw new ArgumentOutOfRangeException(nameof(confidence), confidence, "confidence must be in [0, 1].");
        long setSize = (long)Math.Round(sketchHashes * scaled);
        return SetSizeExactProbability(setSize, scaled, relativeError) >= confidence;
    }

    // sourmash set_size_exact_prob: P(|sketch size * scaled - set size| within relative_error * set size), sketch size ~ Binomial(set_size, 1/scaled).
    private static double SetSizeExactProbability(long setSize, double scaled, double relativeError)
    {
        double p = 1.0 / scaled;
        double upper = setSize / scaled * (relativeError + 1);
        double lowerArg = -setSize / scaled * (relativeError - 1);
        double prob = StatisticsHelper.BinomialCdf(upper, setSize, p) - StatisticsHelper.BinomialCdf(lowerArg, setSize, p);
        if (lowerArg == Math.Truncate(lowerArg))
            prob += Math.Exp(StatisticsHelper.BinomialLogDensity(lowerArg, setSize, p));
        return prob;
    }

    private static void ValidateAniArguments(double value, int ksize, double scaled, long nUniqueKmers, string name)
    {
        if (!(value >= 0.0 && value <= 1.0)) throw new ArgumentOutOfRangeException(name, value, $"{name} must be in [0, 1].");
        if (ksize < 1) throw new ArgumentOutOfRangeException(nameof(ksize), ksize, "k-mer size must be at least 1.");
        if (!(scaled >= 1.0)) throw new ArgumentOutOfRangeException(nameof(scaled), scaled, "scaled must be at least 1.");
        if (nUniqueKmers < 0) throw new ArgumentOutOfRangeException(nameof(nUniqueKmers), nUniqueKmers, "The number of distinct k-mers cannot be negative.");
    }

    // sourmash r1_to_q / exp_n_mutated / var_n_mutated (operation order kept: the formulas are numerically delicate).
    private static double MutationToKmerProbability(int k, double r1) => 1.0 - Math.Pow(1.0 - r1, k);

    private static double ExpectedMutatedKmers(double l, int k, double r1) => l * MutationToKmerProbability(k, r1);

    private static double VarianceMutatedKmers(double l, int k, double r1)
    {
        if (r1 == 0.0) return 0.0;
        double q = MutationToKmerProbability(k, r1);
        double varN =
            l * (1 - q) * (q * (2 * k + 2 / r1 - 1) - 2 * k)
            + k * (k - 1) * Math.Pow(1 - q, 2)
            + 2 * (1 - q) / Math.Pow(r1, 2) * ((1 + (k - 1) * (1 - q)) * r1 - q);
        if (varN < 0.0)
            throw new ArgumentException("varN < 0.0 (the sketch is too small for the variance formula).");
        return varN;
    }

    private static double ExpectedMutatedKmersSquared(double l, int k, double p) => VarianceMutatedKmers(l, k, p) + Math.Pow(ExpectedMutatedKmers(l, k, p), 2);

    // sourmash get_exp_probability_nothing_common.
    private static double ProbabilityNothingInCommon(double mutationRate, int ksize, double scaled, double nUniqueKmers)
    {
        if (mutationRate == 1.0) return 1.0;
        if (mutationRate == 0.0) return 0.0;
        double fScaled = 1.0 / scaled;
        double expNmut = ExpectedMutatedKmers(nUniqueKmers, ksize, mutationRate);
        double logP = (nUniqueKmers - expNmut) * Math.Log(1.0 - fScaled);
        return double.IsNaN(logP) ? 0.0 : Math.Exp(logP);
    }

    private static (double? Low, double? High) ContainmentDistanceInterval(double containment, int ksize, double scaled, double l, double confidence)
    {
        try
        {
            double alpha = 1 - confidence;
            double zAlpha = StatisticsHelper.NormalQuantile(1 - alpha / 2);
            double fScaled = 1.0 / scaled;
            double biasFactor = 1 - Math.Pow(1 - fScaled, l);
            double term1 = (1.0 - fScaled) / (fScaled * Math.Pow(l, 3) * Math.Pow(biasFactor, 2));

            double Term2(double pest) => l * ExpectedMutatedKmers(l, ksize, pest) - ExpectedMutatedKmersSquared(l, ksize, pest);
            double Term3(double pest) => VarianceMutatedKmers(l, ksize, pest) / Math.Pow(l, 2);
            double VarDirect(double pest) => term1 * Term2(pest) + Term3(pest);
            double F1(double pest) => Math.Pow(1 - pest, ksize) + zAlpha * Math.Sqrt(VarDirect(pest)) - containment;
            double F2(double pest) => Math.Pow(1 - pest, ksize) - zAlpha * Math.Sqrt(VarDirect(pest)) - containment;

            double sol1 = StatisticsHelper.BrentRoot(F1, 0.0000001, 0.9999999);
            double sol2 = StatisticsHelper.BrentRoot(F2, 0.0000001, 0.9999999);
            return (sol2, sol1);
        }
        catch (ArgumentException)
        {
            return (null, null); // sourmash: ValueError (no bracket, varN < 0) -> sol1 = sol2 = None
        }
    }

    #endregion
}

/// <summary>
/// Options of the sourmash ANI estimates (<c>containment_ani</c> / <c>max_containment_ani</c> / <c>jaccard_ani</c>).
/// </summary>
/// <param name="EstimateConfidenceInterval">Compute the confidence interval of the containment ANI (<c>estimate_ci</c>; two root searches per estimate).</param>
/// <param name="Confidence">Confidence level of the interval, in (0, 1) (sourmash default 0.95).</param>
/// <param name="ProbabilityThreshold">Threshold of the probability of no shared hash by chance (<c>prob_threshold</c>, default 1e-3).</param>
/// <param name="JaccardErrorThreshold">Threshold of the Jaccard-ANI approximation error (<c>err_threshold</c>, default 1e-4); above it the Jaccard ANI is withheld.</param>
public sealed record AniOptions(
    bool EstimateConfidenceInterval = false, double Confidence = 0.95, double ProbabilityThreshold = 1e-3, double JaccardErrorThreshold = 1e-4);

/// <summary>
/// A sourmash ANI result (<c>ciANIResult</c> / <c>jaccardANIResult</c>): the mutation-rate distance d, its confidence interval when estimated, and the quality flags.
/// </summary>
/// <param name="Distance">Point estimate of the distance d = 1 − C^(1/k) (containment) or 1 − (2J/(1 + J))^(1/k) (Jaccard).</param>
/// <param name="DistanceLow">Lower end of the confidence interval of d (null when not requested or not estimable).</param>
/// <param name="DistanceHigh">Upper end of the confidence interval of d.</param>
/// <param name="ProbabilityNothingInCommon">Probability that two sketches of sequences at this distance share no hash by chance (a false-negative check).</param>
/// <param name="ProbabilityExceedsThreshold">The probability exceeds the threshold: the scaled factor may be too large for this comparison.</param>
/// <param name="SizeIsInaccurate">A sketch size estimate (hashes × scaled) deviates from the true set size by more than 20 % with probability above 5 %; the ANI is withheld.</param>
/// <param name="JaccardError">Jaccard ANI only: lower bound of the error of the point estimate; null for containment results.</param>
/// <param name="JaccardErrorExceedsThreshold">Jaccard ANI only: the error exceeds the threshold; the ANI is withheld.</param>
public sealed record AniEstimate(
    double Distance, double? DistanceLow, double? DistanceHigh, double ProbabilityNothingInCommon, bool ProbabilityExceedsThreshold, bool SizeIsInaccurate,
    double? JaccardError = null, bool JaccardErrorExceedsThreshold = false)
{
    /// <summary>ANI = 1 − d; null when a sketch size is inaccurate or (Jaccard) the error is above the threshold (sourmash returns None).</summary>
    public double? Ani => SizeIsInaccurate || JaccardErrorExceedsThreshold ? null : 1.0 - Distance;

    /// <summary>Lower end of the ANI interval, 1 − <see cref="DistanceHigh"/>; null when not estimated or the size is inaccurate.</summary>
    public double? AniLow => DistanceHigh is null || SizeIsInaccurate ? null : 1.0 - DistanceHigh;

    /// <summary>Upper end of the ANI interval, 1 − <see cref="DistanceLow"/>.</summary>
    public double? AniHigh => DistanceLow is null || SizeIsInaccurate ? null : 1.0 - DistanceLow;
}
