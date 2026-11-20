namespace Seqeron.Genomics.Analysis;

// B09 A2-31 / A2-35: sketch-based similarity (Mash bottom-s MinHash, sourmash FracMinHash), containment, Jaccard distance and the
// ANI estimates derived from them. Every number comes from the canonical KmerAnalyzer (B06): exact k-mer sets, MurmurHash3-x64-128
// sketches (seed 42), Mash distance, sourmash-debiased containment; this file only selects the estimator and converts.
public static partial class GenomicAnalyzer
{
    #region Jaccard distance and string overloads

    /// <summary>
    /// Jaccard <b>distance</b> d_J = 1 − J(A, B) of the k-mer sets of two sequences (a metric on sets), as a fraction in [0, 1]; 1 when both k-mer sets are empty
    /// (J is 0 by the convention of <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int, bool)"/>).
    /// </summary>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length, ≥ 1.</param>
    /// <param name="canonical">Strand-agnostic canonical k-mers (default: literal k-mers).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> &lt; 1.</exception>
    public static double CalculateJaccardDistance(DnaSequence sequence1, DnaSequence sequence2, int kmerSize = DefaultKmerSize, bool canonical = false)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return 1.0 - ExactJaccard(sequence1.Sequence, sequence2.Sequence, kmerSize, canonical);
    }

    /// <summary>
    /// <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int, bool)"/> on nucleotide <b>strings</b> (N runs, IUPAC codes, lower case):
    /// a k-mer window that contains a symbol other than A/C/G/T is skipped (it matches nothing; Mash <c>-n</c> / Jellyfish rule), so the result for an
    /// A/C/G/T string equals the <see cref="DnaSequence"/> overload.
    /// </summary>
    /// <param name="sequence1">First nucleotide string.</param>
    /// <param name="sequence2">Second nucleotide string.</param>
    /// <param name="kmerSize">k-mer length, ≥ 1.</param>
    /// <param name="canonical">Strand-agnostic canonical k-mers.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> &lt; 1.</exception>
    public static double CalculateSimilarity(string sequence1, string sequence2, int kmerSize = DefaultKmerSize, bool canonical = false)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ExactJaccard(sequence1, sequence2, kmerSize, canonical) * PercentScale;
    }

    private static double ExactJaccard(string a, string b, int kmerSize, bool canonical)
    {
        if (kmerSize < 1)
            throw new ArgumentOutOfRangeException(nameof(kmerSize), kmerSize, "k-mer size must be at least 1.");
        return KmerAnalyzer.JaccardSimilarity(a, b, kmerSize, new KmerCountingOptions(Canonical: canonical, AcgtOnly: true));
    }

    /// <summary><see cref="CalculateJaccardDistance(DnaSequence, DnaSequence, int, bool)"/> on nucleotide strings (non-A/C/G/T windows are skipped).</summary>
    /// <param name="sequence1">First nucleotide string.</param>
    /// <param name="sequence2">Second nucleotide string.</param>
    /// <param name="kmerSize">k-mer length, ≥ 1.</param>
    /// <param name="canonical">Strand-agnostic canonical k-mers.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> &lt; 1.</exception>
    public static double CalculateJaccardDistance(string sequence1, string sequence2, int kmerSize = DefaultKmerSize, bool canonical = false)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return 1.0 - ExactJaccard(sequence1, sequence2, kmerSize, canonical);
    }

    #endregion

    #region Sketch estimators

    /// <summary>
    /// k-mer Jaccard similarity (percentage) estimated from a bottom-<paramref name="sketchSize"/> MinHash sketch, exactly as <c>mash dist</c> does
    /// (Ondov et al. 2016): each sequence is reduced to its <paramref name="sketchSize"/> smallest MurmurHash3-x64-128 k-mer hashes (seed 42) and J is
    /// the fraction of shared hashes among the smallest <paramref name="sketchSize"/> hashes of the union.
    /// </summary>
    /// <remarks>
    /// Sketches are built by <see cref="KmerAnalyzer.CreateMinHashSketch(string, int, int, bool, uint)"/> and compared by
    /// <see cref="KmerAnalyzer.CompareMinHashSketches"/> (canonical k-mers by default, as <c>mash sketch</c>; windows with a non-ACGT base are skipped). When
    /// <paramref name="sketchSize"/> is at least the number of distinct k-mers of the union the estimate is the exact Jaccard index. For k ≥ 17 the hashes equal
    /// sourmash <c>MinHash(n=s)</c>; for k ≤ 16 Mash keeps only 32 bits of each hash. Use <see cref="CompareSimilarity(DnaSequence, DnaSequence, int, bool, int?, long?)"/>
    /// for the Mash distance and p-value of the same comparison.
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length, 1…32 (Mash range).</param>
    /// <param name="canonical">Canonical k-mers (Mash default); false = Mash <c>-n</c>.</param>
    /// <param name="sketchSize">Sketch size s ≥ 1 (Mash default 1000).</param>
    /// <returns>The estimated Jaccard index as a percentage in [0, 100].</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> outside 1..32 or <paramref name="sketchSize"/> &lt; 1.</exception>
    public static double CalculateSimilarity(DnaSequence sequence1, DnaSequence sequence2, int kmerSize, bool canonical, int sketchSize)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return CompareSimilarity(sequence1.Sequence, sequence2.Sequence, kmerSize, canonical, sketchSize).Jaccard * PercentScale;
    }

    /// <summary>
    /// <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int, bool, int)"/> on nucleotide <b>strings</b> (windows with a symbol other than A/C/G/T are skipped).
    /// </summary>
    /// <param name="sequence1">First nucleotide string.</param>
    /// <param name="sequence2">Second nucleotide string.</param>
    /// <param name="kmerSize">k-mer length, 1…32.</param>
    /// <param name="canonical">Canonical k-mers (Mash default).</param>
    /// <param name="sketchSize">Sketch size s ≥ 1.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> outside 1..32 or <paramref name="sketchSize"/> &lt; 1.</exception>
    public static double CalculateSimilarity(string sequence1, string sequence2, int kmerSize, bool canonical, int sketchSize)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return CompareSimilarityCore(sequence1, sequence2, kmerSize, canonical, sketchSize, null).Jaccard * PercentScale;
    }

    /// <summary>
    /// Compares two sequences by k-mer sets and returns every set-similarity quantity of one estimator: the exact k-mer sets (default), a Mash bottom-s
    /// sketch (<paramref name="sketchSize"/>) or a sourmash FracMinHash "scaled" sketch (<paramref name="scaled"/>).
    /// </summary>
    /// <remarks>
    /// <list type="bullet">
    /// <item><b>Exact</b>: Jaccard, containments and Mash distance from the distinct k-mer sets (<see cref="KmerAnalyzer.JaccardSimilarity(string, string, int, KmerCountingOptions)"/>,
    /// <see cref="KmerAnalyzer.ContainmentIndex(string, string, int, KmerCountingOptions)"/>, <see cref="KmerAnalyzer.MashDistance"/>).</item>
    /// <item><b>MashBottomSketch</b> (<paramref name="sketchSize"/> = s): the <c>mash dist</c> estimate — Jaccard x/s, Mash distance D = −ln(2J/(1+J))/k, p-value; containment is
    /// not defined for bottom-s sketches and is null. k ≤ 32.</item>
    /// <item><b>FracMinHash</b> (<paramref name="scaled"/> = S): sourmash <c>MinHash(n=0, scaled=S)</c> — Jaccard, the bias-corrected containments in both directions
    /// and the max-containment (<see cref="KmerAnalyzer.CompareFracMinHashSketches"/>).</item>
    /// </list>
    /// Derived quantities (sourmash <c>distance_utils</c>, Hera et al. 2023, Genome Res): <see cref="SimilarityEstimate.ContainmentAniAInB"/> = C^(1/k) (point estimate of the
    /// average nucleotide identity from the containment; 0 for C = 0, 1 for C = 1), <see cref="SimilarityEstimate.JaccardAni"/> = (2J/(1+J))^(1/k),
    /// <see cref="SimilarityEstimate.JaccardDistance"/> = 1 − J. Canonical k-mers by default for sketches. Sketch estimators are verified against sourmash 4.9.4
    /// (docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md).
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length ≥ 1 (≤ 32 with <paramref name="sketchSize"/>).</param>
    /// <param name="canonical">Canonical (strand-agnostic) k-mers. Default false for the exact estimator, as <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int)"/>.</param>
    /// <param name="sketchSize">Mash bottom-s sketch size s ≥ 1; null = no Mash sketch. Cannot be combined with <paramref name="scaled"/>.</param>
    /// <param name="scaled">sourmash scaled factor 1…4294967295; null = no FracMinHash. Cannot be combined with <paramref name="sketchSize"/>.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentException">Both <paramref name="sketchSize"/> and <paramref name="scaled"/> are given.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid k, sketch size or scaled.</exception>
    public static SimilarityEstimate CompareSimilarity(
        DnaSequence sequence1, DnaSequence sequence2, int kmerSize = DefaultKmerSize, bool canonical = false, int? sketchSize = null, long? scaled = null, AniOptions? aniOptions = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return CompareSimilarityCore(sequence1.Sequence, sequence2.Sequence, kmerSize, canonical, sketchSize, scaled, aniOptions);
    }

    /// <summary>
    /// <see cref="CompareSimilarity(DnaSequence, DnaSequence, int, bool, int?, long?)"/> on nucleotide <b>strings</b>: k-mer windows with a symbol other than A/C/G/T are skipped.
    /// </summary>
    /// <param name="sequence1">First nucleotide string.</param>
    /// <param name="sequence2">Second nucleotide string.</param>
    /// <param name="kmerSize">k-mer length ≥ 1.</param>
    /// <param name="canonical">Canonical k-mers.</param>
    /// <param name="sketchSize">Mash bottom-s sketch size; null = none.</param>
    /// <param name="scaled">sourmash scaled factor; null = none.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentException">Both <paramref name="sketchSize"/> and <paramref name="scaled"/> are given.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid k, sketch size or scaled.</exception>
    public static SimilarityEstimate CompareSimilarity(
        string sequence1, string sequence2, int kmerSize = DefaultKmerSize, bool canonical = false, int? sketchSize = null, long? scaled = null, AniOptions? aniOptions = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return CompareSimilarityCore(sequence1, sequence2, kmerSize, canonical, sketchSize, scaled, aniOptions);
    }

    private static SimilarityEstimate CompareSimilarityCore(string a, string b, int kmerSize, bool canonical, int? sketchSize, long? scaled, AniOptions? aniOptions = null)
    {
        int k = kmerSize;
        if (k < 1)
            throw new ArgumentOutOfRangeException(nameof(kmerSize), k, "k-mer size must be at least 1.");
        if (sketchSize is not null && scaled is not null)
            throw new ArgumentException("Give either a Mash sketch size or a sourmash scaled factor, not both.", nameof(scaled));

        if (sketchSize is int s)
        {
            var sa = KmerAnalyzer.CreateMinHashSketch(a, k, s, canonical);
            var sb = KmerAnalyzer.CreateMinHashSketch(b, k, s, canonical);
            var mash = KmerAnalyzer.CompareMinHashSketches(sa, sb);
            return new SimilarityEstimate(
                SimilarityEstimator.MashBottomSketch, k, canonical, mash.Jaccard, 1.0 - mash.Jaccard, mash.Distance, JaccardToAni(mash.Jaccard, k),
                null, null, null, null, null, mash.SharedHashes, mash.Denominator, mash.PValue);
        }

        if (scaled is long factor)
        {
            var sa = KmerAnalyzer.CreateFracMinHashSketch(a, k, factor, canonical);
            var sb = KmerAnalyzer.CreateFracMinHashSketch(b, k, factor, canonical);
            var frac = KmerAnalyzer.CompareFracMinHashSketches(sa, sb);
            return new SimilarityEstimate(
                SimilarityEstimator.FracMinHash, k, canonical, frac.Jaccard, 1.0 - frac.Jaccard, KmerAnalyzer.MashDistanceFromJaccard(frac.Jaccard, k), JaccardToAni(frac.Jaccard, k),
                frac.ContainmentAInB, frac.ContainmentBInA, frac.MaxContainment, ContainmentToAni(frac.ContainmentAInB, k), ContainmentToAni(frac.ContainmentBInA, k),
                frac.SharedHashes, frac.UnionHashes, null,
                AniEstimates(sa, sb, frac, k, aniOptions ?? new AniOptions()));
        }

        var options = new KmerCountingOptions(Canonical: canonical, AcgtOnly: true);
        double jaccard = KmerAnalyzer.JaccardSimilarity(a, b, k, options);
        double containmentAB = KmerAnalyzer.ContainmentIndex(a, b, k, options);
        string subject = b, container = a;
        double containmentBA = KmerAnalyzer.ContainmentIndex(subject, container, k, options);
        return new SimilarityEstimate(
            SimilarityEstimator.Exact, k, canonical, jaccard, 1.0 - jaccard, KmerAnalyzer.MashDistance(a, b, k, options), JaccardToAni(jaccard, k),
            containmentAB, containmentBA, Math.Max(containmentAB, containmentBA), ContainmentToAni(containmentAB, k), ContainmentToAni(containmentBA, k),
            null, null, null);
    }

    // sourmash MinHash.containment_ani / max_containment_ani / jaccard_ani: n_unique_kmers = sketch hashes * scaled (not an exact count), size estimates
    // checked with size_is_accurate (a sketch whose size estimate is unreliable withholds its ANI).
    private static SimilarityAniEstimates AniEstimates(FracMinHashSketch sa, FracMinHashSketch sb, FracMinHashComparison frac, int k, AniOptions o)
    {
        int la = sa.Hashes.Count, lb = sb.Hashes.Count;
        double scaled = sa.Scaled;
        bool inaccurate = !SketchSizeIsAccurate(la, scaled) || !SketchSizeIsAccurate(lb, scaled);
        AniEstimate Flag(AniEstimate e) => inaccurate ? e with { SizeIsInaccurate = true } : e;
        AniEstimate Containment(double c, long hashes) => Flag(ContainmentToDistance(c, k, scaled, (long)Math.Round(hashes * scaled), o.Confidence, o.EstimateConfidenceInterval, o.ProbabilityThreshold));
        long avgKmers = (long)Math.Round((la + lb) / 2.0 * scaled, MidpointRounding.ToEven);
        return new SimilarityAniEstimates(
            Containment(frac.ContainmentAInB, la),
            Containment(frac.ContainmentBInA, lb),
            Containment(frac.MaxContainment, Math.Min(la, lb)),
            Flag(JaccardToDistance(frac.Jaccard, k, scaled, avgKmers, o.ProbabilityThreshold, o.JaccardErrorThreshold)));
    }

    // sourmash distance_utils.containment_to_distance point estimate: dist = 1 - C^(1/k); ANI = 1 - dist (C = 0 -> 0, C = 1 -> 1).
    private static double ContainmentToAni(double containment, int k)
    {
        if (containment <= 0.0)
            return 0.0;
        return containment >= 1.0 ? 1.0 : Math.Pow(containment, 1.0 / k);
    }

    // sourmash distance_utils.jaccard_to_distance point estimate: dist = 1 - (2J/(1+J))^(1/k); ANI = 1 - dist (J = 0 -> 0, J = 1 -> 1).
    private static double JaccardToAni(double jaccard, int k)
    {
        if (jaccard <= 0.0)
            return 0.0;
        return jaccard >= 1.0 ? 1.0 : Math.Pow(2.0 * jaccard / (1.0 + jaccard), 1.0 / k);
    }

    #endregion
}

/// <summary>Which set-similarity estimator produced a <see cref="SimilarityEstimate"/>.</summary>
public enum SimilarityEstimator
{
    /// <summary>The exact distinct k-mer sets.</summary>
    Exact = 0,

    /// <summary>A Mash bottom-s MinHash sketch (<c>mash dist</c>).</summary>
    MashBottomSketch = 1,

    /// <summary>A sourmash FracMinHash (scaled) sketch.</summary>
    FracMinHash = 2,
}

/// <summary>
/// Set-similarity of two sequences' k-mers (<see cref="GenomicAnalyzer.CompareSimilarity(DnaSequence, DnaSequence, int, bool, int?, long?)"/>).
/// All fractions are in [0, 1].
/// </summary>
/// <param name="Estimator">The estimator behind the numbers.</param>
/// <param name="KmerSize">k.</param>
/// <param name="Canonical">Canonical (strand-agnostic) k-mers were used.</param>
/// <param name="Jaccard">J = |A ∩ B| / |A ∪ B| (exact, or the sketch estimate).</param>
/// <param name="JaccardDistance">1 − J.</param>
/// <param name="MashDistance">D = −(1/k)·ln(2J/(1+J)) (Ondov et al. 2016, eq. 4; 0 for J = 1, 1 for J = 0).</param>
/// <param name="JaccardAni">sourmash Jaccard ANI point estimate (2J/(1+J))^(1/k) (0 for J = 0, 1 for J = 1).</param>
/// <param name="ContainmentAInB">C(A, B) = |A ∩ B| / |A| (null for a Mash bottom-s sketch).</param>
/// <param name="ContainmentBInA">C(B, A).</param>
/// <param name="MaxContainment">max(C(A, B), C(B, A)).</param>
/// <param name="ContainmentAniAInB">sourmash containment ANI point estimate C(A, B)^(1/k).</param>
/// <param name="ContainmentAniBInA">C(B, A)^(1/k).</param>
/// <param name="SharedHashes">Shared sketch hashes (sketch estimators only).</param>
/// <param name="DenominatorHashes">Union hashes visited (Mash) or the union size (FracMinHash); sketch estimators only.</param>
/// <param name="MashPValue">Probability of at least this many shared hashes by chance (Mash bottom-s sketch only).</param>
/// <param name="AniEstimates">FracMinHash only: the sourmash ANI results with distance, confidence interval (when requested) and quality flags
/// (<see cref="AniOptions"/>); null for the other estimators.</param>
public sealed record SimilarityEstimate(
    SimilarityEstimator Estimator, int KmerSize, bool Canonical, double Jaccard, double JaccardDistance, double MashDistance, double JaccardAni,
    double? ContainmentAInB, double? ContainmentBInA, double? MaxContainment, double? ContainmentAniAInB, double? ContainmentAniBInA,
    int? SharedHashes, int? DenominatorHashes, double? MashPValue, SimilarityAniEstimates? AniEstimates = null);

/// <summary>
/// The sourmash ANI results of a FracMinHash comparison (<c>containment_ani</c> both ways, <c>max_containment_ani</c>, <c>jaccard_ani</c>), each with the point estimate, the
/// confidence interval (<see cref="AniOptions.EstimateConfidenceInterval"/>) and the quality flags of sourmash 4.9.4 (<c>sourmash.distance_utils</c>).
/// </summary>
/// <param name="ContainmentAInB">A.containment_ani(B).</param>
/// <param name="ContainmentBInA">B.containment_ani(A).</param>
/// <param name="MaxContainment">A.max_containment_ani(B).</param>
/// <param name="Jaccard">A.jaccard_ani(B).</param>
public sealed record SimilarityAniEstimates(
    AniEstimate ContainmentAInB, AniEstimate ContainmentBInA, AniEstimate MaxContainment, AniEstimate Jaccard);
