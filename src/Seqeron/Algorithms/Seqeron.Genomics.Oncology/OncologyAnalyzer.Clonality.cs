using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Clonal vs Subclonal Classification (ONCO-CLONAL-001)

    /// <summary>
    /// Cancer-cell-fraction (CCF) threshold separating clonal from subclonal mutations. A mutation is called
    /// clonal when its CCF exceeds this value (with sufficient posterior probability); subclonal otherwise.
    /// Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — "We classified a mutation as clonal if the
    /// CCF harboring it was &gt;0.95 with probability &gt; 0.5, and subclonal otherwise."
    /// </summary>
    public const double ClonalCcfThreshold = 0.95;

    /// <summary>
    /// Minimum posterior probability that the CCF exceeds <see cref="ClonalCcfThreshold"/> required to call a
    /// mutation clonal. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — clonal if P(CCF &gt; 0.95) &gt; 0.5.
    /// </summary>
    public const double ClonalProbabilityThreshold = 0.5;

    /// <summary>
    /// Number of CCF grid points used to evaluate the posterior over the CCF c ∈ [0.01, 1]. Source: Landau et al.
    /// (2013), <i>Cell</i> 152(4):714–726, Extended Experimental Procedures — "calculating these values over a
    /// regular grid of 100 c values and normalizing by dividing them by their sum".
    /// </summary>
    private const int CcfGridPointCount = 100;

    /// <summary>
    /// Lower bound of the CCF grid. Source: Landau et al. (2013) — "with c ∈ [0.01,1]" (a mutation is present in
    /// at least one cancer cell, so the CCF cannot be exactly zero).
    /// </summary>
    private const double CcfGridLowerBound = 0.01;

    /// <summary>Upper bound of the CCF grid (a mutation present in every cancer cell). Source: Landau et al. (2013), c ∈ [0.01,1].</summary>
    private const double CcfGridUpperBound = 1.0;

    /// <summary>Clonal/subclonal status of one somatic mutation.</summary>
    public enum ClonalityStatus
    {
        /// <summary>Present in (essentially) all cancer cells: CCF &gt; 0.95 with posterior probability &gt; 0.5.</summary>
        Clonal,

        /// <summary>Present in only a subpopulation of cancer cells: it does not meet the clonal criterion.</summary>
        Subclonal,
    }

    /// <summary>
    /// Read evidence for one somatic mutation at a locus with known purity-corrected copy-number state, used to
    /// infer its cancer cell fraction. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — the posterior
    /// over CCF is built from <c>a</c> alternate reads out of <c>N</c> total reads at a locus of absolute somatic
    /// copy number <c>q</c>, with mutation multiplicity <c>M</c> (the number of tumour-genome copies carrying the
    /// mutant allele; DeCiFering / Satas et al. 2021, <i>Cell Systems</i> 12(10):1004–1018, Eq. 1).
    /// </summary>
    /// <param name="AltReads">Alternate-allele supporting reads <c>a</c> (≥ 0, ≤ <paramref name="TotalReads"/>).</param>
    /// <param name="TotalReads">Total reads at the locus <c>N</c> (≥ 1).</param>
    /// <param name="LocalCopyNumber">Absolute tumour total copy number at the locus <c>q</c> (≥ 1).</param>
    /// <param name="Multiplicity">Mutation multiplicity <c>M</c> — tumour-genome copies carrying the mutation (≥ 1, ≤ <paramref name="LocalCopyNumber"/>).</param>
    public readonly record struct ClonalityVariant(int AltReads, int TotalReads, int LocalCopyNumber, int Multiplicity)
    {
        /// <summary>Creates a variant at multiplicity 1 (one mutated copy), the default for a heterozygous SNV.</summary>
        public ClonalityVariant(int altReads, int totalReads, int localCopyNumber)
            : this(altReads, totalReads, localCopyNumber, 1)
        {
        }
    }

    /// <summary>The clonality classification of one mutation, with the CCF point estimate and posterior that produced it.</summary>
    /// <param name="Variant">The variant that was classified.</param>
    /// <param name="Ccf">Posterior-mean cancer cell fraction estimate (grid expectation), in [0.01, 1].</param>
    /// <param name="ProbabilityClonal">Posterior probability that the CCF exceeds <see cref="ClonalCcfThreshold"/>.</param>
    /// <param name="Status">Clonal / Subclonal classification.</param>
    public readonly record struct ClonalityCall(ClonalityVariant Variant, double Ccf, double ProbabilityClonal, ClonalityStatus Status);

    /// <summary>Summary of a clonal/subclonal classification over a set of variants.</summary>
    /// <param name="Calls">Per-variant classifications, in input order.</param>
    /// <param name="ClonalCount">Number of variants classified <see cref="ClonalityStatus.Clonal"/>.</param>
    /// <param name="SubclonalCount">Number of variants classified <see cref="ClonalityStatus.Subclonal"/>.</param>
    /// <param name="ClonalFraction">Fraction of variants that are clonal (ClonalCount / total); 0 for an empty set.</param>
    public readonly record struct ClonalityResult(
        IReadOnlyList<ClonalityCall> Calls,
        int ClonalCount,
        int SubclonalCount,
        double ClonalFraction);

    /// <summary>
    /// Classifies each somatic mutation as clonal or subclonal from its read evidence and the tumour purity.
    /// For each variant the posterior over the cancer cell fraction c is built on a regular grid of
    /// <see cref="CcfGridPointCount"/> points c ∈ [<see cref="CcfGridLowerBound"/>, 1] as
    /// <c>P(c) ∝ Binomial(a | N, f(c))</c> with a uniform prior, where the expected alternate-allele fraction is
    /// <c>f(c) = ρ·M·c / (2(1−ρ) + ρ·q)</c> — purity ρ, multiplicity M, local copy number q, normal diploid
    /// contribution 2(1−ρ). A mutation is called <see cref="ClonalityStatus.Clonal"/> when the posterior
    /// probability that c &gt; <see cref="ClonalCcfThreshold"/> exceeds <see cref="ClonalProbabilityThreshold"/>,
    /// and <see cref="ClonalityStatus.Subclonal"/> otherwise. Source: Landau et al. (2013), <i>Cell</i>
    /// 152(4):714–726 (Extended Experimental Procedures); multiplicity generalisation per Satas et al. (2021),
    /// <i>Cell Systems</i> 12(10):1004–1018 (Eq. 1).
    /// </summary>
    /// <param name="variants">Per-variant read evidence and copy-number state.</param>
    /// <param name="purity">Tumor purity ρ ∈ (0, 1].</param>
    /// <returns>Per-variant calls plus clonal/subclonal counts and the clonal fraction.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">purity ∉ (0, 1].</exception>
    /// <exception cref="ArgumentException">A variant has invalid read counts, copy number, or multiplicity.</exception>
    public static ClonalityResult ClassifyClonality(IEnumerable<ClonalityVariant> variants, double purity)
    {
        ArgumentNullException.ThrowIfNull(variants);
        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in the range (0, 1]; the CCF model divides by purity.");
        }

        var calls = new List<ClonalityCall>();
        int clonalCount = 0;
        foreach (ClonalityVariant variant in variants)
        {
            ClonalityCall call = ClassifyOne(variant, purity);
            calls.Add(call);
            if (call.Status == ClonalityStatus.Clonal)
            {
                clonalCount++;
            }
        }

        int total = calls.Count;
        int subclonalCount = total - clonalCount;
        // ClonalFraction is the clonal share; undefined for an empty set so reported as 0.
        double clonalFraction = total == 0 ? 0.0 : (double)clonalCount / total;
        return new ClonalityResult(calls, clonalCount, subclonalCount, clonalFraction);
    }

    /// <summary>
    /// Classifies a set of already-estimated cancer cell fractions (CCF point estimates, e.g. from
    /// <c>EstimateCCF</c>) as clonal vs subclonal and returns the indices of the clonal mutations. A CCF is
    /// clonal when it exceeds <see cref="ClonalCcfThreshold"/> (CCF &gt; 0.95), reflecting a mutation present in
    /// (essentially) all cancer cells. Source: Landau et al. (2013), <i>Cell</i> 152(4):714–726 — "classified a
    /// mutation as clonal if the CCF harboring it was &gt;0.95 … and subclonal otherwise".
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions, each in [0, 1].</param>
    /// <returns>The 0-based indices of the clonal CCF values, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">A CCF value is NaN or outside [0, 1].</exception>
    public static IReadOnlyList<int> IdentifyClonalMutations(IEnumerable<double> ccfValues)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);

        var clonalIndices = new List<int>();
        int index = 0;
        foreach (double ccf in ccfValues)
        {
            if (double.IsNaN(ccf) || ccf < 0.0 || ccf > 1.0)
            {
                throw new ArgumentException(
                    $"Cancer cell fraction must be in [0, 1]; got {ccf} at index {index}.", nameof(ccfValues));
            }

            if (ccf > ClonalCcfThreshold)
            {
                clonalIndices.Add(index);
            }

            index++;
        }

        return clonalIndices;
    }

    /// <summary>
    /// Builds the posterior over CCF for one variant and classifies it (Landau et al. 2013 grid model).
    /// </summary>
    private static ClonalityCall ClassifyOne(ClonalityVariant variant, double purity)
    {
        ValidateClonalityVariant(variant);

        // Expected alt-allele fraction f(c) = ρ·M·c / (2(1−ρ) + ρ·q): mutant copies M·c per cell scaled by
        // purity over the total DNA (normal 2(1−ρ) + tumour ρ·q). Landau 2013 (M=1) generalised by DeCiFering Eq.1.
        double denominator = NormalDiploidCopyNumber * (1.0 - purity) + purity * variant.LocalCopyNumber;
        double alleleFractionPerUnitCcf = purity * variant.Multiplicity / denominator;

        // Posterior P(c) ∝ Binomial(a | N, f(c)) on a uniform grid c ∈ [0.01, 1], uniform prior; normalise by sum.
        // The binomial coefficient C(N, a) is constant in c, so it cancels in normalisation and is omitted.
        double step = (CcfGridUpperBound - CcfGridLowerBound) / (CcfGridPointCount - 1);
        Span<double> weights = stackalloc double[CcfGridPointCount];
        double weightSum = 0.0;
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double c = CcfGridLowerBound + step * i;
            double f = Math.Min(1.0, alleleFractionPerUnitCcf * c);
            double likelihood = BinomialLikelihoodKernel(variant.AltReads, variant.TotalReads, f);
            weights[i] = likelihood;
            weightSum += likelihood;
        }

        double ccfMean = 0.0;
        double probabilityClonal = 0.0;
        // Guard against an all-zero posterior (e.g. f≈0 with a>0): fall back to a flat posterior over the grid.
        bool degenerate = weightSum <= 0.0 || double.IsNaN(weightSum);
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double c = CcfGridLowerBound + step * i;
            double posterior = degenerate ? 1.0 / CcfGridPointCount : weights[i] / weightSum;
            ccfMean += c * posterior;
            if (c > ClonalCcfThreshold)
            {
                probabilityClonal += posterior;
            }
        }

        ClonalityStatus status = probabilityClonal > ClonalProbabilityThreshold
            ? ClonalityStatus.Clonal
            : ClonalityStatus.Subclonal;

        // The reported posterior summaries are bounded by their documented invariants: the CCF point estimate is a
        // grid expectation in [0.01, 1] (INV-03) and the clonal probability is normalised posterior mass in [0, 1]
        // (INV-04). Summing the normalised grid weights can overshoot the bound by one ulp (e.g. 1.0000000000000002),
        // so clamp the reported values to their invariant range. The unclamped probabilityClonal already determined
        // status above, and clamping a near-1 value cannot change the > 0.5 decision.
        double reportedCcf = Math.Clamp(ccfMean, CcfGridLowerBound, CcfGridUpperBound);
        double reportedProbabilityClonal = Math.Clamp(probabilityClonal, 0.0, 1.0);
        return new ClonalityCall(variant, reportedCcf, reportedProbabilityClonal, status);
    }

    /// <summary>
    /// Binomial likelihood kernel L(a | N, p) = p^a · (1−p)^(N−a), without the constant C(N, a) factor (it cancels
    /// under grid normalisation). Computed via log-space to avoid underflow for large N.
    /// </summary>
    private static double BinomialLikelihoodKernel(int altReads, int totalReads, double p)
    {
        int refReads = totalReads - altReads;
        if (p <= 0.0)
        {
            // p = 0 explains zero alternate reads exactly, nothing else.
            return altReads == 0 ? 1.0 : 0.0;
        }

        if (p >= 1.0)
        {
            // p = 1 explains all-alternate reads exactly, nothing else.
            return refReads == 0 ? 1.0 : 0.0;
        }

        double logLikelihood = altReads * Math.Log(p) + refReads * Math.Log(1.0 - p);
        return Math.Exp(logLikelihood);
    }

    /// <summary>Validates read counts, local copy number, and multiplicity of a clonality variant.</summary>
    private static void ValidateClonalityVariant(ClonalityVariant variant)
    {
        if (variant.TotalReads < 1)
        {
            throw new ArgumentException(
                $"Total reads must be at least 1; got {variant.TotalReads}.", nameof(variant));
        }

        if (variant.AltReads < 0 || variant.AltReads > variant.TotalReads)
        {
            throw new ArgumentException(
                $"Alternate reads must be in [0, {variant.TotalReads}]; got {variant.AltReads}.", nameof(variant));
        }

        if (variant.LocalCopyNumber < 1)
        {
            throw new ArgumentException(
                $"Local copy number must be at least 1; got {variant.LocalCopyNumber}.", nameof(variant));
        }

        if (variant.Multiplicity < 1 || variant.Multiplicity > variant.LocalCopyNumber)
        {
            throw new ArgumentException(
                $"Multiplicity must be in [1, {variant.LocalCopyNumber}]; got {variant.Multiplicity}.", nameof(variant));
        }
    }

    #endregion


    #region EstimateCcf

    /// <summary>
    /// Normal locus copy number (diploid) contributing the 2(1−ρ) term in the CCF denominator.
    /// Source: McGranahan et al. (2016), <i>Science</i> 351(6280):1463–1469 — n_mut = VAF·(1/p)·[p·CN_t + CN_n·(1−p)]
    /// with CN_n = 2; Tarabichi et al. (2021), <i>Nat. Methods</i> 18:144–155 (Box 1).
    /// </summary>
    private const double NormalLocusCopyNumber = 2.0;

    /// <summary>Upper bound on a reported cancer cell fraction (a mutation in all cancer cells has CCF = 1).</summary>
    private const double MaxCancerCellFraction = 1.0;

    /// <summary>
    /// A single cancer-cell-fraction point estimate for one somatic mutation.
    /// </summary>
    /// <param name="Ccf">Reported cancer cell fraction, capped to [0, 1] (the registry invariant; a mutation
    /// present in all cancer cells has CCF = 1, per McGranahan et al. 2016).</param>
    /// <param name="RawCcf">Uncapped formula value VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m); may exceed 1 under sampling
    /// noise (CNAqc reports e.g. 1.06).</param>
    public readonly record struct CcfEstimate(double Ccf, double RawCcf);

    /// <summary>
    /// Result of clustering cancer cell fractions into clones/subclones by deterministic 1D k-means.
    /// </summary>
    /// <param name="Centroids">Cluster centroids (means) sorted in ascending order, one per cluster.</param>
    /// <param name="Assignments">For each input CCF (input order), the 0-based index of its assigned cluster
    /// in <paramref name="Centroids"/>.</param>
    /// <param name="ClonalClusterIndex">Index (into <paramref name="Centroids"/>) of the clonal cluster — the
    /// cluster with the highest centroid (Tarabichi et al. 2021: "the cluster with the highest CP … deemed clonal").</param>
    public readonly record struct CcfClustering(
        IReadOnlyList<double> Centroids,
        IReadOnlyList<int> Assignments,
        int ClonalClusterIndex);

    /// <summary>
    /// Estimates the cancer cell fraction (CCF) of a somatic mutation from its variant allele fraction, the tumor
    /// purity, the local tumor copy number, and the mutation multiplicity, using the standard point estimate
    /// <c>CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m)</c> — VAF the variant allele fraction, ρ the purity, N_T the local
    /// tumor copy number, the normal contributing 2(1−ρ), and m the integer mutation multiplicity (number of
    /// mutated copies per cancer cell). Source: McGranahan et al. (2016), <i>Science</i> 351(6280):1463–1469
    /// (n_mut = VAF·(1/p)·[p·CN_t + 2(1−p)], CCF = n_mut/m); Tarabichi et al. (2021), <i>Nat. Methods</i>
    /// 18:144–155 (Box 1); Zheng et al. (2022), <i>Bioinformatics</i> 38(15):3677–3683 (VAF = m·CCF·p/(c·p+2(1−p))).
    /// The reported <see cref="CcfEstimate.Ccf"/> is capped to [0, 1] to honour the 0 ≤ CCF ≤ 1 invariant while the
    /// uncapped value is exposed as <see cref="CcfEstimate.RawCcf"/>.
    /// </summary>
    /// <param name="vaf">Variant allele fraction (mutant read fraction), in [0, 1].</param>
    /// <param name="purity">Tumor purity ρ ∈ (0, 1].</param>
    /// <param name="tumorCopyNumber">Local tumor copy number N_T (≥ 1).</param>
    /// <param name="multiplicity">Mutation multiplicity m (number of mutated copies per cancer cell), in
    /// [1, <paramref name="tumorCopyNumber"/>].</param>
    /// <returns>The capped and raw cancer cell fraction.</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0,1], purity ∉ (0,1], or tumorCopyNumber &lt; 1.</exception>
    /// <exception cref="ArgumentException">multiplicity ∉ [1, tumorCopyNumber].</exception>
    public static CcfEstimate EstimateCcf(double vaf, double purity, int tumorCopyNumber, int multiplicity)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in (0, 1]; the CCF formula divides by purity.");
        }

        if (tumorCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(tumorCopyNumber), tumorCopyNumber, "Tumor copy number must be at least 1.");
        }

        if (multiplicity < 1 || multiplicity > tumorCopyNumber)
        {
            throw new ArgumentException(
                $"Multiplicity must be in [1, {tumorCopyNumber}]; got {multiplicity}.", nameof(multiplicity));
        }

        // CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m): total DNA per cell = tumour ρ·N_T + normal 2(1−ρ); dividing the
        // observed mutant fraction by ρ·m / totalDna recovers the fraction of cancer cells carrying the mutation.
        double totalDnaPerCell = purity * tumorCopyNumber + NormalLocusCopyNumber * (1.0 - purity);
        double rawCcf = vaf * totalDnaPerCell / (purity * multiplicity);
        double cappedCcf = Math.Min(MaxCancerCellFraction, rawCcf);
        return new CcfEstimate(cappedCcf, rawCcf);
    }

    /// <summary>
    /// Clusters cancer cell fractions into <paramref name="clusterCount"/> clones/subclones using a deterministic
    /// one-dimensional Lloyd's k-means (Lloyd 1982, <i>IEEE Trans. Inf. Theory</i> 28(2):129–137): each value is
    /// assigned to the nearest centroid (least squared distance) and each centroid is recomputed as the mean of
    /// its members, iterating to convergence. Determinism is achieved without any RNG by seeding the centroids at
    /// evenly spaced quantiles of the sorted input. The clonal cluster is the one with the highest centroid
    /// (Tarabichi et al. 2021, <i>Nat. Methods</i> 18:144–155: "the cluster with the highest CP can be deemed clonal").
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions to cluster (each finite).</param>
    /// <param name="clusterCount">Number of clusters k, in [1, count of values].</param>
    /// <returns>Ascending-sorted centroids, per-value cluster assignments (input order), and the clonal cluster index.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, or a value is NaN/infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count].</exception>
    public static CcfClustering ClusterCcfValues(IReadOnlyList<double> ccfValues, int clusterCount)
    {
        ArgumentNullException.ThrowIfNull(ccfValues);

        int n = ccfValues.Count;
        if (n == 0)
        {
            throw new ArgumentException("At least one CCF value is required.", nameof(ccfValues));
        }

        for (int i = 0; i < n; i++)
        {
            if (double.IsNaN(ccfValues[i]) || double.IsInfinity(ccfValues[i]))
            {
                throw new ArgumentException($"CCF value must be finite; got {ccfValues[i]} at index {i}.", nameof(ccfValues));
            }
        }

        if (clusterCount < 1 || clusterCount > n)
        {
            throw new ArgumentOutOfRangeException(
                nameof(clusterCount), clusterCount, $"Cluster count must be in [1, {n}].");
        }

        // Sort values (carrying original indices) so seeding and assignment are deterministic and order-independent.
        int[] order = Enumerable.Range(0, n).OrderBy(i => ccfValues[i]).ToArray();
        double[] sorted = order.Select(i => ccfValues[i]).ToArray();

        // Deterministic seeding: place centroid j at the value at quantile (j + 0.5)/k of the sorted data.
        double[] centroids = new double[clusterCount];
        for (int j = 0; j < clusterCount; j++)
        {
            int idx = (int)((j + 0.5) / clusterCount * n);
            if (idx >= n)
            {
                idx = n - 1;
            }

            centroids[j] = sorted[idx];
        }

        int[] sortedAssignments = new int[n];
        // Lloyd iterations: assignment step then update step; converges when no assignment changes.
        // Bounded by n iterations (each iteration strictly reduces WCSS until a fixed point on a finite set).
        for (int iteration = 0; iteration <= n; iteration++)
        {
            bool changed = AssignToNearestCentroid(sorted, centroids, sortedAssignments);
            RecomputeCentroids(sorted, sortedAssignments, centroids);
            if (!changed && iteration > 0)
            {
                break;
            }
        }

        // Map cluster labels to ascending-centroid order so the result is canonical and assignments are stable.
        int[] centroidRank = Enumerable.Range(0, clusterCount).OrderBy(j => centroids[j]).ToArray();
        int[] rankOfCluster = new int[clusterCount];
        double[] sortedCentroids = new double[clusterCount];
        for (int rank = 0; rank < clusterCount; rank++)
        {
            rankOfCluster[centroidRank[rank]] = rank;
            sortedCentroids[rank] = centroids[centroidRank[rank]];
        }

        int[] assignments = new int[n];
        for (int s = 0; s < n; s++)
        {
            assignments[order[s]] = rankOfCluster[sortedAssignments[s]];
        }

        // Clonal cluster = highest centroid; centroids are ascending so it is the last index.
        int clonalClusterIndex = clusterCount - 1;
        return new CcfClustering(sortedCentroids, assignments, clonalClusterIndex);
    }

    /// <summary>Assignment step: assigns each value to the nearest centroid; returns whether any assignment changed.</summary>
    private static bool AssignToNearestCentroid(double[] values, double[] centroids, int[] assignments)
    {
        bool changed = false;
        for (int i = 0; i < values.Length; i++)
        {
            int best = 0;
            double bestDistance = double.PositiveInfinity;
            for (int j = 0; j < centroids.Length; j++)
            {
                double diff = values[i] - centroids[j];
                double distance = diff * diff;
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    best = j;
                }
            }

            if (assignments[i] != best)
            {
                assignments[i] = best;
                changed = true;
            }
        }

        return changed;
    }

    /// <summary>Update step: recomputes each centroid as the mean of its assigned values (empty clusters keep their centroid).</summary>
    private static void RecomputeCentroids(double[] values, int[] assignments, double[] centroids)
    {
        Span<double> sums = stackalloc double[centroids.Length];
        Span<int> counts = stackalloc int[centroids.Length];
        for (int i = 0; i < values.Length; i++)
        {
            sums[assignments[i]] += values[i];
            counts[assignments[i]]++;
        }

        for (int j = 0; j < centroids.Length; j++)
        {
            if (counts[j] > 0)
            {
                centroids[j] = sums[j] / counts[j];
            }
        }
    }

    #endregion


    #region Clonal Hematopoiesis Filtering (ONCO-CHIP-001)

    /// <summary>
    /// Minimum variant allele fraction (VAF) for a driver-gene somatic mutation in blood to meet the
    /// clonal-hematopoiesis-of-indeterminate-potential (CHIP) definition. Source: Steensma et al. (2015),
    /// <i>Blood</i> 126(1):9–16 — "the mutant allele fraction must be ≥2% in the peripheral blood"
    /// (threshold is inclusive: ≥ 0.02).
    /// </summary>
    public const double ChipVafThreshold = 0.02;

    /// <summary>
    /// Default canonical CHIP driver-gene panel (HGNC symbols, upper-case): genes recurrently mutated in
    /// clonal hematopoiesis. Source: Steensma et al. (2015) Fig. 2A and Genovese et al. (2014),
    /// <i>NEJM</i> 371(26):2477–2487 ("Four genes (DNMT3A, TET2, ASXL1, and PPM1D) had disproportionately
    /// high numbers of somatic mutations"; JAK2 V617F and SF3B1 K700E recurrent; SRSF2/TP53 are established
    /// CHIP drivers). This is a labelled canonical set, NOT an invented value — callers may override it via
    /// the <c>chipGenes</c> parameter (the algorithm is gene-panel-agnostic, per Razavi et al. 2019).
    /// </summary>
    public static readonly IReadOnlyCollection<string> DefaultChipGenes = new[]
    {
        "DNMT3A", "TET2", "ASXL1", "TP53", "JAK2", "SF3B1", "SRSF2", "PPM1D"
    };

    /// <summary>
    /// A variant observed in plasma cell-free DNA (cfDNA) for CHIP analysis, carrying its locus, gene
    /// symbol, and plasma VAF, plus optional matched white-blood-cell (WBC) alt-read evidence used by
    /// <see cref="FilterCHIP(IEnumerable{ChipVariant}, IEnumerable{ChipVariant}, IReadOnlyCollection{string}?, double, int)"/>.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier of the locus.</param>
    /// <param name="Position">1-based reference position.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Alternate (mutant) allele.</param>
    /// <param name="Gene">HGNC gene symbol the variant falls in (case-insensitive on comparison).</param>
    /// <param name="Vaf">Plasma variant allele fraction in [0, 1].</param>
    /// <param name="AltReads">Alternate (mutant) supporting reads at this locus (≥ 0); used as WBC evidence.</param>
    public readonly record struct ChipVariant(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        string Gene,
        double Vaf,
        int AltReads = 0);

    /// <summary>
    /// Reports whether a gene symbol belongs to the CHIP driver-gene panel (case-insensitive). Source:
    /// Steensma et al. (2015) / Genovese et al. (2014) canonical driver genes; the panel is caller-supplied
    /// when <paramref name="chipGenes"/> is provided, otherwise <see cref="DefaultChipGenes"/>.
    /// </summary>
    /// <param name="gene">Gene symbol to test (null/empty ⇒ <c>false</c>).</param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <returns><c>true</c> when <paramref name="gene"/> is in the panel; otherwise <c>false</c>.</returns>
    public static bool IsCanonicalChipGene(string? gene, IReadOnlyCollection<string>? chipGenes = null)
    {
        if (string.IsNullOrEmpty(gene))
        {
            return false;
        }

        IReadOnlyCollection<string> panel = chipGenes ?? DefaultChipGenes;
        foreach (string g in panel)
        {
            if (string.Equals(g, gene, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Identifies candidate clonal-hematopoiesis (CHIP) variants by the gene + VAF heuristic: a variant is
    /// flagged CHIP when its gene is in the CHIP driver panel AND its plasma VAF is at or above
    /// <paramref name="minVaf"/> (default <see cref="ChipVafThreshold"/> = 0.02). Source: Steensma et al.
    /// (2015) — a somatic mutation in a gene recurrently mutated in hematologic malignancies at VAF ≥ 2%.
    /// This is a candidate flag; the definitive tumour-vs-CH origin test is matched-WBC subtraction
    /// (<see cref="FilterCHIP(IEnumerable{ChipVariant}, IEnumerable{ChipVariant}, IReadOnlyCollection{string}?, double, int)"/>,
    /// Razavi et al. 2019).
    /// </summary>
    /// <param name="variants">cfDNA variants to screen (non-null).</param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <param name="minVaf">Minimum VAF to meet the CHIP definition (default <see cref="ChipVafThreshold"/>); must be in (0, 1].</param>
    /// <returns>The subset of <paramref name="variants"/> flagged as candidate CHIP, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minVaf"/> ∉ (0, 1].</exception>
    public static IReadOnlyList<ChipVariant> IdentifyCHIPVariants(
        IEnumerable<ChipVariant> variants,
        IReadOnlyCollection<string>? chipGenes = null,
        double minVaf = ChipVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);

        if (!(minVaf > 0.0) || minVaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minVaf), minVaf, "Minimum CHIP VAF must be in the interval (0, 1].");
        }

        var result = new List<ChipVariant>();
        foreach (ChipVariant variant in variants)
        {
            if (IsCanonicalChipGene(variant.Gene, chipGenes) && variant.Vaf >= minVaf)
            {
                result.Add(variant);
            }
        }

        return result;
    }

    /// <summary>
    /// Removes clonal-hematopoiesis (CHIP) confounder variants from a plasma cfDNA call set so that the
    /// retained variants are candidate tumour-derived variants. A cfDNA variant is removed when EITHER
    /// (a) it is also detected in the matched white-blood-cell (WBC) sample at the same locus — the
    /// definitive matched-WBC origin test (Razavi et al. 2019, <i>Nat Med</i> 25:1928–1937: matched
    /// cfDNA–WBC sequencing assigns variant origin tumour-vs-CH) — OR (b) it meets the gene + VAF CHIP
    /// heuristic (<see cref="IdentifyCHIPVariants"/>, Steensma et al. 2015). Rule (a) applies regardless of
    /// gene; rule (b) is the fallback when no matched-WBC evidence exists. Output is a subset of the input
    /// in input order.
    /// </summary>
    /// <param name="variants">Plasma cfDNA variants to filter (non-null).</param>
    /// <param name="whiteBloodCellVariants">
    /// Matched WBC variants; a cfDNA variant sharing a locus (chromosome, position, ref, alt) with a WBC
    /// variant carrying ≥ <paramref name="minWbcAltReads"/> alt reads is treated as WBC/CH-derived.
    /// </param>
    /// <param name="chipGenes">Optional caller-supplied CHIP panel; defaults to <see cref="DefaultChipGenes"/>.</param>
    /// <param name="minVaf">Minimum VAF for the gene+VAF CHIP heuristic (default <see cref="ChipVafThreshold"/>); in (0, 1].</param>
    /// <param name="minWbcAltReads">Minimum alt reads in matched WBC to count the locus as present (default 1; Wan et al. 2020).</param>
    /// <returns>cfDNA variants retained as candidate tumour-derived, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> or <paramref name="whiteBloodCellVariants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minVaf"/> ∉ (0, 1], or <paramref name="minWbcAltReads"/> &lt; 1.</exception>
    public static IReadOnlyList<ChipVariant> FilterCHIP(
        IEnumerable<ChipVariant> variants,
        IEnumerable<ChipVariant> whiteBloodCellVariants,
        IReadOnlyCollection<string>? chipGenes = null,
        double minVaf = ChipVafThreshold,
        int minWbcAltReads = DefaultMrdMinSupportingReads)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ArgumentNullException.ThrowIfNull(whiteBloodCellVariants);

        if (!(minVaf > 0.0) || minVaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minVaf), minVaf, "Minimum CHIP VAF must be in the interval (0, 1].");
        }

        if (minWbcAltReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minWbcAltReads), minWbcAltReads, "Minimum WBC alt reads must be at least 1.");
        }

        // Build the set of loci present in the matched WBC (alt reads >= cutoff) for O(1) lookup.
        var wbcLoci = new HashSet<(string, int, string, string)>();
        foreach (ChipVariant wbc in whiteBloodCellVariants)
        {
            if (wbc.AltReads >= minWbcAltReads)
            {
                wbcLoci.Add(LocusKey(wbc));
            }
        }

        var retained = new List<ChipVariant>();
        foreach (ChipVariant variant in variants)
        {
            bool inMatchedWbc = wbcLoci.Contains(LocusKey(variant));
            bool meetsChipHeuristic = IsCanonicalChipGene(variant.Gene, chipGenes) && variant.Vaf >= minVaf;
            if (!inMatchedWbc && !meetsChipHeuristic)
            {
                retained.Add(variant);
            }
        }

        return retained;
    }

    /// <summary>Locus identity key (chromosome, 1-based position, ref allele, alt allele) for matched-WBC subtraction.</summary>
    private static (string, int, string, string) LocusKey(ChipVariant v) =>
        (v.Chromosome, v.Position, v.ReferenceAllele, v.AlternateAllele);

    /// <summary>Locus identity key (chromosome, 1-based position, ref allele, alt allele) for a matched-WBC observation.</summary>
    private static (string, int, string, string) LocusKey(WbcObservation o) =>
        (o.Chromosome, o.Position, o.ReferenceAllele, o.AlternateAllele);

    /// <summary>
    /// Minimum supporting alternate reads in the matched white-blood-cell (WBC) sample for a variant to count
    /// as a confident clonal-hematopoiesis (CH) call. Source: Bolton et al. (2020), <i>Nat Genet</i>
    /// 52(11):1219–1226 — CH mutations were defined as having "a variant allele fraction of at least 2% and
    /// at least 10 supporting reads."
    /// </summary>
    public const int ChipMinWbcSupportingReads = 10;

    /// <summary>
    /// Default WBC-to-tumour VAF fold ratio at or above which a variant is called WBC / clonal-hematopoiesis
    /// origin rather than tumour-derived. Source: Bolton et al. (2020), <i>Nat Genet</i> 52(11):1219–1226 —
    /// a variant detected in blood with "a VAF of at least twice that in the tumor … were considered" CH
    /// (the ratio "was chosen … through simulations of leukocyte contamination in the tumor").
    /// </summary>
    public const double DefaultWbcVafFold = 2.0;

    /// <summary>
    /// WBC-to-tumour VAF fold ratio for lymph-node tumour biopsy sites, where leukocyte admixture is higher.
    /// Source: Bolton et al. (2020), <i>Nat Genet</i> 52(11):1219–1226 — "1.5 times the VAF if the tumor
    /// biopsy site was a lymph node."
    /// </summary>
    public const double LymphNodeWbcVafFold = 1.5;

    /// <summary>
    /// The called origin of a tumour/plasma variant under strict matched-WBC origin calling.
    /// Source: Bolton et al. (2020), <i>Nat Genet</i> 52(11):1219–1226; Razavi et al. (2019), <i>Nat Med</i>
    /// 25:1928–1937.
    /// </summary>
    public enum VariantOrigin
    {
        /// <summary>The variant is tumour-derived (somatic): not confidently present in the matched WBC.</summary>
        Tumor = 0,

        /// <summary>
        /// The variant is white-blood-cell / clonal-hematopoiesis derived: present in the matched WBC at a
        /// VAF ≥ <see cref="ChipVafThreshold"/> with ≥ <see cref="ChipMinWbcSupportingReads"/> supporting
        /// reads, and a WBC-to-tumour VAF ratio ≥ the configured fold (Bolton et al. 2020).
        /// </summary>
        Chip = 1
    }

    /// <summary>
    /// A variant observation in the matched white-blood-cell (WBC) sample: the locus plus the WBC variant
    /// allele fraction and supporting alt-read count used to call origin in
    /// <see cref="CallVariantOrigin(IEnumerable{ChipVariant}, IEnumerable{WbcObservation}, double, double, int)"/>.
    /// Source: Bolton et al. (2020), <i>Nat Genet</i> 52(11):1219–1226.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier of the locus.</param>
    /// <param name="Position">1-based reference position.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Alternate (mutant) allele.</param>
    /// <param name="Vaf">WBC variant allele fraction in [0, 1].</param>
    /// <param name="AltReads">Alternate (mutant) supporting reads in the WBC sample (≥ 0).</param>
    public readonly record struct WbcObservation(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        double Vaf,
        int AltReads = 0);

    /// <summary>
    /// The called origin of a single tumour/plasma <see cref="ChipVariant"/> together with the matched-WBC
    /// evidence the call was based on.
    /// </summary>
    /// <param name="Variant">The tumour/plasma variant whose origin was called.</param>
    /// <param name="Origin">The called origin (<see cref="VariantOrigin.Chip"/> or <see cref="VariantOrigin.Tumor"/>).</param>
    /// <param name="WbcVaf">The matched-WBC VAF at this locus, or <c>0</c> when the locus is absent from the WBC sample.</param>
    /// <param name="WbcAltReads">The matched-WBC supporting alt reads at this locus, or <c>0</c> when absent.</param>
    public readonly record struct VariantOriginCall(
        ChipVariant Variant,
        VariantOrigin Origin,
        double WbcVaf,
        int WbcAltReads);

    /// <summary>
    /// Performs <b>strict matched-WBC origin calling</b>: given per-variant matched white-blood-cell (WBC)
    /// observations, assigns each tumour/plasma variant an origin instead of using the gene + VAF heuristic.
    /// A variant is called <see cref="VariantOrigin.Chip"/> (white-blood-cell / clonal-hematopoiesis derived)
    /// when a matched-WBC observation exists at the same locus that ALL hold: its WBC VAF is at or above
    /// <paramref name="chipMinWbcVaf"/> (Bolton et al. 2020: ≥ 2%), it carries at least
    /// <paramref name="minWbcAltReads"/> supporting reads (Bolton et al. 2020: ≥ 10), and its WBC VAF is at
    /// least <paramref name="wbcVafFold"/> times the variant's tumour/plasma VAF (Bolton et al. 2020: ≥ 2×,
    /// or use <see cref="LymphNodeWbcVafFold"/> = 1.5× for a lymph-node biopsy). Otherwise the variant is
    /// called <see cref="VariantOrigin.Tumor"/> (tumour-derived / somatic). This is the definitive
    /// origin test (Razavi et al. 2019: matched cfDNA–WBC sequencing assigns variant origin tumour-vs-CH);
    /// it does NOT apply the <see cref="IdentifyCHIPVariants"/> gene + VAF fallback, so it does not
    /// over-remove driver-gene variants that are genuinely absent from the matched WBC.
    /// </summary>
    /// <param name="variants">Tumour/plasma variants whose origin is to be called (non-null).</param>
    /// <param name="whiteBloodCellObservations">Matched-WBC observations carrying per-locus WBC VAF and alt reads (non-null).</param>
    /// <param name="wbcVafFold">
    /// Minimum WBC-to-tumour VAF ratio for a WBC call (default <see cref="DefaultWbcVafFold"/> = 2.0; pass
    /// <see cref="LymphNodeWbcVafFold"/> = 1.5 for a lymph-node biopsy site); must be ≥ 1.
    /// </param>
    /// <param name="chipMinWbcVaf">Minimum WBC VAF for a WBC call (default <see cref="ChipVafThreshold"/> = 0.02); in (0, 1].</param>
    /// <param name="minWbcAltReads">Minimum WBC supporting alt reads for a WBC call (default <see cref="ChipMinWbcSupportingReads"/> = 10); ≥ 1.</param>
    /// <returns>One <see cref="VariantOriginCall"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> or <paramref name="whiteBloodCellObservations"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="wbcVafFold"/> &lt; 1, <paramref name="chipMinWbcVaf"/> ∉ (0, 1], or <paramref name="minWbcAltReads"/> &lt; 1.</exception>
    public static IReadOnlyList<VariantOriginCall> CallVariantOrigin(
        IEnumerable<ChipVariant> variants,
        IEnumerable<WbcObservation> whiteBloodCellObservations,
        double wbcVafFold = DefaultWbcVafFold,
        double chipMinWbcVaf = ChipVafThreshold,
        int minWbcAltReads = ChipMinWbcSupportingReads)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ArgumentNullException.ThrowIfNull(whiteBloodCellObservations);

        if (!(wbcVafFold >= 1.0))
        {
            throw new ArgumentOutOfRangeException(
                nameof(wbcVafFold), wbcVafFold, "WBC-to-tumour VAF fold ratio must be at least 1.");
        }

        if (!(chipMinWbcVaf > 0.0) || chipMinWbcVaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(chipMinWbcVaf), chipMinWbcVaf, "Minimum WBC VAF must be in the interval (0, 1].");
        }

        if (minWbcAltReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minWbcAltReads), minWbcAltReads, "Minimum WBC alt reads must be at least 1.");
        }

        // Index the matched-WBC observations by locus. If several observations map to the same locus, keep
        // the one with the strongest evidence (highest VAF) so a confident WBC call is not missed.
        var wbcByLocus = new Dictionary<(string, int, string, string), WbcObservation>();
        foreach (WbcObservation obs in whiteBloodCellObservations)
        {
            (string, int, string, string) key = LocusKey(obs);
            if (!wbcByLocus.TryGetValue(key, out WbcObservation existing) || obs.Vaf > existing.Vaf)
            {
                wbcByLocus[key] = obs;
            }
        }

        var calls = new List<VariantOriginCall>();
        foreach (ChipVariant variant in variants)
        {
            double wbcVaf = 0.0;
            int wbcAltReads = 0;
            VariantOrigin origin = VariantOrigin.Tumor;

            if (wbcByLocus.TryGetValue(LocusKey(variant), out WbcObservation wbc))
            {
                wbcVaf = wbc.Vaf;
                wbcAltReads = wbc.AltReads;

                // Bolton et al. (2020): WBC/CH origin requires WBC VAF >= 2%, >= 10 supporting reads, and a
                // WBC VAF at least (fold) x the tumour VAF.
                bool meetsVaf = wbc.Vaf >= chipMinWbcVaf;
                bool meetsReads = wbc.AltReads >= minWbcAltReads;
                bool meetsFold = wbc.Vaf >= wbcVafFold * variant.Vaf;
                if (meetsVaf && meetsReads && meetsFold)
                {
                    origin = VariantOrigin.Chip;
                }
            }

            calls.Add(new VariantOriginCall(variant, origin, wbcVaf, wbcAltReads));
        }

        return calls;
    }

    #endregion

}
