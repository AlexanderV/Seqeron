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

        // Posterior P(c) ∝ Binomial(a | N, f(c)) on Landau's regular grid of 100 values c = 0.01, 0.02, …, 1.00,
        // uniform prior, normalised by the sum. The binomial coefficient C(N, a) is constant in c, so it cancels in
        // normalisation and is omitted. The kernel is kept in log space and shifted by its maximum before
        // exponentiation (log-sum-exp): without the C(N, a) factor the raw kernel p^a(1−p)^(N−a) underflows to 0 at
        // every grid point once N ≳ 1100 (e.g. a = 1000, N = 2000), which previously collapsed the posterior to a
        // flat grid (CCF 0.505, subclonal) instead of R dbinom's normalised posterior (CCF 0.985, clonal).
        Span<double> weights = stackalloc double[CcfGridPointCount]; // log weights, then max-shifted weights
        double maxLogWeight = double.NegativeInfinity;
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double f = Math.Min(1.0, alleleFractionPerUnitCcf * CcfGridPoint(i));
            double logLikelihood = BinomialLogLikelihoodKernel(variant.AltReads, variant.TotalReads, f);
            weights[i] = logLikelihood;
            maxLogWeight = Math.Max(maxLogWeight, logLikelihood);
        }

        // f(c) ∈ (0, 1) for every grid point except possibly f(1) = 1 (ρ = 1, M = q), so at least one log weight is
        // finite and the shifted weights sum to ≥ 1: the posterior is always well defined.
        double weightSum = 0.0;
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double weight = Math.Exp(weights[i] - maxLogWeight);
            weights[i] = weight;
            weightSum += weight;
        }

        double ccfMean = 0.0;
        double probabilityClonal = 0.0;
        for (int i = 0; i < CcfGridPointCount; i++)
        {
            double c = CcfGridPoint(i);
            double posterior = weights[i] / weightSum;
            ccfMean += c * posterior;
            // Strict "CCF > 0.95" (Landau 2013). Grid points are exact decimals, so c = 0.95 is excluded.
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
    /// The i-th point (0-based) of Landau's regular CCF grid of <see cref="CcfGridPointCount"/> values over
    /// [<see cref="CcfGridLowerBound"/>, <see cref="CcfGridUpperBound"/>]: c_i = (i + 1) / 100, i.e. 0.01, 0.02, …, 1.00.
    /// Computed as a single correctly-rounded division so every grid point is the double nearest its decimal value;
    /// the accumulated form 0.01 + i·(0.99/99) yields 0.9500000000000001 at i = 94, which made the strict
    /// "CCF &gt; 0.95" test count the c = 0.95 grid point as clonal.
    /// </summary>
    private static double CcfGridPoint(int index) => (double)(index + 1) / CcfGridPointCount;

    /// <summary>
    /// Binomial log-likelihood kernel ln L(a | N, p) = a·ln p + (N−a)·ln(1−p), without the constant ln C(N, a)
    /// (it cancels under grid normalisation). Returns −∞ where the likelihood is exactly zero (p = 0 with a &gt; 0,
    /// p = 1 with a &lt; N).
    /// </summary>
    private static double BinomialLogLikelihoodKernel(int altReads, int totalReads, double p)
    {
        int refReads = totalReads - altReads;
        if (p <= 0.0)
        {
            // p = 0 explains zero alternate reads exactly, nothing else.
            return altReads == 0 ? 0.0 : double.NegativeInfinity;
        }

        if (p >= 1.0)
        {
            // p = 1 explains all-alternate reads exactly, nothing else.
            return refReads == 0 ? 0.0 : double.NegativeInfinity;
        }

        return (altReads * Math.Log(p)) + (refReads * Math.Log(1.0 - p));
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
        // The observed mutation copy number n_mut = VAF·(ρ·N_T + 2(1−ρ))/ρ is the canonical CNAqc purity/copy-number
        // VAF correction (AdjustVAFForPurity, ONCO-VAF-001); CCF = n_mut / m (McGranahan 2016).
        double rawCcf = AdjustVAFForPurity(vaf, purity, tumorCopyNumber) / multiplicity;
        double cappedCcf = Math.Min(MaxCancerCellFraction, rawCcf);
        return new CcfEstimate(cappedCcf, rawCcf);
    }

    /// <summary>
    /// Upper bound on the number of cells (clusters × values) of the Ckmeans.1d.dp backtrack matrix J that
    /// <see cref="ClusterCcfValues"/> allocates (4 bytes per cell, i.e. at most 400 MB).
    /// </summary>
    private const long MaxCcfClusteringDpCells = 100_000_000;

    /// <summary>
    /// Clusters cancer cell fractions into <paramref name="clusterCount"/> clones/subclones by <b>optimal</b>
    /// one-dimensional k-means: the partition minimising the within-cluster sum of squares
    /// Σ_j Σ_{x∈S_j} (x − μ_j)² (the k-means objective, Lloyd 1982), found exactly by the dynamic program of
    /// Wang &amp; Song (2011), <i>The R Journal</i> 3(2):29–33 (R package Ckmeans.1d.dp). This is a line-by-line port
    /// of Ckmeans.1d.dp 4.3.x <c>EWL2::fill_dp_matrix</c> (median-shifted prefix sums, log-linear row fill
    /// <c>fill_row_q_log_linear</c>) and <c>backtrack</c>; unlike Lloyd iterations it cannot stop in a local
    /// optimum and needs no seeding, so the result is deterministic and independent of input order. As in
    /// Ckmeans.1d.dp, when the input has fewer distinct values than <paramref name="clusterCount"/> the number of
    /// clusters is reduced to the number of distinct values (every returned cluster is non-empty). The clonal
    /// cluster is the one with the highest centroid (Tarabichi et al. 2021, <i>Nat. Methods</i> 18:144–155: "the
    /// cluster with the highest CP can be deemed clonal").
    /// </summary>
    /// <param name="ccfValues">Cancer cell fractions to cluster (each finite).</param>
    /// <param name="clusterCount">Number of clusters k, in [1, count of values].</param>
    /// <returns>Ascending centroids (min(k, number of distinct values) of them, each the mean of a non-empty
    /// cluster), per-value cluster assignments (input order), and the clonal cluster index (the last centroid).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="ccfValues"/> is null.</exception>
    /// <exception cref="ArgumentException">no values are supplied, or a value is NaN/infinite.</exception>
    /// <exception cref="ArgumentOutOfRangeException">clusterCount ∉ [1, count], or the dynamic-programming
    /// matrix (effective k × count) would exceed 10⁸ cells.</exception>
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

        // Sort values (carrying original indices; stable) — Ckmeans.1d.dp works on the sorted data.
        int[] order = Enumerable.Range(0, n).OrderBy(i => ccfValues[i]).ToArray();
        double[] x = new double[n];
        for (int s = 0; s < n; s++)
        {
            x[s] = ccfValues[order[s]];
        }

        // Ckmeans.1d.dp: Kmax = min(k, number of unique values).
        int distinct = 1;
        for (int s = 1; s < n; s++)
        {
            if (x[s] != x[s - 1])
            {
                distinct++;
            }
        }

        int k = Math.Min(clusterCount, distinct);
        if (k > 1 && k < n && (long)k * n > MaxCcfClusteringDpCells)
        {
            throw new ArgumentOutOfRangeException(
                nameof(clusterCount), clusterCount,
                $"Optimal 1-D k-means needs a {k} × {n} matrix, above the {MaxCcfClusteringDpCells}-cell limit.");
        }

        int[] clusterStart = CkmeansClusterStarts(x, k);

        // Backtrack (Ckmeans.1d.dp backtrack): centre = arithmetic mean of each contiguous sorted block.
        double[] centroids = new double[k];
        int[] sortedCluster = new int[n];
        for (int q = 0; q < k; q++)
        {
            int left = clusterStart[q];
            int right = q + 1 < k ? clusterStart[q + 1] - 1 : n - 1;
            double sum = 0.0;
            for (int s = left; s <= right; s++)
            {
                sum += x[s];
                sortedCluster[s] = q;
            }

            centroids[q] = sum / (right - left + 1);
        }

        int[] assignments = new int[n];
        for (int s = 0; s < n; s++)
        {
            assignments[order[s]] = sortedCluster[s];
        }

        // Clonal cluster = highest centroid; blocks of sorted data give ascending centroids, so it is the last index.
        return new CcfClustering(centroids, assignments, k - 1);
    }

    /// <summary>
    /// Ckmeans.1d.dp (Wang &amp; Song 2011) optimal 1-D k-means on sorted <paramref name="x"/> with
    /// <paramref name="k"/> ≤ number of distinct values: returns the first sorted index of each of the k clusters.
    /// Port of <c>EWL2::fill_dp_matrix</c> + <c>EWL2::fill_row_q_log_linear</c> + <c>backtrack</c>.
    /// </summary>
    private static int[] CkmeansClusterStarts(double[] x, int k)
    {
        int n = x.Length;
        int[] starts = new int[k];
        if (k == 1)
        {
            return starts;
        }

        if (k == n)
        {
            // Every value distinct and its own cluster: the unique zero-cost optimum.
            for (int q = 0; q < k; q++)
            {
                starts[q] = q;
            }

            return starts;
        }

        // Median-shifted running sums for numerical stability (Ckmeans.1d.dp: shift = x[N/2]).
        double shift = x[n / 2];
        double[] sumX = new double[n];
        double[] sumXSq = new double[n];
        sumX[0] = x[0] - shift;
        sumXSq[0] = (x[0] - shift) * (x[0] - shift);

        // S keeps only rows q−1 and q; J (backtrack) keeps every row.
        double[] sPrev = new double[n];
        double[] sCur = new double[n];
        int[][] j = new int[k][];
        j[0] = new int[n];
        for (int i = 1; i < n; i++)
        {
            sumX[i] = sumX[i - 1] + x[i] - shift;
            sumXSq[i] = sumXSq[i - 1] + (x[i] - shift) * (x[i] - shift);
            sPrev[i] = CkmeansSsq(0, i, sumX, sumXSq);
        }

        for (int q = 1; q < k; q++)
        {
            j[q] = new int[n];
            int imin = q < k - 1 ? Math.Max(1, q) : n - 1;
            CkmeansFillRowLogLinear(imin, n - 1, q, q, n - 1, sPrev, sCur, j[q - 1], j[q], sumX, sumXSq);
            (sPrev, sCur) = (sCur, sPrev);
        }

        int right = n - 1;
        for (int q = k - 1; q >= 0; q--)
        {
            int left = j[q][right];
            starts[q] = left;
            if (q > 0)
            {
                right = left - 1;
            }
        }

        return starts;
    }

    /// <summary>Ckmeans.1d.dp <c>EWL2::fill_row_q_log_linear</c> (divide and conquer over i with monotone J).</summary>
    private static void CkmeansFillRowLogLinear(
        int imin, int imax, int q, int jmin, int jmax,
        double[] sPrev, double[] sCur, int[] jPrev, int[] jCur, double[] sumX, double[] sumXSq)
    {
        if (imin > imax)
        {
            return;
        }

        int n = sCur.Length;
        int i = (imin + imax) / 2;

        sCur[i] = sPrev[i - 1];
        jCur[i] = i;

        int jlow = q;
        if (imin > q)
        {
            jlow = Math.Max(jlow, jmin);
        }

        jlow = Math.Max(jlow, jPrev[i]);

        int jhigh = i - 1;
        if (imax < n - 1)
        {
            jhigh = Math.Min(jhigh, jmax);
        }

        for (int jj = jhigh; jj >= jlow; --jj)
        {
            double sji = CkmeansSsq(jj, i, sumX, sumXSq);
            if (sji + sPrev[jlow - 1] >= sCur[i])
            {
                break;
            }

            double ssqJlow = CkmeansSsq(jlow, i, sumX, sumXSq) + sPrev[jlow - 1];
            if (ssqJlow < sCur[i])
            {
                sCur[i] = ssqJlow;
                jCur[i] = jlow;
            }

            jlow++;

            double ssqJ = sji + sPrev[jj - 1];
            if (ssqJ < sCur[i])
            {
                sCur[i] = ssqJ;
                jCur[i] = jj;
            }
        }

        int leftJmin = imin > q ? jCur[imin - 1] : q;
        CkmeansFillRowLogLinear(imin, i - 1, q, leftJmin, jCur[i], sPrev, sCur, jPrev, jCur, sumX, sumXSq);

        int rightJmax = imax < n - 1 ? jCur[imax + 1] : imax;
        CkmeansFillRowLogLinear(i + 1, imax, q, jCur[i], rightJmax, sPrev, sCur, jPrev, jCur, sumX, sumXSq);
    }

    /// <summary>Ckmeans.1d.dp <c>EWL2::ssq</c>: within-cluster sum of squares of sorted x[j..i] from prefix sums.</summary>
    private static double CkmeansSsq(int j, int i, double[] sumX, double[] sumXSq)
    {
        double sji;
        if (j >= i)
        {
            sji = 0.0;
        }
        else if (j > 0)
        {
            double muji = (sumX[i] - sumX[j - 1]) / (i - j + 1);
            sji = sumXSq[i] - sumXSq[j - 1] - (i - j + 1) * muji * muji;
        }
        else
        {
            sji = sumXSq[i] - sumX[i] * sumX[i] / (i + 1);
        }

        return sji < 0 ? 0 : sji;
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
