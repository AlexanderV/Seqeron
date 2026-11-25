// SPDX-License-Identifier: GPL-3.0-or-later
// Port of CNAqc 1.1.5 (caravagnalab/CNAqc, "License: GPL (>= 3) + file LICENSE"), BMix 1.0.0 (caravagnalab/BMix, "License: GPL (>= 3)") and the R stats kmeans (Hartigan-Wong, MacQueen).
// Part of Seqeron.Genomics.PurityPeakQc, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.PurityPeakQc;

/// <summary>CNAqc peak-based purity QC (port of CNAqc 1.1.5 analyze_peaks).</summary>
public static partial class PurityPeakAnalyzer
{

    /// <summary>
    /// A somatic mutation for CNAqc's peak-based purity QC: its VAF and the clonal allele-specific copy-number state
    /// (Major:minor) of the segment it maps to (CNAqc <c>karyotype = "Major:minor"</c>).
    /// </summary>
    /// <param name="Vaf">Variant allele frequency NV/DP in [0, 1].</param>
    /// <param name="MajorCopyNumber">Major allele copy number (≥ 0).</param>
    /// <param name="MinorCopyNumber">Minor allele copy number (≥ 0).</param>
    public readonly record struct PurityPeakMutation(double Vaf, int MajorCopyNumber, int MinorCopyNumber)
    {
        /// <summary>Alternate-allele read count NV (CNAqc <c>NV</c>); used by the BMix fit (<see cref="PurityPeakOptions.FitMixturePeaks"/>).</summary>
        public int AlternateReads { get; init; }

        /// <summary>
        /// Read depth DP (CNAqc <c>DP</c>; 0 = not supplied). When every analysed mutation has DP ≥ 1, BMix mixture peaks are
        /// fitted by default (<see cref="PurityPeakOptions.FitMixturePeaks"/>).
        /// </summary>
        public int Depth { get; init; }
    }

    /// <summary>How detected VAF peaks are assigned to the expected clonal peaks of a karyotype.</summary>
    public enum PurityPeakMatchingStrategy
    {
        /// <summary>
        /// Each expected peak takes the nearest non-discarded data peak (first on ties) — CNAqc 1.1.5
        /// <c>analyze_peaks_common</c>, the only strategy its <c>analyze_peaks</c> runs.
        /// </summary>
        Closest,

        /// <summary>
        /// Expected peaks (descending) are paired with the highest-VAF non-discarded data peaks (descending), padding
        /// with the rightmost peak when there are fewer data peaks — CNAqc's legacy <c>peak_detector</c>
        /// (<c>matching_strategy = "rightmost"</c>; still documented by <c>analyze_peaks</c> but no longer forwarded to
        /// <c>analyze_peaks_common</c> in 1.1.5).
        /// </summary>
        Rightmost,
    }

    /// <summary>Origin of a data peak.</summary>
    public enum PurityPeakSource
    {
        /// <summary>Gaussian-KDE maximum found by <c>peakPick</c> (CNAqc <c>from = "KDE"</c>).</summary>
        Kde,

        /// <summary>A caller-supplied mixture-component mean (e.g. BMix Binomial <c>B.params</c>; CNAqc <c>from = "BMix"</c>).</summary>
        Mixture,
    }

    /// <summary>
    /// Parameters of <see cref="AnalyzePurityPeaks"/>; defaults are CNAqc 1.1.5 <c>analyze_peaks</c> defaults.
    /// </summary>
    public sealed record PurityPeakOptions
    {
        /// <summary>CNAqc default simple clonal karyotypes <c>c('1:0', '1:1', '2:0', '2:1', '2:2')</c> as (Major, minor).</summary>
        public static IReadOnlyList<(int Major, int Minor)> SimpleClonalKaryotypes { get; } =
            new[] { (1, 0), (1, 1), (2, 0), (2, 1), (2, 2) };

        /// <summary>Default options (CNAqc defaults).</summary>
        public static PurityPeakOptions Default { get; } = new();

        /// <summary>Karyotypes to QC (<c>karyotypes</c>); must be a subset of <see cref="SimpleClonalKaryotypes"/>.</summary>
        public IReadOnlyList<(int Major, int Minor)> Karyotypes { get; init; } = SimpleClonalKaryotypes;

        /// <summary>Minimum share n_k / N of all mutations (<c>min_karyotype_size</c>, default 0), in [0, 1).</summary>
        public double MinKaryotypeSize { get; init; }

        /// <summary>Minimum mutations per karyotype (<c>min_absolute_karyotype_mutations</c>, default 100; n ≥ this).</summary>
        public int MinAbsoluteKaryotypeMutations { get; init; } = 100;

        /// <summary>Purity error ε (<c>purity_error</c>, default 0.05), in (0, 1); sets the VAF bands δ.</summary>
        public double PurityError { get; init; } = 0.05;

        /// <summary>VAF tolerance around a data peak for band overlap (<c>VAF_tolerance</c>, default 0.015), ≥ 0.</summary>
        public double VafTolerance { get; init; } = 0.015;

        /// <summary>KDE bandwidth multiplier (<c>kernel_adjust</c>, default 1), &gt; 0.</summary>
        public double KernelAdjust { get; init; } = 1.0;

        /// <summary>Peak-matching rule (default <see cref="PurityPeakMatchingStrategy.Closest"/>).</summary>
        public PurityPeakMatchingStrategy MatchingStrategy { get; init; } = PurityPeakMatchingStrategy.Closest;

        /// <summary>Only mutations with VAF &gt; this are analysed (<c>min_VAF</c>, default 0).</summary>
        public double MinVaf { get; init; }

        /// <summary>
        /// Optional caller-supplied mixture-component means per karyotype (CNAqc adds BMix Binomial-mixture peaks,
        /// <c>bmixfit(K.Binomials = 1:4)</c>, to the KDE peaks). Each mean is snapped to the nearest KDE grid point.
        /// null = none; to fit them as CNAqc does use <see cref="FitMixturePeaks"/> instead (not both).
        /// </summary>
        public IReadOnlyDictionary<(int Major, int Minor), IReadOnlyList<double>>? MixturePeaks { get; init; }

        /// <summary>true reproduces R ≤ 4.3 <c>density</c> values (<c>old.coords = TRUE</c>); default R ≥ 4.4.</summary>
        public bool LegacyDensityCoordinates { get; init; }

        /// <summary>
        /// Fit CNAqc's BMix Binomial-mixture peaks (<c>mixture_peak_detector</c>: <see cref="FitBinomialMixture"/> on
        /// (NV, DP) of each analysed karyotype, ICL-best component means snapped to the KDE grid). CNAqc's
        /// <c>analyze_peaks</c> always runs them (<c>combined_peak_detector</c> → <c>mixture_peak_detector</c> → <c>bmixfit</c>).
        /// <para>null (default) = automatic, as CNAqc: on when every mutation of the analysed karyotypes carries read counts
        /// (<see cref="PurityPeakMutation.Depth"/> ≥ 1), off when none does (VAF-only input → KDE peaks only) or when
        /// <see cref="MixturePeaks"/> are supplied; a mix of mutations with and without read counts throws
        /// <see cref="ArgumentException"/>. true = always fit (read counts required); false = never (KDE peaks only, plus
        /// any <see cref="MixturePeaks"/>).</para>
        /// </summary>
        public bool? FitMixturePeaks { get; init; }

        /// <summary>
        /// Peak-detection bootstrap replicates (<c>n_bootstrap</c>, default 1 = none), ≥ 1. With n &gt; 1 every analysed
        /// karyotype's VAFs are resampled n times with replacement (<c>sample_n(replace = TRUE)</c>) and the KDE peaks of the
        /// replicates not already found are added (heights re-read from the full-data density); with
        /// <see cref="FitMixturePeaks"/> the BMix fit is bootstrapped likewise.
        /// </summary>
        public int BootstrapCount { get; init; } = 1;

        /// <summary>
        /// R <c>set.seed</c> value for the random steps (bootstrap resampling, BMix k-means starts and jitter); the R
        /// stream is reproduced exactly (Mersenne-Twister, <c>sample.kind = "Rejection"</c>), so a run equals CNAqc's after
        /// <c>set.seed(Seed)</c>. Unused when nothing is random.
        /// </summary>
        public int Seed { get; init; }
    }

    /// <summary>A VAF peak detected in a karyotype's data (CNAqc <c>xy_peaks</c> row).</summary>
    /// <param name="X">Peak VAF (KDE peaks rounded to 2 decimals; mixture peaks at the snapped KDE grid point).</param>
    /// <param name="Y">Density at the peak (KDE peaks rounded to 2 decimals).</param>
    /// <param name="CountsPerBin">Mutations in the 0.01-wide VAF histogram bin <c>round(100·X)</c> (null if outside 1..100).</param>
    /// <param name="Discarded">KDE peak with Y ≤ max(Y)/20 (never used for matching); mixture peaks are never discarded.</param>
    /// <param name="Source">KDE or mixture.</param>
    public readonly record struct PurityDataPeak(double X, double Y, int? CountsPerBin, bool Discarded, PurityPeakSource Source);

    /// <summary>One expected clonal peak and its matched data peak (CNAqc <c>peaks_analysis$matches</c> row).</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="Multiplicity">Mutation multiplicity m (1 or Major).</param>
    /// <param name="ExpectedPeak">m·π / (2(1−π) + π·(Major + minor)).</param>
    /// <param name="DeltaVaf">VAF band half-width δ = 2·m·ε / (2 + π·(ploidy − 2))².</param>
    /// <param name="MatchedPeak">The data peak assigned to this expectation.</param>
    /// <param name="OffsetVaf">ExpectedPeak − MatchedPeak.X.</param>
    /// <param name="Offset">Purity-space offset 2·m·OffsetVaf / (m + X·(2 − ploidy))².</param>
    /// <param name="Weight">Karyotype weight n_k / Σ n over analysed karyotypes.</param>
    /// <param name="Matched">[X ± VafTolerance] overlaps [ExpectedPeak ± δ].</param>
    public readonly record struct PurityPeakMatch(
        int MajorCopyNumber,
        int MinorCopyNumber,
        int Multiplicity,
        double ExpectedPeak,
        double DeltaVaf,
        PurityDataPeak MatchedPeak,
        double OffsetVaf,
        double Offset,
        double Weight,
        bool Matched);

    /// <summary>Per-karyotype result of <see cref="AnalyzePurityPeaks"/>.</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="MutationCount">Mutations pooled for this karyotype (VAF &gt; MinVaf).</param>
    /// <param name="Weight">n_k / Σ n over analysed karyotypes.</param>
    /// <param name="Score">Σ Weight·Offset over this karyotype's expected peaks.</param>
    /// <param name="Pass">QC of the expected peak whose matched data peak holds the most mutations (CountsPerBin).</param>
    /// <param name="Density">The Gaussian KDE of the karyotype's VAFs.</param>
    /// <param name="Peaks">All data peaks (KDE, then mixture).</param>
    /// <param name="Matches">Expected-peak matches (Closest: m = 1 then Major; Rightmost: descending expected VAF).</param>
    public sealed record PurityPeakKaryotype(
        int MajorCopyNumber,
        int MinorCopyNumber,
        int MutationCount,
        double Weight,
        double Score,
        bool Pass,
        KernelDensityEstimate Density,
        IReadOnlyList<PurityDataPeak> Peaks,
        IReadOnlyList<PurityPeakMatch> Matches);

    /// <summary>Result of <see cref="AnalyzePurityPeaks"/> (CNAqc <c>x$peaks_analysis</c> for simple clonal CNAs).</summary>
    /// <param name="Purity">The purity being QC'd.</param>
    /// <param name="Score">Σ Weight·Offset over all matches (CNAqc λ, printed as "Purity correction"); NaN if nothing analysed.</param>
    /// <param name="Pass">Sample QC: the PASS/FAIL class with the larger summed match weight (ties → FAIL); null if no karyotype passed the filters.</param>
    /// <param name="Karyotypes">Analysed karyotypes in (Major, minor) order.</param>
    /// <param name="Matches">All matches, karyotype by karyotype.</param>
    public sealed record PurityPeakAnalysis(
        double Purity,
        double Score,
        bool? Pass,
        IReadOnlyList<PurityPeakKaryotype> Karyotypes,
        IReadOnlyList<PurityPeakMatch> Matches);

    /// <summary>
    /// CNAqc peak-based QC of a tumour purity estimate for simple clonal karyotypes (CNAqc 1.1.5
    /// <c>analyze_peaks</c> → <c>analyze_peaks_common</c>; Antonello et al. 2024, <i>Genome Biology</i> 25:38).
    /// For each karyotype K = Major:minor in <see cref="PurityPeakOptions.Karyotypes"/> with n_K ≥ MinAbsoluteKaryotypeMutations
    /// and n_K / N ≥ MinKaryotypeSize (N = all mutations with VAF &gt; MinVaf, any karyotype):
    /// <list type="number">
    /// <item>pool its VAFs; Gaussian KDE (R <c>density</c>, bw.nrd0 × KernelAdjust); peaks = union of
    /// <c>peakPick::peakpick(neighlim = 1..5)</c> maxima, (x, y) rounded to 2 decimals, distinct x, clamped to [0, 1];
    /// a KDE peak is discarded if y ≤ max(y)/20; optional mixture peaks are appended;</item>
    /// <item>expected peaks for m ∈ {1, Major}: v_m = m·π / (2(1−π) + π·(Major + minor)), bands
    /// δ_m = 2·m·ε / (2 + π·(ploidy − 2))²;</item>
    /// <item>match each expected peak to a non-discarded data peak (<see cref="PurityPeakMatchingStrategy"/>);
    /// offset = 2·m·(v_m − x) / (m + x·(2 − ploidy))² (purity units), weight = n_K / Σ n_analysed;
    /// matched ⇔ [x ± VafTolerance] ∩ [v_m ± δ_m] ≠ ∅;</item>
    /// <item>karyotype PASS ⇔ the expected peak whose matched data peak has the largest histogram count is matched;
    /// sample PASS ⇔ PASS rows carry more total weight than FAIL rows; score λ = Σ weight·offset.</item>
    /// </list>
    /// CNAqc does not threshold λ against ε (ε only sets the bands) and proposes no corrected purity: <c>print</c>
    /// reports λ as "Purity correction". <c>p_binsize_peaks</c> is accepted by CNAqc 1.1.5 but unused. KDE peak
    /// detection is deterministic and ported exactly; CNAqc's BMix mixture peaks are fitted by default whenever the
    /// analysed mutations carry read counts (<see cref="PurityPeakOptions.FitMixturePeaks"/> null = auto;
    /// <see cref="FitBinomialMixture"/>, R random stream reproduced from <see cref="PurityPeakOptions.Seed"/>) or supplied
    /// via <see cref="PurityPeakOptions.MixturePeaks"/>;
    /// <c>n_bootstrap</c> is <see cref="PurityPeakOptions.BootstrapCount"/>. Complex and subclonal karyotypes:
    /// <see cref="AnalyzeComplexKaryotypePeaks"/>, <see cref="AnalyzeSubclonalPurityPeaks"/>.
    /// </summary>
    /// <param name="mutations">Mutations with VAF and karyotype (all karyotypes; non-simple ones count towards N).</param>
    /// <param name="purity">The purity π ∈ (0, 1] to QC.</param>
    /// <param name="options">CNAqc parameters (null = defaults).</param>
    /// <returns>Per-karyotype peaks, matches, scores and the QC verdict.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid purity, option, VAF or copy number.</exception>
    /// <exception cref="ArgumentException">A requested karyotype is not a simple clonal karyotype.</exception>
    /// <exception cref="InvalidOperationException">An analysed karyotype yields no usable VAF peak (CNAqc errors too).</exception>
    public static PurityPeakAnalysis AnalyzePurityPeaks(
        IEnumerable<PurityPeakMutation> mutations,
        double purity,
        PurityPeakOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(mutations);
        options ??= PurityPeakOptions.Default;
        ValidatePurityPeakArguments(purity, options);

        var byKaryotype = new SortedDictionary<(int Major, int Minor), List<PurityPeakMutation>>();
        int total = 0;
        foreach (PurityPeakMutation mutation in mutations)
        {
            if (double.IsNaN(mutation.Vaf) || mutation.Vaf < 0.0 || mutation.Vaf > 1.0)
                throw new ArgumentOutOfRangeException(nameof(mutations), mutation.Vaf, "VAF must be in [0, 1].");
            if (mutation.MajorCopyNumber < 0 || mutation.MinorCopyNumber < 0)
                throw new ArgumentOutOfRangeException(nameof(mutations), "Allele copy numbers must be non-negative.");
            if (!(mutation.Vaf > options.MinVaf)) continue; // CNAqc: filter(VAF > min_VAF)

            var key = (mutation.MajorCopyNumber, mutation.MinorCopyNumber);
            if (!byKaryotype.TryGetValue(key, out List<PurityPeakMutation>? group))
            {
                group = new List<PurityPeakMutation>();
                byKaryotype[key] = group;
            }

            group.Add(mutation);
            total++;
        }

        var analysed = new List<(int Major, int Minor)>();
        int analysedTotal = 0;
        foreach (KeyValuePair<(int Major, int Minor), List<PurityPeakMutation>> entry in byKaryotype)
        {
            int n = entry.Value.Count;
            if (options.Karyotypes.Contains(entry.Key)
                && n >= options.MinAbsoluteKaryotypeMutations
                && (double)n / total >= options.MinKaryotypeSize)
            {
                analysed.Add(entry.Key);
                analysedTotal += n;
            }
        }

        if (analysed.Count == 0)
        {
            return new PurityPeakAnalysis(purity, double.NaN, null, Array.Empty<PurityPeakKaryotype>(), Array.Empty<PurityPeakMatch>());
        }

        bool fitMixture = ResolveFitMixturePeaks(options, analysed.SelectMany(k => byKaryotype[k]), nameof(mutations));

        // One R random stream for the whole call (set.seed before analyze_peaks_common), consumed karyotype by karyotype.
        RMersenneTwister? rng = fitMixture || options.BootstrapCount > 1 ? new RMersenneTwister(options.Seed) : null;
        var karyotypes = new List<PurityPeakKaryotype>(analysed.Count);
        var allMatches = new List<PurityPeakMatch>();
        foreach ((int major, int minor) in analysed)
        {
            List<PurityPeakMutation> group = byKaryotype[(major, minor)];
            double weight = (double)group.Count / analysedTotal;
            IReadOnlyList<double>? mixturePeaks = null;
            if (options.MixturePeaks is not null && options.MixturePeaks.TryGetValue((major, minor), out IReadOnlyList<double>? m))
                mixturePeaks = m;

            PurityPeakKaryotype result = AnalyzePurityPeaksKaryotype(major, minor, group, weight, purity, options, fitMixture, mixturePeaks, rng);
            karyotypes.Add(result);
            allMatches.AddRange(result.Matches);
        }

        double score = 0.0, passWeight = 0.0, failWeight = 0.0;
        foreach (PurityPeakKaryotype k in karyotypes)
        {
            foreach (PurityPeakMatch match in k.Matches)
            {
                score += match.Weight * match.Offset;
                if (k.Pass) passWeight += match.Weight; else failWeight += match.Weight;
            }
        }

        // dplyr group_by(QC) orders "FAIL" < "PASS"; arrange(desc(prop)) is stable, so a tie resolves to FAIL.
        return new PurityPeakAnalysis(purity, score, passWeight > failWeight, karyotypes, allMatches);
    }

    // PurityPeakOptions.FitMixturePeaks: explicit true/false as given; null (auto) = CNAqc's always-on BMix whenever every
    // analysed mutation carries read counts (DP ≥ 1), off for VAF-only input or caller-supplied MixturePeaks. Mixed input
    // (some DP = 0, some ≥ 1) is rejected: CNAqc requires NV/DP on every mutation, and silently dropping BMix would hide
    // partially missing read counts.
    private static bool ResolveFitMixturePeaks(PurityPeakOptions options, IEnumerable<PurityPeakMutation> analysed, string paramName)
    {
        if (options.FitMixturePeaks == false || (options.FitMixturePeaks is null && options.MixturePeaks is not null))
            return false;

        int withReads = 0, withoutReads = 0;
        foreach (PurityPeakMutation m in analysed)
        {
            if (m.Depth >= 1) withReads++; else withoutReads++;
            if (m.Depth >= 1 && (m.AlternateReads < 0 || m.AlternateReads > m.Depth))
                throw new ArgumentException("BMix mixture peaks need 0 ≤ AlternateReads ≤ Depth on every analysed mutation.", paramName);
        }

        if (options.FitMixturePeaks == true)
        {
            if (withoutReads > 0)
                throw new ArgumentException(
                    "FitMixturePeaks needs read counts: every analysed mutation must have Depth ≥ 1 and 0 ≤ AlternateReads ≤ Depth.",
                    paramName);
            return true;
        }

        if (withReads > 0 && withoutReads > 0)
            throw new ArgumentException(
                $"Read counts supplied for {withReads} of {withReads + withoutReads} analysed mutations: supply Depth ≥ 1 for all " +
                "(BMix mixture peaks, CNAqc default) or for none (KDE peaks only), or set FitMixturePeaks explicitly.",
                paramName);
        return withReads > 0;
    }

    private static void ValidatePurityPeakArguments(double purity, PurityPeakOptions options)
    {
        if (!(purity > 0.0 && purity <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        if (!(options.PurityError > 0.0 && options.PurityError < 1.0))
            throw new ArgumentOutOfRangeException(nameof(options), options.PurityError, "PurityError must be in (0, 1).");
        if (!(options.MinKaryotypeSize >= 0.0 && options.MinKaryotypeSize < 1.0))
            throw new ArgumentOutOfRangeException(nameof(options), options.MinKaryotypeSize, "MinKaryotypeSize must be in [0, 1).");
        if (options.MinAbsoluteKaryotypeMutations < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.MinAbsoluteKaryotypeMutations, "MinAbsoluteKaryotypeMutations must be ≥ 0.");
        if (!(options.VafTolerance >= 0.0) || double.IsPositiveInfinity(options.VafTolerance))
            throw new ArgumentOutOfRangeException(nameof(options), options.VafTolerance, "VafTolerance must be finite and ≥ 0.");
        if (!(options.KernelAdjust > 0.0) || double.IsPositiveInfinity(options.KernelAdjust))
            throw new ArgumentOutOfRangeException(nameof(options), options.KernelAdjust, "KernelAdjust must be finite and > 0.");
        if (double.IsNaN(options.MinVaf))
            throw new ArgumentOutOfRangeException(nameof(options), options.MinVaf, "MinVaf must not be NaN.");
        if (options.BootstrapCount < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.BootstrapCount, "BootstrapCount must be ≥ 1.");
        if (options.FitMixturePeaks == true && options.MixturePeaks is not null)
            throw new ArgumentException("Use either FitMixturePeaks or caller-supplied MixturePeaks, not both.", nameof(options));
        if (options.Karyotypes is null)
            throw new ArgumentException("Karyotypes must not be null.", nameof(options));
        foreach ((int Major, int Minor) k in options.Karyotypes)
        {
            if (!PurityPeakOptions.SimpleClonalKaryotypes.Contains(k))
                throw new ArgumentException(
                    $"Karyotype {k.Major}:{k.Minor} is not a simple clonal karyotype (1:0, 1:1, 2:0, 2:1, 2:2).", nameof(options));
        }
    }

    private static PurityPeakKaryotype AnalyzePurityPeaksKaryotype(
        int major, int minor, List<PurityPeakMutation> group, double weight, double purity, PurityPeakOptions options,
        bool fitMixture, IReadOnlyList<double>? mixturePeaks, RMersenneTwister? rng)
    {
        int ploidy = major + minor;
        var vafs = group.Select(g => g.Vaf).ToList();
        // combined_peak_detector: KDE peaks (simple_peak_detector, with bootstrap), then BMix peaks (mixture_peak_detector).
        (KernelDensityEstimate density, List<PurityDataPeak> peaks) = SimplePeakDetector(
            vafs, options.KernelAdjust, options.LegacyDensityCoordinates, $"karyotype {major}:{minor}", options.BootstrapCount, rng);
        int[] histogram = VafHistogram(vafs);
        if (fitMixture)
            peaks.AddRange(MixturePeakDetector(group, options.KernelAdjust, options.LegacyDensityCoordinates, options.BootstrapCount, rng!));

        if (mixturePeaks is not null)
        {
            foreach (double p in mixturePeaks)
            {
                // CNAqc mixture_peak_detector: w_den = which.min(abs(den$x − p)); x, y from the grid (unrounded).
                int w = 0;
                for (int i = 1; i < density.X.Count; i++)
                {
                    if (Math.Abs(density.X[i] - p) < Math.Abs(density.X[w] - p)) w = i;
                }

                double x = density.X[w];
                peaks.Add(new PurityDataPeak(x, density.Y[w], HistogramCount(histogram, Math.Round(x * 100, MidpointRounding.ToEven)), false, PurityPeakSource.Mixture));
            }
        }

        var candidates = peaks.Where(p => !p.Discarded).ToList();
        if (candidates.Count == 0)
            throw new InvalidOperationException($"No non-discarded VAF peak for karyotype {major}:{minor}.");

        // expected_vaf_peak: multiplicities unique(c(1, Major)).
        var multiplicities = major == 1 ? new[] { 1 } : new[] { 1, major };
        var expected = multiplicities
            .Select(m => (Multiplicity: m, Peak: CnaqcExpectedVaf(m, purity, ploidy)))
            .ToList();

        var assignment = new List<((int Multiplicity, double Peak) Expectation, PurityDataPeak Peak)>();
        if (options.MatchingStrategy == PurityPeakMatchingStrategy.Closest)
        {
            foreach (var e in expected)
            {
                PurityDataPeak best = candidates[0];
                foreach (PurityDataPeak c in candidates)
                {
                    if (Math.Abs(c.X - e.Peak) < Math.Abs(best.X - e.Peak)) best = c;
                }

                assignment.Add((e, best));
            }
        }
        else
        {
            var expectedDesc = expected.OrderByDescending(e => e.Peak).ToList(); // stable
            var peaksDesc = candidates.OrderByDescending(p => p.X).ToList();      // stable
            for (int i = 0; i < expectedDesc.Count; i++)
            {
                assignment.Add((expectedDesc[i], i < peaksDesc.Count ? peaksDesc[i] : peaksDesc[0]));
            }
        }

        var matches = new List<PurityPeakMatch>(assignment.Count);
        double karyotypeScore = 0.0;
        foreach (((int m, double peak), PurityDataPeak data) in assignment)
        {
            double band = 2 * m * options.PurityError;
            double spread = 2 + (purity * (ploidy - 2));
            double deltaVaf = band / (spread * spread);
            double offsetVaf = peak - data.X;
            double denom = m + (data.X * (2 - ploidy));
            double offset = 2 * m * offsetVaf / (denom * denom);
            bool matched = Math.Max(data.X - options.VafTolerance, peak - deltaVaf)
                           <= Math.Min(data.X + options.VafTolerance, peak + deltaVaf);
            matches.Add(new PurityPeakMatch(major, minor, m, peak, deltaVaf, data, offsetVaf, offset, weight, matched));
            karyotypeScore += weight * offset;
        }

        // QC per karyotype: arrange(desc(counts_per_bin)) (stable, NA last) → first row's `matched`.
        PurityPeakMatch lead = matches[0];
        foreach (PurityPeakMatch match in matches)
        {
            if ((match.MatchedPeak.CountsPerBin ?? int.MinValue) > (lead.MatchedPeak.CountsPerBin ?? int.MinValue)) lead = match;
        }

        return new PurityPeakKaryotype(major, minor, vafs.Count, weight, karyotypeScore, lead.Matched, density, peaks, matches);
    }

    // CNAqc simple_peak_detector: union of peakpick(neighlim = 1..5) on the density's y; arrange(x); round(x, 2),
    // round(y, 2); distinct(x) keeping the first; x ∈ (1, 1.01) → 1, x ∈ (−0.01, 0) → 0; keep 0 ≤ x ≤ 1; counts from
    // hist(breaks = seq(0, 1, 0.01)) at round(100·x) (first index 0 → 1); discarded ⇔ y ≤ max(y)·(1/20).
    private static List<PurityDataPeak> DetectKdePeaks(KernelDensityEstimate density, int[] histogram)
    {
        var picked = new SortedSet<int>();
        for (int neighlim = 1; neighlim <= 5; neighlim++)
        {
            bool[] mask = StatisticsHelper.PeakPick(density.Y, neighlim);
            for (int i = 0; i < mask.Length; i++)
            {
                if (mask[i]) picked.Add(i);
            }
        }

        var rows = new List<(double X, double Y)>();
        var seen = new HashSet<double>();
        foreach (int i in picked) // grid x increases with the index, so this is arrange(x)
        {
            double x = RRound(density.X[i], 2);
            if (!seen.Add(x)) continue;
            double y = RRound(density.Y[i], 2);
            if (x > 1 && x < 1.01) x = 1;
            else if (x < 0 && x > -0.01) x = 0;
            if (x <= 1 && x >= 0) rows.Add((x, y));
        }

        var peaks = new List<PurityDataPeak>(rows.Count);
        if (rows.Count == 0) return peaks;
        double maxY = rows.Max(r => r.Y);
        for (int k = 0; k < rows.Count; k++)
        {
            double index = Math.Round(rows[k].X * 100, MidpointRounding.ToEven);
            if (k == 0 && index == 0) index = 1;
            peaks.Add(new PurityDataPeak(rows[k].X, rows[k].Y, HistogramCount(histogram, index),
                rows[k].Y <= maxY * (1.0 / 20), PurityPeakSource.Kde));
        }

        return peaks;
    }

    // R hist(VAF, breaks = seq(0, 1, 0.01), plot = FALSE)$counts: right-closed bins, include.lowest, breaks fuzzed by
    // 1e-7·median(diff(breaks)) (first break down, the others up) before R's C_BinCount bisection.
    private static int[] VafHistogram(IReadOnlyList<double> vafs)
    {
        const int Bins = 100;
        var breaks = new double[Bins + 1];
        for (int k = 0; k <= Bins; k++) breaks[k] = Math.Min(k * 0.01, 1.0); // seq(0, 1, 0.01) = from + (0:n)·by, pmin(to)
        var widths = new double[Bins];
        for (int k = 0; k < Bins; k++) widths[k] = breaks[k + 1] - breaks[k];
        double diddle = 1e-7 * StatisticsHelper.Median(widths);
        var fuzzy = new double[Bins + 1];
        fuzzy[0] = breaks[0] - diddle;
        for (int k = 1; k <= Bins; k++) fuzzy[k] = breaks[k] + diddle;

        var counts = new int[Bins];
        foreach (double v in vafs)
        {
            if (!(fuzzy[0] <= v && v <= fuzzy[Bins])) continue;
            int lo = 0, hi = Bins;
            while (hi - lo >= 2)
            {
                int mid = (hi + lo) / 2;
                if (v > fuzzy[mid]) lo = mid; else hi = mid;
            }

            counts[lo]++;
        }

        return counts;
    }

    // R 1-based hst[index]; indices outside 1..100 give NA (null).
    private static int? HistogramCount(int[] histogram, double index) =>
        index >= 1 && index <= histogram.Length ? histogram[(int)index - 1] : null;

    // R ≥ 4.0 round(x, digits) (src/nmath/fround.c): pick the closer of floor/ceil(x·10^d)/10^d, even on ties.
    private static double RRound(double x, int digits)
    {
        if (double.IsNaN(x) || double.IsInfinity(x) || x == 0.0) return x;
        double sgn = 1.0;
        if (x < 0)
        {
            sgn = -1.0;
            x = -x;
        }

        double l10x = 0.30102999566398119521 * (0.5 + Math.Floor(Math.Log2(x))); // M_LOG10_2·(0.5 + logb(x))
        if (l10x + digits > 15) return sgn * x; // DBL_DIG
        double pow10 = Math.Pow(10, digits);
        double x10 = x * pow10, i10 = Math.Floor(x10);
        double xd = i10 / pow10, xu = Math.Ceiling(x10) / pow10;
        double du = xu - x, dd = x - xd;
        return sgn * ((du < dd || (i10 % 2 == 1 && du == dd)) ? xu : xd);
    }

    // ---- CNAqc analyze_peaks_general: complex clonal karyotypes (FIN-B24 F62) ----

    /// <summary>One expected clonal VAF peak of a complex karyotype (CNAqc <c>peaks_analysis$general$expected_peaks</c> row).</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="Multiplicity">Mutation multiplicity m ∈ 1..max(Major, minor, 1).</param>
    /// <param name="ExpectedPeak">m·π / (2(1−π) + π·(Major + minor)) (CNAqc <c>expected_vaf_fun</c>).</param>
    /// <param name="Matched">Some data peak (discarded ones included) lies strictly within ε of the expected peak.</param>
    public readonly record struct ComplexKaryotypeExpectedPeak(
        int MajorCopyNumber, int MinorCopyNumber, int Multiplicity, double ExpectedPeak, bool Matched);

    /// <summary>Per-karyotype result of <see cref="AnalyzeComplexKaryotypePeaks"/> (CNAqc <c>general$summary</c> row + its fit).</summary>
    /// <param name="MajorCopyNumber">Karyotype Major.</param>
    /// <param name="MinorCopyNumber">Karyotype minor.</param>
    /// <param name="MutationCount">Mutations of this karyotype (VAF &gt; MinVaf) — CNAqc <c>n</c>.</param>
    /// <param name="MatchedPeaks">Expected peaks with a data peak within ε.</param>
    /// <param name="MismatchedPeaks">Expected peaks without one.</param>
    /// <param name="MatchedProportion">MatchedPeaks / (MatchedPeaks + MismatchedPeaks) — CNAqc <c>prop</c>.</param>
    /// <param name="Pass">CNAqc <c>analyze_peaks</c> segment/mutation <c>QC_PASS</c>: <c>prop ≥ 0.5</c>.</param>
    /// <param name="Density">Gaussian KDE of the karyotype's VAFs.</param>
    /// <param name="Peaks">KDE data peaks (same detector as the simple karyotypes; no mixture peaks).</param>
    /// <param name="ExpectedPeaks">Expected peaks, multiplicity ascending.</param>
    public sealed record ComplexKaryotypePeakResult(
        int MajorCopyNumber,
        int MinorCopyNumber,
        int MutationCount,
        int MatchedPeaks,
        int MismatchedPeaks,
        double MatchedProportion,
        bool Pass,
        KernelDensityEstimate Density,
        IReadOnlyList<PurityDataPeak> Peaks,
        IReadOnlyList<ComplexKaryotypeExpectedPeak> ExpectedPeaks);

    /// <summary>Result of <see cref="AnalyzeComplexKaryotypePeaks"/> (CNAqc <c>x$peaks_analysis$general</c>).</summary>
    /// <param name="Purity">The purity being QC'd.</param>
    /// <param name="Ran">false when CNAqc's <c>analyze_peaks</c> gate skips the analysis (no karyotype outside
    /// <see cref="PurityPeakOptions.Karyotypes"/> has more than MinAbsoluteKaryotypeMutations mutations).</param>
    /// <param name="Karyotypes">Analysed karyotypes in CNAqc summary order: descending matched proportion, ties by
    /// the "Major:minor" string (C collation).</param>
    /// <param name="ExpectedPeaks">All expected peaks, karyotypes in first-appearance order of the input mutations.</param>
    public sealed record ComplexKaryotypePeakAnalysis(
        double Purity,
        bool Ran,
        IReadOnlyList<ComplexKaryotypePeakResult> Karyotypes,
        IReadOnlyList<ComplexKaryotypeExpectedPeak> ExpectedPeaks);

    /// <summary>
    /// CNAqc peak QC of complex clonal karyotypes (CNAqc 1.1.5 <c>analyze_peaks</c> → <c>analyze_peaks_general</c>,
    /// called with <c>n_min = min_absolute_karyotype_mutations</c>, <c>epsilon = purity_error</c>, <c>kernel_adjust</c>).
    /// <list type="number">
    /// <item>Gate (<c>analyze_peaks</c>): run only if some karyotype not in <see cref="PurityPeakOptions.Karyotypes"/> has
    /// n &gt; MinAbsoluteKaryotypeMutations (strict).</item>
    /// <item>Analysed karyotypes: every karyotype other than the five simple ones (1:0, 1:1, 2:0, 2:1, 2:2 — fixed, not
    /// the option) with n ≥ MinAbsoluteKaryotypeMutations, in first-appearance order of the mutations (VAF &gt; MinVaf).</item>
    /// <item>Per karyotype: KDE peaks exactly as for simple karyotypes (no BMix); expected peaks for m = 1..max(Major,
    /// minor, 1): m·π / (2(1−π) + π·ploidy); an expected peak is matched iff |x − peak| &lt; ε for some data peak
    /// (discarded peaks included).</item>
    /// <item>Summary: matched / mismatched counts, prop = matched / total; <c>analyze_peaks</c> sets QC PASS iff
    /// prop ≥ 0.5. No score and no sample-level verdict (CNAqc: "a complex segment with many matched peaks is likely to
    /// be correct").</item>
    /// </list>
    /// </summary>
    /// <param name="mutations">Mutations with VAF and karyotype.</param>
    /// <param name="purity">The purity π ∈ (0, 1] to QC.</param>
    /// <param name="options">CNAqc parameters (null = defaults; uses Karyotypes, MinAbsoluteKaryotypeMutations, PurityError,
    /// KernelAdjust, MinVaf, LegacyDensityCoordinates).</param>
    /// <returns>Per-karyotype matches and proportions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid purity, option, VAF or copy number.</exception>
    /// <exception cref="InvalidOperationException">The gate opens but no complex karyotype has n ≥ MinAbsoluteKaryotypeMutations,
    /// or a karyotype yields no KDE peak (CNAqc stops with an R error in both cases).</exception>
    public static ComplexKaryotypePeakAnalysis AnalyzeComplexKaryotypePeaks(
        IEnumerable<PurityPeakMutation> mutations,
        double purity,
        PurityPeakOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(mutations);
        options ??= PurityPeakOptions.Default;
        ValidatePurityPeakArguments(purity, options);

        // x$mutations after filter(VAF > min_VAF); karyotypes kept in first-appearance order (unique / setdiff).
        var order = new List<(int Major, int Minor)>();
        var byKaryotype = new Dictionary<(int Major, int Minor), List<double>>();
        foreach (PurityPeakMutation mutation in mutations)
        {
            ValidatePurityPeakMutation(mutation);
            if (!(mutation.Vaf > options.MinVaf)) continue;
            var key = (mutation.MajorCopyNumber, mutation.MinorCopyNumber);
            if (!byKaryotype.TryGetValue(key, out List<double>? vafs))
            {
                vafs = new List<double>();
                byKaryotype[key] = vafs;
                order.Add(key);
            }

            vafs.Add(mutation.Vaf);
        }

        bool gate = order.Any(k => !options.Karyotypes.Contains(k) && byKaryotype[k].Count > options.MinAbsoluteKaryotypeMutations);
        if (!gate)
        {
            return new ComplexKaryotypePeakAnalysis(purity, false, Array.Empty<ComplexKaryotypePeakResult>(), Array.Empty<ComplexKaryotypeExpectedPeak>());
        }

        var analysis = order
            .Where(k => !PurityPeakOptions.SimpleClonalKaryotypes.Contains(k) && byKaryotype[k].Count >= options.MinAbsoluteKaryotypeMutations)
            .ToList();
        if (analysis.Count == 0)
        {
            // R: expected_peaks = Reduce(bind_rows, list()) = NULL, then `for (e in 1:nrow(NULL))` errors.
            throw new InvalidOperationException(
                "CNAqc analyze_peaks_general has no complex karyotype with n ≥ MinAbsoluteKaryotypeMutations (R error).");
        }

        // Data peaks in group_split(karyotype) order (C collation) — the order the R random stream is consumed in.
        RMersenneTwister? rng = options.BootstrapCount > 1 ? new RMersenneTwister(options.Seed) : null;
        var fits = new Dictionary<(int Major, int Minor), (KernelDensityEstimate Density, List<PurityDataPeak> Peaks)>();
        foreach ((int major, int minor) in analysis.OrderBy(k => $"{k.Major}:{k.Minor}", StringComparer.Ordinal))
        {
            fits[(major, minor)] = SimplePeakDetector(byKaryotype[(major, minor)], options.KernelAdjust, options.LegacyDensityCoordinates,
                $"{major}:{minor}", options.BootstrapCount, rng);
        }

        var results = new List<ComplexKaryotypePeakResult>(analysis.Count);
        var allExpected = new List<ComplexKaryotypeExpectedPeak>();
        foreach ((int major, int minor) in analysis)
        {
            List<double> vafs = byKaryotype[(major, minor)];
            (KernelDensityEstimate density, List<PurityDataPeak> peaks) = fits[(major, minor)];
            int ploidy = major + minor;
            int maxMultiplicity = Math.Max(Math.Max(major, minor), 1);
            var expected = new List<ComplexKaryotypeExpectedPeak>(maxMultiplicity);
            int matched = 0;
            for (int m = 1; m <= maxMultiplicity; m++)
            {
                double peak = CnaqcExpectedVaf(m, purity, ploidy); // expected_vaf_fun
                bool hit = peaks.Any(p => Math.Abs(p.X - peak) < options.PurityError);
                if (hit) matched++;
                expected.Add(new ComplexKaryotypeExpectedPeak(major, minor, m, peak, hit));
            }

            allExpected.AddRange(expected);
            double prop = (double)matched / (matched + (maxMultiplicity - matched));
            results.Add(new ComplexKaryotypePeakResult(major, minor, vafs.Count, matched, maxMultiplicity - matched, prop,
                prop >= 0.5, density, peaks, expected));
        }

        // summary: group_by(karyotype) (C collation) … arrange(desc(prop)) (stable).
        var summary = results
            .OrderBy(r => $"{r.MajorCopyNumber}:{r.MinorCopyNumber}", StringComparer.Ordinal)
            .OrderByDescending(r => r.MatchedProportion)
            .ToList();
        return new ComplexKaryotypePeakAnalysis(purity, true, summary, allExpected);
    }

    private static void ValidatePurityPeakMutation(PurityPeakMutation mutation)
    {
        if (double.IsNaN(mutation.Vaf) || mutation.Vaf < 0.0 || mutation.Vaf > 1.0)
            throw new ArgumentOutOfRangeException(nameof(mutation), mutation.Vaf, "VAF must be in [0, 1].");
        if (mutation.MajorCopyNumber < 0 || mutation.MinorCopyNumber < 0)
            throw new ArgumentOutOfRangeException(nameof(mutation), "Allele copy numbers must be non-negative.");
    }

    // CNAqc simple_peak_detector: density + peakPick peaks + histogram counts; with n_bootstrap > 1 the peaks of n
    // resamples (sample_n(replace = TRUE)) are pooled, distinct(x), heights re-read from the full density
    // (phase_to_density), and the new x values appended (distinct(x) after the full-data peaks).
    private static (KernelDensityEstimate Density, List<PurityDataPeak> Peaks) SimplePeakDetector(
        IReadOnlyList<double> vafs, double kernelAdjust, bool legacyCoordinates, string label,
        int bootstrapCount = 1, RMersenneTwister? rng = null)
    {
        (KernelDensityEstimate density, List<PurityDataPeak> peaks) = KdeSingleRun(vafs, kernelAdjust, legacyCoordinates, label);
        if (bootstrapCount > 1)
        {
            var pooled = new List<PurityDataPeak>();
            for (int b = 0; b < bootstrapCount; b++)
                pooled.AddRange(KdeSingleRun(BootstrapResample(vafs, rng!), kernelAdjust, legacyCoordinates, label).Peaks);
            AppendDistinctX(peaks, PhaseToDensity(DistinctX(pooled), density));
        }

        return (density, peaks);
    }

    private static (KernelDensityEstimate Density, List<PurityDataPeak> Peaks) KdeSingleRun(
        IReadOnlyList<double> vafs, double kernelAdjust, bool legacyCoordinates, string label)
    {
        KernelDensityEstimate density = StatisticsHelper.GaussianKernelDensity(vafs, kernelAdjust, legacyCoordinates: legacyCoordinates);
        List<PurityDataPeak> peaks = DetectKdePeaks(density, VafHistogram(vafs));
        if (peaks.Count == 0)
        {
            // simple_peak_detector: `if (indexes[1] == 0)` on an empty peak set is an R error.
            throw new InvalidOperationException($"Cannot find KDE peaks for {label}.");
        }

        return (density, peaks);
    }

    // dplyr sample_n(x, nrow(x), replace = TRUE) = x[sample.int(n, n, replace = TRUE)]: n draws R_unif_index(n).
    private static List<T> BootstrapResample<T>(IReadOnlyList<T> rows, RMersenneTwister rng)
    {
        var sample = new List<T>(rows.Count);
        for (int i = 0; i < rows.Count; i++) sample.Add(rows[rng.UnifIndex(rows.Count)]);
        return sample;
    }

    // dplyr distinct(x, .keep_all = TRUE): first row of each x.
    private static List<PurityDataPeak> DistinctX(IEnumerable<PurityDataPeak> peaks)
    {
        var seen = new HashSet<double>();
        return peaks.Where(p => seen.Add(p.X)).ToList();
    }

    private static void AppendDistinctX(List<PurityDataPeak> peaks, IEnumerable<PurityDataPeak> extra)
    {
        var seen = new HashSet<double>(peaks.Select(p => p.X));
        peaks.AddRange(extra.Where(p => seen.Add(p.X)));
    }

    // CNAqc phase_to_density: y = density$y[which.min(|density$x − x|)] (other columns kept).
    private static IEnumerable<PurityDataPeak> PhaseToDensity(IEnumerable<PurityDataPeak> peaks, KernelDensityEstimate density) =>
        peaks.Select(p => p with { Y = density.Y[NearestGridIndex(density, p.X)] });

    private static int NearestGridIndex(KernelDensityEstimate density, double x)
    {
        int w = 0;
        for (int i = 1; i < density.X.Count; i++)
        {
            if (Math.Abs(density.X[i] - x) < Math.Abs(density.X[w] - x)) w = i;
        }

        return w;
    }

    // CNAqc mixture_peak_detector: BMix means snapped to the KDE grid (x, y unrounded; counts from the same data's
    // histogram; never discarded). n_bootstrap > 1: the full-data fit still runs (R draws its random numbers) but the
    // peaks are those of n resampled fits (each snapped to its own resample's density), distinct(x), with heights read
    // from the density of one further resample (single_run_kde).
    private static List<PurityDataPeak> MixturePeakDetector(
        List<PurityPeakMutation> group, double kernelAdjust, bool legacyCoordinates, int bootstrapCount, RMersenneTwister rng)
    {
        if (group.Select(g => g.Vaf).Distinct().Count() == 1) return new List<PurityDataPeak>(); // returns NULL

        List<PurityDataPeak> SingleRun(IReadOnlyList<PurityPeakMutation> rows)
        {
            var fit = BMixFit(rows.Select(r => r.AlternateReads).ToList(), rows.Select(r => r.Depth).ToList(), new[] { 1, 2, 3, 4 }, rng).Best;
            var vafs = rows.Select(r => r.Vaf).ToList();
            KernelDensityEstimate density = StatisticsHelper.GaussianKernelDensity(vafs, kernelAdjust, legacyCoordinates: legacyCoordinates);
            int[] histogram = VafHistogram(vafs);
            return fit.Means.Select(b =>
            {
                int w = NearestGridIndex(density, b);
                double x = density.X[w];
                return new PurityDataPeak(x, density.Y[w], HistogramCount(histogram, Math.Round(x * 100, MidpointRounding.ToEven)), false, PurityPeakSource.Mixture);
            }).ToList();
        }

        List<PurityDataPeak> full = SingleRun(group);
        if (bootstrapCount <= 1) return full;

        var pooled = new List<PurityDataPeak>();
        for (int b = 0; b < bootstrapCount; b++) pooled.AddRange(SingleRun(BootstrapResample(group, rng)));
        var kdeSample = BootstrapResample(group, rng).Select(g => g.Vaf).ToList();
        KernelDensityEstimate kdeDensity = StatisticsHelper.GaussianKernelDensity(kdeSample, kernelAdjust, legacyCoordinates: legacyCoordinates);
        // single_run_kde: `hst[round(pks$x * 100)]` without the index-0 guard is an R error for a KDE peak at x < 0.005.
        if (DetectKdePeaks(kdeDensity, VafHistogram(kdeSample)).Any(p => Math.Round(p.X * 100, MidpointRounding.ToEven) == 0))
            throw new InvalidOperationException("CNAqc single_run_kde fails on a KDE peak at VAF 0 (R error).");
        return DistinctX(PhaseToDensity(DistinctX(pooled), kdeDensity));
    }

    // ---- CNAqc analyze_peaks_subclonal: subclonal simple segments (FIN-B24 F62) ----
}
