// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Subclonal copy-number states (Battenberg, port of Battenberg 3.0.0).</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{
    /// <summary>
    /// One integer allele-specific copy-number state of a (possibly sub-clonal) segment, present in a given
    /// fraction of tumour cells. Mirrors the Battenberg output (Nik-Zainal et al. 2012, <i>Cell</i> 149:994:
    /// <c>nMaj1_A, nMin1_A, frac1_A</c>): a state is a (major, minor) integer pair plus the cellular fraction.
    /// </summary>
    /// <param name="MajorCopyNumber">Major-allele integer copy number of this state (≥ 0).</param>
    /// <param name="MinorCopyNumber">Minor-allele integer copy number of this state (≥ 0).</param>
    /// <param name="CellFraction">Fraction of tumour cells carrying this state (Battenberg <c>frac1_A</c>/<c>frac2_A</c>; 1 for a
    /// clonal state). Reported unclamped, as Battenberg does, so it can leave [0, 1] when the segment BAF lies beyond
    /// the nearest edge's corners.</param>
    public readonly record struct SubclonalCopyNumberState(
        int MajorCopyNumber,
        int MinorCopyNumber,
        double CellFraction)
    {
        /// <summary>Total integer copy number of this state (major + minor).</summary>
        public int TotalCopyNumber => MajorCopyNumber + MinorCopyNumber;
    }

    /// <summary>
    /// Sub-clonal copy-number fit of one segment under the Battenberg two-population model (Nik-Zainal et al. 2012,
    /// <i>Cell</i> 149:994): a segment is either <b>clonal</b> (one integer state in all tumour cells) or
    /// <b>sub-clonal</b> (a mixture of two adjacent integer states, fractions summing to 1).
    /// </summary>
    /// <param name="Segment">The segment these states describe.</param>
    /// <param name="PrimaryState">State 1 (Battenberg <c>nMaj1_A, nMin1_A, frac1_A</c>).</param>
    /// <param name="SecondaryState">State 2 (Battenberg <c>frac2</c>), or <c>null</c> for a clonal segment.</param>
    /// <param name="IsSubclonal">True when two states were needed (the observed CN was not (near-)integer).</param>
    public readonly record struct SubclonalSegmentFit(
        AlleleSpecificSegmentSummary Segment,
        SubclonalCopyNumberState PrimaryState,
        SubclonalCopyNumberState? SecondaryState,
        bool IsSubclonal);

    /// <summary>
    /// Maximum distance of an allele-specific copy number from the nearest integer below which the pre-2026-09
    /// implementation called a segment clonal. <b>No longer used</b> by <see cref="FitSubclonalCopyNumber"/>, which now
    /// applies Battenberg's own BAF-space tolerance <see cref="BattenbergMaxBafDistance"/>; retained for API compatibility.
    /// </summary>
    [Obsolete("Unused since FIN-B24/F14: Battenberg determine_copynumber uses maxdist; kept for binary compatibility.")]
    public const double SubclonalIntegerTolerance = 0.05;

    /// <summary>
    /// Battenberg <c>maxdist = 0.01</c> (<c>callSubclones</c> default): a segment whose BAF lies within 0.01 of the BAF of
    /// the closest clonal corner of its nearest edge is called clonal (<c>if (abs(l − test.level) &lt; maxdist) pval = 1</c>).
    /// </summary>
    public const double BattenbergMaxBafDistance = 0.01;

    /// <summary>Battenberg <c>cn_upper_limit = 1000</c>: major copy number used when the minor allele is negative at BAF = 1.</summary>
    private const double BattenbergCopyNumberUpperLimit = 1000.0;

    /// <summary>Battenberg: a negative minor copy number is raised to <c>nMinor = 0.01</c> (major raised along the BAF line).</summary>
    private const double BattenbergMinimumMinorCopyNumber = 0.01;

    /// <summary>
    /// Fits each segment's allele-specific copy number to one integer state (clonal) or a mixture of two integer states
    /// with a sub-clonal cellular fraction, porting Battenberg's <c>determine_copynumber</c> (Wedge-lab/battenberg,
    /// R/fitcopynumber.R; <c>orderEdges</c> in R/orderEdges.R; Nik-Zainal et al. 2012, <i>Cell</i> 149:994):
    /// <list type="number">
    /// <item>l = max(BAF, 1 − BAF); ψ_all = ρ·ψ + 2(1 − ρ);
    /// nMajor = (ρ − 1 + l·ψ_all·2^(logR/γ))/ρ, nMinor = (ρ − 1 + (1 − l)·ψ_all·2^(logR/γ))/ρ (the ASCAT equations).
    /// A negative nMinor is raised to 0.01 with nMajor moved along the BAF line (or set to 1000 when l = 1).</item>
    /// <item>The four corners (⌊nMaj⌋/⌈nMaj⌉ × ⌊nMin⌋/⌈nMin⌉) and their BAF levels (1 − ρ + ρ·nMaj)/(2 − 2ρ + ρ·(nMaj + nMin))
    /// select the nearest edge of the copy-number square (<c>orderEdges</c> option 1; the two states differ in one
    /// allele), using the total-copy-number (logR) priority <c>ntot &lt; x + y + 1</c>.</item>
    /// <item>Clonal when the closest corner of that edge has |l − level| &lt; <see cref="BattenbergMaxBafDistance"/>
    /// (reported as that corner with fraction 1); otherwise sub-clonal with state 1 at fraction
    /// τ = (1 − ρ + ρ·M₂ − 2l(1 − ρ) − lρ(m₂ + M₂)) / (lρ(m₁ + M₁) − lρ(m₂ + M₂) − ρM₁ + ρM₂) and state 2 at 1 − τ, the
    /// fraction that reproduces the segment BAF as a mixture of the two states (Battenberg reports τ unclamped).</item>
    /// </list>
    /// A segment summary carries a single BAF value, i.e. its SNP BAFs have zero spread; Battenberg's per-SNP
    /// t-test then returns <c>pval = 0</c> (<c>if (is.na(sd(BAFke)) || sd(BAFke) == 0) pval = 0</c>), so only the
    /// maxdist rule can make a segment clonal — exactly the path ported here (with per-SNP BAFs, Battenberg's t-test is
    /// run by <see cref="FitSubclonalCopyNumberWithSnpTest"/>; the alternative solutions B–F, SDfrac and the seeded
    /// bootstrap confidence intervals by <see cref="FitSubclonalCopyNumberWithBootstrap"/>).
    /// </summary>
    /// <param name="segments">Segment summaries (e.g. from <see cref="SegmentAlleleSpecificAspcf"/>). Non-null; finite logR,
    /// BAF in [0, 1].</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (the ASCAT/Battenberg <c>psit</c>, &gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>Per-segment clonal/sub-clonal copy-number fits in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException">a segment has a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0,1], ψ ≤ 0, or γ ≤ 0.</exception>
    public static IReadOnlyList<SubclonalSegmentFit> FitSubclonalCopyNumber(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        double psiAll = MixtureCopiesPerCell(purity, ploidy); // psi = rho*psit + 2*(1-rho)
        var fits = new List<SubclonalSegmentFit>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            BattenbergSegmentGeometry g = BattenbergSegmentEdge(s, purity, psiAll, gamma, nameof(segments));
            // Constant (summary) BAF ⇒ Battenberg pval = 0; only maxdist can make the segment clonal.
            fits.Add(BattenbergSegmentFit(s, g, purity, g.ClosestDistance < BattenbergMaxBafDistance));
        }

        return fits;
    }

    /// <summary>
    /// Battenberg's significance level <c>siglevel = 0.05</c> (<c>callSubclones</c> default): a segment is sub-clonal when
    /// its SNP-BAF t-test p-value is ≤ siglevel (<c>if (pval[i] &lt;= siglevel)</c>).
    /// </summary>
    public const double BattenbergSignificanceLevel = 0.05;

    /// <summary>
    /// A segment summary together with the phased BAFs of its heterozygous SNPs (Battenberg <c>BAFphased</c> — the raw
    /// SNP BAF flipped to the segment's side, <c>ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>; used as given, not mirrored
    /// again), the input of Battenberg's per-segment clonality t-test.
    /// </summary>
    /// <param name="Segment">The segment (its <see cref="AlleleSpecificSegmentSummary.MeanBAF"/> is Battenberg's
    /// segment level <c>BAFseg</c>, its <see cref="AlleleSpecificSegmentSummary.MeanLogR"/> the segment logR).</param>
    /// <param name="PhasedSnpBafs">The segment's phased SNP BAFs (Battenberg <c>BAFke</c>), each in [0, 1]; may be empty.</param>
    public readonly record struct SubclonalSegmentSnpBafs(
        AlleleSpecificSegmentSummary Segment,
        IReadOnlyList<double> PhasedSnpBafs);

    /// <summary>
    /// A Battenberg segment fit with the clonality-test p-value Battenberg reports (<c>subcloneres$pval</c>).
    /// </summary>
    /// <param name="Fit">The clonal / sub-clonal copy-number fit.</param>
    /// <param name="PValue">Battenberg <c>pval</c>: 1 when the segment BAF is within maxdist of the closest corner level,
    /// 0 when the SNP BAFs have no spread (fewer than two SNPs, or constant), else the two-sided one-sample t-test p-value
    /// of the SNP BAFs against that corner level.</param>
    public readonly record struct SubclonalSegmentTestedFit(SubclonalSegmentFit Fit, double PValue);

    /// <summary>
    /// <see cref="FitSubclonalCopyNumber"/> with Battenberg's per-SNP clonality test (Wedge-lab/battenberg
    /// R/fitcopynumber.R, <c>determine_copynumber</c>; B24 F39). After the nearest edge and its closest corner level
    /// <c>test.level</c> are found exactly as in <see cref="FitSubclonalCopyNumber"/>:
    /// <code>
    /// if (is.na(sd(BAFke)) || sd(BAFke) == 0) pval = 0
    /// else pval = t.test(BAFke, alternative = "two.sided", mu = test.level)$p.value
    /// if (abs(l − test.level) &lt; maxdist) pval = 1
    /// sub-clonal ⇔ pval ≤ siglevel
    /// </code>
    /// so a segment whose BAF is beyond maxdist but whose SNP BAFs are too noisy to reject the clonal level is called
    /// clonal (and a segment with ≥ 2 distinct SNP BAFs that rejects it is sub-clonal). The t-test is
    /// <see cref="StatisticsHelper.OneSampleTTestPValue"/>; where R's <c>t.test</c> would stop ("data are essentially
    /// constant", Battenberg would abort) the p-value is taken as 0, as for exactly constant BAFs. With every
    /// segment's SNP BAFs constant this reproduces <see cref="FitSubclonalCopyNumber"/> exactly.
    /// </summary>
    /// <param name="segments">Segments with their phased SNP BAFs. Non-null; finite logR, BAF in [0, 1], SNP BAFs in [0, 1].</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <param name="significanceLevel">Battenberg <c>siglevel</c> ∈ [0, 1] (default <see cref="BattenbergSignificanceLevel"/>).</param>
    /// <param name="maxBafDistance">Battenberg <c>maxdist</c> ≥ 0 (default <see cref="BattenbergMaxBafDistance"/>).</param>
    /// <returns>Per-segment fits with Battenberg's p-value, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or a SNP-BAF list is null.</exception>
    /// <exception cref="ArgumentException">a segment has a non-finite logR, a BAF outside [0, 1], or a SNP BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0,1], ψ ≤ 0, γ ≤ 0, siglevel ∉ [0, 1] or maxdist &lt; 0 / not finite.</exception>
    public static IReadOnlyList<SubclonalSegmentTestedFit> FitSubclonalCopyNumberWithSnpTest(
        IReadOnlyList<SubclonalSegmentSnpBafs> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma,
        double significanceLevel = BattenbergSignificanceLevel,
        double maxBafDistance = BattenbergMaxBafDistance)
    {
        ArgumentNullException.ThrowIfNull(segments);
        ValidateAscatModelParameters(purity, ploidy, gamma);
        if (!(significanceLevel >= 0.0 && significanceLevel <= 1.0))
        {
            throw new ArgumentOutOfRangeException(nameof(significanceLevel), significanceLevel, "siglevel must lie in [0, 1].");
        }

        if (!(maxBafDistance >= 0.0) || double.IsPositiveInfinity(maxBafDistance))
        {
            throw new ArgumentOutOfRangeException(nameof(maxBafDistance), maxBafDistance, "maxdist must be finite and ≥ 0.");
        }

        double psiAll = MixtureCopiesPerCell(purity, ploidy);
        var fits = new List<SubclonalSegmentTestedFit>(segments.Count);
        foreach (SubclonalSegmentSnpBafs item in segments)
        {
            if (item.PhasedSnpBafs is null)
            {
                throw new ArgumentNullException(nameof(segments), "Every segment needs a (possibly empty) phased SNP-BAF list.");
            }

            foreach (double b in item.PhasedSnpBafs)
            {
                if (!(b >= 0.0 && b <= 1.0))
                {
                    throw new ArgumentException("Every phased SNP BAF must lie in [0, 1].", nameof(segments));
                }
            }

            BattenbergSegmentGeometry g = BattenbergSegmentEdge(item.Segment, purity, psiAll, gamma, nameof(segments));
            double pValue = StatisticsHelper.OneSampleTTestPValue(item.PhasedSnpBafs, g.ClosestTestLevel);
            if (double.IsNaN(pValue))
            {
                pValue = 0.0; // sd NA (< 2 SNPs) or 0 (constant BAF): Battenberg pval = 0
            }

            if (g.ClosestDistance < maxBafDistance)
            {
                pValue = 1.0;
            }

            fits.Add(new SubclonalSegmentTestedFit(
                BattenbergSegmentFit(item.Segment, g, purity, isClonal: !(pValue <= significanceLevel)),
                pValue));
        }

        return fits;
    }

    /// <summary>Battenberg per-segment geometry: mirrored BAF l, nearest edge, and the closest corner of that edge.</summary>
    private readonly record struct BattenbergSegmentGeometry(
        double L, double Maj1, double Min1, double Maj2, double Min2, bool FirstClosest, double ClosestDistance)
    {
        public double ClosestTestLevel => FirstClosest ? Level1 : Level2;

        public double Level1 { get; init; }

        public double Level2 { get; init; }

        /// <summary>Battenberg <c>levels[2]</c>/<c>levels[3]</c> (corner BAFs of the square, 0/0 ⇒ 0.5), x, y and ntot.</summary>
        public double SquareLevel2 { get; init; }

        public double SquareLevel3 { get; init; }

        public double X { get; init; }

        public double Y { get; init; }

        public double Ntot { get; init; }
    }

    /// <summary>The <c>determine_copynumber</c> body up to <c>whichclosestlevel.test</c> (shared by both clonality rules).</summary>
    private static BattenbergSegmentGeometry BattenbergSegmentEdge(
        AlleleSpecificSegmentSummary s, double rho, double psiAll, double gamma, string paramName)
    {
        if (!IsValidAlleleSignal(s.MeanLogR, s.MeanBAF))
        {
            throw new ArgumentException("Every segment needs a finite mean logR and a mean BAF in [0, 1].", paramName);
        }

        double l = MirrorBaf(s.MeanBAF); // max(BAF, 1 − BAF)
        double scaled = psiAll * Math.Pow(2.0, s.MeanLogR / gamma);
        double nMajor = (rho - 1.0 + l * scaled) / rho;
        double nMinor = (rho - 1.0 + (1.0 - l) * scaled) / rho;

        // Increase nMajor and nMinor together, to avoid impossible combinations (negative sub-clonal fractions).
        if (nMinor < 0.0)
        {
            nMajor = l == 1.0
                ? BattenbergCopyNumberUpperLimit
                : nMajor + l * (BattenbergMinimumMinorCopyNumber - nMinor) / (1.0 - l);
            nMinor = BattenbergMinimumMinorCopyNumber;
        }

        double x = Math.Floor(nMinor), y = Math.Floor(nMajor);
        double ntot = nMajor + nMinor;
        // Corners, sorted in the order of ascending BAF: (⌊M⌋,⌈m⌉), (⌈M⌉,⌈m⌉), (⌊M⌋,⌊m⌋), (⌈M⌉,⌊m⌋).
        double level2 = BattenbergCornerLevel(Math.Ceiling(nMajor), Math.Ceiling(nMinor), rho, zeroCornerIsBalanced: true);
        double level3 = BattenbergCornerLevel(Math.Floor(nMajor), Math.Floor(nMinor), rho, zeroCornerIsBalanced: true);
        (double maj1, double min1, double maj2, double min2) = BattenbergNearestEdge(level2, level3, l, ntot, x, y);

        // Clonality test on the corners of the nearest edge (test.levels carry no 0/0 correction in Battenberg).
        double testLevel1 = BattenbergCornerLevel(maj1, min1, rho, zeroCornerIsBalanced: false);
        double testLevel2 = BattenbergCornerLevel(maj2, min2, rho, zeroCornerIsBalanced: false);
        double dist1 = double.IsNaN(testLevel1) ? double.PositiveInfinity : Math.Abs(testLevel1 - l);
        double dist2 = double.IsNaN(testLevel2) ? double.PositiveInfinity : Math.Abs(testLevel2 - l);
        bool firstClosest = dist1 <= dist2; // which.min: first index on a tie
        return new BattenbergSegmentGeometry(l, maj1, min1, maj2, min2, firstClosest, firstClosest ? dist1 : dist2)
        {
            Level1 = testLevel1,
            Level2 = testLevel2,
            SquareLevel2 = level2,
            SquareLevel3 = level3,
            X = x,
            Y = y,
            Ntot = ntot,
        };
    }

    /// <summary>
    /// Battenberg output row: the closest corner with fraction 1 when clonal, else state 1 at τ and state 2 at 1 − τ.
    /// </summary>
    private static SubclonalSegmentFit BattenbergSegmentFit(
        AlleleSpecificSegmentSummary s, BattenbergSegmentGeometry g, double rho, bool isClonal)
    {
        double l = g.L, maj1 = g.Maj1, min1 = g.Min1, maj2 = g.Maj2, min2 = g.Min2;
        if (isClonal)
        {
            var clonal = g.FirstClosest
                ? new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), 1.0)
                : new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0);
            return new SubclonalSegmentFit(s, clonal, SecondaryState: null, IsSubclonal: false);
        }

        double tau = BattenbergTau(l, rho, maj1, min1, maj2, min2);
        return new SubclonalSegmentFit(
            s,
            new SubclonalCopyNumberState(AscatCopyNumberToInt(maj1), AscatCopyNumberToInt(min1), tau),
            new SubclonalCopyNumberState(AscatCopyNumberToInt(maj2), AscatCopyNumberToInt(min2), 1.0 - tau),
            IsSubclonal: true);
    }

    /// <summary>
    /// Battenberg corner BAF level <c>(1 − ρ + ρ·nMaj)/(2 − 2ρ + ρ·(nMaj + nMin))</c>; for the corner (0, 0) at ρ = 1
    /// (0/0) the square's <c>levels</c> are set to 0.5 while the edge <c>test.levels</c> stay NaN.
    /// </summary>
    private static double BattenbergCornerLevel(double nMaj, double nMin, double rho, bool zeroCornerIsBalanced)
    {
        if (zeroCornerIsBalanced && nMaj == 0.0 && nMin == 0.0)
        {
            return BalancedBaf;
        }

        return (1.0 - rho + rho * nMaj) / (2.0 - 2.0 * rho + rho * (nMaj + nMin));
    }

    /// <summary>
    /// Battenberg <c>orderEdges(levels, l, ntot, x, y)</c>, first option (the nearest edge): case 1/2a when
    /// l &gt; levels[3], case 2c when l &gt; levels[2], else case 2b; within each case the logR criterion
    /// <c>ntot &lt; x + y + 1</c> picks the lower- or higher-total edge. Returns (nMaj1, nMin1, nMaj2, nMin2).
    /// </summary>
    private static (double Maj1, double Min1, double Maj2, double Min2) BattenbergNearestEdge(
        double level2, double level3, double l, double ntot, double x, double y)
    {
        bool lowerTotal = ntot < x + y + 1.0;
        if (l > level3)
        {
            // case 1 or 2a
            return lowerTotal ? (y, x, y + 1.0, x) : (y + 1.0, x, y + 1.0, x + 1.0);
        }

        if (l > level2)
        {
            // case 2c
            return lowerTotal ? (y, x, y, x + 1.0) : (y + 1.0, x, y + 1.0, x + 1.0);
        }

        // case 2b
        return lowerTotal ? (y, x, y, x + 1.0) : (y, x + 1.0, y + 1.0, x + 1.0);
    }

    /// <summary>
    /// One of Battenberg's six candidate sub-clonal solutions A–F of a segment (<c>callSubclones</c> extended output
    /// columns <c>nMaj1_X, nMin1_X, frac1_X, nMaj2_X, nMin2_X, frac2_X, SDfrac_X, SDfrac_X_BS, frac1_X_0.025, frac1_X_0.975</c>).
    /// A solution whose edge would need a negative copy number is <c>NA</c> in Battenberg: its copy numbers are
    /// <c>null</c> and its numbers NaN (such solutions are listed last).
    /// </summary>
    /// <param name="MajorCopyNumber1">nMaj1 (state 1 major), or null (NA).</param>
    /// <param name="MinorCopyNumber1">nMin1, or null.</param>
    /// <param name="Fraction1">frac1 = τ (unclamped), NaN when NA.</param>
    /// <param name="MajorCopyNumber2">nMaj2 (state 2 major), or null.</param>
    /// <param name="MinorCopyNumber2">nMin2, or null.</param>
    /// <param name="Fraction2">frac2 = 1 − τ.</param>
    /// <param name="FractionSd">SDfrac: |τ(l + s) − τ|/2 + |τ(l − s) − τ|/2 with s = sd(BAFke)/√n (NaN for &lt; 2 SNPs).</param>
    /// <param name="FractionBootstrapSd">SDfrac_BS: sd of the bootstrap τ values (NaN when any is NaN).</param>
    /// <param name="Fraction1Lower">frac1_0.025: the 25th smallest bootstrap τ (NaN when fewer than 25 non-NaN values).</param>
    /// <param name="Fraction1Upper">frac1_0.975: the 975th smallest bootstrap τ (NaN when fewer than 975).</param>
    public readonly record struct BattenbergSubclonalSolution(
        int? MajorCopyNumber1,
        int? MinorCopyNumber1,
        double Fraction1,
        int? MajorCopyNumber2,
        int? MinorCopyNumber2,
        double Fraction2,
        double FractionSd,
        double FractionBootstrapSd,
        double Fraction1Lower,
        double Fraction1Upper);

    /// <summary>
    /// A Battenberg <c>determine_copynumber</c> output row with its alternative solutions and confidence intervals.
    /// </summary>
    /// <param name="Fit">Solution A as a clonal / sub-clonal fit (identical to <see cref="FitSubclonalCopyNumberWithSnpTest"/>).</param>
    /// <param name="PValue">Battenberg <c>pval</c>.</param>
    /// <param name="Baf">Battenberg <c>BAF</c> column: the mirrored segment level l = max(BAFseg, 1 − BAFseg).</param>
    /// <param name="TotalCopyNumber">Battenberg <c>ntot</c> = nMajor + nMinor (raw, after the negative-minor correction).</param>
    /// <param name="Solutions">Solutions A–F (six entries) for a sub-clonal segment; empty for a clonal one (Battenberg NA).</param>
    public sealed record BattenbergSegmentCall(
        SubclonalSegmentFit Fit,
        double PValue,
        double Baf,
        double TotalCopyNumber,
        IReadOnlyList<BattenbergSubclonalSolution> Solutions)
    {
        /// <summary>
        /// True when Battenberg <c>mask_high_cn_segments</c> masked this row (B24 F60): its <c>nMaj1_A</c>, <c>nMin1_A</c>,
        /// <c>nMaj2_A</c>, <c>nMin2_A</c> are <c>NA</c> in Battenberg's output (solution A's copy numbers are then null here;
        /// <see cref="Fit"/> keeps the pre-mask states, the fractions and solutions B–F are unchanged, as in R).
        /// </summary>
        public bool IsMasked { get; init; }
    }

    /// <summary>Battenberg <c>noperms = 1000</c> (<c>callSubclones</c> default).</summary>
    public const int BattenbergDefaultPermutations = 1000;

    /// <summary>
    /// <see cref="FitSubclonalCopyNumberWithSnpTest"/> plus the rest of Battenberg's <c>determine_copynumber</c> row
    /// (Wedge-lab/battenberg R/fitcopynumber.R; B24 F41): for a sub-clonal segment all six <c>orderEdges</c> solutions
    /// A–F (rows with a negative copy number become NA and move to the end), each with τ, the delta-method
    /// <c>SDfrac</c> and the bootstrap <c>SDfrac_BS</c> / <c>frac1_0.025</c> / <c>frac1_0.975</c>:
    /// <code>
    /// for (option in 1:6) for (j in 1:noperms) {
    ///   permBAFs = sample(BAFke, length(BAFke), replace = T); permMeanBAF = mean(permBAFs)
    ///   permFraction[j] = τ(permMeanBAF; option) }
    /// SDfrac_BS = sd(permFraction); frac1_0.025 = sort(permFraction)[25]; frac1_0.975 = sort(permFraction)[975]
    /// </code>
    /// Note the bounds are fixed order statistics, not <c>quantile()</c>: they are the 2.5 %/97.5 % points only for
    /// noperms = 1000 and NA for fewer resamples. Resampling uses a port of R's default generator — Mersenne-Twister
    /// seeded by <c>set.seed(seed)</c> (initial scrambling + LCG fill) with R ≥ 3.6 "Rejection" <c>sample</c>
    /// (<c>R_unif_index</c>/<c>rbits</c>) — consuming draws in Battenberg's order (segments in input order, options
    /// A–F, NA options included), so the resamples are R's: <c>set.seed(seed); determine_copynumber(…)</c> is reproduced
    /// (the bootstrap columns agree to ~1e−15; R's long-double <c>mean</c>/<c>sd</c> differ only in the last bits).
    /// <c>callSubclones</c> calls <c>set.seed</c> once and then <c>determine_copynumber</c> twice (before/after
    /// <c>merge_segments</c>), so its final table matches this method only for the RNG state at the second call —
    /// <see cref="CallBattenbergSubclones"/> runs the whole sequence on one stream (B24 F61).
    /// </summary>
    /// <param name="segments">Segments with their phased SNP BAFs (e.g. from <see cref="BuildBattenbergSegments"/>).</param>
    /// <param name="purity">Fitted tumour purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Fitted tumour ploidy ψ (&gt; 0).</param>
    /// <param name="seed">R <c>set.seed</c> value (Battenberg <c>seed</c>; any 32-bit integer).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <param name="significanceLevel">Battenberg <c>siglevel</c> ∈ [0, 1].</param>
    /// <param name="maxBafDistance">Battenberg <c>maxdist</c> ≥ 0.</param>
    /// <param name="permutations">Battenberg <c>noperms</c> ≥ 1 (default 1000).</param>
    /// <returns>One call per segment, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or a SNP-BAF list is null.</exception>
    /// <exception cref="ArgumentException">invalid segment or SNP BAF (as <see cref="FitSubclonalCopyNumberWithSnpTest"/>).</exception>
    /// <exception cref="ArgumentOutOfRangeException">invalid ρ, ψ, γ, siglevel, maxdist, or permutations &lt; 1.</exception>
    public static IReadOnlyList<BattenbergSegmentCall> FitSubclonalCopyNumberWithBootstrap(
        IReadOnlyList<SubclonalSegmentSnpBafs> segments,
        double purity,
        double ploidy,
        int seed,
        double gamma = AscatSequencingGamma,
        double significanceLevel = BattenbergSignificanceLevel,
        double maxBafDistance = BattenbergMaxBafDistance,
        int permutations = BattenbergDefaultPermutations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(permutations, 1);
        return BattenbergDetermineCopyNumber(
            segments, purity, ploidy, new RMersenneTwister(seed), gamma, significanceLevel, maxBafDistance, permutations);
    }

    /// <summary>
    /// <c>determine_copynumber</c> drawing its bootstrap resamples from <paramref name="rng"/> — shared by
    /// <see cref="FitSubclonalCopyNumberWithBootstrap"/> (fresh <c>set.seed</c>) and <see cref="CallBattenbergSubclones"/>
    /// (one stream across both calls, as <c>callSubclones</c>).
    /// </summary>
    private static List<BattenbergSegmentCall> BattenbergDetermineCopyNumber(
        IReadOnlyList<SubclonalSegmentSnpBafs> segments,
        double purity,
        double ploidy,
        RMersenneTwister rng,
        double gamma,
        double significanceLevel,
        double maxBafDistance,
        int permutations)
    {
        IReadOnlyList<SubclonalSegmentTestedFit> tested = FitSubclonalCopyNumberWithSnpTest(
            segments, purity, ploidy, gamma, significanceLevel, maxBafDistance);
        double psiAll = MixtureCopiesPerCell(purity, ploidy);
        var calls = new List<BattenbergSegmentCall>(tested.Count);
        for (int i = 0; i < tested.Count; i++)
        {
            SubclonalSegmentTestedFit t = tested[i];
            BattenbergSegmentGeometry g = BattenbergSegmentEdge(segments[i].Segment, purity, psiAll, gamma, nameof(segments));
            IReadOnlyList<BattenbergSubclonalSolution> solutions = t.Fit.IsSubclonal
                ? BattenbergSolutions(g, purity, segments[i].PhasedSnpBafs, permutations, rng)
                : Array.Empty<BattenbergSubclonalSolution>();
            calls.Add(new BattenbergSegmentCall(t.Fit, t.PValue, g.L, g.Ntot, solutions));
        }

        return calls;
    }

    /// <summary>The sub-clonal branch of <c>determine_copynumber</c>: six solutions with SDfrac and bootstrap columns.</summary>
    private static BattenbergSubclonalSolution[] BattenbergSolutions(
        BattenbergSegmentGeometry g, double rho, IReadOnlyList<double> bafke, int noperms, RMersenneTwister rng)
    {
        (double Maj1, double Min1, double Maj2, double Min2)[] edges = BattenbergAllEdges(g.SquareLevel2, g.SquareLevel3, g.L, g.Ntot, g.X, g.Y);
        double l = g.L;
        int n = bafke.Count;
        double sdl = n < 2 ? double.NaN : Math.Sqrt(StatisticsHelper.RSampleVariance(bafke)) / Math.Sqrt(n);
        var result = new BattenbergSubclonalSolution[edges.Length];
        var perm = new double[n];
        var permFraction = new double[noperms];
        for (int o = 0; o < edges.Length; o++)
        {
            (double maj1, double min1, double maj2, double min2) = edges[o];
            double tau = BattenbergTau(l, rho, maj1, min1, maj2, min2);
            double sdtau = Math.Abs(BattenbergTau(l + sdl, rho, maj1, min1, maj2, min2) - tau) / 2
                           + Math.Abs(BattenbergTau(l - sdl, rho, maj1, min1, maj2, min2) - tau) / 2;
            for (int j = 0; j < noperms; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    perm[k] = bafke[rng.UnifIndex(n)];
                }

                double permMean = n == 0 ? double.NaN : StatisticsHelper.RMean(perm);
                permFraction[j] = BattenbergTau(permMean, rho, maj1, min1, maj2, min2);
            }

            double sdBoot = Math.Sqrt(StatisticsHelper.RSampleVariance(permFraction));
            var ordered = permFraction.Where(v => !double.IsNaN(v)).ToArray(); // sort() drops NA/NaN
            Array.Sort(ordered);
            bool na = double.IsNaN(maj1);
            result[o] = new BattenbergSubclonalSolution(
                na ? null : AscatCopyNumberToInt(maj1),
                na ? null : AscatCopyNumberToInt(min1),
                tau,
                na ? null : AscatCopyNumberToInt(maj2),
                na ? null : AscatCopyNumberToInt(min2),
                1.0 - tau,
                sdtau,
                sdBoot,
                ordered.Length >= 25 ? ordered[24] : double.NaN,
                ordered.Length >= 975 ? ordered[974] : double.NaN);
        }

        return result;
    }

    /// <summary>
    /// Battenberg sub-clonal fraction of state 1 at BAF <paramref name="l"/>:
    /// τ = (1 − ρ + ρM₂ − 2l(1 − ρ) − lρ(m₂ + M₂)) / (lρ(m₁ + M₁) − lρ(m₂ + M₂) − ρM₁ + ρM₂) (R operation order).
    /// </summary>
    private static double BattenbergTau(double l, double rho, double maj1, double min1, double maj2, double min2) =>
        (1.0 - rho + rho * maj2 - 2.0 * l * (1.0 - rho) - l * rho * (min2 + maj2))
        / (l * rho * (min1 + maj1) - l * rho * (min2 + maj2) - rho * maj1 + rho * maj2);

    /// <summary>
    /// Battenberg <c>orderEdges(levels, l, ntot, x, y)</c>, all six options (R/orderEdges.R), negative-copy-number rows
    /// set to NaN and moved to the end (<c>determine_copynumber</c>: <c>rbind(all.edges[-na.indices,], all.edges[na.indices,])</c>).
    /// </summary>
    private static (double Maj1, double Min1, double Maj2, double Min2)[] BattenbergAllEdges(
        double level2, double level3, double l, double ntot, double x, double y)
    {
        double[] maj1, min1, maj2, min2;
        bool lower = ntot < x + y + 1.0;
        if (l > level3)
        {
            // case 1 or 2a
            if (lower)
            {
                maj1 = new[] { y, y - 1, y, y + 1, y + 1, y + 1 };
                min1 = new[] { x, x, x, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 2, y + 1, y + 1, y + 1 };
                min2 = new[] { x, x, x, x + 1, x + 1, x + 2 };
            }
            else
            {
                maj1 = new[] { y + 1, y + 1, y + 1, y, y - 1, y };
                min1 = new[] { x, x - 1, x, x, x, x };
                maj2 = new[] { y + 1, y + 1, y + 1, y + 1, y + 1, y + 2 };
                min2 = new[] { x + 1, x + 1, x + 2, x, x, x };
            }
        }
        else if (l > level2)
        {
            // case 2c
            if (lower)
            {
                maj1 = new[] { y, y, y, y + 1, y + 1, y + 1 };
                min1 = new[] { x, x - 1, x, x, x - 1, x };
                maj2 = new[] { y, y, y, y + 1, y + 1, y + 1 };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 2 };
            }
            else
            {
                maj1 = new[] { y + 1, y + 1, y + 1, y, y, y };
                min1 = new[] { x, x - 1, x, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 1, y, y, y };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 2 };
            }
        }
        else
        {
            // case 2b
            if (lower)
            {
                maj1 = new[] { y, y, y, y, y - 1, y };
                min1 = new[] { x, x - 1, x, x + 1, x + 1, x + 1 };
                maj2 = new[] { y, y, y, y + 1, y + 1, y + 2 };
                min2 = new[] { x + 1, x + 1, x + 2, x + 1, x + 1, x + 1 };
            }
            else
            {
                maj1 = new[] { y, y - 1, y, y, y, y };
                min1 = new[] { x + 1, x + 1, x + 1, x, x - 1, x };
                maj2 = new[] { y + 1, y + 1, y + 2, y, y, y };
                min2 = new[] { x + 1, x + 1, x + 1, x + 1, x + 1, x + 2 };
            }
        }

        var valid = new List<(double, double, double, double)>(6);
        var invalid = new List<(double, double, double, double)>(6);
        for (int i = 0; i < 6; i++)
        {
            if (maj1[i] < 0 || min1[i] < 0 || maj2[i] < 0 || min2[i] < 0)
            {
                invalid.Add((double.NaN, double.NaN, double.NaN, double.NaN));
            }
            else
            {
                valid.Add((maj1[i], min1[i], maj2[i], min2[i]));
            }
        }

        valid.AddRange(invalid);
        return valid.ToArray();
    }
}
