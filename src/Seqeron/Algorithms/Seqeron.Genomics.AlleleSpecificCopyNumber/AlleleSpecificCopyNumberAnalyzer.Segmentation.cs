// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Allele-specific piecewise-constant segmentation (ASCAT ascat.aspcf / asmultipcf, port of ASCAT 3.2.0).</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{
    /// <summary>
    /// Default ASPCF penalty used when a caller does not supply one: ASCAT's <c>ascat.aspcf(..., penalty = 70)</c>
    /// (VanLoo-lab/ascat, ascat.aspcf.R; Ross et al. 2021, <i>Bioinformatics</i> 37:1909). The penalty is charged per
    /// breakpoint on the <b>standardised</b> joint cost (each track's squared error divided by its MAD-based variance),
    /// so it is scale-free. The copynumber package's single-track default is γ = 40 (Nilsen et al. 2012).
    /// </summary>
    public const double AspcfDefaultPenalty = 70.0;

    /// <summary>ASCAT/copynumber minimum ASPCF segment length <c>kmin = 6</c> (<c>fastAspcf(..., 6, ...)</c>).</summary>
    private const int AspcfMinSegmentLength = 6;

    /// <summary>ASCAT <c>madWins(x, 2.5, 25)</c>: winsorisation at τ = 2.5 MAD-standard deviations.</summary>
    private const double AspcfWinsorTau = 2.5;

    /// <summary>ASCAT <c>madWins</c>/<c>getMad</c> running-median half-window k = 25 (filter width 2k + 1 = 51).</summary>
    private const int AspcfMedianHalfWindow = 25;

    /// <summary>ASCAT <c>fastAspcf</c> window size <c>w = 1000</c> and overlap <c>d = 100</c>.</summary>
    private const int AspcfWindowSize = 1000;

    /// <summary>ASCAT <c>fastAspcf</c> window overlap <c>d = 100</c>.</summary>
    private const int AspcfWindowOverlap = 100;

    /// <summary>ASCAT: re-run segmentation with the next larger penalty while ≥ 800 distinct logR levels remain.</summary>
    private const int AspcfMaxSegmentLevels = 800;

    /// <summary>ASCAT penalty ladder <c>segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140))</c>, kept where ≥ penalty.</summary>
    private static readonly double[] AspcfPenaltyLadder = { 35.0, 50.0, 70.0, 100.0, 140.0 };

    /// <summary>
    /// Allele-Specific Piecewise Constant Fitting (ASPCF) as run by ASCAT (<c>ascat.aspcf</c>, VanLoo-lab/ascat
    /// ascat.aspcf.R; the bivariate PCF of the copynumber package, Nilsen et al. 2012, <i>BMC Genomics</i> 13:591;
    /// Ross et al. 2021, <i>Bioinformatics</i> 37:1909), ported per chromosome on the heterozygous loci:
    /// <list type="number">
    /// <item>logR and mirrored BAF are MAD-winsorised (<c>madWins(x, 2.5, 25)</c>: running median of width 51 with
    /// Tukey end rule, residuals clipped at ±2.5·MAD).</item>
    /// <item>Chromosomes with fewer than 6 loci form one segment (mean winsorised mirrored BAF).</item>
    /// <item>Otherwise <c>fastAspcf</c>: overlapping windows (1000 loci, overlap 100); in each window the joint
    /// penalised cost <c>Σ_segments [SSE_logR/sd₁² + SSE_BAF/sd₂²] + γ·(#breakpoints)</c> is minimised exactly by the
    /// PCF dynamic program with minimum segment length kmin = 6 (<c>aspcfpart</c>), sd₁ / sd₂ being the MAD of the
    /// residuals from a running median (<c>getMad</c>) of logR and of the flipped BAF min(b, 1 − b); a window whose
    /// sd is 0 or undefined (e.g. noise-free data) contributes no breakpoint.</item>
    /// <item>Each segment's BAF is 0.5 + mean|b − 0.5|, shrunk to exactly 0.5 when
    /// <c>sqrt(sd₂² + μ²) &lt; 2·sd₂</c> (balanced); its logR is the mean of the raw (unwinsorised) logR.</item>
    /// <item>While ≥ 800 distinct segment logR levels remain, segmentation is repeated with the next larger penalty
    /// of the ASCAT ladder (35, 50, 70, 100, 140).</item>
    /// </list>
    /// Breakpoints never cross a contig boundary (loci are grouped into contiguous same-chromosome runs, input order).
    /// Every supplied locus is treated as a germline-heterozygous SNP; for mixed heterozygous / homozygous input (logR
    /// averaging over homozygous probes, homozygous-stretch resegmentation) use the germline-aware overload
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/> (B24 F36).
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome. LogR must be
    /// finite and BAF in [0, 1].</param>
    /// <param name="penalty">ASCAT/copynumber penalty γ &gt; 0 per breakpoint on the standardised cost
    /// (see <see cref="AspcfDefaultPenalty"/>).</param>
    /// <returns>The segment summaries (mean logR, ASPCF mirrored BAF) in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException">a locus has a null chromosome label, a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static IReadOnlyList<AlleleSpecificSegmentSummary> SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(loci);
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        var ordered = new List<AlleleSpecificLocus>();
        foreach (AlleleSpecificLocus locus in loci)
        {
            if (locus.Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            if (!IsValidAlleleSignal(locus.LogR, locus.BAF))
            {
                throw new ArgumentException("Every locus needs a finite logR and a BAF in [0, 1].", nameof(loci));
            }

            ordered.Add(locus);
        }

        // Contiguous same-chromosome runs (ASCATobj$chr), in input order.
        var runs = new List<(int Lo, int Hi)>();
        int start = 0;
        while (start < ordered.Count)
        {
            int end = start;
            while (end + 1 < ordered.Count && ordered[end + 1].Chromosome == ordered[start].Chromosome)
            {
                end++;
            }

            runs.Add((start, end));
            start = end + 1;
        }

        // segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140)); segmentlengths[segmentlengths >= penalty]
        var ladder = new List<double> { penalty };
        foreach (double p in AspcfPenaltyLadder)
        {
            if (p > penalty)
            {
                ladder.Add(p);
            }
        }

        List<AlleleSpecificSegmentSummary> result = new();
        foreach (double segmentPenalty in ladder)
        {
            result = new List<AlleleSpecificSegmentSummary>();
            foreach ((int lo, int hi) in runs)
            {
                SegmentChromosomeAspcf(ordered, lo, hi, segmentPenalty, result);
            }

            var levels = new HashSet<double>();
            foreach (AlleleSpecificSegmentSummary s in result)
            {
                levels.Add(s.MeanLogR);
            }

            if (levels.Count < AspcfMaxSegmentLevels)
            {
                break;
            }
        }

        return result;
    }

    /// <summary>
    /// ASCAT per-chromosome ASPCF (the loop body of <c>ascat.aspcf</c>) on one run of loci [lo, hi]; appends the
    /// segment summaries to <paramref name="output"/>.
    /// </summary>
    private static void SegmentChromosomeAspcf(
        List<AlleleSpecificLocus> loci, int lo, int hi, double penalty,
        List<AlleleSpecificSegmentSummary> output)
    {
        int n = hi - lo + 1;
        var lr = new double[n];
        var baf = new double[n];
        var mirrored = new double[n];
        for (int i = 0; i < n; i++)
        {
            lr[i] = loci[lo + i].LogR;
            baf[i] = loci[lo + i].BAF;
            mirrored[i] = MirrorBaf(baf[i]); // ifelse(bafsel > 0.5, bafsel, 1 - bafsel)
        }

        double[] lrWins = MadWinsorize(lr, AspcfWinsorTau, AspcfMedianHalfWindow);
        double[] bafWinsMirrored = MadWinsorize(mirrored, AspcfWinsorTau, AspcfMedianHalfWindow);
        var bafWins = new double[n];
        for (int i = 0; i < n; i++)
        {
            bafWins[i] = OrientBaf(baf[i], bafWinsMirrored[i]);
        }

        if (n < AspcfMinSegmentLength)
        {
            // logRASPCF = mean(logRaveraged); bafASPCF = mean(bafselwinsmirrored); level re-adapted to mean(lr).
            output.Add(BuildAspcfSummary(loci, lo, 0, n, lr, Mean(bafWinsMirrored, 0, n)));
            return;
        }

        (int[] breakpoints, double[] segmentBaf) = FastAspcf(lrWins, bafWins, AspcfMinSegmentLength, penalty);
        for (int s = 0; s + 1 < breakpoints.Length; s++)
        {
            // Segment s covers the 0-based loci [breakpoints[s], breakpoints[s+1]); logR level = mean raw logR.
            output.Add(BuildAspcfSummary(loci, lo, breakpoints[s], breakpoints[s + 1], lr, segmentBaf[s]));
        }
    }

    /// <summary>Summary of loci [from, to) of the run starting at <paramref name="lo"/>: mean raw logR and the ASPCF BAF.</summary>
    private static AlleleSpecificSegmentSummary BuildAspcfSummary(
        List<AlleleSpecificLocus> loci, int lo, int from, int to, double[] rawLogR, double segmentBaf) =>
        new(
            Chromosome: loci[lo + from].Chromosome,
            Start: loci[lo + from].Position,
            End: loci[lo + to - 1].Position,
            MeanLogR: Mean(rawLogR, from, to),
            MeanBAF: segmentBaf,
            LocusCount: to - from);

    /// <summary>R <c>mean(x[from..to))</c>, bit-identical to R's long-double two-pass <c>real_mean</c>
    /// (<see cref="StatisticsHelper.ExtendedPrecisionMean"/>, B24 F66). Call sites: ascat.aspcf <c>mean(logRaveraged)</c>,
    /// <c>mean(bafselwinsmirrored)</c>, fastAspcf <c>yhat1 = mean(logR[frst:last])</c>, exactPcf <c>mean(y)</c>.</summary>
    private static double Mean(double[] x, int from, int to) =>
        StatisticsHelper.ExtendedPrecisionMean(x.AsSpan(from, to - from));

    /// <summary>R <c>mean(x[from..to), na.rm = TRUE)</c>: NaN entries (R NA) are dropped first, NaN when none remains;
    /// bit-identical to R (<see cref="StatisticsHelper.ExtendedPrecisionMean"/> with skipNaN, B24 F66). Call sites:
    /// ascat.aspcf / ascat.asmultipcf <c>mean(lr[startprobe:endprobe], na.rm=TRUE)</c>, <c>mean(lr, na.rm=TRUE)</c>,
    /// <c>logRaveraged[i] = mean(lrwins[start:end], na.rm=TRUE)</c>, asmultipcf <c>mean(x, na.rm=TRUE)</c>.</summary>
    private static double MeanNaRm(double[] x, int from, int to) =>
        StatisticsHelper.ExtendedPrecisionMean(x.AsSpan(from, to - from), skipNaN: true);

    /// <summary>
    /// ASCAT <c>fastAspcf(logR, allB, kmin, gamma)</c>: windowed exact bivariate PCF. Returns the segment boundaries
    /// (0-based, <c>[b₀ = 0, …, b_S = N]</c>) and each segment's BAF level <c>0.5 + μ</c> (μ shrunk to 0 when
    /// <c>sqrt(sd₂² + μ²) &lt; 2·sd₂</c>).
    /// </summary>
    private static (int[] Breakpoints, double[] SegmentBaf) FastAspcf(double[] logR, double[] allB, int kmin, double gamma)
    {
        int bigN = logR.Length;
        int startw = -AspcfWindowOverlap;
        int stopw = AspcfWindowSize - AspcfWindowOverlap;
        int validWindows = 0;
        double var2 = 0.0;
        var breakpts = new List<int> { 0 };
        while (true)
        {
            int from = Math.Max(1, startw);
            int to = Math.Min(stopw, bigN);
            int len = to - from + 1;
            var logRPart = new double[len];
            var allBFlip = new double[len];
            for (int i = 0; i < len; i++)
            {
                logRPart[i] = logR[from - 1 + i];
                double b = allB[from - 1 + i];
                allBFlip[i] = ToAscatBaf(b); // min(b, 1 − b)
            }

            double sd1 = GetMad(logRPart, AspcfMedianHalfWindow);
            double sd2 = GetMad(allBFlip, AspcfMedianHalfWindow);
            if (!double.IsNaN(sd1) && !double.IsNaN(sd2) && sd1 != 0.0 && sd2 != 0.0)
            {
                List<int> part = AspcfPart(logRPart, allBFlip, startw, stopw, AspcfWindowOverlap, sd1, sd2, bigN, kmin, gamma);
                int last = breakpts[^1];
                foreach (int bp in part)
                {
                    if (bp > last)
                    {
                        breakpts.Add(bp); // breakptspart[larger]
                    }
                }

                var2 += sd2 * sd2;
                validWindows++;
            }

            if (stopw < bigN + AspcfWindowOverlap)
            {
                startw = Math.Min(stopw - 2 * AspcfWindowOverlap + 1, bigN - 2 * AspcfWindowOverlap);
                stopw = startw + AspcfWindowSize;
            }
            else
            {
                break;
            }
        }

        // breakpts <- unique(c(breakpts, N))
        var unique = new List<int>();
        foreach (int bp in breakpts.Append(bigN))
        {
            if (!unique.Contains(bp))
            {
                unique.Add(bp);
            }
        }

        if (validWindows == 0)
        {
            validWindows = 1; // "just in case the sd-test never passes"
        }

        double sd2Pooled = Math.Sqrt(var2 / validWindows);
        var segmentBaf = new double[unique.Count - 1];
        for (int s = 0; s + 1 < unique.Count; s++)
        {
            int first = unique[s], lastExclusive = unique[s + 1];
            // mu = mean(abs(yi2 - 0.5)) (0 for an empty segment), R's long-double mean (B24 F66).
            var deviation = new double[lastExclusive - first];
            for (int i = first; i < lastExclusive; i++)
            {
                deviation[i - first] = Math.Abs(allB[i] - BalancedBaf);
            }

            double mu = deviation.Length > 0 ? StatisticsHelper.ExtendedPrecisionMean(deviation) : 0.0;
            segmentBaf[s] = AscatShrunkBafLevel(mu, sd2Pooled);
        }

        return (unique.ToArray(), segmentBaf);
    }

    /// <summary>
    /// ASCAT segment BAF level <c>0.5 + μ</c> with μ = mean |b − 0.5|, shrunk to 0 when <c>sqrt(sd² + μ²) &lt; 2·sd</c>
    /// (ascat.aspcf <c>fastAspcf</c> and ascat.asmultipcf BAF correction, verbatim).
    /// </summary>
    private static double AscatShrunkBafLevel(double mu, double sd) =>
        (Math.Sqrt(sd * sd + mu * mu) < 2.0 * sd ? 0.0 : mu) + BalancedBaf;

    /// <summary>
    /// ASCAT/copynumber <c>aspcfpart</c>: exact bivariate PCF (minimum segment length <paramref name="kmin"/>) on one
    /// window, cost = SSE₁/sd₁² + SSE₂/sd₂² per segment plus γ per breakpoint. Returns the global breakpoints (1-based
    /// "last index of a segment", as in R) that fall inside the window's use-range [a + d, b − d].
    /// </summary>
    private static List<int> AspcfPart(
        double[] y1, double[] y2, int a, int b, int d, double sd1, double sd2, int bigN, int kmin, double gamma)
    {
        int from = Math.Max(1, a);
        int useFrom = Math.Max(1, a + d);
        int useTo = Math.Min(bigN, b - d);
        int n = y1.Length;
        var result = new List<int>();
        if (n < 2 * kmin)
        {
            // R returns breakpts <- 0 unshifted ("Check that vectors are long enough to run algorithm").
            result.Add(0);
            return result;
        }

        double s1Sq = sd1 * sd1, s2Sq = sd2 * sd2;

        // 1-based arrays of length n + 1 (index 0 unused), mirroring the R code.
        double initSum1 = 0.0, initKvad1 = 0.0, initSum2 = 0.0, initKvad2 = 0.0;
        for (int i = 1; i <= kmin; i++)
        {
            initSum1 += y1[i - 1];
            initKvad1 += y1[i - 1] * y1[i - 1];
            initSum2 += y2[i - 1];
            initKvad2 += y2[i - 1] * y2[i - 1];
        }

        double initAve1 = initSum1 / kmin, initAve2 = initSum2 / kmin;
        var bestCost = new double[n + 1];
        var bestSplit = new int[n + 1];
        bestCost[kmin] = (initKvad1 - initSum1 * initAve1) / s1Sq + (initKvad2 - initSum2 * initAve2) / s2Sq;

        var sum1 = new double[n + 1];
        var sum2 = new double[n + 1];
        var kvad1 = new double[n + 1];
        var kvad2 = new double[n + 1];
        var aver1 = new double[n + 1];
        var aver2 = new double[n + 1];
        var cost = new double[n + 1];
        int kminP1 = kmin + 1;

        for (int k = kminP1; k <= 2 * kmin - 1; k++)
        {
            double yk1 = y1[k - 1], yk2 = y2[k - 1];
            for (int t = kminP1; t <= k; t++)
            {
                sum1[t] += yk1;
                aver1[t] = sum1[t] / (k - t + 1);
                kvad1[t] += yk1 * yk1;
                sum2[t] += yk2;
                aver2[t] = sum2[t] / (k - t + 1);
                kvad2[t] += yk2 * yk2;
            }

            double bestAver1 = (initSum1 + sum1[kminP1]) / k;
            double bestAver2 = (initSum2 + sum2[kminP1]) / k;
            double cost1 = ((initKvad1 + kvad1[kminP1]) - k * bestAver1 * bestAver1) / s1Sq;
            double cost2 = ((initKvad2 + kvad2[kminP1]) - k * bestAver2 * bestAver2) / s2Sq;
            bestCost[k] = cost1 + cost2;
        }

        for (int m = 2 * kmin; m <= n; m++)
        {
            int nMkminP1 = m - kmin + 1;
            double ym1 = y1[m - 1], ym2 = y2[m - 1];
            for (int t = kminP1; t <= m; t++)
            {
                sum1[t] += ym1;
                aver1[t] = sum1[t] / (m - t + 1);
                kvad1[t] += ym1 * ym1;
                sum2[t] += ym2;
                aver2[t] = sum2[t] / (m - t + 1);
                kvad2[t] += ym2 * ym2;
            }

            // Cost[kminP1:nMkminP1] <- bestCost[kmin:(n-kmin)] + cost1 + cost2; Pos <- which.min(...) + kmin
            int pos = -1;
            double minCost = double.PositiveInfinity;
            for (int t = kminP1; t <= nMkminP1; t++)
            {
                double c1 = (kvad1[t] - sum1[t] * aver1[t]) / s1Sq;
                double c2 = (kvad2[t] - sum2[t] * aver2[t]) / s2Sq;
                cost[t] = bestCost[t - 1] + c1 + c2;
                if (pos < 0 || cost[t] < minCost)
                {
                    minCost = cost[t];
                    pos = t;
                }
            }

            double best = cost[pos] + gamma;
            double totAver1 = (sum1[kminP1] + initSum1) / m;
            double totCost1 = ((kvad1[kminP1] + initKvad1) - m * totAver1 * totAver1) / s1Sq;
            double totAver2 = (sum2[kminP1] + initSum2) / m;
            double totCost2 = ((kvad2[kminP1] + initKvad2) - m * totAver2 * totAver2) / s2Sq;
            double totCost = totCost1 + totCost2;
            if (totCost < best)
            {
                pos = 1;
                best = totCost;
            }

            bestCost[m] = best;
            bestSplit[m] = pos - 1;
        }

        // Trace back: breakpts <- c(bestSplit[n], breakpts) until n == 0 (the leading 0 is included).
        var trace = new List<int> { n };
        int cursor = n;
        while (cursor > 0)
        {
            cursor = bestSplit[cursor];
            trace.Insert(0, cursor);
        }

        foreach (int bp in trace)
        {
            int global = bp + from - 1;
            if (global >= useFrom && global <= useTo)
            {
                result.Add(global);
            }
        }

        return result;
    }

    /// <summary>
    /// ASCAT <c>getMad(x, k)</c>: zeros removed (likely imputed), then <c>mad(x − medianFilter(x, k))</c>. NaN when no
    /// value remains (R <c>mad(numeric(0))</c> = NA).
    /// </summary>
    private static double GetMad(double[] x, int k)
    {
        var nonZero = new List<double>(x.Length);
        foreach (double v in x)
        {
            if (v != 0.0)
            {
                nonZero.Add(v);
            }
        }

        if (nonZero.Count == 0)
        {
            return double.NaN;
        }

        double[] values = nonZero.ToArray();
        double[] runMedian = MedianFilter(values, k);
        var dif = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            dif[i] = values[i] - runMedian[i];
        }

        return RMad(dif);
    }

    /// <summary>ASCAT <c>madWins(x, tau, k)</c>: <c>xhat = medianFilter(x, k)</c>; residuals clipped at ±tau·mad(x − xhat).</summary>
    private static double[] MadWinsorize(double[] x, double tau, int k)
    {
        if (x.Length == 0)
        {
            return Array.Empty<double>();
        }

        double[] xhat = MedianFilter(x, k);
        var d = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            d[i] = x[i] - xhat[i];
        }

        double z = tau * RMad(d);
        var win = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            win[i] = xhat[i] + Math.Clamp(d[i], -z, z); // psi(d, z)
        }

        return win;
    }

    /// <summary>
    /// ASCAT's NA-aware winsorisation (ascat.aspcf <c>lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)], 2.5, 25)$ywin</c>;
    /// ascat.asmultipcf <c>madWinsMatrixWithNA</c> per column): <see cref="MadWinsorize"/> of the non-NaN values in order,
    /// NaN kept at the missing positions. Equals <see cref="MadWinsorize"/> on complete data.
    /// </summary>
    private static double[] MadWinsorizeWithNa(double[] x, double tau, int k)
    {
        int present = 0;
        foreach (double v in x)
        {
            if (!double.IsNaN(v))
            {
                present++;
            }
        }

        if (present == x.Length)
        {
            return MadWinsorize(x, tau, k);
        }

        var values = new double[present];
        for (int i = 0, j = 0; i < x.Length; i++)
        {
            if (!double.IsNaN(x[i]))
            {
                values[j++] = x[i];
            }
        }

        double[] wins = MadWinsorize(values, tau, k);
        var result = new double[x.Length];
        for (int i = 0, j = 0; i < x.Length; i++)
        {
            result[i] = double.IsNaN(x[i]) ? double.NaN : wins[j++];
        }

        return result;
    }

    /// <summary>
    /// R <c>mad(x)</c> = 1.4826 · median(|x − median(x)|) (constant <see cref="MadConsistencyConstant"/>, raw MAD from
    /// <see cref="RawMedianAbsoluteDeviation"/>, both shared with the MATH score).
    /// </summary>
    private static double RMad(double[] x) =>
        MadConsistencyConstant * RawMedianAbsoluteDeviation(x, StatisticsHelper.Median(x));

    /// <summary>
    /// Unscaled median absolute deviation about <paramref name="center"/>: median(|xᵢ − center|) (canonical
    /// <see cref="StatisticsHelper.Median"/>). Shared by R <c>mad</c> (<see cref="RMad"/>, ASCAT ASPCF) and the maftools
    /// MATH score (<see cref="CalculateITH(IReadOnlyList{double})"/>), which apply the 1.4826 factor in different operation orders.
    /// </summary>
    private static double RawMedianAbsoluteDeviation(double[] x, double center)
    {
        var dev = new double[x.Length];
        for (int i = 0; i < x.Length; i++)
        {
            dev[i] = Math.Abs(x[i] - center);
        }

        return StatisticsHelper.Median(dev);
    }

    /// <summary>
    /// ASCAT <c>medianFilter(x, k)</c>: <c>runmed(x, 2k + 1, endrule = "median")</c>, the width reduced to n (odd n) or
    /// n − 1 (even n) when it exceeds the series length.
    /// </summary>
    private static double[] MedianFilter(double[] x, int k)
    {
        int n = x.Length;
        int width = 2 * k + 1;
        if (width > n)
        {
            if (n == 0)
            {
                width = 1;
            }
            else
            {
                width = n % 2 == 0 ? n - 1 : n;
            }
        }

        return RunningMedian(x, width);
    }

    /// <summary>
    /// R <c>runmed(x, k, endrule = "median")</c> for odd k ≤ n: centred running medians of width k in the interior,
    /// the first/last k%/%2 values kept, then <c>smoothEnds(res, k)</c> (Tukey's end-point rule).
    /// </summary>
    private static double[] RunningMedian(double[] x, int k)
    {
        int n = x.Length;
        var res = (double[])x.Clone();
        int half = k / 2;
        if (half < 1)
        {
            return res;
        }

        var window = new double[k];
        for (int i = half; i < n - half; i++)
        {
            Array.Copy(x, i - half, window, 0, k);
            res[i] = StatisticsHelper.Median(window);
        }

        return SmoothEnds(res, k);
    }

    /// <summary>R <c>smoothEnds(y, k)</c> (stats package), verbatim.</summary>
    private static double[] SmoothEnds(double[] y, int k)
    {
        int half = k / 2;
        if (half < 1)
        {
            return y;
        }

        int n = y.Length;
        var sm = (double[])y.Clone();
        // 0-based: R index i ↔ C# index i − 1.
        if (half >= 2)
        {
            sm[1] = Med3(y[0], y[1], y[2]);
            sm[n - 2] = Med3(y[n - 1], y[n - 2], y[n - 3]);
            for (int i = 3; i <= half; i++)
            {
                int j = 2 * i - 1;
                sm[i - 1] = MedianOfRange(y, 0, j);
                sm[n - i] = MedianOfRange(y, n - j, j);
            }
        }

        sm[0] = Med3(y[0], sm[1], sm[1] - 2.0 * (sm[2] - sm[1]));
        sm[n - 1] = Med3(y[n - 1], sm[n - 2], sm[n - 2] - 2.0 * (sm[n - 3] - sm[n - 2]));
        return sm;
    }

    /// <summary>Median of y[start .. start + count) (odd count, R <c>med.odd</c>).</summary>
    private static double MedianOfRange(double[] y, int start, int count)
    {
        var tmp = new double[count];
        Array.Copy(y, start, tmp, 0, count);
        return StatisticsHelper.Median(tmp);
    }

    /// <summary>R <c>smoothEnds</c> helper <c>med3(a, b, c)</c>, verbatim.</summary>
    private static double Med3(double a, double b, double c)
    {
        double m = b;
        if (a < b)
        {
            if (c < b)
            {
                m = a >= c ? a : c;
            }
        }
        else
        {
            if (c > b)
            {
                m = a <= c ? a : c;
            }
        }

        return m;
    }

    // ---- ASCAT ascat.aspcf with germline-homozygous probes (B24 F36) ----

    /// <summary>ASCAT <c>predictGermlineHomozygousStretches</c>: hard-coded p-value 0.001 of a homozygous run.</summary>
    private const double AscatHomozygousStretchPValue = 0.001;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: PCF context of ±100 probes around the stretch.</summary>
    private const int AscatHomStretchContext = 100;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: replacement margin of ±5 probes around the stretch.</summary>
    private const int AscatHomStretchMargin = 5;

    /// <summary>ascat.aspcf homozygous-stretch resegmentation: a probe is replaced when |Δ level| &gt; 0.3 …</summary>
    private const double AscatHomStretchLevelDifference = 0.3;

    /// <summary>… and only when more than 5 probes of the stretch differ (<c>sum(dif &gt; 0.3) &gt; 5</c>).</summary>
    private const int AscatHomStretchMinDifferingProbes = 5;

    /// <summary>
    /// One ASCAT logR segment of a germline-aware ASPCF segmentation (<see cref="AspcfSegmentation"/>): a maximal run of
    /// loci with the same segmented logR inside one chromosome — exactly the segments runASCAT builds from
    /// <c>rle(lrrsegmented)</c> ∪ chromosome ends (ascat.runAscat.R, <c>tlrstart</c>/<c>tlrend</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Start">Position of the first locus of the segment.</param>
    /// <param name="End">Position of the last locus of the segment.</param>
    /// <param name="MeanLogR">Segmented logR level (mean raw logR of the segment's loci, heterozygous and homozygous).</param>
    /// <param name="MeanBAF">Mirrored (≥ 0.5) segmented BAF of the segment's <b>first</b> germline-heterozygous locus
    /// (runASCAT <c>bafke = bafsegmented[bafpos][1]</c>); <see cref="double.NaN"/> when the segment has no heterozygous
    /// locus (a germline-homozygous stretch: ASCAT then sets <c>bafke = 0</c> and only nA + nB is meaningful).</param>
    /// <param name="LocusCount">Number of loci (all germline genotypes) in the segment.</param>
    /// <param name="HeterozygousLocusCount">Number of germline-heterozygous loci in the segment.</param>
    public readonly record struct AspcfSegment(
        string Chromosome,
        long Start,
        long End,
        double MeanLogR,
        double MeanBAF,
        int LocusCount,
        int HeterozygousLocusCount)
    {
        /// <summary>True when the segment has a segmented BAF (at least one germline-heterozygous locus).</summary>
        public bool HasBaf => !double.IsNaN(MeanBAF);
    }

    /// <summary>
    /// Result of ASCAT allele-specific segmentation with germline genotypes
    /// (<see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/>): the
    /// per-locus ASCAT <c>Tumor_LogR_segmented</c> (every locus) and <c>Tumor_BAF_segmented</c> (heterozygous loci only),
    /// plus the logR segments runASCAT builds from them. Consumed by <see cref="FitPurityPloidyFromAspcf"/>.
    /// </summary>
    public sealed class AspcfSegmentation
    {
        internal AspcfSegmentation(
            AlleleSpecificLocus[] loci, bool[] heterozygous, double[] segmentedLogR, double[] segmentedBaf,
            AspcfSegment[] segments)
        {
            Loci = Array.AsReadOnly(loci);
            GermlineHeterozygous = Array.AsReadOnly(heterozygous);
            SegmentedLogR = Array.AsReadOnly(segmentedLogR);
            SegmentedBaf = Array.AsReadOnly(segmentedBaf);
            Segments = Array.AsReadOnly(segments);
        }

        /// <summary>The input loci, in input order.</summary>
        public IReadOnlyList<AlleleSpecificLocus> Loci { get; }

        /// <summary>The germline genotype of each locus (true = heterozygous), in input order.</summary>
        public IReadOnlyList<bool> GermlineHeterozygous { get; }

        /// <summary>ASCAT <c>Tumor_LogR_segmented</c>: the segmented logR level of every locus.</summary>
        public IReadOnlyList<double> SegmentedLogR { get; }

        /// <summary>Mirrored (≥ 0.5) segmented BAF of each heterozygous locus (ASCAT <c>Tumor_BAF_segmented = 1 − value</c>);
        /// <see cref="double.NaN"/> at germline-homozygous loci (ASCAT segments BAF on heterozygous probes only) and at
        /// heterozygous loci that ASCAT leaves out of <c>Tumor_BAF_segmented</c> because of missing data (B24 F59).</summary>
        public IReadOnlyList<double> SegmentedBaf { get; }

        /// <summary>The logR segments (runs of equal segmented logR within a chromosome), in input order.</summary>
        public IReadOnlyList<AspcfSegment> Segments { get; }
    }

    /// <summary>
    /// ASCAT allele-specific segmentation with germline genotypes — a port of the complete <c>ascat.aspcf</c>
    /// (VanLoo-lab/ascat ascat.aspcf.R; Ross et al. 2021, <i>Bioinformatics</i> 37:1909) for mixed heterozygous /
    /// homozygous input, as ASCAT runs it with <c>ascat.gg$germlinegenotypes</c> (TRUE = homozygous):
    /// <list type="number">
    /// <item>Germline-homozygous stretches (<c>predictGermlineHomozygousStretches</c>): with h = the fraction of homozygous
    /// loci, a run of ≥ <c>ceiling(log(0.001, h))</c> consecutive homozygous loci is a stretch.</item>
    /// <item>Per chromosome: logR of <b>all</b> loci is MAD-winsorised; BAF of the heterozygous loci only. Each
    /// heterozygous locus gets the mean winsorised logR of the loci between the midpoints to its heterozygous
    /// neighbours; ASPCF (<c>fastAspcf</c>, kmin 6; &lt; 6 heterozygous loci ⇒ one segment) runs on these averages and
    /// the heterozygous BAF. Between consecutive heterozygous loci of different levels the logR breakpoint is placed
    /// in the homozygous gap at the minimum absolute deviation (ascat.aspcf "find best breakpoint", verbatim
    /// including its <c>1:0</c> indexing at bp = 0); levels are re-estimated as the mean raw logR. A chromosome without
    /// heterozygous loci becomes one logR segment (mean raw logR) with no BAF.</item>
    /// <item>Each homozygous stretch is re-segmented: exact PCF (<c>exactPcf</c>, kmin 6, penalty floor(γ/4)) of the
    /// winsorised raw logR over the stretch ± 100 loci; the stretch ± 5 loci takes the PCF level where it differs by
    /// more than 0.3, if more than 5 loci differ.</item>
    /// <item>Genome-wide: zero levels are filled from their neighbours (<c>fillNA(zeroIsNA = TRUE)</c>) and every run of
    /// equal levels is re-estimated as its mean raw logR; the penalty ladder (35, 50, 70, 100, 140) is climbed while
    /// ≥ 800 distinct levels remain.</item>
    /// </list>
    /// With every locus heterozygous the logR levels and BAF are those of
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, double)"/> (which is unchanged), apart from
    /// two R genome-wide steps that overload omits — a level exactly 0 is replaced by <c>fillNA</c>, and equal adjacent
    /// levels on consecutive chromosomes are re-averaged together.
    /// Loci of one chromosome must be contiguous in the input (each contiguous same-label run is one ASCAT <c>chr</c>
    /// part); positions are not re-sorted.
    /// </summary>
    /// <param name="loci">Per-locus measurements in genome order. LogR must be finite or NaN; BAF must be in [0, 1] or NaN
    /// at heterozygous loci and is ignored at homozygous loci (it may be NaN there, e.g. a copy-number-only probe). NaN is
    /// R's <c>NA</c> and follows ascat.aspcf's NA path (B24 F59): a heterozygous locus with a missing BAF or logR is left
    /// out of the ASPCF (<c>Select_het</c>) and gets no segmented BAF; missing logR is skipped in winsorisation, averages,
    /// the gap breakpoint and level means; all-missing levels are filled by <c>fillNA</c> / the previous level.</param>
    /// <param name="germlineHeterozygous">Germline genotype per locus: true = heterozygous (ASCAT
    /// <c>germlinegenotypes == FALSE</c>), false = homozygous. Same length as <paramref name="loci"/>.</param>
    /// <param name="penalty">ASPCF penalty (ASCAT <c>penalty</c>, default <see cref="AspcfDefaultPenalty"/> = 70).</param>
    /// <returns>The per-locus segmented logR/BAF and the runASCAT logR segments.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Lengths differ, a chromosome label is null, a logR is infinite, or a
    /// heterozygous locus has a BAF outside [0, 1] (and not NaN).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static AspcfSegmentation SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        IReadOnlyList<bool> germlineHeterozygous,
        double penalty = AspcfDefaultPenalty)
        => SegmentAlleleSpecificAspcfCore(loci, germlineHeterozygous, null, penalty);

    /// <summary>
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/> with ASCAT's
    /// male X non-PAR germline re-genotyping (<c>ascat.aspcf(…, seed)</c> with <c>gender = "XY"</c> and <c>X_nonPAR</c>,
    /// B24 F58; see <see cref="AscatMaleXGenotyping"/>). As in ascat.aspcf, the germline-homozygous stretches are predicted
    /// from the caller's genotypes <i>before</i> the re-genotyping, and the segmentation then uses the re-genotyped flags
    /// (returned in <see cref="AspcfSegmentation.GermlineHeterozygous"/>; a re-genotyped heterozygous locus needs a BAF in
    /// [0, 1]). With an inactive <paramref name="maleX"/> (female, or no X non-PAR interval) the result equals the
    /// three-argument overload.
    /// </summary>
    /// <param name="loci">See the three-argument overload.</param>
    /// <param name="germlineHeterozygous">The caller's germline genotypes (true = heterozygous), before re-genotyping.</param>
    /// <param name="maleX">Sex model, seed and optional germline BAF.</param>
    /// <param name="penalty">ASPCF penalty (default <see cref="AspcfDefaultPenalty"/> = 70).</param>
    /// <returns>The per-locus segmented logR/BAF (with the re-genotyped genotypes) and the runASCAT logR segments.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">As the three-argument overload; also a germline-BAF count other than one per
    /// locus, or (without germline BAF) no heterozygous autosomal locus while re-genotyping (R's <c>sample()</c> fails).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static AspcfSegmentation SegmentAlleleSpecificAspcf(
        IEnumerable<AlleleSpecificLocus> loci,
        IReadOnlyList<bool> germlineHeterozygous,
        AscatMaleXGenotyping maleX,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(maleX);
        return SegmentAlleleSpecificAspcfCore(loci, germlineHeterozygous, maleX, penalty);
    }

    /// <summary>One sample of a multi-sample <c>ascat.aspcf</c> run (<see cref="SegmentAlleleSpecificAspcfSamples"/>).</summary>
    /// <param name="Loci">The sample's per-locus logR/BAF in genome order (as the single-sample overloads).</param>
    /// <param name="GermlineHeterozygous">The caller's germline genotypes (true = heterozygous), before re-genotyping.</param>
    public sealed record AspcfSampleInput(IReadOnlyList<AlleleSpecificLocus> Loci, IReadOnlyList<bool> GermlineHeterozygous)
    {
        /// <summary>
        /// The sample's sex model (ASCAT <c>gender[sample]</c> and <c>X_nonPAR</c>); null = no re-genotyping (as a female
        /// sample or a run without <c>X_nonPAR</c>).
        /// </summary>
        public AscatSexModel? SexModel { get; init; }

        /// <summary>The sample's germline BAF per locus (ASCAT <c>Germline_BAF[, sample]</c>), or null (<c>sample()</c> branch).</summary>
        public IReadOnlyList<double>? GermlineBaf { get; init; }
    }

    /// <summary>
    /// <c>ascat.aspcf</c> over several samples in one call (VanLoo-lab/ascat ascat.aspcf.R: <c>set.seed(seed)</c> once,
    /// then <c>for (sample in selectsamples)</c>; B24 F65). Each sample is segmented exactly as
    /// <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, AscatMaleXGenotyping, double)"/>,
    /// but the male X non-PAR re-genotyping (<see cref="AscatMaleXGenotyping"/>) of all samples draws from ONE R random
    /// stream, in sample order: a later male sample continues where the previous one stopped (samples without active
    /// re-genotyping — female, no X non-PAR interval, or ≤ 5 known non-PAR probes — draw nothing). The first male sample
    /// therefore equals the single-sample call with the same seed; later ones generally differ from it.
    /// </summary>
    /// <param name="samples">The samples, in ASCAT's <c>selectsamples</c> order.</param>
    /// <param name="seed">R <c>set.seed</c> value (ASCAT <c>seed</c>); null = <see cref="AscatMaleXGenotyping.DefaultSeed"/>.</param>
    /// <param name="penalty">ASPCF penalty (default <see cref="AspcfDefaultPenalty"/> = 70), shared by all samples.</param>
    /// <returns>One segmentation per sample, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="samples"/> is null.</exception>
    /// <exception cref="ArgumentException">A sample is null, or a sample fails the single-sample contract.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite.</exception>
    public static IReadOnlyList<AspcfSegmentation> SegmentAlleleSpecificAspcfSamples(
        IReadOnlyList<AspcfSampleInput> samples,
        int? seed = null,
        double penalty = AspcfDefaultPenalty)
    {
        ArgumentNullException.ThrowIfNull(samples);
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        int rSeed = seed ?? AscatMaleXGenotyping.DefaultSeed();
        var rng = new RMersenneTwister(rSeed); // ascat.aspcf: set.seed(seed) before the sample loop.
        var results = new AspcfSegmentation[samples.Count];
        for (int s = 0; s < samples.Count; s++)
        {
            AspcfSampleInput sample = samples[s]
                ?? throw new ArgumentException($"Sample {s} is null.", nameof(samples));
            AscatMaleXGenotyping? maleX = sample.SexModel is null
                ? null
                : new AscatMaleXGenotyping(sample.SexModel, rSeed, sample.GermlineBaf);
            results[s] = SegmentAlleleSpecificAspcfCore(sample.Loci, sample.GermlineHeterozygous, maleX, penalty, rng);
        }

        return results;
    }

    private static AspcfSegmentation SegmentAlleleSpecificAspcfCore(
        IEnumerable<AlleleSpecificLocus> loci,
        IReadOnlyList<bool> germlineHeterozygous,
        AscatMaleXGenotyping? maleX,
        double penalty,
        RMersenneTwister? rng = null)
    {
        ArgumentNullException.ThrowIfNull(loci);
        ArgumentNullException.ThrowIfNull(germlineHeterozygous);
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        AlleleSpecificLocus[] probes = loci.ToArray();
        if (probes.Length != germlineHeterozygous.Count)
        {
            throw new ArgumentException(
                "Exactly one germline genotype is required per locus.", nameof(germlineHeterozygous));
        }

        var het = new bool[probes.Length];
        for (int i = 0; i < probes.Length; i++)
        {
            het[i] = germlineHeterozygous[i];
            if (probes[i].Chromosome is null)
            {
                throw new ArgumentException("A locus has a null chromosome label.", nameof(loci));
            }

            if (double.IsInfinity(probes[i].LogR))
            {
                throw new ArgumentException("A logR must be finite or NaN (missing).", nameof(loci));
            }
        }

        int n = probes.Length;
        var logR = new double[n];
        for (int i = 0; i < n; i++)
        {
            logR[i] = probes[i].LogR;
        }

        List<(int Lo, int Hi)> runs = ContiguousChromosomeRuns(probes);

        // ascat.aspcf: ghs = predictGermlineHomozygousStretches(chr, gg) precedes the male X non-PAR re-genotyping.
        List<(int Run, int Start, int End)> stretches = PredictGermlineHomozygousStretches(runs, het);
        RegenotypeMaleXNonPar(probes, het, maleX, nameof(maleX), rng);
        for (int i = 0; i < n; i++)
        {
            if (het[i] && !IsAscatBafOrMissing(probes[i].BAF))
            {
                throw new ArgumentException("Every heterozygous locus needs a BAF in [0, 1] or NaN (missing).", nameof(loci));
            }
        }

        var ladder = new List<double> { penalty };
        foreach (double p in AspcfPenaltyLadder)
        {
            if (p > penalty)
            {
                ladder.Add(p);
            }
        }

        double[] logRPcfed = Array.Empty<double>();
        double[] bafPcfed = new double[n];
        foreach (double segmentPenalty in ladder)
        {
            logRPcfed = new double[n];
            Array.Fill(bafPcfed, double.NaN);
            for (int r = 0; r < runs.Count; r++)
            {
                SegmentChromosomeAspcfGermline(logR, probes, het, runs[r].Lo, runs[r].Hi, segmentPenalty, logRPcfed, bafPcfed);
                ResegmentHomozygousStretches(logR, runs[r].Lo, runs[r].Hi, stretches, r, segmentPenalty, logRPcfed);
            }

            AscatFillNa(logRPcfed);
            ReadaptLevels(logRPcfed, logR, 0, n);

            if (logRPcfed.Distinct().Count() < AspcfMaxSegmentLevels)
            {
                break;
            }
        }

        return new AspcfSegmentation(probes, het, logRPcfed, bafPcfed, BuildAspcfSegments(probes, het, runs, logRPcfed, bafPcfed));
    }

    /// <summary>Contiguous same-chromosome runs of the loci (ASCAT <c>chr</c> parts), in input order.</summary>
    private static List<(int Lo, int Hi)> ContiguousChromosomeRuns(AlleleSpecificLocus[] loci)
    {
        var runs = new List<(int Lo, int Hi)>();
        int start = 0;
        while (start < loci.Length)
        {
            int end = start;
            while (end + 1 < loci.Length && loci[end + 1].Chromosome == loci[start].Chromosome)
            {
                end++;
            }

            runs.Add((start, end));
            start = end + 1;
        }

        return runs;
    }

    /// <summary>
    /// ASCAT <c>predictGermlineHomozygousStretches</c>: runs of ≥ <c>ceiling(log(0.001, h))</c> consecutive homozygous loci
    /// per chromosome, h = genome-wide homozygous fraction (R <c>log(x, base) = log(x)/log(base)</c>). Returns 0-based
    /// global (first, last) indices.
    /// </summary>
    private static List<(int Run, int Start, int End)> PredictGermlineHomozygousStretches(
        List<(int Lo, int Hi)> runs, bool[] het)
    {
        var stretches = new List<(int Run, int Start, int End)>();
        if (het.Length == 0)
        {
            return stretches;
        }

        int homCount = 0;
        foreach (bool h in het)
        {
            if (!h)
            {
                homCount++;
            }
        }

        double fractionHom = (double)homCount / het.Length;
        double threshold = Math.Ceiling(Math.Log(AscatHomozygousStretchPValue) / Math.Log(fractionHom));
        for (int r = 0; r < runs.Count; r++)
        {
            int count = 0, first = -1;
            for (int i = runs[r].Lo; i <= runs[r].Hi; i++)
            {
                if (!het[i])
                {
                    if (count == 0)
                    {
                        first = i;
                    }

                    count++;
                }
                else
                {
                    // An empty run gives R an NA row (min of an empty vector), which ascat.aspcf skips.
                    if (count > 0 && count >= threshold)
                    {
                        stretches.Add((r, first, i - 1));
                    }

                    count = 0;
                }
            }

            if (count > 0 && count >= threshold)
            {
                stretches.Add((r, first, runs[r].Hi));
            }
        }

        return stretches;
    }

    /// <summary>
    /// The per-chromosome body of <c>ascat.aspcf</c> with germline genotypes, on loci [lo, hi]: writes the segmented logR
    /// of every locus to <paramref name="logRPcfed"/> and the mirrored segmented BAF of each heterozygous locus to
    /// <paramref name="bafPcfed"/>.
    /// </summary>
    private static void SegmentChromosomeAspcfGermline(
        double[] logR, AlleleSpecificLocus[] loci, bool[] het, int lo, int hi, double penalty,
        double[] logRPcfed, double[] bafPcfed)
    {
        int len = hi - lo + 1;
        var lr = new double[len];
        Array.Copy(logR, lo, lr, 0, len);
        double[] lrWins = MadWinsorizeWithNa(lr, AspcfWinsorTau, AspcfMedianHalfWindow);

        // R "indices" = which(Select_het), Select_het = !homo & !is.na(baf) & !is.na(lr) (1-based local indices).
        var indices = new List<int>();
        for (int i = 0; i < len; i++)
        {
            if (het[lo + i] && !double.IsNaN(loci[lo + i].BAF) && !double.IsNaN(lr[i]))
            {
                indices.Add(i + 1);
            }
        }

        int h = indices.Count;
        if (h == 0)
        {
            // No usable heterozygous probe: a single logR segment (mean(lr, na.rm = TRUE)), no BAF.
            double level = MeanNaRm(lr, 0, len);
            Array.Fill(logRPcfed, level, lo, len);
            return;
        }

        var bafSel = new double[h];
        var mirrored = new double[h];
        for (int k = 0; k < h; k++)
        {
            bafSel[k] = loci[lo + indices[k] - 1].BAF;
            mirrored[k] = MirrorBaf(bafSel[k]);
        }

        double[] bafWinsMirrored = MadWinsorize(mirrored, AspcfWinsorTau, AspcfMedianHalfWindow);
        var bafWins = new double[h];
        for (int k = 0; k < h; k++)
        {
            bafWins[k] = OrientBaf(bafSel[k], bafWinsMirrored[k]);
        }

        // averageIndices = c(1, (indices[-h] + indices[-1])/2, length(lr) + 0.01); start = ceiling, end = floor(· − 0.01).
        var localIndices = indices.ConvertAll(i => i - 1);
        var logRAveraged = new double[h];
        for (int k = 0; k < h; k++)
        {
            (int startIndex, int endIndex) = AscatAveragingWindow(localIndices, k, len);
            logRAveraged[k] = RColonMeanNaRm(lrWins, startIndex, endIndex);
        }

        var levelPerHet = new double[h];
        var bafPerHet = new double[h];
        if (h < AspcfMinSegmentLength)
        {
            Array.Fill(levelPerHet, Mean(logRAveraged, 0, h));
            Array.Fill(bafPerHet, Mean(bafWinsMirrored, 0, h));
        }
        else
        {
            (int[] breakpoints, double[] segmentBaf) = FastAspcf(logRAveraged, bafWins, AspcfMinSegmentLength, penalty);
            for (int s = 0; s + 1 < breakpoints.Length; s++)
            {
                double level = Mean(logRAveraged, breakpoints[s], breakpoints[s + 1]); // fastAspcf yhat1
                for (int k = breakpoints[s]; k < breakpoints[s + 1]; k++)
                {
                    levelPerHet[k] = level;
                    bafPerHet[k] = segmentBaf[s];
                }
            }
        }

        // logRc: extend the heterozygous levels over the homozygous gaps ("find best breakpoint").
        var logRc = new double[len];
        for (int p = 0; p < h; p++)
        {
            int at = indices[p]; // 1-based
            if (p == 0)
            {
                Array.Fill(logRc, levelPerHet[0], 0, at);
            }

            if (p == h - 1)
            {
                Array.Fill(logRc, levelPerHet[p], at, len - at);
            }
            else if (levelPerHet[p] == levelPerHet[p + 1])
            {
                Array.Fill(logRc, levelPerHet[p], at, indices[p + 1] - at);
            }
            else
            {
                int total = indices[p + 1] - at;
                int breakpoint = AscatBestGapBreakpoint(lr, at, total, levelPerHet[p], levelPerHet[p + 1]);
                Array.Fill(logRc, levelPerHet[p], at, breakpoint);
                Array.Fill(logRc, levelPerHet[p + 1], at + breakpoint, total - breakpoint);
            }
        }

        // 2nd step: adapt levels (rle(logRc) runs → mean raw logR).
        int runStart = 0;
        for (int i = 1; i <= len; i++)
        {
            if (i == len || logRc[i] != logRc[runStart])
            {
                double level = MeanNaRm(lr, runStart, i);
                Array.Fill(logRPcfed, level, lo + runStart, i - runStart);
                runStart = i;
            }
        }

        for (int k = 0; k < h; k++)
        {
            bafPcfed[lo + indices[k] - 1] = bafPerHet[k];
        }
    }

    /// <summary>
    /// R <c>mean(x[start:end])</c> for 1-based <paramref name="start"/>, <paramref name="end"/>. When end &lt; start, R's
    /// <c>:</c> counts down (<c>x[n:(n-1)]</c> = two elements) — this happens in ascat.aspcf / ascat.asmultipcf for the last
    /// averaging window, because <c>floor(n + 0.01 − 0.01)</c> is n − 1 in IEEE doubles for n = 2, 32, 128, 16384, 65536, ….
    /// Missing (NaN) values are skipped (<c>na.rm = TRUE</c>, as both R functions average <c>logRaveraged</c>).
    /// </summary>
    private static double RColonMeanNaRm(double[] x, int start, int end) =>
        end >= start ? MeanNaRm(x, start - 1, end) : MeanNaRm(x, end - 1, start);

    /// <summary>
    /// ascat.aspcf "find best breakpoint" between heterozygous loci at 1-based local positions <paramref name="at"/> and
    /// at + <paramref name="total"/>: for bp ∈ 0..total−1, d(bp) = Σ|lr[(1:bp) + at] − left| + Σ|lr[((bp+1):total) + at] −
    /// right|; returns <c>which.min(d) − 1</c>. R's <c>1:0 = c(1, 0)</c> makes bp = 0 compare lr[at + 1] and lr[at] with
    /// the left level — kept verbatim.
    /// </summary>
    private static int AscatBestGapBreakpoint(double[] lr, int at, int total, double left, double right)
    {
        // sum(abs(lr[...] − level), na.rm = TRUE): a missing logR contributes nothing.
        static double AbsDeviationNaRm(double x, double level) => double.IsNaN(x) ? 0.0 : Math.Abs(x - level);

        int best = 0;
        double bestDistance = double.NaN;
        for (int bp = 0; bp < total; bp++)
        {
            double leftSum = 0.0;
            if (bp == 0)
            {
                leftSum += AbsDeviationNaRm(lr[at], left);     // lr[1 + at] (R 1-based) = C# lr[at]
                leftSum += AbsDeviationNaRm(lr[at - 1], left); // lr[0 + at]
            }
            else
            {
                for (int q = 1; q <= bp; q++)
                {
                    leftSum += AbsDeviationNaRm(lr[q + at - 1], left);
                }
            }

            double rightSum = 0.0;
            for (int q = bp + 1; q <= total; q++)
            {
                rightSum += AbsDeviationNaRm(lr[q + at - 1], right);
            }

            double distance = leftSum + rightSum;
            if (bp == 0 || distance < bestDistance)
            {
                bestDistance = distance;
                best = bp;
            }
        }

        return best;
    }

    /// <summary>
    /// ascat.aspcf "correct wrong segments in germline homozygous stretches" for the stretches of chromosome run
    /// <paramref name="run"/> [lo, hi]: exact PCF (kmin 6, penalty floor(γ/4)) of the winsorised raw logR over the stretch
    /// ± 100 loci; the stretch ± 5 loci takes the PCF level where it differs by &gt; 0.3, if &gt; 5 loci differ.
    /// </summary>
    private static void ResegmentHomozygousStretches(
        double[] logR, int lo, int hi, List<(int Run, int Start, int End)> stretches, int run, double penalty,
        double[] logRPcfed)
    {
        double pcfPenalty = Math.Floor(penalty / 4.0);
        foreach ((int stretchRun, int start, int end) in stretches)
        {
            if (stretchRun != run)
            {
                continue;
            }

            int start2 = Math.Max(start - AscatHomStretchContext, lo);
            int end2 = Math.Min(end + AscatHomStretchContext, hi);
            int start3 = Math.Max(start - AscatHomStretchMargin, lo);
            int end3 = Math.Min(end + AscatHomStretchMargin, hi);

            // towins; winsed = madWins(towins[!is.na(towins)]); pcfed[!is.na(towins)] = exactPcf(winsed), 0 elsewhere.
            var window = new double[end2 - start2 + 1];
            Array.Copy(logR, start2, window, 0, window.Length);
            var present = new List<int>(window.Length);
            for (int i = 0; i < window.Length; i++)
            {
                if (!double.IsNaN(window[i]))
                {
                    present.Add(i);
                }
            }

            double[] pcfed;
            if (present.Count == window.Length)
            {
                pcfed = ExactPcf(MadWinsorize(window, AspcfWinsorTau, AspcfMedianHalfWindow), AspcfMinSegmentLength, pcfPenalty);
            }
            else
            {
                var towins = new double[present.Count];
                for (int j = 0; j < towins.Length; j++)
                {
                    towins[j] = window[present[j]];
                }

                double[] fitted = ExactPcf(MadWinsorize(towins, AspcfWinsorTau, AspcfMedianHalfWindow), AspcfMinSegmentLength, pcfPenalty);
                pcfed = new double[window.Length];
                for (int j = 0; j < towins.Length; j++)
                {
                    pcfed[present[j]] = fitted[j];
                }
            }

            int count = end3 - start3 + 1;
            var dif = new double[count];
            int differing = 0;
            bool anyMissing = false;
            for (int i = 0; i < count; i++)
            {
                dif[i] = Math.Abs(pcfed[start3 - start2 + i] - logRPcfed[start3 + i]);
                anyMissing |= double.IsNaN(dif[i]);
                if (dif[i] > AscatHomStretchLevelDifference)
                {
                    differing++;
                }
            }

            if (!anyMissing && differing > AscatHomStretchMinDifferingProbes) // !anyNA(dif) && sum(dif > 0.3) > 5
            {
                for (int i = 0; i < count; i++)
                {
                    if (dif[i] > AscatHomStretchLevelDifference)
                    {
                        logRPcfed[start3 + i] = pcfed[start3 - start2 + i];
                    }
                }
            }
        }
    }

    /// <summary>
    /// ASCAT <c>exactPcf(y, kmin, gamma)</c> (ascat.aspcf.R): exact univariate piecewise-constant fit by Potts filtering,
    /// minimum segment length kmin, penalty γ per discontinuity on the raw (unstandardised) squared error.
    /// </summary>
    private static double[] ExactPcf(double[] y, int kmin, double gamma)
    {
        int n = y.Length;
        var yhat = new double[n];
        if (n < 2 * kmin)
        {
            Array.Fill(yhat, Mean(y, 0, n));
            return yhat;
        }

        // 1-based arrays (index 0 unused), mirroring the R code.
        double initSum = 0.0, initKvad = 0.0;
        for (int i = 1; i <= kmin; i++)
        {
            initSum += y[i - 1];
            initKvad += y[i - 1] * y[i - 1];
        }

        double initAve = initSum / kmin;
        var bestCost = new double[n + 1];
        var bestSplit = new int[n + 1];
        var bestAver = new double[n + 1];
        var sum = new double[n + 1];
        var kvad = new double[n + 1];
        var aver = new double[n + 1];
        var cost = new double[n + 1];
        bestCost[kmin] = initKvad - initSum * initAve;
        bestAver[kmin] = initAve;
        int kminP1 = kmin + 1;
        for (int k = kminP1; k <= 2 * kmin - 1; k++)
        {
            double yk = y[k - 1];
            for (int t = kminP1; t <= k; t++)
            {
                sum[t] += yk;
                aver[t] = sum[t] / (k - t + 1);
                kvad[t] += yk * yk;
            }

            bestAver[k] = (initSum + sum[kminP1]) / k;
            bestCost[k] = (initKvad + kvad[kminP1]) - k * (bestAver[k] * bestAver[k]);
        }

        for (int m = 2 * kmin; m <= n; m++)
        {
            double yn = y[m - 1];
            double yn2 = yn * yn;
            for (int t = kminP1; t <= m; t++)
            {
                sum[t] += yn;
                aver[t] = sum[t] / (m - t + 1);
                kvad[t] += yn2;
            }

            int nMkminP1 = m - kmin + 1;
            int pos = -1;
            for (int t = kminP1; t <= nMkminP1; t++)
            {
                cost[t] = bestCost[t - 1] + kvad[t] - sum[t] * aver[t] + gamma;
                if (pos < 0 || cost[t] < cost[pos])
                {
                    pos = t; // which.min: first minimum
                }
            }

            double best = cost[pos];
            double bestLevel = aver[pos];
            double totAver = (sum[kminP1] + initSum) / m;
            double totCost = (kvad[kminP1] + initKvad) - m * totAver * totAver;
            if (totCost < best)
            {
                pos = 1;
                best = totCost;
                bestLevel = totAver;
            }

            bestCost[m] = best;
            bestAver[m] = bestLevel;
            bestSplit[m] = pos - 1;
        }

        int cursor = n;
        while (cursor > 0)
        {
            for (int i = bestSplit[cursor]; i < cursor; i++)
            {
                yhat[i] = bestAver[cursor];
            }

            cursor = bestSplit[cursor];
        }

        return yhat;
    }

    /// <summary>
    /// ASCAT <c>fillNA(vec, zeroIsNA = TRUE)</c> (ascat.aspcf.R), verbatim: levels exactly 0 become missing; a missing run
    /// at the start/end takes the next/previous value; an interior run takes the previous value up to
    /// <c>start + ceiling(N/2)</c> and the next value after it (for N = 2, 3 R's descending <c>(end+1):end</c> sets only
    /// the last element to the next value).
    /// </summary>
    private static void AscatFillNa(double[] vec)
    {
        int n = vec.Length;
        var missing = new List<int>();
        for (int i = 0; i < n; i++)
        {
            if (vec[i] == 0.0 || double.IsNaN(vec[i]))
            {
                vec[i] = double.NaN;
                missing.Add(i);
            }
        }

        if (missing.Count == 0)
        {
            return;
        }

        var starts = new List<int> { missing[0] };
        var ends = new List<int>();
        for (int k = 1; k < missing.Count; k++)
        {
            if (missing[k] - missing[k - 1] > 1)
            {
                ends.Add(missing[k - 1]);
                starts.Add(missing[k]);
            }
        }

        ends.Add(missing[^1]);
        double At(int index) => index < n ? vec[index] : double.NaN; // R: out-of-range index → NA

        int startAt = 0;
        if (starts[0] == 0)
        {
            double fill = At(ends[0] + 1);
            Array.Fill(vec, fill, 0, ends[0] + 1);
            startAt = 1;
        }

        if (startAt >= starts.Count)
        {
            return;
        }

        int endAt = starts.Count - 1;
        if (double.IsNaN(vec[n - 1]))
        {
            Array.Fill(vec, vec[starts[endAt] - 1], starts[endAt], ends[endAt] - starts[endAt] + 1);
            endAt--;
        }

        for (int k = startAt; k <= endAt; k++)
        {
            int start = starts[k], end = ends[k];
            int count = 1 + end - start;
            if (count == 1)
            {
                vec[start] = vec[start - 1];
                continue;
            }

            int midpoint = start + (count + 1) / 2; // start + ceiling(N/2)
            double previous = vec[start - 1], next = vec[end + 1];
            Array.Fill(vec, previous, start, midpoint - start + 1);
            if (midpoint < end)
            {
                Array.Fill(vec, next, midpoint + 1, end - midpoint);
            }
            else
            {
                vec[end] = next; // vec[(end+1):end] = vec[end+1]
            }
        }
    }

    /// <summary>
    /// ascat.aspcf / ascat.asmultipcf "adapt levels again": every run of equal values of <paramref name="levels"/> in
    /// [from, to) — across chromosome boundaries, as R's genome-wide <c>rle</c> — is set to the mean raw logR of the run
    /// (<c>na.rm = TRUE</c>); a run without any non-missing raw logR takes the previous run's level ("making sure no NA's
    /// get filled in", <c>prevlevel</c> starting at 0). R's <c>rle</c> splits every NaN into its own run where this merges
    /// adjacent NaN levels; both give each such probe the mean of its (all-missing) raw logR, i.e. <c>prevlevel</c>.
    /// </summary>
    private static void ReadaptLevels(double[] levels, double[] rawLogR, int from, int to)
    {
        int runStart = from;
        double previousLevel = 0.0;
        for (int i = from + 1; i <= to; i++)
        {
            if (i == to || !levels[i].Equals(levels[runStart]))
            {
                double level = MeanNaRm(rawLogR, runStart, i);
                if (double.IsNaN(level))
                {
                    level = previousLevel;
                }
                else
                {
                    previousLevel = level;
                }

                Array.Fill(levels, level, runStart, i - runStart);
                runStart = i;
            }
        }
    }

    /// <summary>
    /// runASCAT segments: runs of equal segmented logR split at chromosome ends (<c>union(tlrend, tlrend.chr)</c>),
    /// each carrying the mirrored BAF of its first heterozygous locus with a segmented BAF (NaN if none).
    /// </summary>
    private static AspcfSegment[] BuildAspcfSegments(
        AlleleSpecificLocus[] loci, bool[] het, List<(int Lo, int Hi)> runs, double[] logRPcfed, double[] bafPcfed)
    {
        var segments = new List<AspcfSegment>();
        foreach ((int lo, int hi) in runs)
        {
            int start = lo;
            for (int i = lo + 1; i <= hi + 1; i++)
            {
                if (i == hi + 1 || !logRPcfed[i].Equals(logRPcfed[start]))
                {
                    double baf = double.NaN;
                    int hetCount = 0;
                    for (int k = start; k < i; k++)
                    {
                        if (het[k])
                        {
                            // runASCAT bafke = bafsegmented[bafpos][1]: the first locus that HAS a segmented BAF (a missing
                            // tumour BAF/logR leaves a heterozygous locus out of Tumor_BAF_segmented, B24 F59).
                            if (double.IsNaN(baf))
                            {
                                baf = bafPcfed[k];
                            }

                            hetCount++;
                        }
                    }

                    segments.Add(new AspcfSegment(
                        loci[start].Chromosome, loci[start].Position, loci[i - 1].Position, logRPcfed[start], baf,
                        i - start, hetCount));
                    start = i;
                }
            }
        }

        return segments.ToArray();
    }

    // ---- ASCAT ascat.asmultipcf: joint multi-sample allele-specific PCF (B24 F53) ----

    /// <summary>ascat.asmultipcf penalty ladder <c>unique(c(penalty, 25, 50, 100, 200, 400, 800))</c>, kept where ≥ penalty.</summary>
    private static readonly double[] AsMultiPcfPenaltyLadder = { 25.0, 50.0, 100.0, 200.0, 400.0, 800.0 };

    /// <summary>ascat.asmultipcf with <c>nrow(logRaveraged) &lt; 6</c>: one segment per chromosome part (no PCF).</summary>
    private const int AsMultiPcfMinProbes = 6;

    /// <summary><c>runFastASMultiPCF</c> window size <c>subsize = 5000</c> (windows hold subsize + 1 probes).</summary>
    private const int AsMultiPcfFastSubsize = 5000;

    /// <summary><c>runFastASMultiPCF</c> window step <c>4·subsize/5</c>.</summary>
    private const int AsMultiPcfFastStep = 4 * AsMultiPcfFastSubsize / 5;

    /// <summary>ascat.asmultipcf <c>selectAlg</c>: the joint segmentation algorithm.</summary>
    public enum AsMultiPcfAlgorithm
    {
        /// <summary><c>"exact"</c> (ASCAT default): <c>ASmultiPCFcompact</c> on every probe of the chromosome part.</summary>
        Exact,

        /// <summary><c>"fast"</c>: <c>runFastASMultiPCF</c> — exact segmentation of overlapping 5001-probe windows (step 4000)
        /// collects candidate breakpoints, then the exact algorithm runs on the compacted data. Identical to
        /// <see cref="Exact"/> on parts of ≤ 5001 probes only in the absence of sample weights.</summary>
        Fast,
    }

    /// <summary>Expert options of <see cref="SegmentAlleleSpecificAsMultiPcf"/> (ascat.asmultipcf arguments).</summary>
    public sealed record AsMultiPcfOptions
    {
        /// <summary>The ASCAT defaults: no sample weights, exact algorithm, per-sample refinement on.</summary>
        public static AsMultiPcfOptions Default { get; } = new();

        /// <summary>ASCAT <c>wsample</c>: positive finite weights, either one per sample (applied to its logR and BAF
        /// tracks) or one per track (2·S values: the S logR tracks, then the S BAF tracks). Null = all 1 (default).</summary>
        public IReadOnlyList<double>? SampleWeights { get; init; }

        /// <summary>ASCAT <c>selectAlg</c> (default <see cref="AsMultiPcfAlgorithm.Exact"/>).</summary>
        public AsMultiPcfAlgorithm Algorithm { get; init; } = AsMultiPcfAlgorithm.Exact;

        /// <summary>ASCAT <c>refine</c> (default true): re-segment each sample on the joint breakpoints with penalty γ/S, so a
        /// breakpoint not supported by that sample is removed from it.</summary>
        public bool Refine { get; init; } = true;

        /// <summary>ASCAT's male X non-PAR germline re-genotyping (<c>gender[1] == "XY"</c> with <c>X_nonPAR</c>, <c>seed</c>,
        /// <c>Germline_BAF[, 1]</c>; B24 F58, see <see cref="AscatMaleXGenotyping"/>), applied to the one germline before
        /// segmentation. Null (default) = the caller's genotypes are used as given.</summary>
        public AscatMaleXGenotyping? MaleXGenotyping { get; init; }
    }

    /// <summary>
    /// ASCAT joint multi-sample allele-specific segmentation — a port of <c>ascat.asmultipcf</c> (VanLoo-lab/ascat
    /// ASCAT/R/ascat.asmultipcf.R; Ross et al. 2021, <i>Bioinformatics</i> 37:1909), for samples that share part of their
    /// breakpoints (common ancestry, e.g. multi-region sequencing of one tumour), on one probe set and one germline:
    /// <list type="number">
    /// <item>Per chromosome part and sample: logR of every probe is MAD-winsorised (<c>madWins(x, 2.5, 25)</c>); the BAF of
    /// the germline-heterozygous probes is mirrored to ≥ 0.5, winsorised and mirrored again; homozygous probes have no BAF
    /// (weight 0). Each probe keeps its own winsorised logR (the last probe averages <c>lr[n:floor(n + 0.01 − 0.01)]</c>,
    /// i.e. probes n and n − 1 for n = 2, 32, 128, …, as R does).</item>
    /// <item>Parts with fewer than 6 probes form one segment per sample (logR = mean, BAF = mean mirrored BAF).</item>
    /// <item>Otherwise the 2·S tracks (S logR, S mirrored BAF) are segmented jointly by <c>ASmultiPCFcompact</c>: exact
    /// dynamic programming minimising <c>Σ_segments Σ_tracks −(Σ w·y)²/Σ w + γ·(#segments − 1)</c> on the raw
    /// (unstandardised) values, track weights = 0 for missing BAF times <see cref="AsMultiPcfOptions.SampleWeights"/>;
    /// <see cref="AsMultiPcfAlgorithm.Fast"/> runs <c>runFastASMultiPCF</c> instead.</item>
    /// <item>Refinement (<see cref="AsMultiPcfOptions.Refine"/>): each sample's two tracks are compacted on the joint
    /// breakpoints and re-segmented exactly with penalty γ/S; if that removes segments, the sample takes the refined
    /// segmentation.</item>
    /// <item>BAF per segment (runs of equal value over the heterozygous probes): 0.5 + |b − 0.5|, shrunk to 0.5 when
    /// <c>sqrt(sd² + μ²) &lt; 2·sd</c>, sd = <c>getMadwithNA</c> of the sample's winsorised heterozygous BAF.</item>
    /// <item>LogR levels exactly 0 are replaced by the closest non-zero level of the part; then, genome-wide per sample,
    /// every run of equal levels takes the mean raw logR; the ladder (penalty, 25, 50, 100, 200, 400, 800 — values ≥
    /// penalty) is climbed while any sample has ≥ 800 distinct levels.</item>
    /// </list>
    /// Unlike <c>ascat.aspcf</c>, homozygous stretches are not re-segmented (they take part in the joint logR segmentation
    /// from the start). The male-only <c>X_nonPAR</c> random re-genotyping runs when
    /// <see cref="AsMultiPcfOptions.MaleXGenotyping"/> is set (B24 F58; otherwise the caller's genotypes are used as given).
    /// Missing values (NaN = R <c>NA</c>, B24 F59) follow R's NA path: only probes with a heterozygous BAF or a logR in at
    /// least one sample enter the segmentation (<c>Select_sites</c>), each missing value with weight 0; winsorisation and
    /// the logR averages skip NA; probes with both in at least one sample get a segmented BAF (<c>Select_sites2</c>); a part
    /// without such probes gets <c>mean(logR, na.rm = TRUE)</c>, and an all-missing level run takes the previous level.
    /// <para>R fails (<c>bafna[homo, ] &lt;- NA</c>: "incorrect number of subscripts on matrix") with a single sample or a
    /// single-probe chromosome part, because <c>Tumor_LogR[chr[[k]], ]</c> drops to a vector; both are rejected here.
    /// For one sample use <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/>.</para>
    /// </summary>
    /// <param name="samples">Per sample, the loci of the common probe set in genome order (same chromosome and position at
    /// each index in every sample). LogR must be finite or NaN (missing); BAF must be in [0, 1] or NaN at heterozygous
    /// probes and is ignored at homozygous probes. At least two samples; each contiguous same-chromosome run needs ≥ 2
    /// probes.</param>
    /// <param name="germlineHeterozygous">Germline genotype per probe (true = heterozygous, ASCAT
    /// <c>germlinegenotypes == FALSE</c>); one germline for all samples, as ascat.asmultipcf uses. Null = every probe
    /// heterozygous.</param>
    /// <param name="penalty">ascat.asmultipcf <c>penalty</c> (default 70 = <see cref="AspcfDefaultPenalty"/>), applied to the
    /// raw joint cost.</param>
    /// <param name="options">Sample weights, algorithm and refinement (null = <see cref="AsMultiPcfOptions.Default"/>).</param>
    /// <returns>One <see cref="AspcfSegmentation"/> per sample (input order): ASCAT <c>Tumor_LogR_segmented[, s]</c>,
    /// mirrored <c>Tumor_BAF_segmented[[s]]</c> (NaN at homozygous probes and where no BAF is segmented) and the runASCAT
    /// logR segments — ready for
    /// <see cref="FitPurityPloidyFromAspcf(AspcfSegmentation, double)"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="samples"/> or one of its entries is null.</exception>
    /// <exception cref="ArgumentException">Fewer than two samples, no probes, differing probe sets, a wrong genotype or
    /// weight count, a null chromosome label, an infinite logR, a heterozygous BAF outside [0, 1] (and not NaN), or a
    /// single-probe chromosome part.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="penalty"/> ≤ 0, NaN or infinite, or a weight that is
    /// not positive and finite.</exception>
    public static IReadOnlyList<AspcfSegmentation> SegmentAlleleSpecificAsMultiPcf(
        IReadOnlyList<IReadOnlyList<AlleleSpecificLocus>> samples,
        IReadOnlyList<bool>? germlineHeterozygous = null,
        double penalty = AspcfDefaultPenalty,
        AsMultiPcfOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(samples);
        options ??= AsMultiPcfOptions.Default;
        if (!double.IsFinite(penalty) || penalty <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(penalty), penalty, "The ASPCF penalty γ must be positive and finite.");
        }

        int sampleCount = samples.Count;
        if (sampleCount < 2)
        {
            throw new ArgumentException(
                "ascat.asmultipcf needs at least two samples (R fails on one); use SegmentAlleleSpecificAspcf for a single sample.",
                nameof(samples));
        }

        AlleleSpecificLocus[][] probes = new AlleleSpecificLocus[sampleCount][];
        for (int s = 0; s < sampleCount; s++)
        {
            if (samples[s] is null)
            {
                throw new ArgumentNullException(nameof(samples), "A sample is null.");
            }

            probes[s] = samples[s].ToArray();
        }

        int n = probes[0].Length;
        if (n == 0)
        {
            throw new ArgumentException("At least one probe is required.", nameof(samples));
        }

        var het = new bool[n];
        if (germlineHeterozygous is null)
        {
            Array.Fill(het, true);
        }
        else if (germlineHeterozygous.Count != n)
        {
            throw new ArgumentException("Exactly one germline genotype is required per probe.", nameof(germlineHeterozygous));
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                het[i] = germlineHeterozygous[i];
            }
        }

        if (options.MaleXGenotyping is not null)
        {
            foreach (AlleleSpecificLocus locus in probes[0])
            {
                if (locus.Chromosome is null)
                {
                    throw new ArgumentException("A locus has a null chromosome label.", nameof(samples));
                }
            }

            RegenotypeMaleXNonPar(probes[0], het, options.MaleXGenotyping, nameof(options));
        }

        ValidateAsMultiPcfProbes(probes, het);
        double[] trackWeights = AsMultiPcfTrackWeights(options.SampleWeights, sampleCount);

        List<(int Lo, int Hi)> runs = ContiguousChromosomeRuns(probes[0]);
        foreach ((int lo, int hi) in runs)
        {
            if (hi == lo)
            {
                throw new ArgumentException(
                    "Every chromosome part needs at least two probes (R's Tumor_LogR[chr[[k]], ] drops to a vector and ascat.asmultipcf fails).",
                    nameof(samples));
            }
        }

        var rawLogR = new double[sampleCount][];
        for (int s = 0; s < sampleCount; s++)
        {
            rawLogR[s] = new double[n];
            for (int i = 0; i < n; i++)
            {
                rawLogR[s][i] = probes[s][i].LogR;
            }
        }

        var ladder = new List<double> { penalty };
        foreach (double p in AsMultiPcfPenaltyLadder)
        {
            if (p > penalty)
            {
                ladder.Add(p);
            }
        }

        double[][] logRPcfed = Array.Empty<double[]>();
        double[][] bafPcfed = Array.Empty<double[]>();
        foreach (double segmentLength in ladder)
        {
            logRPcfed = new double[sampleCount][];
            bafPcfed = new double[sampleCount][];
            for (int s = 0; s < sampleCount; s++)
            {
                logRPcfed[s] = new double[n];
                bafPcfed[s] = new double[n];
                Array.Fill(bafPcfed[s], double.NaN);
            }

            foreach ((int lo, int hi) in runs)
            {
                SegmentPartAsMultiPcf(probes, rawLogR, het, lo, hi, segmentLength, trackWeights, options, logRPcfed, bafPcfed);
            }

            bool fewLevels = true;
            for (int s = 0; s < sampleCount; s++)
            {
                ReadaptLevels(logRPcfed[s], rawLogR[s], 0, n);
                if (logRPcfed[s].Distinct().Count() >= AspcfMaxSegmentLevels)
                {
                    fewLevels = false;
                }
            }

            if (fewLevels)
            {
                break;
            }
        }

        var result = new AspcfSegmentation[sampleCount];
        for (int s = 0; s < sampleCount; s++)
        {
            result[s] = new AspcfSegmentation(
                probes[s], (bool[])het.Clone(), logRPcfed[s], bafPcfed[s],
                BuildAspcfSegments(probes[s], het, runs, logRPcfed[s], bafPcfed[s]));
        }

        return result;
    }

    /// <summary>Same probe set in every sample; logR finite or NaN (R NA); heterozygous BAF in [0, 1] or NaN (B24 F59).</summary>
    private static void ValidateAsMultiPcfProbes(AlleleSpecificLocus[][] samples, bool[] het)
    {
        int n = samples[0].Length;
        foreach (AlleleSpecificLocus[] sample in samples)
        {
            if (sample.Length != n)
            {
                throw new ArgumentException("Every sample must have the same probes.", nameof(samples));
            }

            for (int i = 0; i < n; i++)
            {
                if (sample[i].Chromosome is null)
                {
                    throw new ArgumentException("A locus has a null chromosome label.", nameof(samples));
                }

                if (!string.Equals(sample[i].Chromosome, samples[0][i].Chromosome, StringComparison.Ordinal)
                    || sample[i].Position != samples[0][i].Position)
                {
                    throw new ArgumentException("Every sample must have the same probes (chromosome and position) in the same order.", nameof(samples));
                }

                if (double.IsInfinity(sample[i].LogR))
                {
                    throw new ArgumentException("A logR must be finite or NaN (missing).", nameof(samples));
                }

                if (het[i] && !IsAscatBafOrMissing(sample[i].BAF))
                {
                    throw new ArgumentException("Every heterozygous locus needs a BAF in [0, 1] or NaN (missing).", nameof(samples));
                }
            }
        }
    }

    /// <summary>ascat.asmultipcf <c>w * wsample</c>: per-track weights (S logR tracks, then S BAF tracks).</summary>
    private static double[] AsMultiPcfTrackWeights(IReadOnlyList<double>? sampleWeights, int sampleCount)
    {
        var weights = new double[2 * sampleCount];
        if (sampleWeights is null)
        {
            Array.Fill(weights, 1.0);
            return weights;
        }

        if (sampleWeights.Count == sampleCount)
        {
            for (int s = 0; s < sampleCount; s++)
            {
                weights[s] = sampleWeights[s];
                weights[s + sampleCount] = sampleWeights[s]; // wsample <- c(wsample, wsample)
            }
        }
        else if (sampleWeights.Count == 2 * sampleCount)
        {
            for (int t = 0; t < weights.Length; t++)
            {
                weights[t] = sampleWeights[t];
            }
        }
        else
        {
            throw new ArgumentException(
                $"Length of the sample weights is {sampleWeights.Count} but should be {sampleCount} or {2 * sampleCount}.",
                nameof(sampleWeights));
        }

        foreach (double w in weights)
        {
            if (!double.IsFinite(w) || w <= 0.0)
            {
                throw new ArgumentOutOfRangeException(nameof(sampleWeights), w, "Sample weights must be positive and finite.");
            }
        }

        return weights;
    }

    /// <summary>
    /// The per-chromosome-part body of <c>ascat.asmultipcf</c> on probes [lo, hi], including R's NA path (B24 F59; NaN =
    /// R <c>NA</c>): <c>Select_sites</c> = probes with a BAF (heterozygous, non-missing) or a logR in at least one sample;
    /// <c>Select_sites2</c> = probes with both in at least one sample (they get a segmented BAF). The joint segmentation runs
    /// on the selected sites with weight 0 at every missing value; missing logR is skipped in the winsorisation and the
    /// per-site logR averages (<c>na.rm = TRUE</c>). With complete data every probe is selected and the heterozygous ones
    /// get a BAF.
    /// </summary>
    private static void SegmentPartAsMultiPcf(
        AlleleSpecificLocus[][] probes, double[][] rawLogR, bool[] het, int lo, int hi, double segmentLength,
        double[] trackWeights, AsMultiPcfOptions options, double[][] logRPcfed, double[][] bafPcfed)
    {
        int sampleCount = probes.Length;
        int tracks = 2 * sampleCount;
        int len = hi - lo + 1;

        // bafna (NA at homozygous probes), Select_sites ("indices", 0-based local) and Select_sites2.
        var indices = new List<int>();
        var outputBaf = new List<bool>(); // Select_sites2 per selected site
        for (int i = 0; i < len; i++)
        {
            bool anyBaf = false, anyLogR = false;
            for (int s = 0; s < sampleCount; s++)
            {
                anyBaf |= het[lo + i] && !double.IsNaN(probes[s][lo + i].BAF);
                anyLogR |= !double.IsNaN(rawLogR[s][lo + i]);
            }

            if (anyBaf || anyLogR)
            {
                indices.Add(i);
                outputBaf.Add(anyBaf && anyLogR);
            }
        }

        int m = indices.Count;
        if (m == 0)
        {
            // "if there are no probes in the segment, don't do anything, except add a LogR segment" (mean, na.rm = TRUE).
            for (int s = 0; s < sampleCount; s++)
            {
                Array.Fill(logRPcfed[s], MeanNaRm(rawLogR[s], lo, hi + 1), lo, len);
            }

            return;
        }

        // bafwins = mirror(madWinsMatrixWithNA(mirror(bafsel))): per sample over its non-missing values; with a single
        // selected site R's madWinsMatrixWithNA winsorises that row ACROSS the samples (nrow(x) == 1 branch).
        var bafWins = new double[sampleCount][];
        for (int s = 0; s < sampleCount; s++)
        {
            bafWins[s] = new double[m];
            for (int k = 0; k < m; k++)
            {
                int i = lo + indices[k];
                bafWins[s][k] = het[i] ? MirrorBaf(probes[s][i].BAF) : double.NaN; // MirrorBaf(NaN) = NaN
            }
        }

        if (m == 1)
        {
            var row = new double[sampleCount];
            for (int s = 0; s < sampleCount; s++)
            {
                row[s] = bafWins[s][0];
            }

            double[] rowWins = MadWinsorizeWithNa(row, AspcfWinsorTau, AspcfMedianHalfWindow);
            for (int s = 0; s < sampleCount; s++)
            {
                bafWins[s][0] = MirrorBaf(rowWins[s]);
            }
        }
        else
        {
            for (int s = 0; s < sampleCount; s++)
            {
                double[] wins = MadWinsorizeWithNa(bafWins[s], AspcfWinsorTau, AspcfMedianHalfWindow);
                for (int k = 0; k < m; k++)
                {
                    bafWins[s][k] = MirrorBaf(wins[k]);
                }
            }
        }

        // lrwins (all probes of the part) and logRaveraged: mean(lrwins[start:end], na.rm = TRUE) around each selected site.
        var logRAveraged = new double[sampleCount][];
        for (int s = 0; s < sampleCount; s++)
        {
            var lr = new double[len];
            Array.Copy(rawLogR[s], lo, lr, 0, len);
            double[] lrWins = MadWinsorizeWithNa(lr, AspcfWinsorTau, AspcfMedianHalfWindow);
            logRAveraged[s] = new double[m];
            for (int k = 0; k < m; k++)
            {
                (int start, int end) = AscatAveragingWindow(indices, k, len);
                logRAveraged[s][k] = RColonMeanNaRm(lrWins, start, end);
            }
        }

        // logRASPCF / bafASPCF: per sample, per selected site.
        var logRAspcf = new double[sampleCount][];
        var bafAspcf = new double[sampleCount][];
        if (m < AsMultiPcfMinProbes)
        {
            // (nrow == 1 keeps the site's own values, which the means below reproduce.)
            for (int s = 0; s < sampleCount; s++)
            {
                double bafMean = MeanNaRm(bafWins[s], 0, m);
                logRAspcf[s] = Filled(m, MeanNaRm(logRAveraged[s], 0, m));
                bafAspcf[s] = Filled(m, MirrorBaf(bafMean)); // R bafMean >= 0.5 ? … : 1 − …: equal at 0.5; NaN stays NaN
            }
        }
        else
        {
            // lrANDbaf with NA → 0 and the weight matrix w (0 at NA, times wsample), tracks in rows.
            var value = new double[tracks, m];
            var weight = new double[tracks, m];
            for (int s = 0; s < sampleCount; s++)
            {
                for (int k = 0; k < m; k++)
                {
                    double r = logRAveraged[s][k];
                    value[s, k] = double.IsNaN(r) ? 0.0 : r;
                    weight[s, k] = double.IsNaN(r) ? 0.0 : trackWeights[s];
                    double b = bafWins[s][k];
                    value[s + sampleCount, k] = double.IsNaN(b) ? 0.0 : b;
                    weight[s + sampleCount, k] = double.IsNaN(b) ? 0.0 : trackWeights[s + sampleCount];
                }
            }

            AsMultiPcfFit joint = options.Algorithm == AsMultiPcfAlgorithm.Exact
                ? AsMultiPcfCompact(weight, WeightedValues(value, weight), segmentLength, Ones(m))
                : RunFastAsMultiPcf(value, weight, segmentLength);
            double[,] yhat = ExpandMulti(m, tracks, joint.Lengths, joint.Means);
            for (int s = 0; s < sampleCount; s++)
            {
                logRAspcf[s] = new double[m];
                bafAspcf[s] = new double[m];
                for (int k = 0; k < m; k++)
                {
                    logRAspcf[s][k] = yhat[s, k];
                    bafAspcf[s][k] = yhat[s + sampleCount, k];
                }
            }

            if (options.Refine && joint.Starts.Length > 1)
            {
                RefineAsMultiPcf(value, weight, joint, segmentLength, logRAspcf, bafAspcf);
            }
        }

        // BAF output rows (Select_sites2) and the getMadwithNA input bafwins[het selected sites, s].
        var outRows = new List<int>();
        var hetRows = new List<int>();
        for (int k = 0; k < m; k++)
        {
            if (outputBaf[k])
            {
                outRows.Add(k);
            }

            if (het[lo + indices[k]])
            {
                hetRows.Add(k);
            }
        }

        var logRC = new double[len];
        for (int s = 0; s < sampleCount; s++)
        {
            var sdInput = new double[hetRows.Count];
            for (int q = 0; q < sdInput.Length; q++)
            {
                sdInput[q] = bafWins[s][hetRows[q]];
            }

            CorrectAsMultiPcfBaf(outRows, indices, bafAspcf[s], GetMadWithNa(sdInput, AspcfMedianHalfWindow), lo, bafPcfed[s]);

            // Expand to every probe: site k covers (indices[k−1], indices[k]]; the first from the part start, the last to
            // its end.
            int from = 0;
            for (int k = 0; k < m; k++)
            {
                int to = k == m - 1 ? len - 1 : indices[k];
                for (int i = from; i <= to; i++)
                {
                    logRC[i] = logRAspcf[s][k];
                }

                from = to + 1;
            }

            FillZeroLevelsFromClosest(logRC);
            Array.Copy(logRC, 0, logRPcfed[s], lo, len);
        }
    }

    /// <summary>
    /// ascat.aspcf / ascat.asmultipcf logR averaging window of selected site <paramref name="k"/> (1-based [start, end]
    /// into the part of <paramref name="len"/> probes): <c>averageIndices = c(1, (indices[-n] + indices[-1])/2, len + 0.01)</c>,
    /// start = ceiling, end = floor(· − 0.01); one site ⇒ the whole part. <paramref name="indices"/> are 0-based local.
    /// </summary>
    private static (int Start, int End) AscatAveragingWindow(List<int> indices, int k, int len)
    {
        int h = indices.Count;
        if (h == 1)
        {
            return (1, len);
        }

        double lower = k == 0 ? 1.0 : ((indices[k - 1] + 1) + (indices[k] + 1)) / 2.0;
        double upper = k == h - 1 ? len + 0.01 : ((indices[k] + 1) + (indices[k + 1] + 1)) / 2.0;
        return ((int)Math.Ceiling(lower), (int)Math.Floor(upper - 0.01));
    }

    /// <summary>ASCAT <c>getMadwithNA(x, k)</c>: <see cref="GetMad"/> after dropping missing (NaN) values.</summary>
    private static double GetMadWithNa(double[] x, int k) =>
        GetMad(Array.FindAll(x, v => !double.IsNaN(v)), k);

    /// <summary>
    /// ascat.asmultipcf refinement: per sample, its logR and BAF tracks are compacted on the joint breakpoints
    /// (<c>compactASMulti</c>) and re-segmented by <c>ASmultiPCFcompact</c> with γ/S; fewer segments ⇒ the sample takes them.
    /// </summary>
    private static void RefineAsMultiPcf(
        double[,] value, double[,] weight, AsMultiPcfFit joint, double segmentLength, double[][] logRAspcf, double[][] bafAspcf)
    {
        int sampleCount = logRAspcf.Length;
        int len = value.GetLength(1);
        var mark = new bool[len];
        for (int k = 1; k < joint.Starts.Length; k++)
        {
            mark[joint.Starts[k] - 2] = true; // start0[-1] - 1 (1-based) = last probe of the previous segment
        }

        mark[len - 1] = true;
        for (int s = 0; s < sampleCount; s++)
        {
            int[] rows = { s, s + sampleCount };
            var y = new double[2, len];
            var w = new double[2, len];
            for (int r = 0; r < 2; r++)
            {
                for (int i = 0; i < len; i++)
                {
                    y[r, i] = value[rows[r], i];
                    w[r, i] = weight[rows[r], i];
                }
            }

            (double[,] wSum, double[,] compWeight, int[] probesPerSeg) = CompactAsMulti(y, mark, w);
            AsMultiPcfFit refined = AsMultiPcfCompact(compWeight, wSum, segmentLength / sampleCount, probesPerSeg);
            if (refined.Starts.Length < joint.Starts.Length)
            {
                double[,] potts = ExpandMulti(len, 2, refined.Lengths, refined.Means);
                for (int i = 0; i < len; i++)
                {
                    logRAspcf[s][i] = potts[0, i];
                    bafAspcf[s][i] = potts[1, i];
                }
            }
        }
    }

    /// <summary>
    /// ascat.asmultipcf BAF correction of one sample on the BAF output sites of a part (<c>Select_sites2</c>,
    /// <paramref name="outRows"/> = positions among the selected sites): each <c>rle</c> run of equal segmented BAF gets
    /// 0.5 + μ, μ = |b − 0.5| (the mean of a constant run), shrunk to 0 when <c>sqrt(sd² + μ²) &lt; 2·sd</c>; a missing run
    /// stays missing. <paramref name="sd"/> = <c>getMadwithNA</c> of the sample's winsorised mirrored heterozygous BAF.
    /// </summary>
    private static void CorrectAsMultiPcfBaf(
        List<int> outRows, List<int> indices, double[] bafAspcf, double sd, int lo, double[] bafPcfed)
    {
        int h = outRows.Count;
        int runStart = 0;
        for (int k = 1; k <= h; k++)
        {
            if (k < h && bafAspcf[outRows[k]].Equals(bafAspcf[outRows[runStart]]))
            {
                continue;
            }

            double yi = bafAspcf[outRows[runStart]];
            double level = double.NaN;
            if (!double.IsNaN(yi))
            {
                level = AscatShrunkBafLevel(Math.Abs(yi - BalancedBaf), sd);
            }

            for (int q = runStart; q < k; q++)
            {
                bafPcfed[lo + indices[outRows[q]]] = level;
            }

            runStart = k;
        }
    }

    /// <summary>
    /// ascat.asmultipcf "fill in NAs in logR data": levels that are NaN or exactly 0 take the level of the closest other
    /// probe of the part (<c>which.min</c>: the lower index on a tie); all 0 when none remains.
    /// </summary>
    private static void FillZeroLevelsFromClosest(double[] levels)
    {
        var keep = new List<int>();
        for (int i = 0; i < levels.Length; i++)
        {
            if (!double.IsNaN(levels[i]) && levels[i] != 0.0)
            {
                keep.Add(i);
            }
        }

        if (keep.Count == levels.Length)
        {
            return;
        }

        var original = (double[])levels.Clone();
        for (int i = 0; i < levels.Length; i++)
        {
            if (!double.IsNaN(original[i]) && original[i] != 0.0)
            {
                continue;
            }

            if (keep.Count == 0)
            {
                levels[i] = 0.0;
                continue;
            }

            int best = keep[0];
            foreach (int k in keep)
            {
                if (Math.Abs(k - i) < Math.Abs(best - i))
                {
                    best = k;
                }
            }

            levels[i] = original[best];
        }
    }

    /// <summary>Result of <c>ASmultiPCFcompact</c>: segment lengths in probes (<c>Lengde</c>), 1-based starts (<c>sta</c>),
    /// per-track segment means [track, segment] and the number of segments.</summary>
    private sealed record AsMultiPcfFit(int[] Lengths, int[] Starts, double[,] Means);

    /// <summary>
    /// ASCAT <c>ASmultiPCFcompact(nr, wSum, gamma, wsample, nProbesPerSeg)</c>, verbatim: exact multi-track Potts DP on
    /// (compacted) columns with weights <paramref name="nr"/> and weighted sums <paramref name="wSum"/> [track, column].
    /// Cost[j] = Σ_tracks −Sum²/Nevner (Nevner 0 ⇒ 1) + bestCost[j − 1] + γ; <c>which.min</c> skips NaN (bestCost[1] is NaN
    /// when the first column has a zero weight, exactly as in R). <c>wsample</c> is unused by R and not taken.
    /// </summary>
    private static AsMultiPcfFit AsMultiPcfCompact(double[,] nr, double[,] wSum, double gamma, int[] nProbesPerSeg)
    {
        int tracks = wSum.GetLength(0);
        int bigN = nr.GetLength(1);
        var bestCost = new double[bigN + 1];   // 1-based
        var bestSplit = new int[bigN + 2];     // 1-based
        var bestAver = new double[tracks, bigN + 1];
        var sum = new double[tracks, bigN + 1];
        var nevner = new double[tracks, bigN + 1];
        var cost = new double[bigN + 1];

        double first = 0.0;
        for (int t = 0; t < tracks; t++)
        {
            sum[t, 1] = wSum[t, 0];
            nevner[t, 1] = nr[t, 0];
            bestAver[t, 1] = wSum[t, 0] / nr[t, 0];
            first += -sum[t, 1] * bestAver[t, 1]; // helper %*% (-Sum[, 1] * bestAver[, 1])
        }

        bestCost[1] = first;
        for (int m = 2; m <= bigN; m++)
        {
            for (int j = 1; j <= m; j++)
            {
                double c = 0.0;
                for (int t = 0; t < tracks; t++)
                {
                    sum[t, j] += wSum[t, m - 1];
                    nevner[t, j] += nr[t, m - 1];
                    double denominator = nevner[t, j] == 0.0 ? 1.0 : nevner[t, j];
                    c += -(sum[t, j] * sum[t, j]) / denominator;
                }

                cost[j] = j >= 2 ? c + bestCost[j - 1] + gamma : c;
            }

            int pos = -1;
            for (int j = 1; j <= m; j++)
            {
                if (!double.IsNaN(cost[j]) && (pos < 0 || cost[j] < cost[pos]))
                {
                    pos = j; // which.min: first minimum, NA/NaN skipped
                }
            }

            bestCost[m] = cost[pos];
            for (int t = 0; t < tracks; t++)
            {
                bestAver[t, m] = sum[t, pos] / nevner[t, pos];
            }

            bestSplit[m] = pos - 1;
        }

        var lengths = new List<int>();
        var means = new List<double[]>();
        int cursor = bigN;
        while (cursor > 0)
        {
            int probesInSegment = 0;
            for (int k = bestSplit[cursor] + 1; k <= cursor; k++)
            {
                probesInSegment += nProbesPerSeg[k - 1];
            }

            lengths.Add(probesInSegment);
            var mean = new double[tracks];
            for (int t = 0; t < tracks; t++)
            {
                mean[t] = bestAver[t, cursor];
            }

            means.Add(mean);
            cursor = bestSplit[cursor];
        }

        lengths.Reverse();
        means.Reverse();
        var starts = new int[lengths.Count];
        starts[0] = 1;
        for (int k = 1; k < starts.Length; k++)
        {
            starts[k] = starts[k - 1] + lengths[k - 1];
        }

        var meanMatrix = new double[tracks, lengths.Count];
        for (int k = 0; k < lengths.Count; k++)
        {
            for (int t = 0; t < tracks; t++)
            {
                meanMatrix[t, k] = means[k][t];
            }
        }

        return new AsMultiPcfFit(lengths.ToArray(), starts, meanMatrix);
    }

    /// <summary>
    /// ASCAT <c>compactASMulti(y, mark, w)</c>: weighted sums, weights and probe counts of the blocks ending at each marked
    /// column (the last column must be marked).
    /// </summary>
    private static (double[,] WSum, double[,] Weight, int[] ProbesPerSeg) CompactAsMulti(double[,] y, bool[] mark, double[,] w)
    {
        int tracks = y.GetLength(0);
        int len = y.GetLength(1);
        int blocks = mark.Count(m => m);
        var wSum = new double[tracks, blocks];
        var weight = new double[tracks, blocks];
        var probesPerSeg = new int[blocks];
        var delSum = new double[tracks];
        var delWeight = new double[tracks];
        int count = 0, oldPos = -1;
        for (int pos = 0; pos < len; pos++)
        {
            for (int t = 0; t < tracks; t++)
            {
                double contribution = y[t, pos] * w[t, pos];
                if (mark[pos])
                {
                    wSum[t, count] = delSum[t] + contribution;
                    weight[t, count] = delWeight[t] + w[t, pos];
                    delSum[t] = 0.0;
                    delWeight[t] = 0.0;
                }
                else
                {
                    delSum[t] += contribution;
                    delWeight[t] += w[t, pos];
                }
            }

            if (mark[pos])
            {
                probesPerSeg[count] = pos - oldPos;
                oldPos = pos;
                count++;
            }
        }

        return (wSum, weight, probesPerSeg);
    }

    /// <summary>
    /// ASCAT <c>runFastASMultiPCF(x, w, gamma, yest = TRUE, subsize = 5000)</c>, verbatim: <c>ASmultiPCFcompact</c> on
    /// windows of subsize + 1 probes (step 4·subsize/5) with the <b>unweighted</b> values as <c>wSum</c> (as R passes
    /// <c>t(x[…])</c>), marking each segment start; all probes from the last window start on are marked; the exact
    /// algorithm then runs on the compacted data.
    /// </summary>
    private static AsMultiPcfFit RunFastAsMultiPcf(double[,] value, double[,] weight, double gamma)
    {
        int tracks = value.GetLength(0);
        int antGen = value.GetLength(1);
        var mark = new bool[antGen];
        int start0 = 1;
        while (start0 + AsMultiPcfFastSubsize < antGen)
        {
            int width = AsMultiPcfFastSubsize + 1;
            var nr = new double[tracks, width];
            var x = new double[tracks, width];
            for (int t = 0; t < tracks; t++)
            {
                for (int i = 0; i < width; i++)
                {
                    nr[t, i] = weight[t, start0 - 1 + i];
                    x[t, i] = value[t, start0 - 1 + i];
                }
            }

            AsMultiPcfFit res = AsMultiPcfCompact(nr, x, gamma, Ones(width));
            foreach (int sta in res.Starts)
            {
                mark[start0 + sta - 2] = true; // mark[start0 + res$sta - 1] (1-based)
            }

            start0 += AsMultiPcfFastStep;
        }

        for (int i = start0 - 1; i < antGen; i++)
        {
            mark[i] = true;
        }

        (double[,] wSum, double[,] compWeight, int[] probesPerSeg) = CompactAsMulti(value, mark, weight);
        return AsMultiPcfCompact(compWeight, wSum, gamma, probesPerSeg);
    }

    /// <summary>ASCAT <c>expandMulti(nProbes, nSamples, lengthInt, mean)</c>: [track, probe] piecewise-constant expansion.</summary>
    private static double[,] ExpandMulti(int nProbes, int tracks, int[] lengths, double[,] means)
    {
        var potts = new double[tracks, nProbes];
        int k = 0;
        for (int seg = 0; seg < lengths.Length; seg++)
        {
            for (int j = 0; j < lengths[seg]; j++)
            {
                for (int t = 0; t < tracks; t++)
                {
                    potts[t, k] = means[t, seg];
                }

                k++;
            }
        }

        return potts;
    }

    /// <summary><c>t(lrANDbaf) * w</c>.</summary>
    private static double[,] WeightedValues(double[,] value, double[,] weight)
    {
        int tracks = value.GetLength(0), len = value.GetLength(1);
        var result = new double[tracks, len];
        for (int t = 0; t < tracks; t++)
        {
            for (int i = 0; i < len; i++)
            {
                result[t, i] = value[t, i] * weight[t, i];
            }
        }

        return result;
    }

    /// <summary><c>rep(1, n)</c> as probe counts.</summary>
    private static int[] Ones(int n)
    {
        var ones = new int[n];
        Array.Fill(ones, 1);
        return ones;
    }
}
