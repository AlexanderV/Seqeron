// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Battenberg phased-BAF segmentation (port of Battenberg 3.0.0).</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{

    /// <summary>
    /// One germline-heterozygous SNP with its tumour BAF, expressed on a <b>caller-supplied</b> haplotype phase (Battenberg
    /// <c>combine.baf.files</c> output, columns <c>Chromosome, Position, BAF</c>; the BAF is that of haplotype 1 after
    /// IMPUTE2/Beagle5 phasing). <see cref="double.NaN"/> marks a missing BAF (dropped, as Battenberg drops <c>NA</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">SNP position (the SNPs of a chromosome are taken in input order, as in Battenberg).</param>
    /// <param name="Baf">Haplotype-phased BAF in [0, 1], or NaN.</param>
    public readonly record struct PhasedBafSnp(string Chromosome, long Position, double Baf);

    /// <summary>A prior breakpoint (e.g. a structural-variant junction): Battenberg <c>prior_breakpoints_file</c> row.</summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">Breakpoint position.</param>
    public readonly record struct BattenbergPriorBreakpoint(string Chromosome, long Position);

    /// <summary>Battenberg <c>calc_seg_baf_option</c>: how the BAF of a segment is recalculated after segmentation.</summary>
    public enum BattenbergSegmentBafOption
    {
        /// <summary>1 — median of the segment's phased BAFs (<c>adjustSegmValues</c>).</summary>
        Median = 1,

        /// <summary>2 — the PCF segment mean, unchanged.</summary>
        Mean = 2,

        /// <summary>3 (Battenberg default) — the median, unless it is exactly 0 or 1, then the mean.</summary>
        MedianUnlessExtreme = 3,
    }

    /// <summary>
    /// Parameters of Battenberg <c>segment.baf.phased</c> (Wedge-lab/battenberg R/segmentation.R); defaults are
    /// Battenberg's (<c>gamma = 10, phasegamma = 3, kmin = 3, phasekmin = 3, no_segmentation = F,
    /// calc_seg_baf_option = 3</c>, no prior breakpoints).
    /// </summary>
    public sealed record BattenbergPhasedSegmentationOptions
    {
        /// <summary>Battenberg defaults.</summary>
        public static BattenbergPhasedSegmentationOptions Default { get; } = new();

        /// <summary><c>gamma</c>: PCF penalty factor of the second (copy-number) segmentation (multiplied by the MAD sd).</summary>
        public double Gamma { get; init; } = 10.0;

        /// <summary><c>phasegamma</c>: PCF penalty factor of the first (phase-correcting) segmentation.</summary>
        public double PhaseGamma { get; init; } = 3.0;

        /// <summary><c>kmin</c>: minimum SNPs per segment of the second segmentation, in [1, <see cref="BattenbergMaxKmin"/>].</summary>
        public int Kmin { get; init; } = 3;

        /// <summary><c>phasekmin</c>: minimum SNPs per segment of the first segmentation, in [1, <see cref="BattenbergMaxKmin"/>].</summary>
        public int PhaseKmin { get; init; } = 3;

        /// <summary><c>no_segmentation</c>: switch the haplotype blocks but take the mean phased BAF as the segment BAF.</summary>
        public bool NoSegmentation { get; init; }

        /// <summary><c>calc_seg_baf_option</c> (default <see cref="BattenbergSegmentBafOption.MedianUnlessExtreme"/>).</summary>
        public BattenbergSegmentBafOption SegmentBafOption { get; init; } = BattenbergSegmentBafOption.MedianUnlessExtreme;

        /// <summary>Prior breakpoints (<c>prior_breakpoints_file</c>), in file order; empty = none.</summary>
        public IReadOnlyList<BattenbergPriorBreakpoint> PriorBreakpoints { get; init; } = Array.Empty<BattenbergPriorBreakpoint>();
    }

    /// <summary>
    /// Largest <c>kmin</c>/<c>phasekmin</c> accepted: Battenberg's <c>filterMarkS4</c> indexes <c>3·kmin + 6</c> SNPs,
    /// so larger values leave R's index ranges for the minimum PCF input of 50 SNPs (R then misbehaves or stops).
    /// </summary>
    public const int BattenbergMaxKmin = 14;

    /// <summary>
    /// One output row of Battenberg <c>segment.baf.phased</c> (<c>Chromosome, Position, BAF, BAFphased, BAFseg</c>).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">SNP position.</param>
    /// <param name="Baf">The input phased BAF.</param>
    /// <param name="BafPhased">Haplotype-block-corrected BAF <c>ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>.</param>
    /// <param name="BafSegment">The segment BAF <c>BAFseg</c> (per <see cref="BattenbergSegmentBafOption"/>).</param>
    public readonly record struct PhasedBafSegmentedSnp(
        string Chromosome, long Position, double Baf, double BafPhased, double BafSegment);

    /// <summary>A raw logR probe (Battenberg <c>LogRvals</c>: chromosome, position, logR). NaN/±∞ logR are ignored.</summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Position">Probe position.</param>
    /// <param name="LogR">Raw logR.</param>
    public readonly record struct LogRProbe(string Chromosome, long Position, double LogR);

    /// <summary>Battenberg <c>bkps_to_presegment_breakpoints</c>: <c>maxsnpdist = 3000000</c> (a gap ≥ this splits).</summary>
    private const long BattenbergMaxSnpDistance = 3_000_000;

    /// <summary>Battenberg <c>run_pcf</c>: presegments with fewer SNPs are not segmented (their mean BAF is used).</summary>
    private const int BattenbergMinPcfSnps = 50;

    /// <summary>Battenberg <c>run_pcf</c>: MAD sd floor 0.09 ("binomial around 0.5 at depth 30").</summary>
    private const double BattenbergMinBafSd = 0.09;

    /// <summary>
    /// Segments caller-phased SNP BAFs exactly as Battenberg <c>segment.baf.phased</c> (Wedge-lab/battenberg
    /// R/segmentation.R, default path; B24 F40), the step between haplotype phasing and <c>callSubclones</c>:
    /// <list type="number">
    /// <item>Per chromosome (first-appearance order), SNPs with a missing BAF are dropped and the chromosome is
    /// pre-segmented at prior breakpoints and at every gap ≥ 3 Mb between consecutive SNPs
    /// (<c>bkps_to_presegment_breakpoints</c>/<c>addin_bigholes</c>, including R's index quirks).</item>
    /// <item>Per presegment (<c>run_pcf</c>): sd = <c>getMad(min(BAF, 1 − BAF), k = 25)</c> (NA → 0, floored at 0.09);
    /// a first PCF <c>selectFastPcf(BAF, phasekmin, phasegamma·sd)</c> finds the switched haplotype blocks,
    /// <c>BAFphased = ifelse(BAFsegm &gt; 0.5, BAF, 1 − BAF)</c>; a second PCF <c>selectFastPcf(BAFphased, kmin,
    /// gamma·sd)</c> gives <c>BAFseg</c> (presegments of &lt; 50 SNPs, or <c>no_segmentation</c>, use the mean); the
    /// segment BAF is then recalculated per <c>calc_seg_baf_option</c>.</item>
    /// </list>
    /// <c>selectFastPcf</c> is the copynumber-package fast PCF that Battenberg ships in R/fastPCF.R (<c>filterMarkS4</c>
    /// candidate breakpoints + <c>PottsCompact</c>, and the 5000-SNP windowed <c>runPcfSubset</c> above 15 000 SNPs) —
    /// a different routine from ASCAT's <c>exactPcf</c>/<c>fastAspcf</c>, ported here once; <c>getMad</c>/<c>runmed</c>
    /// are the shared ASCAT helpers (identical code). Haplotype phasing itself (IMPUTE2/Beagle5 against the 1000 Genomes
    /// reference panel) is <b>not</b> run: the caller supplies phased BAFs. Agreement with Battenberg (R 4.3.3) is to
    /// ≲ 1e−12; R accumulates <c>cumsum</c>/<c>mean</c> in 80-bit long double, which can only matter at exact cost ties.
    /// </summary>
    /// <param name="snps">Phased SNP BAFs. Non-null; non-null chromosome; BAF in [0, 1] or NaN.</param>
    /// <param name="options">Battenberg parameters (default <see cref="BattenbergPhasedSegmentationOptions.Default"/>).</param>
    /// <returns>Battenberg's per-SNP rows (chromosomes in first-appearance order, presegments in order).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="snps"/>, a chromosome, or the breakpoint list is null.</exception>
    /// <exception cref="ArgumentException">a BAF outside [0, 1] (and not NaN).</exception>
    /// <exception cref="ArgumentOutOfRangeException">gamma/phasegamma not finite or &lt; 0; kmin/phasekmin outside [1, 14];
    /// an undefined <see cref="BattenbergSegmentBafOption"/>.</exception>
    public static IReadOnlyList<PhasedBafSegmentedSnp> SegmentPhasedBaf(
        IReadOnlyList<PhasedBafSnp> snps,
        BattenbergPhasedSegmentationOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(snps);
        options ??= BattenbergPhasedSegmentationOptions.Default;
        ValidateBattenbergSegmentationOptions(options);

        var chromosomes = new List<string>();
        var byChromosome = new Dictionary<string, (List<long> Pos, List<double> Baf)>(StringComparer.Ordinal);
        foreach (PhasedBafSnp snp in snps)
        {
            if (snp.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(snps), "Every SNP needs a chromosome.");
            }

            if (!(double.IsNaN(snp.Baf) || (snp.Baf >= 0.0 && snp.Baf <= 1.0)))
            {
                throw new ArgumentException("Every phased BAF must lie in [0, 1] (or be NaN = missing).", nameof(snps));
            }

            if (!byChromosome.TryGetValue(snp.Chromosome, out var lists))
            {
                lists = (new List<long>(), new List<double>());
                byChromosome.Add(snp.Chromosome, lists);
                chromosomes.Add(snp.Chromosome);
            }

            if (!double.IsNaN(snp.Baf))
            {
                lists.Pos.Add(snp.Position);
                lists.Baf.Add(snp.Baf);
            }
        }

        var output = new List<PhasedBafSegmentedSnp>(snps.Count);
        foreach (string chr in chromosomes)
        {
            (List<long> posList, List<double> bafList) = byChromosome[chr];
            long[] pos = posList.ToArray();
            double[] baf = bafList.ToArray();
            var bkps = new List<long>();
            foreach (BattenbergPriorBreakpoint b in options.PriorBreakpoints)
            {
                if (b.Chromosome == chr)
                {
                    bkps.Add(b.Position);
                }
            }

            foreach ((long? start, long? end) in BattenbergPresegments(pos, bkps))
            {
                BattenbergRunPcf(chr, pos, baf, start, end, options, output);
            }
        }

        return output;
    }

    /// <summary>
    /// Groups <see cref="SegmentPhasedBaf"/> rows into the segments Battenberg <c>determine_copynumber</c> fits — a new
    /// segment starts wherever <c>BAFseg</c> or the chromosome changes (<c>switchpoints</c>) — with each segment's logR
    /// = mean of the finite logR probes of the same chromosome within [first SNP position, last SNP position]
    /// (<c>mean(LogRvals[LogRpos &gt;= startpos &amp; LogRpos &lt;= endpos &amp; !is.infinite(.), 3], na.rm = T)</c>, 0 when
    /// none). The result feeds <see cref="FitSubclonalCopyNumberWithSnpTest"/>: segment BAF = <c>BAFseg</c> (the level
    /// <c>l</c> is mirrored there), SNP BAFs = <c>BAFphased</c> (<c>BAFke</c>), Start/End = min/max SNP position.
    /// </summary>
    /// <param name="segmentedSnps">Rows of <see cref="SegmentPhasedBaf"/> (BAFs in [0, 1]). Non-null.</param>
    /// <param name="logR">Raw logR probes. Non-null (may be empty; NaN/±∞ ignored).</param>
    /// <returns>One entry per Battenberg segment, in row order.</returns>
    /// <exception cref="ArgumentNullException">an argument or a chromosome is null.</exception>
    /// <exception cref="ArgumentException">a BAF or phased BAF outside [0, 1].</exception>
    public static IReadOnlyList<SubclonalSegmentSnpBafs> BuildBattenbergSegments(
        IReadOnlyList<PhasedBafSegmentedSnp> segmentedSnps,
        IReadOnlyList<LogRProbe> logR)
    {
        ArgumentNullException.ThrowIfNull(segmentedSnps);
        ArgumentNullException.ThrowIfNull(logR);
        var probes = new Dictionary<string, List<(long Pos, double LogR)>>(StringComparer.Ordinal);
        foreach (LogRProbe p in logR)
        {
            if (p.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(logR), "Every logR probe needs a chromosome.");
            }

            if (!double.IsFinite(p.LogR))
            {
                continue; // !is.infinite(.) and na.rm = T
            }

            if (!probes.TryGetValue(p.Chromosome, out var list))
            {
                list = new List<(long, double)>();
                probes.Add(p.Chromosome, list);
            }

            list.Add((p.Position, p.LogR));
        }

        var result = new List<SubclonalSegmentSnpBafs>();
        int n = segmentedSnps.Count;
        int from = 0;
        for (int i = 0; i < n; i++)
        {
            PhasedBafSegmentedSnp row = segmentedSnps[i];
            if (row.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(segmentedSnps), "Every row needs a chromosome.");
            }

            if (!(row.BafSegment >= 0.0 && row.BafSegment <= 1.0) || !(row.BafPhased >= 0.0 && row.BafPhased <= 1.0))
            {
                throw new ArgumentException("Every BAFseg and BAFphased must lie in [0, 1].", nameof(segmentedSnps));
            }

            bool last = i == n - 1
                        || segmentedSnps[i + 1].BafSegment != row.BafSegment
                        || segmentedSnps[i + 1].Chromosome != row.Chromosome;
            if (!last)
            {
                continue;
            }

            long start = long.MaxValue, end = long.MinValue;
            var phased = new double[i - from + 1];
            for (int j = from; j <= i; j++)
            {
                start = Math.Min(start, segmentedSnps[j].Position);
                end = Math.Max(end, segmentedSnps[j].Position);
                phased[j - from] = segmentedSnps[j].BafPhased;
            }

            var inRange = new List<double>();
            if (probes.TryGetValue(row.Chromosome, out var chrProbes))
            {
                foreach ((long pos, double value) in chrProbes)
                {
                    if (pos >= start && pos <= end)
                    {
                        inRange.Add(value);
                    }
                }
            }

            double segLogR = inRange.Count == 0 ? 0.0 : StatisticsHelper.RMean(inRange); // is.na(LogR) ⇒ 0
            result.Add(new SubclonalSegmentSnpBafs(
                new AlleleSpecificSegmentSummary(row.Chromosome, start, end, segLogR, row.BafSegment, phased.Length),
                phased));
            from = i + 1;
        }

        return result;
    }

    private static void ValidateBattenbergSegmentationOptions(BattenbergPhasedSegmentationOptions options)
    {
        BattenbergPhasedSegmentationOptions o = options;
        if (!(o.Gamma >= 0.0) || double.IsPositiveInfinity(o.Gamma))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.Gamma, "gamma must be finite and ≥ 0.");
        }

        if (!(o.PhaseGamma >= 0.0) || double.IsPositiveInfinity(o.PhaseGamma))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.PhaseGamma, "phasegamma must be finite and ≥ 0.");
        }

        if (o.Kmin < 1 || o.Kmin > BattenbergMaxKmin)
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.Kmin, "kmin must lie in [1, 14].");
        }

        if (o.PhaseKmin < 1 || o.PhaseKmin > BattenbergMaxKmin)
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.PhaseKmin, "phasekmin must lie in [1, 14].");
        }

        if (!Enum.IsDefined(o.SegmentBafOption))
        {
            throw new ArgumentOutOfRangeException(nameof(options), o.SegmentBafOption, "Unknown calc_seg_baf_option.");
        }

        ArgumentNullException.ThrowIfNull(o.PriorBreakpoints);
        foreach (BattenbergPriorBreakpoint b in o.PriorBreakpoints)
        {
            if (b.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(options), "Every prior breakpoint needs a chromosome.");
            }
        }
    }

    /// <summary>
    /// Battenberg <c>bkps_to_presegment_breakpoints(…, addin_bigholes = T)</c> on one chromosome; R <c>NA</c> = null
    /// (a presegment with an <c>NA</c> bound selects no SNP).
    /// </summary>
    private static List<(long? Start, long? End)> BattenbergPresegments(long[] pos, List<long> bkps)
    {
        int n = pos.Length;
        long? At(int oneBased) => oneBased >= 1 && oneBased <= n ? pos[oneBased - 1] : null; // R x[i], NA beyond
        var result = new List<(long?, long?)>();

        // addin_bigholes over SNP positions sel (1-based indices into pos)
        long? AddBigHoles(List<int> sel, long? startpos)
        {
            for (int k = 0; k + 1 < sel.Count; k++)
            {
                if (pos[sel[k + 1] - 1] - pos[sel[k] - 1] >= BattenbergMaxSnpDistance)
                {
                    result.Add((startpos, pos[sel[k] - 1]));
                    startpos = pos[sel[k + 1] - 1];
                }
            }

            return startpos;
        }

        if (bkps.Count > 0)
        {
            long? startpos;
            int startFromSv;
            long? first = At(1);
            if (first.HasValue && first.Value < bkps[0])
            {
                startpos = first;
                startFromSv = 1;
            }
            else
            {
                // R: NA < bkp is NA ⇒ if() would stop; an empty chromosome never reaches here with data.
                startpos = bkps[0];
                startFromSv = 2;
            }

            foreach (int idx in RColonRange(startFromSv, bkps.Count))
            {
                long? sv = idx >= 1 && idx <= bkps.Count ? bkps[idx - 1] : null;
                var selected = new List<int>();
                if (startpos.HasValue && sv.HasValue)
                {
                    for (int j = 1; j <= n; j++)
                    {
                        if (pos[j - 1] >= startpos.Value && pos[j - 1] <= sv.Value)
                        {
                            selected.Add(j);
                        }
                    }
                }

                if (selected.Count > 0)
                {
                    startpos = AddBigHoles(selected, startpos);
                    int endIndex = selected[^1];
                    result.Add((startpos, pos[endIndex - 1]));
                    startpos = At(endIndex + 1);
                }
            }

            if (n > 0 && pos[n - 1] > bkps[^1])
            {
                result.Add((startpos, pos[n - 1]));
            }
        }
        else
        {
            long? startpos = At(1);
            var all = new List<int>(n);
            for (int j = 1; j <= n; j++)
            {
                all.Add(j);
            }

            startpos = AddBigHoles(all, startpos);
            result.Add((startpos, At(n)));
        }

        return result;
    }

    /// <summary>R <c>a:b</c> for integers (descending when a &gt; b).</summary>
    private static IEnumerable<int> RColonRange(int a, int b)
    {
        int step = a <= b ? 1 : -1;
        int count = Math.Abs(b - a) + 1;
        for (int k = 0; k < count; k++)
        {
            yield return a + k * step;
        }
    }

    /// <summary>Battenberg <c>run_pcf</c> on one presegment; appends its rows to <paramref name="output"/>.</summary>
    private static void BattenbergRunPcf(
        string chr, long[] pos, double[] bafAll, long? start, long? end,
        BattenbergPhasedSegmentationOptions o, List<PhasedBafSegmentedSnp> output)
    {
        var idx = new List<int>();
        if (start.HasValue && end.HasValue)
        {
            for (int j = 0; j < pos.Length; j++)
            {
                if (pos[j] >= start.Value && pos[j] <= end.Value)
                {
                    idx.Add(j);
                }
            }
        }

        int n = idx.Count;
        if (n == 0)
        {
            return;
        }

        var baf = new double[n];
        var folded = new double[n];
        for (int k = 0; k < n; k++)
        {
            baf[k] = bafAll[idx[k]];
            folded[k] = ToAscatBaf(baf[k]); // min(b, 1 − b); equal at 0.5
        }

        double sdev = GetMad(folded, 25);
        if (double.IsNaN(sdev))
        {
            sdev = 0.0;
        }

        if (sdev < BattenbergMinBafSd)
        {
            sdev = BattenbergMinBafSd;
        }

        double[] bafSegm = n < BattenbergMinPcfSnps
            ? Filled(n, StatisticsHelper.RMean(baf))
            : BattenbergSelectFastPcf(baf, o.PhaseKmin, o.PhaseGamma * sdev);

        var phased = new double[n];
        for (int k = 0; k < n; k++)
        {
            phased[k] = OrientBaf(bafSegm[k], baf[k]);
        }

        double[] phSeg = n < BattenbergMinPcfSnps || o.NoSegmentation
            ? Filled(n, StatisticsHelper.RMean(phased))
            : BattenbergSelectFastPcf(phased, o.Kmin, o.Gamma * sdev);

        if (o.SegmentBafOption != BattenbergSegmentBafOption.Mean)
        {
            double[] median = BattenbergAdjustSegmValues(phased, phSeg);
            for (int k = 0; k < n; k++)
            {
                bool keepMean = o.SegmentBafOption == BattenbergSegmentBafOption.MedianUnlessExtreme
                                && (median[k] == 0.0 || median[k] == 1.0);
                phSeg[k] = keepMean ? phSeg[k] : median[k];
            }
        }

        for (int k = 0; k < n; k++)
        {
            output.Add(new PhasedBafSegmentedSnp(chr, pos[idx[k]], baf[k], phased[k], phSeg[k]));
        }
    }

    /// <summary>Battenberg <c>adjustSegmValues</c>: each run of equal <c>BAFseg</c> (R <c>rle</c>) → median of its BAFphased.</summary>
    private static double[] BattenbergAdjustSegmValues(double[] phased, double[] seg)
    {
        int n = seg.Length;
        var result = new double[n];
        int from = 0;
        for (int i = 0; i < n; i++)
        {
            if (i < n - 1 && seg[i + 1] == seg[i])
            {
                continue;
            }

            double median = StatisticsHelper.Median(new ArraySegment<double>(phased, from, i - from + 1));
            for (int k = from; k <= i; k++)
            {
                result[k] = median;
            }

            from = i + 1;
        }

        return result;
    }

    private static double[] Filled(int n, double value)
    {
        var a = new double[n];
        Array.Fill(a, value);
        return a;
    }


    /// <summary>
    /// Battenberg/copynumber <c>selectFastPcf(x, kmin, gamma, yest = T)$yhat</c> (R/fastPCF.R): <c>runFastPcf</c> with
    /// filter fractions (0.15, 0.15) below 1000 values, (0.12, 0.05) below 15 000, else <c>runPcfSubset</c>.
    /// </summary>
    private static double[] BattenbergSelectFastPcf(double[] x, int kmin, double gamma)
    {
        int n = x.Length;
        if (n < 1000)
        {
            return BattenbergRunFastPcf(x, kmin, gamma, 0.15, 0.15);
        }

        return n < 15000
            ? BattenbergRunFastPcf(x, kmin, gamma, 0.12, 0.05)
            : BattenbergRunPcfSubset(x, kmin, gamma, 0.12, 0.05);
    }

    private static double[] BattenbergRunFastPcf(double[] x, int kmin, double gamma, double frac1, double frac2)
    {
        bool[] mark = BattenbergFilterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9);
        mark[^1] = true;
        (int[] nr, double[] sum, double[] sq) = BattenbergCompact(x, x.Length, mark);
        return BattenbergPottsCompact(kmin, gamma, nr, sum, sq);
    }

    /// <summary>copynumber <c>runPcfSubset</c>: Potts marks refined over 5000-value windows advancing by 4000.</summary>
    private static double[] BattenbergRunPcfSubset(double[] x, int kmin, double gamma, double frac1, double frac2)
    {
        const int subSize = 5000;
        int antGen = x.Length;
        bool[] mark = BattenbergFilterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9);
        var markInit = new bool[subSize];
        Array.Copy(mark, markInit, subSize - 1);
        markInit[subSize - 1] = true;
        (int[] nr, double[] sum, double[] sq) = BattenbergCompact(x, subSize, markInit);
        var mark2 = new bool[antGen];
        Array.Copy(BattenbergMarkWithPotts(kmin, gamma, nr, sum, sq, subSize), mark2, subSize);
        mark2[(4 * subSize / 5) - 1] = true;
        int start = (4 * subSize / 5) + 1; // 1-based
        while (start + subSize < antGen)
        {
            int slutt = start + subSize - 1;
            var markSub = new bool[slutt];
            Array.Copy(mark2, markSub, start - 1);
            Array.Copy(mark, start - 1, markSub, start - 1, slutt - start + 1);
            markSub[slutt - 1] = true;
            (nr, sum, sq) = BattenbergCompact(x, slutt, markSub);
            Array.Copy(BattenbergMarkWithPotts(kmin, gamma, nr, sum, sq, slutt), mark2, slutt);
            start += 4 * subSize / 5;
            mark2[start - 2] = true;
        }

        var finalMark = new bool[antGen];
        Array.Copy(mark2, finalMark, start - 1);
        Array.Copy(mark, start - 1, finalMark, start - 1, antGen - start + 1);
        (nr, sum, sq) = BattenbergCompact(x, antGen, finalMark);
        return BattenbergPottsCompact(kmin, gamma, nr, sum, sq);
    }

    /// <summary>
    /// copynumber <c>compact(y[1:length], mark)</c>: counts, sums and sums of squares between marked positions, from the
    /// running <c>cumsum</c>s (differences of the cumulative sums at the marks, as R does).
    /// </summary>
    private static (int[] Nr, double[] Sum, double[] Sq) BattenbergCompact(double[] y, int length, bool[] mark)
    {
        var nr = new List<int>();
        var sum = new List<double>();
        var sq = new List<double>();
        double cy = 0.0, cy2 = 0.0, lowCy = 0.0, lowCy2 = 0.0;
        int low = 0;
        for (int i = 0; i < length; i++)
        {
            cy += y[i];
            cy2 += y[i] * y[i];
            if (mark[i])
            {
                nr.Add(i + 1 - low);
                sum.Add(cy - lowCy);
                sq.Add(cy2 - lowCy2);
                low = i + 1;
                lowCy = cy;
                lowCy2 = cy2;
            }
        }

        return (nr.ToArray(), sum.ToArray(), sq.ToArray());
    }

    /// <summary>
    /// copynumber <c>PottsCompact(kmin, gamma, nr, res, sq, yest = T)$yhat</c>: exact Potts filtering on the compacted
    /// array (1-based arrays below mirror the R code line by line), then <c>findEst</c>.
    /// </summary>
    private static double[] BattenbergPottsCompact(int kmin, double gamma, int[] nr, double[] res, double[] sq)
    {
        int bigN = nr.Length;
        long total = 0;
        double totalSum = 0.0;
        for (int i = 0; i < bigN; i++)
        {
            total += nr[i];
            totalSum += res[i];
        }

        if (total < 2 * kmin)
        {
            return Filled((int)total, totalSum / total); // R returns the scalar estimate
        }

        int[] bestSplit = BattenbergPottsSplits(kmin, gamma, nr, res, sq, markSub: null);
        return BattenbergFindEst(bestSplit, bigN, nr, res);
    }

    /// <summary>
    /// The Potts recursion shared by copynumber <c>PottsCompact</c> and <c>markWithPotts</c>; returns 1-based
    /// <c>bestSplit</c> and, when <paramref name="markSub"/> is given, sets <c>markSub[Pos − 1]</c> as markWithPotts does.
    /// </summary>
    private static int[] BattenbergPottsSplits(int kmin, double gamma, int[] nr, double[] res, double[] sq, bool[]? markSub)
    {
        int bigN = nr.Length;
        var ant = new double[bigN + 1];
        var sum = new double[bigN + 1];
        var kvad = new double[bigN + 1];
        var cost = new double[bigN + 1];
        var bestCost = new double[bigN + 1];
        var bestSplit = new int[bigN + 1];
        double initAnt = nr[0], initSum = res[0], initKvad = sq[0];
        double initAve = initSum / initAnt;
        bestCost[1] = initKvad - initSum * initAve;
        int k = 2;
        double cum = bigN > 1 ? nr[0] + nr[1] : nr[0];
        while (cum < 2 * kmin && k < bigN)
        {
            for (int j = 2; j <= k; j++)
            {
                ant[j] += nr[k - 1];
                sum[j] += res[k - 1];
                kvad[j] += sq[k - 1];
            }

            double s2 = initSum + sum[2];
            bestCost[k] = (initKvad + kvad[2]) - s2 * s2 / (initAnt + ant[2]);
            k++;
            cum += nr[k - 1];
        }

        for (int n = k; n <= bigN; n++)
        {
            for (int j = 2; j <= n; j++)
            {
                ant[j] += nr[n - 1];
                sum[j] += res[n - 1];
                kvad[j] += sq[n - 1];
            }

            int limit = n;
            while (limit > 2 && ant[limit] < kmin)
            {
                limit--;
            }

            int pos = -1;
            double best = double.NaN;
            for (int j = 2; j <= limit; j++)
            {
                cost[j] = bestCost[j - 1] + kvad[j] - sum[j] * sum[j] / ant[j];
                if (!double.IsNaN(cost[j]) && (pos < 0 || cost[j] < best))
                {
                    pos = j;
                    best = cost[j];
                }
            }

            double c = best + gamma;
            double st = sum[2] + initSum;
            double totCost = (kvad[2] + initKvad) - st * st / (ant[2] + initAnt);
            if (totCost < c)
            {
                pos = 1;
                c = totCost;
            }

            bestCost[n] = c;
            bestSplit[n] = pos - 1;
            if (markSub is not null && pos - 1 >= 1)
            {
                markSub[pos - 2] = true; // R markSub[Pos − 1] (1-based)
            }
        }

        return bestSplit;
    }

    /// <summary>copynumber <c>findEst(bestSplit, N, Nr, Sum, yest = T)$yhat</c>.</summary>
    private static double[] BattenbergFindEst(int[] bestSplit, int bigN, int[] nr, double[] sum)
    {
        var lengths = new List<int>();
        int n = bigN;
        while (n > 0)
        {
            lengths.Add(n - bestSplit[n]);
            n = bestSplit[n];
        }

        lengths.Reverse();
        var yhat = new List<double>();
        int start = 0; // 0-based into the compact arrays
        foreach (int len in lengths)
        {
            int lengdeOrig = 0;
            double s = 0.0;
            for (int i = start; i < start + len; i++)
            {
                lengdeOrig += nr[i];
                s += sum[i];
            }

            double verdi = s / lengdeOrig;
            for (int i = 0; i < lengdeOrig; i++)
            {
                yhat.Add(verdi);
            }

            start += len;
        }

        return yhat.ToArray();
    }

    /// <summary>copynumber <c>markWithPotts</c> + <c>findMarks</c>: Potts split points mapped back to original indices.</summary>
    private static bool[] BattenbergMarkWithPotts(int kmin, double gamma, int[] nr, double[] res, double[] sq, int subSize)
    {
        var markSub = new bool[nr.Length];
        BattenbergPottsSplits(kmin, gamma, nr, res, sq, markSub);
        var mark = new bool[subSize];
        int orig = 0;
        for (int i = 0; i < markSub.Length; i++)
        {
            orig += nr[i];
            if (markSub[i])
            {
                mark[orig - 1] = true; // findMarks: mark[startOrig − 1] at the end of compact cell i
            }
        }

        return mark;
    }

    /// <summary>
    /// copynumber <c>filterMarkS4(x, kmin, L, L2, frac1, frac2, frac3, thres)</c>: candidate breakpoints from two
    /// high-pass filters (widths 6L, 6L2) and a kmin-segment filter, with type-7 quantile limits
    /// (<see cref="StatisticsHelper.SampleQuantileType7"/>). Arrays are 1-based (index 0 unused) as in R.
    /// </summary>
    private static bool[] BattenbergFilterMarkS4(
        double[] x, int kmin, int bigL, int bigL2, double frac1, double frac2, double frac3, double thres)
    {
        int n = x.Length;
        var xc = new double[n + 1]; // R xc = c(0, cumsum(x)); xc[j] (1-based R) = xc[j − 1] here
        for (int i = 0; i < n; i++)
        {
            xc[i + 1] = xc[i] + x[i];
        }

        double Xc(int rIndex) => xc[rIndex - 1];

        // cost1 (length n, 1-based)
        var cost1 = new double[n + 1];
        for (int t = 1; t <= n - 6 * bigL + 1; t++)
        {
            cost1[3 * bigL - 1 + t] = Math.Abs(
                4 * Xc(t + 3 * bigL) - Xc(t) - Xc(t + bigL) - Xc(t + 5 * bigL) - Xc(t + 6 * bigL));
        }

        double[] test = SevenMax(cost1, n);
        var cost1B = new List<double>();
        for (int j = 1; j <= n; j++)
        {
            if (cost1[j] >= thres * test[j])
            {
                cost1B.Add(cost1[j]);
            }
        }

        double frac1B = Math.Min(0.8, frac1 * n / cost1B.Count);
        double limit = StatisticsHelper.SampleQuantileType7(cost1B, 1 - frac1B);
        var mark = new bool[n + 1];
        for (int j = 1; j <= n; j++)
        {
            mark[j] = cost1[j] > limit && cost1[j] > 0.9 * test[j];
        }

        int m2 = n - 6 * bigL2 + 1;
        var cost2 = new double[m2];
        for (int t = 1; t <= m2; t++)
        {
            cost2[t - 1] = Math.Abs(
                4 * Xc(t + 3 * bigL2) - Xc(t) - Xc(t + bigL2) - Xc(t + 5 * bigL2) - Xc(t + 6 * bigL2));
        }

        double limit2 = StatisticsHelper.SampleQuantileType7(cost2, 1 - frac2);
        var mark2 = new bool[n + 1];
        for (int t = 1; t <= m2; t++)
        {
            mark2[3 * bigL2 - 1 + t] = cost2[t - 1] > limit2;
        }

        if (3 * bigL > kmin)
        {
            SetRange(mark, kmin, 3 * bigL - 1, true);
            SetRange(mark, n - 3 * bigL + 1, n - kmin, true);
        }
        else
        {
            SetRange(mark, kmin, kmin, true);
            SetRange(mark, n - kmin, n - kmin, true);
        }

        if (kmin > 1)
        {
            int m = n - 3 * kmin + 1;
            var shortAb = new double[m + 1];
            for (int t = 1; t <= m; t++)
            {
                shortAb[t] = Math.Abs(3 * (Xc(t + 2 * kmin) - Xc(t + kmin)) - (Xc(t + 3 * kmin) - Xc(t)));
            }

            double[] test2 = SevenMax(shortAb, m);
            var cost1C = new List<double>();
            for (int t = 1; t <= m; t++)
            {
                if (shortAb[t] >= thres * test2[t])
                {
                    cost1C.Add(shortAb[t]);
                }
            }

            double frac1C = Math.Min(0.8, frac3 * m / cost1C.Count);
            double limit3 = StatisticsHelper.SampleQuantileType7(cost1C, 1 - frac1C);
            for (int t = 1; t <= m; t++)
            {
                if (shortAb[t] > limit3 && shortAb[t] > thres * test2[t])
                {
                    mark[kmin - 1 + t] = true;     // markH2
                    mark[2 * kmin - 1 + t] = true; // markH3
                }
            }
        }

        for (int j = 1; j <= n; j++)
        {
            mark[j] |= mark2[j];
        }

        if (3 * bigL > kmin)
        {
            SetRange(mark, 1, kmin - 1, false);
            SetRange(mark, kmin, 3 * bigL - 1, true);
            SetRange(mark, n - 3 * bigL + 1, n - kmin, true);
            SetRange(mark, n - kmin + 1, n - 1, false);
            mark[n] = true;
        }
        else
        {
            SetRange(mark, 1, kmin - 1, false);
            SetRange(mark, n - kmin + 1, n - 1, false);
            mark[n] = true;
            mark[kmin] = true;
            mark[n - kmin] = true;
        }

        var result = new bool[n];
        Array.Copy(mark, 1, result, 0, n);
        return result;
    }

    /// <summary>
    /// filterMarkS4's <c>c(rep(0, 3), pmax(v[i..i+6]), rep(0, 3))</c> over a 1-based vector v of length len.
    /// </summary>
    private static double[] SevenMax(double[] v, int len)
    {
        var test = new double[len + 1];
        for (int i = 1; i <= len - 6; i++)
        {
            double mx = v[i];
            for (int d = 1; d <= 6; d++)
            {
                mx = Math.Max(mx, v[i + d]);
            }

            test[3 + i] = mx;
        }

        return test;
    }

    /// <summary>R <c>mark[a:b] &lt;- value</c> on a 1-based array (descending when a &gt; b; index 0 ignored).</summary>
    private static void SetRange(bool[] mark, int a, int b, bool value)
    {
        foreach (int i in RColonRange(a, b))
        {
            if (i >= 1)
            {
                mark[i] = value;
            }
        }
    }
}
