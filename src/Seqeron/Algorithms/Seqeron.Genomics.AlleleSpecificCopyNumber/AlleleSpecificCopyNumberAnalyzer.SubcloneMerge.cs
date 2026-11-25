// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Battenberg segment merging, high copy-number masking and the callSubclones driver (port of Battenberg 3.0.0).</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{

    /// <summary>Battenberg <c>callSubclones</c> default <c>max_allowed_state = 250</c> (<c>mask_high_cn_segments</c>).</summary>
    public const int BattenbergMaxAllowedState = 250;

    /// <summary>Battenberg <c>merge_segments</c>: segments farther apart than 3 Mb (<c>GenomicRanges::distance</c>) are never merged.</summary>
    private const long BattenbergMergeMaxDistance = 3_000_000;

    /// <summary>Battenberg <c>merge_segments</c>: both segments need more than 10 non-NA logR and BAF values for the t-tests.</summary>
    private const int BattenbergMergeMinValues = 10;

    /// <summary>Battenberg <c>merge_segments</c>: significance level of the two Welch t-tests (logR and BAF).</summary>
    private const double BattenbergMergeSignificance = 0.05;

    /// <summary>
    /// One row of Battenberg's merged <c>subclones</c> table (<c>merge_segments</c> output): the segment extent (the larger
    /// neighbour extended over the absorbed segments, gaps included), its <c>BAF</c> and <c>LogR</c> (recomputed after a
    /// merge, else the <c>determine_copynumber</c> values) and the <c>determine_copynumber</c> row whose remaining
    /// columns it keeps (Battenberg only updates <c>start</c>/<c>end</c>/<c>BAF</c>/<c>LogR</c> of the surviving row).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Start">startpos.</param>
    /// <param name="End">endpos.</param>
    /// <param name="Baf">Battenberg <c>BAF</c>: the mirrored level l of an unmerged row, or the merged segment's
    /// <c>calc_seg_baf_option</c> summary of its <c>BAFphased</c> (not mirrored, as in Battenberg).</param>
    /// <param name="LogR">Battenberg <c>LogR</c>: segment mean logR (merged: mean of the non-NA probes in [Start, End], 0 if
    /// there are none).</param>
    /// <param name="Call">The <c>determine_copynumber</c> row this one descends from (copy-number states, fractions, p-value).</param>
    /// <param name="WasMerged">True when at least one neighbour was merged into this segment.</param>
    public sealed record BattenbergMergedSegment(
        string Chromosome, long Start, long End, double Baf, double LogR, BattenbergSegmentCall Call, bool WasMerged);

    /// <summary>Battenberg <c>merge_segments</c> output: the merged <c>subclones</c> rows and the updated <c>BAFsegmented</c>.</summary>
    /// <param name="Segments">Merged segments, chromosomes in <c>GenomeInfoDb</c> seqlevel order (see
    /// <see cref="MergeBattenbergSegments"/>), segments in genome order.</param>
    /// <param name="SegmentedSnps">The <c>BAFsegmented</c> rows, chromosomes in the same order (rows of a chromosome in input
    /// order), with <c>BAFseg</c> of every merged segment replaced by its new BAF.</param>
    public sealed record BattenbergSegmentMerge(
        IReadOnlyList<BattenbergMergedSegment> Segments,
        IReadOnlyList<PhasedBafSegmentedSnp> SegmentedSnps);

    /// <summary>Battenberg <c>mask_high_cn_segments</c> output.</summary>
    /// <param name="Calls">The calls, masked rows flagged <see cref="BattenbergSegmentCall.IsMasked"/>.</param>
    /// <param name="SegmentedSnps">The <c>BAFsegmented</c> rows with <c>BAFseg</c> = NaN (R <c>NA</c>) for the SNPs of masked
    /// segments with <c>startpos &lt; Position ≤ endpos</c> (Battenberg's half-open test: a masked segment's first SNP keeps
    /// its BAFseg).</param>
    /// <param name="MaskedCount">Battenberg <c>masked_count</c>.</param>
    /// <param name="MaskedSize">Battenberg <c>masked_size</c> = Σ (endpos − startpos) over masked segments.</param>
    public sealed record BattenbergHighCopyNumberMask(
        IReadOnlyList<BattenbergSegmentCall> Calls,
        IReadOnlyList<PhasedBafSegmentedSnp> SegmentedSnps,
        int MaskedCount,
        long MaskedSize);

    /// <summary>
    /// Result of <see cref="CallBattenbergSubclones"/>: the three stages of Battenberg <c>callSubclones</c>.
    /// </summary>
    /// <param name="InitialCalls">First <c>determine_copynumber</c> (Battenberg's <c>…_1.txt</c> table).</param>
    /// <param name="Merge"><c>merge_segments</c> of the initial calls.</param>
    /// <param name="Mask"><c>mask_high_cn_segments</c> of the second <c>determine_copynumber</c> on the merged
    /// <c>BAFsegmented</c>; <see cref="Calls"/> is Battenberg's final (extended) copy-number table.</param>
    public sealed record BattenbergSubcloneCalls(
        IReadOnlyList<BattenbergSegmentCall> InitialCalls,
        BattenbergSegmentMerge Merge,
        BattenbergHighCopyNumberMask Mask)
    {
        /// <summary>The final per-segment calls (second <c>determine_copynumber</c>, masked).</summary>
        public IReadOnlyList<BattenbergSegmentCall> Calls => Mask.Calls;
    }

    /// <summary>
    /// Battenberg <c>merge_segments</c> (Wedge-lab/battenberg R/fitcopynumber.R; B24 F60): merges adjacent segments of a
    /// fitted copy-number profile when there is not enough evidence for them to be separate. Per chromosome, repeatedly
    /// the <b>smallest</b> segment (width end − start + 1; first on a tie) that still has an unchecked neighbour is taken,
    /// and its neighbours are visited closest first (previous one first on a tie):
    /// <list type="number">
    /// <item>already checked → skipped; farther than 3 Mb (<c>GenomicRanges::distance</c> = start₂ − end₁ − 1) → not merged;</item>
    /// <item>same clonal solution (<c>nMaj1_A</c>, <c>nMin1_A</c> equal and <c>frac1_A == 1</c> in both) → merged;</item>
    /// <item>else, when <c>round(nmin)</c> or <c>round(nmaj)</c> agree, with
    /// nmin = (ρ − 1 − (BAF − 1)·2^(LogR/γ)·((1 − ρ)·2 + ρ·ψ))/ρ and nmaj = (ρ − 1 + BAF·2^(LogR/γ)·((1 − ρ)·2 + ρ·ψ))/ρ
    /// (ψ = <c>callSubclones</c>' psi of all cells, ρ·psit + 2(1 − ρ) — Battenberg's formula, used verbatim), and both
    /// segments have &gt; 10 non-NA logR probes and SNP BAFs in their ranges: merged when neither Welch
    /// <c>t.test</c> (logR vs logR, <c>BAFphased</c> vs <c>BAFphased</c>; <see cref="StatisticsHelper.WelchTTestPValue"/>)
    /// has p &lt; 0.05; otherwise not merged.</item>
    /// </list>
    /// A merge extends the neighbour over the smaller segment, recomputes its BAF from the <c>BAFphased</c> in the new range
    /// (<paramref name="bafOption"/>: median, mean, or median unless it is 0 or 1) and LogR (mean of the non-NA probes, 0
    /// when none), writes the BAF into those rows' <c>BAFseg</c>, keeps the neighbour's copy-number columns, and re-opens
    /// both sides of the merged segment for checking. Rounding is R's (half to even); the input BAF/LogR/frac1 pass
    /// through R's 15-significant-digit <c>as.character</c> round trip of <c>determine_copynumber</c>'s character matrix.
    /// Output chromosomes are ordered as <c>makeGRangesFromDataFrame</c> orders seqlevels (<c>GenomeInfoDb::rankSeqlevels</c>:
    /// optional <c>CHR</c>/<c>chr</c>/<c>CH</c>/<c>ch</c> prefix, then numbers, X, Y, M, MT, other names), which changes
    /// the row order of <c>BAFsegmented</c> for unsorted input (and with it the order of the second
    /// <c>determine_copynumber</c>).
    /// </summary>
    /// <param name="calls">The fitted profile (<c>determine_copynumber</c> rows, e.g. from
    /// <see cref="FitSubclonalCopyNumberWithBootstrap"/>); rows on chromosomes absent from
    /// <paramref name="segmentedSnps"/> are dropped, as in Battenberg.</param>
    /// <param name="segmentedSnps">The <c>BAFsegmented</c> rows the profile was fitted from (BAFphased, BAFseg in [0, 1]).</param>
    /// <param name="logR">Raw logR probes; NaN = R <c>NA</c> (ignored by the means and tests, not counted).</param>
    /// <param name="purity">ρ ∈ (0, 1] the profile was fitted with.</param>
    /// <param name="ploidy">Tumour ploidy psit (&gt; 0) the profile was fitted with.</param>
    /// <param name="gamma">Platform γ (&gt; 0).</param>
    /// <param name="bafOption">Battenberg <c>calc_seg_baf_option</c> (default 3, median unless extreme).</param>
    /// <returns>The merged segments and the updated <c>BAFsegmented</c> rows.</returns>
    /// <exception cref="ArgumentNullException">an argument, a chromosome or a call is null.</exception>
    /// <exception cref="ArgumentException">a BAF outside [0, 1]; an infinite logR (R's tests and means fail); a chromosome
    /// of <paramref name="segmentedSnps"/> without any call or logR probe (R <c>stopifnot</c>).</exception>
    /// <exception cref="ArgumentOutOfRangeException">invalid ρ, ψ, γ or an undefined <paramref name="bafOption"/>.</exception>
    /// <exception cref="InvalidOperationException">a merge t-test whose data are essentially constant (R <c>t.test</c> stops).</exception>
    public static BattenbergSegmentMerge MergeBattenbergSegments(
        IReadOnlyList<BattenbergSegmentCall> calls,
        IReadOnlyList<PhasedBafSegmentedSnp> segmentedSnps,
        IReadOnlyList<LogRProbe> logR,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma,
        BattenbergSegmentBafOption bafOption = BattenbergSegmentBafOption.MedianUnlessExtreme)
    {
        ArgumentNullException.ThrowIfNull(calls);
        ArgumentNullException.ThrowIfNull(segmentedSnps);
        ArgumentNullException.ThrowIfNull(logR);
        ValidateAscatModelParameters(purity, ploidy, gamma);
        if (!Enum.IsDefined(bafOption))
        {
            throw new ArgumentOutOfRangeException(nameof(bafOption), bafOption, "Undefined calc_seg_baf_option.");
        }

        var snpsByChr = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        var firstAppearance = new List<string>();
        for (int i = 0; i < segmentedSnps.Count; i++)
        {
            PhasedBafSegmentedSnp row = segmentedSnps[i];
            if (row.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(segmentedSnps), "Every row needs a chromosome.");
            }

            if (!(row.BafPhased >= 0.0 && row.BafPhased <= 1.0) || !(row.BafSegment >= 0.0 && row.BafSegment <= 1.0))
            {
                throw new ArgumentException("Every BAFseg and BAFphased must lie in [0, 1].", nameof(segmentedSnps));
            }

            if (!snpsByChr.TryGetValue(row.Chromosome, out var list))
            {
                list = new List<int>();
                snpsByChr.Add(row.Chromosome, list);
                firstAppearance.Add(row.Chromosome);
            }

            list.Add(i);
        }

        var probesByChr = new Dictionary<string, List<(long Pos, double LogR)>>(StringComparer.Ordinal);
        foreach (LogRProbe p in logR)
        {
            if (p.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(logR), "Every logR probe needs a chromosome.");
            }

            if (double.IsInfinity(p.LogR))
            {
                throw new ArgumentException("merge_segments needs finite logR values (NaN = NA is allowed).", nameof(logR));
            }

            if (!probesByChr.TryGetValue(p.Chromosome, out var list))
            {
                list = new List<(long, double)>();
                probesByChr.Add(p.Chromosome, list);
            }

            list.Add((p.Position, p.LogR));
        }

        var segsByChr = new Dictionary<string, List<BattenbergMergeSegment>>(StringComparer.Ordinal);
        foreach (BattenbergSegmentCall call in calls)
        {
            if (call is null)
            {
                throw new ArgumentNullException(nameof(calls), "Every call must be non-null.");
            }

            AlleleSpecificSegmentSummary s = call.Fit.Segment;
            if (s.Chromosome is null)
            {
                throw new ArgumentNullException(nameof(calls), "Every call needs a chromosome.");
            }

            if (!segsByChr.TryGetValue(s.Chromosome, out var list))
            {
                list = new List<BattenbergMergeSegment>();
                segsByChr.Add(s.Chromosome, list);
            }

            list.Add(new BattenbergMergeSegment(call)
            {
                Start = s.Start,
                End = s.End,
                Baf = RAsCharacterRoundTrip(call.Baf),
                LogR = RAsCharacterRoundTrip(s.MeanLogR),
                Frac1 = RAsCharacterRoundTrip(call.Fit.PrimaryState.CellFraction),
            });
        }

        double psiAll = MixtureCopiesPerCell(purity, ploidy);
        double cellCopies = MixtureCopiesPerCell(purity, psiAll); // Battenberg ((1-rho)*2+rho*psi), psi of all cells
        var bafSegOut = new double[segmentedSnps.Count];
        for (int i = 0; i < bafSegOut.Length; i++)
        {
            bafSegOut[i] = segmentedSnps[i].BafSegment;
        }

        var outSegments = new List<BattenbergMergedSegment>();
        var outRows = new List<PhasedBafSegmentedSnp>(segmentedSnps.Count);
        foreach (string chr in BattenbergSeqlevelOrder(firstAppearance))
        {
            if (!segsByChr.TryGetValue(chr, out var segs) || !probesByChr.TryGetValue(chr, out var probes))
            {
                throw new ArgumentException(
                    $"Chromosome '{chr}' needs at least one call and one logR probe (Battenberg merge_segments stopifnot).",
                    nameof(segmentedSnps));
            }

            var ctx = new BattenbergMergeChromosome(segmentedSnps, snpsByChr[chr], probes, bafSegOut, bafOption);
            BattenbergMergeChromosomeSegments(segs, ctx, purity, cellCopies, gamma);
            foreach (BattenbergMergeSegment m in segs)
            {
                outSegments.Add(new BattenbergMergedSegment(chr, m.Start, m.End, m.Baf, m.LogR, m.Call, m.WasMerged));
            }

            foreach (int i in snpsByChr[chr])
            {
                outRows.Add(segmentedSnps[i] with { BafSegment = bafSegOut[i] });
            }
        }

        return new BattenbergSegmentMerge(outSegments, outRows);
    }

    /// <summary>
    /// Battenberg <c>mask_high_cn_segments</c> (R/fitcopynumber.R; B24 F60): a segment whose solution-A state 1 has
    /// <c>nMaj1_A &gt; max_allowed_state</c> or <c>nMin1_A &gt; max_allowed_state</c> is masked — <c>nMaj1_A</c>,
    /// <c>nMin1_A</c>, <c>nMaj2_A</c>, <c>nMin2_A</c> become NA (flag <see cref="BattenbergSegmentCall.IsMasked"/>; solution A's
    /// copy numbers null) and the <c>BAFseg</c> of its SNPs with <c>startpos &lt; Position ≤ endpos</c> becomes NA; the
    /// count and total size Σ(endpos − startpos) are reported. "In part an artifact of small segments" (Battenberg).
    /// </summary>
    /// <param name="calls">The fitted profile (<c>determine_copynumber</c> rows).</param>
    /// <param name="segmentedSnps">The <c>BAFsegmented</c> rows (BAFseg masked in the output).</param>
    /// <param name="maxAllowedState">Battenberg <c>max_allowed_state</c> (default 250).</param>
    /// <returns>The masked calls and rows, <c>masked_count</c> and <c>masked_size</c>.</returns>
    /// <exception cref="ArgumentNullException">an argument, a call or a chromosome is null.</exception>
    public static BattenbergHighCopyNumberMask MaskHighCopyNumberSegments(
        IReadOnlyList<BattenbergSegmentCall> calls,
        IReadOnlyList<PhasedBafSegmentedSnp> segmentedSnps,
        int maxAllowedState = BattenbergMaxAllowedState)
    {
        ArgumentNullException.ThrowIfNull(calls);
        ArgumentNullException.ThrowIfNull(segmentedSnps);
        var rows = segmentedSnps.ToArray();
        var outCalls = new List<BattenbergSegmentCall>(calls.Count);
        int count = 0;
        long size = 0;
        foreach (BattenbergSegmentCall call in calls)
        {
            if (call is null)
            {
                throw new ArgumentNullException(nameof(calls), "Every call must be non-null.");
            }

            SubclonalCopyNumberState a = call.Fit.PrimaryState;
            if (!(a.MajorCopyNumber > maxAllowedState || a.MinorCopyNumber > maxAllowedState))
            {
                outCalls.Add(call);
                continue;
            }

            AlleleSpecificSegmentSummary s = call.Fit.Segment;
            for (int i = 0; i < rows.Length; i++)
            {
                if (rows[i].Chromosome == s.Chromosome && s.Start < rows[i].Position && s.End >= rows[i].Position)
                {
                    rows[i] = rows[i] with { BafSegment = double.NaN };
                }
            }

            IReadOnlyList<BattenbergSubclonalSolution> solutions = call.Solutions;
            if (solutions.Count > 0)
            {
                var masked = solutions.ToArray();
                masked[0] = masked[0] with
                {
                    MajorCopyNumber1 = null, MinorCopyNumber1 = null, MajorCopyNumber2 = null, MinorCopyNumber2 = null,
                };
                solutions = masked;
            }

            outCalls.Add(call with { Solutions = solutions, IsMasked = true });
            count++;
            size += s.End - s.Start;
        }

        return new BattenbergHighCopyNumberMask(outCalls, rows, count, size);
    }

    /// <summary>
    /// Battenberg <c>callSubclones</c> copy-number driver (R/fitcopynumber.R, B24 F61): <c>set.seed(seed)</c> →
    /// <c>determine_copynumber</c> (<see cref="BuildBattenbergSegments"/> + <see cref="FitSubclonalCopyNumberWithBootstrap"/>)
    /// → <see cref="MergeBattenbergSegments"/> → <c>determine_copynumber</c> on the merged <c>BAFsegmented</c> →
    /// <see cref="MaskHighCopyNumberSegments"/>. Both <c>determine_copynumber</c> calls draw their bootstrap resamples
    /// from <b>one</b> R Mersenne-Twister stream seeded once, as in Battenberg, so the final table (including
    /// <c>SDfrac_BS</c> and the bootstrap bounds) is Battenberg's for the same seed. Plots, file output, the PGA-clonal
    /// statistic and the ploidy recalculation of <c>callSubclones</c> are not part of this method. ψ of all cells is
    /// derived from <paramref name="ploidy"/> (psit) as in <c>callSubclones</c>: <c>psi = rho·psit + 2(1 − rho)</c>.
    /// </summary>
    /// <param name="segmentedSnps">Battenberg <c>BAFsegmented</c> rows (e.g. <see cref="SegmentPhasedBaf"/>).</param>
    /// <param name="logR">Raw logR probes (NaN = NA; ±∞ rejected by the merge step).</param>
    /// <param name="purity">ρ ∈ (0, 1] (<c>rho.psi.file</c> rho).</param>
    /// <param name="ploidy">Tumour ploidy psit (&gt; 0).</param>
    /// <param name="seed">R <c>set.seed</c> value (Battenberg <c>seed</c>).</param>
    /// <param name="gamma">Platform γ (&gt; 0).</param>
    /// <param name="significanceLevel">Battenberg <c>siglevel</c> (default 0.05).</param>
    /// <param name="maxBafDistance">Battenberg <c>maxdist</c> (default 0.01).</param>
    /// <param name="permutations">Battenberg <c>noperms</c> (default 1000).</param>
    /// <param name="maxAllowedState">Battenberg <c>max_allowed_state</c> (default 250).</param>
    /// <param name="bafOption">Battenberg <c>calc_seg_baf_option</c> used by the merge (default 3).</param>
    /// <returns>The initial calls, the merge, and the masked final calls.</returns>
    /// <exception cref="ArgumentNullException">an argument is null.</exception>
    /// <exception cref="ArgumentException">as <see cref="BuildBattenbergSegments"/> / <see cref="MergeBattenbergSegments"/>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">invalid ρ, ψ, γ, siglevel, maxdist, permutations &lt; 1 or bafOption.</exception>
    /// <exception cref="InvalidOperationException">a merge t-test on essentially constant data (R stops).</exception>
    public static BattenbergSubcloneCalls CallBattenbergSubclones(
        IReadOnlyList<PhasedBafSegmentedSnp> segmentedSnps,
        IReadOnlyList<LogRProbe> logR,
        double purity,
        double ploidy,
        int seed,
        double gamma = AscatSequencingGamma,
        double significanceLevel = BattenbergSignificanceLevel,
        double maxBafDistance = BattenbergMaxBafDistance,
        int permutations = BattenbergDefaultPermutations,
        int maxAllowedState = BattenbergMaxAllowedState,
        BattenbergSegmentBafOption bafOption = BattenbergSegmentBafOption.MedianUnlessExtreme)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(permutations, 1);
        var rng = new RMersenneTwister(seed);
        List<BattenbergSegmentCall> initial = BattenbergDetermineCopyNumber(
            BuildBattenbergSegments(segmentedSnps, logR), purity, ploidy, rng, gamma, significanceLevel, maxBafDistance, permutations);
        BattenbergSegmentMerge merge = MergeBattenbergSegments(initial, segmentedSnps, logR, purity, ploidy, gamma, bafOption);
        List<BattenbergSegmentCall> final = BattenbergDetermineCopyNumber(
            BuildBattenbergSegments(merge.SegmentedSnps, logR), purity, ploidy, rng, gamma, significanceLevel, maxBafDistance, permutations);
        return new BattenbergSubcloneCalls(initial, merge, MaskHighCopyNumberSegments(final, merge.SegmentedSnps, maxAllowedState));
    }

    /// <summary>A mutable <c>merge_segments</c> GRanges row (<c>ID</c> = list index, <c>Prev_checked</c>, <c>Next_checked</c>).</summary>
    private sealed class BattenbergMergeSegment(BattenbergSegmentCall call)
    {
        public BattenbergSegmentCall Call { get; } = call;

        public long Start { get; set; }

        public long End { get; set; }

        public double Baf { get; set; }

        public double LogR { get; set; }

        public double Frac1 { get; init; }

        public bool PrevChecked { get; set; }

        public bool NextChecked { get; set; }

        public bool WasMerged { get; set; }

        public int Major1 => Call.Fit.PrimaryState.MajorCopyNumber;

        public int Minor1 => Call.Fit.PrimaryState.MinorCopyNumber;
    }

    /// <summary>Per-chromosome <c>bafsegmented</c>/<c>logR</c> data of <c>merge_segments</c> (<c>findOverlaps</c> by position).</summary>
    private sealed class BattenbergMergeChromosome(
        IReadOnlyList<PhasedBafSegmentedSnp> rows,
        List<int> rowIndices,
        List<(long Pos, double LogR)> probes,
        double[] bafSegOut,
        BattenbergSegmentBafOption bafOption)
    {
        public List<double> PhasedBafs(long start, long end)
        {
            var values = new List<double>();
            foreach (int i in rowIndices)
            {
                if (rows[i].Position >= start && rows[i].Position <= end)
                {
                    values.Add(rows[i].BafPhased);
                }
            }

            return values;
        }

        /// <summary>The probes overlapping [start, end]: their count (NA included) and the non-NA values.</summary>
        public (int Count, List<double> Values) LogRs(long start, long end)
        {
            int count = 0;
            var values = new List<double>();
            foreach ((long pos, double value) in probes)
            {
                if (pos >= start && pos <= end)
                {
                    count++;
                    if (!double.IsNaN(value))
                    {
                        values.Add(value);
                    }
                }
            }

            return (count, values);
        }

        /// <summary><c>merge_seg</c> BAF (<c>calc_seg_baf_option</c>) and LogR of [start, end]; writes BAFseg.</summary>
        public (double Baf, double LogR) Recompute(long start, long end)
        {
            List<double> baf = PhasedBafs(start, end);
            double newBaf;
            switch (bafOption)
            {
                case BattenbergSegmentBafOption.Median:
                    newBaf = StatisticsHelper.Median(baf);
                    break;
                case BattenbergSegmentBafOption.Mean:
                    newBaf = StatisticsHelper.RMean(baf);
                    break;
                default:
                    double median = StatisticsHelper.Median(baf);
                    newBaf = median != 0.0 && median != 1.0 ? median : StatisticsHelper.RMean(baf);
                    break;
            }

            (int count, List<double> values) = LogRs(start, end);
            double newLogR = count == 0 ? 0.0 : StatisticsHelper.RMean(values); // mean(numeric(0)) = NaN when every probe is NA
            foreach (int i in rowIndices)
            {
                if (rows[i].Position >= start && rows[i].Position <= end)
                {
                    bafSegOut[i] = newBaf;
                }
            }

            return (newBaf, newLogR);
        }
    }

    /// <summary>The <c>merge_segments</c> loop over one chromosome (segments in place).</summary>
    private static void BattenbergMergeChromosomeSegments(
        List<BattenbergMergeSegment> segs, BattenbergMergeChromosome data, double rho, double cellCopies, double gamma)
    {
        segs[0].PrevChecked = true;
        segs[^1].NextChecked = true;
        while (true)
        {
            // Smallest segment (width = end − start + 1, first on a tie) among those with an unchecked side.
            int index = -1;
            long bestWidth = long.MaxValue;
            for (int k = 0; k < segs.Count; k++)
            {
                if ((!segs[k].PrevChecked || !segs[k].NextChecked) && segs[k].End - segs[k].Start + 1 < bestWidth)
                {
                    bestWidth = segs[k].End - segs[k].Start + 1;
                    index = k;
                }
            }

            if (index < 0)
            {
                return;
            }

            // Neighbours ordered by distance (order() is stable: the previous neighbour first on a tie).
            var neighbours = new List<int>(2);
            if (index > 0)
            {
                neighbours.Add(index - 1);
            }

            if (index < segs.Count - 1)
            {
                neighbours.Add(index + 1);
            }

            if (neighbours.Count == 2 && GRangesDistance(segs[index], segs[index + 1]) < GRangesDistance(segs[index], segs[index - 1]))
            {
                neighbours.Reverse();
            }

            foreach (int other in neighbours)
            {
                BattenbergMergeSegment a = segs[index], b = segs[other];
                bool checkedPair = other > index ? a.NextChecked && b.PrevChecked : a.PrevChecked && b.NextChecked;
                if (checkedPair)
                {
                    continue;
                }

                bool merge;
                if (GRangesDistance(a, b) > BattenbergMergeMaxDistance)
                {
                    merge = false;
                }
                else if (a.Major1 == b.Major1 && a.Minor1 == b.Minor1 && a.Frac1 == 1.0 && b.Frac1 == 1.0)
                {
                    merge = true; // same clonal copy-number solution
                }
                else
                {
                    merge = BattenbergSameSquareAndNotDifferent(a, b, data, rho, cellCopies, gamma);
                }

                if (merge)
                {
                    BattenbergMergeInto(segs, index, other, data);
                    break;
                }

                if (other > index)
                {
                    a.NextChecked = true;
                    b.PrevChecked = true;
                }
                else
                {
                    a.PrevChecked = true;
                    b.NextChecked = true;
                }
            }
        }
    }

    /// <summary>
    /// <c>merge_segments</c> for two different solutions: same "square" (rounded nmin or nmaj equal) and neither the
    /// logR nor the BAF Welch t-test significant (each side needs &gt; 10 non-NA values).
    /// </summary>
    private static bool BattenbergSameSquareAndNotDifferent(
        BattenbergMergeSegment a, BattenbergMergeSegment b, BattenbergMergeChromosome data, double rho, double cellCopies, double gamma)
    {
        // calc_nmin / calc_nmax in R's left-to-right operation order.
        double powA = Math.Pow(2.0, a.LogR / gamma), powB = Math.Pow(2.0, b.LogR / gamma);
        double nminA = RRound((rho - 1.0 - (a.Baf - 1.0) * powA * cellCopies) / rho);
        double nmajA = RRound((rho - 1.0 + a.Baf * powA * cellCopies) / rho);
        double nminB = RRound((rho - 1.0 - (b.Baf - 1.0) * powB * cellCopies) / rho);
        double nmajB = RRound((rho - 1.0 + b.Baf * powB * cellCopies) / rho);
        if (!(nminA == nminB || nmajA == nmajB))
        {
            return false; // different squares
        }

        List<double> logRA = data.LogRs(a.Start, a.End).Values, logRB = data.LogRs(b.Start, b.End).Values;
        List<double> bafA = data.PhasedBafs(a.Start, a.End), bafB = data.PhasedBafs(b.Start, b.End);
        if (!(logRA.Count > BattenbergMergeMinValues && logRB.Count > BattenbergMergeMinValues
              && bafA.Count > BattenbergMergeMinValues && bafB.Count > BattenbergMergeMinValues))
        {
            return false; // too few values
        }

        double pLogR = StatisticsHelper.WelchTTestPValue(logRA, logRB);
        double pBaf = StatisticsHelper.WelchTTestPValue(bafA, bafB);
        if (double.IsNaN(pLogR) || double.IsNaN(pBaf))
        {
            throw new InvalidOperationException(
                "Battenberg merge_segments: t.test data are essentially constant (R stops with an error).");
        }

        return !(pLogR < BattenbergMergeSignificance) && !(pBaf < BattenbergMergeSignificance);
    }

    /// <summary><c>merge_seg</c>: extend <paramref name="into"/> over <paramref name="index"/>, drop it, recompute, re-open.</summary>
    private static void BattenbergMergeInto(List<BattenbergMergeSegment> segs, int index, int into, BattenbergMergeChromosome data)
    {
        BattenbergMergeSegment target = segs[into];
        if (into < index)
        {
            target.End = segs[index].End;
        }
        else
        {
            target.Start = segs[index].Start;
        }

        segs.RemoveAt(index);
        int k = into < index ? index - 1 : index;
        // updateAround: the merged segment and its neighbours must be re-checked.
        if (k > 0)
        {
            segs[k].PrevChecked = false;
            segs[k - 1].NextChecked = false;
        }
        else
        {
            segs[k].PrevChecked = true;
        }

        if (k < segs.Count - 1)
        {
            segs[k].NextChecked = false;
            segs[k + 1].PrevChecked = false;
        }
        else
        {
            segs[k].NextChecked = true;
        }

        (target.Baf, target.LogR) = data.Recompute(target.Start, target.End);
        target.WasMerged = true;
    }

    /// <summary><c>GenomicRanges::distance</c> of two ranges on one chromosome: 0 if they overlap or are adjacent.</summary>
    private static long GRangesDistance(BattenbergMergeSegment a, BattenbergMergeSegment b) =>
        Math.Max(0L, Math.Max(b.Start - a.End - 1, a.Start - b.End - 1));

    /// <summary>R <c>as.numeric(as.character(x))</c> (15 significant digits), as <c>determine_copynumber</c>'s character matrix.</summary>
    private static double RAsCharacterRoundTrip(double x) =>
        double.IsFinite(x) ? double.Parse(x.ToString("G15", System.Globalization.CultureInfo.InvariantCulture), System.Globalization.CultureInfo.InvariantCulture) : x;

    /// <summary>
    /// Seqlevel order of <c>makeGRangesFromDataFrame</c> for character seqnames: <c>seqlevels[rankSeqlevels(seqlevels)] &lt;-
    /// seqlevels</c> (GenomeInfoDb 1.38 <c>rankSeqlevels</c>, C collation). Super-groups by prefix <c>CHR</c>, <c>chr</c>,
    /// <c>CH</c>, <c>ch</c>, none; inside each: Roman numerals, numbers (≤ 6 digits; suffix A/a/B/b/L/R), W, Z, X (sex
    /// chromosome unless it is the only kind of name), Y, U, M, MT, number+text, Wxxx, Zxxx, Xxxx, Yxxx, Uxxx, Mxxx, MTxxx,
    /// then everything else alphabetically. Ties (only for equal Roman values such as IV/IIII) keep ordinal order.
    /// </summary>
    internal static IReadOnlyList<string> BattenbergSeqlevelOrder(IReadOnlyList<string> names)
    {
        string[] levels = names.Distinct(StringComparer.Ordinal).OrderBy(n => n, StringComparer.Ordinal).ToArray();
        var prov = new int?[levels.Length];
        int last = 0;

        void Assign(List<int> idx, List<int>? ints = null)
        {
            if (idx.Count == 0)
            {
                return;
            }

            ints ??= Enumerable.Repeat(0, idx.Count).ToList();
            int min = ints.Min(), max = int.MinValue;
            for (int k = 0; k < idx.Count; k++)
            {
                int id = last + ints[k] - min + 1;
                prov[idx[k]] = id;
                max = Math.Max(max, id);
            }

            last = max;
        }

        foreach (string prefix in new[] { "CHR", "chr", "CH", "ch", string.Empty })
        {
            var sg = new List<int>();
            for (int i = 0; i < levels.Length; i++)
            {
                if (prov[i] is null && levels[i].StartsWith(prefix, StringComparison.Ordinal))
                {
                    sg.Add(i);
                }
            }

            if (sg.Count == 0)
            {
                continue;
            }

            string[] suf = sg.Select(i => levels[i][prefix.Length..]).ToArray();
            int n = suf.Length;
            var nbKind = new int[n]; // −1 none; 0 plain; 1 A; 2 a; 3 B; 4 b; 5 L; 6 R
            var nbxxx = new bool[n];
            for (int j = 0; j < n; j++)
            {
                nbKind[j] = ShortNumberKind(suf[j]);
                nbxxx[j] = nbKind[j] < 0 && LeadingDigits(suf[j]) is { Length: > 0 and <= 6 } d && d[0] != '0';
            }

            bool anyX = suf.Contains("X");
            bool xSexual = !anyX || suf.Contains("Y") || nbKind.Any(k => k >= 0) || nbxxx.Any(v => v);
            var roman = new int?[n];
            for (int j = 0; j < n; j++)
            {
                bool seXual = suf[j] == "X" && xSexual;
                roman[j] = seXual ? null : RomanValue(suf[j]);
            }

            bool Has(int j, string p) => suf[j].StartsWith(p, StringComparison.Ordinal);
            List<int> Where(Func<int, bool> f) => Enumerable.Range(0, n).Where(f).ToList();
            List<int> Map(List<int> js) => js.Select(j => sg[j]).ToList();
            List<int> FactorCodes(List<int> js, int drop) =>
                FactorInts(js.Select(j => suf[j][drop..]).ToList());

            var isRoman = Where(j => roman[j] is not null);
            Assign(Map(isRoman), isRoman.Select(j => roman[j]!.Value).ToList());
            var isNb = Where(j => nbKind[j] >= 0);
            Assign(Map(isNb), isNb.Select(j =>
            {
                string digits = nbKind[j] == 0 ? suf[j] : suf[j][..^1];
                return 7 * int.Parse(digits, System.Globalization.CultureInfo.InvariantCulture) + nbKind[j];
            }).ToList());
            Assign(Map(Where(j => suf[j] == "W")));
            Assign(Map(Where(j => suf[j] == "Z")));
            Assign(Map(Where(j => suf[j] == "X" && xSexual)));
            Assign(Map(Where(j => suf[j] == "Y")));
            Assign(Map(Where(j => suf[j] == "U")));
            Assign(Map(Where(j => suf[j] == "M")));
            Assign(Map(Where(j => suf[j] == "MT")));
            var isNbxxx = Where(j => nbxxx[j]);
            if (isNbxxx.Count > 0)
            {
                var ints1 = isNbxxx.Select(j => int.Parse(LeadingDigits(suf[j]), System.Globalization.CultureInfo.InvariantCulture)).ToList();
                var ints2 = FactorInts(isNbxxx.Select(j => suf[j][LeadingDigits(suf[j]).Length..]).ToList());
                int m2 = ints2.Max() + 1;
                Assign(Map(isNbxxx), ints1.Select((v, k) => m2 * v + ints2[k]).ToList());
            }

            bool Special(int j) => roman[j] is null && !(suf[j] == "X" && xSexual);
            var wx = Where(j => Has(j, "W") && suf[j] != "W" && roman[j] is null);
            var zx = Where(j => Has(j, "Z") && suf[j] != "Z" && roman[j] is null);
            var xx = Where(j => Has(j, "X") && suf[j] != "X" && roman[j] is null);
            var yx = Where(j => Has(j, "Y") && suf[j] != "Y" && roman[j] is null);
            var ux = Where(j => Has(j, "U") && suf[j] != "U" && roman[j] is null);
            var mtx = Where(j => Has(j, "MT") && suf[j] != "MT" && roman[j] is null);
            var mx = Where(j => Has(j, "M") && suf[j] != "M" && suf[j] != "MT" && !mtx.Contains(j) && roman[j] is null);
            Assign(Map(wx), FactorCodes(wx, 1));
            Assign(Map(zx), FactorCodes(zx, 1));
            Assign(Map(xx), FactorCodes(xx, 1));
            Assign(Map(yx), FactorCodes(yx, 1));
            Assign(Map(ux), FactorCodes(ux, 1));
            Assign(Map(mx), FactorCodes(mx, 1));
            Assign(Map(mtx), FactorCodes(mtx, 2));
            var named = new HashSet<string>(StringComparer.Ordinal) { "W", "Z", "Y", "U", "M", "MT" };
            var rest = Where(j => Special(j) && nbKind[j] < 0 && !nbxxx[j] && !named.Contains(suf[j]) && suf[j] != "X"
                                  && !wx.Contains(j) && !zx.Contains(j) && !xx.Contains(j) && !yx.Contains(j)
                                  && !ux.Contains(j) && !mx.Contains(j) && !mtx.Contains(j));
            Assign(Map(rest), FactorCodes(rest, 0));
        }

        return Enumerable.Range(0, levels.Length).OrderBy(i => prov[i]!.Value).ThenBy(i => i).Select(i => levels[i]).ToArray();

        static string LeadingDigits(string s)
        {
            int k = 0;
            while (k < s.Length && s[k] >= '0' && s[k] <= '9')
            {
                k++;
            }

            return s[..k];
        }

        // .isShortNb(x, abc): ^[1-9][0-9]*abc$ with at most 6 leading digits; abc ∈ {"", A, a, B, b, L, R}.
        static int ShortNumberKind(string s)
        {
            string d = LeadingDigits(s);
            if (d.Length == 0 || d.Length > 6 || d[0] == '0')
            {
                return -1;
            }

            if (d.Length == s.Length)
            {
                return 0;
            }

            return s.Length == d.Length + 1 ? "AaBbLR".IndexOf(s[^1], StringComparison.Ordinal) switch { -1 => -1, int i => i + 1 } : -1;
        }

        // as.integer(factor(x)) in C collation: 1-based rank among the distinct values.
        static List<int> FactorInts(List<string> values)
        {
            string[] distinct = values.Distinct(StringComparer.Ordinal).OrderBy(v => v, StringComparer.Ordinal).ToArray();
            return values.Select(v => Array.BinarySearch(distinct, v, StringComparer.Ordinal) + 1).ToList();
        }

        // GenomeInfoDb isRoman / utils:::.roman2numeric: upper-case valid Roman numeral → its value, else null.
        static int? RomanValue(string s)
        {
            if (s.Length == 0 || s != s.ToUpperInvariant())
            {
                return null;
            }

            string y = s.Replace("CM", "DCCCC", StringComparison.Ordinal).Replace("CD", "CCCC", StringComparison.Ordinal)
                .Replace("XC", "LXXXX", StringComparison.Ordinal).Replace("XL", "XXXX", StringComparison.Ordinal)
                .Replace("IX", "VIIII", StringComparison.Ordinal).Replace("IV", "IIII", StringComparison.Ordinal);
            if (!System.Text.RegularExpressions.Regex.IsMatch(y, "^M{0,3}D?C{0,4}L?X{0,4}V?I{0,4}$"))
            {
                return null;
            }

            int v = 0;
            foreach (char c in y)
            {
                v += c switch { 'M' => 1000, 'D' => 500, 'C' => 100, 'L' => 50, 'X' => 10, 'V' => 5, _ => 1 };
            }

            return v;
        }
    }
}
