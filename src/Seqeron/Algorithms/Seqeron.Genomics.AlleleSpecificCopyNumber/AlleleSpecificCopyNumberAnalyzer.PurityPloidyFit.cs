// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Allele-specific copy-number fit: purity / ploidy by the ASCAT distance minimisation, multiplicity (port of ASCAT 3.2.0, VanLoo-lab).</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{

    /// <summary>
    /// Platform/technology parameter γ in the ASCAT logR model. For massively parallel sequencing data
    /// (WGS/WES/TS) γ = 1; the SNP-array default 0.55 does not apply. Source: ASCAT README / Van Loo lab
    /// (VanLoo-lab/ascat): "For massively parallel sequencing data, gamma should always be set to 1."
    /// </summary>
    public const double AscatSequencingGamma = 1.0;

    /// <summary>
    /// Worst-case squared distance of a value to the nearest integer, (1/2)² = 0.25; the per-segment term in
    /// the ASCAT theoretical-maximum-distance used to normalise goodness of fit to a percentage. Source:
    /// ascat.runAscat.R — <c>TheoretMaxdist = sum(rep(0.25, n) * length * ...)</c>.
    /// </summary>
    private const double AscatWorstCaseIntegerDistance = 0.25;

    /// <summary>
    /// Down-weight applied to balanced (BAF = 0.5) segments in the ASCAT goodness-of-fit, because such
    /// segments carry little allele-specific information. Source: ascat.runAscat.R —
    /// <c>ifelse(b == 0.5, 0.05, 1)</c>.
    /// </summary>
    private const double AscatBalancedSegmentWeight = 0.05;

    /// <summary>BAF value of a perfectly balanced (1:1) heterozygous segment; the down-weight pivot in the GoF.</summary>
    private const double BalancedBaf = 0.5;

    /// <summary>
    /// A single per-locus allele-specific measurement at a germline-heterozygous SNP: the log-R ratio (total
    /// signal, "r") and the B-allele frequency (allelic contrast, "b"). These are the two ASCAT input tracks
    /// (Van Loo et al. 2010, PNAS) and are <b>observed measurements</b> supplied by the caller — they are the
    /// raw data, not a derived quantity.
    /// </summary>
    /// <param name="Chromosome">Contig label (used to group loci into per-chromosome segments).</param>
    /// <param name="Position">0-based genomic coordinate of the SNP.</param>
    /// <param name="LogR">Log-R ratio r (log2 total-signal ratio vs the reference baseline).</param>
    /// <param name="BAF">B-allele frequency b ∈ [0, 1] at the germline-heterozygous SNP.</param>
    public readonly record struct AlleleSpecificLocus(string Chromosome, long Position, double LogR, double BAF);

    /// <summary>
    /// A genomic segment summarised by its mean logR and mean BAF, produced by allele-specific segmentation
    /// of per-locus <see cref="AlleleSpecificLocus"/> data. A single fitted logR value and a BAF value are
    /// obtained per segment (Van Loo et al. 2010, ASCAT).
    /// </summary>
    /// <param name="Chromosome">Contig label.</param>
    /// <param name="Start">0-based start coordinate (first locus position in the segment).</param>
    /// <param name="End">End coordinate (last locus position in the segment; End ≥ Start).</param>
    /// <param name="MeanLogR">Length-unweighted mean logR over the segment's loci.</param>
    /// <param name="MeanBAF">Mean "folded" BAF (distance from 0.5, re-centred) over the segment's loci.</param>
    /// <param name="LocusCount">Number of loci summarised by the segment.</param>
    public readonly record struct AlleleSpecificSegmentSummary(
        string Chromosome,
        long Start,
        long End,
        double MeanLogR,
        double MeanBAF,
        int LocusCount)
    {
        /// <summary>Segment length in base pairs (End − Start).</summary>
        public long Length => End - Start;
    }

    /// <summary>
    /// Result of the joint ASCAT purity/ploidy fit: the recovered purity ρ, the ASCAT output ploidy, the goodness of
    /// fit, the allele-specific integer copy-number segments, and the model ploidy parameter ψ.
    /// </summary>
    /// <param name="Purity">Recovered tumour purity ρ (aberrant cell fraction) ∈ (0, 1] (ASCAT <c>purity</c>).</param>
    /// <param name="Ploidy">ASCAT output <c>ploidy</c>: the mean integer total copy number (major + minor) over the
    /// heterozygous probes, i.e. the <see cref="AlleleSpecificSegmentSummary.LocusCount"/>-weighted mean of the emitted
    /// segments' total copy number. It equals ψ for an exact integer-copy-number genome. From the germline-aware path
    /// (<see cref="FitPurityPloidyFromAspcf"/>, B24 F37) it is ASCAT's <c>mean(nA + nB)</c> over <b>all</b> probes,
    /// heterozygous and homozygous.</param>
    /// <param name="GoodnessOfFit">Percentage goodness of fit (1 − distance/TheoretMaxdist)·100, in (−∞, 100].</param>
    /// <param name="Segments">The allele-specific integer copy-number segments (major/minor CN) implied by (ρ, ψ).</param>
    public readonly record struct PurityPloidyFit(
        double Purity,
        double Ploidy,
        double GoodnessOfFit,
        IReadOnlyList<AlleleSpecificSegment> Segments)
    {
        /// <summary>The selected ploidy parameter ψ of the ASCAT model (ASCAT <c>psi</c>, a grid value).</summary>
        public double Psi { get; init; }

        /// <summary>ASCAT <c>nonaberrant</c> flag: ≤ 3 % of the probes are allelically imbalanced and no imbalanced
        /// segment exceeds 0.5 % of the probes (MINABB / MINABBREGION).</summary>
        public bool IsNonAberrant { get; init; }
    }

    /// <summary>ASCAT sample sex (<c>gender</c> argument of <c>ascat.loadData</c> / <c>runASCAT</c>).</summary>
    public enum AscatGender
    {
        /// <summary><c>"XX"</c> (ASCAT default when <c>gender = NULL</c>): X and Y are emitted with the diploid model.</summary>
        XX,

        /// <summary><c>"XY"</c>: X (non-PAR) and Y are haploid — <c>nA = (ρ − 1 + (2(1 − ρ) + ρψ)·2^(r/γ))/ρ</c>, <c>nB = 0</c>.</summary>
        XY,
    }

    /// <summary>
    /// Sex-chromosome model of the ASCAT copy-number output (runASCAT, ascat.runAscat.R, VanLoo-lab/ascat). It does
    /// <b>not</b> change the purity/ploidy fit — ASCAT always fits on autosomes only
    /// (<c>autoprobes = !(SNPposhet[,1] %in% sexchromosomes)</c>) — only the integer segments emitted on X and Y and
    /// hence the reported ploidy:
    /// <list type="bullet">
    /// <item><see cref="AscatGender.XX"/> (<see cref="Female"/>, the default): <c>haploidchrs</c> is empty, so X and Y
    /// use the diploid equations (Y is not special-cased: the <c>nullprobes</c> branch is reached only for non-diploid
    /// probes). <see cref="XNonPar"/> is ignored, as in ASCAT (<c>!is.null(X_nonPAR) &amp;&amp; gender == "XY"</c>).</item>
    /// <item><see cref="AscatGender.XY"/>: <c>haploidchrs = c("X", "Y")</c>; a haploid segment gets
    /// <c>nAraw = (rho − 1 + ((1 − rho)·2 + rho·psi)·2^(logR/gamma))/rho</c>, <c>nBraw = 0</c> (normal cells carry one
    /// copy), followed by the usual negative-value correction and R rounding. Without <see cref="XNonPar"/>
    /// (ASCAT <c>genomeVersion = NULL</c>) the whole of X is haploid. With it, ASCAT's
    /// <c>diploidprobes_fixnonPAR</c> applies: an X segment is haploid only when its overlap with the closed interval
    /// [Start, End] exceeds 50 % of the segment's closed span [first probe position, last probe position]
    /// (IRanges widths), otherwise it is diploid (pseudo-autosomal). Y is always haploid.</item>
    /// </list>
    /// Chromosome labels are matched case-insensitively with an optional "chr" prefix (as for the autosomal filter).
    /// The interval is compared with the segment coordinates exactly as given (ASCAT compares it with <c>SNPpos</c>,
    /// which is 1-based; the presets are ASCAT's 1-based constants).
    /// </summary>
    public sealed record AscatSexModel
    {
        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "hg19"</c> (GRCh37): <c>c(2699521, 154931043)</c>.</summary>
        public static (long Start, long End) XNonParHg19 { get; } = (2_699_521, 154_931_043);

        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "hg38"</c> (GRCh38): <c>c(2781480, 155701382)</c>.</summary>
        public static (long Start, long End) XNonParHg38 { get; } = (2_781_480, 155_701_382);

        /// <summary>ASCAT <c>X_nonPAR</c> for <c>genomeVersion = "CHM13"</c> (T2T-CHM13): <c>c(2394411, 153925834)</c>.</summary>
        public static (long Start, long End) XNonParChm13 { get; } = (2_394_411, 153_925_834);

        /// <summary>ASCAT default (<c>gender = "XX"</c>): X and Y diploid. Bit-identical to the gender-less entry points.</summary>
        public static AscatSexModel Female { get; } = new(AscatGender.XX, null);

        /// <summary><c>gender = "XY"</c>, <c>X_nonPAR = NULL</c> (ASCAT <c>genomeVersion = NULL</c>): all of X and Y haploid.</summary>
        public static AscatSexModel Male { get; } = new(AscatGender.XY, null);

        /// <summary>The sample sex.</summary>
        public AscatGender Gender { get; }

        /// <summary>Closed X non-PAR interval (ASCAT <c>X_nonPAR</c>); <c>null</c> = whole X haploid in a male.</summary>
        public (long Start, long End)? XNonPar { get; }

        /// <summary>Creates a sex model.</summary>
        /// <param name="gender">The sample sex.</param>
        /// <param name="xNonPar">Optional closed X non-PAR interval (Start ≤ End); used only for <see cref="AscatGender.XY"/>.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="gender"/> is undefined or the interval has Start &gt; End.</exception>
        public AscatSexModel(AscatGender gender, (long Start, long End)? xNonPar = null)
        {
            if (gender is not (AscatGender.XX or AscatGender.XY))
            {
                throw new ArgumentOutOfRangeException(nameof(gender), gender, "Gender must be XX or XY.");
            }

            if (xNonPar is { } interval && interval.Start > interval.End)
            {
                throw new ArgumentOutOfRangeException(nameof(xNonPar), xNonPar, "The X non-PAR interval needs Start ≤ End.");
            }

            Gender = gender;
            XNonPar = xNonPar;
        }

        /// <summary>A male (<c>"XY"</c>) model with ASCAT's X non-PAR interval of <paramref name="genome"/>
        /// (<c>ascat.loadData(genomeVersion = "hg19" | "hg38")</c>).</summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="genome"/> is undefined.</exception>
        public static AscatSexModel MaleWithXNonPar(ReferenceGenome genome) => genome switch
        {
            ReferenceGenome.GRCh37 => new AscatSexModel(AscatGender.XY, XNonParHg19),
            ReferenceGenome.GRCh38 => new AscatSexModel(AscatGender.XY, XNonParHg38),
            _ => throw new ArgumentOutOfRangeException(nameof(genome), genome, "Unknown reference genome."),
        };

        /// <summary>
        /// True when a segment is haploid under this model (runASCAT <c>!diploidprobes</c> at the segment start; ASCAT's
        /// <c>nullprobes</c> branch is unreachable with the default <c>sexchromosomes = c("X", "Y")</c>).
        /// </summary>
        internal bool IsHaploid(string chromosome, long start, long end)
        {
            if (Gender != AscatGender.XY)
            {
                return false;
            }

            ReadOnlySpan<char> name = StripChrPrefix(chromosome.AsSpan().Trim());
            if (name.Equals("Y", StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }

            if (!name.Equals("X", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }

            if (XNonPar is not { } nonPar)
            {
                return true;
            }

            // diploidprobes_fixnonPAR: findOverlaps(nonPAR, SEGMENTS) and width(pintersect)/width(segment) > 0.5.
            double overlap = (double)Math.Min(end, nonPar.End) - Math.Max(start, nonPar.Start) + 1.0;
            return overlap > 0.0 && overlap / ((double)end - start + 1.0) > 0.5;
        }
    }

    /// <summary>
    /// Male X non-PAR germline re-genotyping of <c>ascat.aspcf</c> / <c>ascat.asmultipcf</c> (VanLoo-lab/ascat
    /// ascat.aspcf.R, ascat.asmultipcf.R, B24 F58). It applies only when <see cref="SexModel"/> is
    /// <see cref="AscatGender.XY"/> with an <see cref="AscatSexModel.XNonPar"/> interval
    /// (<c>!is.null(X_nonPAR) &amp;&amp; gender == "XY"</c>); otherwise the caller's genotypes are used unchanged. ASCAT then
    /// takes the X probes (label "X", optional "chr" prefix, case-insensitive as in <see cref="AscatSexModel"/>) whose
    /// position lies in the closed non-PAR interval and whose genotype is known; if there are more than 5, it marks all of
    /// them homozygous and re-marks k = <c>round(m · h_auto)</c> of the m probes heterozygous, h_auto = the heterozygous
    /// fraction of the autosomal probes (every label other than X and Y; R <c>table(gg[autosomes])</c>), R half-even
    /// rounding:
    /// <list type="bullet">
    /// <item>with <see cref="GermlineBaf"/> (ASCAT <c>Germline_BAF</c>): the k probes of smallest
    /// <c>DIST = 1 − max(b, 1 − b)</c> (distance of the germline BAF to 0/1), ties broken at random —
    /// <c>rank(DIST, ties.method = "random") ≤ k</c>, which R computes as <c>sort.list(order(DIST, runif(m)))</c>;</item>
    /// <item>without it (ASCAT <c>Germline_BAF = NULL</c>, genotypes from <c>ascat.gg</c>): <c>sample(nonPAR_index, k)</c>.</item>
    /// </list>
    /// Draws come from R's default generator seeded by <c>set.seed(<see cref="Seed"/>)</c> (Mersenne-Twister, R ≥ 3.6
    /// "Rejection" sampling — the same port as the Battenberg bootstrap). ASCAT's default seed is
    /// <c>as.integer(Sys.time())</c> (the Unix time in whole seconds), so an R run is reproducible only with an explicit
    /// seed; a null seed here mirrors that default. One call consumes the stream of one ASCAT call on one male sample;
    /// <c>ascat.aspcf</c> over several samples (one <c>set.seed</c>, the stream continuing from male sample to male
    /// sample) is <see cref="SegmentAlleleSpecificAspcfSamples"/> (B24 F65).
    /// </summary>
    public sealed record AscatMaleXGenotyping
    {
        /// <summary>Creates the re-genotyping settings.</summary>
        /// <param name="sexModel">Sex model; re-genotyping runs only for <see cref="AscatGender.XY"/> with an X non-PAR interval.</param>
        /// <param name="seed">R <c>set.seed</c> value (ASCAT <c>seed</c>); null = <see cref="DefaultSeed"/> (ASCAT's
        /// <c>as.integer(Sys.time())</c>).</param>
        /// <param name="germlineBaf">Optional germline BAF per locus (ASCAT <c>Germline_BAF</c>, first column), each in
        /// [0, 1] or NaN (NaN = genotype unknown, as ASCAT's <c>gg = Germline_BAF &lt; 0.3 | Germline_BAF &gt; 0.7</c> is NA
        /// there: the locus is left out of the non-PAR selection and of the autosomal fraction). Null selects ASCAT's
        /// <c>sample()</c> branch.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sexModel"/> is null.</exception>
        /// <exception cref="ArgumentException">A germline BAF is outside [0, 1] (and not NaN).</exception>
        public AscatMaleXGenotyping(AscatSexModel sexModel, int? seed = null, IReadOnlyList<double>? germlineBaf = null)
        {
            ArgumentNullException.ThrowIfNull(sexModel);
            if (germlineBaf is not null)
            {
                foreach (double b in germlineBaf)
                {
                    if (!double.IsNaN(b) && b is not (>= 0.0 and <= 1.0))
                    {
                        throw new ArgumentException("Every germline BAF must be in [0, 1] or NaN.", nameof(germlineBaf));
                    }
                }
            }

            SexModel = sexModel;
            Seed = seed ?? DefaultSeed();
            GermlineBaf = germlineBaf;
        }

        /// <summary>The sex model (gender and X non-PAR interval).</summary>
        public AscatSexModel SexModel { get; }

        /// <summary>The R <c>set.seed</c> value.</summary>
        public int Seed { get; }

        /// <summary>The germline BAF per locus, or null (ASCAT's <c>sample()</c> branch).</summary>
        public IReadOnlyList<double>? GermlineBaf { get; }

        /// <summary>ASCAT's default <c>seed = as.integer(Sys.time())</c>: the current Unix time in whole seconds.</summary>
        public static int DefaultSeed() => unchecked((int)DateTimeOffset.UtcNow.ToUnixTimeSeconds());

        /// <summary>True when ASCAT re-genotypes (<c>!is.null(X_nonPAR) &amp;&amp; gender == "XY"</c>).</summary>
        internal bool IsActive => SexModel.Gender == AscatGender.XY && SexModel.XNonPar is not null;
    }

    /// <summary>
    /// Applies <see cref="AscatMaleXGenotyping"/> to <paramref name="het"/> in place (ascat.aspcf.R ll. 51–70,
    /// ascat.asmultipcf.R ll. 44–63, verbatim order of random draws). No-op when <paramref name="spec"/> is null or inactive.
    /// Draws from <paramref name="rng"/> when given (a stream shared across samples), else from a fresh
    /// <c>set.seed(spec.Seed)</c> stream.
    /// </summary>
    private static void RegenotypeMaleXNonPar(
        AlleleSpecificLocus[] loci, bool[] het, AscatMaleXGenotyping? spec, string paramName, RMersenneTwister? rng = null)
    {
        if (spec is null || !spec.IsActive)
        {
            return;
        }

        IReadOnlyList<double>? gbaf = spec.GermlineBaf;
        if (gbaf is not null && gbaf.Count != loci.Length)
        {
            throw new ArgumentException("Exactly one germline BAF is required per locus.", paramName);
        }

        (long start, long end) = spec.SexModel.XNonPar!.Value;
        var nonPar = new List<int>();
        int autosomal = 0;
        int autosomalHet = 0;
        for (int i = 0; i < loci.Length; i++)
        {
            bool known = gbaf is null || !double.IsNaN(gbaf[i]);
            ReadOnlySpan<char> name = StripChrPrefix(loci[i].Chromosome.AsSpan().Trim());
            bool isX = name.Equals("X", StringComparison.OrdinalIgnoreCase);
            if (isX)
            {
                if (known && loci[i].Position >= start && loci[i].Position <= end)
                {
                    nonPar.Add(i);
                }
            }
            else if (known && !name.Equals("Y", StringComparison.OrdinalIgnoreCase))
            {
                autosomal++;
                if (het[i])
                {
                    autosomalHet++;
                }
            }
        }

        int m = nonPar.Count;
        if (m <= 5)
        {
            return;
        }

        foreach (int i in nonPar)
        {
            het[i] = false;
        }

        rng ??= new RMersenneTwister(spec.Seed);
        if (gbaf is not null)
        {
            // DIST = 1 - max(x, 1 - x); rank(DIST, ties.method = "random") = sort.list(order(DIST, runif(m))).
            var dist = new double[m];
            var tieBreak = new double[m];
            for (int t = 0; t < m; t++)
            {
                double x = gbaf[nonPar[t]];
                dist[t] = 1.0 - MirrorBaf(x);
            }

            for (int t = 0; t < m; t++)
            {
                tieBreak[t] = rng.UnifRand();
            }

            if (autosomalHet == 0)
            {
                return; // autosomes_info["FALSE"] is NA: rank(...) <= NA selects nothing.
            }

            int k = (int)Math.Round(m * ((double)autosomalHet / autosomal), MidpointRounding.ToEven);
            int[] order = Enumerable.Range(0, m)
                .OrderBy(t => dist[t]).ThenBy(t => tieBreak[t]).ThenBy(t => t)
                .ToArray();
            for (int r = 0; r < k; r++)
            {
                het[nonPar[order[r]]] = true;
            }
        }
        else
        {
            if (autosomalHet == 0)
            {
                throw new ArgumentException(
                    "No heterozygous autosomal locus: ASCAT's sample(nonPAR_index, NA) fails (invalid 'size' argument).", paramName);
            }

            int k = (int)Math.Round(m * ((double)autosomalHet / autosomal), MidpointRounding.ToEven);
            foreach (int t in rng.SampleWithoutReplacement(m, k))
            {
                het[nonPar[t]] = true;
            }
        }
    }

    /// <summary>
    /// Segments per-locus allele-specific signal (logR, BAF) into contiguous regions, producing one
    /// (mean logR, mirrored BAF) summary per segment, by running the published ASCAT allele-specific segmentation
    /// (<c>ascat.aspcf</c>, VanLoo-lab/ascat ascat.aspcf.R; Nilsen et al. 2012, <i>BMC Genomics</i> 13:591; Ross et al.
    /// 2021, <i>Bioinformatics</i> 37:1909) with ASCAT's default penalty
    /// (<see cref="AspcfDefaultPenalty"/> = 70). It is exactly
    /// <c><see cref="SegmentAlleleSpecificAspcf"/>(loci, <see cref="AspcfDefaultPenalty"/>)</c>.
    /// <para>
    /// Compatibility (B24 F35): this name formerly ran an unsourced greedy left-to-right mean-shift heuristic with
    /// caller-chosen absolute thresholds. That heuristic has no published reference and was replaced; the threshold
    /// parameters below have no ASPCF meaning (ASPCF's cost is MAD-standardised and penalised, with the fixed minimum
    /// segment length kmin = 6) and are <b>ignored</b> — they are kept, and still range-checked, only so existing call
    /// sites keep compiling and keep their argument-validation behaviour. Pass a penalty to
    /// <see cref="SegmentAlleleSpecificAspcf"/> to tune the segmentation.
    /// </para>
    /// </summary>
    /// <param name="loci">Per-locus measurements; processed in input order within each chromosome. LogR must be
    /// finite and BAF in [0, 1] (as <see cref="SegmentAlleleSpecificAspcf"/>).</param>
    /// <param name="logRChangeThreshold">Ignored (former greedy-heuristic logR threshold; no ASPCF equivalent). Must be &gt; 0.</param>
    /// <param name="bafChangeThreshold">Ignored (former greedy-heuristic BAF threshold; no ASPCF equivalent). Must be &gt; 0.</param>
    /// <param name="minLociPerSegment">Ignored (ASPCF uses ASCAT's fixed kmin = 6). Must be ≥ 1.</param>
    /// <returns>The ASPCF segment summaries (mean raw logR, ASPCF mirrored BAF) in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException">a locus has a null chromosome label, a non-finite logR or a BAF outside [0, 1].</exception>
    /// <exception cref="ArgumentOutOfRangeException">a threshold ≤ 0 or NaN, or minLociPerSegment &lt; 1.</exception>
    public static IReadOnlyList<AlleleSpecificSegmentSummary> SegmentAlleleSpecific(
        IEnumerable<AlleleSpecificLocus> loci,
        double logRChangeThreshold,
        double bafChangeThreshold = 0.1,
        int minLociPerSegment = 1)
    {
        ArgumentNullException.ThrowIfNull(loci);

        if (double.IsNaN(logRChangeThreshold) || logRChangeThreshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(logRChangeThreshold), logRChangeThreshold, "The logR change threshold must be positive.");
        }

        if (double.IsNaN(bafChangeThreshold) || bafChangeThreshold <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(bafChangeThreshold), bafChangeThreshold, "The BAF change threshold must be positive.");
        }

        if (minLociPerSegment < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minLociPerSegment), minLociPerSegment, "At least one locus per segment is required.");
        }

        return SegmentAlleleSpecificAspcf(loci, AspcfDefaultPenalty);
    }

    /// <summary>
    /// Raw (real-valued) ASCAT allele-specific copy numbers (nA, nB) for one segment given (r, b, ρ, ψ, γ).
    /// Source: ascat.runAscat.R (VanLoo-lab/ascat), verbatim:
    /// <code>
    /// nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
    /// nB = (rho-1 +  b   *2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
    /// </code>
    /// With ASCAT's segmented BAF convention b ≤ 0.5 (<c>Tumor_BAF_segmented = 1 − bafPCFed</c>) nA is the major
    /// and nB the minor allele.
    /// </summary>
    private static (double NA, double NB) AscatRawCopyNumbers(double r, double b, double rho, double psi, double gamma)
    {
        double scaledTotal = Math.Pow(2.0, r / gamma) * MixtureCopiesPerCell(rho, psi);
        double nA = (rho - 1.0 - (b - 1.0) * scaledTotal) / rho;
        double nB = (rho - 1.0 + b * scaledTotal) / rho;
        return (nA, nB);
    }

    /// <summary>
    /// Raw ASCAT copy number of a haploid (male X non-PAR / Y) segment, ascat.runAscat.R <c>seg_raw</c>, verbatim:
    /// <c>nAraw = (rho-1 + ((1-rho)*2+rho*psi)*2^(logR/gamma))/rho</c> (the normal cells contribute one copy).
    /// </summary>
    private static double AscatRawHaploidCopyNumber(double r, double rho, double psi, double gamma) =>
        (rho - 1.0 + Math.Pow(2.0, r / gamma) * MixtureCopiesPerCell(rho, psi)) / rho;

    /// <summary>
    /// Average number of copies of a locus per cell in a tumour sample of purity ρ whose tumour cells carry
    /// <paramref name="tumorCopies"/> copies and whose normal cells are diploid: <c>2(1 − ρ) + ρ·n</c>. This mixture
    /// denominator is shared by the ASCAT logR model (<c>(1 − rho)·2 + rho·psi</c>, ascat.runAscat.R), Battenberg
    /// (<c>psi = rho·psit + 2(1 − rho)</c>) and the Landau/CNAqc expected-VAF model (<c>f = ρ·M·c / (2(1 − ρ) + ρ·q)</c>,
    /// <see cref="ClassifyClonality"/>). IEEE addition and multiplication are commutative, so every caller is
    /// bit-identical to its former inline form.
    /// </summary>
    public static double MixtureCopiesPerCell(double purity, double tumorCopies) =>
        NormalDiploidCopyNumber * (1.0 - purity) + purity * tumorCopies;

    /// <summary>
    /// CNAqc <c>expected_vaf_fun(m, purity, ploidy)</c> = <c>m·ρ / (2(1 − ρ) + ρ·n)</c>: the clonal VAF peak of a mutation
    /// on <paramref name="multiplicity"/> of the <paramref name="tumorCopies"/> tumour copies (shared by
    /// <see cref="AnalyzePurityPeaks"/> and <see cref="AnalyzeComplexKaryotypePeaks"/>; R operation order).
    /// </summary>
    public static double CnaqcExpectedVaf(int multiplicity, double purity, double tumorCopies) =>
        multiplicity * purity / MixtureCopiesPerCell(purity, tumorCopies);

    // ---- ASCAT runASCAT solution-selection constants (ascat.runAscat.R, verbatim) ----

    /// <summary>ASCAT <c>MINRHO = 0.2</c>: minimum aberrant-cell fraction of an accepted optimum.</summary>
    private const double AscatMinRho = 0.2;

    /// <summary>ASCAT <c>MINGOODNESSOFFIT = 80</c> (%): minimum goodness of fit of an accepted optimum.</summary>
    private const double AscatMinGoodnessOfFit = 80.0;

    /// <summary>ASCAT <c>MINPERCZERO = 0.02</c>: minimum fraction of allele copies rounded to 0 (pass 1).</summary>
    private const double AscatMinPercentZero = 0.02;

    /// <summary>ASCAT <c>MINPERCZEROABB = 0.1</c>: minimum zero-allele fraction over aberrant segments (passes 2–3).</summary>
    private const double AscatMinPercentZeroAberrant = 0.1;

    /// <summary>ASCAT <c>MINPERCODDEVEN = 0.05</c>: minimum fraction of odd/even allele pairs (pass 3).</summary>
    private const double AscatMinPercentOddEven = 0.05;

    /// <summary>ASCAT <c>MINPLOIDYSTRICT = 1.7</c>: strict lower ploidy bound of the fallback passes 2 and 4.</summary>
    private const double AscatMinPloidyStrict = 1.7;

    /// <summary>ASCAT <c>MAXPLOIDYSTRICT = 2.3</c>: strict upper ploidy bound of the fallback passes 2 and 4.</summary>
    private const double AscatMaxPloidyStrict = 2.3;

    /// <summary>ASCAT <c>MINABB = 0.03</c>: a sample with ≤ 3 % aberrant (BAF ≠ 0.5) probes may be non-aberrant.</summary>
    private const double AscatMinAberrantFraction = 0.03;

    /// <summary>ASCAT <c>MINABBREGION = 0.005</c>: … and no aberrant segment larger than 0.5 % of the probes.</summary>
    private const double AscatMinAberrantRegionFraction = 0.005;

    /// <summary>ASCAT local-minimum window half-width: <c>seld = d[(i-3):(i+3), (j-3):(j+3)]</c> (a 7 × 7 window).</summary>
    private const int AscatLocalMinimumHalfWindow = 3;

    /// <summary>ASCAT pass-3 mask value for grid columns with ρ &gt; 1: <c>d[, cold] = 1E20</c>.</summary>
    private const double AscatMaskedDistance = 1e20;

    /// <summary>ASCAT <c>limitround = 0.5</c>: odd-total evidence threshold for balanced (BAF = 0.5) segments.</summary>
    private const double AscatLimitRound = 0.5;

    /// <summary>
    /// One autosomal ASCAT fitting segment (<c>make_segments</c> row): segmented logR r, segmented BAF b in ASCAT's
    /// ≤ 0.5 orientation, and <c>length</c> = the number of germline-heterozygous probes of the segment.
    /// </summary>
    private readonly record struct AscatFitSegment(double R, double B, double Length);

    /// <summary>
    /// True when <paramref name="chromosome"/> is one of ASCAT's default <c>sexchromosomes = c("X", "Y")</c>
    /// (an optional "chr" prefix is ignored, case-insensitive). ASCAT excludes these probes from the purity/ploidy
    /// fit: <c>autoprobes = !(SNPposhet[,1] %in% sexchromosomes)</c>.
    /// </summary>
    private static bool IsAscatSexChromosome(string chromosome)
    {
        ReadOnlySpan<char> name = StripChrPrefix(chromosome.AsSpan().Trim());
        return name.Equals("X", StringComparison.OrdinalIgnoreCase) || name.Equals("Y", StringComparison.OrdinalIgnoreCase);
    }

    /// <summary>
    /// Input contract of the ASCAT / ASPCF / Battenberg ports for one (logR, BAF) pair — per locus or per segment
    /// summary: a finite logR and a BAF in [0, 1].
    /// </summary>
    private static bool IsValidAlleleSignal(double logR, double baf) =>
        double.IsFinite(logR) && !double.IsNaN(baf) && baf >= 0.0 && baf <= 1.0;

    /// <summary>ascat.aspcf / ascat.asmultipcf input contract for one BAF: in [0, 1] or NaN = R <c>NA</c> (B24 F59).</summary>
    private static bool IsAscatBafOrMissing(double baf) => double.IsNaN(baf) || (baf >= 0.0 && baf <= 1.0);

    /// <summary>Mirrored BAF max(b, 1 − b) ∈ [0.5, 1] (ascat.aspcf <c>ifelse(b &gt; 0.5, b, 1 − b)</c>; Battenberg <c>l</c>).</summary>
    private static double MirrorBaf(double baf) => OrientBaf(baf, baf);

    /// <summary>
    /// <paramref name="value"/> when <paramref name="reference"/> &gt; 0.5, otherwise 1 − value (R
    /// <c>ifelse(reference &gt; 0.5, value, 1 − value)</c>): mirrors a BAF onto the upper half
    /// (<see cref="MirrorBaf"/> = <c>OrientBaf(b, b)</c>) or restores the orientation of a mirrored value from the original
    /// BAF (ascat.aspcf / ascat.asmultipcf <c>bafwins</c>, Battenberg <c>BAFphased</c>). NaN reference ⇒ 1 − value.
    /// </summary>
    private static double OrientBaf(double reference, double value) => reference > BalancedBaf ? value : 1.0 - value;

    /// <summary>ASCAT segmented-BAF orientation: <c>Tumor_BAF_segmented = 1 − mirroredBAF</c> ∈ [0, 0.5].</summary>
    private static double ToAscatBaf(double baf) => baf > BalancedBaf ? 1.0 - baf : baf;

    /// <summary>R <c>.Machine$double.eps</c> = 2⁻⁵² (used by <c>seq.default</c>).</summary>
    private const double RMachineEpsilon = 2.220446049250313e-16;

    /// <summary>R <c>round()</c> (IEC 60559 round-half-to-even), as used throughout ascat.runAscat.R.</summary>
    public static double RRound(double x) => Math.Round(x, MidpointRounding.ToEven);

    /// <summary>R floored modulus <c>x %% 2</c> for an integral double (0 or 1, also for negative x).</summary>
    private static double RMod2(double x)
    {
        double m = x % 2.0;
        return m < 0.0 ? m + 2.0 : m;
    }

    /// <summary>
    /// R <c>seq(from, to, by)</c> for by &gt; 0 (<c>seq.default</c>): <c>n = as.integer((to−from)/by + 1e−10)</c>,
    /// <c>x = from + (0:n)·by</c>, clamped by <c>pmin(x, to)</c>; <c>from</c> alone when the span is negligible.
    /// </summary>
    private static double[] RSeq(double from, double to, double by)
    {
        double del = to - from;
        double scale = Math.Max(Math.Abs(to), Math.Abs(from));
        if (del == 0.0 || (scale > 0.0 && Math.Abs(del) / scale < 100.0 * RMachineEpsilon))
        {
            return new[] { from };
        }

        int n = (int)(del / by + 1e-10);
        var values = new double[n + 1];
        for (int k = 0; k <= n; k++)
        {
            values[k] = Math.Min(from + k * by, to);
        }

        return values;
    }

    /// <summary>
    /// The value R obtains from <c>as.numeric(rownames(d)[i])</c>: the grid value printed with 15 significant digits
    /// and parsed back (ASCAT reads the candidate ψ/ρ from the distance-matrix dimnames).
    /// </summary>
    private static double RDimnameValue(double value) =>
        double.Parse(value.ToString("G15", System.Globalization.CultureInfo.InvariantCulture),
            System.Globalization.CultureInfo.InvariantCulture);

    /// <summary>
    /// ASCAT <c>create_distance_matrix</c> cell: nA/nB at (ρ, ψ), the minor allele chosen genome-wide
    /// (<c>if (sum(nA) &lt; sum(nB)) nMinor = nA else nMinor = nB</c>), and
    /// <c>d = sum(abs(nMinor − pmax(round(nMinor), 0))^2 · length · ifelse(b == 0.5, 0.05, 1), na.rm = TRUE)</c>.
    /// </summary>
    private static double AscatDistance(AscatFitSegment[] segs, double rho, double psi, double gamma)
    {
        double sumA = 0.0, sumB = 0.0;
        var nA = new double[segs.Length];
        var nB = new double[segs.Length];
        for (int i = 0; i < segs.Length; i++)
        {
            (nA[i], nB[i]) = AscatRawCopyNumbers(segs[i].R, segs[i].B, rho, psi, gamma);
            if (!double.IsNaN(nA[i])) sumA += nA[i];
            if (!double.IsNaN(nB[i])) sumB += nB[i];
        }

        double[] nMinor = sumA < sumB ? nA : nB;
        double d = 0.0;
        for (int i = 0; i < segs.Length; i++)
        {
            double dev = Math.Abs(nMinor[i] - Math.Max(RRound(nMinor[i]), 0.0));
            double term = dev * dev * segs[i].Length * (segs[i].B == BalancedBaf ? AscatBalancedSegmentWeight : 1.0);
            if (!double.IsNaN(term))
            {
                d += term; // na.rm = TRUE
            }
        }

        return d;
    }

    /// <summary>
    /// Per-candidate statistics of one ASCAT local minimum (ascat.runAscat.R): the recomputed ploidy
    /// <c>sum((nA+nB)·length)/sum(length)</c>, <c>percentzero</c>, <c>perczeroAbb</c> (NaN → 0) and <c>percOddEven</c>.
    /// </summary>
    private static (double Ploidy, double PercentZero, double PercentZeroAberrant, double PercentOddEven) AscatCandidateStatistics(
        AscatFitSegment[] segs, double rho, double psi, double gamma)
    {
        double totalLength = 0.0, ploidyNumerator = 0.0, zero = 0.0, zeroAbb = 0.0, abbLength = 0.0, oddEven = 0.0;
        foreach (AscatFitSegment s in segs)
        {
            (double nA, double nB) = AscatRawCopyNumbers(s.R, s.B, rho, psi, gamma);
            double rA = RRound(nA), rB = RRound(nB);
            double aberrant = s.B == BalancedBaf ? 0.0 : 1.0;
            totalLength += s.Length;
            ploidyNumerator += (nA + nB) * s.Length;
            double zeros = (rA == 0.0 ? 1.0 : 0.0) + (rB == 0.0 ? 1.0 : 0.0);
            zero += zeros * s.Length;
            zeroAbb += zeros * s.Length * aberrant;
            abbLength += s.Length * aberrant;
            double modA = RMod2(rA), modB = RMod2(rB);
            if ((modA == 0.0 && modB == 1.0) || (modA == 1.0 && modB == 0.0))
            {
                oddEven += s.Length;
            }
        }

        double percentZeroAberrant = zeroAbb / abbLength;
        if (double.IsNaN(percentZeroAberrant))
        {
            percentZeroAberrant = 0.0; // "the next can happen if BAF is a flat line at 0.5"
        }

        return (ploidyNumerator / totalLength, zero / totalLength, percentZeroAberrant, oddEven / totalLength);
    }

    /// <summary>
    /// Integer allele-specific copy number of one segment at the chosen (ρ, ψ), ported verbatim from the
    /// <c>seg_raw</c> construction of ascat.runAscat.R: raw nA/nB, the negative-value correction
    /// (<c>nA+nB &lt; 0 ⇒ 0,0</c>; a negative allele is folded into the other), R half-to-even rounding, and the
    /// balanced-segment odd-total rule (<c>limitround = 0.5</c>: for BAF = 0.5, if nA+nB exceeds the rounded total by
    /// more than 0.5 nA is raised by one; if it falls short by more than 0.5 nB is lowered by one).
    /// </summary>
    /// <param name="r">Segmented logR.</param>
    /// <param name="bAscat">Segmented BAF in ASCAT's ≤ 0.5 orientation (0 when the segment has no heterozygous probe).</param>
    /// <param name="rho">Purity ρ.</param>
    /// <param name="psi">Ploidy parameter ψ.</param>
    /// <param name="gamma">Platform parameter γ.</param>
    /// <param name="haploid">Male X (non-PAR) / Y segment (<see cref="AscatSexModel"/>): runASCAT's non-diploid branch
    /// <c>nAraw = (rho − 1 + ((1 − rho)·2 + rho·psi)·2^(logR/gamma))/rho</c>, <c>nBraw = 0</c>.</param>
    private static (double Major, double Minor) AscatRoundSegment(
        double r, double bAscat, double rho, double psi, double gamma, bool haploid = false)
    {
        (double nAraw, double nBraw) = haploid
            ? (AscatRawHaploidCopyNumber(r, rho, psi, gamma), 0.0)
            : AscatRawCopyNumbers(r, bAscat, rho, psi, gamma);
        if (nAraw + nBraw < 0.0)
        {
            nAraw = 0.0;
            nBraw = 0.0;
        }
        else if (nAraw < 0.0)
        {
            nBraw = nAraw + nBraw;
            nAraw = 0.0;
        }
        else if (nBraw < 0.0)
        {
            nAraw = nAraw + nBraw;
            nBraw = 0.0;
        }

        double rA = RRound(nAraw), rB = RRound(nBraw);
        double nA = rA, nB = rB;
        if (bAscat == BalancedBaf)
        {
            if (nAraw + nBraw > rA + rB + AscatLimitRound)
            {
                nA = rA + 1.0;
            }
            else if (nAraw + nBraw < rA + rB - AscatLimitRound)
            {
                nB = rB - 1.0;
            }
        }

        return (nA, nB);
    }

    /// <summary>Converts an ASCAT integer copy number (≥ 0) to <see cref="int"/>, saturating at Int32.MaxValue (NaN → 0).</summary>
    private static int AscatCopyNumberToInt(double value)
    {
        if (double.IsNaN(value) || value <= 0.0)
        {
            return 0;
        }

        return value >= int.MaxValue ? int.MaxValue : (int)value;
    }

    /// <summary>
    /// Builds the fit result at the selected (ρ, ψ): the integer major/minor segments (ASCAT <c>seg_raw</c> nMajor/nMinor,
    /// one per input summary, all chromosomes), the ASCAT output ploidy — the mean integer total copy number over the
    /// heterozygous probes (<c>ploidy = mean(nA + nB)</c>, here weighted by <c>LocusCount</c>) — and the given GoF.
    /// </summary>
    private static PurityPloidyFit BuildAscatFit(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments, double rho, double psi, double gamma, double goodnessOfFit,
        bool nonAberrant, AscatSexModel sexModel)
    {
        var result = new List<AlleleSpecificSegment>(segments.Count);
        var probeCounts = new int[segments.Count];
        for (int i = 0; i < segments.Count; i++)
        {
            AlleleSpecificSegmentSummary s = segments[i];
            (double major, double minor) = AscatRoundSegment(s.MeanLogR, ToAscatBaf(s.MeanBAF), rho, psi, gamma,
                sexModel.IsHaploid(s.Chromosome, s.Start, s.End));
            int majorInt = AscatCopyNumberToInt(major);
            int minorInt = AscatCopyNumberToInt(minor);
            // Segments with End == Start (single-position) get a 1 bp span so AlleleSpecificSegment.Length > 0.
            long end = s.End > s.Start ? s.End : s.Start + 1;
            result.Add(new AlleleSpecificSegment(s.Chromosome, s.Start, end, majorInt, minorInt));
            probeCounts[i] = s.LocusCount;
        }

        // ASCAT ploidy = mean(nA + nB) over probes = probe-count-weighted mean (canonical overload).
        return new PurityPloidyFit(rho, EstimatePloidy(result, probeCounts), goodnessOfFit, result)
        {
            Psi = psi,
            IsNonAberrant = nonAberrant,
        };
    }

    /// <summary>Validates the segment summaries consumed by the ASCAT fit and returns the autosomal fitting segments.</summary>
    private static AscatFitSegment[] PrepareAscatSegments(IReadOnlyList<AlleleSpecificSegmentSummary> segments)
    {
        ArgumentNullException.ThrowIfNull(segments);
        if (segments.Count == 0)
        {
            throw new ArgumentException("At least one segment is required to fit purity and ploidy.", nameof(segments));
        }

        var autosomal = new List<AscatFitSegment>(segments.Count);
        foreach (AlleleSpecificSegmentSummary s in segments)
        {
            if (s.Chromosome is null)
            {
                throw new ArgumentException("A segment has a null chromosome label.", nameof(segments));
            }

            if (!IsValidAlleleSignal(s.MeanLogR, s.MeanBAF))
            {
                throw new ArgumentException(
                    "Every segment needs a finite mean logR and a mean BAF in [0, 1].", nameof(segments));
            }

            if (s.LocusCount < 1)
            {
                throw new ArgumentException(
                    "Every segment must summarise at least one heterozygous locus (ASCAT weights segments by probe count).",
                    nameof(segments));
            }

            if (!IsAscatSexChromosome(s.Chromosome))
            {
                autosomal.Add(new AscatFitSegment(s.MeanLogR, ToAscatBaf(s.MeanBAF), s.LocusCount));
            }
        }

        if (autosomal.Count == 0)
        {
            throw new ArgumentException(
                "ASCAT fits purity and ploidy on autosomal segments only (sex chromosomes X/Y are excluded); none was supplied.",
                nameof(segments));
        }

        return autosomal.ToArray();
    }

    /// <summary>
    /// Jointly estimates tumour purity ρ and ploidy ψ from segment-level (logR, BAF) summaries with the ASCAT
    /// algorithm (Van Loo et al. 2010, <i>PNAS</i> 107:16910), ported from <c>runASCAT</c> in ascat.runAscat.R
    /// (VanLoo-lab/ascat):
    /// <list type="number">
    /// <item><b>Distance matrix</b> (<c>create_distance_matrix</c>): for ψ ∈ seq(ploidyMin − 0.5, ploidyMax + 0.5,
    /// ploidyStep) × ρ ∈ seq(purityMin, purityMax, purityStep), d(ψ, ρ) = Σ (nMinor − max(round(nMinor), 0))² ·
    /// length · w_b over the <b>autosomal</b> segments, where length = number of heterozygous probes
    /// (<see cref="AlleleSpecificSegmentSummary.LocusCount"/>, ASCAT <c>make_segments</c>) and w_b = 0.05 for BAF = 0.5.</item>
    /// <item><b>Local minima</b>: a grid cell is a candidate when it is the strict minimum of its 7 × 7 neighbourhood.</item>
    /// <item><b>Filter cascade</b> (first pass that yields a candidate wins): (1) non-aberrant sample excluded,
    /// ploidyMin &lt; ploidy &lt; ploidyMax, ρ ≥ 0.2, GoF &gt; 80 %, percentzero &gt; 0.02; (2) 1.7 &lt; ploidy &lt; 2.3 and
    /// perczeroAbb &gt; 0.1; (3) grid columns with ρ &gt; 1 masked (so ρ = 1 can be a border optimum), percentzero /
    /// perczeroAbb / percOddEven alternatives; (4) 1.7 &lt; ploidy &lt; 2.3 only. Here ploidy = Σ(nA+nB)·length/Σlength.</item>
    /// <item><b>Selection</b>: the candidate with the smallest distance (the last one in grid order on an exact tie);
    /// ρ &gt; 1 is reported as 1. GoF = (1 − d/TheoretMaxdist)·100 with TheoretMaxdist = Σ 0.25·length·w_b.</item>
    /// <item><b>Integer segments</b> (ASCAT <c>seg_raw</c>): negative-value correction, R half-to-even rounding and the
    /// balanced odd-total rule, see <see cref="AscatRoundSegment"/>.</item>
    /// </list>
    /// Sex-chromosome (X/Y) segments are excluded from the fit (as in ASCAT, for either sex) and emitted with the diploid
    /// model (ASCAT gender "XX"); pass an <see cref="AscatSexModel"/> (overload, B24 F38) for a male sample.
    /// </summary>
    /// <param name="segments">Segment summaries (from <see cref="SegmentAlleleSpecificAspcf"/> or a caller's segmenter). Non-empty;
    /// each needs a finite mean logR, a mean BAF in [0, 1] and LocusCount ≥ 1; at least one autosomal segment.</param>
    /// <param name="purityMin">Lower bound of the purity grid, in (0, 1] (ASCAT <c>min_purity</c> = 0.1).</param>
    /// <param name="purityMax">Upper bound of the purity grid, finite and ≥ purityMin (ASCAT <c>max_purity</c> = 1.05;
    /// grid points above 1 let ρ = 1 be an interior optimum and are reported as ρ = 1).</param>
    /// <param name="purityStep">Purity grid step (&gt; 0; ASCAT 0.01).</param>
    /// <param name="ploidyMin">ASCAT <c>min_ploidy</c> (&gt; 0, default 1.5): lower ploidy filter; the ψ grid starts at ploidyMin − 0.5.</param>
    /// <param name="ploidyMax">ASCAT <c>max_ploidy</c> (≥ ploidyMin, default 5.5): upper ploidy filter; the ψ grid ends at ploidyMax + 0.5.</param>
    /// <param name="ploidyStep">Ploidy grid step (&gt; 0; ASCAT 0.05).</param>
    /// <param name="gamma">Platform parameter γ (sequencing = <see cref="AscatSequencingGamma"/> = 1).</param>
    /// <returns>The recovered ρ, the ASCAT output ploidy, the percentage GoF, the integer segments, and ψ
    /// (<see cref="PurityPloidyFit.Psi"/>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum ("ASCAT could not find an optimal
    /// ploidy and purity value"); use <see cref="TryFitPurityPloidy"/> to test without an exception.</exception>
    public static PurityPloidyFit FitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => FitPurityPloidy(segments, AscatSexModel.Female, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// ASCAT purity/ploidy fit (see <see cref="FitPurityPloidy(IReadOnlyList{AlleleSpecificSegmentSummary}, double, double, double, double, double, double, double)"/>)
    /// with an explicit sex-chromosome model (runASCAT <c>gender</c> / <c>X_nonPAR</c>, B24 F38). The fit (ρ, ψ, GoF,
    /// non-aberrant flag) is unchanged — sex chromosomes never enter it; only the X/Y integer segments and therefore the
    /// reported ploidy depend on <paramref name="sexModel"/> (see <see cref="AscatSexModel"/>). Each summary is one
    /// runASCAT segment; its [Start, End] is the span tested against the X non-PAR interval.
    /// </summary>
    /// <param name="segments">Segment summaries (as the gender-less overload).</param>
    /// <param name="sexModel">Sex-chromosome model; <see cref="AscatSexModel.Female"/> reproduces the gender-less overload.</param>
    /// <param name="purityMin">See the gender-less overload.</param>
    /// <param name="purityMax">See the gender-less overload.</param>
    /// <param name="purityStep">See the gender-less overload.</param>
    /// <param name="ploidyMin">See the gender-less overload.</param>
    /// <param name="ploidyMax">See the gender-less overload.</param>
    /// <param name="ploidyStep">See the gender-less overload.</param>
    /// <param name="gamma">See the gender-less overload.</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">a grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum.</exception>
    public static PurityPloidyFit FitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        if (!TryFitPurityPloidy(segments, sexModel, out PurityPloidyFit fit, purityMin, purityMax, purityStep,
                ploidyMin, ploidyMax, ploidyStep, gamma))
        {
            throw new InvalidOperationException(
                "ASCAT could not find an optimal ploidy and purity value: no local minimum of the distance matrix passes " +
                "the ASCAT solution filters (ascat.runAscat.R).");
        }

        return fit;
    }

    /// <summary>
    /// ASCAT purity/ploidy fit (see <see cref="FitPurityPloidy"/>) that reports failure instead of throwing:
    /// returns <c>false</c> (and a default <paramref name="fit"/>) when no local minimum of the distance matrix passes
    /// the ASCAT filters — ASCAT's <c>rho = NA</c> outcome. Argument errors still throw.
    /// </summary>
    public static bool TryFitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => TryFitPurityPloidy(segments, AscatSexModel.Female, out fit, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidy(IReadOnlyList{AlleleSpecificSegmentSummary}, AscatSexModel, double, double, double, double, double, double, double)"/>
    /// reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing (B24 F38).
    /// </summary>
    public static bool TryFitPurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);
        if (!TryFindAscatOptimum(s, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma,
                out double rhoOpt, out double psiOpt, out double goodnessOfFit, out bool nonAberrant))
        {
            fit = default;
            return false;
        }

        fit = BuildAscatFit(segments, rhoOpt, psiOpt, gamma, goodnessOfFit, nonAberrant, sexModel);
        return true;
    }

    /// <summary>
    /// The runASCAT grid search (ascat.runAscat.R) on validated autosomal fitting segments: distance matrix, strict
    /// 7 × 7 local minima, the four-pass filter cascade and the optimum selection (see <see cref="FitPurityPloidy"/>).
    /// Returns <c>false</c> when no candidate passes (ASCAT <c>rho = NA</c>); ρ &gt; 1 is reported as 1.
    /// </summary>
    private static bool TryFindAscatOptimum(
        AscatFitSegment[] s, double purityMin, double purityMax, double purityStep,
        double ploidyMin, double ploidyMax, double ploidyStep, double gamma,
        out double rhoOpt, out double psiOpt, out double goodnessOfFitOpt, out bool nonAberrant)
    {
        double[] psiPos = RSeq(ploidyMin - 0.5, ploidyMax + 0.5, ploidyStep);
        double[] rhoPos = RSeq(purityMin, purityMax, purityStep);
        int rows = psiPos.Length, cols = rhoPos.Length;
        var d = new double[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                d[i, j] = AscatDistance(s, rhoPos[j], psiPos[i], gamma);
            }
        }

        (double theoreticalMaxDistance, nonAberrant) = AscatSampleSummary(s);
        bool sampleNonAberrant = nonAberrant;

        var candidates = new List<(double M, int I, int J, double GoodnessOfFit)>();
        bool strictPloidyWindowReachable = ploidyMin < AscatMaxPloidyStrict && ploidyMax > AscatMinPloidyStrict;

        // Pass 1: all filters.
        CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
            !sampleNonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
            && st.GoodnessOfFit > AscatMinGoodnessOfFit && st.PercentZero > AscatMinPercentZero);

        // Pass 2: drop percentzero (allow non-aberrant solutions) with strict ploidy borders.
        if (candidates.Count == 0 && strictPloidyWindowReachable)
        {
            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                st.Ploidy > AscatMinPloidyStrict && st.Ploidy < AscatMaxPloidyStrict && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit && st.PercentZeroAberrant > AscatMinPercentZeroAberrant);
        }

        // Pass 3: allow 100 % aberrant cells — mask the rho > 1 columns so rho = 1 can be a (border) minimum.
        if (candidates.Count == 0)
        {
            for (int j = 0; j < cols; j++)
            {
                if (RDimnameValue(rhoPos[j]) > 1.0)
                {
                    for (int i = 0; i < rows; i++)
                    {
                        d[i, j] = AscatMaskedDistance;
                    }
                }
            }

            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                !sampleNonAberrant && st.Ploidy > ploidyMin && st.Ploidy < ploidyMax && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit
                && (st.PercentZeroAberrant > AscatMinPercentZeroAberrant || st.PercentZero > AscatMinPercentZero
                    || st.PercentOddEven > AscatMinPercentOddEven));
        }

        // Pass 4: drop the percentzero filters, strict ploidy borders.
        if (candidates.Count == 0 && strictPloidyWindowReachable)
        {
            CollectAscatOptima(d, psiPos, rhoPos, s, gamma, theoreticalMaxDistance, candidates, st =>
                st.Ploidy > AscatMinPloidyStrict && st.Ploidy < AscatMaxPloidyStrict && st.Rho >= AscatMinRho
                && st.GoodnessOfFit > AscatMinGoodnessOfFit);
        }

        if (candidates.Count == 0)
        {
            rhoOpt = psiOpt = goodnessOfFitOpt = double.NaN;
            return false;
        }

        // optlim = sort(localmin)[1]; the loop keeps the LAST optimum whose distance equals optlim.
        double optimum = double.PositiveInfinity;
        foreach (var c in candidates)
        {
            optimum = Math.Min(optimum, c.M);
        }

        (double M, int I, int J, double GoodnessOfFit) best = default;
        foreach (var c in candidates)
        {
            if (c.M == optimum)
            {
                best = c;
            }
        }

        psiOpt = RDimnameValue(psiPos[best.I]);
        rhoOpt = Math.Min(1.0, RDimnameValue(rhoPos[best.J])); // if (rho_opt1 > 1) rho_opt1 = 1
        goodnessOfFitOpt = best.GoodnessOfFit;
        return true;
    }

    /// <summary>
    /// Sample-level ASCAT quantities shared by the grid fit (<see cref="TryFitPurityPloidy"/>) and the manual fit
    /// (<see cref="EvaluatePurityPloidy"/>), ascat.runAscat.R: <c>TheoretMaxdist = sum(rep(0.25, n) · length ·
    /// ifelse(b == 0.5, 0.05, 1))</c> (the GoF normaliser) and the <c>nonaberrant</c> flag (MINABB: aberrant probes ≤ 3 %,
    /// MINABBREGION: no aberrant segment above 0.5 % of the probes).
    /// </summary>
    private static (double TheoreticalMaxDistance, bool NonAberrant) AscatSampleSummary(AscatFitSegment[] s)
    {
        double theoreticalMaxDistance = 0.0, totalLength = 0.0, aberrantLength = 0.0, maxAberrantSegment = 0.0;
        foreach (AscatFitSegment seg in s)
        {
            bool balanced = seg.B == BalancedBaf;
            theoreticalMaxDistance += AscatWorstCaseIntegerDistance * seg.Length * (balanced ? AscatBalancedSegmentWeight : 1.0);
            totalLength += seg.Length;
            if (!balanced)
            {
                aberrantLength += seg.Length;
                maxAberrantSegment = Math.Max(maxAberrantSegment, seg.Length);
            }
        }

        bool nonAberrant = aberrantLength / totalLength <= AscatMinAberrantFraction
                           && maxAberrantSegment / totalLength <= AscatMinAberrantRegionFraction;
        return (theoreticalMaxDistance, nonAberrant);
    }

    /// <summary>Statistics of one ASCAT local-minimum candidate, as tested by the runASCAT filter passes.</summary>
    private readonly record struct AscatCandidate(
        double Rho, double Ploidy, double GoodnessOfFit, double PercentZero, double PercentZeroAberrant, double PercentOddEven);

    /// <summary>
    /// One runASCAT filter pass: scans the interior of the distance matrix (i ∈ 4..nrow−3, j ∈ 4..ncol−3 in R's 1-based
    /// indices) for cells that are the strict minimum of their 7 × 7 window (<c>seld[4,4] = max(seld); min(seld) &gt; m</c>)
    /// and appends those whose statistics pass <paramref name="accept"/>.
    /// </summary>
    private static void CollectAscatOptima(
        double[,] d, double[] psiPos, double[] rhoPos, AscatFitSegment[] segs, double gamma, double theoreticalMaxDistance,
        List<(double M, int I, int J, double GoodnessOfFit)> candidates, Func<AscatCandidate, bool> accept)
    {
        int rows = d.GetLength(0), cols = d.GetLength(1);
        const int h = AscatLocalMinimumHalfWindow;
        for (int i = h; i < rows - h; i++)
        {
            for (int j = h; j < cols - h; j++)
            {
                double m = d[i, j];
                double windowMax = double.NegativeInfinity;
                for (int di = -h; di <= h; di++)
                {
                    for (int dj = -h; dj <= h; dj++)
                    {
                        windowMax = Math.Max(windowMax, d[i + di, j + dj]);
                    }
                }

                bool strictMinimum = true;
                for (int di = -h; di <= h && strictMinimum; di++)
                {
                    for (int dj = -h; dj <= h; dj++)
                    {
                        double v = di == 0 && dj == 0 ? windowMax : d[i + di, j + dj];
                        if (!(v > m))
                        {
                            strictMinimum = false;
                            break;
                        }
                    }
                }

                if (!strictMinimum)
                {
                    continue;
                }

                double psi = RDimnameValue(psiPos[i]);
                double rho = RDimnameValue(rhoPos[j]);
                var stats = AscatCandidateStatistics(segs, rho, psi, gamma);
                double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
                var candidate = new AscatCandidate(
                    rho, stats.Ploidy, goodnessOfFit, stats.PercentZero, stats.PercentZeroAberrant, stats.PercentOddEven);
                if (accept(candidate))
                {
                    candidates.Add((m, i, j, goodnessOfFit));
                }
            }
        }
    }

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (<c>rho_manual</c> / <c>psi_manual</c> in ascat.runAscat.R): skips the
    /// grid search and returns the goodness of fit — d computed as in <c>create_distance_matrix</c> over the autosomal
    /// segments, GoF = (1 − d/TheoretMaxdist)·100 — together with the integer segments and ASCAT ploidy at (ρ, ψ).
    /// </summary>
    /// <param name="segments">Segment summaries (same requirements as <see cref="FitPurityPloidy"/>).</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
        => EvaluatePurityPloidy(segments, AscatSexModel.Female, purity, ploidy, gamma);

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (as the gender-less overload) and an explicit sex-chromosome model
    /// for the emitted X/Y segments (<see cref="AscatSexModel"/>, B24 F38).
    /// </summary>
    /// <param name="segments">Segment summaries.</param>
    /// <param name="sexModel">Sex-chromosome model.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segments"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="segments"/> is empty, malformed, or has no autosomal segment.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidy(
        IReadOnlyList<AlleleSpecificSegmentSummary> segments,
        AscatSexModel sexModel,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = PrepareAscatSegments(segments);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);
        double m = AscatDistance(s, purity, ploidy, gamma);
        double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
        return BuildAscatFit(segments, purity, ploidy, gamma, goodnessOfFit, nonAberrant, sexModel);
    }

    /// <summary>
    /// ASCAT purity/ploidy fit (runASCAT, ascat.runAscat.R) on a germline-aware ASPCF segmentation
    /// (<see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, double)"/>), i.e.
    /// with germline-homozygous probes, exactly as runASCAT treats them:
    /// <list type="bullet">
    /// <item><b>Fit</b>: the distance matrix, filters, goodness of fit and non-aberrant flag use the heterozygous
    /// autosomal probes only (<c>r = lrrsegmented[names(bafsegmented)]</c>), grouped by <c>make_segments</c> into runs of
    /// identical (segmented logR, segmented BAF) with length = number of heterozygous probes. Homozygous probes enter
    /// only through the logR levels; the grid search is the one of <see cref="FitPurityPloidy"/>.</item>
    /// <item><b>Segments</b> (<c>seg_raw</c>): one per <see cref="AspcfSegmentation.Segments"/> entry (runs of equal
    /// segmented logR within a chromosome), using the BAF of the run's first heterozygous probe; a run without
    /// heterozygous probes (BAF NA) uses ASCAT's <c>bafke = 0</c>, which after the negative-value correction puts the
    /// whole total copy number on the major allele (only nA + nB is meaningful there).</item>
    /// <item><b>Ploidy</b>: ASCAT's reported <c>ploidy = mean(nA + nB)</c> over <b>all</b> probes (heterozygous and
    /// homozygous, sex chromosomes included), i.e. the <see cref="AspcfSegment.LocusCount"/>-weighted mean total copy
    /// number of the segments.</item>
    /// </list>
    /// Parameters and the no-optimum behaviour are those of <see cref="FitPurityPloidy"/>.
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="purityMin">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="purityMax">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="purityStep">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyMin">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyMax">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="ploidyStep">See <see cref="FitPurityPloidy"/>.</param>
    /// <param name="gamma">See <see cref="FitPurityPloidy"/>.</param>
    /// <returns>ρ, the all-probe ASCAT ploidy, GoF, the <c>seg_raw</c> integer segments and ψ.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum (use
    /// <see cref="TryFitPurityPloidyFromAspcf"/>).</exception>
    public static PurityPloidyFit FitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => FitPurityPloidyFromAspcf(segmentation, AscatSexModel.Female, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf(AspcfSegmentation, double, double, double, double, double, double, double)"/>
    /// with an explicit sex-chromosome model (runASCAT <c>gender</c> / <c>X_nonPAR</c>, B24 F38): the fit is unchanged
    /// (sex chromosomes never enter it); the X/Y <c>seg_raw</c> segments — and so the all-probe ploidy — follow
    /// <paramref name="sexModel"/> (see <see cref="AscatSexModel"/>). For the X non-PAR rule a segment's span is its
    /// first/last probe position (<see cref="AspcfSegment.Start"/>/<see cref="AspcfSegment.End"/>, i.e. the run of
    /// equal segmented logR that <c>diploidprobes_fixnonPAR</c> builds with <c>rle</c>).
    /// <para>For a male with <c>X_nonPAR</c>, <c>ascat.aspcf</c> additionally re-labels the germline genotypes of non-PAR X
    /// probes (all homozygous, then a random autosome-matched fraction heterozygous); that step belongs to the
    /// segmentation — run <see cref="SegmentAlleleSpecificAspcf(IEnumerable{AlleleSpecificLocus}, IReadOnlyList{bool}, AscatMaleXGenotyping, double)"/>
    /// (B24 F58). The fit uses <see cref="AspcfSegmentation.GermlineHeterozygous"/> as given.</para>
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="sexModel">Sex-chromosome model; <see cref="AscatSexModel.Female"/> reproduces the gender-less overload.</param>
    /// <param name="purityMin">See the gender-less overload.</param>
    /// <param name="purityMax">See the gender-less overload.</param>
    /// <param name="purityStep">See the gender-less overload.</param>
    /// <param name="ploidyMin">See the gender-less overload.</param>
    /// <param name="ploidyMax">See the gender-less overload.</param>
    /// <param name="ploidyStep">See the gender-less overload.</param>
    /// <param name="gamma">See the gender-less overload.</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A grid bound or step is out of range.</exception>
    /// <exception cref="InvalidOperationException">ASCAT finds no acceptable optimum.</exception>
    public static PurityPloidyFit FitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        if (!TryFitPurityPloidyFromAspcf(segmentation, sexModel, out PurityPloidyFit fit, purityMin, purityMax, purityStep,
                ploidyMin, ploidyMax, ploidyStep, gamma))
        {
            throw new InvalidOperationException(
                "ASCAT could not find an optimal ploidy and purity value: no local minimum of the distance matrix passes " +
                "the ASCAT solution filters (ascat.runAscat.R).");
        }

        return fit;
    }

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf"/> reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing.
    /// </summary>
    public static bool TryFitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
        => TryFitPurityPloidyFromAspcf(segmentation, AscatSexModel.Female, out fit, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);

    /// <summary>
    /// <see cref="FitPurityPloidyFromAspcf(AspcfSegmentation, AscatSexModel, double, double, double, double, double, double, double)"/>
    /// reporting ASCAT's <c>rho = NA</c> outcome as <c>false</c> instead of throwing (B24 F38).
    /// </summary>
    public static bool TryFitPurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        out PurityPloidyFit fit,
        double purityMin = 0.1,
        double purityMax = 1.05,
        double purityStep = 0.01,
        double ploidyMin = 1.5,
        double ploidyMax = 5.5,
        double ploidyStep = 0.05,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = AscatMakeSegments(segmentation);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateGrid(purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma);
        if (!TryFindAscatOptimum(s, purityMin, purityMax, purityStep, ploidyMin, ploidyMax, ploidyStep, gamma,
                out double rhoOpt, out double psiOpt, out double goodnessOfFit, out bool nonAberrant))
        {
            fit = default;
            return false;
        }

        fit = BuildAscatFitFromAspcf(segmentation, rhoOpt, psiOpt, gamma, goodnessOfFit, nonAberrant, sexModel);
        return true;
    }

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy (<c>rho_manual</c>/<c>psi_manual</c>) on a germline-aware ASPCF
    /// segmentation: GoF over the heterozygous autosomal <c>make_segments</c> runs, <c>seg_raw</c> segments and the
    /// all-probe ploidy as in <see cref="FitPurityPloidyFromAspcf"/>.
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
        => EvaluatePurityPloidyFromAspcf(segmentation, AscatSexModel.Female, purity, ploidy, gamma);

    /// <summary>
    /// ASCAT with a user-supplied purity and ploidy on a germline-aware ASPCF segmentation (as the gender-less overload)
    /// with an explicit sex-chromosome model for the X/Y <c>seg_raw</c> segments (<see cref="AscatSexModel"/>, B24 F38).
    /// </summary>
    /// <param name="segmentation">The germline-aware ASPCF segmentation.</param>
    /// <param name="sexModel">Sex-chromosome model.</param>
    /// <param name="purity">Purity ρ ∈ (0, 1].</param>
    /// <param name="ploidy">Ploidy parameter ψ (&gt; 0).</param>
    /// <param name="gamma">Platform parameter γ (&gt; 0).</param>
    /// <exception cref="ArgumentNullException"><paramref name="segmentation"/> or <paramref name="sexModel"/> is null.</exception>
    /// <exception cref="ArgumentException">The segmentation has no heterozygous autosomal probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException">ρ ∉ (0, 1], ψ ≤ 0 or γ ≤ 0 (or any is non-finite).</exception>
    public static PurityPloidyFit EvaluatePurityPloidyFromAspcf(
        AspcfSegmentation segmentation,
        AscatSexModel sexModel,
        double purity,
        double ploidy,
        double gamma = AscatSequencingGamma)
    {
        AscatFitSegment[] s = AscatMakeSegments(segmentation);
        ArgumentNullException.ThrowIfNull(sexModel);
        ValidateAscatModelParameters(purity, ploidy, gamma);

        (double theoreticalMaxDistance, bool nonAberrant) = AscatSampleSummary(s);
        double m = AscatDistance(s, purity, ploidy, gamma);
        double goodnessOfFit = (1.0 - m / theoreticalMaxDistance) * 100.0;
        return BuildAscatFitFromAspcf(segmentation, purity, ploidy, gamma, goodnessOfFit, nonAberrant, sexModel);
    }

    /// <summary>
    /// runASCAT fitting segments from a germline-aware segmentation: the heterozygous autosomal probes' (segmented logR,
    /// 1 − mirrored segmented BAF) pairs, grouped by ASCAT <c>make_segments</c> into runs of identical pairs (the runs
    /// ignore chromosome boundaries, as in R); length = number of probes in the run.
    /// </summary>
    private static AscatFitSegment[] AscatMakeSegments(AspcfSegmentation segmentation)
    {
        ArgumentNullException.ThrowIfNull(segmentation);
        var result = new List<AscatFitSegment>();
        double previousR = 1e10, previousB = -1.0;
        int count = 0;
        for (int i = 0; i < segmentation.Loci.Count; i++)
        {
            // SNPposhet = SNPpos[names(bafsegmented), ]: heterozygous probes WITH a segmented BAF (a missing tumour
            // BAF/logR drops a heterozygous probe from Tumor_BAF_segmented, B24 F59).
            if (!segmentation.GermlineHeterozygous[i] || double.IsNaN(segmentation.SegmentedBaf[i])
                || IsAscatSexChromosome(segmentation.Loci[i].Chromosome))
            {
                continue;
            }

            double r = segmentation.SegmentedLogR[i];
            double b = 1.0 - segmentation.SegmentedBaf[i]; // Tumor_BAF_segmented = 1 - bafPCFed
            if (b != previousB || r != previousR)
            {
                if (count > 0)
                {
                    result.Add(new AscatFitSegment(previousR, previousB, count));
                }

                count = 0;
            }

            count++;
            previousR = r;
            previousB = b;
        }

        if (count == 0)
        {
            throw new ArgumentException(
                "ASCAT fits purity and ploidy on heterozygous autosomal probes; the segmentation has none.",
                nameof(segmentation));
        }

        result.Add(new AscatFitSegment(previousR, previousB, count));
        return result.ToArray();
    }

    /// <summary>
    /// runASCAT output at the selected (ρ, ψ) for a germline-aware segmentation: <c>seg_raw</c> integer segments (one per
    /// logR segment; <c>bafke</c> = first heterozygous BAF, 0 when none) and <c>ploidy = mean(nA + nB, na.rm = TRUE)</c> over
    /// all probes, i.e. without those whose raw logR is missing or that have a segmented BAF but a missing raw BAF.
    /// </summary>
    private static PurityPloidyFit BuildAscatFitFromAspcf(
        AspcfSegmentation segmentation, double rho, double psi, double gamma, double goodnessOfFit, bool nonAberrant,
        AscatSexModel sexModel)
    {
        var result = new List<AlleleSpecificSegment>(segmentation.Segments.Count);
        var counted = new List<AlleleSpecificSegment>(segmentation.Segments.Count);
        var probeCounts = new List<int>(segmentation.Segments.Count);
        int locus = 0;
        for (int i = 0; i < segmentation.Segments.Count; i++)
        {
            AspcfSegment s = segmentation.Segments[i];

            // ploidy = mean(nA + nB, na.rm = TRUE) over the probes (B24 F59): n1all/n2all are NA where the raw logR is NA, and
            // at a probe with a segmented BAF whose raw BAF is NA (ifelse(baf <= 0.5, …) on NA — ascat.asmultipcf gives
            // such a probe a BAF from the other samples).
            int observed = 0;
            for (int k = 0; k < s.LocusCount; k++, locus++)
            {
                AlleleSpecificLocus l = segmentation.Loci[locus];
                bool bafless = !double.IsNaN(segmentation.SegmentedBaf[locus]) && double.IsNaN(l.BAF);
                if (!double.IsNaN(l.LogR) && !bafless)
                {
                    observed++;
                }
            }

            double bafke = s.HasBaf ? 1.0 - s.MeanBAF : 0.0; // "if (is.na(bafke)) bafke = 0"
            (double major, double minor) = AscatRoundSegment(s.MeanLogR, bafke, rho, psi, gamma,
                sexModel.IsHaploid(s.Chromosome, s.Start, s.End));
            long end = s.End > s.Start ? s.End : s.Start + 1;
            var segment = new AlleleSpecificSegment(
                s.Chromosome, s.Start, end, AscatCopyNumberToInt(major), AscatCopyNumberToInt(minor));
            result.Add(segment);
            if (observed > 0)
            {
                counted.Add(segment);
                probeCounts.Add(observed);
            }
        }

        return new PurityPloidyFit(rho, EstimatePloidy(counted, probeCounts), goodnessOfFit, result)
        {
            Psi = psi,
            IsNonAberrant = nonAberrant,
        };
    }

    /// <summary>
    /// Validates a fixed ASCAT/Battenberg model point (ρ, ψ, γ) — shared by <see cref="EvaluatePurityPloidy"/> and
    /// <see cref="FitSubclonalCopyNumber"/>: ρ ∈ (0, 1], ψ &gt; 0 and γ &gt; 0, all finite (the ASCAT equations divide by ρ
    /// and by γ).
    /// </summary>
    private static void ValidateAscatModelParameters(double purity, double ploidy, double gamma)
    {
        if (!double.IsFinite(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity ρ must be in (0, 1].");
        }

        if (!double.IsFinite(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy ψ must be positive and finite.");
        }

        if (!double.IsFinite(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive and finite.");
        }
    }

    private static void ValidateGrid(
        double purityMin, double purityMax, double purityStep,
        double ploidyMin, double ploidyMax, double ploidyStep, double gamma)
    {
        if (!double.IsFinite(purityMin) || purityMin <= 0.0 || purityMin > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMin), purityMin, "purityMin must be in (0, 1].");
        }

        if (!double.IsFinite(purityMax) || purityMax < purityMin)
        {
            throw new ArgumentOutOfRangeException(nameof(purityMax), purityMax, "purityMax must be finite and ≥ purityMin.");
        }

        if (!double.IsFinite(purityStep) || purityStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purityStep), purityStep, "purityStep must be positive.");
        }

        if (!double.IsFinite(ploidyMin) || ploidyMin <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMin), ploidyMin, "ploidyMin must be positive.");
        }

        if (!double.IsFinite(ploidyMax) || ploidyMax < ploidyMin)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyMax), ploidyMax, "ploidyMax must be finite and ≥ ploidyMin.");
        }

        if (!double.IsFinite(ploidyStep) || ploidyStep <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidyStep), ploidyStep, "ploidyStep must be positive.");
        }

        if (!double.IsFinite(gamma) || gamma <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be positive.");
        }

        // Bound the grid so a pathological step cannot allocate an unbounded distance matrix.
        double cells = ((ploidyMax - ploidyMin + 1.0) / ploidyStep + 1.0) * ((purityMax - purityMin) / purityStep + 1.0);
        if (cells > AscatMaxGridCells)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purityStep), purityStep, $"The (ρ, ψ) grid would have {cells:G3} cells (limit {AscatMaxGridCells:G3}).");
        }
    }

    /// <summary>Upper bound on the number of (ρ, ψ) grid cells evaluated (ASCAT default grid: 101 × 96 = 9 696).</summary>
    private const double AscatMaxGridCells = 4_000_000;

    /// <summary>
    /// Derives the integer mutation multiplicity m (number of mutated copies per cancer cell) of a somatic
    /// variant from its VAF, the tumour purity ρ, and the local total / major copy number, so that
    /// <see cref="EstimateCcf"/> can be driven without a caller-supplied multiplicity. The expected number of
    /// mutated copies for a clonal mutation is n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)] (McGranahan et al. 2016,
    /// <i>Science</i> 351:1463; equivalently the inversion of the PICTograph model VAF = m·CCF·ρ /
    /// (N_T·ρ + 2(1−ρ)) at CCF = 1, Zheng et al. 2022, <i>Bioinformatics</i> 38:3677). The result is rounded to
    /// the nearest integer with ties to even and floored at 1, exactly as facets-suite
    /// <c>expected_mutant_copies</c> (mskcc/facets-suite <c>R/ccf-annotate-maf.R</c>, "Based on PMID 28270531":
    /// <c>mu &lt; 1 → 1</c>, then R <c>round</c> = IEC 60559 half-to-even, so n_mut = 2.5 → 2 and 1.5 → 2). In
    /// addition (not in facets-suite) the result is capped at majorCopyNumber (a variant present on at least one
    /// copy cannot exceed the major-allele copy number).
    /// </summary>
    /// <param name="vaf">Observed variant allele fraction ∈ [0, 1].</param>
    /// <param name="purity">Tumour purity ρ ∈ (0, 1].</param>
    /// <param name="totalCopyNumber">Local tumour total copy number N_T (≥ 1).</param>
    /// <param name="majorCopyNumber">Local major-allele copy number, the upper bound on multiplicity (in [1, N_T]).</param>
    /// <returns>The integer mutation multiplicity m ∈ [1, majorCopyNumber].</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0,1], purity ∉ (0,1], totalCopyNumber &lt; 1, or majorCopyNumber ∉ [1, totalCopyNumber].</exception>
    public static int DeriveMultiplicity(double vaf, double purity, int totalCopyNumber, int majorCopyNumber)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(purity), purity, "Purity must be in (0, 1].");
        }

        if (totalCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(totalCopyNumber), totalCopyNumber, "Total copy number must be ≥ 1.");
        }

        if (majorCopyNumber < 1 || majorCopyNumber > totalCopyNumber)
        {
            throw new ArgumentOutOfRangeException(
                nameof(majorCopyNumber), majorCopyNumber, $"Major copy number must be in [1, {totalCopyNumber}].");
        }

        // n_mut = VAF·(1/ρ)·[ρ·N_T + 2(1−ρ)] — McGranahan 2016 observed mutation copy number (CCF=1 ⇒ m = n_mut).
        // This is exactly the canonical CNAqc purity/copy-number VAF correction (AdjustVAFForPurity, ONCO-VAF-001).
        double rawMultiplicity = AdjustVAFForPurity(vaf, purity, totalCopyNumber);
        // facets-suite expected_mutant_copies: R round() is IEC 60559 half-to-even (2.5 → 2), and mu < 1 → 1.
        int rounded = (int)Math.Round(rawMultiplicity, MidpointRounding.ToEven);
        // Clamp to [1, major CN]: the floor at 1 is facets-suite's mu < 1 → 1 (round(mu) ≤ 1 there); the cap at the
        // major-allele copy number is a documented extra (facets-suite leaves expected_alt_copies uncapped).
        return Math.Clamp(rounded, 1, majorCopyNumber);
    }
}
