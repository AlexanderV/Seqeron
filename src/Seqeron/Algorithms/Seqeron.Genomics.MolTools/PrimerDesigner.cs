using System.Text;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs PCR primers for DNA sequences with various quality criteria.
/// </summary>
public static partial class PrimerDesigner
{
    /// <summary>
    /// Default primer design parameters (library conventions — <c>MaxDinucleotideRepeats</c> = 4 is an unsourced
    /// library screen, see <see cref="FindLongestDinucleotideRepeat"/>; the Primer3 3′-end checks PRIMER_GC_CLAMP,
    /// PRIMER_MAX_END_GC and PRIMER_MAX_END_STABILITY are at their Primer3 defaults, i.e. inactive;
    /// <c>Check3PrimeStability</c> is deprecated and has no effect, see <see cref="PrimerParameters"/>).
    /// </summary>
    public static readonly PrimerParameters DefaultParameters = new(
        MinLength: 18,
        MaxLength: 25,
        OptimalLength: 20,
        MinGcContent: 40,
        MaxGcContent: 60,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 4,
        MaxDinucleotideRepeats: 4,
        Avoid3PrimeGC: false,
        Check3PrimeStability: true
    );

    /// <summary>
    /// Primer3's default per-primer settings (<c>libprimer3.cc</c> <c>pr_set_default_global_args_1/_2</c>,
    /// primer3-py 2.3.1 <c>design_primers</c> defaults): PRIMER_MIN/OPT/MAX_SIZE = 18/20/27,
    /// PRIMER_MIN/MAX_GC = 20/80 %, PRIMER_MIN/OPT/MAX_TM = 57/60/63 °C, PRIMER_MAX_POLY_X = 5, no
    /// dinucleotide-repeat limit (Primer3 has none), PRIMER_GC_CLAMP = 0, PRIMER_MAX_END_GC = 5 and
    /// PRIMER_MAX_END_STABILITY = 100 (<see cref="PrimerParameters.GcClamp"/>, <see cref="PrimerParameters.MaxEndGc"/>,
    /// <see cref="PrimerParameters.MaxEndStability"/> left at their Primer3 defaults; 100 cannot be exceeded by an ACGT
    /// 3′ pentamer), thermodynamic structure screen with every limit 47 °C. With <see cref="PrimerPairOptions.Primer3Defaults"/> this makes
    /// <see cref="DesignPrimerPairs"/> reproduce <c>primer3.design_primers</c> run with only
    /// SEQUENCE_TEMPLATE / SEQUENCE_TARGET set. (<see cref="DefaultParameters"/> keeps the library's
    /// stricter conventions: 18–25 nt, 40–60 % GC, poly-X ≤ 4, dinucleotide repeat ≤ 4.)
    /// </summary>
    public static readonly PrimerParameters Primer3DefaultParameters = new(
        MinLength: 18,
        MaxLength: 27,
        OptimalLength: 20,
        MinGcContent: 20,
        MaxGcContent: 80,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 5,
        MaxDinucleotideRepeats: int.MaxValue,
        Avoid3PrimeGC: false,
        Check3PrimeStability: false
    );

    /// <summary>
    /// Designs the best forward/reverse primer pair for a target region with Primer3's pair-selection
    /// algorithm (Untergasser et al. 2012; <c>libprimer3.cc</c> <c>make_detection_primer_lists</c>,
    /// <c>choose_pair_or_triple</c>, <c>characterize_pair</c>, <c>obj_fn</c>, <c>compare_primer_pair</c>)
    /// — the first pair (rank 0) of <see cref="DesignPrimerPairs"/>:
    /// <list type="number">
    /// <item><b>Search region</b> (Primer3 <c>SEQUENCE_TARGET</c> / <c>SEQUENCE_INCLUDED_REGION</c> /
    /// <c>PRIMER_PRODUCT_SIZE_RANGE</c>): forward candidates end at or before <paramref name="targetStart"/>,
    /// reverse candidates start at or after <paramref name="targetEnd"/> (primers never overlap the target),
    /// both inside <see cref="PrimerPairOptions.IncludedRegion"/> (default: the whole template); the
    /// product-size ranges (default 100–300 bp, Primer3's default) bound the product. Every candidate is
    /// evaluated by <see cref="EvaluatePrimer"/> (per-primer limits of <paramref name="parameters"/>) and
    /// carries its Primer3 per-primer penalty (<see cref="PrimerCandidate.Penalty"/>).</item>
    /// <item><b>Pair search</b>: candidates sorted as Primer3's <c>sort_primer_array</c>; pairs are examined
    /// with Primer3's pruning and must satisfy, in <c>characterize_pair</c> order, the product size range
    /// (ranges are tried in order: the next range is used only when no pair fits the current one), the
    /// product Tm limits (<see cref="PrimerPairOptions.ProductMinTm"/>/<see cref="PrimerPairOptions.ProductMaxTm"/>,
    /// Primer3 <c>long_seq_tm</c>), |Tm_f − Tm_r| ≤ <see cref="PrimerPairOptions.MaxTmDifference"/>
    /// (PRIMER_PAIR_MAX_DIFF_TM), the per-primer secondary-structure screen and the pair complementarity
    /// screen — with the default <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> Primer3's ntthal
    /// limits (self-any/self-end/hairpin and pair compl-any/compl-end Tm ≤ 47 °C), evaluated lazily as
    /// <c>characterize_pair</c> does; with <see cref="PrimerStructureScreen.Heuristic"/>
    /// <see cref="HasHairpinPotential"/> and <see cref="HasPrimerDimer"/> — and, with
    /// <see cref="PrimerPairOptions.PickInternalOligo"/>, a hybridization oligo between the primers
    /// (<c>choose_internal_oligo</c>).</item>
    /// <item><b>Objective</b>: the pair penalty is Primer3's <c>obj_fn</c> with
    /// <see cref="PrimerPairOptions.Weights"/> (default: PRIMER_PAIR_WT_PR_PENALTY = 1, all other pair weights
    /// 0, i.e. the sum of the two primer penalties); ties within 1e-6 are broken as
    /// <c>compare_primer_pair</c> (left primer further 3′, then right primer further 5′, then shorter left,
    /// then shorter right).</item>
    /// </list>
    /// Coordinates: the target is the half-open interval [<paramref name="targetStart"/>,
    /// <paramref name="targetEnd"/>) (0-based). Reverse-primer <see cref="PrimerCandidate.Position"/> is
    /// the leftmost template coordinate of its binding site; the product size is
    /// <c>reverse.Position + reverse.Length − forward.Position</c> (Primer3 PRIMER_PAIR_PRODUCT_SIZE).
    /// When candidates exist but no pair satisfies the pair constraints, the individually lowest-penalty
    /// forward and reverse candidates are returned with <c>IsValid = false</c> and a message naming the
    /// violated constraint. Library defaults that differ from Primer3 (documented): the per-primer limits of
    /// <see cref="DefaultParameters"/> and PRIMER_PAIR_MAX_DIFF_TM = 5 °C (Primer3: 100); pass
    /// <see cref="Primer3DefaultParameters"/> and <see cref="PrimerPairOptions.Primer3Defaults"/> for
    /// Primer3's defaults. Reaction conditions: the primer Tm, the ntthal structure / pair complementarity values and
    /// the product Tm use the primer conditions of <paramref name="parameters"/> (PRIMER_SALT_MONOVALENT,
    /// PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC, PRIMER_DNA_CONC; Primer3 defaults 50 mM / 1.5 mM / 0.6 mM / 50 nM), the
    /// internal oligo those of <see cref="PrimerPairOptions.InternalOligo"/> (PRIMER_INTERNAL_*; defaults 50 mM / 0 / 0 /
    /// 50 nM). Verified against primer3-py 2.3.1 <c>design_primers</c>, including random reaction conditions.
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>Primer pair result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">Invalid target or options (see <see cref="PrimerPairOptions"/>).</exception>
    public static PrimerPairResult DesignPrimers(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        return DesignPrimersCore(template.Sequence, null, targetStart, targetEnd, parameters, pairOptions);
    }

    /// <summary>
    /// <see cref="DesignPrimers(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)"/> on a case-preserving
    /// template string (Primer3 SEQUENCE_TEMPLATE as given): with <see cref="PrimerPairOptions.LowercaseMasking"/>
    /// (PRIMER_LOWERCASE_MASKING) or <see cref="PrimerPairOptions.MaskTemplate"/> a primer or internal oligo whose 3′-terminal
    /// template base is lower case is rejected (Primer3 <c>is_lowercase_masked</c>); otherwise case is ignored and the result
    /// equals the <see cref="DnaSequence"/> overload. All other computations use the upper-cased template.
    /// </summary>
    /// <param name="template">Template (A/C/G/T, any case; lower case marks masked bases).</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>Primer pair result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">A character other than A/C/G/T, an invalid target or options.</exception>
    public static PrimerPairResult DesignPrimers(
        string template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        return DesignPrimersCore(new DnaSequence(template).Sequence, template, targetStart, targetEnd, parameters, pairOptions);
    }

    private static PrimerPairResult DesignPrimersCore(string seq, string? caseTemplate, int targetStart, int targetEnd,
        PrimerParameters? parameters, PrimerPairOptions? pairOptions)
    {
        var search = new PrimerPairSearch(seq, caseTemplate, targetStart, targetEnd,
            parameters ?? DefaultParameters, pairOptions ?? PrimerPairOptions.Default);
        var pairs = search.Run(1);
        return pairs.Count > 0 ? pairs[0] : search.Failure();
    }

    /// <summary>
    /// Returns up to <see cref="PrimerPairOptions.NumReturn"/> primer pairs (Primer3 PRIMER_NUM_RETURN,
    /// default 5), best first, exactly as Primer3's <c>choose_pair_or_triple</c> picks them (see
    /// <see cref="DesignPrimers"/>): after a pair is selected it is removed and the search is repeated;
    /// by default primers may be reused in later pairs (Primer3 PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE = −1),
    /// otherwise a selected pair excludes every left / right primer whose 3′ end lies closer than
    /// <see cref="PrimerPairOptions.MinLeftThreePrimeDistance"/> / <see cref="PrimerPairOptions.MinRightThreePrimeDistance"/>
    /// to its own (0: only the identical primer) from all later pairs; product-size ranges are tried in order. Each result is valid and carries
    /// <see cref="PrimerPairResult.PairPenalty"/> (PRIMER_PAIR_k_PENALTY),
    /// <see cref="PrimerPairResult.ProductTm"/> (PRIMER_PAIR_k_PRODUCT_TM), the pair complementarity Tm
    /// values and, with <see cref="PrimerPairOptions.PickInternalOligo"/>, the internal oligo
    /// (PRIMER_INTERNAL_k_*). An empty list means no pair satisfies the constraints
    /// (<see cref="DesignPrimers"/> then reports why).
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>The selected pairs, rank 0 first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">Invalid target or options.</exception>
    public static IReadOnlyList<PrimerPairResult> DesignPrimerPairs(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var opts = pairOptions ?? PrimerPairOptions.Default;
        return new PrimerPairSearch(template.Sequence, null, targetStart, targetEnd, parameters ?? DefaultParameters, opts)
            .Run(opts.NumReturn);
    }

    /// <summary>
    /// <see cref="DesignPrimerPairs(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)"/> on a case-preserving
    /// template string: with <see cref="PrimerPairOptions.LowercaseMasking"/> (PRIMER_LOWERCASE_MASKING) or
    /// <see cref="PrimerPairOptions.MaskTemplate"/> primers and internal oligos whose 3′-terminal template base is lower case
    /// are rejected (Primer3 <c>is_lowercase_masked</c>); otherwise case is ignored.
    /// </summary>
    /// <param name="template">Template (A/C/G/T, any case; lower case marks masked bases).</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>The selected pairs, rank 0 first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">A character other than A/C/G/T, an invalid target or options.</exception>
    public static IReadOnlyList<PrimerPairResult> DesignPrimerPairs(
        string template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var opts = pairOptions ?? PrimerPairOptions.Default;
        return new PrimerPairSearch(new DnaSequence(template).Sequence, template, targetStart, targetEnd,
            parameters ?? DefaultParameters, opts).Run(opts.NumReturn);
    }

    /// <summary>
    /// Former fixed search flank (bp). <see cref="DesignPrimers"/> now uses Primer3's search region
    /// (included region + <see cref="PrimerPairOptions.ProductSizeRanges"/>); kept for source compatibility.
    /// </summary>
    [Obsolete("DesignPrimers now uses Primer3's search region: PrimerPairOptions.IncludedRegion and ProductSizeRanges.")]
    public const int PrimerSearchFlank = 200;

    /// <summary>
    /// Library default maximum |Tm_forward − Tm_reverse| (°C) for a primer pair (Addgene "within 5 °C";
    /// stricter than Primer3's PRIMER_PAIR_MAX_DIFF_TM default of 100) — the default of
    /// <see cref="PrimerPairOptions.MaxTmDifference"/>.
    /// </summary>
    public const double MaxPairTmDifference = 5.0;

    /// <summary>Primer3's default PRIMER_PAIR_MAX_DIFF_TM (°C), <c>pr_set_default_global_args_1</c>.</summary>
    public const double Primer3MaxPairTmDifference = 100.0;

    /// <summary>
    /// Primer3 product melting temperature (<c>oligotm.c</c> <c>long_seq_tm</c>, PRIMER_PAIR_k_PRODUCT_TM):
    /// Tm = 81.5 + 16.6·log10([Mon]_eq/1000) + 41·(G+C)/N − 600/N, with [Mon]_eq = [Mon] +
    /// 120·√([Mg²⁺] − [dNTP]) (mM) and N the product length (G/C counted case-insensitively; every
    /// character counts towards N). Defaults are Primer3's primer conditions (50 mM, 1.5 mM Mg²⁺, 0.6 mM dNTP).
    /// </summary>
    /// <param name="product">Product (amplicon) sequence.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>Product Tm in °C.</returns>
    /// <exception cref="ArgumentException">Empty product.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative concentrations or zero total cation.</exception>
    public static double CalculateProductMeltingTemperaturePrimer3(
        string product,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar)
    {
        if (string.IsNullOrEmpty(product))
            throw new ArgumentException("Product cannot be null or empty.", nameof(product));
        if (!(monovalentMillimolar >= 0) || !(divalentMillimolar >= 0) || !(dntpMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), "Concentrations must be ≥ 0 mM.");
        double monovalentEq = Primer3MonovalentEquivalent(monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        if (!(monovalentEq > 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), "Total monovalent-equivalent cation concentration must be > 0 mM.");
        int gc = 0;
        foreach (char c in product)
            if (c is 'G' or 'C' or 'g' or 'c') gc++;
        return LongSeqTm(gc, product.Length, monovalentEq);
    }

    // long_seq_tm (oligotm.c) with DMSO = formamide = 0; the canonical salt-adjusted GC formula.
    private static double LongSeqTm(int gc, int length, double monovalentEqMillimolar) =>
        ThermoConstants.CalculateSaltAdjustedTm((double)gc / length, length, monovalentEqMillimolar / 1000.0);

    // Primer3 pair search (make_detection_primer_lists + choose_pair_or_triple + characterize_pair + obj_fn).
    private sealed class PrimerPairSearch
    {
        private readonly int _targetStart, _targetEnd, _incStart, _incEnd;
        private readonly PrimerParameters _param;
        private readonly PrimerPairOptions _opt;
        private readonly Primer3PairWeights _w;
        private readonly List<(PrimerCandidate C, double Tm)> _fwd = new(), _rev = new();
        private readonly bool?[] _fwdOk, _revOk;
        // Primer3 primer_rec.overlaps: excluded by PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE after a pair was selected.
        private readonly bool[] _fwdUsed, _revUsed;
        private readonly Dictionary<string, bool> _structureBySequence = new(StringComparer.Ordinal);
        private readonly Dictionary<(string F, string R), (bool Fails, double? Any, double? End)> _dimer = new();
        private readonly int[] _gcPrefix;
        // long_seq_tm salt for the product Tm: Primer3 passes the primer conditions (p_args) — "skewed" by its own comment.
        private readonly double _productMonovalentEq;
        private readonly ProbeDesigner.Primer3ProbeSettings? _intlSettings;
        private readonly List<ProbeDesigner.Primer3Probe>? _intl;
        private readonly ProbeDesigner.Primer3Probe?[]? _intlChecked;
        private readonly bool?[]? _intlOk;
        private bool _sawCharacterized, _sawTm, _sawDimer, _sawProductTm, _sawInternal, _sawLibrary;
        // PRIMER_MISPRIMING_LIBRARY (null when absent / empty) and each candidate's repeat_sim scores.
        private readonly PrimerMisprimingLibrary? _library;
        private readonly Dictionary<PrimerCandidate, LibraryMispriming> _libraryByCandidate = new(ReferenceEqualityComparer.Instance);
        // Template mispriming (null when Primer3 never computes it: no limit / weight set in the active mode).
        private readonly TemplateContext? _template;
        private readonly bool _needPairTemplate;
        private readonly Dictionary<PrimerCandidate, TemplateMisprimingScore> _templateByCandidate = new(ReferenceEqualityComparer.Instance);
        private bool _sawTemplate;
        // Non-default PRIMER_INSIDE/OUTSIDE_PENALTY: primers are scored against the target (compute_position_penalty)
        // instead of being kept off it; null under Primer3's defaults.
        private readonly TargetPosition? _position;
        // SEQUENCE_QUALITY (null when absent).
        private readonly QualityContext? _quality;
        // PRIMER_MASK_TEMPLATE: the masked copies of the included region (null without masking).
        private readonly MaskContext? _mask;
        // PRIMER_LOWERCASE_MASKING (forced on by PRIMER_MASK_TEMPLATE): the case-preserving template (trimmed_orig_seq);
        // null when lower-case masking is off or the template has no case information.
        private readonly string? _lowercase;

        // seq: the upper-case template; caseTemplate: the same template as given (null = upper case only).
        public PrimerPairSearch(string seq, string? caseTemplate, int targetStart, int targetEnd,
            PrimerParameters param, PrimerPairOptions opt)
        {
            ArgumentNullException.ThrowIfNull(opt);
            int n = seq.Length;
            if (targetStart < 0 || targetEnd >= n || targetStart >= targetEnd)
                throw new ArgumentException("Invalid target region.");
            _targetStart = targetStart;
            _targetEnd = targetEnd;
            _param = param;
            _opt = opt;
            _w = opt.Weights ?? throw new ArgumentException("Pair weights cannot be null.", nameof(opt));
            _library = param.ActiveLibrary;
            ValidateOptions(n);
            _quality = QualityContext.Create(opt.SequenceQuality, param.EffectiveQualityRangeMax);
            // _pr_need_pair_template_mispriming[_thermod] / _pr_need_template_mispriming[_thermod].
            bool thermoTemplate = param.ThermodynamicTemplateAlignment;
            _needPairTemplate = thermoTemplate
                ? opt.MaxTemplateMisprimingTh >= 0 || _w.TemplateMisprimingTh > 0
                : opt.MaxTemplateMispriming >= 0 || _w.TemplateMispriming > 0;
            if (_needPairTemplate || param.ActiveMaxTemplateMispriming >= 0 || param.ActiveTemplateMisprimingWeight > 0)
                _template = new TemplateContext(seq, thermoTemplate, param.EffectiveMonovalentMillimolar,
                    param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar, param.EffectiveDnaConcentrationNanomolar);
            _productMonovalentEq = Primer3MonovalentEquivalent(param.EffectiveMonovalentMillimolar,
                param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar);
            (_incStart, _incEnd) = opt.IncludedRegion is { } inc ? (inc.Start, inc.Start + inc.Length) : (0, n);
            // primer3-py sets lowercase_masking = mask_template when masking (thermoanalysis.pyx), so PRIMER_MASK_TEMPLATE implies
            // PRIMER_LOWERCASE_MASKING.
            if ((opt.LowercaseMasking || opt.MaskTemplate) && caseTemplate is not null)
                _lowercase = caseTemplate;
            // Primer3 masks the included region (trimmed_orig_seq, original case: input lower case stays masked) once per
            // design, before the primer lists are built.
            if (opt.MaskTemplate)
                _mask = new MaskContext(_lowercase ?? seq, _incStart, _incEnd, opt);

            _gcPrefix = new int[n + 1];
            for (int i = 0; i < n; i++)
                _gcPrefix[i + 1] = _gcPrefix[i] + (seq[i] is 'G' or 'C' ? 1 : 0);

            int minProduct = int.MaxValue, maxProduct = 0;
            foreach (var r in opt.ProductSizeRanges)
            {
                minProduct = Math.Min(minProduct, r.Min);
                maxProduct = Math.Max(maxProduct, r.Max);
            }

            // Default position penalties: a primer may not overlap the target (oligo_overlaps_interval). Otherwise
            // (compute_position_penalty) a left primer's 3' end may lie anywhere up to the target's last base and a
            // right primer's from the target's first base on (make_detection_primer_lists then searches the whole
            // included region); the pair must still span the target (pair_spans_target).
            if (!opt.DefaultPositionPenalties)
                _position = new TargetPosition(targetStart, targetEnd, opt.InsidePenalty, opt.OutsidePenalty);
            int leftLimit = _position is null ? targetStart : targetEnd;   // left primer end (exclusive) ≤ leftLimit
            int rightLimit = _position is null ? targetEnd : targetStart;  // right primer start ≥ rightLimit

            // Forward candidates: 3' end before the target; Primer3 drops a left primer that starts past
            // n − min product (pick_primer_range); starts that cannot reach any product range are skipped.
            int fStartLo = Math.Max(_incStart, rightLimit + param.MinLength - maxProduct);
            int fStartHi = Math.Min(leftLimit - param.MinLength, _incEnd - minProduct);
            for (int start = fStartLo; start <= fStartHi; start++)
            {
                for (int len = param.MinLength; len <= param.MaxLength && start + len <= leftLimit; len++)
                {
                    var (candidate, tm, lib, tmp) = EvaluatePrimerCore(seq.Substring(start, len), start, true, param, evaluateStructure: false, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase);
                    if (candidate.IsValid)
                    {
                        _fwd.Add((candidate, tm));
                        if (lib is not null) _libraryByCandidate[candidate] = lib;
                        if (tmp is { } t) _templateByCandidate[candidate] = t;
                    }
                }
            }

            // Reverse candidates (evaluated as the reverse complement): 5' end on the top strand = end − 1.
            int rEndLo = Math.Max(rightLimit + param.MinLength, _incStart + minProduct);
            int rEndHi = Math.Min(_incEnd, leftLimit - param.MinLength + maxProduct);
            for (int end = rEndLo; end <= rEndHi; end++)
            {
                for (int len = param.MinLength; len <= param.MaxLength && end - len >= rightLimit; len++)
                {
                    int start = end - len;
                    var revComp = DnaSequence.GetReverseComplementString(seq.Substring(start, len));
                    var (candidate, tm, lib, tmp) = EvaluatePrimerCore(revComp, start, false, param, evaluateStructure: false, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase);
                    if (candidate.IsValid)
                    {
                        _rev.Add((candidate, tm));
                        if (lib is not null) _libraryByCandidate[candidate] = lib;
                        if (tmp is { } t) _templateByCandidate[candidate] = t;
                    }
                }
            }

            // Primer3 examines primers in increasing penalty order (sort_primer_array).
            _fwd.Sort((a, b) => CompareLeft(a.C, b.C));
            _rev.Sort((a, b) => CompareRight(a.C, b.C));
            _fwdOk = new bool?[_fwd.Count];
            _revOk = new bool?[_rev.Count];
            _fwdUsed = new bool[_fwd.Count];
            _revUsed = new bool[_rev.Count];

            if (opt.PickInternalOligo)
            {
                // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT is global in Primer3: the internal oligo follows the primer screen.
                _intlSettings = (opt.InternalOligo ?? new ProbeDesigner.Primer3ProbeSettings()) with
                {
                    ThermodynamicOligoAlignment = param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
                    // PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS and PRIMER_ANNEALING_TEMP are global too.
                    LibraryAmbiguityCodesConsensus = param.EffectiveLibraryAmbiguityCodesConsensus,
                    AnnealingTemperature = param.EffectiveAnnealingTemperature,
                    // PRIMER_QUALITY_RANGE_MIN/MAX are global as well.
                    QualityRangeMin = param.EffectiveQualityRangeMin,
                    QualityRangeMax = param.EffectiveQualityRangeMax,
                };
                // make_internal_oligo_list over the included region; choose_internal_oligo takes the
                // lowest-penalty oligo (first in enumeration order among equals) → stable sort.
                var list = ProbeDesigner.EnumeratePrimer3InternalOligos(seq, _incStart, _incEnd, _intlSettings, screenStructure: false, _quality?.Values, _lowercase);
                _intl = list
                    .Select((p, i) => (p, i))
                    .OrderBy(t => t.p.Penalty).ThenBy(t => t.i)
                    .Select(t => t.p)
                    .ToList();
                _intlChecked = new ProbeDesigner.Primer3Probe?[_intl.Count];
                _intlOk = new bool?[_intl.Count];
            }
        }

        private void ValidateOptions(int n)
        {
            var o = _opt;
            if (o.ProductSizeRanges is null || o.ProductSizeRanges.Count == 0)
                throw new ArgumentException("At least one product size range is required (PRIMER_PRODUCT_SIZE_RANGE).");
            int minProduct = int.MaxValue;
            foreach (var r in o.ProductSizeRanges)
            {
                if (r.Min < 1 || r.Max < r.Min)
                    throw new ArgumentException($"Invalid product size range {r.Min}-{r.Max}: need 1 ≤ min ≤ max.");
                minProduct = Math.Min(minProduct, r.Min);
            }
            _param.ValidateConditions("parameters");
            if (_param.MinLength < 1 || _param.MaxLength < _param.MinLength)
                throw new ArgumentException("Primer sizes must satisfy 1 ≤ MinLength ≤ MaxLength.");
            if (_param.MaxLength > minProduct)
                throw new ArgumentException("PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE (Primer3 _pr_data_control).");
            if (!(o.MaxTmDifference >= 0))
                throw new ArgumentException("MaxTmDifference must be ≥ 0 °C.");
            if (o.NumReturn < 1)
                throw new ArgumentException("PRIMER_NUM_RETURN < 1 (Primer3 _pr_data_control).");
            if (o.MinLeftThreePrimeDistance < -1 || o.MinRightThreePrimeDistance < -1)
                throw new ArgumentException("Minimum 3' distance must be >= -1 (min_*_three_prime_distance) (Primer3 _pr_data_control).");
            // _pr_data_control sequence-quality checks (they precede the objective-function checks; the PRIMER_INTERNAL_*
            // values count whether or not an internal oligo is picked).
            var io = o.InternalOligo ?? new ProbeDesigner.Primer3ProbeSettings();
            ValidatePrimer3Quality(o.SequenceQuality, n, _param.MinQuality, io.MinQuality, _param.EffectiveQualityRangeMin,
                _param.EffectiveQualityRangeMax, _param.PenaltyWeights?.SequenceQuality ?? 0.0, io.WeightSequenceQuality, nameof(o.SequenceQuality));
            if ((_w.ProductTmLt != 0 || _w.ProductTmGt != 0) && o.ProductOptTm is null)
                throw new ArgumentException("Product temperature is part of objective function while optimum temperature is not defined (Primer3 _pr_data_control).");
            if ((_w.ProductSizeLt != 0 || _w.ProductSizeGt != 0) && o.ProductOptSize is null)
                throw new ArgumentException("Product size is part of objective function while optimum size is not defined (Primer3 _pr_data_control).");
            // _pr_data_control: GC weights without PRIMER_[INTERNAL_]OPT_GC_PERCENT (the internal-oligo check applies whether or
            // not an internal oligo is picked).
            _param.ValidateGcOptimum("parameters");
            ProbeDesigner.ValidatePrimer3ProbeGcOptimum(io, nameof(o.InternalOligo));
            // primer3-py: PRIMER_MASK_TEMPLATE needs the k-mer lists (PRIMER_MASK_KMERLIST_PATH).
            if (o.MaskTemplate && o.MaskKmerLists is null)
                throw new ArgumentException("masking template chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified (primer3-py).");
            ValidateMaskSettings(o.MaskFailureRate, o.MaskFivePrimeDirection, o.MaskThreePrimeDirection, nameof(o.MaskFailureRate));
            if (double.IsNaN(_param.PenaltyWeights?.MaskFailureRate ?? 0.0))
                throw new ArgumentException("PRIMER_WT_MASK_FAILURE_RATE must not be NaN.");
            // _pr_data_control: PRIMER_PAIR_WT_IO_PENALTY without PRIMER_PICK_INTERNAL_OLIGO.
            if (_w.InternalOligoPenalty != 0 && !o.PickInternalOligo)
                throw new ArgumentException("Internal oligo quality is part of objective function while internal oligo choice is not required (Primer3 _pr_data_control).");
            if (!(o.MaxComplAny >= 0 && o.MaxComplAny <= short.MaxValue && o.MaxComplEnd >= 0 && o.MaxComplEnd <= short.MaxValue
                  && _param.EffectiveMaxSelfAny >= 0 && _param.EffectiveMaxSelfAny <= short.MaxValue
                  && _param.EffectiveMaxSelfEnd >= 0 && _param.EffectiveMaxSelfEnd <= short.MaxValue))
                throw new ArgumentException("Illegal value for primer complementarity restrictions (Primer3 _pr_data_control: 0 ≤ limit ≤ 32767).");
            if (double.IsNaN(o.MaxLibraryMispriming)
                || (o.MaxLibraryMispriming > short.MaxValue && _param.StructureScreen == PrimerStructureScreen.Primer3Alignment))
                throw new ArgumentException("Value too large at tag PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (Primer3 _pr_data_control).");
            if (_w.LibraryMispriming != 0 && _library is null)
                throw new ArgumentException("Mispriming score is part of objective function, but mispriming library is not defined (Primer3 _pr_data_control).");
            if (double.IsNaN(o.MaxTemplateMispriming) || double.IsNaN(o.MaxTemplateMisprimingTh)
                || (o.MaxTemplateMispriming > short.MaxValue && !_param.ThermodynamicTemplateAlignment))
                throw new ArgumentException("Value too large at tag PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING (Primer3 _pr_data_control).");
            if (_w.PrimerPenalty < 0 || _w.InternalOligoPenalty < 0 || _w.DiffTm < 0 || _w.ComplAnyTh < 0 || _w.ComplEndTh < 0
                || _w.ComplAny < 0 || _w.ComplEnd < 0 || _w.LibraryMispriming < 0
                || !(_w.TemplateMispriming >= 0) || !(_w.TemplateMisprimingTh >= 0)
                || _w.ProductTmLt < 0 || _w.ProductTmGt < 0 || _w.ProductSizeLt < 0 || _w.ProductSizeGt < 0)
                throw new ArgumentException("Pair weights must be ≥ 0.");
            if (o.IncludedRegion is { } inc)
            {
                if (inc.Start < 0 || inc.Length < 1 || inc.Start + inc.Length > n)
                    throw new ArgumentException("Included region outside the template.");
                if (_targetStart < inc.Start || _targetEnd > inc.Start + inc.Length)
                    throw new ArgumentException("TARGET outside of INCLUDED_REGION (Primer3 _check_and_adjust_1_interval).");
            }
            if (o.PickInternalOligo)
            {
                var s = (o.InternalOligo ?? new ProbeDesigner.Primer3ProbeSettings()) with
                {
                    ThermodynamicOligoAlignment = _param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
                    AnnealingTemperature = _param.EffectiveAnnealingTemperature,
                };
                ProbeDesigner.ValidatePrimer3ProbeSettings(s, nameof(o.InternalOligo));
                if (s.MaxSize > minProduct)
                    throw new ArgumentException("PRIMER_INTERNAL_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE (Primer3 _pr_data_control).");
            }
        }

        // Primer3 primer_rec_comp (sort_primer_array): penalty ascending, then Primer3 "start"
        // descending, then shorter first. A left primer's start is its 5' end (Position); a right
        // primer's start is its 5' end on the top strand (Position + Length − 1).
        private static int CompareLeft(PrimerCandidate a, PrimerCandidate b)
        {
            int c = a.Penalty.CompareTo(b.Penalty);
            if (c != 0) return c;
            c = b.Position.CompareTo(a.Position);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        }

        private static int CompareRight(PrimerCandidate a, PrimerCandidate b)
        {
            int c = a.Penalty.CompareTo(b.Penalty);
            if (c != 0) return c;
            c = (b.Position + b.Length).CompareTo(a.Position + a.Length);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        }

        // Primer3 compare_primer_pair (libprimer3.cc): quality (±1e-6), then left start descending,
        // right start (its 5' end) ascending, left length ascending, right length ascending.
        private static int ComparePair(
            double q1, PrimerCandidate l1, PrimerCandidate r1,
            double q2, PrimerCandidate l2, PrimerCandidate r2)
        {
            if (q1 + PairQualityEpsilon < q2) return -1;
            if (q1 > q2 + PairQualityEpsilon) return 1;
            int c = l2.Position.CompareTo(l1.Position);
            if (c != 0) return c;
            c = (r1.Position + r1.Length).CompareTo(r2.Position + r2.Length);
            if (c != 0) return c;
            c = l1.Length.CompareTo(l2.Length);
            return c != 0 ? c : r1.Length.CompareTo(r2.Length);
        }

        // Pair-level complementarity screen used by DesignPrimers (Primer3 characterize_pair): whether the
        // pair fails, and (thermodynamic screen) PRIMER_PAIR_COMPL_ANY_TH / _COMPL_END_TH.
        private static (bool Fails, double? Any, double? End) PairScreen(
            string forward, string reverse, PrimerParameters param, PrimerPairOptions opt)
        {
            if (param.StructureScreen == PrimerStructureScreen.Heuristic)
                return (HasPrimerDimer(forward, reverse), null, null);
            if (param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
                return AlignmentPairScreen(forward, reverse, param, opt);
            if (!IsAcgtOnly(forward) || !IsAcgtOnly(reverse))
                return (false, null, null);
            // Same values as CalculatePrimer3PairComplementarity(...), stopping at the first alignment over the
            // limit (characterize_pair also fails the pair on compl_any first); complete when the pair passes.
            double max = param.EffectiveMaxStructureTm;
            // characterize_pair uses the primer conditions (thal_arg_to_use = create_thal_arg_holder(p_args)).
            var (any, end) = Primer3PairTms(forward.ToUpperInvariant(), reverse.ToUpperInvariant(),
                param.EffectiveMonovalentMillimolar / 1000.0, param.EffectiveDivalentMillimolar / 1000.0,
                param.EffectiveDntpMillimolar / 1000.0, param.EffectiveDnaConcentrationNanomolar * 1e-9, max);
            return (any > max || end > max, any, end);
        }

        // characterize_pair, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0 (s1 = left, s2 = revcomp(right),
        // s1_rev = revcomp(left), s2_rev = right): compl_any = align(s1, s2, LOCAL) > PRIMER_PAIR_MAX_COMPL_ANY fails;
        // compl_end = align(s1, s2, GLOBAL_END) > PRIMER_PAIR_MAX_COMPL_END fails; then, when
        // align(s2_rev, s1_rev, GLOBAL_END) is larger, it fails above PRIMER_MAX_SELF_END (Primer3 uses the
        // per-primer limit there) and becomes compl_end.
        private static (bool Fails, double? Any, double? End) AlignmentPairScreen(
            string forward, string reverse, PrimerParameters param, PrimerPairOptions opt)
        {
            string l = forward.ToUpperInvariant(), r = reverse.ToUpperInvariant();
            string rcL = ReverseComplementPrimer3(l), rcR = ReverseComplementPrimer3(r);
            double any = DpalLocalScore(l, rcR);
            if (any > opt.MaxComplAny)
                return (true, any, null);
            double end = DpalGlobalEndScore(l, rcR);
            if (end > opt.MaxComplEnd)
                return (true, any, end);
            double end2 = DpalGlobalEndScore(r, rcL);
            if (end2 > end)
            {
                if (end2 > param.EffectiveMaxSelfEnd)
                    return (true, any, end2);
                end = end2;
            }
            return (false, any, end);
        }

        private sealed record PairEval(
            int Fi, int Ri, double Penalty, int ProductSize, double ProductTm,
            double? ComplAny, double? ComplEnd, ProbeDesigner.Primer3Probe? Internal,
            int? LibraryScore = null, string? LibraryName = null, double? TemplateScore = null);

        public List<PrimerPairResult> Run(int numReturn)
        {
            var results = new List<PrimerPairResult>();
            if (_fwd.Count == 0 || _rev.Count == 0)
                return results;
            // PRIMER_PICK_INTERNAL_OLIGO with no acceptable internal oligo: make_internal_oligo_list fails and Primer3 returns
            // before examining any pair.
            if (_intl is { Count: 0 })
            {
                _sawCharacterized = _sawInternal = true;
                return results;
            }

            var ranges = _opt.ProductSizeRanges;
            var cache = new Dictionary<(int Ri, int Fi), PairEval?>();
            double wq = _w.PrimerPenalty;
            int rangeIndex = 0;
            while (true)
            {
                PairEval? best = null;
                for (int i = 0; i < _rev.Count; i++)
                {
                    // Only primers that are (still) legal; their expensive checks run in characterize_pair.
                    if (_revOk[i] == false)
                        continue;
                    // No pair with this or any later reverse primer can beat the best pair.
                    if (wq * (_rev[i].C.Penalty + _fwd[0].C.Penalty) > BestQuality(best))
                        break;
                    if (_revUsed[i])
                        continue; // 3′ end too close to a right primer of a selected pair

                    for (int j = 0; j < _fwd.Count; j++)
                    {
                        if (_revOk[i] == false)
                            break;
                        if (_fwdOk[j] == false)
                            continue;
                        if (wq * (_fwd[j].C.Penalty + _rev[i].C.Penalty) > BestQuality(best))
                            break;
                        if (_fwdUsed[j])
                            continue;

                        int product = ProductSize(j, i);
                        var range = ranges[rangeIndex];
                        if (product < range.Min || product > range.Max)
                            continue;

                        if (!cache.TryGetValue((i, j), out var e))
                        {
                            e = Characterize(j, i, product);
                            cache[(i, j)] = e;
                        }
                        if (e is null)
                            continue; // illegal, or already selected

                        if (best is null || ComparePair(e.Penalty, _fwd[j].C, _rev[i].C,
                                best.Penalty, _fwd[best.Fi].C, _rev[best.Ri].C) < 0)
                            best = e;
                        if (best.Penalty == 0)
                            break; // there cannot be a better pair
                    }
                    if (best is { Penalty: 0 })
                        break;
                }

                if (best is null)
                {
                    // No pair in this product-size range: try the next one, if any.
                    if (++rangeIndex >= ranges.Count)
                        break;
                    continue;
                }

                results.Add(ToResult(best));
                cache[(best.Ri, best.Fi)] = null; // mark as selected
                MarkUsedPrimers(best);
                if (results.Count == numReturn)
                    break;
            }
            return results;
        }

        private static double BestQuality(PairEval? best) => best?.Penalty ?? double.MaxValue;

        // Primer3 left/right_oligo_in_pair_overlaps_used_oligo after a pair is selected: min distance −1 → nothing;
        // 0 → the identical primer (same start and length); d > 0 → every primer whose 3′ end is < d bases from the
        // selected primer's 3′ end (left 3′ end = Position + Length − 1; right 3′ end on the top strand = Position).
        private void MarkUsedPrimers(PairEval best)
        {
            MarkUsed(_fwd, _fwdUsed, _fwd[best.Fi].C, _opt.MinLeftThreePrimeDistance, c => c.Position + c.Length - 1);
            MarkUsed(_rev, _revUsed, _rev[best.Ri].C, _opt.MinRightThreePrimeDistance, c => c.Position);

            static void MarkUsed(List<(PrimerCandidate C, double Tm)> list, bool[] used, PrimerCandidate chosen,
                int minDistance, Func<PrimerCandidate, int> threePrime)
            {
                if (minDistance == -1)
                    return;
                int chosenEnd = threePrime(chosen);
                for (int k = 0; k < list.Count; k++)
                {
                    var c = list[k].C;
                    if (minDistance == 0
                        ? c.Position == chosen.Position && c.Length == chosen.Length
                        : Math.Abs(threePrime(c) - chosenEnd) < minDistance)
                        used[k] = true;
                }
            }
        }

        private int ProductSize(int fi, int ri) =>
            _rev[ri].C.Position + _rev[ri].C.Length - _fwd[fi].C.Position;

        // characterize_pair (+ choose_internal_oligo and obj_fn for a legal pair).
        private PairEval? Characterize(int fi, int ri, int product)
        {
            _sawCharacterized = true;
            var f = _fwd[fi];
            var r = _rev[ri];

            // pair_spans_target: the left primer's 3' end must precede the right primer's (always true when neither
            // may overlap the target).
            if (_position is not null && f.C.Position + f.C.Length - 1 >= r.C.Position)
                return null;

            int fStart = f.C.Position;
            int gc = _gcPrefix[fStart + product] - _gcPrefix[fStart];
            double productTm = LongSeqTm(gc, product, _productMonovalentEq);
            if ((_opt.ProductMinTm is { } minTm && productTm < minTm)
                || (_opt.ProductMaxTm is { } maxTm && productTm > maxTm))
            {
                _sawProductTm = true;
                return null;
            }

            double diffTm = Math.Abs(f.Tm - r.Tm);
            if (diffTm > _opt.MaxTmDifference)
            {
                _sawTm = true;
                return null;
            }

            if (!StructureOk(_fwd, _fwdOk, fi) || !StructureOk(_rev, _revOk, ri))
                return null;

            var (fails, any, end) = PairComplementarity(f.C.Sequence, r.C.Sequence);
            if (fails)
            {
                _sawDimer = true;
                return null;
            }

            // pair_repeat_sim > PRIMER_PAIR_MAX_LIBRARY_MISPRIMING fails the pair.
            int? libScore = null;
            string? libName = null;
            if (_library is not null)
            {
                (int score, libName) = PairLibraryMispriming(LibraryOf(f.C), LibraryOf(r.C), _library);
                if (score > _opt.MaxLibraryMispriming)
                {
                    _sawLibrary = true;
                    return null;
                }
                libScore = score;
            }

            // Pair template mispriming (end of characterize_pair): max(left.T + right.T_r, left.T_r + right.T); fails
            // above PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING when that is ≥ 0, or (thermodynamic mode) above
            // PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH when that is non-zero.
            double? templateScore = null;
            if (_needPairTemplate)
            {
                if (!_templateByCandidate.TryGetValue(f.C, out var lt) || !_templateByCandidate.TryGetValue(r.C, out var rt))
                    throw new InvalidOperationException(
                        "Primer3 PR_ASSERT(template_mispriming != ALIGN_SCORE_UNDEF): with a mispriming library scored at pick time (PRIMER_WT_LIBRARY_MISPRIMING ≠ 0) Primer3 never scores template mispriming for the pair unless PRIMER_WT_TEMPLATE_MISPRIMING[_TH] ≠ 0.");
                double v = lt.SameStrand + rt.OtherStrand;
                if (lt.OtherStrand + rt.SameStrand > v)
                    v = lt.OtherStrand + rt.SameStrand;
                bool templateFails = _param.ThermodynamicTemplateAlignment
                    ? _opt.MaxTemplateMisprimingTh != 0 && v > _opt.MaxTemplateMisprimingTh
                    : _opt.MaxTemplateMispriming >= 0 && v > _opt.MaxTemplateMispriming;
                if (templateFails)
                {
                    _sawTemplate = true;
                    return null;
                }
                templateScore = v;
            }

            ProbeDesigner.Primer3Probe? intl = null;
            if (_intl is not null)
            {
                intl = ChooseInternalOligo(f.C, r.C);
                if (intl is null)
                {
                    _sawInternal = true;
                    return null;
                }
            }

            double penalty = ObjectiveFunction(f, r, diffTm, any, end, productTm, product, intl, libScore ?? 0, templateScore ?? 0);
            return new PairEval(fi, ri, penalty, product, productTm, any, end, intl, libScore, libName, templateScore);
        }

        // The candidate's repeat_sim scores (computed at pick time when weighted, otherwise now).
        private LibraryMispriming LibraryOf(PrimerCandidate c)
        {
            if (!_libraryByCandidate.TryGetValue(c, out var lib))
            {
                lib = ComputeLibraryMispriming(c.Sequence, c.IsForward, _library!,
                    _param.EffectiveLibraryAmbiguityCodesConsensus, _param.EffectiveMaxLibraryMispriming);
                _libraryByCandidate[c] = lib;
            }
            return lib;
        }

        // obj_fn (libprimer3.cc): thermodynamic mode (ComplAnyTh/ComplEndTh terms), alignment mode (linear
        // ComplAny/ComplEnd terms) or, under the heuristic screen, without complementarity terms.
        private double ObjectiveFunction(
            (PrimerCandidate C, double Tm) f, (PrimerCandidate C, double Tm) r, double diffTm,
            double? complAny, double? complEnd, double productTm, int product, ProbeDesigner.Primer3Probe? intl,
            int libraryScore, double templateScore)
        {
            double sum = 0.0;
            double lowerTm = r.Tm;
            if (f.Tm < r.Tm) lowerTm = f.Tm;

            if (_w.PrimerPenalty != 0)
                sum += _w.PrimerPenalty * (f.C.Penalty + r.C.Penalty);
            if (_w.InternalOligoPenalty != 0 && intl is { } io)
                sum += _w.InternalOligoPenalty * io.Penalty;
            if (_w.DiffTm != 0)
                sum += _w.DiffTm * diffTm;
            if (_param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
            {
                // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0: linear PRIMER_PAIR_WT_COMPL_ANY / _COMPL_END terms.
                if (_w.ComplAny != 0)
                    sum += _w.ComplAny * complAny!.Value;
                if (_w.ComplEnd != 0)
                    sum += _w.ComplEnd * complEnd!.Value;
            }
            else
            {
                if (complAny is { } a)
                    sum += ThermodynamicStructurePenalty(_w.ComplAnyTh, lowerTm, a);
                if (complEnd is { } e)
                    sum += ThermodynamicStructurePenalty(_w.ComplEndTh, lowerTm, e);
            }
            if (_w.ProductTmLt != 0 && productTm < _opt.ProductOptTm!.Value)
                sum += _w.ProductTmLt * (_opt.ProductOptTm.Value - productTm);
            if (_w.ProductTmGt != 0 && productTm > _opt.ProductOptTm!.Value)
                sum += _w.ProductTmGt * (productTm - _opt.ProductOptTm.Value);
            if (_w.ProductSizeLt != 0 && product < _opt.ProductOptSize!.Value)
                sum += _w.ProductSizeLt * (_opt.ProductOptSize.Value - product);
            if (_w.ProductSizeGt != 0 && product > _opt.ProductOptSize!.Value)
                sum += _w.ProductSizeGt * (product - _opt.ProductOptSize.Value);
            if (_w.LibraryMispriming != 0)
                sum += _w.LibraryMispriming * libraryScore;
            if (!_param.ThermodynamicTemplateAlignment)
            {
                if (_w.TemplateMispriming != 0)
                    sum += _w.TemplateMispriming * templateScore;
            }
            else
                sum += ThermodynamicStructurePenalty(_w.TemplateMisprimingTh, lowerTm, templateScore);
            // obj_fn ends with PR_ASSERT(sum >= 0.0): Primer3 aborts on a negative pair penalty, which only negative
            // per-primer terms produce (e.g. the default PRIMER_INSIDE_PENALTY −1 with a changed PRIMER_OUTSIDE_PENALTY).
            if (sum < 0.0)
                throw new InvalidOperationException(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Negative primer-pair penalty {sum} (Primer3 obj_fn PR_ASSERT(sum >= 0.0) aborts): negative per-primer penalties — e.g. PRIMER_INSIDE_PENALTY < 0 (default −1) with a non-default PRIMER_OUTSIDE_PENALTY — are not usable."));
            return sum;
        }

        // choose_internal_oligo: the lowest-penalty internal oligo lying strictly between the primers
        // (start after the left primer's 3' end, end before the right primer's 5'-most top-strand base)
        // whose postponed self-any / self-end / hairpin checks pass; a failing oligo stays rejected.
        private ProbeDesigner.Primer3Probe? ChooseInternalOligo(PrimerCandidate left, PrimerCandidate right)
        {
            int leftEnd = left.Position + left.Length - 1;
            for (int k = 0; k < _intl!.Count; k++)
            {
                var h = _intl[k];
                if (h.Start <= leftEnd || h.Start + h.Length - 1 >= right.Position)
                    continue;
                if (_intlOk![k] is null)
                {
                    _intlChecked![k] = ProbeDesigner.PassesPrimer3ProbeStructure(h, _intlSettings!);
                    _intlOk[k] = _intlChecked[k] is not null;
                }
                if (_intlOk[k] == true)
                    return _intlChecked![k];
            }
            return null;
        }

        private bool StructureOk(List<(PrimerCandidate C, double Tm)> list, bool?[] cache, int i)
        {
            if (cache[i] is { } known)
                return known;
            string seq = list[i].C.Sequence;
            if (!_structureBySequence.TryGetValue(seq, out bool ok))
            {
                var issues = new List<string>();
                AddStructureIssues(seq, _param, issues);
                ok = issues.Count == 0;
                _structureBySequence[seq] = ok;
            }
            // characterize_pair then checks the primer against the mispriming library (oligo_repeat_library_mispriming).
            if (ok && _library is not null && LibraryOf(list[i].C).Exceeds)
            {
                ok = false;
                _sawLibrary = true;
            }
            // ... and, when its library scores were not computed at pick time (repeat_sim.score == NULL: no library, or
            // PRIMER_WT_LIBRARY_MISPRIMING = 0), against the template (oligo_template_mispriming).
            if (ok && _template is not null
                && (_library is null || EffectivePenaltyWeights(_param).LibraryMispriming == 0)
                && TemplateMisprimingExceeds(TemplateOf(list[i].C), _param))
            {
                ok = false;
                _sawTemplate = true;
            }
            cache[i] = ok;
            return ok;
        }

        private (bool Fails, double? Any, double? End) PairComplementarity(string f, string r)
        {
            if (!_dimer.TryGetValue((f, r), out var v))
            {
                v = PairScreen(f, r, _param, _opt);
                _dimer[(f, r)] = v;
            }
            return v;
        }

        // The candidate's template mispriming scores (computed at pick time when weighted, otherwise now).
        private TemplateMisprimingScore TemplateOf(PrimerCandidate c)
        {
            if (!_templateByCandidate.TryGetValue(c, out var t))
            {
                t = _template!.Score(c.Sequence, c.Position, c.IsForward);
                _templateByCandidate[c] = t;
            }
            return t;
        }

        private PrimerPairResult ToResult(PairEval e)
        {
            var forward = WithTemplate(Reevaluate(_fwd[e.Fi].C), _fwd[e.Fi].C);
            var reverse = WithTemplate(Reevaluate(_rev[e.Ri].C), _rev[e.Ri].C);
            bool alignment = _param.StructureScreen == PrimerStructureScreen.Primer3Alignment;
            return new PrimerPairResult(
                Forward: forward,
                Reverse: reverse,
                IsValid: true,
                Message: "Valid primer pair found.",
                ProductSize: e.ProductSize,
                PairPenalty: e.Penalty,
                ProductTm: e.ProductTm,
                ComplAnyTh: alignment ? null : e.ComplAny,
                ComplEndTh: alignment ? null : e.ComplEnd,
                InternalOligo: e.Internal)
            {
                ComplAny = alignment ? e.ComplAny : null,
                ComplEnd = alignment ? e.ComplEnd : null,
                LibraryMispriming = e.LibraryScore,
                LibraryMisprimingName = e.LibraryName,
                TemplateMispriming = e.TemplateScore,
            };
        }

        private PrimerCandidate WithTemplate(PrimerCandidate full, PrimerCandidate picked) =>
            _templateByCandidate.TryGetValue(picked, out var t) ? full with { TemplateMispriming = t.Max } : full;

        // Full evaluation (including the structure values) of a chosen primer.
        private PrimerCandidate Reevaluate(PrimerCandidate c) =>
            EvaluatePrimerCore(c.Sequence, c.Position, c.IsForward, _param, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase).Candidate;

        // Result when no pair qualifies: the individually lowest-penalty primers that pass their own
        // (structure) constraints, with the violated pair constraint.
        public PrimerPairResult Failure()
        {
            const string noPrimers = "Could not find valid primers for the target region.";
            string noProduct = $"No primer pair with a product size in {string.Join(", ", _opt.ProductSizeRanges.Select(x => $"{x.Min}-{x.Max}"))} bp (PRIMER_PRODUCT_SIZE_RANGE).";
            if (_fwd.Count == 0 || _rev.Count == 0)
            {
                // A search region shorter than the smallest product cannot hold any pair (Primer3
                // _pr_data_control per-sequence error).
                int minProduct = _opt.ProductSizeRanges.Min(x => x.Min);
                return new PrimerPairResult(null, null, false,
                    _incEnd - _incStart < minProduct ? "SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE" : noPrimers, 0);
            }
            int f0i = FirstStructurallyValid(_fwd, _fwdOk);
            int r0i = FirstStructurallyValid(_rev, _revOk);
            if (f0i < 0 || r0i < 0)
                return new PrimerPairResult(null, null, false, noPrimers, 0);

            var f0 = Reevaluate(_fwd[f0i].C);
            var r0 = Reevaluate(_rev[r0i].C);
            double maxDiff = _opt.MaxTmDifference;
            string reason;
            if (!_sawCharacterized)
                reason = noProduct;
            else if (_sawProductTm || _sawInternal || _sawLibrary || _sawTemplate)
            {
                var parts = new List<string>();
                if (_sawLibrary) parts.Add("the mispriming-library limits");
                if (_sawTemplate) parts.Add("the template-mispriming limits");
                if (_sawProductTm) parts.Add("product Tm limits");
                if (_sawTm) parts.Add($"the {maxDiff:0.##}°C Tm-difference limit");
                if (_sawDimer) parts.Add("primer-dimer avoidance");
                if (_sawInternal) parts.Add("an acceptable internal oligo between the primers");
                reason = $"No primer pair satisfies {string.Join(", ", parts)}.";
            }
            else if (_sawTm && !_sawDimer)
                reason = $"No primer pair within the {maxDiff:0.##}°C Tm-difference limit (best primers: Tm {f0.MeltingTemperature:F1}/{r0.MeltingTemperature:F1}°C).";
            else if (_sawDimer && !_sawTm)
                reason = "Every candidate primer pair forms a primer-dimer.";
            else
                reason = $"No primer pair satisfies both the {maxDiff:0.##}°C Tm-difference limit and primer-dimer avoidance.";
            return new PrimerPairResult(
                Forward: f0,
                Reverse: r0,
                IsValid: false,
                Message: reason,
                ProductSize: r0.Position + r0.Length - f0.Position);

            int FirstStructurallyValid(List<(PrimerCandidate C, double Tm)> list, bool?[] cache)
            {
                for (int i = 0; i < list.Count; i++)
                    if (StructureOk(list, cache, i))
                        return i;
                return -1;
            }
        }
    }

    // Primer3 compare_primer_pair epsilon on pair_quality.
    private const double PairQualityEpsilon = 1e-6;

    /// <summary>
    /// Evaluates a single primer candidate against the per-primer constraints of
    /// <paramref name="parameters"/> (length, GC%, Tm, homopolymer, dinucleotide repeat, secondary
    /// structure, Primer3's 3′-end checks PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY
    /// (<see cref="PrimerParameters.GcClamp"/>, <see cref="PrimerParameters.MaxEndGc"/>,
    /// <see cref="PrimerParameters.MaxEndStability"/>), the deprecated library rule <c>Avoid3PrimeGC</c>). The secondary-structure screen is
    /// <see cref="PrimerParameters.StructureScreen"/>: by default Primer3's thermodynamic limits
    /// (ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ 47 °C, reported in
    /// <see cref="PrimerCandidate.SelfAnyTh"/>/<see cref="PrimerCandidate.SelfEndTh"/>/<see cref="PrimerCandidate.HairpinTh"/>),
    /// Primer3's alignment-mode limits (<see cref="PrimerStructureScreen.Primer3Alignment"/>: dpal self_any ≤ 8,
    /// self_end ≤ 3, reported in <see cref="PrimerCandidate.SelfAny"/>/<see cref="PrimerCandidate.SelfEnd"/>),
    /// or the sequence-only <see cref="HasHairpinPotential"/>. The Tm is Primer3's primer Tm
    /// (<see cref="CalculateMeltingTemperaturePrimer3"/>: SantaLucia 1998 nearest-neighbour,
    /// SantaLucia salt correction) at the reaction conditions of <paramref name="parameters"/>
    /// (<see cref="PrimerParameters.MonovalentMillimolar"/>, <see cref="PrimerParameters.DivalentMillimolar"/>,
    /// <see cref="PrimerParameters.DntpMillimolar"/>, <see cref="PrimerParameters.DnaConcentrationNanomolar"/>; Primer3's
    /// defaults 50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo), the scale on which the Primer3-sourced Tm
    /// window 57–63 °C (opt 60) is defined; the ntthal structure values use the same conditions (Primer3
    /// <c>create_thal_arg_holder(p_args)</c>). A sequence
    /// containing a non-ACGT base has no computable Tm (Primer3 PRIMER_MAX_NS_ACCEPTED = 0): it is
    /// reported with Tm 0 and an issue. <see cref="PrimerCandidate.Penalty"/> is the Primer3
    /// per-primer penalty (<see cref="CalculatePrimer3Penalty"/>, <see cref="PrimerParameters.PenaltyWeights"/>, optima
    /// OptimalTm / OptimalLength / <see cref="PrimerParameters.OptimalGcPercent"/>) used by <see cref="DesignPrimers"/> for ranking;
    /// <see cref="PrimerCandidate.Score"/> is an informational additive quality score (0–100,
    /// higher is better) that does not drive selection.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">Illegal reaction conditions (Primer3 <c>_pr_data_control</c>:
    /// salt or DNA concentration ≤ 0, negative divalent / dNTP concentration, NaN / ∞).</exception>
    /// <exception cref="ArgumentException">PRIMER_MIN_QUALITY / a quality weight without a template, or a GC weight without
    /// <see cref="PrimerParameters.OptimalGcPercent"/> (Primer3 <c>_pr_data_control</c>).</exception>
    public static PrimerCandidate EvaluatePrimer(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;
        param.ValidateConditions(nameof(parameters));
        // No template, hence no SEQUENCE_QUALITY: Primer3 rejects PRIMER_MIN_QUALITY / PRIMER_WT_SEQ_QUAL without it.
        ValidatePrimer3Quality(null, 0, param.MinQuality, 0, param.EffectiveQualityRangeMin, param.EffectiveQualityRangeMax,
            param.PenaltyWeights?.SequenceQuality ?? 0.0, 0.0, nameof(parameters));
        param.ValidateGcOptimum(nameof(parameters));
        return EvaluatePrimerCore(sequence, position, isForward, param).Candidate;
    }

    // Primer3 is_lowercase_masked: only a / c / g / t count as masked.
    internal static bool IsLowercaseMaskedBase(char c) => c is 'a' or 'c' or 'g' or 't';

    // Evaluates a candidate and also returns its unrounded Tm (Primer3 compares unrounded Tm values).
    // With evaluateStructure = false the secondary-structure screen is skipped (DesignPrimers runs it
    // lazily, like Primer3's characterize_pair, and re-evaluates the chosen primers in full).
    // template (pair search only): Primer3 scores template mispriming at pick time when the active template weight is
    // non-zero (calc_and_check_oligo_features), for an oligo that passed every earlier check.
    private static (PrimerCandidate Candidate, double Tm, LibraryMispriming? Library, TemplateMisprimingScore? Template) EvaluatePrimerCore(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters param,
        bool evaluateStructure = true,
        TemplateContext? template = null,
        TargetPosition? target = null,
        QualityContext? quality = null,
        MaskContext? mask = null,
        string? lowercase = null)
    {
        var seq = sequence.ToUpperInvariant();

        double gcContent = CalculateGcContent(seq);
        // seqtm (MAX_NN_TM_LENGTH 36) with the fraction bound at PRIMER_ANNEALING_TEMP.
        double annealing = param.EffectiveAnnealingTemperature;
        var (tmRaw, boundRaw) = Primer3SeqTm(seq, param.EffectiveDnaConcentrationNanomolar,
            Primer3MonovalentEquivalent(param.EffectiveMonovalentMillimolar, param.EffectiveDivalentMillimolar,
                param.EffectiveDntpMillimolar), annealing);
        bool tmComputable = !double.IsNaN(tmRaw);
        double tm = tmComputable ? tmRaw : 0.0;
        int homopolymer = FindLongestHomopolymer(seq);
        int dinucRepeat = FindLongestDinucleotideRepeat(seq);
        double stability3Prime = Calculate3PrimeStability(seq);

        var issues = new List<string>();

        // Validate against parameters
        if (seq.Length < param.MinLength || seq.Length > param.MaxLength)
            issues.Add($"Length {seq.Length} outside range [{param.MinLength}-{param.MaxLength}]");

        if (gcContent < param.MinGcContent || gcContent > param.MaxGcContent)
            issues.Add($"GC content {gcContent:F1}% outside range [{param.MinGcContent}-{param.MaxGcContent}]%");

        if (!tmComputable)
            issues.Add("Tm not computable: sequence contains a non-ACGT base");
        else if (tm < param.MinTm || tm > param.MaxTm)
            issues.Add($"Tm {tm:F1}°C outside range [{param.MinTm}-{param.MaxTm}]°C");

        // Primer3 fraction bound (calc_and_check_oligo_features, after the Tm checks): only when PRIMER_ANNEALING_TEMP > 0.
        if (annealing > 0.0)
        {
            if (boundRaw < param.EffectiveMinBound)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Fraction bound {BoundText(boundRaw)} below {param.EffectiveMinBound:0.##}% (Primer3 PRIMER_MIN_BOUND)"));
            if (boundRaw > param.EffectiveMaxBound)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Fraction bound {BoundText(boundRaw)} above {param.EffectiveMaxBound:0.##}% (Primer3 PRIMER_MAX_BOUND)"));
        }

        // PRIMER_MASK_TEMPLATE (is_lowercase_masked, among the first checks of calc_and_check_oligo_features): the primer's
        // 3' base must not be masked on its strand's copy of the template.
        // With masking the masked copies (built from the case-preserving template) carry the input's lower case too.
        if (mask is not null && mask.ThreePrimeEndMasked(position, seq.Length, isForward))
            issues.Add("3' end overlaps masked sequence (Primer3 PRIMER_MASK_TEMPLATE)");
        else if (mask is null && lowercase is not null
                 && IsLowercaseMaskedBase(lowercase[isForward ? position + seq.Length - 1 : position]))
            issues.Add("3' end overlaps lower-case masked sequence (Primer3 PRIMER_LOWERCASE_MASKING)");

        // Primer3 position penalty relative to the target (non-default PRIMER_INSIDE/OUTSIDE_PENALTY, pair search only):
        // a 3' end past the target is an infinite position penalty (OP_OVERLAPS_TARGET).
        double? positionPenalty = null;
        if (target is { } tp)
        {
            positionPenalty = CalculatePositionPenaltyPrimer3(position, seq.Length, isForward, tp.Start, tp.End,
                tp.InsidePenalty, tp.OutsidePenalty);
            if (positionPenalty is null)
                issues.Add("3' end beyond the target (Primer3 infinite position penalty, overlaps target)");
        }

        if (homopolymer > param.MaxHomopolymer)
            issues.Add($"Homopolymer run of {homopolymer} exceeds max {param.MaxHomopolymer}");

        if (dinucRepeat > param.MaxDinucleotideRepeats)
            issues.Add($"Dinucleotide repeat of {dinucRepeat} exceeds max {param.MaxDinucleotideRepeats}");

        bool hasHairpin = false;
        StructureValues structure = default;
        if (evaluateStructure)
            (hasHairpin, structure) = AddStructureIssues(seq, param, issues);

        // Primer3 3′-end checks of left/right primers (calc_and_check_oligo_features): PRIMER_GC_CLAMP — the
        // gc_clamp 3′-most bases must all be G/C; PRIMER_MAX_END_GC — checked only when < 5, at most max_end_gc
        // G/C among the five 3′-most bases (a shorter primer counts its own bases); PRIMER_MAX_END_STABILITY —
        // end_stability = end_oligodg(seq, 5) = −ΔG of the 3′ pentamer must not exceed the limit.
        int gcClamp = param.GcClamp;
        for (int i = 0; i < gcClamp; i++)
        {
            int k = seq.Length - 1 - i;
            if (k < 0 || seq[k] is not ('G' or 'C'))
            {
                issues.Add($"No GC clamp: the {gcClamp} 3'-most bases must be G/C (PRIMER_GC_CLAMP)");
                break;
            }
        }

        int maxEndGc = param.EffectiveMaxEndGc;
        if (maxEndGc < Primer3MaxEndGc)
        {
            int endGc = 0;
            for (int k = Math.Max(0, seq.Length - 5); k < seq.Length; k++)
                if (seq[k] is 'G' or 'C') endGc++;
            if (endGc > maxEndGc)
                issues.Add($"{endGc} G/C in the last 5 bases exceeds max {maxEndGc} (PRIMER_MAX_END_GC)");
        }

        // Primer3 sequence_quality_is_ok (after the end-GC check; only with SEQUENCE_QUALITY): the minimum base quality,
        // then the minimum over the five 3'-most bases.
        int? minQuality = null;
        if (quality is not null)
        {
            var (qMin, qEnd) = CalculateSequenceQualityPrimer3(quality.Values, position, seq.Length, isForward, quality.RangeMax);
            minQuality = qMin;
            if (qMin < param.MinQuality)
                issues.Add($"Minimum base quality {qMin} below {param.MinQuality} (Primer3 PRIMER_MIN_QUALITY)");
            else if (qEnd < param.MinEndQuality)
                issues.Add($"Minimum 3'-end base quality {qEnd} below {param.MinEndQuality} (Primer3 PRIMER_MIN_END_QUALITY)");
        }

        double endStability = -stability3Prime;
        if (endStability > param.EffectiveMaxEndStability)
            issues.Add($"3' end stability {endStability:F2} kcal/mol exceeds max {param.EffectiveMaxEndStability} (PRIMER_MAX_END_STABILITY)");

        // Library rule kept for source compatibility (not Primer3; superseded by GcClamp / MaxEndGc).
        if (param.Avoid3PrimeGC && seq.Length >= 2)
        {
            string last2 = seq.Substring(seq.Length - 2);
            int gcCount = last2.Count(c => c == 'G' || c == 'C');
            if (gcCount == 0)
                issues.Add("No GC clamp at 3' end");
        }

        // Primer3 mispriming library (oligo_repeat_library_mispriming): at pick time only when weighted
        // (PRIMER_WT_LIBRARY_MISPRIMING ≠ 0), otherwise in characterize_pair — the pair search does that lazily.
        var weights = EffectivePenaltyWeights(param);
        LibraryMispriming? library = null;
        if (param.ActiveLibrary is { } lib && (evaluateStructure || weights.LibraryMispriming != 0))
        {
            library = ComputeLibraryMispriming(seq, isForward, lib, param.EffectiveLibraryAmbiguityCodesConsensus,
                param.EffectiveMaxLibraryMispriming);
            if (library.Exceeds)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Library mispriming score {library.Scores.Max():0.00} exceeds {param.EffectiveMaxLibraryMispriming:0.00} (Primer3 PRIMER_MAX_LIBRARY_MISPRIMING)"));
        }

        TemplateMisprimingScore? templateScore = null;
        if (template is not null && issues.Count == 0 && param.ActiveTemplateMisprimingWeight != 0)
        {
            templateScore = template.Score(seq, position, isForward);
            if (TemplateMisprimingExceeds(templateScore.Value, param))
                issues.Add(TemplateMisprimingIssue(templateScore.Value, param));
        }
        // Without a template (EvaluatePrimer) the template mispriming score is undefined: no template term
        // (callers combine CalculateTemplateMispriming with CalculatePrimer3Penalty themselves).
        if (templateScore is null)
            weights = weights with { TemplateMispriming = 0.0, TemplateMisprimingTh = 0.0 };

        bool isValid = issues.Count == 0;

        // Masker failure rate (end of calc_and_check_oligo_features, 16-nt window): only with PRIMER_MASK_TEMPLATE.
        double? failureRate = mask is null ? null : CalculateMaskFailureRatePrimer3(seq, mask.Lists);

        // Informational heuristic score and the Primer3 ranking penalty.
        double score = CalculatePrimerScore(seq, gcContent, tm, homopolymer, param);
        if (!evaluateStructure && UsesStructureTerms(weights, param.StructureScreen))
            structure = ComputeStructureValues(seq, param); // Primer3 computes them at pick time when weighted
        double penalty = CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(tm, seq.Length, gcContent,
                SelfAny: (weights.ThermodynamicOligoAlignment ? structure.Thermo?.SelfAnyTh : structure.SelfAny) ?? 0.0,
                SelfEnd: (weights.ThermodynamicOligoAlignment ? structure.Thermo?.SelfEndTh : structure.SelfEnd) ?? 0.0,
                HairpinTh: structure.Thermo?.HairpinTh ?? 0.0,
                // p_obj_fn end_stability term (left/right primers only): h->end_stability =
                // end_oligodg(seq, 5, santalucia) = −ΔG of the 3′ pentamer (positive magnitude).
                EndStability: double.IsNaN(stability3Prime) ? 0.0 : -stability3Prime)
            {
                LibraryMispriming = library?.MaxScore ?? 0.0,
                TemplateMispriming = templateScore?.Max ?? 0.0,
                Bound = annealing > 0.0 ? boundRaw : null,
                PositionPenalty = positionPenalty ?? 0.0,
                SequenceQuality = minQuality,
                QualityRangeMax = param.EffectiveQualityRangeMax,
                MaskFailureRate = failureRate ?? 0.0,
            },
            weights,
            new Primer3Optima(param.OptimalTm, param.OptimalLength, param.OptimalGcPercent) { OptBound = param.EffectiveOptBound });

        var candidate = new PrimerCandidate(
            Sequence: seq,
            Position: position,
            IsForward: isForward,
            Length: seq.Length,
            GcContent: Math.Round(gcContent, 1),
            MeltingTemperature: Math.Round(tm, 1),
            HomopolymerLength: homopolymer,
            HasHairpin: hasHairpin,
            Stability3Prime: Math.Round(stability3Prime, 1),
            IsValid: isValid,
            Issues: issues.AsReadOnly(),
            Score: Math.Round(score, 2),
            Penalty: penalty,
            SelfAnyTh: structure.Thermo?.SelfAnyTh,
            SelfEndTh: structure.Thermo?.SelfEndTh,
            HairpinTh: structure.Thermo?.HairpinTh
        )
        {
            SelfAny = structure.SelfAny,
            SelfEnd = structure.SelfEnd,
            LibraryMispriming = library?.MaxScore,
            LibraryMisprimingName = library?.Name,
            TemplateMispriming = templateScore?.Max,
            Bound = annealing > 0.0 && boundRaw != Primer3OligoTmError ? boundRaw : null,
            PositionPenalty = positionPenalty,
            MinSequenceQuality = minQuality,
            MaskFailureRate = failureRate,
        };
        return (candidate, tm, library, templateScore);

        static string BoundText(double b) => b == Primer3OligoTmError
            ? "undefined (Primer3 OLIGOTM_ERROR: > 36 bases or non-ACGT)"
            : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{b:0.##}%");
    }

    // The single target of a pair search when the position penalties are not Primer3's defaults
    // (compute_position_penalty): [Start, End) 0-based, PRIMER_INSIDE_PENALTY / PRIMER_OUTSIDE_PENALTY.
    private readonly record struct TargetPosition(int Start, int End, double InsidePenalty, double OutsidePenalty);

    // Primer3 secondary-structure values of one primer: ntthal Tm values (thermodynamic screen) or dpal
    // self_any / self_end scores (alignment screen).
    private readonly record struct StructureValues(Primer3OligoStructure? Thermo, double? SelfAny, double? SelfEnd);

    // The per-primer weights with the structure mode of the screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT).
    private static Primer3PenaltyWeights EffectivePenaltyWeights(PrimerParameters param) =>
        (param.PenaltyWeights ?? DefaultPrimer3Weights) with
        {
            ThermodynamicOligoAlignment = param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
            ThermodynamicTemplateAlignment = param.ThermodynamicTemplateAlignment,
        };

    // primer_mispriming_to_template[_thermod]: max_template_mispriming[_th] ≥ 0 and the max score strictly greater.
    private static bool TemplateMisprimingExceeds(TemplateMisprimingScore score, PrimerParameters param)
    {
        double max = param.ActiveMaxTemplateMispriming;
        return max >= 0 && score.Max > max;
    }

    private static string TemplateMisprimingIssue(TemplateMisprimingScore score, PrimerParameters param) =>
        string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"Template mispriming score {score.Max:0.00} exceeds {param.ActiveMaxTemplateMispriming:0.00} (Primer3 PRIMER_MAX_TEMPLATE_MISPRIMING{(param.ThermodynamicTemplateAlignment ? "_TH" : "")})");

    private static bool UsesStructureTerms(Primer3PenaltyWeights w, PrimerStructureScreen screen) => screen switch
    {
        PrimerStructureScreen.Primer3Alignment => w.SelfAny != 0 || w.SelfEnd != 0,
        PrimerStructureScreen.Primer3Thermodynamic => w.SelfAnyTh != 0 || w.SelfEndTh != 0 || w.HairpinTh != 0,
        _ => false,
    };

    private static StructureValues ComputeStructureValues(string seq, PrimerParameters param) => param.StructureScreen switch
    {
        PrimerStructureScreen.Primer3Alignment => new StructureValues(null,
            CalculatePrimerSelfAnyComplementarity(seq), CalculatePrimerSelfEndComplementarity(seq)),
        PrimerStructureScreen.Primer3Thermodynamic => new StructureValues(CalculatePrimer3OligoStructure(seq,
            param.EffectiveMonovalentMillimolar, param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar,
            param.EffectiveDnaConcentrationNanomolar), null, null),
        _ => default,
    };

    // Per-primer secondary-structure screen. Primer3Thermodynamic: Primer3's default
    // (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1) ntthal limits PRIMER_MAX_SELF_ANY_TH / _SELF_END_TH /
    // _HAIRPIN_TH = 47 °C (a non-ACGT primer has no ntthal structure; it is already invalid by Tm).
    // Primer3Alignment: PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0 oligo_compl — dpal self_any > PRIMER_MAX_SELF_ANY,
    // then self_end > PRIMER_MAX_SELF_END (no hairpin value in this mode).
    // Heuristic: the sequence-only library stem-loop screen HasHairpinPotential (default stem 4 = unsourced library
    // threshold, loop 3 = thal.c min_hrpn_loop); not Primer3's hairpin screen.
    private static (bool HasHairpin, StructureValues Structure) AddStructureIssues(
        string seq, PrimerParameters param, List<string> issues)
    {
        if (param.StructureScreen == PrimerStructureScreen.Heuristic)
        {
            bool hp = HasHairpinPotential(seq);
            if (hp)
                issues.Add("Potential hairpin structure detected");
            return (hp, default);
        }

        var values = ComputeStructureValues(seq, param);
        if (param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
        {
            double maxAny = param.EffectiveMaxSelfAny, maxEnd = param.EffectiveMaxSelfEnd;
            if (values.SelfAny > maxAny)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Self-complementarity {values.SelfAny:0.00} exceeds {maxAny:0.00} (Primer3 PRIMER_MAX_SELF_ANY)"));
            if (values.SelfEnd > maxEnd)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"3' self-complementarity {values.SelfEnd:0.00} exceeds {maxEnd:0.00} (Primer3 PRIMER_MAX_SELF_END)"));
            return (false, values);
        }

        if (values.Thermo is not { } v)
            return (false, values);
        double max = param.EffectiveMaxStructureTm;
        bool hasHairpin = v.HairpinTh > max;
        if (hasHairpin)
            issues.Add($"Hairpin melting temperature {v.HairpinTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_HAIRPIN_TH)");
        if (v.SelfAnyTh > max)
            issues.Add($"Self-dimer melting temperature {v.SelfAnyTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_ANY_TH)");
        if (v.SelfEndTh > max)
            issues.Add($"3' self-dimer melting temperature {v.SelfEndTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_END_TH)");
        return (hasHairpin, values);
    }

    /// <summary>
    /// Calculates the "basic" melting temperature for DNA primers (OligoCalc basic Tm,
    /// Kibbe 2007, NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for
    /// &lt; 14 counted bases, and Tm = 64.9 + 41·(G+C − 16.4)/N (customarily attributed to
    /// Marmur &amp; Doty 1962) for ≥ 14 counted bases. Both formulas assume fixed standard conditions
    /// (50 nM primer, 50 mM Na+, pH 7.0); use <see cref="CalculateMeltingTemperatureWithSalt"/>
    /// for another [Na+], or <see cref="CalculateMeltingTemperatureNN"/> for a nearest-neighbor Tm.
    /// A, C, G, T and U are counted (case-insensitive; U read as T, as Biopython <c>MeltingTemp._check</c>
    /// back-transcribes RNA for <c>Tm_Wallace</c>/<c>Tm_GC</c>); all other characters are ignored.
    /// Delegates to the canonical <see cref="ThermoConstants.CalculateBasicTm"/>.
    /// </summary>
    public static double CalculateMeltingTemperature(string primer) => ThermoConstants.CalculateBasicTm(primer);

    /// <summary>
    /// Calculates the OligoCalc "salt adjusted" melting temperature (Kibbe 2007, NAR 35:W43),
    /// the [Na+]-aware counterpart of <see cref="CalculateMeltingTemperature(string)"/>:
    /// <list type="bullet">
    /// <item>&lt; 14 valid bases: Tm = 2(A+T) + 4(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+]) —
    /// the Wallace rule (defined at 50 mM Na+) shifted by the Schildkraut–Lifson relative correction;</item>
    /// <item>≥ 14 valid bases: Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+]).</item>
    /// </list>
    /// [Na+] enters in mol/L (converted from the mM argument). The basic formulas already assume
    /// 50 mM Na+, so the salt term is never added on top of them. Result rounded to one decimal.
    /// A/C/G/T/U are counted as in <see cref="CalculateMeltingTemperature(string)"/> (U read as T);
    /// returns 0 for null/empty input or no counted bases.
    /// </summary>
    /// <param name="primer">Primer sequence.</param>
    /// <param name="naConcentration">Na+ concentration in mM (default: 50). Must be &gt; 0.</param>
    /// <returns>Salt-adjusted melting temperature in °C.</returns>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="naConcentration"/> is not a positive finite number.</exception>
    public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
    {
        if (!(naConcentration > 0) || double.IsInfinity(naConcentration))
            throw new ArgumentOutOfRangeException(nameof(naConcentration), naConcentration,
                "Na+ concentration must be a positive finite value in mM.");

        if (string.IsNullOrEmpty(primer))
            return 0;

        var (at, gc) = ThermoConstants.CountBasicTmBases(primer);
        if (at + gc == 0)
            return 0;

        return Math.Round(
            ThermoConstants.CalculateOligoCalcSaltAdjustedTm(at, gc, naConcentration / 1000.0), 1);
    }

    // ---- Primer3 default primer Tm (PRIMER-DESIGN-001) -------------------------------------
    // Reproduces Primer3's seqtm()/oligotm() with PRIMER_TM_FORMULA = 1 (SantaLucia 1998) and
    // PRIMER_SALT_CORRECTIONS = 1 (SantaLucia 1998), the Primer3 ≥ 2.0 defaults, including the
    // divalent→monovalent equivalence of von Ahsen et al. (2001) and the long_seq_tm fallback
    // for oligos longer than MAX_NN_TM_LENGTH (36). Source: primer3 src/oligotm.c (oligotm,
    // seqtm, long_seq_tm, divalent_to_monovalent, symmetry; SantaLucia_1998_dH/dS tables);
    // SantaLucia (1998) PNAS 95:1460 Table 2 (NN + terminal-initiation terms); von Ahsen et al.
    // (2001) Clin Chem 47:1956. Cross-checked against primer3-py 2.3.1 primer3.calc_tm.

    /// <summary>Primer3 MAX_NN_TM_LENGTH: longer oligos use the long_seq_tm GC/length formula.</summary>
    public const int Primer3MaxNnTmLength = 36;

    // SantaLucia (1998) Table 2 unified NN parameters exactly as tabulated in Primer3 oligotm.c
    // (SantaLucia_1998_dH in units of −100 cal/mol, SantaLucia_1998_dS in units of −0.1 cal/(K·mol);
    // index order A, C, G, T). Primer3 sums these as integers, so primers with the same
    // nearest-neighbour multiset get bit-identical Tm values (this matters for tie-breaking in
    // DesignPrimers). Differs from the SantaLucia & Hicks (2004) set (ThermoConstants NnParameterSet.SantaLuciaHicks2004) in
    // AA/TT (−7.9/−22.2 vs −7.6/−21.3).
    private static readonly int[,] Primer3SantaLucia1998Dh =
    {
        { 79, 84, 78, 72 },
        { 85, 80, 106, 78 },
        { 82, 98, 80, 84 },
        { 72, 82, 85, 79 }
    };
    private static readonly int[,] Primer3SantaLucia1998Ds =
    {
        { 222, 224, 210, 204 },
        { 227, 199, 272, 210 },
        { 222, 244, 199, 224 },
        { 213, 222, 227, 222 }
    };

    // SantaLucia (1998) Table 2 initiation per duplex end, in the same integer units
    // (oligotm.c "Terminal penalty"): terminal A·T ΔH +2.3 / ΔS +4.1 → (−23, −41);
    // terminal G·C ΔH +0.1 / ΔS −2.8 → (−1, +28); symmetry ΔS −1.4 → +14.
    private const int Primer3TerminalAtDh = -23, Primer3TerminalAtDs = -41;
    private const int Primer3TerminalGcDh = -1, Primer3TerminalGcDs = 28;
    private const int Primer3SymmetryDs = 14;

    // oligotm.c constants: R = 1.987 cal/(K·mol); divalent_to_monovalent factor 120 (mM units).
    private const double Primer3GasConstant = 1.987;
    private const double Primer3DivalentFactor = 120.0;

    /// <summary>
    /// Primer3's default primer melting temperature (°C) — Primer3 <c>seqtm</c> with
    /// PRIMER_TM_FORMULA = SantaLucia 1998 and PRIMER_SALT_CORRECTIONS = SantaLucia 1998:
    /// <para>
    /// [Mon]_eq = [Mon] + 120·√([Mg²⁺] − [dNTP]) (mM; von Ahsen 2001, 0 when [Mg²⁺] ≤ [dNTP]);
    /// ΔS = ΔS°(1 M) + 0.368·(N − 1)·ln([Mon]_eq / 1000);
    /// Tm = ΔH / (ΔS + R·ln(C/4)) − 273.15, with C/1 for a self-complementary oligo,
    /// R = 1.987, C = DNA concentration in mol/L;
    /// for N &gt; <see cref="Primer3MaxNnTmLength"/>: Tm = 81.5 + 16.6·log10([Mon]_eq/1000) + 41·GC/N − 600/N.
    /// </para>
    /// Defaults are Primer3's PRIMER_DNA_CONC = 50 nM, PRIMER_SALT_MONOVALENT = 50 mM,
    /// PRIMER_SALT_DIVALENT = 1.5 mM, PRIMER_DNTP_CONC = 0.6 mM. Matches primer3-py
    /// <c>calc_tm</c> (e.g. AGCTAGCTAGCTAGCTAGCT → 58.101 °C).
    /// </summary>
    /// <param name="primer">Primer sequence (case-insensitive).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when the sequence is null/shorter than 2 bases or contains a
    /// non-ACGT character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">For a non-positive DNA concentration, negative ion
    /// concentrations, or a zero total monovalent-equivalent concentration.</exception>
    public static double CalculateMeltingTemperaturePrimer3(
        string primer,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar) =>
        CalculateMeltingTemperaturePrimer3(primer, dnaConcentrationNanomolar, monovalentMillimolar,
            divalentMillimolar, dntpMillimolar, Primer3MaxNnTmLength);

    /// <summary>
    /// Primer3 <c>seqtm</c> (<see cref="CalculateMeltingTemperaturePrimer3(string, double, double, double, double)"/>)
    /// with the nearest-neighbour length limit as an argument: oligos of at most
    /// <paramref name="maxNearestNeighborLength"/> nt get the SantaLucia 1998 nearest-neighbour Tm, longer ones
    /// Primer3's <c>long_seq_tm</c> (<c>oligotm.c</c> <c>seqtm</c>: <c>if (len &gt; nn_max_len) long_seq_tm … else
    /// oligotm</c>; primer3-py <c>calc_tm(…, max_nn_length)</c>). Primer3's own primer/probe picker uses
    /// <see cref="Primer3MaxNnTmLength"/> = 36; e.g. OligoArray-style microarray probes of 50–60 nt need the
    /// nearest-neighbour Tm over the whole oligo (limit ≥ 60).
    /// </summary>
    /// <param name="primer">Oligo sequence (case-insensitive).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <param name="maxNearestNeighborLength">Longest oligo (nt) that gets the nearest-neighbour Tm (Primer3
    /// <c>nn_max_len</c>, ≥ 0).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when the sequence is null/shorter than 2 bases or contains a
    /// non-ACGT character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">For a non-positive DNA concentration, negative ion
    /// concentrations, a zero total monovalent-equivalent concentration or a negative
    /// <paramref name="maxNearestNeighborLength"/>.</exception>
    public static double CalculateMeltingTemperaturePrimer3(
        string primer,
        double dnaConcentrationNanomolar,
        double monovalentMillimolar,
        double divalentMillimolar,
        double dntpMillimolar,
        int maxNearestNeighborLength)
    {
        if (maxNearestNeighborLength < 0)
            throw new ArgumentOutOfRangeException(nameof(maxNearestNeighborLength), maxNearestNeighborLength,
                "The nearest-neighbour length limit must be ≥ 0 nt.");
        if (!(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(nameof(dnaConcentrationNanomolar), dnaConcentrationNanomolar,
                "DNA concentration must be a positive finite value in nM.");
        if (!(monovalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Monovalent cation concentration must be ≥ 0 mM.");
        if (!(divalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(divalentMillimolar), divalentMillimolar,
                "Divalent cation concentration must be ≥ 0 mM.");
        if (!(dntpMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMillimolar), dntpMillimolar,
                "dNTP concentration must be ≥ 0 mM.");

        double monovalentEq = Primer3MonovalentEquivalent(monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        if (!(monovalentEq > 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Total monovalent-equivalent cation concentration must be > 0 mM.");

        return Primer3SeqTm(primer, dnaConcentrationNanomolar, monovalentEq, Primer3DefaultAnnealingTemperature,
            maxNearestNeighborLength).Tm;
    }

    // Primer3 seqtm() / oligotm() (SantaLucia 1998 Tm, SantaLucia salt correction, nn_max_len = maxNnLength, Primer3's
    // MAX_NN_TM_LENGTH = 36 by default) with the
    // fraction bound at the annealing temperature (oligotm.c, PRIMER_ANNEALING_TEMP): Tm is NaN when Primer3 reports
    // OLIGOTM_ERROR (fewer than 2 bases, a non-ACGT base); Bound is Primer3OligoTmError unless annealingTemperature > 0
    // and the nearest-neighbour branch applies (long_seq_tm leaves bound = OLIGOTM_ERROR).
    // Input: validated conditions, monovalentEq = [Mon] + 120·√([Mg²⁺] − [dNTP]) in mM.
    internal static (double Tm, double Bound) Primer3SeqTm(
        string? primer, double dnaConcentrationNanomolar, double monovalentEq, double annealingTemperature,
        int maxNnLength = Primer3MaxNnTmLength)
    {
        if (string.IsNullOrEmpty(primer) || primer.Length < 2)
            return (double.NaN, Primer3OligoTmError);

        string seq = primer.ToUpperInvariant();
        int n = seq.Length;
        int gc = 0;
        foreach (char c in seq)
        {
            if (c is 'G' or 'C') gc++;
            else if (c is not ('A' or 'T')) return (double.NaN, Primer3OligoTmError);
        }

        if (n > maxNnLength)
            return (ThermoConstants.CalculateSaltAdjustedTm((double)gc / n, n, monovalentEq / 1000.0), Primer3OligoTmError);

        // oligotm(): integer accumulation, then ΔH = dh·(−100) cal/mol, ΔS = ds·(−0.1) cal/(K·mol).
        bool symmetric = IsSelfComplementary(seq);
        int dh = 0, ds = symmetric ? Primer3SymmetryDs : 0;
        foreach (char end in new[] { seq[0], seq[^1] })
        {
            if (end is 'A' or 'T') { dh += Primer3TerminalAtDh; ds += Primer3TerminalAtDs; }
            else { dh += Primer3TerminalGcDh; ds += Primer3TerminalGcDs; }
        }
        for (int i = 0; i < n - 1; i++)
        {
            int x = BaseIndex(seq[i]), y = BaseIndex(seq[i + 1]);
            dh += Primer3SantaLucia1998Dh[x, y];
            ds += Primer3SantaLucia1998Ds[x, y];
        }

        double deltaH = dh * -100.0;
        double deltaS = ds * -0.1;
        deltaS += SantaLuciaEntropySaltCoefficient * (n - 1) * Math.Log(monovalentEq / 1000.0);
        // Equation A (self-complementary, C_T/1) or Equation B (C_T/4) of oligotm.c.
        double strandDivisor = symmetric ? 1000000000.0 : 4000000000.0;
        double tm = deltaH / (deltaS + Primer3GasConstant * Math.Log(dnaConcentrationNanomolar / strandDivisor)) - KelvinOffset;
        double bound = Primer3OligoTmError;
        if (annealingTemperature > 0.0)
        {
            // oligotm.c (santalucia salt correction): ddG = ΔH − (Ta + 273.15)·ΔS; Ka = exp(−ddG / (R·(Ta + 273.15)));
            // bound = 100 / (1 + √(1 / ((C/x)·Ka))), x = 1e9 (Equation A) or 4e9 (Equation B).
            double ddG = deltaH - (annealingTemperature + KelvinOffset) * deltaS;
            double ka = Math.Exp(-ddG / (Primer3GasConstant * (annealingTemperature + KelvinOffset)));
            bound = (1 / (1 + Math.Sqrt(1 / ((dnaConcentrationNanomolar / strandDivisor) * ka)))) * 100;
        }
        return (tm, bound);

        static int BaseIndex(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };
    }

    /// <summary>
    /// Calculates GC content as a percentage.
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcContentFast();

    /// <summary>
    /// Finds the longest homopolymer run (consecutive identical nucleotides, case-insensitive) — the
    /// quantity Primer3 limits with <c>PRIMER_MAX_POLY_X</c> (a primer fails when the run exceeds
    /// it; <c>libprimer3.cc</c> <c>_pr_violates_poly_x</c>). As in Primer3, N is a wildcard that takes
    /// the worst case: a forward scan assigns each N to the preceding non-N base, a reverse scan to
    /// the following one, and the longer run is reported (e.g. ANA → 3, GNGNG → 5, ANGNG → 4).
    /// </summary>
    /// <returns>0 for null/empty input, otherwise ≥ 1.</returns>
    public static int FindLongestHomopolymer(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string seq = sequence.ToUpperInvariant();
        int len = seq.Length;

        // Forward scan (N counts as the last seen non-N base; leading Ns as the first non-N base).
        char lastNonN = seq[0];
        if (lastNonN == 'N')
        {
            for (int i = 1; i < len; i++)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        int run = 1, maxRun = 1;
        bool hasN = false;
        for (int i = 1; i < len; i++)
        {
            if (seq[i] == 'N') { hasN = true; run++; }
            else if (seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        if (!hasN)
            return maxRun;

        // Reverse scan (N counts as the next non-N base).
        lastNonN = seq[len - 1];
        if (lastNonN == 'N')
        {
            for (int i = len - 2; i >= 0; i--)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        run = 1;
        for (int i = len - 2; i >= 0; i--)
        {
            if (seq[i] == 'N' || seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        return maxRun;
    }

    /// <summary>
    /// Library screen (unsourced heuristic, no Primer3 counterpart): the number of consecutive copies of the
    /// most-repeated 2-mer (e.g. ATATAT = 3 units of AT), case-insensitive; inputs shorter than 4 nt return 0.
    /// Any 2-mer counts, including homo-dinucleotides (AAAA = 2 units of AA), and no base is treated as a
    /// wildcard. <see cref="PrimerParameters.MaxDinucleotideRepeats"/> limits this count in
    /// <see cref="EvaluatePrimer"/> (library default 4 in <see cref="DefaultParameters"/>). No authoritative
    /// published definition of this screen or of the limit 4 was found (audit round 3, A3-8): Primer3 has no
    /// dinucleotide-repeat setting (<c>primer3_manual.htm</c> lists only PRIMER_MAX_POLY_X, the mononucleotide run =
    /// <see cref="FindLongestHomopolymer"/>, and repeat-library mispriming PRIMER_MISPRIMING_LIBRARY), so
    /// <see cref="Primer3DefaultParameters"/> disables it (<see cref="int.MaxValue"/>). For Primer3's sourced
    /// alternatives use PRIMER_MAX_POLY_X (<see cref="PrimerParameters.MaxHomopolymer"/>) and the mispriming
    /// library PRIMER_MISPRIMING_LIBRARY (<see cref="PrimerParameters.MisprimingLibrary"/>).
    /// </summary>
    /// <param name="sequence">Nucleotide sequence.</param>
    /// <returns>Unit count of the longest dinucleotide tandem repeat (≥ 1 for inputs of ≥ 4 nt), 0 below 4 nt.</returns>
    public static int FindLongestDinucleotideRepeat(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 4)
            return 0;

        var seq = sequence.ToUpperInvariant();
        int maxRepeats = 0;

        for (int i = 0; i < seq.Length - 3; i++)
        {
            string dinuc = seq.Substring(i, 2);
            int repeats = 1;
            int j = i + 2;

            while (j + 1 < seq.Length && seq.Substring(j, 2) == dinuc)
            {
                repeats++;
                j += 2;
            }

            maxRepeats = Math.Max(maxRepeats, repeats);
        }

        return maxRepeats;
    }

    /// <summary>
    /// Library screen (sequence-only heuristic, no Primer3 counterpart): <c>true</c> when the sequence contains
    /// two non-overlapping segments of <paramref name="minStemLength"/> bases that are exact Watson–Crick reverse
    /// complements of each other (an antiparallel stem) separated by at least
    /// <paramref name="minLoopLength"/> unpaired bases. Case-insensitive; no G·T wobble, mismatches or energies.
    /// The default minimum loop of 3 nt is Primer3's hairpin minimum (<c>thal.c</c> <c>min_hrpn_loop = 3</c>); the
    /// default minimum stem of 4 bp is an unsourced library threshold (no published definition of a
    /// "≥ 4-bp stem + ≥ 3-nt loop" rule was found, audit round 3, A3-8). It is not equivalent to Primer3's
    /// hairpin screen: e.g. <c>AAAACCCTTTT</c> is flagged although Primer3 <c>calc_hairpin</c> finds no structure,
    /// and <c>CAGTAAAACCCTTTTGCAGC</c> is flagged although its hairpin Tm is 37.65 °C (primer3-py 2.3.1,
    /// 50 mM / 1.5 mM / 0.6 mM / 50 nM), below PRIMER_MAX_HAIRPIN_TH = 47 °C. The sourced Primer3 hairpin screen
    /// (PRIMER_HAIRPIN_TH) is <see cref="CalculatePrimer3OligoStructure"/> /
    /// <see cref="CalculateHairpinThermodynamicsNtthal(string, double)"/>, which
    /// <see cref="EvaluatePrimer"/> uses by default; this screen is used only with
    /// <see cref="PrimerStructureScreen.Heuristic"/>. Uses an O(n²) scan below 100 nt and a suffix tree
    /// at ≥ 100 nt (identical results).
    /// </summary>
    /// <param name="sequence">DNA sequence to check.</param>
    /// <param name="minStemLength">Minimum stem length (default 4).</param>
    /// <param name="minLoopLength">Minimum loop length (default 3).</param>
    /// <returns>True if hairpin potential detected.</returns>
    public static bool HasHairpinPotential(string sequence, int minStemLength = 4, int minLoopLength = 3)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minStemLength * 2 + minLoopLength)
            return false;

        var seq = sequence.ToUpperInvariant();

        // For short sequences (typical primers), use simple O(n²) approach
        // Break-even point is ~100bp based on suffix tree construction overhead
        if (seq.Length < 100)
        {
            return HasHairpinPotentialSimple(seq, minStemLength, minLoopLength);
        }

        // For longer sequences, use suffix tree for O(n) lookup
        return HasHairpinPotentialWithSuffixTree(seq, minStemLength, minLoopLength);
    }

    /// <summary>
    /// Simple O(n²) hairpin detection for short sequences.
    /// </summary>
    private static bool HasHairpinPotentialSimple(string seq, int minStemLength, int minLoopLength)
    {
        // Check for self-complementary regions
        for (int i = 0; i <= seq.Length - minStemLength; i++)
        {
            string fragment = seq.Substring(i, minStemLength);
            // Look for complementary sequence at least minLoopLength positions away
            for (int j = i + minStemLength + minLoopLength; j <= seq.Length - minStemLength; j++)
            {
                string target = seq.Substring(j, minStemLength);
                if (AreComplementary(fragment, Reverse(target)))
                    return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Suffix tree-based O(n) hairpin detection for long sequences.
    /// 
    /// Algorithm: A hairpin forms when a substring S at position i is complementary
    /// to a substring at position j (in reverse). This is equivalent to:
    /// - seq[i..i+k] being present in revComp at some position p
    /// - The positions must satisfy: j = n - p - k, and j >= i + k + minLoopLength
    /// 
    /// We build a suffix tree on seq and search for all substrings of revComp,
    /// checking if any match satisfies the loop constraint.
    /// </summary>
    private static bool HasHairpinPotentialWithSuffixTree(string seq, int minStemLength, int minLoopLength)
    {
        var revComp = DnaSequence.GetReverseComplementString(seq);
        var tree = global::SuffixTree.SuffixTree.Build(seq);

        // For each position in revComp, find matches in seq via suffix tree
        // and check if they form valid hairpin (sufficient loop distance)
        int n = seq.Length;

        // Slide through revComp looking for stems
        for (int p = 0; p <= n - minStemLength; p++)
        {
            var pattern = revComp.AsSpan(p, minStemLength);
            var matches = tree.FindAllOccurrences(pattern);

            foreach (int i in matches)
            {
                // Position in revComp p corresponds to position (n - p - minStemLength) in seq
                // when we reverse complement back
                int j = n - p - minStemLength;

                // Check if positions form valid hairpin: j >= i + minStemLength + minLoopLength
                // Also check i and j don't overlap with the stem itself
                if (j >= i + minStemLength + minLoopLength && j + minStemLength <= n)
                {
                    return true;
                }
                // Also check the reverse case: i is the 3' stem, j is the 5' stem
                if (i >= j + minStemLength + minLoopLength && i + minStemLength <= n)
                {
                    return true;
                }
            }
        }

        return false;
    }

    /// <summary>
    /// Checks whether two primers can form a 3′-end primer-dimer using Primer3's alignment-based
    /// (non-thermodynamic, <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>) pair 3′-complementarity:
    /// returns <c>true</c> when <see cref="CalculatePrimerDimerEndComplementarity"/> ≥
    /// <paramref name="minComplementarity"/>. For ACGT primers the score is an integer, so the
    /// default threshold 4 is exactly Primer3's default <c>PRIMER_PAIR_MAX_COMPL_END = 3.00</c>
    /// (a pair fails when <c>compl_end &gt; 3</c>). The thermodynamic counterpart (Primer3's default
    /// mode) is <see cref="CalculatePrimer3PairComplementarity"/>.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′).</param>
    /// <param name="primer2">Second primer (5′→3′).</param>
    /// <param name="minComplementarity">Minimum 3′-anchored complementarity score that flags a dimer (default 4).</param>
    /// <returns><c>false</c> for null/empty primers.</returns>
    public static bool HasPrimerDimer(string primer1, string primer2, int minComplementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return false;
        return CalculatePrimerDimerEndComplementarity(primer1, primer2) >= minComplementarity;
    }

    /// <summary>
    /// Primer3 alignment-mode pair 3′-complementarity <c>compl_end</c> (Primer3 <c>libprimer3.cc</c>
    /// <c>characterize_pair</c>, <c>PRIMER_PAIR_COMPL_END</c> with
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; Rozen &amp; Skaletsky 2000; Untergasser et al. 2012):
    /// the maximum of <c>align(p1, revcomp(p2))</c> and <c>align(p2, revcomp(p1))</c>, where
    /// <c>align</c> is the <c>dpal</c> end-anchored (<c>DPAL_GLOBAL_END</c>) alignment — it must end
    /// at the 3′-terminal base of the first sequence — scored +1 per complementary base pair, −1 per
    /// mismatch, −0.25 against N (any non-ACGT character is scored as N), −2 per single-base gap
    /// (max gap 1), and floored at 0. Two primers whose 3′-terminal k bases are reverse complements
    /// score k; a homopolymer against itself (e.g. A₈/A₈) scores 0 because it cannot pair.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′), case-insensitive.</param>
    /// <param name="primer2">Second primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>compl_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerDimerEndComplementarity(string primer1, string primer2)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return 0;
        string p1 = primer1.ToUpperInvariant();
        string p2 = primer2.ToUpperInvariant();
        // Primer3 compares s1 with s2 taken from the same (top) strand, i.e. the right primer's
        // template-strand copy = revcomp(right primer); then also s2_rev against s1_rev.
        double a = DpalGlobalEndScore(p1, ReverseComplementPrimer3(p2));
        double b = DpalGlobalEndScore(p2, ReverseComplementPrimer3(p1));
        return Math.Max(a, b);
    }

    /// <summary>
    /// Primer3 alignment-mode self 3′-complementarity <c>self_end</c>
    /// (<c>oligo_compl</c>: <c>align(oligo, revcomp(oligo), DPAL_GLOBAL_END)</c>, Primer3
    /// <c>PRIMER_LEFT/RIGHT_SELF_END</c>; default limit <c>PRIMER_MAX_SELF_END = 3.00</c>).
    /// </summary>
    /// <param name="primer">Primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>self_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerSelfEndComplementarity(string primer)
    {
        if (string.IsNullOrEmpty(primer))
            return 0;
        string p = primer.ToUpperInvariant();
        return DpalGlobalEndScore(p, ReverseComplementPrimer3(p));
    }

    /// <summary>
    /// Primer3 alignment-mode self-complementarity <c>self_any</c> (<c>libprimer3.cc</c> <c>oligo_compl</c>:
    /// <c>align(oligo, revcomp(oligo), DPAL_LOCAL)</c>; Primer3 <c>PRIMER_LEFT/RIGHT/INTERNAL_n_SELF_ANY</c> with
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; default limits <c>PRIMER_MAX_SELF_ANY = 8.00</c>,
    /// <c>PRIMER_INTERNAL_MAX_SELF_ANY = 12.00</c>): the best <c>dpal</c> local alignment of the oligo with its
    /// reverse complement, scored +1 per complementary base pair, −1 per mismatch, −0.25 against N (any non-ACGT
    /// character is scored as N), −2 per single-base gap (max gap 1), never below 0 — a line-by-line port of
    /// dpal.c <c>_dpal_long_nopath_maxgap1_local</c> (the routine Primer3 runs in DPM_FAST mode) with Primer3's
    /// <c>align</c> rule that a second sequence shorter than 3 scores its length. No length limit (dpal's
    /// score-only routine is linear in memory). Bit-exact to dpal.c / primer3-py 2.3.1.
    /// </summary>
    /// <param name="oligo">Oligo (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>self_any</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerSelfAnyComplementarity(string oligo)
    {
        if (string.IsNullOrEmpty(oligo))
            return 0;
        string p = oligo.ToUpperInvariant();
        return DpalLocalScore(p, ReverseComplementPrimer3(p));
    }

    /// <summary>
    /// Primer3 alignment-mode pair complementarity <c>compl_any</c> (<c>libprimer3.cc</c> <c>characterize_pair</c>,
    /// <c>PRIMER_PAIR_COMPL_ANY</c> with <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; default limit
    /// <c>PRIMER_PAIR_MAX_COMPL_ANY = 8.00</c>): <c>align(s1, s2, DPAL_LOCAL)</c> where s1 is the left primer and s2 the
    /// right primer's top-strand copy, i.e. the dpal local alignment of <paramref name="leftPrimer"/> with
    /// revcomp(<paramref name="rightPrimer"/>) (scoring as <see cref="CalculatePrimerSelfAnyComplementarity"/>).
    /// Primer3 evaluates this one orientation only.
    /// </summary>
    /// <param name="leftPrimer">Left (forward) primer, 5′→3′, case-insensitive.</param>
    /// <param name="rightPrimer">Right (reverse) primer, 5′→3′, case-insensitive.</param>
    /// <returns>The Primer3 <c>compl_any</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerDimerAnyComplementarity(string leftPrimer, string rightPrimer)
    {
        if (string.IsNullOrEmpty(leftPrimer) || string.IsNullOrEmpty(rightPrimer))
            return 0;
        return DpalLocalScore(leftPrimer.ToUpperInvariant(), ReverseComplementPrimer3(rightPrimer.ToUpperInvariant()));
    }

    /// <summary>Primer3 default PRIMER_MAX_SELF_ANY (alignment mode, <c>pr_set_default_global_args_1</c>).</summary>
    public const double Primer3MaxSelfAny = 8.0;

    /// <summary>Primer3 default PRIMER_MAX_SELF_END (alignment mode).</summary>
    public const double Primer3MaxSelfEnd = 3.0;

    /// <summary>Primer3 default PRIMER_PAIR_MAX_COMPL_ANY (alignment mode).</summary>
    public const double Primer3MaxPairComplAny = 8.0;

    /// <summary>Primer3 default PRIMER_PAIR_MAX_COMPL_END (alignment mode).</summary>
    public const double Primer3MaxPairComplEnd = 3.0;

    /// <summary>Primer3 default PRIMER_INTERNAL_MAX_SELF_ANY and PRIMER_INTERNAL_MAX_SELF_END (alignment mode, 12.00).</summary>
    public const double Primer3InternalMaxSelfComplementarity = 12.0;

    // Primer3 oligo_compl decision (self_any > maxSelfAny, else self_end > maxSelfEnd) without computing the full
    // scores: the dpal recurrences stop as soon as the running optimum exceeds the limit (the decision is identical
    // to comparing the full scores, since the optimum never decreases). Used by the probe fallback screen, where
    // probes of several hundred nt make the full O(n²) scores the dominant cost.
    internal static bool ExceedsPrimer3SelfComplementarity(string oligo, double maxSelfAny, double maxSelfEnd,
        out bool selfAnyExceeded)
    {
        selfAnyExceeded = false;
        if (string.IsNullOrEmpty(oligo))
            return false;
        string x = oligo.ToUpperInvariant();
        string y = ReverseComplementPrimer3(x);
        if (y.Length < 3 ? y.Length > maxSelfAny : DpalLocalFast(x, y, maxSelfAny * 100.0) / 100.0 > maxSelfAny)
            return selfAnyExceeded = true;
        return DpalGlobalEndScore(x, y) > maxSelfEnd;
    }

    // Primer3 align() with DPAL_LOCAL: a second sequence shorter than 3 "maxes out" the score to its length
    // (not divided by 100); otherwise dpal score / 100, floored at 0.
    private static double DpalLocalScore(string x, string y)
    {
        if (y.Length < 3)
            return y.Length;
        int smax = DpalLocalFast(x, y);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    // dpal.c _dpal_long_nopath_maxgap1_local (score only, max gap 1, cells floored at 0). With a finite
    // stopAbove the scan returns as soon as the running optimum exceeds it (the value is then a lower bound).
    // matrix: a dpal substitution matrix (MatrixSsm); null = DpalSsm (non-ACGT scored as N).
    private static int DpalLocalFast(string x, string y, double stopAbove = double.PositiveInfinity, int[]? matrix = null)
    {
        int Ssm(char a, char b) => matrix is null ? DpalSsm(a, b) : MatrixSsm(matrix, a, b);
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[ylen];
        var s1 = new int[ylen];
        var s2 = new int[ylen];
        int smax = 0, score, a;

        // Row 0.
        for (int j = 0; j < ylen; j++)
        {
            score = Ssm(x[0], y[j]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s0[j] = score;
        }
        // Row 1 (for |X| = 1 the C code reads X[1] = NUL, whose ssm row is INT_MIN: every cell is 0).
        if (xlen == 1)
            return smax;
        score = Ssm(x[1], y[0]);
        if (score < 0) score = 0;
        else if (score > smax) smax = score;
        s1[0] = score;
        for (int j = 1; j < ylen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += Ssm(x[1], y[j]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s1[j] = score;
        }

        for (int i = 2; i < xlen; i++)
        {
            score = Ssm(x[i], y[0]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s2[0] = score;
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += Ssm(x[i], y[1]);
            if (score < 0) score = 0;
            else if (score > smax) smax = score;
            s2[1] = score;
            for (int j = 2; j < ylen; j++)
            {
                score = s0[j - 1];
                if ((a = s1[j - 2]) > score) score = a;
                score += gap;
                if ((a = s1[j - 1]) > score) score = a;
                score += Ssm(x[i], y[j]);
                if (score < 0) score = 0;
                else if (score > smax) smax = score;
                s2[j] = score;
            }
            if (smax > stopAbove)
                return smax;
            (s0, s1, s2) = (s1, s2, s0);
        }
        return smax;
    }

    // Primer3 p3_reverse_complement: ACGT complemented, every other character becomes N.
    private static string ReverseComplementPrimer3(string seq)
    {
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            chars[seq.Length - 1 - i] = seq[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => 'N',
            };
        }
        return new string(chars);
    }

    // dpal default primer-picking scoring (dpal.c set_dpal_args): identity matrix ×100.
    private const int DpalMatch = 100, DpalMismatch = -100, DpalN = -25, DpalGap = -200;

    private static int DpalSsm(char x, char y)
    {
        bool xn = x is not ('A' or 'C' or 'G' or 'T');
        bool yn = y is not ('A' or 'C' or 'G' or 'T');
        if (xn || yn) return DpalN;
        return x == y ? DpalMatch : DpalMismatch;
    }

    /// <summary>
    /// Primer3 <c>align(X, Y, DPAL_GLOBAL_END)</c> score / 100 floored at 0 (libprimer3.cc
    /// <c>align</c>), computed by a line-by-line port of dpal.c
    /// <c>_dpal_long_nopath_maxgap1_global_end</c> (the routine Primer3 runs in DPM_FAST mode with
    /// max_gap = 1); inputs too short for that routine (|X| ≤ 3 or |Y| = 1, where the C code reads
    /// past the sequence end) use dpal.c's <c>_dpal_generic</c> recurrence for GLOBAL_END.
    /// </summary>
    private static double DpalGlobalEndScore(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        // The fast routine reads Y[ylen] (C's NUL terminator) when |X| ≤ 3 and Y[1] when |Y| = 1.
        int smax = xlen < 4 || ylen < 2 ? DpalGlobalEndGeneric(x, y) : DpalGlobalEndFast(x, y);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    private static int DpalGlobalEndFast(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[xlen];
        var s1 = new int[xlen];
        var s2 = new int[xlen];
        int score, a;

        int smax = DpalSsm(x[xlen - 1], y[0]);
        for (int j = 0; j < xlen; j++) s0[j] = DpalSsm(x[j], y[0]);

        s1[0] = DpalSsm(x[0], y[1]);
        for (int j = 1; j < xlen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += DpalSsm(x[j], y[1]);
            if (score > smax && j == xlen - 1) smax = score;
            s1[j] = score;
        }

        int k = ylen - xlen / 2 + 1;
        if (k < 1) k = 1;

        // Rectangular part.
        for (int j = 2; j < k + 1; j++)
        {
            s2[0] = DpalSsm(x[0], y[j]);
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += DpalSsm(x[1], y[j]);
            s2[1] = score;
            for (int i = 2; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }

        // Triangular part (cells left of the band are not recomputed, exactly as dpal.c).
        int t = 2;
        for (int j = k + 1; j < ylen; j++)
        {
            for (int i = t; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            t += 2;
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }
        return smax;
    }

    // dpal.c _dpal_generic with flag DPAL_GLOBAL_END and max_gap = 1 (score only).
    private static int DpalGlobalEndGeneric(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        var sm = new int[xlen, ylen];
        for (int i = 0; i < xlen; i++) sm[i, 0] = DpalSsm(x[i], y[0]);
        int smax = sm[xlen - 1, 0];
        for (int j = 0; j < ylen; j++) sm[0, j] = DpalSsm(x[0], y[j]);
        for (int i = 1; i < xlen; i++)
        {
            for (int j = 1; j < ylen; j++)
            {
                long a = sm[i - 1, j - 1];
                long b = i > 1 ? (long)sm[i - 2, j - 1] + DpalGap : long.MinValue;
                long c = j > 1 ? (long)sm[i - 1, j - 2] + DpalGap : long.MinValue;
                long best;
                if (a >= b && a >= c) best = a;
                else if (b > a && b >= c) best = b;
                else best = c;
                int score = (int)(best + DpalSsm(x[i], y[j]));
                if (score >= smax && i == xlen - 1) smax = score;
                sm[i, j] = score;
            }
        }
        return smax;
    }

    /// <summary>
    /// 3′-end stability: the duplex ΔG°37 (kcal/mol, 1 M NaCl) of the last five bases (the whole
    /// primer when it is shorter than five), computed exactly as Primer3's
    /// <c>end_oligodg(seq, 5, santalucia)</c> (<c>oligotm.c</c> <c>oligodg</c>; the value Primer3 reports
    /// as <c>PRIMER_{LEFT,RIGHT}_n_END_STABILITY</c> and limits with <c>PRIMER_MAX_END_STABILITY</c>):
    /// SantaLucia (1998) Table 1 unified NN ΔG°37 values summed over the steps, plus initiation
    /// +1.96, +0.05 per terminal A·T base pair and +0.43 for a self-complementary (even-length)
    /// sequence. For a 5-mer this is identical to SantaLucia (1998)'s "initiation with terminal G·C
    /// +0.98 / terminal A·T +1.03" form (0.98+0.98 = 1.96, 1.03−0.98 = 0.05). N is accepted with
    /// Primer3's N-row values (NA 0.58, NC 1.30, NG 1.28, NT 0.88, NN 0.58 negated, etc.) and no A·T
    /// penalty.
    /// <para><b>Sign convention:</b> this method returns the physical ΔG (negative = stable); Primer3's
    /// END_STABILITY is the same quantity with the opposite sign (GCGCG → −6.86 here, 6.86 in
    /// Primer3; TATAT → −0.86 / 0.86).</para>
    /// </summary>
    /// <param name="sequence">Primer sequence (5′→3′), case-insensitive.</param>
    /// <returns>ΔG°37 in kcal/mol; 0 for null/empty input; <c>double.NaN</c> when the 3′ window
    /// contains a character other than A, C, G, T, N (Primer3 <c>OLIGOTM_ERROR</c>).</returns>
    public static double Calculate3PrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string window = (sequence.Length > 5 ? sequence[^5..] : sequence).ToUpperInvariant();
        var idx = new int[window.Length];
        for (int i = 0; i < window.Length; i++)
        {
            int b = window[i] switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4, _ => -1 };
            if (b < 0) return double.NaN;
            idx[i] = b;
        }

        // oligodg (santalucia): dg = −1960 [−430 if symmetric] [−50 per terminal A/T] + Σ table, in cal/mol
        // of −ΔG; the method returns ΔG = −dg/1000.
        int dg = -1960;
        if (IsPrimer3Symmetric(window)) dg += -430;
        if (window[0] is 'A' or 'T') dg += -50;
        for (int i = 0; i + 1 < idx.Length; i++)
            dg += Primer3SantaLucia1998Dg[idx[i], idx[i + 1]];
        if (window[^1] is 'A' or 'T') dg += -50;

        return -dg / 1000.0;
    }

    // oligotm.c SantaLucia_1998_dG (−ΔG°37, cal/mol), rows/cols A, C, G, T, N.
    private static readonly int[,] Primer3SantaLucia1998Dg =
    {
        { 1000, 1440, 1280,  880,  880 },
        { 1450, 1840, 2170, 1280, 1450 },
        { 1300, 2240, 1840, 1440, 1300 },
        {  580, 1300, 1450, 1000,  580 },
        {  580, 1300, 1280,  880,  580 },
    };

    // oligotm.c symmetry(): even length and every A/T and C/G position Watson-Crick paired with its mirror.
    private static bool IsPrimer3Symmetric(string seq)
    {
        int n = seq.Length;
        if (n % 2 == 1) return false;
        for (int i = 0; i < n / 2; i++)
        {
            char s = seq[i], e = seq[n - 1 - i];
            if ((s == 'A' && e != 'T') || (s == 'T' && e != 'A') || (e == 'A' && s != 'T') || (e == 'T' && s != 'A'))
                return false;
            if ((s == 'C' && e != 'G') || (s == 'G' && e != 'C') || (e == 'C' && s != 'G') || (e == 'G' && s != 'C'))
                return false;
        }
        return true;
    }

    // ---- Nearest-neighbour salt-corrected Tm (PRIMER-NNTM-001, opt-in) ----------
    // SantaLucia & Hicks (2004) Watson-Crick NN ΔH°/ΔS° (1 M NaCl; Biopython DNA_NN4 — NOT the
    // SantaLucia 1998 / Allawi & SantaLucia 1997 set: AA/TT is −7.6/−21.3 here, −7.9/−22.2 in 1998,
    // and 2004 uses one duplex-initiation term plus a terminal A·T penalty instead of per-end
    // initiation) with the bimolecular Tm equation and published salt corrections. This is an OPT-IN
    // design Tm: the default CalculateMeltingTemperature (Wallace / Marmur-Doty) is unchanged.
    // All arithmetic is delegated to the canonical NN core ThermoConstants.CalculateNearestNeighborDuplex
    // (a line-by-line port of Biopython Tm_NN); this class only maps its API onto it:
    //   parameter set = SantaLuciaHicks2004 (DNA_NN4); R = 1.9872 (SantaLucia & Hicks 2004 Eq. 3;
    //   Biopython uses 1.987); C_T total → dnac1 = dnac2 = C_T/2 (k = C_T/4), self-complementary →
    //   dnac1 = C_T (k = C_T); self-complementarity detected from the sequence; salt mode → Biopython
    //   salt_correction method 0 / 5 / 6 / 7 with [Na⁺] (and, for method 7 only, [Mg²⁺], [dNTPs]).
    // Sources: SantaLucia J, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415, Table 1 + Eq. 3/5;
    //   Owczarzy R et al. (2004) Biochemistry 43:3537 (monovalent); Owczarzy R et al. (2008)
    //   Biochemistry 47:5336 (Mg²⁺/dNTP); Biopython 1.88 Bio.SeqUtils.MeltingTemp (reference).

    /// <summary>R of SantaLucia &amp; Hicks (2004) Eq. 3, used by the NN Tm methods of this class.</summary>
    private const double GasConstant = ThermoConstants.NnGasConstantSantaLuciaHicks2004;

    // Strand-concentration divisor x in Tm = ΔH°/(ΔS° + R·ln(C_T/x)) (SantaLucia & Hicks 2004 Eq. 3):
    // x = 4 for non-self-complementary, x = 1 for self-complementary duplexes (used by the dimer Tm).
    private const double NonSelfComplementaryFactor = 4.0;
    private const double SelfComplementaryFactor = 1.0;

    // Kelvin-to-Celsius offset.
    private const double KelvinOffset = 273.15;

    // Default total strand concentration C_T = 0.5 µM (a common PCR primer working
    // concentration). Exposed as a parameter; the caller may override.
    private const double DefaultStrandConcentrationMolar = 0.5e-6;

    // Default total strand concentration for the intermolecular dimer Tm. Primer3 / ntthal
    // uses dna_conc = 50 nM (thal.c default a->dna_conc = 50, lines 829/844); the opt-in dimer
    // methods adopt the same convention so they reproduce the ntthal reference out of the box.
    private const double DefaultDimerStrandConcentrationMolar = 50e-9;

    // The NN parameter set of this class's NN Tm, hairpin and dimer helpers (SantaLucia & Hicks 2004).
    private const NnParameterSet DesignNnParameterSet = NnParameterSet.SantaLuciaHicks2004;

    // SantaLucia & Hicks (2004) Table 1 initiation-type terms of DesignNnParameterSet
    // (duplex initiation +0.2/−5.7, terminal A·T +2.2/+6.9, symmetry 0/−1.4).
    private static readonly NnInitiationTerms DesignNnInitiation =
        ThermoConstants.GetNearestNeighborInitiation(DesignNnParameterSet);

    // Watson–Crick stack of the top-strand dinucleotide (upper-case ACGT) from DesignNnParameterSet.
    private static bool TryGetDesignStack(string dinucleotide, out (double DeltaH, double DeltaS) p) =>
        ThermoConstants.TryGetNearestNeighborStack(DesignNnParameterSet, dinucleotide, out p);

    /// <summary>Salt-correction mode for <see cref="CalculateMeltingTemperatureNN"/>.</summary>
    public enum SaltCorrectionMode
    {
        /// <summary>No correction — Tm at the SantaLucia 1 M NaCl reference state (Biopython saltcorr 0).</summary>
        None,

        /// <summary>
        /// SantaLucia (1998) / SantaLucia &amp; Hicks (2004) Eq. 5 entropy correction (Biopython saltcorr 5):
        /// ΔS°[Na] = ΔS°[1 M] + 0.368·(N − 1)·ln[Na⁺], N = oligo length (N − 1 = half the duplex phosphates);
        /// applied to ΔS° before the Tm equation.
        /// </summary>
        SantaLuciaEntropy,

        /// <summary>
        /// Owczarzy et al. (2004) monovalent 1/Tm correction (Biochemistry 43:3537; Biopython saltcorr 6):
        /// 1/Tm = 1/Tm(1 M) + (4.29·f(GC) − 3.95)·10⁻⁵·ln[Na⁺] + 9.40·10⁻⁶·ln²[Na⁺]. [Mg²⁺] and dNTPs are ignored.
        /// </summary>
        Owczarzy2004Monovalent,

        /// <summary>
        /// Owczarzy et al. (2008) divalent Mg²⁺ (and dNTP-adjusted) correction (Biochemistry 47:5336; Biopython
        /// saltcorr 7); reduces to the 2004 monovalent form when √[Mg²⁺]/[Na⁺] &lt; 0.22.
        /// </summary>
        Owczarzy2008Divalent
    }

    // SantaLucia (1998) / SantaLucia & Hicks (2004) Eq. 5 entropy salt-correction coefficient (0.368),
    // used by the Primer3 seqtm port and the dimer Tm.
    private const double SantaLuciaEntropySaltCoefficient = 0.368;

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a DNA oligonucleotide and its
    /// perfect complement with the SantaLucia &amp; Hicks (2004) nearest-neighbour parameters (Biopython
    /// <c>DNA_NN4</c>): duplex initiation, a terminal A·T penalty per A·T-closed end, the stacks and (for a
    /// self-complementary sequence) the symmetry correction. Equals Biopython
    /// <c>Tm_NN(seq, nn_table=DNA_NN4, selfcomp=…)</c>'s ΔH°/ΔS°; computed by
    /// <see cref="ThermoConstants.CalculateNearestNeighborThermodynamics"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3'; case-insensitive).</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the sequence is empty,
    /// shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        if (!IsAcgtOnly(seq))
            return null; // non-ACGT base present

        bool selfComp = IsSelfComplementary(seq);
        var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
            seq, parameterSet: DesignNnParameterSet, selfComplementary: selfComp, check: false);
        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a primer/oligonucleotide using the
    /// SantaLucia &amp; Hicks (2004) nearest-neighbour thermodynamics (Biopython <c>DNA_NN4</c>) and the
    /// bimolecular Tm equation, with an optional published salt correction. <b>Opt-in</b>: the default
    /// <see cref="CalculateMeltingTemperature(string)"/> (Wallace / Marmur-Doty) is unchanged.
    /// <para>
    /// Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15, with R = 1.9872 cal/(K·mol),
    /// x = 4 for a non-self-complementary duplex and x = 1 for a self-complementary one
    /// (SantaLucia &amp; Hicks 2004, Eq. 3); salt corrections per <paramref name="saltMode"/>.
    /// Identical to Biopython <c>Tm_NN(seq, nn_table=DNA_NN4, dnac1=dnac2=C_T/2 (self-complementary:
    /// dnac1=C_T), selfcomp, Na, Mg, dNTPs, saltcorr=0/5/6/7)</c> with Biopython's R = 1.987 replaced by
    /// 1.9872 (≈ +0.005 °C); computed by the canonical <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/>.
    /// </para>
    /// </summary>
    /// <param name="primer">DNA primer sequence (5'→3'). Must be ≥ 2 ACGT bases.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation ([Na⁺]+[K⁺]+[Tris]/2) concentration in
    /// mol/L (default 0.05 M = 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the
    /// <see cref="SaltCorrectionMode.Owczarzy2008Divalent"/> mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0; sequesters Mg²⁺
    /// in the divalent mode).</param>
    /// <param name="saltMode">Which salt correction to apply (default
    /// <see cref="SaltCorrectionMode.Owczarzy2004Monovalent"/>).</param>
    /// <returns>The nearest-neighbour Tm in °C, or <c>double.NaN</c> if the sequence is
    /// empty, shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static double CalculateMeltingTemperatureNN(
        string primer,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamics(primer);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, primer.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // Parameter-domain guards: the Tm equation takes R·ln(C_T/x) and the salt corrections take
    // ln([Na⁺]) / ln([Mg²⁺]); a zero or negative concentration would make ln undefined. Reject these
    // out-of-domain inputs explicitly (Biopython raises ValueError for them too).
    private static void ValidateNnConditions(
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar)
    {
        if (!(strandConcentrationMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(strandConcentrationMolar),
                strandConcentrationMolar, "Strand concentration C_T must be > 0 mol/L.");
        if (!(sodiumMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(sodiumMolar),
                sodiumMolar, "Monovalent cation concentration [Na⁺] must be > 0 mol/L.");
        if (!(magnesiumMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(magnesiumMolar),
                magnesiumMolar, "[Mg²⁺] must be ≥ 0 mol/L.");
        if (!(dntpMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMolar),
                dntpMolar, "Total dNTP concentration must be ≥ 0 mol/L.");
    }

    // Maps this class's NN Tm API onto the canonical Biopython Tm_NN core (see the region comment):
    // C_T → dnac1 = dnac2 = C_T/2 (k = C_T/4) or dnac1 = C_T for a self-complementary duplex (k = C_T);
    // mol/L → mM; salt mode → method 0/5/6/7 ([Mg²⁺]/[dNTPs] only for method 7); R = 1.9872.
    private static double NnTm(
        double deltaH, double deltaS, string saltSequence, bool selfComp,
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar,
        SaltCorrectionMode saltMode, double gasConstant = GasConstant)
    {
        double dnac = selfComp ? strandConcentrationMolar * 1e9 : strandConcentrationMolar * 1e9 / 2.0;
        var method = saltMode switch
        {
            SaltCorrectionMode.None => NnSaltCorrection.None,
            SaltCorrectionMode.SantaLuciaEntropy => NnSaltCorrection.SantaLucia1998Entropy,
            SaltCorrectionMode.Owczarzy2004Monovalent => NnSaltCorrection.Owczarzy2004,
            SaltCorrectionMode.Owczarzy2008Divalent => NnSaltCorrection.Owczarzy2008,
            _ => throw new ArgumentOutOfRangeException(nameof(saltMode), saltMode, "Unknown salt-correction mode.")
        };
        bool divalent = method == NnSaltCorrection.Owczarzy2008;
        try
        {
            return ThermoConstants.CalculateNearestNeighborTmFromThermodynamics(
                deltaH, deltaS, saltSequence, dnac1: dnac, dnac2: dnac, selfComplementary: selfComp,
                sodium: sodiumMolar * 1000.0,
                magnesium: divalent ? magnesiumMolar * 1000.0 : 0,
                dntps: divalent ? dntpMolar * 1000.0 : 0,
                saltCorrection: method, gasConstant: gasConstant).MeltingTemperature;
        }
        catch (ArgumentException)
        {
            return double.NaN; // degenerate duplex the reference rejects (Tm_NN ZeroDivisionError)
        }
    }

    // ---- NN internal-mismatch + dangling-end Tm (PRIMER-NNTM-001, opt-in extension) -------
    // Extends the perfect-match NN model to a probe–target duplex with internal mismatches,
    // terminal mismatches and/or a single unpaired dangling base at either end, exactly as
    // Biopython Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4, tmm_table=DNA_TMM1, imm_table=DNA_IMM1,
    // de_table=DNA_DE1) does (the core ThermoConstants.CalculateNearestNeighborDuplex is a port of it).
    //
    // Convention of this API: the top strand is 5'→3'; the bottom strand is supplied 3'→5' (the
    // complement of the top read in the SAME left-to-right order, NOT the reverse complement), so
    // column i pairs top[i] with bottom[i]. A '.' as the first or last character of a strand marks
    // the missing partner of a single dangling base on the other strand. It maps onto Tm_NN as
    // seq = top without '.', c_seq = bottom without '.', shift = (leading '.' of top) − (leading '.'
    // of bottom). Tm_NN then scores, in this order: dangling ends (Bommarito et al. 2000), terminal
    // mismatches (SantaLucia & Peyret 2001, Biopython DNA_TMM1), the initiation terms — the terminal
    // A·T penalty is taken from the first and last base of the TOP strand without '.' (Tm_NN's
    // `ends = seq[0] + seq[-1]`) — and every remaining stack from the internal-mismatch table
    // (Allawi & SantaLucia 1997/1998, Peyret et al. 1999) or the Watson-Crick table.

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a probe–target DNA duplex
    /// that may contain internal single mismatches (Allawi &amp; SantaLucia 1997/1998; Peyret et al. 1999),
    /// terminal mismatches (SantaLucia &amp; Peyret 2001) and/or a single dangling end at either end
    /// (Bommarito et al. 2000), with the SantaLucia &amp; Hicks (2004) Watson–Crick parameters. Identical to
    /// Biopython <c>Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4)</c>'s ΔH°/ΔS° (see the mapping in the
    /// remarks of <see cref="CalculateMeltingTemperatureNNMismatch"/>).
    /// </summary>
    /// <param name="topStrand">Top strand 5'→3'. May start/end with a single '.' marking a
    /// dangling base of the bottom strand.</param>
    /// <param name="bottomStrand">Bottom strand written 3'→5' (the complement of the top
    /// read left-to-right, NOT the reverse complement), so base i pairs with top base i.
    /// May start/end with a single '.' marking a dangling base of the top strand.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the strands are null,
    /// unequal length, shorter than two columns, contain a character other than A/C/G/T (case-insensitive)
    /// or a '.' that is not a single terminal marker facing a base, or contain a stack with no parameter
    /// (e.g. adjacent mismatches other than the tandem G·T motifs).</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsMismatch(
        string topStrand, string bottomStrand)
    {
        var mapped = MapMismatchDuplex(topStrand, bottomStrand);
        if (mapped is null)
            return null;

        var (seq, cSeq, shift, selfComp) = mapped.Value;
        try
        {
            var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
                seq, cSeq, shift, DesignNnParameterSet, selfComp, check: false, strict: true);
            return (dH, dS, selfComp);
        }
        catch (ArgumentException)
        {
            return null; // a neighbour pair with no thermodynamic parameter
        }
    }

    // Validates the column-aligned (top, bottom) pair and maps it to Tm_NN's (seq, c_seq, shift).
    private static (string Seq, string CSeq, int Shift, bool SelfComp)? MapMismatchDuplex(string top, string bottom)
    {
        if (top is null || bottom is null)
            return null;
        string t = top.ToUpperInvariant();
        string b = bottom.ToUpperInvariant();
        int n = t.Length;
        if (n != b.Length || n < 2)
            return null;

        for (int i = 0; i < n; i++)
        {
            char x = t[i], y = b[i];
            bool xDot = x == '.', yDot = y == '.';
            if ((xDot || yDot) && i != 0 && i != n - 1)
                return null; // a dangling marker must be terminal
            if (xDot && yDot)
                return null; // a column with no base at all
            if ((!xDot && x is not ('A' or 'C' or 'G' or 'T')) || (!yDot && y is not ('A' or 'C' or 'G' or 'T')))
                return null;
        }

        string seq = t.Replace(".", string.Empty, StringComparison.Ordinal);
        string cSeq = b.Replace(".", string.Empty, StringComparison.Ordinal);
        if (seq.Length < 1 || cSeq.Length < 1)
            return null;
        int shift = (t[0] == '.' ? 1 : 0) - (b[0] == '.' ? 1 : 0);

        // Symmetry term only for a fully paired, self-complementary duplex (Tm_NN's selfcomp flag).
        bool hasDangling = t.Contains('.') || b.Contains('.');
        bool selfComp = !hasDangling && IsSelfComplementary(t)
                        && string.Equals(b, Complement(t), StringComparison.Ordinal);
        return (seq, cSeq, shift, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) for a probe–target DNA duplex that may contain
    /// internal and terminal mismatches and/or a single dangling end at either end, with the
    /// SantaLucia &amp; Hicks (2004) Watson–Crick parameters plus the internal-mismatch, terminal-mismatch and
    /// Bommarito (2000) dangling-end terms (<see cref="CalculateNearestNeighborThermodynamicsMismatch"/>), the
    /// same bimolecular Tm equation and salt corrections as <see cref="CalculateMeltingTemperatureNN"/>.
    /// <b>Opt-in extension</b>: a fully paired duplex through this path equals
    /// <see cref="CalculateMeltingTemperatureNN"/>.
    /// </summary>
    /// <remarks>
    /// Equals Biopython <c>Tm_NN(seq, c_seq=…, shift=…, nn_table=DNA_NN4, dnac1=dnac2=C_T/2, Na, Mg, dNTPs,
    /// saltcorr=0/5/6/7)</c> with R = 1.9872, where seq/c_seq are the strands without '.', shift = (leading '.'
    /// of the top) − (leading '.' of the bottom). The salt correction uses the top strand without '.'
    /// (its length and GC fraction), as Tm_NN does.
    /// </remarks>
    /// <param name="topStrand">Top strand 5'→3' (may carry a leading/trailing '.' dangling-end marker).</param>
    /// <param name="bottomStrand">Bottom strand 3'→5', aligned base-for-base under the top
    /// (complement direction, NOT reverse complement; may carry a '.' dangling-end marker).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (see <see cref="CalculateNearestNeighborThermodynamicsMismatch"/>).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A non-positive C_T or [Na⁺], or a negative [Mg²⁺]/[dNTPs].</exception>
    public static double CalculateMeltingTemperatureNNMismatch(
        string topStrand,
        string bottomStrand,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamicsMismatch(topStrand, bottomStrand);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        string topPaired = topStrand.ToUpperInvariant().Replace(".", string.Empty, StringComparison.Ordinal);
        return NnTm(dH, dS, topPaired, selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // ---- LNA (locked nucleic acid)-modified NN thermodynamics + Tm (PROBE-LNATM-001) ----------
    // A DNA oligonucleotide with LNA monomers at given INTERNAL positions, hybridised to a DNA strand
    // (perfect complement or a supplied 3'→5' target), modelled by one of two published LNA·DNA
    // nearest-neighbour models, exactly as the MELTING 5 reference implementation (Dumousseau et al.
    // 2012, BMC Bioinformatics 13:101; melting5.jar 5.2.0 as shipped in Bioconductor rmelting:
    // LockedAcidNNMethod / McTigue04LockedAcid / Owczarzy11LockedAcid / Owczarzy11TandemLockedAcid /
    // Owczarzy11SingleMismatchLockedAcid, NearestNeighborMode.getAppropriatePatternModel) realises them:
    //
    //  * Base DNA model (both LNA models): SantaLucia (1998) unified = Allawi & SantaLucia (1997) stacks,
    //    initiation per terminal A·T (+2.3/+4.1) or G·C (+0.1/−2.8) pair (Biopython DNA_NN3, MELTING
    //    "all97", the DNA/DNA default of MELTING). McTigue et al. (2004) derived their LNA increments on
    //    this unified set; Owczarzy et al. (2011)'s single-LNA ΔH° equal those increments added to it
    //    (e.g. TTL/AA −5.574 = −7.9 + 2.326). Steps without LNA: Watson–Crick or internal-mismatch
    //    (DNA_IMM1) value, looked up as Biopython Tm_NN (ThermoConstants.TryGetNearestNeighborDuplexStep).
    //  * McTigue, Peterson & Kahn (2004) Biochemistry 43:5388 (MELTING "mct04"): for an ISOLATED LNA, each
    //    of its two flanking steps = DNA stack + ΔΔH°/ΔΔS° increment (32 increments, McTigue2004lockedmn.xml).
    //    McTigue parameterised single internal LNAs only; a run of ≥ 2 consecutive LNAs is handled — as in
    //    MELTING (getAppropriatePatternModel → tandem/single-mismatch LNA model whatever -lck says) — by the
    //    Owczarzy (2011) tables below.
    //  * Owczarzy, You, Groth & Tataurov (2011) Biochemistry 50:9352 (MELTING "owc11", its default): every
    //    step containing an LNA takes a COMPLETE NN value (no DNA base added): one LNA in the step → 32
    //    single-LNA parameters (Owczarzy2011lockedmn.xml); two LNAs, both pairs Watson–Crick → 16
    //    consecutive-LNA parameters (Owczarzy2011lockedTandemmn.xml); two LNAs, a mismatch opposite one of
    //    them → LNA-mismatch parameters (Owczarzy2011lockedmmn.xml). A mismatched LNA is therefore
    //    computable only when both of its neighbours are LNAs (the paper's +X+M+Y triplets).
    //  * Terminal LNAs and terminal mismatches are not parameterised (MELTING isApplicable → not computable).
    //  * Tm = 1000·ΔH° / (ΔS° + R·ln(C_T/x)) − 273.15 (x = 4, 1 if self-complementary) with the salt
    //    correction of saltMode (MELTING default for Na⁺ only = Owczarzy et al. 2004 Eq. 22 = Biopython
    //    method 6). MELTING's code uses R = 1.99 (NearestNeighborMode.computesMeltingTemperature); this
    //    class uses R = 1.9872 by default — pass gasConstant: 1.99 for bit-exact MELTING parity.
    //
    // Data anomalies kept exactly as in the reference implementation (the primary tables are behind
    // the ACS paywall; see docs/Validation/review-2026-09/B07.md F26): the mismatch file has no
    // "GLAL/CA" entry but a double-mismatch key "CLAL/CA", and lists "GLTL/TA" twice (MELTING's
    // HashMap keeps the last: −14.213/−40.041); an independent transcription (iCarrin/Bio_dpt
    // lna_tm.py) has ΔS° −21.535 for consecutive +T+C where MELTING has −21.735 (MELTING kept).
    // Reference cross-check (melting5.jar, Na = 1 M): CCATT(L)GCTACC, 1e-4 M — mct04 63.61426 °C,
    // owc11 63.48299 °C; GA(L)C(L)C 12.94323 °C (rmelting test-method.locked.R) — reproduced bit-exactly.

    /// <summary>Published LNA·DNA nearest-neighbour models for <see cref="CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?, double, double, double, double, SaltCorrectionMode, double)"/>.</summary>
    public enum LnaNearestNeighborModel
    {
        /// <summary>
        /// Owczarzy, You, Groth &amp; Tataurov (2011) Biochemistry 50:9352 — complete NN parameters for steps with
        /// one LNA, two consecutive LNAs, and an LNA·DNA mismatch inside an LNA triplet (MELTING "owc11", its default).
        /// </summary>
        Owczarzy2011,

        /// <summary>
        /// McTigue, Peterson &amp; Kahn (2004) Biochemistry 43:5388 — ΔΔH°/ΔΔS° increments for an isolated internal
        /// LNA added to the SantaLucia (1998) unified DNA stacks (MELTING "mct04"); runs of consecutive LNAs use
        /// the <see cref="Owczarzy2011"/> parameters, as MELTING does.
        /// </summary>
        McTigue2004
    }

    // The base DNA NN set of both LNA models: SantaLucia (1998) unified / Allawi & SantaLucia (1997).
    private const NnParameterSet LnaBaseParameterSet = NnParameterSet.AllawiSantaLucia1997;

    // LNA parameter tables (ΔH° kcal/mol, ΔS° cal/(K·mol), 1 M Na⁺), transcribed from the MELTING 5.2.0
    // data files (cal/mol ÷ 1000). Key = MELTING notation: top strand 5'→3' with "L" after each locked
    // base / the two opposite bottom-strand bases 3'→5' (e.g. "TTL/AA": step TT, 3' base locked).
    // McTigue (2004): increments ΔΔH°/ΔΔS° (McTigue2004lockedmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> McTigue2004LnaIncrements = new(StringComparer.Ordinal)
    {
        ["ALA/TT"] = (0.707, 2.5), ["ALT/TA"] = (2.282, 7.5),
        ["ALG/TC"] = (0.264, 2.6), ["ALC/TG"] = (1.131, 4.1),
        ["TLA/AT"] = (-0.046, 1.6), ["TLT/AA"] = (1.528, 5.3),
        ["TLG/AC"] = (-1.540, -3), ["TLC/AG"] = (1.893, 6.7),
        ["GLA/CT"] = (3.162, 10.5), ["GLT/CA"] = (-0.212, 0.1),
        ["GLG/CC"] = (-2.844, -6.7), ["GLC/CG"] = (-0.360, -0.3),
        ["CLA/GT"] = (1.049, 4.3), ["CLT/GA"] = (0.708, 4.2),
        ["CLG/GC"] = (0.785, 3.7), ["CLC/GG"] = (2.096, 8),
        ["AAL/TT"] = (0.992, 4.1), ["ATL/TA"] = (1.816, 6.9),
        ["AGL/TC"] = (-1.200, -1.8), ["ACL/TG"] = (2.890, 10.6),
        ["TAL/AT"] = (1.591, 5.3), ["TTL/AA"] = (2.326, 8.1),
        ["TGL/AC"] = (2.165, 7.2), ["TCL/AG"] = (0.609, 3.2),
        ["GAL/CT"] = (0.444, 2.9), ["GTL/CA"] = (-0.635, -0.3),
        ["GGL/CC"] = (-0.943, -0.9), ["GCL/CG"] = (-0.925, -1.1),
        ["CAL/GT"] = (1.358, 4.4), ["CTL/GA"] = (-1.671, -4.1),
        ["CGL/GC"] = (-0.276, -0.7), ["CCL/GG"] = (2.063, 7.6)
    };

    // Owczarzy (2011): complete NN parameters of a step with one LNA (Owczarzy2011lockedmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaSingle = new(StringComparer.Ordinal)
    {
        ["ALA/TT"] = (-7.193, -19.723), ["ALC/TG"] = (-7.269, -18.336),
        ["ALG/TC"] = (-7.536, -18.387), ["ALT/TA"] = (-4.918, -12.943),
        ["CLA/GT"] = (-7.451, -18.38), ["CLC/GG"] = (-5.904, -11.904),
        ["CLG/GC"] = (-9.815, -23.491), ["CLT/GA"] = (-7.092, -16.825),
        ["GLA/CT"] = (-5.038, -11.656), ["GLC/CG"] = (-10.160, -24.651),
        ["GLG/CC"] = (-10.844, -26.58), ["GLT/CA"] = (-8.612, -22.327),
        ["TLA/AT"] = (-7.246, -19.738), ["TLC/AG"] = (-6.307, -15.515),
        ["TLG/AC"] = (-10.040, -25.744), ["TLT/AA"] = (-6.372, -16.902),
        ["AAL/TT"] = (-6.908, -18.135), ["ACL/TG"] = (-5.510, -11.824),
        ["AGL/TC"] = (-9.000, -22.826), ["ATL/TA"] = (-5.384, -13.537),
        ["CAL/GT"] = (-7.142, -18.333), ["CCL/GG"] = (-5.937, -12.335),
        ["CGL/GC"] = (-10.876, -27.918), ["CTL/GA"] = (-9.471, -25.07),
        ["GAL/CT"] = (-7.756, -19.302), ["GCL/CG"] = (-10.725, -25.511),
        ["GGL/CC"] = (-8.943, -20.833), ["GTL/CA"] = (-9.035, -22.742),
        ["TAL/AT"] = (-5.609, -16.019), ["TCL/AG"] = (-7.591, -19.031),
        ["TGL/AC"] = (-6.335, -15.537), ["TTL/AA"] = (-5.574, -14.149)
    };

    // Owczarzy (2011): two consecutive LNAs, Watson–Crick (Owczarzy2011lockedTandemmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaConsecutive = new(StringComparer.Ordinal)
    {
        ["ALAL/TT"] = (-9.991, -27.175), ["ALCL/TG"] = (-11.389, -28.963),
        ["ALGL/TC"] = (-12.793, -31.607), ["ALTL/TA"] = (-14.703, -40.75),
        ["CLAL/GT"] = (-14.177, -35.498), ["CLCL/GG"] = (-15.399, -36.375),
        ["CLGL/GC"] = (-14.558, -35.239), ["CLTL/GA"] = (-15.737, -41.218),
        ["GLAL/CT"] = (-13.959, -35.097), ["GLCL/CG"] = (-16.109, -40.738),
        ["GLGL/CC"] = (-13.022, -29.673), ["GLTL/CA"] = (-17.361, -45.858),
        ["TLAL/AT"] = (-10.318, -26.108), ["TLCL/AG"] = (-9.166, -21.735),
        ["TLGL/AC"] = (-10.046, -22.591), ["TLTL/AA"] = (-10.419, -27.683)
    };

    // Owczarzy (2011): two consecutive LNAs with a mismatch opposite one of them (Owczarzy2011lockedmmn.xml;
    // 98 entries, "GLTL/TA" listed twice — the later value is kept, as MELTING's HashMap does).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaMismatch = new(StringComparer.Ordinal)
    {
        ["ALAL/AT"] = (-3.826, -13.109), ["ALCL/AG"] = (-2.367, -7.322),
        ["ALGL/AC"] = (-4.849, -13.007), ["ALTL/AA"] = (-5.049, -17.514),
        ["ALAL/TA"] = (-4.229, -15.16), ["CLAL/GA"] = (-5.878, -17.663),
        ["CLAL/CA"] = (-8.558, -23.976), ["TLAL/AA"] = (2.074, 3.446),
        ["CLAL/CT"] = (2.218, 4.75), ["CLCL/CG"] = (1.127, 1.826),
        ["CLGL/CC"] = (-10.903, -32.025), ["CLTL/CA"] = (-2.053, -10.517),
        ["ALCL/TC"] = (1.065, -1.403), ["CLCL/GC"] = (-9.522, -27.024),
        ["GLCL/CC"] = (-4.767, -14.897), ["TLCL/AC"] = (4.114, 9.258),
        ["GLAL/GT"] = (-2.920, -9.387), ["GLCL/GG"] = (-8.139, -21.784),
        ["GLGL/GC"] = (-5.149, -12.508), ["GLTL/GA"] = (-8.991, -27.311),
        ["ALGL/TG"] = (-4.980, -15.426), ["CLGL/GG"] = (-4.441, -12.158),
        ["GLGL/CG"] = (-13.505, -36.021), ["TLGL/AG"] = (-2.775, -9.286),
        ["TLAL/TT"] = (-3.744, -12.149), ["TLCL/TG"] = (-4.387, -13.52),
        ["TLGL/TC"] = (-6.346, -16.629), ["TLTL/TA"] = (-7.697, -25.049),
        ["ALTL/TT"] = (-4.207, -14.307), ["CLTL/GT"] = (-8.176, -22.962),
        ["GLTL/CT"] = (-7.241, -20.622), ["TLTL/AT"] = (-2.051, -7.055),
        ["ALAL/CT"] = (-1.362, -5.551), ["ALCL/CG"] = (-1.759, -6.511),
        ["ALGL/CC"] = (-6.549, -18.073), ["ALTL/CA"] = (-3.563, -14.105),
        ["ALAL/TC"] = (-2.078, -10.088), ["CLAL/GC"] = (-5.868, -16.952),
        ["GLAL/CC"] = (-8.477, -24.565), ["TLAL/AC"] = (2.690, 4.965),
        ["CLAL/AT"] = (-9.844, -29.673), ["CLCL/AG"] = (-3.761, -11.204),
        ["CLGL/AC"] = (-9.845, -27.316), ["CLTL/AA"] = (-3.389, -12.517),
        ["ALCL/TA"] = (0.753, -0.503), ["CLCL/GA"] = (-12.714, -35.555),
        ["GLCL/CA"] = (-12.658, -35.729), ["TLCL/AA"] = (-1.719, -7.023),
        ["ALAL/GT"] = (2.193, 4.374), ["ALCL/GG"] = (-8.453, -22.672),
        ["ALGL/GC"] = (-1.164, -2.532), ["ALTL/GA"] = (-7.418, -24.066),
        ["ALAL/TG"] = (-1.963, -9.013), ["CLAL/GG"] = (-8.712, -23.779),
        ["GLAL/CG"] = (-7.875, -21.661), ["TLAL/AG"] = (3.207, 7.156),
        ["GLAL/AT"] = (-2.914, -9.402), ["GLCL/AG"] = (-9.131, -25.347),
        ["GLGL/AC"] = (-2.154, -3.871), ["GLTL/AA"] = (-8.515, -26.313),
        ["ALGL/TA"] = (-6.691, -21.148), ["CLGL/GA"] = (-3.960, -10.588),
        ["GLGL/CA"] = (-12.898, -34.656), ["TLGL/AA"] = (0.334, -0.44),
        ["CLAL/TT"] = (0.382, -0.579), ["CLCL/TG"] = (-2.716, -8),
        ["CLGL/TC"] = (-10.363, -29.315), ["CLTL/TA"] = (-5.783, -20.173),
        ["ALCL/TT"] = (-0.692, -5.278), ["CLCL/GT"] = (-10.288, -28.503),
        ["GLCL/CT"] = (-9.062, -26.356), ["TLCL/AT"] = (2.073, 3.968),
        ["TLAL/CT"] = (-5.485, -17.347), ["TLCL/CG"] = (1.451, 1.556),
        ["TLGL/CC"] = (-7.213, -20.128), ["TLTL/CA"] = (-2.397, -11.371),
        ["ALTL/TC"] = (-0.633, -5.801), ["CLTL/GC"] = (-6.868, -21),
        ["GLTL/CC"] = (-5.853, -16.643), ["TLTL/AC"] = (0.211, -1.446),
        ["GLAL/TT"] = (-5.551, -15.398), ["GLCL/TG"] = (-14.943, -40.148),
        ["GLGL/TC"] = (-8.110, -18.349), ["GLTL/TA"] = (-14.213, -40.041),
        ["ALGL/TT"] = (-7.130, -20.786), ["CLGL/GT"] = (-14.862, -39.43),
        ["GLGL/CT"] = (-14.622, -37.51), ["TLGL/AT"] = (-6.703, -18.111),
        ["TLAL/GT"] = (-4.612, -14.039), ["TLCL/GG"] = (-9.798, -26.406),
        ["TLGL/GC"] = (-4.519, -11.065), ["TLTL/GA"] = (-4.523, -15.693),
        ["ALTL/TG"] = (-2.364, -8.834), ["CLTL/GG"] = (-11.396, -30.732),
        ["GLTL/CG"] = (-6.233, -15.933), ["TLTL/AG"] = (-2.960, -9.305)
    };

    /// <summary>
    /// Duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M Na⁺) of a DNA oligonucleotide with <b>internal</b> LNA
    /// monomers at <paramref name="lnaPositions"/> paired with its perfect DNA complement, by the default
    /// <see cref="LnaNearestNeighborModel.Owczarzy2011"/> LNA·DNA nearest-neighbour model on the SantaLucia (1998)
    /// unified DNA parameters (MELTING 5 default). See
    /// <see cref="CalculateNearestNeighborThermodynamicsLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'; case-insensitive), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based LNA positions (order/duplicates tolerated); terminal or out-of-range → not computable.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> when not computable.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions) =>
        CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions, LnaNearestNeighborModel.Owczarzy2011);

    /// <summary>
    /// Duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M Na⁺) of a DNA oligonucleotide carrying <b>internal</b> LNA
    /// (locked nucleic acid) monomers, hybridised to a DNA strand, by a published LNA·DNA nearest-neighbour model
    /// (see <see cref="LnaNearestNeighborModel"/>) on the SantaLucia (1998) unified DNA parameters (Allawi &amp;
    /// SantaLucia 1997; initiation per terminal A·T / G·C pair) — the MELTING 5 implementation, reproduced
    /// bit-exactly. Steps without an LNA use the Watson–Crick or (for an internal DNA mismatch) the Allawi /
    /// SantaLucia / Peyret internal-mismatch parameters.
    /// </summary>
    /// <param name="sequence">LNA-modified oligonucleotide (DNA letters, 5'→3'; case-insensitive), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the LNA monomers in <paramref name="sequence"/>; order and
    /// duplicates are tolerated. A terminal (0 or length − 1) or out-of-range position is not parameterised.</param>
    /// <param name="model">LNA nearest-neighbour model.</param>
    /// <param name="target">The DNA strand opposite <paramref name="sequence"/>, written 3'→5' (base i pairs with
    /// sequence[i]; same length); <c>null</c> = the perfect complement. Internal mismatches are allowed opposite a DNA
    /// base (DNA_IMM1) or opposite the central LNA of three consecutive LNAs (Owczarzy 2011); the two terminal pairs
    /// must be Watson–Crick.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> when not computable: empty / &lt; 2 bases / non-ACGT
    /// strand, target of another length, terminal or out-of-range LNA, terminal mismatch, or a step with no
    /// published parameter (e.g. a mismatch opposite an isolated LNA, any mismatch with
    /// <see cref="LnaNearestNeighborModel.McTigue2004"/> at an isolated LNA).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="lnaPositions"/> is null.</exception>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        LnaNearestNeighborModel model,
        string? target = null)
    {
        ArgumentNullException.ThrowIfNull(lnaPositions);
        if (model is not (LnaNearestNeighborModel.Owczarzy2011 or LnaNearestNeighborModel.McTigue2004))
            throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown LNA nearest-neighbour model.");
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        if (!IsAcgtOnly(seq))
            return null;
        string bottom = target is null ? Complement(seq) : target.ToUpperInvariant();
        if (bottom.Length != seq.Length || !IsAcgtOnly(bottom))
            return null;

        int n = seq.Length;
        var locked = new bool[n];
        foreach (int pos in lnaPositions)
        {
            // Terminal LNAs are not parameterised by either model (MELTING isApplicable).
            if (pos <= 0 || pos >= n - 1)
                return null;
            locked[pos] = true;
        }

        // Terminal pairs must be Watson–Crick (no terminal-mismatch LNA parameters).
        if (!IsWatsonCrickPair(seq[0], bottom[0]) || !IsWatsonCrickPair(seq[n - 1], bottom[n - 1]))
            return null;

        // An LNA-modified strand paired with an unmodified DNA strand is never a symmetric duplex (the two
        // strands differ chemically), so only an LNA-free self-complementary sequence with its own complement
        // is treated as self-complementary (symmetry term, x = 1); MELTING likewise rejects -self with LNAs.
        bool selfComp = target is null && !locked.Contains(true) && IsSelfComplementary(seq);
        var init = ThermoConstants.GetNearestNeighborInitiation(LnaBaseParameterSet);
        double dH = init.Initiation.DeltaH, dS = init.Initiation.DeltaS;
        var oneOrAll = seq.Any(c => c is 'G' or 'C') ? init.OneGC : init.AllAT;
        dH += oneOrAll.DeltaH; dS += oneOrAll.DeltaS;
        foreach (char end in new[] { seq[0], seq[n - 1] })
        {
            var t = end is 'A' or 'T' ? init.TerminalAT : init.TerminalGC;
            dH += t.DeltaH; dS += t.DeltaS;
        }
        if (seq[0] == 'T') { dH += init.FiveTerminalTA.DeltaH; dS += init.FiveTerminalTA.DeltaS; }
        if (seq[n - 1] == 'A') { dH += init.FiveTerminalTA.DeltaH; dS += init.FiveTerminalTA.DeltaS; }

        for (int i = 0; i < n - 1; i++)
        {
            string top = seq.Substring(i, 2);
            string bot = bottom.Substring(i, 2);
            int lockedCount = (locked[i] ? 1 : 0) + (locked[i + 1] ? 1 : 0);

            if (lockedCount == 0)
            {
                if (!ThermoConstants.TryGetNearestNeighborDuplexStep(LnaBaseParameterSet, top, bot, out var p))
                    return null;
                dH += p.DeltaH; dS += p.DeltaS;
                continue;
            }

            bool isolated = (locked[i] && IsIsolatedLna(locked, i)) || (locked[i + 1] && IsIsolatedLna(locked, i + 1));
            if (model == LnaNearestNeighborModel.McTigue2004 && isolated)
            {
                // McTigue (2004): DNA stack + increment (Watson–Crick steps only).
                if (!IsWatsonCrickPair(top[0], bot[0]) || !IsWatsonCrickPair(top[1], bot[1])
                    || !ThermoConstants.TryGetNearestNeighborStack(LnaBaseParameterSet, top, out var stack)
                    || !McTigue2004LnaIncrements.TryGetValue(LnaKey(top, bot, locked[i], locked[i + 1]), out var inc))
                    return null;
                dH += stack.DeltaH + inc.DeltaH;
                dS += stack.DeltaS + inc.DeltaS;
                continue;
            }

            // Owczarzy (2011): complete parameter of the LNA-containing step.
            bool watsonCrick = IsWatsonCrickPair(top[0], bot[0]) && IsWatsonCrickPair(top[1], bot[1]);
            var table = Owczarzy2011LnaSingle;
            if (lockedCount == 2)
                table = watsonCrick ? Owczarzy2011LnaConsecutive : Owczarzy2011LnaMismatch;
            if (!table.TryGetValue(LnaKey(top, bot, locked[i], locked[i + 1]), out var v))
                return null;
            dH += v.DeltaH; dS += v.DeltaS;
        }

        if (selfComp)
        {
            dH += init.Symmetry.DeltaH; dS += init.Symmetry.DeltaS;
        }
        return (dH, dS, selfComp);

        static bool IsIsolatedLna(bool[] l, int p) =>
            (p == 0 || !l[p - 1]) && (p == l.Length - 1 || !l[p + 1]);

        static string LnaKey(string top, string bot, bool lock0, bool lock1) =>
            string.Concat(top[0].ToString(), lock0 ? "L" : "", top[1].ToString(), lock1 ? "L" : "", "/", bot);
    }

    /// <summary>
    /// LNA-adjusted nearest-neighbour Tm (°C) of a DNA oligonucleotide with <b>internal</b> LNA monomers paired with its
    /// perfect DNA complement, by the default <see cref="LnaNearestNeighborModel.Owczarzy2011"/> model (MELTING 5
    /// default) at the stated conditions. See
    /// <see cref="CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?, double, double, double, double, SaltCorrectionMode, double)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based internal LNA positions.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction (default Owczarzy2004Monovalent, as MELTING for Na⁺ only).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when not computable.</returns>
    /// <exception cref="ArgumentOutOfRangeException">C_T ≤ 0, [Na⁺] ≤ 0, [Mg²⁺] &lt; 0 or [dNTP] &lt; 0.</exception>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent) =>
        CalculateMeltingTemperatureNNLna(sequence, lnaPositions, LnaNearestNeighborModel.Owczarzy2011, null,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);

    /// <summary>
    /// LNA-adjusted nearest-neighbour melting temperature (°C) of an LNA-modified DNA oligonucleotide hybridised to a
    /// DNA strand: ΔH°/ΔS° from
    /// <see cref="CalculateNearestNeighborThermodynamicsLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?)"/>,
    /// Tm = 1000·ΔH° / (ΔS° + R·ln(C_T/x)) − 273.15 (x = 4; 1 for a self-complementary duplex) with the salt correction
    /// of <paramref name="saltMode"/> evaluated on the DNA letters of <paramref name="sequence"/> (the canonical
    /// <see cref="ThermoConstants.CalculateNearestNeighborTmFromThermodynamics"/>). With
    /// <c>gasConstant: 1.99</c> and <see cref="SaltCorrectionMode.Owczarzy2004Monovalent"/> this equals MELTING 5
    /// (<c>-H dnadna -P C_T -E Na=…</c>, <c>-lck mct04|owc11</c>) exactly.
    /// </summary>
    /// <param name="sequence">LNA-modified oligonucleotide (DNA letters, 5'→3'), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based internal LNA positions.</param>
    /// <param name="model">LNA nearest-neighbour model.</param>
    /// <param name="target">Opposite DNA strand 3'→5' (same length), or <c>null</c> for the perfect complement.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction (default Owczarzy2004Monovalent).</param>
    /// <param name="gasConstant">R in cal/(K·mol) (default 1.9872, SantaLucia &amp; Hicks 2004; MELTING uses 1.99).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when not computable.</returns>
    /// <exception cref="ArgumentOutOfRangeException">C_T ≤ 0, [Na⁺] ≤ 0, [Mg²⁺] &lt; 0, [dNTP] &lt; 0, R ≤ 0 or an unknown model.</exception>
    /// <exception cref="ArgumentNullException"><paramref name="lnaPositions"/> is null.</exception>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        LnaNearestNeighborModel model,
        string? target = null,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent,
        double gasConstant = GasConstant)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);
        if (!(gasConstant > 0) || double.IsInfinity(gasConstant))
            throw new ArgumentOutOfRangeException(nameof(gasConstant), gasConstant, "R must be a positive finite number.");

        var thermo = CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions, model, target);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, sequence.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode, gasConstant);
    }

    // ---- DNA hairpin folding + secondary-structure (hairpin) Tm (PRIMER-TM-001, opt-in) ----
    // Finds the most stable intramolecular hairpin (a single stem closing one hairpin loop)
    // of a DNA oligo and computes its ΔH°/ΔS°/ΔG°37 and unimolecular melting temperature.
    // The perfect-match CalculateMeltingTemperatureNN and the default Wallace/Marmur-Doty Tm
    // are UNCHANGED; this is a new opt-in capability.
    //
    // Model (SantaLucia & Hicks 2004, Annu Rev Biophys 33:415, "Hairpin Loops", Eqs 8–11):
    //   ΔG°37(hairpin) = Σ stem NN stacks (Table 1)  +  ΔG°37(hairpin loop of N) (Table 4)
    //   "To compute the stability of a complete hairpin + stem, one simply adds the
    //    salt-corrected base pair NN contributions (Table 1; Equation 3) to the loop energy
    //    from Equations 8–10." (paper, p.428).
    //   The stem contributes only its nearest-neighbour STACKS (Table 1 propagation terms) —
    //   the bimolecular duplex-initiation term (+0.2/−5.7) is a TWO-strand nucleation cost and
    //   does NOT apply to a unimolecular hairpin (the loop-initiation term is the nucleation
    //   cost instead). Loop ΔH° = 0 for every loop size, and the loop ΔS° increment is
    //   ΔS° = −ΔG°37 × 1000 / 310.15 (Table 4 note: "ΔS° = ΔG°37 × 1000/310.15"; the loop is
    //   destabilising so ΔG°37 > 0 → ΔS° < 0).
    //   Two-state hairpin Tm (Eq 11) is UNIMOLECULAR/concentration-independent:
    //       Tm = ΔH° × 1000 / ΔS° − 273.15        (NO R·ln(C_T/x) strand-concentration term).
    //
    // NOT bundled (honest residual): the supplementary triloop/tetraloop bonus tables (length-3
    // and length-4 special loops) and the terminal-mismatch increment (the first mismatch stack
    // closing loops of length ≥4) are separate Annual-Reviews supplementary tables not embedded
    // here. They are exposed as an OPT-IN caller-supplied additive ΔG°37/ΔH° adjustment
    // (default 0) so a caller who has those tables can supply the increment; without it the
    // result is the stem-stack + loop-initiation core, which is exact and fully sourced.

    /// <summary>Hairpin loop ΔG°37 increment (kcal/mol, 1 M NaCl) by loop size (number of
    /// unpaired loop nucleotides). SantaLucia &amp; Hicks (2004) Table 4 "Hairpin loops" column;
    /// sizes 3–30 are tabulated. ΔH° = 0 for all sizes; ΔS° = −ΔG°37·1000/310.15.</summary>
    private static readonly Dictionary<int, double> HairpinLoopInitiationDeltaG = new()
    {
        [3] = 3.5, [4] = 3.5, [5] = 3.3, [6] = 4.0, [7] = 4.2, [8] = 4.3, [9] = 4.5,
        [10] = 4.6, [12] = 5.0, [14] = 5.1, [16] = 5.3, [18] = 5.5, [20] = 5.7,
        [25] = 6.1, [30] = 6.3
    };

    // SantaLucia & Hicks (2004): minimum sterically allowed hairpin loop size is 3 nt
    // ("Hairpin loops with lengths shorter than 3 are sterically prohibited.").
    private const int MinHairpinLoopSize = 3;

    // Reference temperature for the ΔG°37 ↔ ΔS° conversion (310.15 K = 37 °C).
    // SantaLucia & Hicks (2004) Table 4 note: ΔS° = ΔG°37 × 1000/310.15.
    private const double ReferenceTemperatureKelvin = 310.15;

    // Jacobson-Stockmayer entropic extrapolation coefficient for loop sizes beyond the
    // tabulated lengths. SantaLucia & Hicks (2004) Eq. 7: ΔG°37(loop-n) =
    // ΔG°37(loop-x) + 2.44·R·310.15·ln(n/x); the 2.44 coefficient is from recent DNA
    // kinetics measurements (ref 22), preferred over the older 1.75.
    private const double JacobsonStockmayerCoefficient = 2.44;

    /// <summary>
    /// Hairpin loop ΔG°37 (kcal/mol) for a loop of <paramref name="loopSize"/> unpaired
    /// nucleotides. Tabulated sizes return Table 4 directly; non-tabulated sizes are filled by
    /// the Jacobson-Stockmayer extrapolation from the largest tabulated size ≤ n
    /// (SantaLucia &amp; Hicks 2004 Eq. 7). Loop sizes &lt; 3 are sterically prohibited.
    /// </summary>
    private static double HairpinLoopDeltaG(int loopSize)
    {
        if (HairpinLoopInitiationDeltaG.TryGetValue(loopSize, out double dg))
            return dg;

        // Jacobson-Stockmayer from the largest tabulated x ≤ loopSize.
        int x = 0;
        foreach (int size in HairpinLoopInitiationDeltaG.Keys)
            if (size <= loopSize && size > x) x = size;

        return HairpinLoopInitiationDeltaG[x]
               + JacobsonStockmayerCoefficient * GasConstant * ReferenceTemperatureKelvin
                 * 1e-3 * Math.Log((double)loopSize / x);
    }

    private static bool IsWatsonCrickPair(char a, char b) =>
        (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
        (a == 'G' && b == 'C') || (a == 'C' && b == 'G');

    /// <summary>
    /// The most stable intramolecular DNA hairpin found in an oligo: the closing stem span,
    /// stem length (base pairs), loop size, and the hairpin ΔH° (kcal/mol), ΔS° (cal/(K·mol)),
    /// and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="StemStart">5'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemEnd">3'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemLength">Number of base pairs in the stem.</param>
    /// <param name="LoopSize">Number of unpaired loop nucleotides closed by the stem.</param>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol.</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol).</param>
    /// <param name="DeltaG37">Hairpin ΔG°37 in kcal/mol (negative = stable).</param>
    public readonly record struct HairpinResult(
        int StemStart, int StemEnd, int StemLength, int LoopSize,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (minimum ΔG°37) intramolecular DNA <b>hairpin</b> — a single
    /// Watson-Crick stem closing one hairpin loop — in <paramref name="sequence"/>, using the
    /// SantaLucia &amp; Hicks (2004) Table 1 nearest-neighbour stem stacks and their
    /// Table 4 hairpin-loop initiation increments. <b>Opt-in</b>: the duplex Tm methods
    /// are unchanged. Returns <c>null</c> when the sequence is empty, contains a non-ACGT
    /// character, or admits no hairpin at all (no stem of ≥ 2 bp can close a loop of ≥ 3 nt,
    /// e.g. a homopolymer such as poly-A).
    /// <para>
    /// Model: ΔG°37 = Σ stem NN stacks (Table 1) + ΔG°37(loop of N) (Table 4); the bimolecular
    /// duplex-initiation term is intentionally excluded for this unimolecular structure. Loop
    /// ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. The supplementary triloop/tetraloop and
    /// terminal-mismatch increments are not bundled (see <paramref name="loopBonusDeltaG37"/>).
    /// </para>
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2 → at least one
    /// NN stack). Must be ≥ 2.</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied additive ΔG°37 increment
    /// (kcal/mol) for the terminal-mismatch / special triloop-tetraloop bonus that is NOT
    /// bundled (default 0). Added to the loop free energy; its ΔS° contribution follows the
    /// same −ΔG·1000/310.15 rule (ΔH° contribution 0), consistent with the Table 4 loop model.</param>
    /// <returns>The most stable hairpin, or <c>null</c> if none exists / invalid input.</returns>
    public static HairpinResult? FindMostStableHairpin(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        if (string.IsNullOrEmpty(sequence) || minStemLength < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        int n = seq.Length;
        foreach (char c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T'))
                return null; // non-ACGT base present

        HairpinResult? best = null;

        // For every candidate outermost closing pair (i, j), extend the stem inward as far as
        // Watson-Crick pairing allows, then close the remaining inner bases as a hairpin loop.
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (!IsWatsonCrickPair(seq[i], seq[j]))
                    continue;

                double dH = 0.0, dS = 0.0;
                // Extend the stem from the outermost pair (i, j) inward.
                int a = i, b = j;
                int stemPairs = 0;
                while (a < b && IsWatsonCrickPair(seq[a], seq[b]))
                {
                    if (stemPairs > 0)
                    {
                        // NN stack between pair (a-1, b+1) and (a, b): key = 5'-strand dinucleotide seq[a-1..a].
                        string step = seq.Substring(a - 1, 2);
                        if (!TryGetDesignStack(step, out var p))
                            break;
                        dH += p.DeltaH;
                        dS += p.DeltaS;
                    }
                    stemPairs++;
                    a++;
                    b--;
                }

                if (stemPairs < minStemLength)
                    continue;

                // Innermost pair is (a-1, b+1); loop is the bases strictly between them.
                int innerLeft = a - 1;
                int innerRight = b + 1;
                int loopSize = innerRight - innerLeft - 1;
                if (loopSize < MinHairpinLoopSize)
                    continue;

                double loopDg = HairpinLoopDeltaG(loopSize) + loopBonusDeltaG37;
                // Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (destabilising loop).
                double loopDs = -loopDg * 1000.0 / ReferenceTemperatureKelvin;

                double totalDh = dH;                       // loop ΔH° contribution is 0
                double totalDs = dS + loopDs;
                double dG37 = totalDh - ReferenceTemperatureKelvin * totalDs / 1000.0;

                if (best is null || dG37 < best.Value.DeltaG37)
                    best = new HairpinResult(i, j, stemPairs, loopSize, totalDh, totalDs, dG37);
            }
        }

        return best;
    }

    /// <summary>
    /// Computes the secondary-structure (hairpin) melting temperature (°C) of a DNA oligo:
    /// finds its most stable intramolecular hairpin (<see cref="FindMostStableHairpin"/>) and
    /// returns the <b>unimolecular</b> two-state Tm = ΔH°·1000/ΔS° − 273.15
    /// (SantaLucia &amp; Hicks 2004, Eq. 11). A hairpin is intramolecular, so the Tm is
    /// concentration-independent: there is <b>no</b> R·ln(C_T/x) strand-concentration term.
    /// <b>Opt-in</b>: the duplex Tm methods (<see cref="CalculateMeltingTemperatureNN"/>) and the
    /// default <see cref="CalculateMeltingTemperature(string)"/> are unchanged.
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2).</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied terminal-mismatch / special-loop
    /// ΔG°37 increment (default 0; not bundled — see <see cref="FindMostStableHairpin"/>).</param>
    /// <returns>The hairpin Tm in °C, or <c>double.NaN</c> if no hairpin exists / invalid input.</returns>
    public static double CalculateHairpinMeltingTemperature(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        var hairpin = FindMostStableHairpin(sequence, minStemLength, loopBonusDeltaG37);
        if (hairpin is null)
            return double.NaN;

        var h = hairpin.Value;
        // Unimolecular: NO concentration term (Eq. 11).
        return (h.DeltaH * 1000.0) / h.DeltaS - KelvinOffset;
    }

    // ---- Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment ----
    // PRIMER-TM-001, opt-in. Finds the most stable INTERMOLECULAR antiparallel duplex
    // between two oligonucleotides (self-dimer = an oligo against a second copy of itself;
    // hetero/cross-dimer = two different oligos) and returns its NN ΔH°/ΔS° and the
    // bimolecular Tm. Reuses the SantaLucia & Hicks (2004) Table 1 NN stacking table
    // (ThermoConstants NnParameterSet.SantaLuciaHicks2004), the terminal-A·T penalty, the duplex-initiation term and the
    // 0.368 entropy salt coefficient already used by CalculateMeltingTemperatureNN.
    // The duplex Tm / hairpin Tm / default Tm methods and their defaults are UNCHANGED.
    //
    // Model (gapless contiguous-WC scorer on the ntthal duplex terms — SantaLucia & Hicks 2004 unified NN;
    // the full ntthal DP with mismatches/loops/bulges/overhangs is NtthalDimer via CalculateDimerThermodynamicsNtthal):
    //   For each gapless antiparallel offset of strand2 (read 3'→5') under strand1 (5'→3'),
    //   each maximal contiguous run of Watson-Crick pairs (≥ 1 NN stack) is a candidate
    //   duplex with
    //     ΔH° = ΔH°_init + Σ stacks ΔH° + ΔH°_AT-penalty(per A·T-closed end),
    //     ΔS° = ΔS°_init + Σ stacks ΔS° + ΔS°_AT-penalty(per A·T-closed end)
    //           + 0.368·N_stacks·ln[Na⁺]   (salt correction baked into ΔS°, ntthal saltCorrectS),
    //   and bimolecular
    //     Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15,
    //   with x = 1 when BOTH oligos are reverse-complement palindromes (ntthal symmetry_thermo),
    //   else x = 4. ntthal keeps the candidate with the highest Tm — so does this method.
    //
    // Sources (retrieved & extracted this session, 2026-06-25):
    //   SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415-440 — unified NN parameters
    //     (Table 1, ThermoConstants NnParameterSet.SantaLuciaHicks2004) + the bimolecular Tm Eq. 3 + Eq. 5 entropy
    //     salt correction (0.368 coefficient).
    //   Untergasser A et al. (2012) Nucleic Acids Res 40:e115 (Primer3 2.0) — the ntthal
    //     thermodynamic-alignment engine for oligo dimers.
    //   Primer3 `thal.c` (primer3-py vendored libprimer3, retrieved
    //     https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c):
    //     dplx_init_H=200 cal, dplx_init_S=−5.7 (lines 588-589); AT_H=2200, AT_S=6.9
    //     (lines 128-129); saltCorrectS = 0.368·ln((mv+120·√max(0,dv−dntp))/1000) per stack
    //     (lines 623-624, 1042); RC = R·ln(dna_conc/1e9) when both strands symmetric else
    //     R·ln(dna_conc/4e9) (lines 590-593); symmetry_thermo = reverse-complement palindrome
    //     (line 2771). dna_conc is in nM, so /1e9 → mol/L with x=1 and /4e9 → x=4.
    //   Cross-checked against primer3-py 2.3.0 calc_homodimer / calc_heterodimer
    //     (mv=50, dv=0, dntp=0, dna_conc=50 nM): this method reproduces ntthal's ΔH°, ΔS°
    //     and Tm to machine precision for every case whose optimal structure is a contiguous
    //     Watson-Crick duplex (e.g. GCGCGCGC, ACGTACGTACGT, ATCGATCGATCG/CGATCGATCGAT,
    //     CGATCGATCG self-dimer, GCATGC, GGGGCCCC). By design this record keeps the contiguous-WC
    //     optimum (DNA_Dimer_Tm.md §5.3/§5.4): ntthal's mismatch / loop / bulge / terminal-overhang
    //     (tstack2 / dangle) terms are in the full port NtthalDimer — CalculateDimerThermodynamicsNtthal,
    //     CalculateSelfDimerMeltingTemperature and CalculateDimerMeltingTemperature use it (bit-exact to
    //     primer3-py 2.3.1, incl. dv/dntp, temp_c, max_loop; PRIMER-DIMER-001).

    /// <summary>
    /// The most stable intermolecular DNA duplex (self- or hetero-dimer) found between two
    /// oligonucleotides: the aligned spans on each strand, the number of base pairs, and the
    /// dimer ΔH° (kcal/mol), ΔS° (cal/(K·mol)) and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="Strand1Start">5'-most aligned index (0-based) on strand 1.</param>
    /// <param name="Strand2Start">5'-most aligned index (0-based) on strand 2.</param>
    /// <param name="BasePairs">Number of contiguous Watson-Crick base pairs in the duplex.</param>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the 0.368 salt correction.</param>
    /// <param name="DeltaG37">Dimer ΔG°37 = ΔH° − 310.15·ΔS°/1000 in kcal/mol (negative = stable).</param>
    public readonly record struct DimerResult(
        int Strand1Start, int Strand2Start, int BasePairs,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (highest-Tm) gapless, contiguous Watson–Crick intermolecular DNA duplex between two
    /// oligonucleotides, scored with the <c>ntthal</c> duplex terms over the SantaLucia &amp; Hicks (2004) unified
    /// nearest-neighbour model (the full <c>ntthal</c> alignment with mismatches, loops, bulges and terminal
    /// overhangs is <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double)"/>). A <b>self-dimer</b> is obtained by passing the same
    /// sequence as both strands; a <b>hetero/cross-dimer</b> by passing two different sequences.
    /// <b>Opt-in</b>: the duplex (<see cref="CalculateMeltingTemperatureNN"/>) and hairpin Tm
    /// methods, and the default <see cref="CalculateMeltingTemperature(string)"/>, are unchanged.
    /// <para>
    /// The two strands are aligned antiparallel (strand 2 read 3'→5' under strand 1 5'→3') over
    /// every gapless offset; each maximal contiguous Watson-Crick run of ≥ 2 bp is scored as a
    /// duplex with ΔH° = init + Σ stacks + terminal-A·T penalty per A·T-closed end, and ΔS° the
    /// same plus the 0.368·N<sub>stacks</sub>·ln[Na⁺] salt correction. The candidate with the
    /// highest bimolecular Tm is returned (the ΔS° in the result already includes the salt term
    /// for the supplied <paramref name="sodiumMolar"/>).
    /// </para>
    /// </summary>
    /// <param name="strand1">First DNA oligo (5'→3'). Must contain ≥ 2 ACGT bases.</param>
    /// <param name="strand2">Second DNA oligo (5'→3'); the same string as
    /// <paramref name="strand1"/> for a self-dimer. Must contain ≥ 2 ACGT bases.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM); only
    /// the entropy salt correction depends on it.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L for the
    /// bimolecular Tm used to rank candidates (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The most stable dimer, or <c>null</c> if either strand is null/&lt; 2 bases/contains
    /// a non-ACGT character, or no duplex of ≥ 2 contiguous base pairs exists between them.</returns>
    public static DimerResult? FindMostStableDimer(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        if (string.IsNullOrEmpty(strand1) || string.IsNullOrEmpty(strand2))
            return null;

        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();
        if (s1.Length < 2 || s2.Length < 2)
            return null;
        foreach (char c in s1)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;
        foreach (char c in s2)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        // x = 1 only when both strands are reverse-complement palindromes (ntthal symmetry_thermo).
        bool symmetric = IsSelfComplementary(s1) && IsSelfComplementary(s2);
        double x = symmetric ? SelfComplementaryFactor : NonSelfComplementaryFactor;
        // Strand-concentration term R·ln(C_T / x) for the bimolecular Tm (constant over candidates).
        double rcTerm = GasConstant * Math.Log(strandConcentrationMolar / x);
        double saltPerStack = SantaLuciaEntropySaltCoefficient * Math.Log(sodiumMolar);

        // Strand 2 read 3'→5' so its index i pairs base-for-base under strand 1 read 5'→3'.
        string s2Rev = Reverse(s2);
        int n = s1.Length, m = s2.Length;

        DimerResult? best = null;
        double bestTm = double.NegativeInfinity;

        // Slide strand 2 across strand 1 over every gapless antiparallel offset.
        for (int offset = -(m - 1); offset < n; offset++)
        {
            // Within this offset, walk the overlap and split it into maximal contiguous
            // Watson-Crick runs; each run of ≥ 2 bp is a candidate duplex.
            int runStart = -1; // strand-1 index where the current WC run started
            for (int i = 0; i <= n; i++)
            {
                int j = i - offset; // index into s2Rev paired with s1[i]
                bool paired = i < n && j >= 0 && j < m && IsWatsonCrickPair(s1[i], s2Rev[j]);

                if (paired && runStart < 0)
                    runStart = i;

                if (!paired && runStart >= 0)
                {
                    int runEnd = i - 1; // inclusive strand-1 index of the run end
                    EvaluateRun(s1, s2Rev, offset, runStart, runEnd, saltPerStack, rcTerm,
                                ref best, ref bestTm);
                    runStart = -1;
                }
            }
        }

        return best;

        // Scores one contiguous Watson-Crick run [runStart..runEnd] on strand 1 and keeps it
        // if its bimolecular Tm exceeds the best found so far.
        static void EvaluateRun(
            string s1, string s2Rev, int offset, int runStart, int runEnd,
            double saltPerStack, double rcTerm, ref DimerResult? best, ref double bestTm)
        {
            int basePairs = runEnd - runStart + 1;
            if (basePairs < MinDimerBasePairs)
                return;

            double dH = DesignNnInitiation.Initiation.DeltaH;
            double dS = DesignNnInitiation.Initiation.DeltaS;
            for (int k = runStart; k < runEnd; k++)
            {
                string step = s1.Substring(k, 2);
                // A contiguous WC run pairs every column; the stack is the perfect-match NN
                // keyed by the strand-1 dinucleotide (its complement is the strand-2 stack).
                if (!TryGetDesignStack(step, out var p))
                    return;
                dH += p.DeltaH;
                dS += p.DeltaS;
            }

            // Terminal A·T penalty per duplex end closing with an A·T pair.
            if (s1[runStart] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }
            if (s1[runEnd] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }

            int stacks = basePairs - 1;
            dS += stacks * saltPerStack; // 0.368·N_stacks·ln[Na⁺] (ntthal saltCorrectS)

            double tmKelvin = (dH * 1000.0) / (dS + rcTerm);
            double tmCelsius = tmKelvin - KelvinOffset;

            if (tmCelsius > bestTm)
            {
                bestTm = tmCelsius;
                int strand2Start5 = s2Rev.Length - 1 - (runEnd - offset); // → 5'→3' index on strand 2
                double dG37 = dH - ReferenceTemperatureKelvin * dS / 1000.0;
                best = new DimerResult(runStart, strand2Start5, basePairs, dH, dS, dG37);
            }
        }
    }

    // Minimum base pairs for a dimer duplex (at least one NN stack). ntthal requires a paired
    // region; a single base pair has no stacking energy and is not a duplex.
    private const int MinDimerBasePairs = 2;

    /// <summary>
    /// Computes the intermolecular <b>self-dimer</b> melting temperature (°C) of a DNA oligo: the
    /// bimolecular Tm of the most stable duplex it forms with a second copy of itself, via the
    /// Primer3 / <c>ntthal</c> thermodynamic alignment (SantaLucia &amp; Hicks 2004 unified NN).
    /// <b>Opt-in</b>: the perfect-match duplex / hairpin / default Tm methods are unchanged.
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The self-dimer Tm in °C (ntthal also reports weak / single-pair structures, whose Tm
    /// may be negative), or <c>double.NaN</c> if ntthal finds no structure / the sequence is
    /// null, empty or non-ACGT.</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    public static double CalculateSelfDimerMeltingTemperature(
        string sequence,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar) =>
        CalculateDimerMeltingTemperature(sequence, sequence, sodiumMolar, strandConcentrationMolar);

    /// <summary>
    /// Computes the intermolecular <b>dimer</b> melting temperature (°C) between two DNA oligos
    /// (a self-dimer when both arguments are the same sequence; a hetero/cross-dimer otherwise),
    /// as the bimolecular Tm of the most stable duplex found by the Primer3 / <c>ntthal</c>
    /// thermodynamic alignment (SantaLucia &amp; Hicks 2004 unified NN):
    /// Tm = ΔH°·1000/(ΔS° + R·ln(C_T/x)) − 273.15, x = 1 if both oligos are reverse-complement
    /// palindromes else x = 4 (C_T = the total strand concentration, default 50 nM). <b>Opt-in</b>: existing Tm methods are unchanged.
    /// </summary>
    /// <param name="strand1">First DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="strand2">Second DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The dimer Tm in °C (ntthal also reports weak / single-pair structures, whose Tm may
    /// be negative), or <c>double.NaN</c> if ntthal finds no structure / either sequence is null,
    /// empty or non-ACGT.</returns>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    public static double CalculateDimerMeltingTemperature(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        var dimer = CalculateDimerThermodynamicsNtthal(strand1, strand2, sodiumMolar, strandConcentrationMolar);
        return dimer is null ? double.NaN : dimer.Value.TmCelsius;
    }

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics (ΔH°, ΔS°, ΔG°37 and bimolecular Tm) for the most
    /// stable intermolecular DNA duplex between two oligos, computed by the complete Primer3
    /// <c>ntthal</c> dynamic program (mode ANY): matched nearest-neighbour stacks, single internal
    /// mismatches, internal loops, single- and multi-base bulges, and terminal overhangs /
    /// dangling ends (the <c>tstack2</c> terminal table + 5′/3′ dangling-end tables + interior /
    /// bulge loop-length parameters). Unlike <see cref="FindMostStableDimer"/> (which scores only
    /// the best contiguous Watson–Crick run), this reproduces primer3-py's
    /// <c>calc_homodimer</c>/<c>calc_heterodimer</c> for dimers whose optimum is <b>non-contiguous</b>.
    /// This overload uses monovalent salt only (primer3-py with <c>dv_conc=0, dntp_conc=0</c>; note the
    /// primer3-py defaults are 1.5 mM Mg²⁺ / 0.6 mM dNTP — use the overload taking
    /// <see cref="NtthalAlignmentMode"/> for those).
    /// <b>Opt-in</b>: all other Tm methods and defaults are unchanged.
    /// </summary>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    /// <param name="strand1">First DNA oligo (5′→3′); ≥ 1 ACGT base.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′); the same string for a self-dimer.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The most stable dimer's thermodynamics, or <c>null</c> if either strand is
    /// null/empty/contains a non-ACGT character, or no duplex can be formed (ntthal
    /// <c>no_structure</c>).</returns>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        if (string.IsNullOrEmpty(strand1) || string.IsNullOrEmpty(strand2))
            return null;
        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();
        foreach (char c in s1)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;
        foreach (char c in s2)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        var r = NtthalDimer.Run(s1, s2, sodiumMolar, strandConcentrationMolar);
        if (r is null)
            return null;

        var v = r.Value;
        // Convert ntthal native cal/mol → the library's kcal/mol convention for ΔH/ΔG.
        return new DimerThermodynamics(
            DeltaH: v.DeltaH / 1000.0,
            DeltaS: v.DeltaS,
            DeltaG37: v.DeltaG37 / 1000.0,
            TmCelsius: v.TmCelsius,
            BasePairs: v.BasePairs);
    }

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics of the most stable intermolecular duplex.
    /// </summary>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the N·saltCorrection term.</param>
    /// <param name="DeltaG37">Dimer ΔG° = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at the
    /// analysis temperature T (310.15 K = 37 °C unless an overload with <c>temperatureCelsius</c> is used;
    /// ntthal / primer3-py <c>temp_c</c>). ntthal may report a positive ΔG for its optimal structure.</param>
    /// <param name="TmCelsius">Bimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">Number of paired bases in the optimal structure.</param>
    public readonly record struct DimerThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

    /// <summary>
    /// Computes the full <b>ntthal</b> intramolecular-hairpin thermodynamics (ΔH°, ΔS°, ΔG°37, Tm)
    /// of a DNA oligo, reproducing primer3-py's <c>calc_hairpin</c>. This runs the complete
    /// Primer3 <c>ntthal</c> monomer dynamic program (a single stem with internal
    /// mismatches/loops, bulges, terminal mismatch / dangling-end terminal contributions and the
    /// size-keyed hairpin-loop initiation) and — unlike <see cref="FindMostStableHairpin"/>'s
    /// SantaLucia &amp; Hicks (2004) Table 4 model — <b>automatically applies the bundled
    /// sequence-specific special triloop / tetraloop stability bonuses</b> (the primer3
    /// <c>triloop.dh/.ds</c> + <c>tetraloop.dh/.ds</c> tables, keyed on the full loop string
    /// including the closing base pair). No caller-supplied loop bonus is required.
    /// <b>Opt-in</b>: <see cref="FindMostStableHairpin"/>, the duplex/dimer Tm methods and all
    /// defaults are unchanged.
    /// </summary>
    /// <param name="sequence">The DNA oligo (5′→3'); ≥ 1 ACGT base.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM, the
    /// primer3 <c>calc_hairpin</c> default of <c>mv=50</c>). This overload has no divalent cations or
    /// dNTPs: it equals <c>calc_hairpin(seq, mv_conc, dv_conc=0, dntp_conc=0)</c> at 37 °C, max loop 30.</param>
    /// <returns>The most stable hairpin's thermodynamics, or <c>null</c> if the sequence is
    /// null/empty/contains a non-ACGT character, or no hairpin can form (ntthal
    /// <c>no_structure</c>, e.g. a homopolymer).</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c <c>THAL_MAX_ALIGN</c>).</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;
        string seq = sequence.ToUpperInvariant();
        foreach (char c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        var r = NtthalHairpin.Run(seq, sodiumMolar);
        if (r is null)
            return null;

        var v = r.Value;
        // Convert ntthal native cal/mol → the library's kcal/mol convention for ΔH/ΔG.
        return new HairpinThermodynamics(
            DeltaH: v.DeltaH / 1000.0,
            DeltaS: v.DeltaS,
            DeltaG37: v.DeltaG37 / 1000.0,
            TmCelsius: v.TmCelsius,
            BasePairs: v.BasePairs);
    }

    /// <summary>
    /// Full <c>ntthal</c> intramolecular-hairpin thermodynamics of the most stable hairpin
    /// (reproduces primer3-py <c>calc_hairpin</c>; special tri/tetraloop bonuses applied).
    /// </summary>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol), including the (N/2−1)·saltCorrection term.</param>
    /// <param name="DeltaG37">Hairpin ΔG = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at T = 310.15 K
    /// (37 °C), or at <c>temperatureCelsius</c> + 273.15 for the overloads taking it (primer3-py <c>temp_c</c>).</param>
    /// <param name="TmCelsius">Unimolecular melting temperature in °C (no strand-concentration term).</param>
    /// <param name="BasePairs">ntthal N/2: half the number of paired positions among bases 1..len−1 of the
    /// optimal structure (the count thal.c uses in the (N/2 − 1)·saltCorrection term).</param>
    public readonly record struct HairpinThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

    /// <summary>
    /// ntthal dimer alignment type (Primer3 <c>thal_alignment_type</c>, ntthal <c>-a</c>).
    /// </summary>
    public enum NtthalAlignmentMode
    {
        /// <summary>THAL_ANY: the most stable duplex anywhere (primer3-py <c>calc_heterodimer</c>).</summary>
        Any,
        /// <summary>THAL_END1: the duplex must contain the 3′-terminal base of strand 1
        /// (primer3-py <c>calc_end_stability(strand1, strand2)</c>).</summary>
        End1,
        /// <summary>THAL_END2: the duplex must contain the 3′-terminal base of strand 2
        /// (= END1 with the strands swapped).</summary>
        End2,
    }

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics with an explicit alignment type and the complete
    /// ntthal salt model (<c>saltCorrectS</c>: 0.368·ln((mv + 120·√max(0, dv − dntp))/1000), mM),
    /// reproducing primer3-py <c>calc_heterodimer</c> (mode <see cref="NtthalAlignmentMode.Any"/>)
    /// and <c>calc_end_stability</c> (mode <see cref="NtthalAlignmentMode.End1"/>) at any
    /// mv/dv/dntp/dna_conc (temperature 37 °C, max loop 30 — the primer3-py defaults).
    /// </summary>
    /// <param name="strand1">First DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="mode">ntthal alignment type.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="strandConcentrationMolar">Oligo concentration, mol/L (ntthal dna_conc).</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no duplex forms.</returns>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt, or either is longer
    /// than 10 000 nt (thal.c <c>THAL_MAX_ALIGN</c> / <c>THAL_MAX_SEQ</c>; primer3-py raises).</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar) =>
        CalculateDimerThermodynamicsNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, NtthalDefaultTemperatureCelsius, NtthalDefaultMaxLoop);

    /// <summary>primer3-py / ntthal default analysis temperature (°C) at which ΔG is reported.</summary>
    public const double NtthalDefaultTemperatureCelsius = 37.0;

    /// <summary>primer3-py / ntthal default (and maximum) internal-loop / bulge size.</summary>
    public const int NtthalDefaultMaxLoop = 30;

    /// <summary>
    /// thal.h <c>THAL_MAX_ALIGN</c> = 60: Primer3 / primer3-py refuse a dimer whose two strands are both longer, and a
    /// hairpin of a longer oligo ("At least one sequence must be equal to or shorter than 60bp …"). It is a
    /// compile-time constant of thal.c (<c>#ifndef THAL_MAX_ALIGN</c>, used only by <c>thal_check_errors</c>; the DP
    /// tables are allocated from the actual lengths) chosen by the Primer3 authors as "the maximum reasonable length
    /// for nearest neighbor models … only two states of melting" (thal.h). The ntthal overloads taking
    /// <c>maxAlignLength</c> raise it (opt-in) and then reproduce ntthal built with <c>-DTHAL_MAX_ALIGN=…</c>.
    /// </summary>
    public const int NtthalMaxAlignLength = 60;

    /// <summary>thal.h <c>THAL_MAX_SEQ</c> = 10 000: maximum length of either strand (and the largest
    /// <c>maxAlignLength</c> accepted).</summary>
    public const int NtthalMaxSequenceLength = 10000;

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics with every primer3-py <c>calc_heterodimer</c> /
    /// <c>calc_end_stability</c> argument: alignment type, mv/dv/dntp/dna_conc, <c>temp_c</c> (the
    /// temperature at which ΔG is evaluated: ΔG = ΔH − (temp_c + 273.15)·ΔS, thal.c <c>calcDimer</c>;
    /// the DP ranking and Tm do not depend on it) and <c>max_loop</c> (largest internal loop / bulge
    /// considered, 0–30). Bit-faithful port of primer3-py 2.3.1 <c>thal.c</c>.
    /// </summary>
    /// <param name="strand1">First DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="mode">ntthal alignment type.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="strandConcentrationMolar">Oligo concentration, mol/L (ntthal dna_conc).</param>
    /// <param name="temperatureCelsius">primer3-py <c>temp_c</c>; the returned
    /// <see cref="DimerThermodynamics.DeltaG37"/> is ΔG at this temperature.</param>
    /// <param name="maxLoop">primer3-py <c>max_loop</c>, 0–30.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no duplex forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt, or either is longer
    /// than 10 000 nt (thal.c <c>THAL_MAX_ALIGN</c> / <c>THAL_MAX_SEQ</c>).</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: false,
            NtthalMaxAlignLength)?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>
    /// with thal.h <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see
    /// <see cref="NtthalMaxAlignLength"/>): the unchanged ntthal recursions run when at least one strand is at most
    /// <paramref name="maxAlignLength"/> nt — identical to ntthal / thal.c compiled with
    /// <c>-DTHAL_MAX_ALIGN=maxAlignLength</c>. With 60 it is the 9-argument overload. Cost is O(len1·len2·maxLoop²).
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, <see cref="NtthalMaxAlignLength"/> (60) …
    /// <see cref="NtthalMaxSequenceLength"/> (10 000) — smaller values would refuse what Primer3 accepts.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30 or
    /// <paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">Both strands are longer than <paramref name="maxAlignLength"/>, or either
    /// is longer than 10 000 nt.</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: false,
            CheckedMaxAlignLength(maxAlignLength))?.Thermodynamics;

    // Public THAL_MAX_ALIGN override range: 60 (Primer3) … 10 000 (THAL_MAX_SEQ).
    private static int CheckedMaxAlignLength(int maxAlignLength)
    {
        if (maxAlignLength < NtthalMaxAlignLength || maxAlignLength > NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(nameof(maxAlignLength), maxAlignLength,
                "THAL_MAX_ALIGN override must be in 60..10000 (Primer3 default .. THAL_MAX_SEQ).");
        return maxAlignLength;
    }

    /// <summary>
    /// As <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>,
    /// plus the optimal duplex drawn exactly as thal.c <c>drawDimer</c> / primer3-py
    /// <c>ThermoResult.ascii_structure_lines</c> (<c>output_structure=True</c>): four lines
    /// "SEQ	…" (unpaired strand-1 bases), "SEQ	…" (paired strand-1 bases), "STR	…" (paired
    /// strand-2 bases), "STR	…" (unpaired strand-2 bases); strand 2 runs 3′→5′ and '-' pads the
    /// shorter side of a loop. In mode <see cref="NtthalAlignmentMode.End2"/> the strands are drawn
    /// swapped, as ntthal does.
    /// </summary>
    /// <returns>The thermodynamics and structure lines, or <c>null</c> for invalid input or when no
    /// duplex forms.</returns>
    public static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode = NtthalAlignmentMode.Any,
        double sodiumMolar = 0.05,
        double divalentMolar = 0.0015,
        double dntpMolar = 0.0006,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar,
        double temperatureCelsius = NtthalDefaultTemperatureCelsius,
        int maxLoop = NtthalDefaultMaxLoop) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: true, NtthalMaxAlignLength);

    /// <summary>
    /// As <see cref="CalculateDimerStructureNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>
    /// with THAL_MAX_ALIGN raised to <paramref name="maxAlignLength"/> (opt-in, 60–10 000; see
    /// <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int, int)"/>).
    /// </summary>
    public static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: true,
            CheckedMaxAlignLength(maxAlignLength));

    private static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1, string strand2, NtthalAlignmentMode mode, double sodiumMolar, double divalentMolar,
        double dntpMolar, double strandConcentrationMolar, double temperatureCelsius, int maxLoop, bool withStructure,
        int maxAlign)
    {
        if (!IsAcgtOnly(strand1) || !IsAcgtOnly(strand2))
            return null;
        var type = mode switch
        {
            NtthalAlignmentMode.Any => NtthalDimer.AlignmentType.Any,
            NtthalAlignmentMode.End1 => NtthalDimer.AlignmentType.End1,
            NtthalAlignmentMode.End2 => NtthalDimer.AlignmentType.End2,
            _ => throw new ArgumentOutOfRangeException(nameof(mode)),
        };
        var r = NtthalDimer.Run(strand1.ToUpperInvariant(), strand2.ToUpperInvariant(),
            sodiumMolar, strandConcentrationMolar, type, divalentMolar, dntpMolar,
            temperatureCelsius + KelvinOffset, maxLoop, withStructure, maxAlign);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalDimerStructure(
            new DimerThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
    }

    /// <summary>
    /// ntthal dimer thermodynamics plus the thal.c <c>drawDimer</c> ASCII duplex
    /// (primer3-py <c>ThermoResult.ascii_structure_lines</c>).
    /// </summary>
    /// <param name="Thermodynamics">ΔH/ΔS/ΔG/Tm of the optimal duplex.</param>
    /// <param name="AsciiStructureLines">The four "SEQ	"/"SEQ	"/"STR	"/"STR	" lines.</param>
    public sealed record NtthalDimerStructure(DimerThermodynamics Thermodynamics, IReadOnlyList<string> AsciiStructureLines);

    /// <summary>
    /// Full <c>ntthal</c> hairpin thermodynamics with the complete ntthal salt model (divalent
    /// cations and dNTPs enter through <c>saltCorrectS</c>), reproducing primer3-py
    /// <c>calc_hairpin</c> at any mv/dv/dntp (temperature 37 °C, max loop 30 — the primer3-py defaults).
    /// </summary>
    /// <param name="sequence">DNA oligo (5′→3′), ACGT only, at most 60 nt.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no hairpin forms.</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c
    /// <c>THAL_MAX_ALIGN</c>; primer3-py raises).</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar) =>
        CalculateHairpinThermodynamicsNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            NtthalDefaultTemperatureCelsius, NtthalDefaultMaxLoop);

    /// <summary>
    /// Full <c>ntthal</c> hairpin thermodynamics with every primer3-py <c>calc_hairpin</c> argument:
    /// mv/dv/dntp, <c>temp_c</c> (the analysis temperature: it sets the reported ΔG = ΔH − (temp_c +
    /// 273.15)·ΔS and, as in thal.c <c>calc_terminal_bp</c>, the exterior-loop acceptance test
    /// ΔH − T·ΔS &lt; 0, so it can change the selected structure) and <c>max_loop</c> (largest internal
    /// loop / bulge considered, 0–30). Bit-faithful port of primer3-py 2.3.1 <c>thal.c</c> (type 4).
    /// </summary>
    /// <param name="sequence">DNA oligo (5′→3′), ACGT only (case-insensitive), at most 60 nt.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="temperatureCelsius">primer3-py <c>temp_c</c>; the returned
    /// <see cref="HairpinThermodynamics.DeltaG37"/> is ΔG at this temperature.</param>
    /// <param name="maxLoop">primer3-py <c>max_loop</c>, 0–30.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no hairpin forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt.</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: false, NtthalMaxAlignLength)?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int)"/> with thal.h
    /// <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see <see cref="NtthalMaxAlignLength"/>):
    /// oligos up to <paramref name="maxAlignLength"/> nt are folded by the unchanged ntthal hairpin recursions —
    /// identical to ntthal / thal.c compiled with <c>-DTHAL_MAX_ALIGN=maxAlignLength</c> (verified on 61–120-mers
    /// against thal.c from primer3-py 2.3.1 built with a larger THAL_MAX_ALIGN). With 60 it is the 6-argument
    /// overload. Cost is O(n²·maxLoop²). Note that the Primer3 authors chose 60 as the limit of the two-state
    /// nearest-neighbour model (thal.h); longer oligos are folded with the same single-structure model.
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, <see cref="NtthalMaxAlignLength"/> (60) …
    /// <see cref="NtthalMaxSequenceLength"/> (10 000).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30 or
    /// <paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than <paramref name="maxAlignLength"/>.</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: false, CheckedMaxAlignLength(maxAlignLength))?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int)"/>,
    /// plus the optimal hairpin drawn exactly as thal.c <c>drawHairpin</c> / primer3-py
    /// <c>ThermoResult.ascii_structure_lines</c> (<c>output_structure=True</c>): two lines,
    /// "SEQ\t" followed by one character per base ('-' unpaired; for each base pair the 5′ partner
    /// is drawn '/' and the 3′ partner '\') and "STR\t" followed by the (upper-case) oligo.
    /// Defaults are the primer3-py <c>calc_hairpin</c> defaults (mv 50 mM, dv 1.5 mM, dNTP 0.6 mM,
    /// 37 °C, max loop 30).
    /// </summary>
    /// <returns>The thermodynamics and structure lines, or <c>null</c> for invalid input or when no
    /// hairpin forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt.</exception>
    public static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence,
        double sodiumMolar = 0.05,
        double divalentMolar = 0.0015,
        double dntpMolar = 0.0006,
        double temperatureCelsius = NtthalDefaultTemperatureCelsius,
        int maxLoop = NtthalDefaultMaxLoop) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: true, NtthalMaxAlignLength);

    /// <summary>
    /// As <see cref="CalculateHairpinStructureNtthal(string, double, double, double, double, int)"/> with THAL_MAX_ALIGN
    /// raised to <paramref name="maxAlignLength"/> (opt-in, 60–10 000; see
    /// <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int, int)"/>).
    /// </summary>
    public static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: true, CheckedMaxAlignLength(maxAlignLength));

    private static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence, double sodiumMolar, double divalentMolar, double dntpMolar,
        double temperatureCelsius, int maxLoop, bool withStructure, int maxAlign)
    {
        if (!IsAcgtOnly(sequence))
            return null;
        var r = NtthalHairpin.Run(sequence.ToUpperInvariant(), sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius + KelvinOffset, maxLoop, withStructure, maxAlign);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalHairpinStructure(
            new HairpinThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
    }

    /// <summary>
    /// ntthal hairpin thermodynamics plus the thal.c <c>drawHairpin</c> ASCII structure
    /// (primer3-py <c>ThermoResult.ascii_structure_lines</c>).
    /// </summary>
    /// <param name="Thermodynamics">ΔH/ΔS/ΔG/Tm of the optimal hairpin.</param>
    /// <param name="AsciiStructureLines">The "SEQ\t…" and "STR\t…" lines.</param>
    public sealed record NtthalHairpinStructure(HairpinThermodynamics Thermodynamics, IReadOnlyList<string> AsciiStructureLines);

    private static bool IsAcgtOnly(string? s)
    {
        if (string.IsNullOrEmpty(s)) return false;
        foreach (char c in s)
            if (c is not ('A' or 'C' or 'G' or 'T' or 'a' or 'c' or 'g' or 't')) return false;
        return true;
    }

    /// <summary>
    /// Primer3's default (thermodynamic) limit, in °C, on every secondary-structure Tm:
    /// <c>PRIMER_MAX_SELF_ANY_TH</c>, <c>PRIMER_MAX_SELF_END_TH</c>, <c>PRIMER_MAX_HAIRPIN_TH</c>,
    /// <c>PRIMER_PAIR_MAX_COMPL_ANY_TH</c>, <c>PRIMER_PAIR_MAX_COMPL_END_TH</c> = 47.0
    /// (<c>libprimer3.cc</c> <c>pr_set_default_global_args_2</c>). A value is a violation when it is
    /// strictly greater than the limit.
    /// </summary>
    public const double Primer3MaxStructureTm = 47.0;

    /// <summary>Primer3 primer default PRIMER_SALT_MONOVALENT = 50 mM (<c>libprimer3.c</c> <c>p_args.salt_conc</c>).</summary>
    public const double Primer3MonovalentMillimolar = 50.0;

    /// <summary>Primer3 primer default PRIMER_SALT_DIVALENT = 1.5 mM (<c>pr_set_default_global_args_2</c> <c>p_args.divalent_conc</c>).</summary>
    public const double Primer3DivalentMillimolar = 1.5;

    /// <summary>Primer3 primer default PRIMER_DNTP_CONC = 0.6 mM (<c>pr_set_default_global_args_2</c> <c>p_args.dntp_conc</c>).</summary>
    public const double Primer3DntpMillimolar = 0.6;

    /// <summary>Primer3 primer default PRIMER_DNA_CONC = 50 nM (<c>p_args.dna_conc</c>).</summary>
    public const double Primer3DnaConcentrationNanomolar = 50.0;

    // Primer3 _pr_data_control messages for a GC weight without PRIMER_[INTERNAL_]OPT_GC_PERCENT.
    internal const string PrimerGcOptimumUndefinedMessage =
        "Primer GC content is part of objective function while optimum gc_content is not defined (Primer3 _pr_data_control).";
    internal const string ProbeGcOptimumUndefinedMessage =
        "Hyb probe GC content is part of objective function while optimum gc_content is not defined (Primer3 _pr_data_control).";

    /// <summary>
    /// The Primer3 manual's documented PRIMER_OPT_GC_PERCENT (50 %). It is <b>not</b> applied by default: Primer3's code keeps
    /// the optimum undefined (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>, <c>pr_set_default_global_args_1</c>)
    /// and its <c>_pr_data_control</c> rejects a non-zero PRIMER_(INTERNAL_)WT_GC_PERCENT_GT/_LT without an explicit optimum,
    /// as <see cref="PrimerParameters.OptimalGcPercent"/> / <see cref="ProbeDesigner.Primer3ProbeSettings.OptGcPercent"/> /
    /// <see cref="Primer3Optima.OptGcPercent"/> = null now do.
    /// </summary>
    [Obsolete("Primer3 leaves PRIMER_OPT_GC_PERCENT undefined by default; set PrimerParameters.OptimalGcPercent explicitly.")]
    public const double Primer3DefaultOptGcPercent = 50.0;

    /// <summary>
    /// Primer3 default PRIMER_MAX_END_STABILITY = 100 kcal/mol (<c>pr_set_default_global_args_1</c>
    /// <c>max_end_stability</c>): a left/right primer fails when its <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="Calculate3PrimeStability"/> is strictly greater. The largest
    /// value an ACGT 3′ pentamer reaches is 6.86 (GCGCG), so the default never rejects a primer.
    /// </summary>
    public const double Primer3MaxEndStability = 100.0;

    /// <summary>Primer3 default PRIMER_GC_CLAMP = 0 (no G/C required at the 3′ end; <c>gc_clamp</c>).</summary>
    public const int Primer3GcClamp = 0;

    /// <summary>
    /// Primer3 default PRIMER_MAX_END_GC = 5 (<c>max_end_gc</c>): the maximum number of G/C among the five 3′-most
    /// bases of a left/right primer; Primer3 only checks it when it is &lt; 5. Legal range 0–5.
    /// </summary>
    public const int Primer3MaxEndGc = 5;

    // Primer3 data control of one oligo-condition set (_pr_data_control: "Illegal value for primer salt or dna
    // concentration" / "… divalent salt or dNTP concentration"): salt and DNA concentration > 0, divalent ≥ 0;
    // a negative dNTP is also rejected here (seqtm has no meaning for it). NaN / ∞ are rejected.
    internal static void ValidatePrimer3Conditions(
        double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar, double dnaConcentrationNanomolar,
        string paramName)
    {
        if (!(monovalentMillimolar > 0) || double.IsInfinity(monovalentMillimolar)
            || !(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for salt or DNA concentration (Primer3: PRIMER_SALT_MONOVALENT and PRIMER_DNA_CONC must be > 0).");
        if (!(divalentMillimolar >= 0) || double.IsInfinity(divalentMillimolar)
            || !(dntpMillimolar >= 0) || double.IsInfinity(dntpMillimolar))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for divalent salt or dNTP concentration (Primer3: PRIMER_SALT_DIVALENT and PRIMER_DNTP_CONC must be ≥ 0).");
    }

    // oligotm.c divalent_to_monovalent added to the monovalent salt: [Mon] + 120·√([Mg²⁺] − [dNTP]) (mM); no divalent ⇒
    // dNTP ignored, Mg ≤ dNTP ⇒ no contribution.
    internal static double Primer3MonovalentEquivalent(double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar)
    {
        double freeDivalent = divalentMillimolar == 0 ? 0 : Math.Max(0, divalentMillimolar - dntpMillimolar);
        return monovalentMillimolar + Primer3DivalentFactor * Math.Sqrt(freeDivalent);
    }

    /// <summary>Primer3 hybridization-probe (internal-oligo) default PRIMER_INTERNAL_DNA_CONC = 50 nM (<c>libprimer3.cc</c> <c>o_args.dna_conc</c>).</summary>
    public const double Primer3InternalDnaConcentrationNanomolar = 50.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_SALT_MONOVALENT = 50 mM (<c>o_args.salt_conc</c>).</summary>
    public const double Primer3InternalMonovalentMillimolar = 50.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_SALT_DIVALENT = 0 mM (<c>o_args.divalent_conc</c>).</summary>
    public const double Primer3InternalDivalentMillimolar = 0.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_DNTP_CONC = 0 mM (<c>o_args.dntp_conc</c>).</summary>
    public const double Primer3InternalDntpMillimolar = 0.0;

    /// <summary>
    /// Primer3 thermodynamic secondary-structure values of one primer (Tm in °C; 0 when ntthal finds
    /// no structure or the Tm is below 0 °C, as <c>align_thermod</c> reports): <c>PRIMER_*_SELF_ANY_TH</c>, <c>PRIMER_*_SELF_END_TH</c>, <c>PRIMER_*_HAIRPIN_TH</c>.
    /// </summary>
    /// <param name="SelfAnyTh">Self-dimer Tm, ntthal ANY of (primer, primer).</param>
    /// <param name="SelfEndTh">3′-anchored self-dimer Tm, ntthal END1 of (primer, primer).</param>
    /// <param name="HairpinTh">Hairpin Tm, ntthal HAIRPIN of the primer.</param>
    public readonly record struct Primer3OligoStructure(double SelfAnyTh, double SelfEndTh, double HairpinTh)
    {
        /// <summary>True when any value exceeds <paramref name="maxTm"/> (Primer3 rejects the primer).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) =>
            SelfAnyTh > maxTm || SelfEndTh > maxTm || HairpinTh > maxTm;
    }

    /// <summary>
    /// Primer3 thermodynamic pair complementarity (Tm in °C; 0 when no structure or Tm &lt; 0 °C):
    /// <c>PRIMER_PAIR_COMPL_ANY_TH</c> and <c>PRIMER_PAIR_COMPL_END_TH</c>.
    /// </summary>
    /// <param name="ComplAnyTh">Hetero-dimer Tm, ntthal ANY of (left, right).</param>
    /// <param name="ComplEndTh">Max of ntthal END1/END2 of (left, right) and of (rc(right), rc(left)).</param>
    public readonly record struct Primer3PairComplementarity(double ComplAnyTh, double ComplEndTh)
    {
        /// <summary>True when either value exceeds <paramref name="maxTm"/> (Primer3 rejects the pair).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) => ComplAnyTh > maxTm || ComplEndTh > maxTm;
    }

    /// <summary>
    /// Computes a primer's Primer3 thermodynamic secondary-structure Tm values exactly as Primer3
    /// (default <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1</c>) does in <c>oligo_compl_thermod</c> /
    /// <c>oligo_hairpin</c>: self_any = ntthal ANY(primer, primer), self_end = ntthal END1(primer,
    /// primer), hairpin = ntthal HAIRPIN(primer), each the Tm (°C) of the most stable structure, 0
    /// when none forms (ntthal <c>no_structure</c>) or the Tm is negative (<c>align_thermod</c>). Conditions default to Primer3's primer
    /// conditions (50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo). Cross-checked against
    /// primer3-py 2.3.1 <c>calc_homodimer</c>/<c>calc_end_stability</c>/<c>calc_hairpin</c> and the
    /// <c>PRIMER_LEFT_0_SELF_ANY_TH</c>/<c>_SELF_END_TH</c>/<c>_HAIRPIN_TH</c> values of
    /// <c>design_primers</c>.
    /// </summary>
    /// <param name="primer">Primer (5′→3′), case-insensitive, ACGT only.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM.</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM.</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM.</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM.</param>
    /// <returns>The three Tm values, or <c>null</c> when the primer is null/empty or contains a
    /// non-ACGT character.</returns>
    public static Primer3OligoStructure? CalculatePrimer3OligoStructure(
        string primer,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar) =>
        CalculatePrimer3OligoStructureCore(primer, monovalentMillimolar, divalentMillimolar, dntpMillimolar,
            dnaConcentrationNanomolar, NtthalMaxAlignLength);

    /// <summary>
    /// As <see cref="CalculatePrimer3OligoStructure(string, double, double, double, double)"/> with thal.h
    /// <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see <see cref="NtthalMaxAlignLength"/>):
    /// self_any / self_end / hairpin of oligos up to <paramref name="maxAlignLength"/> nt by the unchanged ntthal
    /// recursions (= thal.c compiled with <c>-DTHAL_MAX_ALIGN=maxAlignLength</c>). Primer3 itself never sees such
    /// oligos (PRIMER_MAX_SIZE ≤ 35 / internal ≤ 36 nt and THAL_MAX_ALIGN = 60).
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, 60 … 10 000.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">The oligo is longer than <paramref name="maxAlignLength"/>.</exception>
    public static Primer3OligoStructure? CalculatePrimer3OligoStructure(
        string primer,
        double monovalentMillimolar,
        double divalentMillimolar,
        double dntpMillimolar,
        double dnaConcentrationNanomolar,
        int maxAlignLength) =>
        CalculatePrimer3OligoStructureCore(primer, monovalentMillimolar, divalentMillimolar, dntpMillimolar,
            dnaConcentrationNanomolar, CheckedMaxAlignLength(maxAlignLength));

    private static Primer3OligoStructure? CalculatePrimer3OligoStructureCore(
        string primer, double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar,
        double dnaConcentrationNanomolar, int maxAlign)
    {
        if (!IsAcgtOnly(primer))
            return null;
        string p = primer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        double any = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.Any, dv, dntp, maxAlign: maxAlign));
        double end = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.End1, dv, dntp, maxAlign: maxAlign));
        var h = NtthalHairpin.Run(p, mv, dv, dntp, maxAlign: maxAlign);
        return new Primer3OligoStructure(any, end, h is null ? 0.0 : Math.Max(0.0, h.Value.TmCelsius));
    }

    /// <summary>
    /// Computes Primer3's thermodynamic pair complementarity between a left (forward) and a right
    /// (reverse) primer, both given 5′→3′ as synthesised, exactly as <c>characterize_pair</c> does in
    /// the default thermodynamic mode: compl_any = ntthal ANY(left, right); compl_end = the maximum of
    /// ntthal END1(left, right), END2(left, right), END1(rc(right), rc(left)) and END2(rc(right), rc(left))
    /// (Primer3 evaluates the last two on the reverse complements, <c>align_thermod(s2, s1_rev, …)</c>).
    /// Tm in °C, 0 when no structure or Tm &lt; 0 °C. Conditions default to Primer3's primer conditions.
    /// </summary>
    /// <param name="leftPrimer">Forward primer (5′→3′), ACGT only.</param>
    /// <param name="rightPrimer">Reverse primer (5′→3′), ACGT only.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM.</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM.</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM.</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM.</param>
    /// <returns>The pair values, or <c>null</c> when either primer is null/empty or contains a
    /// non-ACGT character.</returns>
    public static Primer3PairComplementarity? CalculatePrimer3PairComplementarity(
        string leftPrimer,
        string rightPrimer,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar)
    {
        if (!IsAcgtOnly(leftPrimer) || !IsAcgtOnly(rightPrimer))
            return null;
        string l = leftPrimer.ToUpperInvariant(), r = rightPrimer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        var (any, end) = Primer3PairTms(l, r, mv, dv, dntp, conc, double.PositiveInfinity);
        return new Primer3PairComplementarity(any, end);
    }

    // characterize_pair (thermodynamic mode): s1 = left, s2_rev = right, s2 = revcomp(right),
    // s1_rev = revcomp(left); compl_any = ANY(s1, s2_rev), compl_end = max(END1/END2(s1, s2_rev),
    // END1/END2(s2, s1_rev)). Stops as soon as a value exceeds stopAbove (the remaining values are
    // then irrelevant to a pass/fail decision).
    private static (double ComplAny, double ComplEnd) Primer3PairTms(
        string l, string r, double mv, double dv, double dntp, double conc, double stopAbove)
    {
        double any = TmOrZero(NtthalDimer.Run(l, r, mv, conc, NtthalDimer.AlignmentType.Any, dv, dntp));
        if (any > stopAbove)
            return (any, 0.0);
        string rcL = DnaSequence.GetReverseComplementString(l), rcR = DnaSequence.GetReverseComplementString(r);
        double end = 0.0;
        foreach (var (a, b, t) in new[]
                 {
                     (l, r, NtthalDimer.AlignmentType.End1), (l, r, NtthalDimer.AlignmentType.End2),
                     (rcR, rcL, NtthalDimer.AlignmentType.End1), (rcR, rcL, NtthalDimer.AlignmentType.End2),
                 })
        {
            end = Math.Max(end, TmOrZero(NtthalDimer.Run(a, b, mv, conc, t, dv, dntp)));
            if (end > stopAbove)
                break;
        }
        return (any, end);
    }

    // libprimer3.cc align_thermod: Tm of the structure, 0 when none forms or when Tm < 0 °C.
    private static double TmOrZero(NtthalDimer.Result? r) => r is null ? 0.0 : Math.Max(0.0, r.Value.TmCelsius);

    // Watson-Crick complement (same left-to-right order) via the canonical Core per-base complement
    // (SequenceExtensions.TryGetComplement); identical to the former local ACGT switch on the
    // validated ACGT duplexes it is applied to.
    private static string Complement(string seq) =>
        string.Create(seq.Length, seq, static (dest, src) => src.AsSpan().TryGetComplement(dest));

    /// <summary>True if the sequence equals its own reverse complement (self-complementary).</summary>
    private static bool IsSelfComplementary(string seq)
    {
        int n = seq.Length;
        if (n % 2 != 0) return false; // odd-length cannot be self-complementary
        for (int i = 0; i < n; i++)
        {
            char a = seq[i];
            char b = seq[n - 1 - i];
            bool pair = (a == 'A' && b == 'T') || (a == 'T' && b == 'A')
                     || (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
            if (!pair) return false;
        }
        return true;
    }

    /// <summary>
    /// Generates all possible primers for a region.
    /// </summary>
    public static IEnumerable<PrimerCandidate> GeneratePrimerCandidates(
        DnaSequence template,
        int regionStart,
        int regionEnd,
        bool forward = true,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;

        for (int start = regionStart; start + param.MinLength <= regionEnd; start++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= regionEnd; len++)
            {
                var seq = template.Sequence.Substring(start, len);
                if (!forward)
                    seq = DnaSequence.GetReverseComplementString(seq);

                yield return EvaluatePrimer(seq, start, forward, param);
            }
        }
    }

    // ---- Primer3 weighted penalty objective (PRIMER-TM-001) -------------------
    // Reproduces the per-primer objective function `p_obj_fn` (left/right primer
    // branch) from Primer3's reference source `libprimer3.cc`, with the documented
    // default weights and optima. Lower penalty = better primer, exactly as Primer3.
    // Source: Primer3 source libprimer3.cc p_obj_fn / pr_set_default_global_args_2;
    //         Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE";
    //         Untergasser et al. (2012) NAR 40(15):e115; Koressaar & Remm (2007).

    /// <summary>
    /// Primer3 default per-primer objective weights, taken verbatim from
    /// <c>pr_set_default_global_args_2</c> in Primer3's <c>libprimer3.cc</c>:
    /// length/Tm weights = 1; GC, self-complementarity and N weights = 0.
    /// Source: Primer3 source (branch main), function pr_set_default_global_args_2.
    /// </summary>
    public static readonly Primer3PenaltyWeights DefaultPrimer3Weights = new(
        TmGt: 1.0,           // PRIMER_WT_TM_GT  (libprimer3.cc: weights.temp_gt = 1)
        TmLt: 1.0,           // PRIMER_WT_TM_LT  (weights.temp_lt = 1)
        SizeGt: 1.0,         // PRIMER_WT_SIZE_GT (weights.length_gt = 1)
        SizeLt: 1.0,         // PRIMER_WT_SIZE_LT (weights.length_lt = 1)
        GcGt: 0.0,           // PRIMER_WT_GC_PERCENT_GT (weights.gc_content_gt = 0)
        GcLt: 0.0,           // PRIMER_WT_GC_PERCENT_LT (weights.gc_content_lt = 0)
        SelfAny: 0.0,        // PRIMER_WT_SELF_ANY (weights.compl_any = 0)
        SelfEnd: 0.0,        // PRIMER_WT_SELF_END (weights.compl_end = 0)
        NumNs: 0.0,          // PRIMER_WT_NUM_NS  (weights.num_ns = 0)
        SelfAnyTh: 0.0,      // PRIMER_WT_SELF_ANY_TH (weights.compl_any_th = 0)
        SelfEndTh: 0.0,      // PRIMER_WT_SELF_END_TH (weights.compl_end_th = 0)
        HairpinTh: 0.0,      // PRIMER_WT_HAIRPIN_TH  (weights.hairpin_th = 0)
        EndStability: 0.0,   // PRIMER_WT_END_STABILITY (weights.end_stability = 0)
        ThermodynamicOligoAlignment: true // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1 (pr_set_default_global_args_2)
    );

    /// <summary>
    /// Primer3's fixed <c>weights.temp_cutoff</c> (= 5 °C, <c>pr_set_default_global_args_1</c>; not a
    /// user-settable tag). In thermodynamic mode a secondary-structure Tm within this many degrees of
    /// the primer Tm is penalised linearly, otherwise by a reciprocal term.
    /// </summary>
    public const double Primer3TempCutoff = 5.0;

    /// <summary>
    /// Primer3 default per-primer optima. OPT_TM = 60 °C and OPT_SIZE = 20 bases are
    /// from <c>libprimer3.c</c> (opt_tm = 60.0, opt_size = 20); OPT_GC_PERCENT is undefined (null) as in Primer3's code
    /// (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>; the manual's 50 is not applied): a non-zero GC weight
    /// then needs an explicit optimum (Primer3 <c>_pr_data_control</c>).
    /// <see cref="EvaluatePrimer"/> / <see cref="DesignPrimers"/> take the GC optimum from
    /// <see cref="PrimerParameters.OptimalGcPercent"/>.
    /// </summary>
    public static readonly Primer3Optima DefaultPrimer3Optima = new(
        OptTm: 60.0,         // PRIMER_OPT_TM (libprimer3.cc: opt_tm = 60.0) °C
        OptSize: 20,         // PRIMER_OPT_SIZE (libprimer3.cc: opt_size = 20) bases
        OptGcPercent: null   // PRIMER_OPT_GC_PERCENT (libprimer3.c: PR_UNDEFINED_INT_OPT) %
    );

    /// <summary>
    /// Computes the Primer3 per-primer penalty (objective function value) for a single
    /// primer, faithfully reproducing the left/right-primer branch of Primer3's
    /// <c>p_obj_fn</c>. The penalty is the weighted sum of one-sided deviations of Tm,
    /// length and GC% from their optima, plus the weighted secondary-structure terms
    /// (alignment-score mode: <c>compl_any</c>/<c>compl_end</c>; thermodynamic mode, Primer3's
    /// default: <c>compl_any_th</c>/<c>compl_end_th</c>/<c>hairpin_th</c> with the fixed 5 °C
    /// <c>temp_cutoff</c>), the number-of-Ns term and the 3'-end-stability term. Each term is
    /// added only when its weight is non-zero (and, for Tm/GC/size, the deviation has the
    /// matching sign), so the result is always ≥ 0; <b>lower is better</b>, exactly as Primer3
    /// sorts candidates. Cross-checked against primer3-py 2.3.1 <c>design_primers</c>
    /// PRIMER_LEFT/RIGHT_n_PENALTY in both alignment modes.
    /// The library mispriming term <c>repeat_sim</c> is PRIMER_WT_LIBRARY_MISPRIMING
    /// (<see cref="Primer3PenaltyWeights.LibraryMispriming"/>) × <see cref="Primer3PenaltyInputs.LibraryMispriming"/>.
    /// The fraction-bound terms PRIMER_WT_BOUND_GT / _LT (<see cref="Primer3PenaltyWeights.BoundGt"/> /
    /// <see cref="Primer3PenaltyWeights.BoundLt"/> around <see cref="Primer3Optima.OptBound"/>) apply when
    /// <see cref="Primer3PenaltyInputs.Bound"/> is set (Primer3: PRIMER_ANNEALING_TEMP &gt; 0); the position term is
    /// PRIMER_WT_POS_PENALTY × <see cref="Primer3PenaltyInputs.PositionPenalty"/> (PRIMER_INSIDE/OUTSIDE_PENALTY).
    /// The sequence-quality term is PRIMER_WT_SEQ_QUAL (<see cref="Primer3PenaltyWeights.SequenceQuality"/>) ×
    /// (<see cref="Primer3PenaltyInputs.QualityRangeMax"/> − <see cref="Primer3PenaltyInputs.SequenceQuality"/>).
    /// The masker term is PRIMER_WT_MASK_FAILURE_RATE (<see cref="Primer3PenaltyWeights.MaskFailureRate"/>) ×
    /// <see cref="Primer3PenaltyInputs.MaskFailureRate"/>. <para>The template mispriming terms (PRIMER_WT_TEMPLATE_MISPRIMING / _TH,
    /// <see cref="Primer3PenaltyWeights.TemplateMispriming"/> / <see cref="Primer3PenaltyWeights.TemplateMisprimingTh"/>)
    /// use <see cref="Primer3PenaltyInputs.TemplateMispriming"/>.</para>
    /// </summary>
    /// <param name="inputs">Measured primer properties (Tm in °C, length in bases, GC in
    /// percent 0–100, self/3' local-alignment scores, count of N bases).</param>
    /// <param name="weights">Objective weights; defaults to <see cref="DefaultPrimer3Weights"/>.</param>
    /// <param name="optima">Parameter optima; defaults to <see cref="DefaultPrimer3Optima"/>.</param>
    /// <returns>The Primer3 objective-function value (penalty); 0 means every term is at its optimum.</returns>
    public static double CalculatePrimer3Penalty(
        Primer3PenaltyInputs inputs,
        Primer3PenaltyWeights? weights = null,
        Primer3Optima? optima = null)
    {
        var w = weights ?? DefaultPrimer3Weights;
        var o = optima ?? DefaultPrimer3Optima;

        double sum = 0.0;

        // Tm term: one-sided, separate _gt / _lt weights (p_obj_fn temp_gt / temp_lt).
        if (w.TmGt != 0 && inputs.Tm > o.OptTm)
            sum += w.TmGt * (inputs.Tm - o.OptTm);
        if (w.TmLt != 0 && inputs.Tm < o.OptTm)
            sum += w.TmLt * (o.OptTm - inputs.Tm);

        // Fraction-bound terms (bound_gt / bound_lt around PRIMER_OPT_BOUND); for primers only when
        // PRIMER_ANNEALING_TEMP > 0 (inputs.Bound is then set).
        if (inputs.Bound is { } bound)
        {
            if (w.BoundGt != 0 && bound > o.OptBound)
                sum += w.BoundGt * (bound - o.OptBound);
            if (w.BoundLt != 0 && bound < o.OptBound)
                sum += w.BoundLt * (o.OptBound - bound);
        }

        // GC% term (gc_content is a percentage 0–100 in libprimer3.c); a GC weight needs a defined optimum
        // (_pr_data_control: "Primer GC content is part of objective function while optimum gc_content is not defined").
        if (w.GcGt != 0 || w.GcLt != 0)
        {
            double optGc = o.OptGcPercent ?? throw new ArgumentException(PrimerGcOptimumUndefinedMessage, nameof(optima));
            if (w.GcGt != 0 && inputs.GcPercent > optGc)
                sum += w.GcGt * (inputs.GcPercent - optGc);
            if (w.GcLt != 0 && inputs.GcPercent < optGc)
                sum += w.GcLt * (optGc - inputs.GcPercent);
        }

        // Length/size term (p_obj_fn length_lt / length_gt).
        if (w.SizeLt != 0 && inputs.Length < o.OptSize)
            sum += w.SizeLt * (o.OptSize - inputs.Length);
        if (w.SizeGt != 0 && inputs.Length > o.OptSize)
            sum += w.SizeGt * (inputs.Length - o.OptSize);

        // Masker failure-rate term (failure_rate; left/right primers, PRIMER_MASK_TEMPLATE).
        if (w.MaskFailureRate != 0)
            sum += w.MaskFailureRate * inputs.MaskFailureRate;

        // Secondary-structure terms: p_obj_fn switches on thermodynamic_oligo_alignment.
        if (!w.ThermodynamicOligoAlignment)
        {
            // Mode 0: local-alignment scores, linear weights (compl_any, compl_end).
            if (w.SelfAny != 0)
                sum += w.SelfAny * inputs.SelfAny;
            if (w.SelfEnd != 0)
                sum += w.SelfEnd * inputs.SelfEnd;
        }
        else
        {
            // Mode 1 (Primer3 default): SelfAny / SelfEnd / HairpinTh are structure Tm values (°C)
            // (compl_any_th, compl_end_th, hairpin_th).
            sum += ThermodynamicStructurePenalty(w.SelfAnyTh, inputs.Tm, inputs.SelfAny);
            sum += ThermodynamicStructurePenalty(w.SelfEndTh, inputs.Tm, inputs.SelfEnd);
            sum += ThermodynamicStructurePenalty(w.HairpinTh, inputs.Tm, inputs.HairpinTh);
        }

        // Number-of-Ns term (num_ns).
        if (w.NumNs != 0)
            sum += w.NumNs * inputs.NumNs;

        // Library mispriming term (repeat_sim): weight · repeat_sim.score[repeat_sim.max].
        if (w.LibraryMispriming != 0)
            sum += w.LibraryMispriming * inputs.LibraryMispriming;

        // Position term (pos_penalty): weight · position penalty relative to the target.
        if (w.PositionPenalty != 0)
            sum += w.PositionPenalty * inputs.PositionPenalty;

        // 3'-end stability term (end_stability): weight · ΔG magnitude (kcal/mol, as Primer3 reports it).
        if (w.EndStability != 0)
            sum += w.EndStability * inputs.EndStability;

        // Sequence quality term (seq_quality): weight · (PRIMER_QUALITY_RANGE_MAX − h->seq_quality); 0 without quality
        // data (seq_quality = range max). Internal oligos: the last term of the OT_INTL branch (the end-stability,
        // position and template inputs are 0 there). PRIMER_WT_END_QUAL is never read by p_obj_fn.
        if (w.SequenceQuality != 0)
            sum += w.SequenceQuality * (inputs.QualityRangeMax - (inputs.SequenceQuality ?? inputs.QualityRangeMax));

        // Template mispriming terms (after seq_quality): alignment mode linear
        // (weights.template_mispriming), thermodynamic mode the temp_cutoff rule (weights.template_mispriming_th).
        if (!w.ThermodynamicTemplateAlignment)
        {
            if (w.TemplateMispriming != 0)
                sum += w.TemplateMispriming * inputs.TemplateMispriming;
        }
        else
            sum += ThermodynamicStructurePenalty(w.TemplateMisprimingTh, inputs.Tm, inputs.TemplateMispriming);

        return sum;
    }

    // p_obj_fn thermodynamic secondary-structure term (libprimer3.cc):
    //   if (Tm − temp_cutoff) ≤ s : w · (s − (Tm − temp_cutoff − 1))
    //   else                       : w · 1 / (Tm − temp_cutoff + 1 − s)
    private static double ThermodynamicStructurePenalty(double weight, double primerTm, double structureTm)
    {
        if (weight == 0)
            return 0;
        double threshold = primerTm - Primer3TempCutoff;
        return threshold <= structureTm
            ? weight * (structureTm - (threshold - 1.0))
            : weight * (1.0 / (threshold + 1.0 - structureTm));
    }

    private static double CalculatePrimerScore(string seq, double gc, double tm, int homopolymer, PrimerParameters param)
    {
        double score = 100;

        // Penalize for deviation from optimal length
        score -= Math.Abs(seq.Length - param.OptimalLength) * 2;

        // Penalize for deviation from optimal Tm
        score -= Math.Abs(tm - param.OptimalTm) * 2;

        // Penalize for deviation from 50% GC
        score -= Math.Abs(gc - 50) * 0.5;

        // Penalize for homopolymers
        score -= homopolymer * 5;

        // Bonus for GC clamp at 3' end
        if (seq.Length >= 2)
        {
            char last = seq[^1];
            if (last == 'G' || last == 'C')
                score += 5;
        }

        return Math.Max(0, score);
    }

    private static string Reverse(string s)
    {
        var chars = s.ToCharArray();
        Array.Reverse(chars);
        return new string(chars);
    }

    private static bool AreComplementary(string s1, string s2)
    {
        if (s1.Length != s2.Length) return false;
        for (int i = 0; i < s1.Length; i++)
        {
            if (!IsComplementary(s1[i], s2[i]))
                return false;
        }
        return true;
    }

    private static bool IsComplementary(char c1, char c2) =>
        (c1 == 'A' && c2 == 'T') || (c1 == 'T' && c2 == 'A') ||
        (c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G');
}

/// <summary>
/// Parameters for primer design. <c>MaxDinucleotideRepeats</c> limits the unsourced library screen
/// <see cref="PrimerDesigner.FindLongestDinucleotideRepeat"/> (no Primer3 counterpart; <see cref="int.MaxValue"/> disables
/// it, as in <see cref="PrimerDesigner.Primer3DefaultParameters"/>). Primer3's 3′-end checks are <see cref="GcClamp"/> (PRIMER_GC_CLAMP),
/// <see cref="MaxEndGc"/> (PRIMER_MAX_END_GC) and <see cref="MaxEndStability"/> (PRIMER_MAX_END_STABILITY).
/// Two positional members predate them and are kept only for source compatibility (deprecated):
/// <c>Avoid3PrimeGC</c> is a library rule, not Primer3's — despite its name it <i>requires</i> at least one G/C
/// among the two 3′-most bases (use <see cref="GcClamp"/> = 1 for Primer3's "3′-most base is G/C", or
/// <see cref="MaxEndGc"/> to limit 3′ G/C); <c>Check3PrimeStability</c> no longer has any effect — it gated
/// ΔG(3′ pentamer) &lt; −9 kcal/mol, which no ACGT pentamer reaches (minimum −6.86, GCGCG), so it never rejected a
/// primer; Primer3's end-stability limit is <see cref="MaxEndStability"/>, always applied.
/// </summary>
public readonly record struct PrimerParameters(
    int MinLength,
    int MaxLength,
    int OptimalLength,
    double MinGcContent,
    double MaxGcContent,
    double MinTm,
    double MaxTm,
    double OptimalTm,
    int MaxHomopolymer,
    int MaxDinucleotideRepeats,
    bool Avoid3PrimeGC,
    bool Check3PrimeStability,
    PrimerStructureScreen StructureScreen = PrimerStructureScreen.Primer3Thermodynamic,
    double MaxStructureTm = PrimerDesigner.Primer3MaxStructureTm)
{
    /// <summary>
    /// Limit (°C) for every Primer3 thermodynamic structure value under
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> — Primer3's PRIMER_MAX_SELF_ANY_TH,
    /// PRIMER_MAX_SELF_END_TH, PRIMER_MAX_HAIRPIN_TH, PRIMER_PAIR_MAX_COMPL_ANY_TH and
    /// PRIMER_PAIR_MAX_COMPL_END_TH, all 47 °C by default. A value strictly greater fails.
    /// The value <c>0</c> (e.g. from <c>default(PrimerParameters)</c>) also means 47 °C.
    /// </summary>
    public double EffectiveMaxStructureTm => MaxStructureTm > 0 ? MaxStructureTm : PrimerDesigner.Primer3MaxStructureTm;

    /// <summary>
    /// PRIMER_MAX_SELF_ANY under <see cref="PrimerStructureScreen.Primer3Alignment"/>: maximum Primer3 alignment-mode
    /// <c>self_any</c> (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>); null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxSelfAny"/> (8.00). A value strictly greater fails. Must be in [0, 32767].
    /// </summary>
    public double? MaxSelfAny { get; init; }

    /// <summary>
    /// PRIMER_MAX_SELF_END under <see cref="PrimerStructureScreen.Primer3Alignment"/>: maximum Primer3 alignment-mode
    /// <c>self_end</c> (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>); null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxSelfEnd"/> (3.00). As in Primer3's <c>characterize_pair</c>, it also bounds the
    /// pair's reverse-orientation 3′ complementarity <c>align(right, revcomp(left), DPAL_GLOBAL_END)</c>.
    /// </summary>
    public double? MaxSelfEnd { get; init; }

    /// <summary>
    /// Primer3 per-primer penalty weights (PRIMER_WT_*) used for <see cref="PrimerCandidate.Penalty"/>; null =
    /// <see cref="PrimerDesigner.DefaultPrimer3Weights"/>. The secondary-structure mode of the weights follows
    /// <see cref="PrimerParameters.StructureScreen"/> (Primer3 PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT):
    /// <see cref="PrimerStructureScreen.Primer3Alignment"/> uses PRIMER_WT_SELF_ANY / _SELF_END
    /// (<see cref="Primer3PenaltyWeights.SelfAny"/>/<see cref="Primer3PenaltyWeights.SelfEnd"/>) × the dpal scores,
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> PRIMER_WT_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH
    /// with the ntthal Tm values; the <see cref="PrimerStructureScreen.Heuristic"/> screen has no Primer3 structure
    /// values, so its structure terms are 0. When a structure weight is non-zero the values are computed for every
    /// candidate (as Primer3's <c>calc_and_check_oligo_features</c> does). PRIMER_WT_END_STABILITY
    /// (<see cref="Primer3PenaltyWeights.EndStability"/>) multiplies Primer3's <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="PrimerDesigner.Calculate3PrimeStability"/> (positive kcal/mol), as
    /// <c>p_obj_fn</c> does for left/right primers (internal oligos have no such term).
    /// </summary>
    public Primer3PenaltyWeights? PenaltyWeights { get; init; }

    /// <summary>
    /// PRIMER_SALT_MONOVALENT: monovalent cation concentration (mM, &gt; 0) of the primer reaction; null = Primer3's
    /// default <see cref="PrimerDesigner.Primer3MonovalentMillimolar"/> (50 mM). Primer3 uses the primer conditions
    /// (<c>p_args</c>) for the primer Tm (<c>seqtm</c>), the ntthal self-dimer / 3′ self-dimer / hairpin and pair
    /// complementarity values (<c>create_thal_arg_holder(p_args)</c>) and the product Tm (<c>long_seq_tm</c>).
    /// The probe-side counterpart is <see cref="ProbeDesigner.ProbeParameters.MonovalentMillimolar"/>; the internal oligo
    /// of a pair uses <see cref="ProbeDesigner.Primer3ProbeSettings.MonovalentMillimolar"/> (PRIMER_INTERNAL_SALT_MONOVALENT).
    /// </summary>
    public double? MonovalentMillimolar { get; init; }

    /// <summary>PRIMER_SALT_DIVALENT: Mg²⁺ concentration (mM, ≥ 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DivalentMillimolar"/> (1.5 mM).</summary>
    public double? DivalentMillimolar { get; init; }

    /// <summary>PRIMER_DNTP_CONC: dNTP concentration (mM, ≥ 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DntpMillimolar"/> (0.6 mM).</summary>
    public double? DntpMillimolar { get; init; }

    /// <summary>PRIMER_DNA_CONC: primer (oligo) concentration (nM, &gt; 0); null = Primer3's default <see cref="PrimerDesigner.Primer3DnaConcentrationNanomolar"/> (50 nM).</summary>
    public double? DnaConcentrationNanomolar { get; init; }

    /// <summary>
    /// PRIMER_OPT_GC_PERCENT: the GC optimum of the PRIMER_WT_GC_PERCENT_GT/_LT penalty terms
    /// (<see cref="Primer3PenaltyWeights.GcGt"/>/<see cref="Primer3PenaltyWeights.GcLt"/>); null =
    /// undefined, as in Primer3's code (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>; the manual's 50 is not
    /// applied). It is inert while both GC weights are 0 (Primer3's default); a non-zero GC weight without it is rejected with
    /// Primer3's <c>_pr_data_control</c> error ("Primer GC content is part of objective function while optimum gc_content is
    /// not defined", <see cref="ArgumentException"/>).
    /// </summary>
    public double? OptimalGcPercent { get; init; }

    /// <summary>
    /// PRIMER_MAX_END_STABILITY (kcal/mol, ≥ 0): a primer fails when its Primer3 <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="PrimerDesigner.Calculate3PrimeStability"/> (the positive −ΔG°37 of the
    /// 3′ pentamer) is strictly greater; null = Primer3's default <see cref="PrimerDesigner.Primer3MaxEndStability"/>
    /// (100, never exceeded). Primer3 applies it to left/right primers only (no internal-oligo counterpart).
    /// </summary>
    public double? MaxEndStability { get; init; }

    /// <summary>
    /// PRIMER_GC_CLAMP: number of consecutive G/C required at the 3′ end of a primer (0 = none, Primer3's default);
    /// must not exceed <see cref="MinLength"/> (Primer3 <c>_pr_data_control</c>). A negative value requires nothing,
    /// as in Primer3. Left/right primers only.
    /// </summary>
    public int GcClamp { get; init; }

    /// <summary>
    /// PRIMER_MAX_END_GC (0–5): maximum number of G/C among the five 3′-most bases of a primer; null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxEndGc"/> (5, not checked). Left/right primers only.
    /// </summary>
    public int? MaxEndGc { get; init; }

    /// <summary>
    /// PRIMER_MISPRIMING_LIBRARY (primer3-py <c>misprime_lib</c>): when set and non-empty, every primer gets Primer3's
    /// library mispriming score (<see cref="PrimerDesigner.CalculateLibraryMispriming"/>,
    /// <see cref="PrimerCandidate.LibraryMispriming"/>) and fails when any entry's weighted score exceeds
    /// <see cref="MaxLibraryMispriming"/>; pairs are limited by <see cref="PrimerPairOptions.MaxLibraryMispriming"/>.
    /// As in Primer3's <c>calc_and_check_oligo_features</c> / <c>characterize_pair</c>, the pair search scores primers
    /// when they are picked only if PRIMER_WT_LIBRARY_MISPRIMING (<see cref="Primer3PenaltyWeights.LibraryMispriming"/>)
    /// is non-zero, otherwise lazily when a pair is characterized (same results). Null = no library (Primer3's default).
    /// Not part of the JSON form of the parameters (the MCP tools take the library as a separate name → sequence argument).
    /// </summary>
    [System.Text.Json.Serialization.JsonIgnore]
    public PrimerMisprimingLibrary? MisprimingLibrary { get; init; }

    /// <summary>
    /// PRIMER_MAX_LIBRARY_MISPRIMING: maximum weighted library score of one primer; null = Primer3's default
    /// <see cref="PrimerDesigner.Primer3MaxLibraryMispriming"/> (12.00). Primer3 compares with the value truncated to a
    /// C <c>short</c> (12.9 acts as 12). Under <see cref="PrimerStructureScreen.Primer3Alignment"/> it must not exceed
    /// 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double? MaxLibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS: false (null, Primer3's default 0 — <c>pr_set_default_global_args_2</c>, used by
    /// primer3-py; the "version 1" defaults had 1) = an IUPAC code in a library entry never aligns and N scores −0.25;
    /// true (1) = an IUPAC code matches every base it represents (so runs of N match any primer).
    /// </summary>
    public bool? LibraryAmbiguityCodesConsensus { get; init; }

    /// <summary>
    /// PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT: false (Primer3's default 0) = template mispriming is scored with dpal
    /// (<see cref="MaxTemplateMispriming"/>, PRIMER_WT_TEMPLATE_MISPRIMING); true (1) = with the ntthal THAL_END1 Tm
    /// (<see cref="MaxTemplateMisprimingTh"/>, PRIMER_WT_TEMPLATE_MISPRIMING_TH; template ≤ 10000 nt). Independent of
    /// <see cref="StructureScreen"/> (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT). See
    /// <see cref="PrimerDesigner.CalculateTemplateMispriming"/>.
    /// </summary>
    public bool ThermodynamicTemplateAlignment { get; init; }

    /// <summary>
    /// PRIMER_MAX_TEMPLATE_MISPRIMING (alignment mode): in <see cref="PrimerDesigner.DesignPrimers"/> /
    /// <see cref="PrimerDesigner.DesignPrimerPairs"/> a primer whose template mispriming score
    /// (<see cref="PrimerDesigner.CalculateTemplateMispriming"/>, <see cref="PrimerCandidate.TemplateMispriming"/>) is
    /// strictly greater fails. Null = Primer3's default <see cref="PrimerDesigner.Primer3UndefinedTemplateMispriming"/>
    /// (−100); a negative value is not checked. Must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double? MaxTemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_MAX_TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C): as <see cref="MaxTemplateMispriming"/> for the ntthal
    /// template-mispriming Tm; null = −100 (not checked).
    /// </summary>
    public double? MaxTemplateMisprimingTh { get; init; }

    /// <summary>
    /// PRIMER_ANNEALING_TEMP (°C; null = Primer3's default −10 = off). When &gt; 0 Primer3 computes each primer's
    /// fraction bound at this temperature (<see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/>, reported in
    /// <see cref="PrimerCandidate.Bound"/>), rejects primers outside [<see cref="MinBound"/>, <see cref="MaxBound"/>] and
    /// adds the PRIMER_WT_BOUND_GT / _LT terms (<see cref="Primer3PenaltyWeights.BoundGt"/> /
    /// <see cref="Primer3PenaltyWeights.BoundLt"/>) around <see cref="OptBound"/>. Must not exceed 100 °C. In a primer-pair
    /// design it is also the internal oligo's annealing temperature (one global Primer3 setting).
    /// </summary>
    public double? AnnealingTemperature { get; init; }

    /// <summary>PRIMER_MIN_BOUND (% bound; null = Primer3's −10); checked only when <see cref="AnnealingTemperature"/> &gt; 0.
    /// A primer longer than 36 bases has no bound value (Primer3 OLIGOTM_ERROR = −999999.9999) and fails it.</summary>
    public double? MinBound { get; init; }

    /// <summary>PRIMER_MAX_BOUND (% bound; null = Primer3's 110); checked only when <see cref="AnnealingTemperature"/> &gt; 0.</summary>
    public double? MaxBound { get; init; }

    /// <summary>PRIMER_OPT_BOUND (% bound; null = Primer3's 97): the optimum of the PRIMER_WT_BOUND_GT/_LT terms; must lie
    /// in [<see cref="MinBound"/>, <see cref="MaxBound"/>] (Primer3 <c>_pr_data_control</c>).</summary>
    public double? OptBound { get; init; }

    /// <summary>
    /// PRIMER_MIN_QUALITY (default 0): with SEQUENCE_QUALITY (<see cref="PrimerPairOptions.SequenceQuality"/>) a primer
    /// whose minimum base quality (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>,
    /// <see cref="PrimerCandidate.MinSequenceQuality"/>) is lower fails (Primer3 <c>sequence_quality_is_ok</c>). A non-zero
    /// value requires quality data and must lie in [<see cref="QualityRangeMin"/>, <see cref="QualityRangeMax"/>]
    /// (Primer3 <c>_pr_data_control</c>; <see cref="PrimerDesigner.EvaluatePrimer"/> has no quality data, so it rejects it).
    /// </summary>
    public int MinQuality { get; init; }

    /// <summary>
    /// PRIMER_MIN_END_QUALITY (default 0): with SEQUENCE_QUALITY a primer whose minimum quality over its five 3′-most bases
    /// is lower fails (checked after <see cref="MinQuality"/>; left/right primers only). Not range-checked by Primer3.
    /// </summary>
    public int MinEndQuality { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MIN (null = Primer3's 0): every SEQUENCE_QUALITY value must be ≥ it. One global
    /// Primer3 setting: it also applies to the internal oligo of a pair.</summary>
    public int? QualityRangeMin { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MAX (null = Primer3's 100): every SEQUENCE_QUALITY value must be ≤ it; the PRIMER_WT_SEQ_QUAL
    /// term is weight × (QualityRangeMax − min quality). One global Primer3 setting (also used for the internal oligo).</summary>
    public int? QualityRangeMax { get; init; }

    /// <summary>Effective PRIMER_QUALITY_RANGE_MIN.</summary>
    public int EffectiveQualityRangeMin => QualityRangeMin ?? PrimerDesigner.Primer3QualityRangeMin;

    /// <summary>Effective PRIMER_QUALITY_RANGE_MAX.</summary>
    public int EffectiveQualityRangeMax => QualityRangeMax ?? PrimerDesigner.Primer3QualityRangeMax;

    /// <summary>PRIMER_ANNEALING_TEMP in effect (Primer3 default −10 = off).</summary>
    public double EffectiveAnnealingTemperature => AnnealingTemperature ?? PrimerDesigner.Primer3DefaultAnnealingTemperature;

    /// <summary>PRIMER_MIN_BOUND in effect (default −10).</summary>
    public double EffectiveMinBound => MinBound ?? PrimerDesigner.Primer3MinBound;

    /// <summary>PRIMER_MAX_BOUND in effect (default 110).</summary>
    public double EffectiveMaxBound => MaxBound ?? PrimerDesigner.Primer3MaxBound;

    /// <summary>PRIMER_OPT_BOUND in effect (default 97).</summary>
    public double EffectiveOptBound => OptBound ?? PrimerDesigner.Primer3OptBound;

    /// <summary>Effective PRIMER_MAX_TEMPLATE_MISPRIMING.</summary>
    public double EffectiveMaxTemplateMispriming => MaxTemplateMispriming ?? PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>Effective PRIMER_MAX_TEMPLATE_MISPRIMING_TH.</summary>
    public double EffectiveMaxTemplateMisprimingTh => MaxTemplateMisprimingTh ?? PrimerDesigner.Primer3UndefinedTemplateMispriming;

    // The per-primer limit and weight of the active template alignment mode.
    internal double ActiveMaxTemplateMispriming =>
        ThermodynamicTemplateAlignment ? EffectiveMaxTemplateMisprimingTh : EffectiveMaxTemplateMispriming;

    internal double ActiveTemplateMisprimingWeight =>
        ThermodynamicTemplateAlignment
            ? PenaltyWeights?.TemplateMisprimingTh ?? 0.0
            : PenaltyWeights?.TemplateMispriming ?? 0.0;

    /// <summary>Effective PRIMER_MAX_LIBRARY_MISPRIMING.</summary>
    public double EffectiveMaxLibraryMispriming => MaxLibraryMispriming ?? PrimerDesigner.Primer3MaxLibraryMispriming;

    /// <summary>Effective PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS.</summary>
    public bool EffectiveLibraryAmbiguityCodesConsensus => LibraryAmbiguityCodesConsensus ?? false;

    // The library in use: null when absent or empty (Primer3 seq_lib_num_seq == 0).
    internal PrimerMisprimingLibrary? ActiveLibrary => MisprimingLibrary is { Count: > 0 } lib ? lib : null;

    /// <summary>Effective PRIMER_MAX_END_STABILITY (kcal/mol).</summary>
    public double EffectiveMaxEndStability => MaxEndStability ?? PrimerDesigner.Primer3MaxEndStability;

    /// <summary>Effective PRIMER_MAX_END_GC.</summary>
    public int EffectiveMaxEndGc => MaxEndGc ?? PrimerDesigner.Primer3MaxEndGc;

    /// <summary>Effective PRIMER_SALT_MONOVALENT (mM).</summary>
    public double EffectiveMonovalentMillimolar => MonovalentMillimolar ?? PrimerDesigner.Primer3MonovalentMillimolar;

    /// <summary>Effective PRIMER_SALT_DIVALENT (mM).</summary>
    public double EffectiveDivalentMillimolar => DivalentMillimolar ?? PrimerDesigner.Primer3DivalentMillimolar;

    /// <summary>Effective PRIMER_DNTP_CONC (mM).</summary>
    public double EffectiveDntpMillimolar => DntpMillimolar ?? PrimerDesigner.Primer3DntpMillimolar;

    /// <summary>Effective PRIMER_DNA_CONC (nM).</summary>
    public double EffectiveDnaConcentrationNanomolar => DnaConcentrationNanomolar ?? PrimerDesigner.Primer3DnaConcentrationNanomolar;

    /// <summary>Effective PRIMER_MAX_SELF_ANY (<see cref="MaxSelfAny"/> or 8.00).</summary>
    public double EffectiveMaxSelfAny => MaxSelfAny ?? PrimerDesigner.Primer3MaxSelfAny;

    /// <summary>Effective PRIMER_MAX_SELF_END (<see cref="MaxSelfEnd"/> or 3.00).</summary>
    public double EffectiveMaxSelfEnd => MaxSelfEnd ?? PrimerDesigner.Primer3MaxSelfEnd;

    // Primer3 _pr_data_control: PRIMER_WT_GC_PERCENT_GT/_LT need PRIMER_OPT_GC_PERCENT.
    internal void ValidateGcOptimum(string paramName)
    {
        if ((PenaltyWeights is { } w && (w.GcGt != 0 || w.GcLt != 0)) && OptimalGcPercent is null)
            throw new ArgumentException(PrimerDesigner.PrimerGcOptimumUndefinedMessage, paramName);
    }

    // Primer3 _pr_data_control checks of the primer conditions.
    internal void ValidateConditions(string paramName)
    {
        PrimerDesigner.ValidatePrimer3Conditions(EffectiveMonovalentMillimolar, EffectiveDivalentMillimolar,
            EffectiveDntpMillimolar, EffectiveDnaConcentrationNanomolar, paramName);
        if (OptimalGcPercent is { } optGc && !double.IsFinite(optGc))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_OPT_GC_PERCENT must be finite.");
        // _pr_data_control: PRIMER_MAX_END_GC must be between 0 to 5; PRIMER_MAX_END_STABILITY must be non-negative;
        // PRIMER_GC_CLAMP > PRIMER_MIN_SIZE.
        if (EffectiveMaxEndGc is < 0 or > PrimerDesigner.Primer3MaxEndGc)
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_MAX_END_GC must be between 0 to 5.");
        if (!(EffectiveMaxEndStability >= 0))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_MAX_END_STABILITY must be non-negative.");
        if (GcClamp > MinLength)
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_GC_CLAMP > PRIMER_MIN_SIZE.");
        // _pr_data_control: PRIMER_MAX_LIBRARY_MISPRIMING > SHRT_MAX (alignment mode); a library weight without a library.
        double maxLib = EffectiveMaxLibraryMispriming;
        if (double.IsNaN(maxLib) || Math.Abs(maxLib) >= int.MaxValue
            || (maxLib > short.MaxValue && StructureScreen == PrimerStructureScreen.Primer3Alignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_MAX_LIBRARY_MISPRIMING.");
        if ((PenaltyWeights?.LibraryMispriming ?? 0) != 0 && ActiveLibrary is null)
            throw new ArgumentException("Mispriming score is part of objective function, but mispriming library is not defined (Primer3 _pr_data_control).", paramName);
        // _pr_data_control: PRIMER_MAX_TEMPLATE_MISPRIMING > SHRT_MAX in alignment mode.
        if (double.IsNaN(EffectiveMaxTemplateMispriming) || double.IsNaN(EffectiveMaxTemplateMisprimingTh)
            || (EffectiveMaxTemplateMispriming > short.MaxValue && !ThermodynamicTemplateAlignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_MAX_TEMPLATE_MISPRIMING.");
        // A negative template weight makes Primer3 skip the score (_pr_need_template_mispriming: weight > 0) but still
        // add the term (p_obj_fn: weight ≠ 0), which fails its PR_ASSERT.
        if (!((PenaltyWeights?.TemplateMispriming ?? 0) >= 0) || !((PenaltyWeights?.TemplateMisprimingTh ?? 0) >= 0))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_WT_TEMPLATE_MISPRIMING[_TH] must be ≥ 0.");
        // _pr_data_control: PRIMER_OPT_BOUND within [PRIMER_MIN_BOUND, PRIMER_MAX_BOUND]; PRIMER_ANNEALING_TEMP ≤ 100.
        PrimerDesigner.ValidatePrimer3Bound(EffectiveAnnealingTemperature, EffectiveMinBound, EffectiveMaxBound,
            EffectiveOptBound, internalOligo: false, paramName);
    }
}

/// <summary>
/// Secondary-structure screen applied by <see cref="PrimerDesigner.EvaluatePrimer"/> and
/// <see cref="PrimerDesigner.DesignPrimers"/>.
/// </summary>
public enum PrimerStructureScreen
{
    /// <summary>
    /// Primer3's default thermodynamic screen (<c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1</c>): a
    /// primer is rejected when its ntthal self-dimer, 3′ self-dimer or hairpin Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>), a pair when its ntthal
    /// hetero-dimer or 3′ hetero-dimer Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3PairComplementarity"/>).
    /// </summary>
    Primer3Thermodynamic = 0,

    /// <summary>
    /// Sequence-only library screen: <see cref="PrimerDesigner.HasHairpinPotential"/> (a ≥ 4-bp
    /// Watson–Crick stem closing a ≥ 3-nt loop; unsourced library heuristic, not Primer3's hairpin screen) per
    /// primer and <see cref="PrimerDesigner.HasPrimerDimer"/>
    /// (Primer3 alignment-mode pair 3′ complementarity ≥ 4) per pair.
    /// </summary>
    Heuristic = 1,

    /// <summary>
    /// Primer3's alignment-score screen (<c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>, <c>libprimer3.cc</c>
    /// <c>oligo_compl</c> / <c>characterize_pair</c>): a primer is rejected when its dpal self-complementarity
    /// <c>self_any</c> (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>) exceeds
    /// <see cref="PrimerParameters.MaxSelfAny"/> (PRIMER_MAX_SELF_ANY, default 8.00) or its 3′ self-complementarity
    /// <c>self_end</c> (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>) exceeds
    /// <see cref="PrimerParameters.MaxSelfEnd"/> (PRIMER_MAX_SELF_END, 3.00); a pair when <c>compl_any</c>
    /// (<see cref="PrimerDesigner.CalculatePrimerDimerAnyComplementarity"/>) exceeds
    /// <see cref="PrimerPairOptions.MaxComplAny"/> (PRIMER_PAIR_MAX_COMPL_ANY, 8.00) or <c>align(left, rc(right),
    /// GLOBAL_END)</c> exceeds <see cref="PrimerPairOptions.MaxComplEnd"/> (PRIMER_PAIR_MAX_COMPL_END, 3.00) or the
    /// reverse orientation <c>align(right, rc(left), GLOBAL_END)</c>, when larger, exceeds PRIMER_MAX_SELF_END
    /// (Primer3 compares that one with the per-primer limit). No hairpin value exists in this mode. An internal
    /// oligo (<see cref="PrimerPairOptions.PickInternalOligo"/>) is screened with PRIMER_INTERNAL_MAX_SELF_ANY /
    /// _SELF_END (12.00). Verified against primer3-py 2.3.1 <c>design_primers</c> with
    /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0.
    /// </summary>
    Primer3Alignment = 2,
}

/// <summary>
/// A primer candidate with quality metrics. <see cref="MeltingTemperature"/> is the Primer3-default
/// Tm rounded to 0.1 °C; <see cref="Score"/> is an informational heuristic (higher is better);
/// <see cref="Penalty"/> is the unrounded Primer3 per-primer penalty (lower is better) that
/// <see cref="PrimerDesigner.DesignPrimers"/> ranks by. Under the default
/// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> screen <see cref="SelfAnyTh"/>,
/// <see cref="SelfEndTh"/> and <see cref="HairpinTh"/> carry Primer3's PRIMER_*_SELF_ANY_TH /
/// _SELF_END_TH / _HAIRPIN_TH (°C) and <see cref="HasHairpin"/> means HairpinTh &gt; 47 °C; under the
/// heuristic screen they are <c>null</c> and <see cref="HasHairpin"/> is
/// <see cref="PrimerDesigner.HasHairpinPotential"/>.
/// </summary>
public sealed record PrimerCandidate(
    string Sequence,
    int Position,
    bool IsForward,
    int Length,
    double GcContent,
    double MeltingTemperature,
    int HomopolymerLength,
    bool HasHairpin,
    double Stability3Prime,
    bool IsValid,
    IReadOnlyList<string> Issues,
    double Score,
    double Penalty = 0.0,
    double? SelfAnyTh = null,
    double? SelfEndTh = null,
    double? HairpinTh = null)
{
    /// <summary>
    /// Primer3 alignment-mode PRIMER_*_SELF_ANY (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>);
    /// set under <see cref="PrimerStructureScreen.Primer3Alignment"/>, otherwise <c>null</c>.
    /// </summary>
    public double? SelfAny { get; init; }

    /// <summary>
    /// Primer3 alignment-mode PRIMER_*_SELF_END (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>);
    /// set under <see cref="PrimerStructureScreen.Primer3Alignment"/>, otherwise <c>null</c>.
    /// </summary>
    public double? SelfEnd { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING score (<see cref="PrimerDesigner.CalculateLibraryMispriming"/>) when
    /// <see cref="PrimerParameters.MisprimingLibrary"/> is set, otherwise <c>null</c>.
    /// </summary>
    public double? LibraryMispriming { get; init; }

    /// <summary>The library entry named in PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING (with <see cref="LibraryMispriming"/>).</summary>
    public string? LibraryMisprimingName { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING (alignment mode) or _TEMPLATE_MISPRIMING_TH (thermodynamic mode,
    /// °C) — <see cref="TemplateMisprimingScore.Max"/> of <see cref="PrimerDesigner.CalculateTemplateMispriming"/> —
    /// for a primer of a pair returned by <see cref="PrimerDesigner.DesignPrimers"/> / <see cref="PrimerDesigner.DesignPrimerPairs"/>
    /// when Primer3 computes it (a template-mispriming limit or weight is set); otherwise <c>null</c>.
    /// </summary>
    public double? TemplateMispriming { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_BOUND: the fraction (%) bound at <see cref="PrimerParameters.AnnealingTemperature"/>
    /// (<see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/>) when that is &gt; 0 and the value is defined
    /// (≤ 36 ACGT bases); otherwise <c>null</c>.
    /// </summary>
    public double? Bound { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_POSITION_PENALTY (<see cref="PrimerDesigner.CalculatePositionPenaltyPrimer3"/>) of a
    /// designed primer when <see cref="PrimerPairOptions.InsidePenalty"/> / <see cref="PrimerPairOptions.OutsidePenalty"/>
    /// are not Primer3's defaults (−1 / 0); otherwise <c>null</c>.
    /// </summary>
    public double? PositionPenalty { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_MIN_SEQ_QUALITY: the primer's minimum base quality
    /// (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>) when SEQUENCE_QUALITY
    /// (<see cref="PrimerPairOptions.SequenceQuality"/>) is given; otherwise <c>null</c>.
    /// </summary>
    public int? MinSequenceQuality { get; init; }

    /// <summary>
    /// The primer's predicted PCR failure rate (<see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>, the value Primer3
    /// weights with PRIMER_WT_MASK_FAILURE_RATE) when <see cref="PrimerPairOptions.MaskTemplate"/> is set; otherwise <c>null</c>.
    /// </summary>
    public double? MaskFailureRate { get; init; }
}

/// <summary>
/// Measured properties of a single primer used as input to the Primer3 penalty
/// objective (<see cref="PrimerDesigner.CalculatePrimer3Penalty"/>). Units mirror
/// Primer3's <c>p_obj_fn</c>: Tm in °C, length in bases, GC in percent (0–100),
/// self/3'-complementarity as local-alignment scores, and the count of N bases.
/// </summary>
/// <param name="Tm">Primer melting temperature in °C.</param>
/// <param name="Length">Primer length in bases.</param>
/// <param name="GcPercent">GC content as a percentage in [0, 100].</param>
/// <param name="SelfAny">Self-complementarity: the local-alignment score PRIMER_SELF_ANY when
/// <see cref="Primer3PenaltyWeights.ThermodynamicOligoAlignment"/> is false, or the self-dimer Tm in °C
/// (PRIMER_SELF_ANY_TH) when it is true — Primer3 stores both in the same field.</param>
/// <param name="SelfEnd">3'-self-complementarity: PRIMER_SELF_END score (alignment mode) or
/// PRIMER_SELF_END_TH Tm in °C (thermodynamic mode).</param>
/// <param name="NumNs">Number of ambiguous N bases in the primer.</param>
/// <param name="HairpinTh">Hairpin Tm in °C (PRIMER_HAIRPIN_TH); used only in thermodynamic mode.</param>
/// <param name="EndStability">3'-end stability ΔG magnitude in kcal/mol (PRIMER_END_STABILITY).</param>
public readonly record struct Primer3PenaltyInputs(
    double Tm,
    int Length,
    double GcPercent,
    double SelfAny = 0.0,
    double SelfEnd = 0.0,
    int NumNs = 0,
    double HairpinTh = 0.0,
    double EndStability = 0.0)
{
    /// <summary>
    /// Library mispriming score (Primer3 <c>repeat_sim.score[repeat_sim.max]</c>, PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING;
    /// <see cref="PrimerDesigner.CalculateLibraryMispriming"/>); 0 without a mispriming library.
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// Template mispriming score (Primer3 <c>oligo_max_template_mispriming</c>, PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING[_TH];
    /// <see cref="PrimerDesigner.CalculateTemplateMispriming"/>); 0 when not computed.
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// Fraction bound in % (Primer3 <c>h->bound</c>, <see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/>), or
    /// <c>null</c> when the bound terms do not apply: for a left/right primer Primer3 adds the PRIMER_WT_BOUND_GT/_LT
    /// terms only when PRIMER_ANNEALING_TEMP &gt; 0 (pass null otherwise). The internal-oligo branch of <c>p_obj_fn</c>
    /// has no such gate: an internal oligo without a bound value carries Primer3's OLIGOTM_ERROR (−999999.9999).
    /// </summary>
    public double? Bound { get; init; }

    /// <summary>
    /// Position penalty (Primer3 <c>h->position_penalty</c>, <see cref="PrimerDesigner.CalculatePositionPenaltyPrimer3"/>);
    /// 0 under Primer3's default inside / outside penalties.
    /// </summary>
    public double PositionPenalty { get; init; }

    /// <summary>
    /// Minimum base quality of the oligo (Primer3 <c>h->seq_quality</c>, <see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>),
    /// or <c>null</c> without SEQUENCE_QUALITY (Primer3 then sets it to <see cref="QualityRangeMax"/>, so the
    /// PRIMER_WT_SEQ_QUAL term is 0).
    /// </summary>
    public int? SequenceQuality { get; init; }

    /// <summary>PRIMER_QUALITY_RANGE_MAX of the PRIMER_WT_SEQ_QUAL term (default 100).</summary>
    public int QualityRangeMax { get; init; } = PrimerDesigner.Primer3QualityRangeMax;

    /// <summary>Predicted PCR failure rate of the primer (Primer3 <c>h->failure_rate</c>,
    /// <see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>); 0 without template masking.</summary>
    public double MaskFailureRate { get; init; }
}

/// <summary>
/// Weights for the Primer3 per-primer penalty objective (the <c>PRIMER_WT_*</c>
/// parameters). Tm, GC and size each have separate "greater-than" (_gt) and
/// "less-than" (_lt) weights, applied one-sidedly relative to the optimum.
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Weights"/>.
/// </summary>
/// <param name="TmGt">PRIMER_WT_TM_GT.</param>
/// <param name="TmLt">PRIMER_WT_TM_LT.</param>
/// <param name="SizeGt">PRIMER_WT_SIZE_GT.</param>
/// <param name="SizeLt">PRIMER_WT_SIZE_LT.</param>
/// <param name="GcGt">PRIMER_WT_GC_PERCENT_GT.</param>
/// <param name="GcLt">PRIMER_WT_GC_PERCENT_LT.</param>
/// <param name="SelfAny">PRIMER_WT_SELF_ANY (alignment mode only).</param>
/// <param name="SelfEnd">PRIMER_WT_SELF_END (alignment mode only).</param>
/// <param name="NumNs">PRIMER_WT_NUM_NS.</param>
/// <param name="SelfAnyTh">PRIMER_WT_SELF_ANY_TH (thermodynamic mode only).</param>
/// <param name="SelfEndTh">PRIMER_WT_SELF_END_TH (thermodynamic mode only).</param>
/// <param name="HairpinTh">PRIMER_WT_HAIRPIN_TH (thermodynamic mode only).</param>
/// <param name="EndStability">PRIMER_WT_END_STABILITY.</param>
/// <param name="ThermodynamicOligoAlignment">PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: false = alignment-score
/// secondary-structure terms (SelfAny/SelfEnd weights), true = thermodynamic terms (…Th weights). Defaults to
/// false for an explicitly constructed weight set; <see cref="PrimerDesigner.DefaultPrimer3Weights"/> uses
/// Primer3's default (true).</param>
public readonly record struct Primer3PenaltyWeights(
    double TmGt,
    double TmLt,
    double SizeGt,
    double SizeLt,
    double GcGt,
    double GcLt,
    double SelfAny,
    double SelfEnd,
    double NumNs,
    double SelfAnyTh = 0.0,
    double SelfEndTh = 0.0,
    double HairpinTh = 0.0,
    double EndStability = 0.0,
    bool ThermodynamicOligoAlignment = false)
{
    /// <summary>
    /// PRIMER_WT_LIBRARY_MISPRIMING (Primer3 <c>weights.repeat_sim</c>, default 0): × the library mispriming score
    /// (<see cref="Primer3PenaltyInputs.LibraryMispriming"/>). Non-zero requires
    /// <see cref="PrimerParameters.MisprimingLibrary"/> in <see cref="PrimerDesigner.EvaluatePrimer"/> /
    /// <see cref="PrimerDesigner.DesignPrimers"/> (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_WT_TEMPLATE_MISPRIMING (Primer3 <c>weights.template_mispriming</c>, default 0): × the template mispriming
    /// score (<see cref="Primer3PenaltyInputs.TemplateMispriming"/>) when <see cref="ThermodynamicTemplateAlignment"/> is false.
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_WT_TEMPLATE_MISPRIMING_TH (<c>weights.template_mispriming_th</c>, default 0): thermodynamic template term
    /// (when <see cref="ThermodynamicTemplateAlignment"/> is true) with the 5 °C <c>temp_cutoff</c> rule of the other
    /// thermodynamic terms: s ≥ Tm − 5 → w·(s − (Tm − 6)), else w / (Tm − 4 − s).
    /// </summary>
    public double TemplateMisprimingTh { get; init; }

    /// <summary>
    /// PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT for the template terms (default false = 0, Primer3's default);
    /// <see cref="PrimerDesigner.EvaluatePrimer"/> / <see cref="PrimerDesigner.DesignPrimers"/> take it from
    /// <see cref="PrimerParameters.ThermodynamicTemplateAlignment"/>.
    /// </summary>
    public bool ThermodynamicTemplateAlignment { get; init; }

    /// <summary>PRIMER_WT_BOUND_GT (Primer3 <c>weights.bound_gt</c>, default 0): × (bound − <see cref="Primer3Optima.OptBound"/>)
    /// when the fraction bound (<see cref="Primer3PenaltyInputs.Bound"/>) is above the optimum.</summary>
    public double BoundGt { get; init; }

    /// <summary>PRIMER_WT_BOUND_LT (<c>weights.bound_lt</c>, default 0): × (<see cref="Primer3Optima.OptBound"/> − bound)
    /// when the fraction bound is below the optimum.</summary>
    public double BoundLt { get; init; }

    /// <summary>PRIMER_WT_POS_PENALTY (<c>weights.pos_penalty</c>, Primer3 default 1): × the position penalty
    /// (<see cref="Primer3PenaltyInputs.PositionPenalty"/>).</summary>
    public double PositionPenalty { get; init; } = PrimerDesigner.Primer3WeightPositionPenalty;

    /// <summary>
    /// PRIMER_WT_SEQ_QUAL / PRIMER_INTERNAL_WT_SEQ_QUAL (<c>weights.seq_quality</c>, default 0): × (PRIMER_QUALITY_RANGE_MAX −
    /// the oligo's minimum base quality) (<see cref="Primer3PenaltyInputs.SequenceQuality"/>,
    /// <see cref="Primer3PenaltyInputs.QualityRangeMax"/>). Non-zero requires SEQUENCE_QUALITY (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double SequenceQuality { get; init; }

    /// <summary>
    /// PRIMER_WT_END_QUAL / PRIMER_INTERNAL_WT_END_QUAL (<c>weights.end_quality</c>, default 0). Primer3 2.3.1 parses and stores
    /// this weight but its <c>p_obj_fn</c> never reads it, so it has <b>no effect</b> on the penalty — reproduced: it is
    /// accepted and ignored.
    /// </summary>
    public double EndQuality { get; init; }

    /// <summary>
    /// PRIMER_WT_MASK_FAILURE_RATE (<c>weights.failure_rate</c>, default 0): × the primer's predicted failure rate
    /// (<see cref="Primer3PenaltyInputs.MaskFailureRate"/>, <see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>), which
    /// Primer3 computes only with PRIMER_MASK_TEMPLATE (<see cref="PrimerPairOptions.MaskTemplate"/>; otherwise 0). Left/right
    /// primers only; added after the size terms as in <c>p_obj_fn</c>.
    /// </summary>
    public double MaskFailureRate { get; init; }
}

/// <summary>
/// Optimal parameter values for the Primer3 penalty objective
/// (<c>PRIMER_OPT_TM</c>, <c>PRIMER_OPT_SIZE</c>, <c>PRIMER_OPT_GC_PERCENT</c>).
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Optima"/>.
/// </summary>
public readonly record struct Primer3Optima(
    double OptTm,
    int OptSize,
    double? OptGcPercent)
{
    /// <summary>PRIMER_OPT_BOUND / PRIMER_INTERNAL_OPT_BOUND (% bound; Primer3 default
    /// <see cref="PrimerDesigner.Primer3OptBound"/> = 97): the optimum of the bound terms.</summary>
    public double OptBound { get; init; } = PrimerDesigner.Primer3OptBound;
}

/// <summary>
/// Result of primer pair design. For a valid pair the optional values are Primer3's
/// PRIMER_PAIR_k_PENALTY (<see cref="PairPenalty"/>), PRIMER_PAIR_k_PRODUCT_TM (<see cref="ProductTm"/>),
/// PRIMER_PAIR_k_COMPL_ANY_TH / _COMPL_END_TH (<see cref="ComplAnyTh"/>/<see cref="ComplEndTh"/>, thermodynamic
/// screen only) and PRIMER_INTERNAL_k_* (<see cref="InternalOligo"/>, when an internal oligo was requested).
/// </summary>
/// <param name="Forward">Forward (left) primer.</param>
/// <param name="Reverse">Reverse (right) primer.</param>
/// <param name="IsValid">True when the pair satisfies every constraint.</param>
/// <param name="Message">Human-readable status.</param>
/// <param name="ProductSize">Product size (bp).</param>
/// <param name="PairPenalty">Primer3 pair penalty (<c>obj_fn</c>) of a valid pair.</param>
/// <param name="ProductTm">Primer3 product Tm (<c>long_seq_tm</c>, °C) of a valid pair.</param>
/// <param name="ComplAnyTh">Pair hetero-dimer Tm (ntthal ANY, °C).</param>
/// <param name="ComplEndTh">Pair 3′ hetero-dimer Tm (ntthal END1/END2, °C).</param>
/// <param name="InternalOligo">Internal hybridization oligo picked for the pair (PRIMER_PICK_INTERNAL_OLIGO).</param>
public sealed record PrimerPairResult(
    PrimerCandidate? Forward,
    PrimerCandidate? Reverse,
    bool IsValid,
    string Message,
    int ProductSize,
    double? PairPenalty = null,
    double? ProductTm = null,
    double? ComplAnyTh = null,
    double? ComplEndTh = null,
    ProbeDesigner.Primer3Probe? InternalOligo = null)
{
    /// <summary>
    /// Primer3 alignment-mode PRIMER_PAIR_k_COMPL_ANY (<see cref="PrimerDesigner.CalculatePrimerDimerAnyComplementarity"/>)
    /// of a valid pair under <see cref="PrimerStructureScreen.Primer3Alignment"/>; otherwise <c>null</c>.
    /// </summary>
    public double? ComplAny { get; init; }

    /// <summary>
    /// Primer3 alignment-mode PRIMER_PAIR_k_COMPL_END (<see cref="PrimerDesigner.CalculatePrimerDimerEndComplementarity"/>)
    /// of a valid pair under <see cref="PrimerStructureScreen.Primer3Alignment"/>; otherwise <c>null</c>.
    /// </summary>
    public double? ComplEnd { get; init; }

    /// <summary>
    /// Primer3 PRIMER_PAIR_k_LIBRARY_MISPRIMING score of a valid pair (<c>pair_repeat_sim</c>: the maximum over library
    /// entries of the integer part of the left + right primer scores) when <see cref="PrimerParameters.MisprimingLibrary"/>
    /// is set; otherwise <c>null</c>.
    /// </summary>
    public double? LibraryMispriming { get; init; }

    /// <summary>The library entry named in PRIMER_PAIR_k_LIBRARY_MISPRIMING (with <see cref="LibraryMispriming"/>).</summary>
    public string? LibraryMisprimingName { get; init; }

    /// <summary>
    /// Primer3 PRIMER_PAIR_k_TEMPLATE_MISPRIMING (alignment mode) or _TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C) of a
    /// valid pair — max(left same-strand + right other-strand, left other-strand + right same-strand)
    /// (<see cref="TemplateMisprimingScore"/>) — when a pair template-mispriming limit or weight is set; otherwise <c>null</c>.
    /// </summary>
    public double? TemplateMispriming { get; init; }
}

/// <summary>A Primer3 product-size range (PRIMER_PRODUCT_SIZE_RANGE element), inclusive, in bp.</summary>
/// <param name="Min">Smallest product size.</param>
/// <param name="Max">Largest product size.</param>
public readonly record struct ProductSizeRange(int Min, int Max);

/// <summary>
/// Weights of Primer3's pair objective function <c>obj_fn</c> (<c>libprimer3.cc</c>), the
/// <c>PRIMER_PAIR_WT_*</c> tags; defaults are Primer3's (<c>pr_set_default_global_args_1</c>):
/// PRIMER_PAIR_WT_PR_PENALTY = 1, every other weight 0. The thermodynamic complementarity terms use
/// Primer3's fixed <c>temp_cutoff</c> = 5 °C relative to the lower primer Tm (<see cref="PrimerDesigner.Primer3TempCutoff"/>).
/// </summary>
/// <param name="PrimerPenalty">PRIMER_PAIR_WT_PR_PENALTY (× sum of the two primer penalties).</param>
/// <param name="InternalOligoPenalty">PRIMER_PAIR_WT_IO_PENALTY (× internal-oligo penalty).</param>
/// <param name="DiffTm">PRIMER_PAIR_WT_DIFF_TM (× |Tm_left − Tm_right|).</param>
/// <param name="ComplAnyTh">PRIMER_PAIR_WT_COMPL_ANY_TH (thermodynamic screen only).</param>
/// <param name="ComplEndTh">PRIMER_PAIR_WT_COMPL_END_TH (thermodynamic screen only).</param>
/// <param name="ProductTmLt">PRIMER_PAIR_WT_PRODUCT_TM_LT (needs <see cref="PrimerPairOptions.ProductOptTm"/>).</param>
/// <param name="ProductTmGt">PRIMER_PAIR_WT_PRODUCT_TM_GT (needs <see cref="PrimerPairOptions.ProductOptTm"/>).</param>
/// <param name="ProductSizeLt">PRIMER_PAIR_WT_PRODUCT_SIZE_LT (needs <see cref="PrimerPairOptions.ProductOptSize"/>).</param>
/// <param name="ProductSizeGt">PRIMER_PAIR_WT_PRODUCT_SIZE_GT (needs <see cref="PrimerPairOptions.ProductOptSize"/>).</param>
public sealed record Primer3PairWeights(
    double PrimerPenalty = 1.0,
    double InternalOligoPenalty = 0.0,
    double DiffTm = 0.0,
    double ComplAnyTh = 0.0,
    double ComplEndTh = 0.0,
    double ProductTmLt = 0.0,
    double ProductTmGt = 0.0,
    double ProductSizeLt = 0.0,
    double ProductSizeGt = 0.0)
{
    /// <summary>PRIMER_PAIR_WT_COMPL_ANY (× PRIMER_PAIR_COMPL_ANY; <see cref="PrimerStructureScreen.Primer3Alignment"/> only; Primer3 default 0).</summary>
    public double ComplAny { get; init; }

    /// <summary>PRIMER_PAIR_WT_COMPL_END (× PRIMER_PAIR_COMPL_END; <see cref="PrimerStructureScreen.Primer3Alignment"/> only; Primer3 default 0).</summary>
    public double ComplEnd { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_LIBRARY_MISPRIMING (× PRIMER_PAIR_LIBRARY_MISPRIMING, <see cref="PrimerPairResult.LibraryMispriming"/>;
    /// Primer3 default 0). Non-zero requires <see cref="PrimerParameters.MisprimingLibrary"/>.
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_TEMPLATE_MISPRIMING (× PRIMER_PAIR_k_TEMPLATE_MISPRIMING, <see cref="PrimerPairResult.TemplateMispriming"/>;
    /// alignment mode, <see cref="PrimerParameters.ThermodynamicTemplateAlignment"/> false; Primer3 default 0).
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH (thermodynamic mode; Primer3 default 0): the 5 °C <c>temp_cutoff</c> rule
    /// relative to the lower primer Tm, as for <see cref="ComplAnyTh"/>.
    /// </summary>
    public double TemplateMisprimingTh { get; init; }
}

/// <summary>
/// Pair-level options of <see cref="PrimerDesigner.DesignPrimers"/> / <see cref="PrimerDesigner.DesignPrimerPairs"/>,
/// named after the Primer3 tags they implement. Defaults are Primer3's except
/// <see cref="MaxTmDifference"/> (library default 5 °C; Primer3 100 °C — use <see cref="Primer3Defaults"/>).
/// </summary>
public sealed record PrimerPairOptions
{
    /// <summary>Library defaults (Primer3 defaults with PRIMER_PAIR_MAX_DIFF_TM = 5 °C).</summary>
    public static PrimerPairOptions Default { get; } = new();

    /// <summary>Primer3's own defaults (PRIMER_PAIR_MAX_DIFF_TM = 100 °C).</summary>
    public static PrimerPairOptions Primer3Defaults { get; } = new() { MaxTmDifference = PrimerDesigner.Primer3MaxPairTmDifference };

    /// <summary>
    /// PRIMER_PRODUCT_SIZE_RANGE: product-size ranges in order of preference (default 100–300 bp). A later
    /// range is used only when no further pair fits the earlier ones (Primer3 <c>choose_pair_or_triple</c>).
    /// </summary>
    public IReadOnlyList<ProductSizeRange> ProductSizeRanges { get; init; } = [new ProductSizeRange(100, 300)];

    /// <summary>PRIMER_PAIR_MAX_DIFF_TM: maximum |Tm_left − Tm_right| (°C); library default <see cref="PrimerDesigner.MaxPairTmDifference"/> = 5.</summary>
    public double MaxTmDifference { get; init; } = PrimerDesigner.MaxPairTmDifference;

    /// <summary>PRIMER_NUM_RETURN: number of pairs returned by <see cref="PrimerDesigner.DesignPrimerPairs"/> (default 5, ≥ 1).</summary>
    public int NumReturn { get; init; } = 5;

    /// <summary>SEQUENCE_INCLUDED_REGION (start, length): primers and products must lie inside it (default: whole template).</summary>
    public (int Start, int Length)? IncludedRegion { get; init; }

    /// <summary>PRIMER_PRODUCT_OPT_SIZE (bp); required when a product-size weight is non-zero.</summary>
    public int? ProductOptSize { get; init; }

    /// <summary>PRIMER_PRODUCT_OPT_TM (°C); required when a product-Tm weight is non-zero.</summary>
    public double? ProductOptTm { get; init; }

    /// <summary>PRIMER_PRODUCT_MIN_TM (°C): pairs whose product Tm is lower fail (default: no limit).</summary>
    public double? ProductMinTm { get; init; }

    /// <summary>PRIMER_PRODUCT_MAX_TM (°C): pairs whose product Tm is higher fail (default: no limit).</summary>
    public double? ProductMaxTm { get; init; }

    /// <summary>
    /// PRIMER_PAIR_MAX_COMPL_ANY: maximum Primer3 alignment-mode pair <c>compl_any</c> under
    /// <see cref="PrimerStructureScreen.Primer3Alignment"/> (default 8.00; must be in [0, 32767]).
    /// </summary>
    public double MaxComplAny { get; init; } = PrimerDesigner.Primer3MaxPairComplAny;

    /// <summary>
    /// PRIMER_PAIR_MAX_COMPL_END: maximum Primer3 alignment-mode pair <c>compl_end</c> (left vs. right orientation)
    /// under <see cref="PrimerStructureScreen.Primer3Alignment"/> (default 3.00; must be in [0, 32767]).
    /// </summary>
    public double MaxComplEnd { get; init; } = PrimerDesigner.Primer3MaxPairComplEnd;

    /// <summary>PRIMER_PAIR_WT_* pair objective weights (default Primer3's).</summary>
    public Primer3PairWeights Weights { get; init; } = new();

    /// <summary>
    /// PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (default 24.00): with a <see cref="PrimerParameters.MisprimingLibrary"/>, a pair
    /// fails when its library mispriming score (Primer3 <c>pair_repeat_sim</c>: the maximum over library entries of the
    /// integer part of left score + right score, <see cref="PrimerPairResult.LibraryMispriming"/>) is strictly greater.
    /// Under <see cref="PrimerStructureScreen.Primer3Alignment"/> it must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double MaxLibraryMispriming { get; init; } = PrimerDesigner.Primer3PairMaxLibraryMispriming;

    /// <summary>
    /// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING (alignment mode; default −100 = not checked): a pair whose template mispriming
    /// score (<see cref="PrimerPairResult.TemplateMispriming"/>) is strictly greater fails; a negative value is not checked.
    /// Must not exceed 32767 (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public double MaxTemplateMispriming { get; init; } = PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>
    /// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C; default −100). As in Primer3's
    /// <c>characterize_pair</c> the pair score is computed when this is ≥ 0 or PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH &gt; 0,
    /// and a pair fails when the limit is non-zero and the score exceeds it — so with the default −100 and a non-zero
    /// pair weight every pair fails, and 0 means no limit (Primer3 behaviour, reproduced).
    /// </summary>
    public double MaxTemplateMisprimingTh { get; init; } = PrimerDesigner.Primer3UndefinedTemplateMispriming;

    /// <summary>
    /// PRIMER_INSIDE_PENALTY (default <see cref="PrimerDesigner.Primer3DefaultInsidePenalty"/> = −1). With the default
    /// inside (−1) and outside (0) penalties primers never overlap the target. Any other combination makes Primer3 score
    /// every primer's 3′ end against the (single) target instead (<see cref="PrimerDesigner.CalculatePositionPenaltyPrimer3"/>):
    /// a primer may then extend into the target as long as its 3′ end does not pass the target's far end, a 3′ end
    /// inside the target costs <c>InsidePenalty</c> per base, one outside it <see cref="OutsidePenalty"/> per base of
    /// distance (× PRIMER_WT_POS_PENALTY, <see cref="Primer3PenaltyWeights.PositionPenalty"/>), and a pair must still
    /// span the target (left 3′ end before the right primer's 3′ end). Primer3 applies the value as given, so the default
    /// −1 makes inside positions negative when only <see cref="OutsidePenalty"/> is changed.
    /// </summary>
    public double InsidePenalty { get; init; } = PrimerDesigner.Primer3DefaultInsidePenalty;

    /// <summary>PRIMER_OUTSIDE_PENALTY (default 0): per-base penalty of a 3′ end outside the target when the position
    /// penalties are not Primer3's defaults (see <see cref="InsidePenalty"/>).</summary>
    public double OutsidePenalty { get; init; } = PrimerDesigner.Primer3DefaultOutsidePenalty;

    // Primer3 _PR_DEFAULT_POSITION_PENALTIES.
    internal bool DefaultPositionPenalties => PrimerDesigner.IsDefaultPositionPenalties(InsidePenalty, OutsidePenalty);

    /// <summary>
    /// PRIMER_PICK_INTERNAL_OLIGO: pick a hybridization (internal) oligo for every pair — the lowest-penalty
    /// oligo of the Primer3 internal-oligo list (<see cref="ProbeDesigner.DesignProbesPrimer3"/> rules) lying
    /// strictly between the primers; a pair without one fails (Primer3 <c>choose_internal_oligo</c>).
    /// </summary>
    public bool PickInternalOligo { get; init; }

    /// <summary>PRIMER_INTERNAL_* settings of the internal oligo (default Primer3's).</summary>
    public ProbeDesigner.Primer3ProbeSettings? InternalOligo { get; init; }

    /// <summary>
    /// SEQUENCE_QUALITY: one integer quality per template base (null or empty = none, Primer3's default). With it every
    /// primer (and internal oligo) gets its minimum base quality (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>,
    /// <see cref="PrimerCandidate.MinSequenceQuality"/>), checked against <see cref="PrimerParameters.MinQuality"/> /
    /// <see cref="PrimerParameters.MinEndQuality"/> (<see cref="ProbeDesigner.Primer3ProbeSettings.MinQuality"/>) and
    /// weighted by PRIMER_WT_SEQ_QUAL (<see cref="Primer3PenaltyWeights.SequenceQuality"/>). Its length must equal the
    /// template length and every value must lie in [<see cref="PrimerParameters.QualityRangeMin"/>,
    /// <see cref="PrimerParameters.QualityRangeMax"/>] (Primer3 <c>_pr_data_control</c>).
    /// </summary>
    public IReadOnlyList<int>? SequenceQuality { get; init; }

    /// <summary>
    /// PRIMER_MASK_TEMPLATE: mask the included region with Primer3's k-mer masker (<see cref="PrimerDesigner.MaskTemplatePrimer3"/>,
    /// <see cref="MaskKmerLists"/>, <see cref="MaskFailureRate"/>, <see cref="MaskFivePrimeDirection"/>,
    /// <see cref="MaskThreePrimeDirection"/>): a left primer whose 3′ base is masked on the forward copy, or a right primer whose
    /// 3′ base (its leftmost template base) is masked on the reverse copy, is rejected (Primer3 <c>is_lowercase_masked</c>), and
    /// every primer gets its predicted failure rate (<see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>,
    /// <see cref="PrimerCandidate.MaskFailureRate"/>) for PRIMER_WT_MASK_FAILURE_RATE
    /// (<see cref="Primer3PenaltyWeights.MaskFailureRate"/>). Requires <see cref="MaskKmerLists"/> (primer3-py: "masking template
    /// chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified"). Internal oligos are not masked (Primer3 checks them on the
    /// unmasked, case-preserving template — <see cref="LowercaseMasking"/>, which masking implies).
    /// </summary>
    public bool MaskTemplate { get; init; }

    /// <summary>
    /// PRIMER_LOWERCASE_MASKING (default false): with a case-preserving template (the <c>string</c> overloads of
    /// <see cref="PrimerDesigner.DesignPrimers(string, int, int, PrimerParameters?, PrimerPairOptions?)"/> /
    /// <see cref="PrimerDesigner.DesignPrimerPairs(string, int, int, PrimerParameters?, PrimerPairOptions?)"/>) a left primer
    /// or internal oligo whose 3′ base (rightmost template base), or a right primer whose 3′ base (leftmost template base), is a
    /// lower-case a/c/g/t is rejected (Primer3 <c>calc_and_check_oligo_features</c> → <c>is_lowercase_masked</c> on
    /// <c>trimmed_orig_seq</c>); lower-case bases elsewhere in the oligo are accepted. <see cref="MaskTemplate"/> implies it
    /// (primer3-py sets <c>lowercase_masking = mask_template</c>; primers are then checked on the masked copies, which keep the
    /// input's lower case). No effect with a <see cref="DnaSequence"/> template (upper case).
    /// </summary>
    public bool LowercaseMasking { get; init; }

    /// <summary>The masker's 11-mer / 16-mer genome count lists (PRIMER_MASK_KMERLIST_PATH / _PREFIX; see
    /// <see cref="PrimerMaskingKmerLists.FromGenomeTester4Files"/>). Not part of the JSON form.</summary>
    [System.Text.Json.Serialization.JsonIgnore]
    public PrimerMaskingKmerLists? MaskKmerLists { get; init; }

    /// <summary>PRIMER_MASK_FAILURE_RATE (default 0.1): a template position is masked when the failure rate of a primer ending
    /// there exceeds it (0 = no masking by score).</summary>
    public double MaskFailureRate { get; init; } = PrimerDesigner.Primer3MaskFailureRate;

    /// <summary>PRIMER_MASK_5P_DIRECTION (default 1, ≥ 0): bases masked from a high-failure 3′ end towards the primer's 5′ end
    /// (including the 3′ base).</summary>
    public int MaskFivePrimeDirection { get; init; } = PrimerDesigner.Primer3MaskFivePrimeDirection;

    /// <summary>PRIMER_MASK_3P_DIRECTION (default 0, 0–4984): bases masked beyond a high-failure 3′ end.</summary>
    public int MaskThreePrimeDirection { get; init; } = PrimerDesigner.Primer3MaskThreePrimeDirection;

    /// <summary>
    /// PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE (≥ −1, default −1): once <see cref="PrimerDesigner.DesignPrimerPairs"/> has
    /// selected a pair, no later pair may use a left primer whose 3′ end is fewer than this many bases from the selected
    /// left primer's 3′ end; 0 excludes only the identical left primer; −1 allows reuse (Primer3
    /// <c>choose_pair_or_triple</c>, <c>left_oligo_in_pair_overlaps_used_oligo</c>).
    /// </summary>
    public int MinLeftThreePrimeDistance { get; init; } = -1;

    /// <summary>
    /// PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE (≥ −1, default −1): as <see cref="MinLeftThreePrimeDistance"/> for right
    /// primers (3′ end = the leftmost top-strand base, <see cref="PrimerCandidate.Position"/>).
    /// </summary>
    public int MinRightThreePrimeDistance { get; init; } = -1;

    /// <summary>
    /// PRIMER_MIN_THREE_PRIME_DISTANCE: Primer3's shorthand that sets <see cref="MinLeftThreePrimeDistance"/> and
    /// <see cref="MinRightThreePrimeDistance"/> to the same value (and PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE, which
    /// Primer3 applies only with SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST — not modelled here, so it has no effect).
    /// Reading it returns the common value, or null when the left and right distances differ.
    /// </summary>
    public int? MinThreePrimeDistance
    {
        get => MinLeftThreePrimeDistance == MinRightThreePrimeDistance ? MinLeftThreePrimeDistance : null;
        init
        {
            MinLeftThreePrimeDistance = value ?? -1;
            MinRightThreePrimeDistance = value ?? -1;
        }
    }
}
