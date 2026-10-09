namespace Seqeron.Genomics.MolTools;

// Primer3 sequence quality (SEQUENCE_QUALITY, PRIMER_MIN_QUALITY, PRIMER_MIN_END_QUALITY, PRIMER_QUALITY_RANGE_MIN/MAX,
// PRIMER_WT_SEQ_QUAL / PRIMER_WT_END_QUAL and the PRIMER_INTERNAL_* counterparts) — primer3-py 2.3.1 libprimer3.c
// sequence_quality_is_ok / calc_and_check_oligo_features / p_obj_fn / _pr_data_control.
public static partial class PrimerDesigner
{
    /// <summary>Primer3's default PRIMER_QUALITY_RANGE_MIN (0, <c>pr_set_default_global_args_1</c>).</summary>
    public const int Primer3QualityRangeMin = 0;

    /// <summary>Primer3's default PRIMER_QUALITY_RANGE_MAX (100, <c>pr_set_default_global_args_1</c>).</summary>
    public const int Primer3QualityRangeMax = 100;

    /// <summary>
    /// Primer3's sequence quality of an oligo (<c>libprimer3.c</c> <c>sequence_quality_is_ok</c>): the minimum of
    /// <paramref name="qualityRangeMax"/> and the base qualities (SEQUENCE_QUALITY, indexed by template position) over the
    /// whole oligo (<c>h->seq_quality</c>, primer3-py PRIMER_LEFT/RIGHT/INTERNAL_n_MIN_SEQ_QUALITY), and over its five
    /// 3′-most bases (<c>h->seq_end_quality</c>): the last five template positions of a left primer or internal oligo,
    /// the first five of a right primer (its 3′ end is the leftmost template base); a shorter oligo uses all its bases.
    /// </summary>
    /// <param name="quality">Per-base qualities of the whole template (SEQUENCE_QUALITY).</param>
    /// <param name="position">0-based leftmost template coordinate of the oligo (as <see cref="PrimerCandidate.Position"/>).</param>
    /// <param name="length">Oligo length (≥ 1).</param>
    /// <param name="isForward">True for a left primer or an internal oligo, false for a right primer.</param>
    /// <param name="qualityRangeMax">PRIMER_QUALITY_RANGE_MAX (default 100): the starting value of both minima.</param>
    /// <returns>(<c>Min</c> = seq_quality, <c>MinEnd</c> = seq_end_quality).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="quality"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The oligo lies outside <paramref name="quality"/>.</exception>
    public static (int Min, int MinEnd) CalculateSequenceQualityPrimer3(
        IReadOnlyList<int> quality, int position, int length, bool isForward, int qualityRangeMax = Primer3QualityRangeMax)
    {
        ArgumentNullException.ThrowIfNull(quality);
        if (length < 1 || position < 0 || position > quality.Count - length) // overflow-safe (heavy-tier guard fuzz, F74)
            throw new ArgumentOutOfRangeException(nameof(position), "The oligo must lie inside the quality array.");
        int j = position, k = position + length - 1;
        // The 3' window first (left / internal: k-4..k; right: j..j+4), then the rest with the running minimum.
        int endLo = isForward ? Math.Max(j, k - 4) : j;
        int endHi = isForward ? k : Math.Min(k, j + 4);
        int q = qualityRangeMax;
        for (int i = endLo; i <= endHi; i++)
            q = Math.Min(q, quality[i]);
        int minEnd = q;
        for (int i = j; i <= k; i++)
            q = Math.Min(q, quality[i]);
        return (q, minEnd);
    }

    // Primer3 _pr_data_control checks of the sequence-quality settings (primer p_args and internal-oligo o_args), in
    // Primer3's order: the quality length, quality data required by a non-zero PRIMER_[INTERNAL_]MIN_QUALITY, the minimum
    // qualities inside the quality range, every quality inside the range, and a quality weight without quality data.
    internal static void ValidatePrimer3Quality(IReadOnlyList<int>? quality, int templateLength,
        int primerMinQuality, int internalMinQuality, int rangeMin, int rangeMax,
        double primerWeight, double internalWeight, string paramName)
    {
        bool has = quality is { Count: > 0 };
        if (has && quality!.Count != templateLength)
            throw new ArgumentException("Error in sequence quality data (Primer3 _pr_data_control: SEQUENCE_QUALITY length ≠ template length).", paramName);
        if ((primerMinQuality != 0 || internalMinQuality != 0) && !has)
            throw new ArgumentException("Sequence quality data missing (Primer3 _pr_data_control).", paramName);
        if (double.IsNaN(primerWeight) || double.IsNaN(internalWeight))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_[INTERNAL_]WT_SEQ_QUAL must not be NaN.");
        if (has)
        {
            if (primerMinQuality != 0 && primerMinQuality < rangeMin)
                throw new ArgumentException("PRIMER_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN (Primer3 _pr_data_control).", paramName);
            if (primerMinQuality != 0 && primerMinQuality > rangeMax)
                throw new ArgumentException("PRIMER_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX (Primer3 _pr_data_control).", paramName);
            if (internalMinQuality != 0 && internalMinQuality < rangeMin)
                throw new ArgumentException("PRIMER_INTERNAL_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN (Primer3 _pr_data_control).", paramName);
            if (internalMinQuality != 0 && internalMinQuality > rangeMax)
                throw new ArgumentException("PRIMER_INTERNAL_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX (Primer3 _pr_data_control).", paramName);
            foreach (int v in quality!)
                if (v < rangeMin || v > rangeMax)
                    throw new ArgumentException("Sequence quality score out of range (Primer3 _pr_data_control).", paramName);
        }
        else if (primerWeight != 0 || internalWeight != 0)
            throw new ArgumentException(
                "Sequence quality is part of objective function but sequence quality is not defined (Primer3 _pr_data_control).", paramName);
    }

    // Per-template quality data of a design (null when SEQUENCE_QUALITY is absent).
    internal sealed record QualityContext(IReadOnlyList<int> Values, int RangeMax)
    {
        internal static QualityContext? Create(IReadOnlyList<int>? quality, int rangeMax) =>
            quality is { Count: > 0 } q ? new QualityContext(q, rangeMax) : null;
    }
}
