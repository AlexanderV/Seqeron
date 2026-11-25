namespace Seqeron.Genomics.Annotation;

public static partial class MiRnaAnalyzer
{
    #region Grimson et al. (2007) context score — TargetScan 5.0 (MIRNA-TARGET-001 / WP-6)

    /// <summary>
    /// The fitted TargetScan 5.0 context score of one 8mer / 7mer-m8 / 7mer-A1 site and its four
    /// contributions (Grimson et al. 2007, Mol Cell 27:91; <c>targetscan_50_context_scores.pl</c>).
    /// All values are log2 fold-change units: more negative = stronger predicted repression.
    /// </summary>
    /// <param name="Type">The site type scored (8mer, 7mer-m8 or 7mer-A1).</param>
    /// <param name="SiteTypeContribution">−0.310 (8mer), −0.161 (7mer-m8), −0.099 (7mer-A1).</param>
    /// <param name="ThreePrimePairingContribution">slope × 3' pairing score + intercept (not rounded, as the perl prints it).</param>
    /// <param name="LocalAuContribution">slope × local-AU fraction + intercept, rounded to 3 decimals (<c>sprintf("%.3f")</c>).</param>
    /// <param name="PositionContribution">slope × min(distance to the nearest 3'UTR end, 1500) + intercept, rounded to 3 decimals.</param>
    /// <param name="ThreePrimePairingScore">Raw 3' pairing score (1 point per pair at offset-adjusted positions 4–7, ½ elsewhere, minus ½ per nt of offset beyond 2).</param>
    /// <param name="LocalAuFraction">Position-weighted A/U fraction of the 30 nt flanking the site (unrounded).</param>
    /// <param name="ContextScore">Sum of the four contributions; <c>null</c> ("too_close") when the 1-based site start is &lt; 15.</param>
    public readonly record struct GrimsonContextScore(
        TargetSiteType Type,
        double SiteTypeContribution,
        double ThreePrimePairingContribution,
        double LocalAuContribution,
        double PositionContribution,
        double ThreePrimePairingScore,
        double LocalAuFraction,
        double? ContextScore);

    // Regression constants of targetscan_50_context_scores.pl ("Regression info from Grimson et al., 2007
    // supp data"), indexed by the perl site type 1 = 7mer-A1 ("1a"), 2 = 7mer-m8, 3 = 8mer.
    private static readonly double[] GrimsonSiteTypeContribution = { 0, -0.099, -0.161, -0.310 };
    private static readonly double[] GrimsonAuSlope = { 0, -0.42, -0.5, -0.64 };
    private static readonly double[] GrimsonAuIntercept = { 0, 0.236, 0.269, 0.365 };
    private static readonly double[] GrimsonPositionSlope = { 0, 0.000072, 0.000091, 0.000172 };
    private static readonly double[] GrimsonPositionIntercept = { 0, -0.032, -0.037, -0.07 };
    private static readonly double[] GrimsonThreePrimeSlope = { 0, -0.0211, -0.031, -0.0041 };
    private static readonly double[] GrimsonThreePrimeIntercept = { 0, 0.046, 0.067, 0.011 };
    private const int GrimsonMaxDistToUtrEnd = 1500;   // $maxDistToNearestEndOfUTREnd
    private const int GrimsonAuFlank = 30;             // substr(..., 30)

    // Perl site type 1..3 for the three site classes TargetScan 5.0 scores; 0 for every other class.
    private static int GrimsonPerlType(TargetSiteType type) => type switch
    {
        TargetSiteType.Seed7merA1 => 1,
        TargetSiteType.Seed7merM8 => 2,
        TargetSiteType.Seed8mer => 3,
        _ => 0
    };

    /// <summary>
    /// TargetScan 5.0 / Grimson et al. (2007) context score of a seed site: site type + 3' pairing +
    /// local AU + position, each with the published regression constants. <paramref name="mRnaSequence"/>
    /// is the 3'UTR (index 0 = first nt after the stop codon); <paramref name="siteStart"/> /
    /// <paramref name="siteEnd"/> are the 0-based inclusive coordinates of a <see cref="TargetSite"/>
    /// (an 8mer includes its A1).
    /// </summary>
    /// <returns>
    /// <c>null</c> for sites TargetScan 5.0 has no coefficients for (6mer, offset 6mer, supplementary,
    /// centered), for empty input or for coordinates outside the sequence. For a site with a 1-based
    /// start &lt; 15 the contributions are returned and <see cref="GrimsonContextScore.ContextScore"/>
    /// is <c>null</c> (the perl prints <c>too_close</c>).
    /// </returns>
    /// <remarks>
    /// Faithful port of <c>getLocalAUcontribution</c>, <c>getPositionContribution</c> and
    /// <c>get3primePairingContribution</c> of TargetScan 5.0 (the 3' pairing alignment is shared with the
    /// context++ port). Local AU: 30 nt up- and downstream, the site excluded, weights 1/(i+1) from the
    /// site edge outward (1/(i+2) on the upstream flank of a 7mer-A1 — whose window also skips the
    /// nt adjacent to the site — and on the downstream flank of an 8mer / 7mer-A1; the first two nt
    /// downstream of a 7mer-m8 weigh ½ each). TargetScan 5.0 coordinates exclude the A1 of an 8mer, so the
    /// A1 enters the 8mer's downstream window; this port keeps that (it is part of the fit). Near the 5'
    /// end TargetScan 5.0 passes a negative <c>substr</c> offset (a defect TargetScan fixed in 5.3 by
    /// clipping); the window is clipped at the UTR start here.
    /// </remarks>
    public static GrimsonContextScore? ScoreTargetSiteGrimson(
        string mRnaSequence, string miRnaSequence, int siteStart, int siteEnd, TargetSiteType type)
    {
        int perlType = GrimsonPerlType(type);
        if (perlType == 0 || string.IsNullOrEmpty(mRnaSequence) || string.IsNullOrEmpty(miRnaSequence)
            || siteStart < 0 || siteEnd >= mRnaSequence.Length || siteStart > siteEnd)
            return null;

        string mrna = NormalizeSeed(mRnaSequence);
        string mirna = NormalizeSeed(miRnaSequence);

        double fraction = GrimsonLocalAuFraction(mrna, siteStart, siteEnd, perlType, GrimsonAuFlank);
        double au = RoundPerl3(fraction * GrimsonAuSlope[perlType] + GrimsonAuIntercept[perlType]);
        double position = GrimsonPositionContribution(mrna.Length, siteStart, siteEnd, perlType);
        double threePrimeRaw = ThreePrimePairingRaw(mrna, mirna, siteStart, siteEnd, type);
        double threePrime = threePrimeRaw * GrimsonThreePrimeSlope[perlType] + GrimsonThreePrimeIntercept[perlType];
        double siteTypeTerm = GrimsonSiteTypeContribution[perlType];

        double? total = siteStart + 1 < MinDistToCds
            ? null
            : siteTypeTerm + au + position + threePrime;

        return new GrimsonContextScore(type, siteTypeTerm, threePrime, au, position, threePrimeRaw, fraction, total);
    }

    /// <summary>
    /// <see cref="ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/> for a site returned by
    /// <see cref="FindTargetSites"/>.
    /// </summary>
    public static GrimsonContextScore? ScoreTargetSiteGrimson(string mRnaSequence, MiRna miRna, TargetSite site) =>
        ScoreTargetSiteGrimson(mRnaSequence, miRna.Sequence, site.Start, site.End, site.Type);

    // getPositionContribution (5.0 coordinates: the A1 of an 8mer is not part of [utrStart, utrEnd]).
    private static double GrimsonPositionContribution(int utrLength, int siteStart, int siteEnd, int perlType)
    {
        int utrStart = siteStart + 1;
        int utrEnd = siteEnd + 1 - (perlType == 3 ? 1 : 0);
        int dist = Math.Min(utrStart - 1, utrLength - utrEnd);
        if (dist > GrimsonMaxDistToUtrEnd) dist = GrimsonMaxDistToUtrEnd;
        return RoundPerl3(dist * GrimsonPositionSlope[perlType] + GrimsonPositionIntercept[perlType]);
    }

    // getLocalAUcontribution: Σ weight(A/U) / Σ weight, over up to `flank` nt on each side of the site.
    // Upstream: the flank nt ending at utrStart−1 (utrStart−2 for 7mer-A1: perl substr(utrStart−32, 30)),
    // walked outward from the site; downstream: from 1-based utrEnd+1 (5.0 coordinates).
    private static double GrimsonLocalAuFraction(string mrna, int siteStart, int siteEnd, int perlType, int flank)
    {
        int utrStart = siteStart + 1;
        int utrEnd = siteEnd + 1 - (perlType == 3 ? 1 : 0);

        double sum = 0.0, max = 0.0;

        int upEndExclusive = utrStart - 1 - (perlType == 1 ? 1 : 0);   // 0-based exclusive end of the upstream window
        int upOffset = perlType == 1 ? 2 : 1;
        for (int i = 0; i < flank; i++)
        {
            int idx = upEndExclusive - 1 - i;
            if (idx < 0) break;
            double w = 1.0 / (i + upOffset);
            if (mrna[idx] is 'A' or 'U') sum += w;
            max += w;
        }

        for (int i = 0; i < flank; i++)
        {
            int idx = utrEnd + i;
            if (idx >= mrna.Length) break;
            double w = perlType switch
            {
                2 => i == 0 ? 0.5 : 1.0 / (i + 1),   // `if ($i == 0) { $scoreThisPos = 1 / 2; }`
                _ => 1.0 / (i + 2)
            };
            if (mrna[idx] is 'A' or 'U') sum += w;
            max += w;
        }

        return max == 0.0 ? 0.0 : sum / max;
    }

    #endregion

    #region Centered sites (Shin et al. 2010)

    /// <summary>
    /// Finds centered sites (Shin et al. 2010, Mol Cell 38:789): 3'UTR sites that lack perfect seed pairing
    /// (miRNA nt 2–7 are not all Watson–Crick paired) and instead pair the centre of the miRNA with at
    /// least 11 contiguous Watson–Crick pairs — nt 4–14, nt 5–15 or nt 4–15. G:U wobble does not count.
    /// </summary>
    /// <remarks>
    /// <see cref="TargetSite.Start"/>/<see cref="TargetSite.End"/> span the contiguous paired run (11 or 12 nt,
    /// <see cref="TargetSite.SeedMatchLength"/>); the duplex fields are built like those of
    /// <see cref="FindTargetSites"/>. TargetScan has no coefficients for centered sites, so
    /// <see cref="TargetSite.Score"/> is 0 and <see cref="TargetSite.ContextScore"/> is <c>null</c>.
    /// The additional Shin criterion "lacks 3'-compensatory pairing" is not tested.
    /// </remarks>
    public static IEnumerable<TargetSite> FindCenteredSites(string mRnaSequence, MiRna miRna)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || string.IsNullOrEmpty(miRna.Sequence))
            yield break;

        string mrna = NormalizeSeed(mRnaSequence);
        string mirna = NormalizeSeed(miRna.Sequence);
        const int firstNt = 4, lastNt = 15;
        if (mirna.Length < lastNt) yield break;

        // miRNA nt k (1-based) pairs mrna[i + lastNt − k]: window of 12 nt i..i+11 for nt 15..4.
        for (int i = 0; i + (lastNt - firstNt) < mrna.Length; i++)
        {
            bool Paired(int nt) => IsWatsonCrick(mirna[nt - 1], mrna[i + lastNt - nt]);

            bool seed = true;
            for (int nt = 2; nt <= 7 && seed; nt++)
                seed = i + lastNt - nt < mrna.Length && Paired(nt);
            if (seed) continue;   // a canonical seed site, not a centered one

            bool run4to14 = true, run5to15 = true;
            for (int nt = 4; nt <= 14; nt++) run4to14 &= Paired(nt);
            for (int nt = 5; nt <= 15; nt++) run5to15 &= Paired(nt);
            if (!run4to14 && !run5to15) continue;

            int ntFrom = run4to14 ? 4 : 5;
            int ntTo = run5to15 ? 15 : 14;
            int start = i + lastNt - ntTo;
            int length = ntTo - ntFrom + 1;
            yield return CreateTargetSite(mrna, start, length, i + lastNt - 1, mirna, miRna.Name, TargetSiteType.Centered, length);
        }
    }

    private static bool IsWatsonCrick(char miRnaBase, char targetBase) =>
        (miRnaBase, targetBase) is ('A', 'U') or ('U', 'A') or ('G', 'C') or ('C', 'G');

    #endregion
}
