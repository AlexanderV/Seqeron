using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Bucher (1990) weight-matrix promoter elements: the four RNA polymerase II promoter-element count
/// matrices of Bucher P., J Mol Biol 212:563–578 (502 unrelated promoters, EPD), as distributed in the
/// JASPAR POLII collection (MEDLINE 2329577): POL012.1 TATA-Box, POL002.1 INR (Bucher's cap signal),
/// POL004.1 CCAAT-box, POL003.1 GC-box.
/// </summary>
public static partial class MotifFinder
{
    /// <summary>
    /// The Bucher (1990) promoter-element count matrices (JASPAR POLII; identical in JASPAR 2014, 2016,
    /// 2018 and 2020 — pyjaspar 4.0.0 SQLite releases). Counts are rows A, C, G, T.
    /// </summary>
    public static class BucherPromoterMatrices
    {
        /// <summary>TATA box, JASPAR POL012.1 (389 sites, 15 columns; TSS −39…−23 to −25…−9).</summary>
        public static PromoterElementMatrix TataBox { get; } = new(
            "TATA Box", "POL012.1", orientationIndependent: false, new double[,]
            {
                { 61, 16, 352, 3, 354, 268, 360, 222, 155, 56, 83, 82, 82, 68, 77 },
                { 145, 46, 0, 10, 0, 0, 3, 2, 44, 135, 147, 127, 118, 107, 101 },
                { 152, 18, 2, 2, 5, 0, 10, 44, 157, 150, 128, 128, 128, 139, 140 },
                { 31, 309, 35, 374, 30, 121, 6, 121, 33, 48, 31, 52, 61, 75, 71 },
            });

        /// <summary>Cap signal / initiator (INR), JASPAR POL002.1 (303 sites, 8 columns; TSS −2 to +6).</summary>
        public static PromoterElementMatrix CapSignal { get; } = new(
            "Cap Signal", "POL002.1", orientationIndependent: false, new double[,]
            {
                { 49, 0, 288, 26, 77, 67, 45, 50 },
                { 48, 303, 0, 81, 95, 118, 85, 96 },
                { 69, 0, 0, 116, 0, 46, 73, 56 },
                { 137, 0, 15, 80, 131, 72, 100, 101 },
            });

        /// <summary>CCAAT box, JASPAR POL004.1 (175 sites, 12 columns; orientation-independent — Mantovani 1998).</summary>
        public static PromoterElementMatrix CcaatBox { get; } = new(
            "CCAAT Box", "POL004.1", orientationIndependent: true, new double[,]
            {
                { 56, 32, 25, 102, 51, 0, 0, 175, 119, 17, 23, 116 },
                { 55, 52, 47, 1, 6, 173, 174, 0, 8, 0, 90, 6 },
                { 12, 43, 24, 70, 99, 1, 0, 0, 21, 15, 59, 52 },
                { 52, 48, 79, 2, 19, 1, 1, 0, 27, 143, 3, 1 },
            });

        /// <summary>GC box (Sp1), JASPAR POL003.1 (274 sites, 14 columns; orientation-independent — Gidoni et al. 1985).</summary>
        public static PromoterElementMatrix GcBox { get; } = new(
            "GC Box", "POL003.1", orientationIndependent: true, new double[,]
            {
                { 102, 97, 50, 67, 0, 2, 54, 46, 1, 79, 23, 0, 20, 40 },
                { 40, 31, 6, 1, 0, 0, 170, 1, 3, 0, 17, 166, 86, 24 },
                { 50, 112, 154, 206, 274, 272, 0, 224, 222, 171, 192, 35, 52, 109 },
                { 82, 34, 64, 0, 0, 0, 50, 3, 48, 24, 42, 73, 116, 101 },
            });

        /// <summary>All four matrices in Bucher's order: TATA box, cap signal, CCAAT box, GC box.</summary>
        public static IReadOnlyList<PromoterElementMatrix> All { get; } = new[] { TataBox, CapSignal, CcaatBox, GcBox };
    }

    /// <summary>
    /// Scans a sequence with the four <see cref="BucherPromoterMatrices"/>. Each matrix's
    /// <see cref="PromoterElementMatrix.Pwm"/> (JASPAR pseudocounts, uniform background) is scanned through
    /// the canonical PWM path at the score threshold whose background false-positive rate is
    /// <paramref name="falsePositiveRate"/> (<see cref="PwmScoreDistribution.ThresholdFpr"/>, Biopython
    /// <c>pssm.distribution().threshold_fpr(fpr)</c>, precision 10³).
    /// </summary>
    /// <remarks>
    /// With <paramref name="bothStrands"/> = true the orientation-independent CCAAT and GC boxes are scanned on
    /// both strands (<see cref="ScanWithPwmBothStrands"/>); TATA box and cap signal are strand-specific and
    /// scanned on the given strand only (<see cref="ScanWithPwm"/>). Bucher's own cut-off values (e.g. −8.16
    /// for the TATA box) are defined on his smoothed natural-log weight scale, not on these log2-odds scores,
    /// and are therefore not used. Order: matrix order, then ascending position, '+' before '-'.
    /// </remarks>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <param name="falsePositiveRate">Background false-positive rate per window and strand, in [0, 1].</param>
    /// <param name="bothStrands">Scan the orientation-independent matrices on both strands.</param>
    /// <returns>Promoter-element hits.</returns>
    public static IEnumerable<PromoterMatrixHit> FindPromoterElementsByMatrix(
        DnaSequence sequence, double falsePositiveRate, bool bothStrands = true)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (!(falsePositiveRate >= 0 && falsePositiveRate <= 1))
            throw new ArgumentOutOfRangeException(nameof(falsePositiveRate), falsePositiveRate, "Rate must be in [0, 1].");
        return FindPromoterElementsByMatrixCore(sequence, falsePositiveRate, bothStrands);
    }

    private static IEnumerable<PromoterMatrixHit> FindPromoterElementsByMatrixCore(
        DnaSequence sequence, double falsePositiveRate, bool bothStrands)
    {
        foreach (var element in BucherPromoterMatrices.All)
        {
            double threshold = element.Pwm.ScoreDistribution().ThresholdFpr(falsePositiveRate);
            IEnumerable<PwmStrandMatch> hits = bothStrands && element.OrientationIndependent
                ? ScanWithPwmBothStrands(sequence, element.Pwm, threshold)
                : ScanWithPwm(sequence, element.Pwm, threshold)
                    .Select(m => new PwmStrandMatch(m.Position, m.Position, '+', m.MatchedSequence, m.Pattern, m.Score));

            foreach (var hit in hits)
            {
                yield return new PromoterMatrixHit(element.Name, element.MatrixId, hit.Position, hit.Strand,
                    hit.MatchedSequence, hit.Score, threshold);
            }
        }
    }
}

/// <summary>
/// A published promoter-element count matrix with its source ID and the derived PWM.
/// </summary>
public sealed class PromoterElementMatrix
{
    private readonly double[,] _counts;

    internal PromoterElementMatrix(string name, string matrixId, bool orientationIndependent, double[,] counts)
    {
        Name = name;
        MatrixId = matrixId;
        OrientationIndependent = orientationIndependent;
        _counts = counts;
        Pwm = PositionWeightMatrix.FromCounts(counts, MotifFinder.JasparPseudocounts(counts));
    }

    /// <summary>Element name.</summary>
    public string Name { get; }

    /// <summary>Source matrix ID (JASPAR POLII collection).</summary>
    public string MatrixId { get; }

    /// <summary>Whether the element acts in either orientation (scanned on both strands when requested).</summary>
    public bool OrientationIndependent { get; }

    /// <summary>Returns a copy of the 4 × L count matrix (rows A, C, G, T).</summary>
    public double[,] GetCounts() => (double[,])_counts.Clone();

    /// <summary>
    /// Log2-odds PWM against a uniform background with JASPAR pseudocounts √N̄·0.25 — Biopython
    /// <c>m = motifs.read(f, "jaspar"); m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m); m.pssm</c>.
    /// </summary>
    public PositionWeightMatrix Pwm { get; }
}

/// <summary>A promoter-element PWM hit (<see cref="MotifFinder.FindPromoterElementsByMatrix"/>).</summary>
/// <param name="Name">Element name.</param>
/// <param name="MatrixId">Source matrix ID.</param>
/// <param name="Position">0-based forward-strand window start.</param>
/// <param name="Strand">'+' or '-'.</param>
/// <param name="Sequence">The site read 5'→3' on its own strand.</param>
/// <param name="Score">Log2-odds score.</param>
/// <param name="Threshold">Score threshold used (from the false-positive rate).</param>
public readonly record struct PromoterMatrixHit(
    string Name,
    string MatrixId,
    int Position,
    char Strand,
    string Sequence,
    double Score,
    double Threshold);
