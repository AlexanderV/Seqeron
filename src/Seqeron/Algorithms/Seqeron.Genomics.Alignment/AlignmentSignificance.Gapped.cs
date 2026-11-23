using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentSignificance
{
    private const int UngappedMarker = 32767;     // INT2_MAX of blast_stat.c: the existence / extension of the ungapped row

    // ------------------------------------------------------------------ the NCBI tables

    /// <summary>
    /// The tabulated parameters of a standard protein matrix with gap costs in the BLAST convention (a gap of length k costs <paramref name="gapExistence"/> + k x <paramref name="gapExtension"/>; both positive;
    /// existence and extension of 32767 = the ungapped row). Matrices: BLOSUM45 / 50 / 62 / 80 / 90, PAM30 / 70 / 250 (<c>blast_stat.c</c>).
    /// </summary>
    /// <param name="matrixName">The matrix name (case-insensitive; a leading E of an EMBOSS name is ignored).</param>
    /// <param name="gapExistence">BLAST gap existence cost.</param>
    /// <param name="gapExtension">BLAST gap extension cost.</param>
    /// <param name="parameters">The Karlin-Altschul parameters (alpha and beta included).</param>
    public static bool TryGappedProtein(string matrixName, int gapExistence, int gapExtension, out KarlinAltschulParameters parameters)
    {
        if (TryProteinRow(matrixName, gapExistence, gapExtension, out var row, out _))
        {
            parameters = new KarlinAltschulParameters(row[3], row[4], row[5], row[6], row[7]);
            return true;
        }

        parameters = null!;
        return false;
    }

    private static string CanonicalMatrixName(string name)
    {
        string n = name.Trim().ToUpperInvariant();
        if (!ProteinTables.ContainsKey(n) && n.StartsWith('E') && ProteinTables.ContainsKey(n[1..]))
            n = n[1..];
        return n;
    }

    private static bool TryProteinRow(string matrixName, int gapExistence, int gapExtension, out double[] row, out double[] ungappedRow)
    {
        row = null!;
        ungappedRow = null!;
        if (string.IsNullOrWhiteSpace(matrixName) || !ProteinTables.TryGetValue(CanonicalMatrixName(matrixName), out var rows))
            return false;
        ungappedRow = rows[0];
        foreach (var r in rows)
        {
            if (Math.Round(r[0]) == gapExistence && Math.Round(r[1]) == gapExtension)
            {
                row = r;
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// The tabulated parameters of a blastn reward / penalty pair with gap costs in the BLAST convention (<c>Blast_KarlinBlkNuclGappedCalc</c>: the pair is divided by its gcd, the costs of the table are multiplied back,
    /// lambda and alpha divided; gap costs of 0 / 0 = the linear-gap row; costs larger than the largest tabulated pair take the ungapped parameters, which are computed).
    /// Pairs: 1/-5, 1/-4, 2/-7, 1/-3, 2/-5, 1/-2, 2/-3, 3/-4, 4/-5, 1/-1, 3/-2, 5/-4 (and their multiples).
    /// </summary>
    /// <param name="reward">Match reward (positive).</param>
    /// <param name="penalty">Mismatch penalty (negative).</param>
    /// <param name="gapExistence">BLAST gap existence cost (0 or more).</param>
    /// <param name="gapExtension">BLAST gap extension cost.</param>
    /// <param name="parameters">The parameters.</param>
    /// <param name="roundDown">True when the table is only valid for even scores: an odd score is rounded down to the nearest even number before the E-value is computed.</param>
    public static bool TryGappedNucleotide(int reward, int penalty, int gapExistence, int gapExtension, out KarlinAltschulParameters parameters, out bool roundDown)
    {
        parameters = null!;
        roundDown = false;
        if (reward <= 0 || penalty >= 0)
            return false;
        int divisor = Gcd(reward, penalty);
        int r = reward / divisor, p = penalty / divisor;
        if (!NucleotideTables.TryGetValue((r, p), out var rows))
            return false;
        roundDown = (r, p) is (2, -7) or (2, -5) or (2, -3) or (3, -4);
        // the first row with existence = extension = 0 is the linear-gap row, the others are the affine rows
        bool hasLinear = rows[0][0] == 0 && rows[0][1] == 0;
        var affine = hasLinear ? rows.Skip(1).ToArray() : rows;
        int open = gapExistence / divisor, ext = gapExtension / divisor;
        if (gapExistence % divisor != 0 || gapExtension % divisor != 0)
            return false;
        double[]? hit = null;
        if (gapExistence == 0 && gapExtension == 0 && hasLinear)
        {
            hit = rows[0];
        }
        else
        {
            foreach (var row in affine)
            {
                if (row[0] == open && row[1] == ext)
                {
                    hit = row;
                    break;
                }
            }
        }

        if (hit is not null)
        {
            parameters = new KarlinAltschulParameters(hit[2] / divisor, hit[3], hit[4], hit[5] / divisor, hit[6]);
            return true;
        }

        // beyond the largest tabulated costs the ungapped parameters apply
        int maxOpen = affine.Max(x => (int)x[0]), maxExt = affine.Max(x => (int)x[1]);
        if (open >= maxOpen && ext >= maxExt)
        {
            var ungapped = UngappedNucleotide(reward, penalty);
            parameters = ungapped with { Beta = (r == 1 && p == -1) || (r == 2 && p == -3) ? -2 : 0 };
            return true;
        }

        return false;
    }

    // ------------------------------------------------------------------ E-values

    /// <summary>
    /// The edge-effect length adjustment of BLAST (<c>BLAST_ComputeLengthAdjustment</c>): an integer below the fixed point of ell = beta + (alpha / lambda) (ln K + ln((m - ell)(n - N ell))), subtracted from the
    /// query length and from every database sequence.
    /// </summary>
    /// <param name="parameters">The parameters (alpha and beta; for ungapped scoring 1 / H and 0).</param>
    /// <param name="queryLength">m.</param>
    /// <param name="databaseLength">n, the total length of the database (the subject for a pairwise comparison).</param>
    /// <param name="databaseSequences">N, the number of database sequences (1 for a pairwise comparison).</param>
    public static int LengthAdjustment(KarlinAltschulParameters parameters, long queryLength, long databaseLength, long databaseSequences = 1)
    {
        ArgumentNullException.ThrowIfNull(parameters);
        double alphaDLambda = double.IsNaN(parameters.Alpha) ? 1.0 / parameters.H : parameters.Alpha / parameters.Lambda;
        double beta = double.IsNaN(parameters.Beta) ? 0 : parameters.Beta;
        double logK = parameters.LogK;
        double m = queryLength, n = databaseLength, nn = databaseSequences;
        double a = nn, mb = (m * nn) + n, c = (n * m) - (Math.Max(m, n) / parameters.K);
        if (c < 0)
            return 0;
        double ellMax = 2 * c / (mb + Math.Sqrt((mb * mb) - (4 * a * c)));
        double ellMin = 0, ellNext = 0;
        bool converged = false;
        for (int i = 1; i <= 20; i++)
        {
            double ell = ellNext;
            double ss = (m - ell) * (n - (nn * ell));
            double ellBar = (alphaDLambda * (logK + Math.Log(ss))) + beta;
            if (ellBar >= ell)
            {
                ellMin = ell;
                if (ellBar - ellMin <= 1.0)
                {
                    converged = true;
                    break;
                }

                if (ellMin == ellMax)
                    break;
            }
            else
            {
                ellMax = ell;
            }

            if (ellMin <= ellBar && ellBar <= ellMax)
                ellNext = ellBar;
            else
                ellNext = i == 1 ? ellMax : (ellMin + ellMax) / 2;
        }

        int adjustment = (int)ellMin;
        if (converged)
        {
            double ceil = Math.Ceiling(ellMin);
            if (ceil <= ellMax)
            {
                double ss = (m - ceil) * (n - (nn * ceil));
                if ((alphaDLambda * (logK + Math.Log(ss))) + beta >= ceil)
                    adjustment = (int)ceil;
            }
        }

        return adjustment;
    }

    /// <summary>The E-value of the classic (length-adjusted) Karlin-Altschul statistics of BLAST: (m - ell)(n - N ell) K exp(-lambda S), blastn's formula.</summary>
    /// <param name="score">The raw score.</param>
    /// <param name="parameters">The parameters.</param>
    /// <param name="queryLength">Query length.</param>
    /// <param name="databaseLength">Database (subject) length.</param>
    /// <param name="databaseSequences">Number of database sequences.</param>
    public static AlignmentSignificanceResult KarlinAltschulSignificance(double score, KarlinAltschulParameters parameters, long queryLength, long databaseLength, long databaseSequences = 1)
    {
        int ell = LengthAdjustment(parameters, queryLength, databaseLength, databaseSequences);
        double space = Math.Max(1, queryLength - ell) * (double)Math.Max(1, databaseLength - (databaseSequences * ell));
        double e = parameters.EValue(score, space);
        return new AlignmentSignificanceResult(score, parameters.BitScore(score), e, KarlinAltschulParameters.PValue(e), "LengthAdjustedKarlinAltschul", ell, space, parameters);
    }

    /// <summary>
    /// The finite-size-corrected E-value of Spouge (2002; NCBI <c>BLAST_SpougeStoE</c> with the Gumbel block of the protein tables): what BLAST+ prints for a gapped protein alignment of a standard matrix.
    /// </summary>
    /// <param name="score">The raw score.</param>
    /// <param name="matrixName">A matrix of the tables.</param>
    /// <param name="gapExistence">BLAST gap existence cost.</param>
    /// <param name="gapExtension">BLAST gap extension cost.</param>
    /// <param name="queryLength">Query length m.</param>
    /// <param name="subjectLength">Subject length n.</param>
    /// <param name="databaseLength">The length of the database the subject belongs to (the E-value of a pair is scaled by databaseLength / subjectLength; 0 = the subject only, pairwise comparison).</param>
    /// <exception cref="ArgumentException">The matrix or the gap costs are not in the NCBI tables, or the row is the ungapped one (use <see cref="KarlinAltschulSignificance"/>).</exception>
    public static AlignmentSignificanceResult SpougeSignificance(double score, string matrixName, int gapExistence, int gapExtension, long queryLength, long subjectLength, long databaseLength = 0)
    {
        if (!TryProteinRow(matrixName, gapExistence, gapExtension, out var row, out var ungapped) || row[0] == UngappedMarker)
            throw new ArgumentException($"The gap costs {gapExistence} / {gapExtension} are not a gapped row of the NCBI table of {matrixName}.", nameof(gapExistence));
        double lambda = row[3], k = row[4];
        double g = gapExistence + gapExtension;
        double a = row[6], alphaV = row[9], sigma = row[10];
        double aUn = ungapped[6], alphaUn = ungapped[9];
        double b = 2.0 * g * (aUn - a);
        double beta = 2.0 * g * (alphaUn - alphaV);
        double tau = 2.0 * g * (alphaUn - sigma);
        const double constVal = 0.39894228040143267793994605993438;      // 1 / sqrt(2 pi)
        double dbScale = databaseLength > 0 ? (double)databaseLength / subjectLength : 1.0;

        double miY = queryLength - ((a * score) + b);
        double vY = Math.Max(2.0 * alphaV / lambda, (alphaV * score) + beta);
        double sqrtV = Math.Sqrt(vY);
        double mF = miY / sqrtV;
        double pMF = double.IsNaN(mF) ? 0 : Erfc(-mF / Math.Sqrt(2.0)) / 2.0;
        double p1 = (miY * pMF) + (sqrtV * constVal * Math.Exp(-0.5 * mF * mF));
        double njY = subjectLength - ((a * score) + b);
        double nF = njY / sqrtV;
        double pNF = Erfc(-nF / Math.Sqrt(2.0)) / 2.0;
        double p2 = (njY * pNF) + (sqrtV * constVal * Math.Exp(-0.5 * nF * nF));
        double cY = Math.Max(2.0 * sigma / lambda, (sigma * score) + tau);
        double area = (p1 * p2) + (cY * pMF * pNF);
        double e = area * k * Math.Exp(-lambda * score) * dbScale;
        var parameters = new KarlinAltschulParameters(lambda, k, row[5], row[6], row[7]);
        return new AlignmentSignificanceResult(score, parameters.BitScore(score), e, KarlinAltschulParameters.PValue(e), "SpougeFiniteSizeCorrection", 0, area, parameters);
    }

    // complementary error function to double precision: the Taylor series of erf below 2.5 (the cancellation costs at most 1e-9 of the relative value there), the continued fraction
    // erfc(x) = exp(-x^2) / sqrt(pi) / (x + (1/2) / (x + 1 / (x + (3/2) / (x + 2 / (x + ...))))) above
    internal static double Erfc(double x)
    {
        if (double.IsNaN(x))
            return double.NaN;
        if (x < 0)
            return 2.0 - Erfc(-x);
        if (x < 2.5)
        {
            double term = x, sum = x, x2 = x * x;
            for (int n = 1; n < 200; n++)
            {
                term *= -x2 / n;
                double add = term / ((2 * n) + 1);
                sum += add;
                if (Math.Abs(add) < 1e-18 * Math.Abs(sum))
                    break;
            }

            return 1.0 - (2.0 / Math.Sqrt(Math.PI) * sum);
        }

        if (x > 27)
            return 0.0;
        double tail = x;
        for (int k = 80; k >= 1; k--)
            tail = x + (k / 2.0 / tail);
        return Math.Exp(-x * x) / Math.Sqrt(Math.PI) / tail;
    }

    /// <summary>
    /// The significance of a local alignment score the way BLAST+ reports it: a protein matrix of the NCBI tables with its gap costs gets the Spouge finite-size E-value (<c>blastp -subject</c>), a blastn
    /// reward / penalty pair the length-adjusted Karlin-Altschul E-value (odd scores of the even-only tables are rounded down to an even number for the E-value, not for the bit score), the ungapped row the classic formula.
    /// </summary>
    /// <param name="score">The raw alignment score.</param>
    /// <param name="scoring">The scoring system.</param>
    /// <param name="queryLength">Query length.</param>
    /// <param name="subjectLength">Subject length.</param>
    /// <param name="databaseLength">Database length for the E-value (0 = the subject only).</param>
    /// <param name="databaseSequences">Number of database sequences (nucleotide, default 1).</param>
    /// <exception cref="ArgumentException">The scoring system has no tabulated gapped parameters.</exception>
    public static AlignmentSignificanceResult ForMatrix(double score, MatrixScoring scoring, long queryLength, long subjectLength, long databaseLength = 0, long databaseSequences = 1)
    {
        ArgumentNullException.ThrowIfNull(scoring);
        int existence = (int)Math.Round(-scoring.GapOpen), extension = (int)Math.Round(-scoring.GapExtend);
        // the library convention: a gap of length k scores GapOpen + k x GapExtend (negative); BLAST: existence + k x extension
        if (!TryProteinRow(scoring.Matrix.Name, existence, extension, out var row, out var ungappedRow))
            throw new ArgumentException($"No NCBI parameters for {scoring.Matrix.Name} with gap costs {existence} / {extension}.", nameof(scoring));
        if (row[0] == UngappedMarker)
            throw new ArgumentException("The ungapped row has no gap costs.", nameof(scoring));
        // the NCBI tables belong to the NCBI matrices: a matrix of the same name in other units (the EMBOSS BLOSUM80 has A/A = 7, the NCBI one 5) has another lambda and must not borrow their parameters
        double ungappedLambda = Ungapped(scoring.Matrix).Lambda;
        if (Math.Abs((ungappedLambda / ungappedRow[3]) - 1) > 0.02)
            throw new ArgumentException($"The matrix {scoring.Matrix.Name} is not the NCBI matrix of that name (its ungapped lambda is {ungappedLambda.ToString("F4", CultureInfo.InvariantCulture)}, the NCBI table has {ungappedRow[3].ToString("F4", CultureInfo.InvariantCulture)}): the tabulated gapped parameters do not apply.", nameof(scoring));
        return SpougeSignificance(score, scoring.Matrix.Name, existence, extension, queryLength, subjectLength, databaseLength);
    }

    /// <summary><see cref="ForMatrix"/> for a match / mismatch scheme (blastn): the tabulated gapped parameters with the length adjustment.</summary>
    /// <param name="score">The raw alignment score.</param>
    /// <param name="scoring">Match, mismatch, gap open and gap extension (the library convention: a gap of length k scores GapOpen + k x GapExtend).</param>
    /// <param name="queryLength">Query length.</param>
    /// <param name="subjectLength">Subject (database) length.</param>
    /// <param name="databaseSequences">Number of database sequences (default 1).</param>
    /// <exception cref="ArgumentException">The reward / penalty / gap costs are not in the blastn tables.</exception>
    public static AlignmentSignificanceResult ForNucleotide(double score, ScoringMatrix scoring, long queryLength, long subjectLength, long databaseSequences = 1)
    {
        ArgumentNullException.ThrowIfNull(scoring);
        int existence = -scoring.GapOpen, extension = -scoring.GapExtend;
        if (!TryGappedNucleotide(scoring.Match, scoring.Mismatch, existence, extension, out var p, out bool roundDown))
            throw new ArgumentException($"No blastn parameters for reward {scoring.Match}, penalty {scoring.Mismatch} and gap costs {existence} / {extension}.", nameof(scoring));
        var result = KarlinAltschulSignificance(roundDown && ((long)score & 1) == 1 ? score - 1 : score, p, queryLength, subjectLength, databaseSequences);
        return result with { Score = score, BitScore = p.BitScore(score) };
    }

    // ------------------------------------------------------------------ Z-score

    /// <summary>
    /// The shuffle Z-score of a local alignment (the idea of EMBOSS <c>prss</c> and FASTA <c>ssearch</c> "z-score"): the best local score of the pair is compared with the scores of <paramref name="shuffles"/> random
    /// permutations of <paramref name="sequence2"/> (a seeded Fisher-Yates shuffle, so the value is reproducible): z = (S - mean) / standard deviation of the shuffled scores.
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence (shuffled).</param>
    /// <param name="scoring">Scalar local scoring.</param>
    /// <param name="shuffles">Number of shuffles (at least 2; default 1000).</param>
    /// <param name="seed">Seed of the shuffles.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The score of the pair, the mean and standard deviation of the shuffled scores and the Z-score (0 when every shuffled score is equal).</returns>
    public static (int Score, double Mean, double StandardDeviation, double ZScore) ShuffleZScore(
        string sequence1, string sequence2, ScoringMatrix? scoring = null, int shuffles = 1000, int seed = 1, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentOutOfRangeException.ThrowIfLessThan(shuffles, 2);
        int score = SequenceAligner.LocalAlign(sequence1, sequence2, scoring, cancellationToken).Score;
        var rng = new Random(seed);
        var chars = sequence2.ToCharArray();
        double sum = 0, sumSq = 0;
        for (int i = 0; i < shuffles; i++)
        {
            for (int j = chars.Length - 1; j > 0; j--)
            {
                int k = rng.Next(j + 1);
                (chars[j], chars[k]) = (chars[k], chars[j]);
            }

            int s = SequenceAligner.LocalAlign(sequence1, new string(chars), scoring, cancellationToken).Score;
            sum += s;
            sumSq += (double)s * s;
        }

        double mean = sum / shuffles;
        double variance = Math.Max(0, (sumSq - (shuffles * mean * mean)) / (shuffles - 1));
        double sd = Math.Sqrt(variance);
        return (score, mean, sd, sd > 0 ? (score - mean) / sd : 0);
    }
}
