using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The Karlin-Altschul parameters of a scoring system (audit R14-7, B09 F300): <see cref="Lambda"/>, <see cref="K"/>, the relative entropy <see cref="H"/> and, for gapped scoring from the NCBI tables, the
/// length-correction parameters <see cref="Alpha"/> / <see cref="Beta"/>. Bit score S' = (lambda S - ln K) / ln 2, E = search space x K exp(-lambda S), P = 1 - exp(-E).
/// </summary>
/// <param name="Lambda">The scale of the score distribution (nats per score unit).</param>
/// <param name="K">The prefactor of the extreme-value distribution.</param>
/// <param name="H">The relative entropy of the scoring system in nats per aligned pair.</param>
/// <param name="Alpha">BLAST <c>alpha</c> (the expected alignment length grows as alpha / lambda x ln(search space); for ungapped scoring lambda / H); NaN when not tabulated.</param>
/// <param name="Beta">BLAST <c>beta</c> (the additive length correction; 0 for ungapped scoring); NaN when not tabulated.</param>
public sealed record KarlinAltschulParameters(double Lambda, double K, double H, double Alpha = double.NaN, double Beta = double.NaN)
{
    /// <summary>ln K.</summary>
    public double LogK => Math.Log(K);

    /// <summary>The normalised score S' = (lambda S - ln K) / ln 2 (NCBI <c>bitscore</c>).</summary>
    public double BitScore(double score) => ((Lambda * score) - LogK) / Math.Log(2);

    /// <summary>The expected number of alignments of score at least <paramref name="score"/> in a search space of the given size: space x K exp(-lambda S).</summary>
    public double EValue(double score, double searchSpace) => searchSpace * Math.Exp((-Lambda * score) + LogK);

    /// <summary>The probability of at least one alignment with that E-value: 1 - exp(-E) (NCBI <c>BLAST_KarlinEtoP</c>).</summary>
    public static double PValue(double eValue) => -AlignmentSignificance.ExpM1(-eValue);
}

/// <summary>The significance of one local alignment score (<see cref="AlignmentSignificance"/>).</summary>
/// <param name="Score">The raw score.</param>
/// <param name="BitScore">The bit score of the parameters.</param>
/// <param name="EValue">The expected number of chance alignments with at least this score.</param>
/// <param name="PValue">1 - exp(-E).</param>
/// <param name="Method">Which formula gave the E-value: <c>SpougeFiniteSizeCorrection</c> (BLAST+ for the standard protein matrices) or <c>LengthAdjustedKarlinAltschul</c>.</param>
/// <param name="LengthAdjustment">The edge-effect length adjustment subtracted from both lengths (0 with the finite-size correction).</param>
/// <param name="SearchSpace">The effective search space (query x database length after the adjustment; the area of the Spouge correction).</param>
/// <param name="Parameters">The Karlin-Altschul parameters used.</param>
public sealed record AlignmentSignificanceResult(
    double Score,
    double BitScore,
    double EValue,
    double PValue,
    string Method,
    int LengthAdjustment,
    double SearchSpace,
    KarlinAltschulParameters Parameters);

/// <summary>
/// Alignment significance statistics (audit R14-7, B09 F300): the Karlin-Altschul parameters of an integer scoring system computed from the score distribution (lambda, K, H of
/// Karlin &amp; Altschul 1990, PNAS 87:2264; the algorithms of NCBI <c>blast_stat.c</c>: <c>NlmKarlinLambdaNR</c>, <c>BlastKarlinLtoH</c>, <c>BlastKarlinLHtoK</c>), the NCBI tables of the gapped parameters of the
/// standard protein matrices and of the blastn reward / penalty pairs, the BLAST length adjustment (<c>BLAST_ComputeLengthAdjustment</c>), Spouge's finite-size-corrected E-value
/// (<c>BLAST_SpougeStoE</c>, Spouge 2002; what BLAST+ prints for the protein matrices) and the shuffle Z-score (the EMBOSS <c>prss</c> / FASTA idea).
/// </summary>
/// <remarks>
/// Verified against BLAST+ 2.12.0 (<c>blastp -comp_based_stats 0 -seg no</c> and <c>blastn -task blastn -dust no</c> with <c>-subject</c>): the bit scores and E-values of the same alignment score are equal,
/// and against the ungapped values printed by <c>blastp -ungapped</c> (BLOSUM62 lambda 0.3176, K 0.134, H 0.4012). BLAST+ composition-based statistics (the default of blastp) change lambda per alignment and are not
/// modelled: the values are those of <c>-comp_based_stats 0</c>.
/// </remarks>
public static partial class AlignmentSignificance
{
    // amino acid background frequencies of Robinson and Robinson 1991 as NCBI blast_stat.c lists them (per mille); normalised to 1 when used
    private static readonly (char Residue, double Weight)[] RobinsonFrequencies =
    [
        ('A', 78.05), ('C', 19.25), ('D', 53.64), ('E', 62.95), ('F', 38.56), ('G', 73.77), ('H', 21.99), ('I', 51.42), ('K', 57.44), ('L', 90.19),
        ('M', 22.43), ('N', 44.87), ('P', 52.03), ('Q', 42.64), ('R', 51.29), ('S', 71.20), ('T', 58.41), ('V', 64.41), ('W', 13.30), ('Y', 32.16),
    ];

    /// <summary>The Robinson &amp; Robinson amino-acid background frequencies of NCBI BLAST, normalised to sum 1.</summary>
    public static IReadOnlyDictionary<char, double> RobinsonBackground { get; } = BuildRobinson();

    private static Dictionary<char, double> BuildRobinson()
    {
        double sum = RobinsonFrequencies.Sum(f => f.Weight);
        return RobinsonFrequencies.ToDictionary(f => f.Residue, f => f.Weight / sum);
    }

    // exp(x) - 1 without the cancellation of Math.Exp(x) - 1 for a small x (double.ExpM1 of .NET is exactly that subtraction)
    internal static double ExpM1(double x)
    {
        if (Math.Abs(x) > 0.5)
            return Math.Exp(x) - 1.0;
        double term = x, sum = x;
        for (int n = 2; n < 40; n++)
        {
            term *= x / n;
            sum += term;
            if (Math.Abs(term) < 1e-18 * Math.Abs(sum))
                break;
        }

        return sum;
    }

    // ------------------------------------------------------------------ ungapped parameters

    /// <summary>
    /// The ungapped Karlin-Altschul parameters of a substitution matrix under a residue background (default: <see cref="RobinsonBackground"/> for a protein matrix, uniform A / C / G / T for a nucleotide matrix),
    /// the "ideal" block of BLAST (<c>Blast_ScoreBlkKbpIdealCalc</c>): the score distribution is the sum of p_i p_j over the pairs with that score.
    /// </summary>
    /// <param name="matrix">The matrix (integer scores).</param>
    /// <param name="background">Residue frequencies (they are normalised); null = the default of the alphabet.</param>
    /// <exception cref="ArgumentNullException"><paramref name="matrix"/> is null.</exception>
    /// <exception cref="ArgumentException">The expected score is not negative or no score is positive (the statistics do not exist).</exception>
    public static KarlinAltschulParameters Ungapped(SubstitutionMatrix matrix, IReadOnlyDictionary<char, double>? background = null)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        var freq = background ?? (matrix.IsNucleotide
            ? new Dictionary<char, double> { ['A'] = 0.25, ['C'] = 0.25, ['G'] = 0.25, ['T'] = 0.25 }
            : RobinsonBackground);
        double total = freq.Values.Sum();
        var pairs = new List<(int Score, double Probability)>();
        foreach (var (a, pa) in freq)
        {
            foreach (var (b, pb) in freq)
            {
                if (pa > 0 && pb > 0 && matrix.TryGetScore(a, b, out int s))
                    pairs.Add((s, pa / total * (pb / total)));
            }
        }

        return FromScoreDistribution(pairs);
    }

    /// <summary>The ungapped parameters of a match / mismatch scheme over A / C / G / T (default uniform frequencies), e.g. blastn reward 2 penalty -3.</summary>
    /// <param name="reward">The score of an identical pair (positive).</param>
    /// <param name="penalty">The score of a different pair (negative).</param>
    /// <param name="baseFrequencies">Frequencies of A, C, G, T (normalised); null = 0.25 each.</param>
    public static KarlinAltschulParameters UngappedNucleotide(int reward, int penalty, double[]? baseFrequencies = null)
    {
        double[] f = baseFrequencies ?? [0.25, 0.25, 0.25, 0.25];
        if (f.Length != 4)
            throw new ArgumentException("Four frequencies (A, C, G, T) are required.", nameof(baseFrequencies));
        double total = f.Sum();
        var pairs = new List<(int, double)>();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
                pairs.Add((i == j ? reward : penalty, f[i] / total * (f[j] / total)));
        }

        return FromScoreDistribution(pairs);
    }

    /// <summary>The parameters of an explicit score distribution (the pairs of score and probability; equal scores add up).</summary>
    /// <param name="distribution">Pairs of an integer score and its probability.</param>
    /// <exception cref="ArgumentException">The mean score is not negative or the highest score is not positive.</exception>
    public static KarlinAltschulParameters FromScoreDistribution(IEnumerable<(int Score, double Probability)> distribution)
    {
        ArgumentNullException.ThrowIfNull(distribution);
        var byScore = new SortedDictionary<int, double>();
        foreach (var (score, p) in distribution)
            byScore[score] = byScore.GetValueOrDefault(score) + p;
        byScore = new SortedDictionary<int, double>(byScore.Where(kv => kv.Value > 0).ToDictionary(kv => kv.Key, kv => kv.Value));
        if (byScore.Count == 0)
            throw new ArgumentException("The score distribution is empty.", nameof(distribution));
        int low = byScore.Keys.First(), high = byScore.Keys.Last();
        double sum = byScore.Values.Sum();
        var sprob = new double[high - low + 1];
        double average = 0;
        foreach (var (score, p) in byScore)
        {
            sprob[score - low] = p / sum;
            average += score * (p / sum);
        }

        if (average >= 0 || high <= 0)
            throw new ArgumentException($"Karlin-Altschul statistics need a negative expected score and a positive maximum score (mean {average.ToString("R", CultureInfo.InvariantCulture)}, maximum {high}).", nameof(distribution));
        double lambda = KarlinLambda(sprob, low, high);
        double h = KarlinH(sprob, low, high, lambda);
        double k = KarlinK(sprob, low, high, lambda, h, average);
        return new KarlinAltschulParameters(lambda, k, h, lambda / h, 0);
    }

    private static int Gcd(int a, int b)
    {
        a = Math.Abs(a);
        b = Math.Abs(b);
        while (b != 0)
            (a, b) = (b, a % b);
        return a;
    }

    // probs[i - low] = probability of score i; the positive root of sum_i probs[i] exp(lambda i) = 1 (NCBI NlmKarlinLambdaNR: safeguarded Newton iteration on the polynomial in x = exp(-lambda))
    private static double KarlinLambda(double[] probs, int low, int high)
    {
        int d = -low;
        for (int i = 1; i <= high - low && d > 1; i++)
        {
            if (probs[i] != 0.0)
                d = Gcd(d, i);
        }

        double P(int score) => probs[score - low];
        const double tolx = 1e-5;
        const int itmax = 20 + 17;
        const int maxNewton = 20 + 17;
        double x0 = Math.Exp(-0.5);
        double x = x0 is > 0 and < 1 ? x0 : 0.5;
        double a = 0, b = 1, f = 4;
        bool isNewton = false;
        for (int k = 0; k < itmax; k++)
        {
            double fold = f;
            bool wasNewton = isNewton;
            isNewton = false;
            double g = 0;
            f = P(low);
            int i;
            for (i = low + d; i < 0; i += d)
            {
                g = (x * g) + f;
                f = (f * x) + P(i);
            }

            g = (x * g) + f;
            f = (f * x) + P(0) - 1;
            for (i = d; i <= high; i += d)
            {
                g = (x * g) + f;
                f = (f * x) + P(i);
            }

            if (f > 0)
                a = x;
            else if (f < 0)
                b = x;
            else
                break;
            if (b - a < 2 * a * (1 - b) * tolx)
            {
                x = (a + b) / 2;
                break;
            }

            if (k >= maxNewton || (wasNewton && Math.Abs(f) > 0.9 * Math.Abs(fold)) || g >= 0)
            {
                x = (a + b) / 2;
            }
            else
            {
                double p = -f / g;
                double y = x + p;
                if (y <= a || y >= b)
                {
                    x = (a + b) / 2;
                }
                else
                {
                    isNewton = true;
                    x = y;
                    if (Math.Abs(p) < tolx * x * (1 - x))
                        break;
                }
            }
        }

        return -Math.Log(x) / d;
    }

    // H = lambda sum_i i p_i exp(lambda i) (NCBI BlastKarlinLtoH)
    private static double KarlinH(double[] probs, int low, int high, double lambda)
    {
        double etonlam = Math.Exp(-lambda);
        double sum = low * probs[0];
        for (int score = low + 1; score <= high; score++)
            sum = (score * probs[score - low]) + (etonlam * sum);
        double scale = Math.Pow(etonlam, high);
        return scale > 0.0 ? lambda * sum / scale : lambda * Math.Exp((lambda * high) + Math.Log(sum));
    }

    // K of Karlin & Altschul 1990 (NCBI BlastKarlinLHtoK): the closed forms for steps of one and the iterated convolution of the score distribution otherwise
    private static double KarlinK(double[] sprob, int low0, int high0, double lambda0, double h, double scoreAverage)
    {
        const double sumLimit = 0.0001;
        const int iterLimit = 100;
        int low = low0, high = high0;
        int range0 = high - low;
        int divisor = -low;
        for (int i = 1; i <= range0 && divisor > 1; i++)
        {
            if (sprob[i] != 0.0)
                divisor = Gcd(divisor, i);
        }

        high /= divisor;
        low /= divisor;
        double lambda = lambda0 * divisor;
        int range = high - low;
        double firstTermClosedForm = h / lambda;
        double expMinusLambda = Math.Exp(-lambda);
        if (low == -1 && high == 1)
        {
            double pl = sprob[(low * divisor) - low0], ph = sprob[(high * divisor) - low0];
            return (pl - ph) * (pl - ph) / pl;
        }

        if (low == -1 || high == 1)
        {
            if (high != 1)
            {
                double average = scoreAverage / divisor;
                firstTermClosedForm = average * average / firstTermClosedForm;
            }

            return firstTermClosedForm * (1.0 - expMinusLambda);
        }

        // the scores are multiples of the divisor: work on the lattice scores / divisor (NCBI indexes the unscaled array there, which is only right for a divisor of 1)
        var probArrayStartLow = new double[range + 1];
        for (int j = 0; j <= range; j++)
            probArrayStartLow[j] = sprob[((low + j) * divisor) - low0];
        var alignmentScoreProbabilities = new double[(iterLimit * range) + 1];
        double outerSum = 0, innerSum = 1;
        int lowAlignmentScore = 0, highAlignmentScore = 0;
        alignmentScoreProbabilities[0] = 1;
        int iterCounter = 0;
        while (iterCounter < iterLimit && innerSum > sumLimit)
        {
            int first = range, last = range;
            lowAlignmentScore += low;
            highAlignmentScore += high;
            int ptrP;
            for (ptrP = highAlignmentScore - lowAlignmentScore; ptrP >= 0; ptrP--)
            {
                int ptr1 = ptrP - first;
                int ptr1e = ptrP - last;
                int ptr2 = first;
                innerSum = 0;
                while (ptr1 >= ptr1e)
                {
                    innerSum += alignmentScoreProbabilities[ptr1] * probArrayStartLow[ptr2];
                    ptr1--;
                    ptr2++;
                }

                alignmentScoreProbabilities[ptrP] = innerSum;
                if (first > 0)
                    first--;
                if (ptrP <= range)
                    last--;
            }

            // Horner's rule
            int idx = 0;
            innerSum = alignmentScoreProbabilities[idx];
            int ii;
            for (ii = lowAlignmentScore + 1; ii < 0; ii++)
                innerSum = alignmentScoreProbabilities[++idx] + (innerSum * expMinusLambda);
            innerSum *= expMinusLambda;
            for (; ii <= highAlignmentScore; ++ii)
                innerSum += alignmentScoreProbabilities[++idx];
            iterCounter++;
            innerSum /= iterCounter;
            outerSum += innerSum;
        }

        return -Math.Exp(-2.0 * outerSum) / (firstTermClosedForm * ExpM1(-lambda));
    }
}
