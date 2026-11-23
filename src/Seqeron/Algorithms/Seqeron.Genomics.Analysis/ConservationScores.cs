#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The per-column conservation scores of Capra and Singh (2007), "Predicting functionally important residues from sequence conservation", Bioinformatics 23:1875 (B09 WP67, F308, R14-13): a port of the functions of the
/// authors' <c>score_conservation.py</c> (compbio.cs.princeton.edu/conservation, version 7/23/07) — <see cref="ShannonEntropy"/>, <see cref="PropertyEntropy"/>, <see cref="RelativeEntropy"/>, <see cref="JsDivergence"/> (Jensen-Shannon
/// divergence of the column distribution from a background distribution), <see cref="SumOfPairs"/> (the sum-of-pairs column score with a similarity matrix), <see cref="WindowScore"/> and the sequence weights of Henikoff and Henikoff
/// (<see cref="SequenceWeights"/>, 1994). Every score takes a column as a string, the sequence weights of the rows (null or of another length than the column: every row 1) and the gap penalty switch (true: the score times the
/// weighted fraction of non-gap rows). The alphabet of the column distribution is <see cref="AminoAcids"/> (20 residues in the order of the script, then the gap): a letter outside it is not counted by the entropies and is an error
/// (<see cref="KeyNotFoundException"/>, Python <c>KeyError</c>) in <see cref="SumOfPairs"/> and <see cref="SequenceWeights"/>. The Shannon score is the script's information score <c>1 - h</c> with h the entropy in nats divided by
/// <c>ln(min(21, number of rows))</c> (a column of one row divides by zero as in the script).
/// </summary>
public static class ConservationScores
{
    /// <summary>The pseudo count of the script (1e-7) added to the weighted count of every symbol.</summary>
    public const double PseudoCount = 0.0000001;

    /// <summary>The 20 amino acids in the order of the script followed by the gap symbol.</summary>
    public const string AminoAcids = "ARNDCQEGHILKMFPSTWYV-";

    /// <summary>The BLOSUM62 background distribution of the script (in the order of <see cref="AminoAcids"/>, without the gap).</summary>
    public static IReadOnlyList<double> BlosumBackground { get; } = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072];

    /// <summary>The SwissProt background distribution of the script.</summary>
    public static IReadOnlyList<double> SwissProtBackground { get; } = [.0788, .054, .0414, .0535, .0151, .0395, .0668, .0695, .0229, .0590, .0963, .0594, .0237, .0397, .0482, .0682, .0540, .0113, .0304, .0673];

    // the amino-acid groups of Mirny and Shakhnovich 1999 used by the script's property entropy (the last group is the gap)
    private static readonly string[] PropertyPartition = ["AVLIMC", "FWYH", "STNQ", "KR", "DE", "GP", "-"];

    /// <summary>The 20 x 20 BLOSUM62 matrix in the order of <see cref="AminoAcids"/> (the matrix <c>sim_matrix</c> of the script's default).</summary>
    public static double[,] Blosum62()
    {
        var m = SubstitutionMatrix.Get("BLOSUM62");
        var result = new double[20, 20];
        for (int i = 0; i < 20; i++)
            for (int j = 0; j < 20; j++)
                result[i, j] = m.TryGetScore(AminoAcids[i], AminoAcids[j], out int s) ? s : 0;
        return result;
    }

    private static double[] Weights(string column, IReadOnlyList<double>? weights)
    {
        if (weights is null || weights.Count != column.Length)
            return Enumerable.Repeat(1.0, column.Length).ToArray();
        return [.. weights];
    }

    private static double Sum(IEnumerable<double> values)
    {
        double s = 0;
        foreach (double v in values)
            s += v;
        return s;
    }

    /// <summary><c>weighted_freq_count_pseudocount</c>: the weighted frequency of each of the 21 symbols of <see cref="AminoAcids"/> in the column with the pseudo count.</summary>
    public static double[] WeightedFreqCountPseudocount(string column, IReadOnlyList<double>? weights = null, double pseudoCount = PseudoCount)
    {
        var w = Weights(column, weights);
        var counts = Enumerable.Repeat(pseudoCount, AminoAcids.Length).ToArray();
        for (int a = 0; a < AminoAcids.Length; a++)
            for (int j = 0; j < column.Length; j++)
                if (column[j] == AminoAcids[a])
                    counts[a] += 1 * w[j];
        double denominator = Sum(w) + (AminoAcids.Length * pseudoCount);
        for (int j = 0; j < counts.Length; j++)
            counts[j] /= denominator;
        return counts;
    }

    /// <summary><c>weighted_gap_penalty</c>: 1 minus the weighted fraction of gap rows of the column.</summary>
    public static double WeightedGapPenalty(string column, IReadOnlyList<double>? weights = null)
    {
        var w = Weights(column, weights);
        double gapSum = 0;
        for (int i = 0; i < column.Length; i++)
            if (column[i] == '-')
                gapSum += w[i];
        double total = Sum(w);
        if (total == 0)
            throw new DivideByZeroException("float division by zero");
        return 1 - (gapSum / total);
    }

    private static double PyLog(double x)
    {
        if (x <= 0)
            throw new ArgumentException("math domain error");
        return Math.Log(x);
    }

    private static double Divide(double a, double b) => b == 0 ? throw new DivideByZeroException("float division by zero") : a / b;

    /// <summary><c>shannon_entropy</c>: the information score <c>1 - h</c>, h = sum f ln f / ln(min(21, rows)); 1 is a fully conserved column, 0 the maximal entropy.</summary>
    public static double ShannonEntropy(string column, IReadOnlyList<double>? weights = null, bool gapPenalty = true)
    {
        var fc = WeightedFreqCountPseudocount(column, weights);
        double h = 0;
        for (int i = 0; i < fc.Length; i++)
            if (fc[i] != 0)
                h += fc[i] * Math.Log(fc[i]);
        h = Divide(h, PyLog(Math.Min(fc.Length, column.Length)));
        double info = 1 - (-1 * h);
        return gapPenalty ? info * WeightedGapPenalty(column, weights) : info;
    }

    /// <summary><c>property_entropy</c>: the information score of the column over the seven groups of Mirny and Shakhnovich (1999): AVLIMC, FWYH, STNQ, KR, DE, GP and the gap.</summary>
    public static double PropertyEntropy(string column, IReadOnlyList<double>? weights = null, bool gapPenalty = true)
    {
        var fc = WeightedFreqCountPseudocount(column, weights);
        var prop = new double[PropertyPartition.Length];
        for (int p = 0; p < PropertyPartition.Length; p++)
            foreach (char aa in PropertyPartition[p])
                prop[p] += fc[AminoAcids.IndexOf(aa)];
        double h = 0;
        for (int i = 0; i < prop.Length; i++)
            if (prop[i] != 0)
                h += prop[i] * Math.Log(prop[i]);
        h = Divide(h, PyLog(Math.Min(prop.Length, column.Length)));
        double info = 1 - (-1 * h);
        return gapPenalty ? info * WeightedGapPenalty(column, weights) : info;
    }

    // a background of 20 values has no gap: the gap frequency is removed from the column distribution and the rest renormalised
    private static double[] Distribution(string column, IReadOnlyList<double>? weights, int backgroundLength)
    {
        var fc = WeightedFreqCountPseudocount(column, weights);
        if (backgroundLength == 20)
        {
            var kept = fc[..^1];
            double s = Sum(kept);
            for (int i = 0; i < kept.Length; i++)
                kept[i] /= s;
            fc = kept;
        }

        return fc;
    }

    /// <summary><c>relative_entropy</c>: the Kullback-Leibler divergence of the column distribution from the background, divided by ln(20); -1 when the background has neither 20 nor 21 values.</summary>
    public static double RelativeEntropy(string column, IReadOnlyList<double> background, IReadOnlyList<double>? weights = null, bool gapPenalty = true)
    {
        ArgumentNullException.ThrowIfNull(background);
        var fc = Distribution(column, weights, background.Count);
        if (fc.Length != background.Count)
            return -1;
        double d = 0;
        for (int i = 0; i < fc.Length; i++)
            if (background[i] != 0.0)
                d += fc[i] * PyLog(fc[i] / background[i]);
        d = Divide(d, PyLog(fc.Length));
        return gapPenalty ? d * WeightedGapPenalty(column, weights) : d;
    }

    // math.log(x, 2) of CPython: log(x) / log(2)
    private static double Log2(double x) => PyLog(x) / Math.Log(2);

    /// <summary>
    /// <c>js_divergence</c>: the Jensen-Shannon divergence of the column distribution f from the background b (in bits): with r = (f + b) / 2, half the sum of f log2(f / r) + b log2(b / r) (a zero entry contributes nothing); the
    /// divergence is 0 for a column equal to the background and approaches 1 for a conserved column of a rare residue. -1 when the background has neither 20 nor 21 values.
    /// </summary>
    public static double JsDivergence(string column, IReadOnlyList<double> background, IReadOnlyList<double>? weights = null, bool gapPenalty = true)
    {
        ArgumentNullException.ThrowIfNull(background);
        var fc = Distribution(column, weights, background.Count);
        if (fc.Length != background.Count)
            return -1;
        var r = new double[fc.Length];
        for (int i = 0; i < fc.Length; i++)
            r[i] = (.5 * fc[i]) + (.5 * background[i]);
        double d = 0;
        for (int i = 0; i < fc.Length; i++)
        {
            if (r[i] == 0.0)
                continue;
            if (fc[i] == 0.0)
                d += background[i] * Log2(background[i] / r[i]);
            else if (background[i] == 0.0)
                d += fc[i] * Log2(fc[i] / r[i]);
            else
                d += (fc[i] * Log2(fc[i] / r[i])) + (background[i] * Log2(background[i] / r[i]));
        }

        d /= 2;
        return gapPenalty ? d * WeightedGapPenalty(column, weights) : d;
    }

    /// <summary>
    /// <c>sum_of_pairs</c>: the weighted mean similarity of all pairs of non-gap residues of the column (the matrix is indexed in the order of <see cref="AminoAcids"/>, <see cref="Blosum62"/> by default); 0 when fewer than two residues.
    /// </summary>
    public static double SumOfPairs(string column, double[,]? similarity = null, IReadOnlyList<double>? weights = null, bool gapPenalty = true)
    {
        var m = similarity ?? Blosum62();
        var w = Weights(column, weights);
        double sum = 0, maxSum = 0;
        for (int i = 0; i < column.Length; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (column[i] != '-' && column[j] != '-')
                {
                    maxSum += w[i] * w[j];
                    sum += w[i] * w[j] * m[Index(column[i]), Index(column[j])];
                }
            }
        }

        sum = maxSum != 0 ? sum / maxSum : 0.0;
        return gapPenalty ? sum * WeightedGapPenalty(column, weights) : sum;
    }

    private static int Index(char c)
    {
        int k = AminoAcids.IndexOf(c);
        return k < 0 ? throw new KeyNotFoundException($"'{c}'") : k;
    }

    /// <summary>
    /// <c>window_score</c>: every score of the interior (positions window .. n - window - 1) becomes 0.5 of itself plus 0.5 of the mean of the other non-negative scores within <paramref name="windowLength"/> positions on either side;
    /// negative scores (columns that could not be scored) are left alone and ignored.
    /// </summary>
    public static double[] WindowScore(IReadOnlyList<double> scores, int windowLength)
    {
        var result = new double[scores.Count];
        for (int i = 0; i < result.Length; i++)
            result[i] = scores[i];
        const double lam = .5;
        for (int i = windowLength; i < scores.Count - windowLength; i++)
        {
            if (scores[i] < 0)
                continue;
            double sum = 0, terms = 0;
            for (int j = i - windowLength; j < i + windowLength + 1; j++)
            {
                if (i != j && scores[j] >= 0)
                {
                    terms += 1;
                    sum += scores[j];
                }
            }

            if (terms > 0)
                result[i] = ((1 - lam) * (sum / terms)) + (lam * scores[i]);
        }

        return result;
    }

    /// <summary><c>calculate_sequence_weights</c>: the position-based sequence weights of Henikoff and Henikoff (1994): per column every row gets 1 / (number of different residues x number of rows with its residue), averaged over the columns.</summary>
    public static double[] SequenceWeights(IReadOnlyList<string> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        if (rows.Count == 0)
            throw new ArgumentOutOfRangeException(nameof(rows), "list index out of range");
        var weights = new double[rows.Count];
        int length = rows[0].Length;
        for (int i = 0; i < length; i++)
        {
            var counts = new int[AminoAcids.Length];
            foreach (var row in rows)
                if (row[i] != '-')
                    counts[Index(row[i])]++;
            int types = counts.Count(c => c > 0);
            for (int j = 0; j < rows.Count; j++)
            {
                int d = counts[Index(rows[j][i])] * types;
                if (d > 0)
                    weights[j] += 1.0 / d;
            }
        }

        for (int w = 0; w < weights.Length; w++)
            weights[w] /= length;
        return weights;
    }
}
