using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

public static partial class MultipleAlignmentStatistics
{
    private static string DistanceText(float d)
    {
        // ajFmtPrintF("%.2f") of a float: NaN prints as nan / -nan (the sign bit), infinities as inf
        if (float.IsNaN(d))
            return float.IsNegative(d) ? "-nan" : "nan";
        return CFormat.Fixed(d, 2);
    }

    /// <summary>
    /// The <c>-outfile</c> text of EMBOSS <c>distmat</c> 6.6.0 (<c>distmat_outputDist</c>, EMBOSS-6.6.0 emboss/distmat.c via raw.githubusercontent.com; audit A5-8b, B09 F184): the title, the method line, <c>Using ambiguity codes</c>,
    /// <c>Using base positions N in the codon</c> (nucleotides), <c>Gap weighting is</c>, the column numbers and the upper triangle of the matrix (x100, two decimals, tab-separated, padded by magnitude) with the name and number of each row.
    /// </summary>
    /// <param name="names">The sequence names (<c>ajSeqsetGetseqNameS</c>).</param>
    /// <param name="distances">The matrix of <see cref="CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceOptions)"/> (unscaled fractions).</param>
    /// <param name="method">The method number as the program's menu (0 uncorrected ... 5 Jin-Nei).</param>
    /// <param name="options">The options that were used (<see cref="MsaDistanceOptions.Ambiguous"/>, <see cref="MsaDistanceOptions.GapWeight"/>, <see cref="MsaDistanceOptions.Position"/>).</param>
    /// <param name="nucleotide">The set is nucleotide (the codon-position line is printed).</param>
    /// <exception cref="ArgumentException">The names and the matrix do not agree, or the position is not 1, 2, 3, 12 or 123.</exception>
    public static string DistanceMatrixFile(IReadOnlyList<string> names, double[,] distances, int method, MsaDistanceOptions options, bool nucleotide)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(distances);
        ArgumentNullException.ThrowIfNull(options);
        int n = names.Count;
        if (distances.GetLength(0) != n || distances.GetLength(1) != n)
            throw new ArgumentException("The matrix must be square with one row per name.", nameof(distances));
        int posn = options.Position;
        int incr = 1;
        if (posn is >= 1 and <= 3 && nucleotide)
        {
            incr = 3;
            posn--;
        }
        else if (posn == 123)
            posn = 0;
        else if (posn != 12)
            throw new ArgumentException("Choose base positions 1, 2, 3, 12, or 123", nameof(options));
        if (posn == 0)
            posn = incr == 3 ? 1 : 123;
        else if (posn is 1 or 2)
            posn++;
        var sb = new StringBuilder();
        sb.Append("Distance Matrix\n---------------\n\n");
        string[] titles =
        [
            "Uncorrected for Multiple Substitutions", "Using the Jukes-Cantor correction method", "Using the Kimura correction method", "Using the Tamura correction method",
            "Using the Tajima-Nei correction method", "Using the Jin-Nei correction method",
        ];
        if (method is >= 0 and <= 5)
            sb.Append(titles[method]).Append('\n');
        if (options.Ambiguous)
            sb.Append("Using ambiguity codes\n");
        if (nucleotide)
            sb.Append("Using base positions ").Append(posn.ToString(CultureInfo.InvariantCulture)).Append(" in the codon\n");
        sb.Append("Gap weighting is ").Append(CFormat.Fixed((float)options.GapWeight, 6)).Append("\n\n");
        for (int j = 0; j < n; j++)
            sb.Append("\t    ").Append(j + 1);
        sb.Append('\n');
        for (int j = 0; j < n; j++)
        {
            sb.Append('\t');
            for (int i = 0; i < n; i++)
            {
                if (i < j)
                {
                    sb.Append('\t');
                    continue;
                }

                float d = i == j ? 0f : (float)distances[j, i];
                d *= 100f;
                string text = DistanceText(d);
                if (d < 10f)
                    sb.Append("  ").Append(text).Append('\t');
                else if (d < 100f)
                    sb.Append(' ').Append(text).Append('\t');
                else
                    sb.Append(text).Append('\t');
            }

            sb.Append('\t').Append(names[j]).Append(' ').Append(j + 1).Append('\n');
        }

        return sb.ToString();
    }
}
