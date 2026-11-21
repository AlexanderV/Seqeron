using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>mgaps</c> of MUMmer 3.23 as a text program: reads the output of <c>mummer</c> (a <c>&gt;</c> line per query record and strand and the matches under it as three numbers), clusters the matches of every block by diagonal and separation
/// (<see cref="NucmerClusters"/>, the engine <c>nucmer</c> and <c>promer</c> use) and writes the label of the block and its clusters (<see cref="MummerMgapsText"/>).
/// </summary>
public static class MummerMgaps
{
    private const int MaxLine = 1024;

    /// <summary>Runs <c>mgaps [-d diag] [-f factor] [-l minlen] [-s separation] [-e] &lt; matches</c>.</summary>
    /// <param name="matches">The standard input: the text <c>mummer</c> printed.</param>
    /// <param name="options">The switches (null = the defaults of <c>mgaps</c> itself: <c>-d 5 -f 0.05 -l 200 -s 1000</c>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="FormatException">The input has no <c>&gt;</c> line (the program prints an uninitialised label buffer).</exception>
    public static string Run(string matches, NucmerClusterOptions? options = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(matches);
        var o = options ?? new NucmerClusterOptions { MaxSeparation = 1000, MinClusterLength = 200, DiagonalFraction = 0.05 };
        var sb = new StringBuilder();
        string? label = null;
        var block = new List<(long, long, long)>();
        void Flush()
        {
            if (label is null)
                return;
            if (block.Count == 0)
                sb.Append(label);
            else
                MummerMgapsText.AppendBlock(sb, label.EndsWith('\n') ? label[..^1] : label, NucmerClusters.Cluster(block, o, cancellationToken));
        }

        var v = new long[3];
        foreach (string line in MummerCText.Fgets(matches, MaxLine))
        {
            if (line[0] == '>')
            {
                Flush();
                block.Clear();
                label = line;
            }
            else if (MummerCText.ScanLongs(line, v) == 3)
            {
                block.Add((v[0], v[1], v[2]));
            }
        }

        if (label is null && block.Count > 0)
            throw new FormatException("mgaps: matches before the first '>' line: the program prints an uninitialised label (undefined).");
        Flush();
        return sb.ToString();
    }
}
