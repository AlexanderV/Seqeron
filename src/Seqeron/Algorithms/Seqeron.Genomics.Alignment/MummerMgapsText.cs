using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>The text <c>mgaps</c> (MUMmer 3.23) prints for the clusters of one block of its input: the label line before the first cluster, a <c>#</c> line before every further one, the matches as <c>%8ld %8ld %6ld %7s %6s %6s</c>.</summary>
public static class MummerMgapsText
{
    /// <summary>Appends the clusters of one input block (one query record and strand) to <paramref name="sb"/>; a block without any cluster is just its label line (<c>mgaps</c> prints it).</summary>
    /// <param name="sb">The output.</param>
    /// <param name="headerLine">The label line of the block, e.g. <c>&gt; id</c> or <c>&gt; id Reverse</c> (without the newline).</param>
    /// <param name="clusters">The clusters of the block, in the order <see cref="NucmerClusters.Cluster"/> returns them.</param>
    public static void AppendBlock(StringBuilder sb, string headerLine, IEnumerable<NucmerCluster> clusters)
    {
        ArgumentNullException.ThrowIfNull(sb);
        ArgumentNullException.ThrowIfNull(headerLine);
        ArgumentNullException.ThrowIfNull(clusters);
        var ci = CultureInfo.InvariantCulture;
        bool first = true;
        foreach (var cluster in clusters)
        {
            if (first)
                sb.Append(headerLine).Append('\n');
            else
                sb.Append("#\n");
            first = false;
            bool firstMatch = true;
            long prevEnd1 = 0, prevEnd2 = 0;
            foreach (var m in cluster.Matches)
            {
                sb.Append(m.Start1.ToString(ci).PadLeft(8)).Append(' ').Append(m.Start2.ToString(ci).PadLeft(8)).Append(' ').Append(m.Length.ToString(ci).PadLeft(6));
                if (firstMatch)
                    sb.Append(' ').Append("none".PadLeft(7)).Append(' ').Append("-".PadLeft(6)).Append(' ').Append("-".PadLeft(6)).Append('\n');
                else
                    sb.Append(' ').Append((m.Overlap == 0 ? "none" : (-m.Overlap).ToString(ci)).PadLeft(7)).Append(' ')
                        .Append((m.Start1 - prevEnd1).ToString(ci).PadLeft(6)).Append(' ')
                        .Append((m.Start2 - prevEnd2).ToString(ci).PadLeft(6)).Append('\n');
                firstMatch = false;
                prevEnd1 = m.Start1 + m.Length;
                prevEnd2 = m.Start2 + m.Length;
            }
        }

        if (first)
            sb.Append(headerLine).Append('\n');
    }
}
