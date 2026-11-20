using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class Minimap2Chaining
{
    /// <summary>
    /// The PAF text of <see cref="Map"/>'s result as <c>minimap2</c> 2.26 prints it without a base-level alignment (<c>format.c</c> <c>mm_write_paf3</c> + <c>write_tags</c>): the 12 columns (query name / length / start / end, strand, target name /
    /// length / start / end, residue matches, alignment block length, mapping quality) and the tags <c>tp:A:P|S</c>, <c>cm:i:</c> (minimizers on the chain), <c>s1:i:</c> (chain score), <c>s2:i:</c> (second best chain score, primary hits only),
    /// <c>dv:f:</c> (divergence estimated from the minimizers, <c>0</c> or four decimals, omitted when not estimated) and <c>rl:i:</c> (length of the query covered by repetitive minimizers, on every line). A query without a hit prints nothing unless
    /// <paramref name="printUnmapped"/> (<c>--paf-no-hit</c>): then <c>name len 0 0 * * 0 0 0 0 0 0 rl:i:N</c>. Every line ends with a line feed.
    /// </summary>
    /// <param name="result">The mapping result.</param>
    /// <param name="queryName">The query name (column 1).</param>
    /// <param name="queryLength">The query length (column 2).</param>
    /// <param name="targetNames">The target names, indexed by <see cref="MinimapHit.TargetId"/> (the program prints the number when a target has no name).</param>
    /// <param name="targetLengths">The target lengths.</param>
    /// <param name="printUnmapped"><c>--paf-no-hit</c>.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static string FormatPaf(MinimapMapResult result, string queryName, int queryLength, IReadOnlyList<string?> targetNames, IReadOnlyList<int> targetLengths, bool printUnmapped = false)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(queryName);
        ArgumentNullException.ThrowIfNull(targetNames);
        ArgumentNullException.ThrowIfNull(targetLengths);
        var sb = new StringBuilder();
        if (result.Hits.Count == 0)
        {
            if (printUnmapped)
                sb.Append(queryName).Append('\t').Append(queryLength).Append("\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\trl:i:").Append(result.RepeatLength).Append('\n');
            return sb.ToString();
        }

        foreach (var h in result.Hits)
        {
            sb.Append(queryName).Append('\t').Append(queryLength).Append('\t').Append(h.QueryStart).Append('\t').Append(h.QueryEnd).Append('\t').Append(h.Reverse ? '-' : '+').Append('\t');
            sb.Append(h.TargetId < targetNames.Count && targetNames[h.TargetId] is { } name ? name : h.TargetId.ToString(System.Globalization.CultureInfo.InvariantCulture));
            sb.Append('\t').Append(targetLengths[h.TargetId]).Append('\t').Append(h.TargetStart).Append('\t').Append(h.TargetEnd);
            sb.Append('\t').Append(h.MatchLength).Append('\t').Append(h.BlockLength).Append('\t').Append(h.MapQ);
            sb.Append("\ttp:A:").Append(h.IsPrimary ? 'P' : 'S').Append("\tcm:i:").Append(h.AnchorCount).Append("\ts1:i:").Append(h.Score);
            if (h.IsPrimary)
                sb.Append("\ts2:i:").Append(h.SecondaryScore);
            if (h.Divergence >= 0f && h.Divergence <= 1f)
                sb.Append("\tdv:f:").Append(h.Divergence == 0f ? "0" : CFormat.Fixed(h.Divergence, 4));
            sb.Append("\trl:i:").Append(result.RepeatLength).Append('\n');
        }

        return sb.ToString();
    }
}
