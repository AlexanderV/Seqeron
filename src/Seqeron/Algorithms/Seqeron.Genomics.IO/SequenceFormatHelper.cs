using System.Globalization;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Shared parsing logic for GenBank and EMBL format parsers.
/// </summary>
internal static partial class SequenceFormatHelper
{
    /// <summary>
    /// Splits a FASTA/FASTQ title line (the text after the leading '&gt;' or '@') into the record id and
    /// the free-text description. The id is the first whitespace-delimited word and the description is the
    /// remainder after that single whitespace character — Biopython <c>SeqIO</c> semantics
    /// (<c>title.split(None, 1)[0]</c>; any whitespace, not only a space). Shared by
    /// <see cref="FastaParser"/> and <see cref="FastqParser"/> so both formats derive ids identically.
    /// </summary>
    /// <param name="title">The title line without its leading marker; expected to be already trimmed.</param>
    /// <returns>The id and the description, or <c>null</c> description when the title has no whitespace.</returns>
    internal static (string Id, string? Description) SplitTitle(string title)
    {
        for (int i = 0; i < title.Length; i++)
        {
            if (char.IsWhiteSpace(title[i]))
                return (title.Substring(0, i), title.Substring(i + 1));
        }
        return (title, null);
    }

    /// <summary>
    /// Parses feature location parts from a location string.
    /// </summary>
    /// <param name="locationStr">Raw location string (e.g., "complement(join(1..100,200..300))").</param>
    /// <param name="useStartsWithForComplement">
    /// If true, detects complement via StartsWith (GenBank convention).
    /// If false, detects complement via Contains (EMBL convention).
    /// </param>
    /// <returns>Parsed location components.</returns>
    internal static (int Start, int End, bool IsComplement, bool IsJoin, bool IsOrder,
        bool Is5PrimePartial, bool Is3PrimePartial, IReadOnlyList<(int Start, int End)> Parts)
        ParseLocationParts(string locationStr, bool useStartsWithForComplement)
    {
        bool isComplement = useStartsWithForComplement
            ? locationStr.StartsWith("complement(", StringComparison.OrdinalIgnoreCase)
            : locationStr.Contains("complement(", StringComparison.OrdinalIgnoreCase);
        bool isJoin = locationStr.Contains("join(", StringComparison.OrdinalIgnoreCase);
        bool isOrder = locationStr.Contains("order(", StringComparison.OrdinalIgnoreCase);

        // INSDC 3.4.2.1: '<' indicates partial on 5' end, '>' indicates partial on 3' end.
        // These characters appear exclusively as partial indicators in INSDC location syntax.
        bool is5PrimePartial = locationStr.Contains('<');
        bool is3PrimePartial = locationStr.Contains('>');

        var parts = new List<(int Start, int End)>();

        var rangeMatches = LocationRangeRegex().Matches(locationStr);
        foreach (Match match in rangeMatches)
        {
            int start = int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
            int end = match.Groups[2].Success
                ? int.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture)
                : start;
            parts.Add((start, end));
        }

        int overallStart = parts.Count > 0 ? parts.Min(p => p.Start) : 0;
        int overallEnd = parts.Count > 0 ? parts.Max(p => p.End) : 0;

        return (overallStart, overallEnd, isComplement, isJoin, isOrder,
            is5PrimePartial, is3PrimePartial, parts);
    }

    [GeneratedRegex(@"(\d+)(?:\.\.(\d+))?")]
    internal static partial Regex LocationRangeRegex();
}
