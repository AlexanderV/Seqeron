using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Shared parsing logic for GenBank and EMBL format parsers.
/// </summary>
internal static partial class SequenceFormatHelper
{
    /// <summary>
    /// Parses feature location parts from a location string.
    /// </summary>
    /// <param name="locationStr">Raw location string (e.g., "complement(join(1..100,200..300))").</param>
    /// <param name="useStartsWithForComplement">
    /// If true, detects complement via StartsWith (GenBank convention).
    /// If false, detects complement via Contains (EMBL convention).
    /// </param>
    /// <returns>Parsed location components.</returns>
    internal static (int Start, int End, bool IsComplement, bool IsJoin, IReadOnlyList<(int Start, int End)> Parts)
        ParseLocationParts(string locationStr, bool useStartsWithForComplement)
    {
        bool isComplement = useStartsWithForComplement
            ? locationStr.StartsWith("complement(", StringComparison.OrdinalIgnoreCase)
            : locationStr.Contains("complement(", StringComparison.OrdinalIgnoreCase);
        bool isJoin = locationStr.Contains("join(", StringComparison.OrdinalIgnoreCase);

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

        return (overallStart, overallEnd, isComplement, isJoin, parts);
    }

    [GeneratedRegex(@"(\d+)(?:\.\.(\d+))?")]
    internal static partial Regex LocationRangeRegex();
}
