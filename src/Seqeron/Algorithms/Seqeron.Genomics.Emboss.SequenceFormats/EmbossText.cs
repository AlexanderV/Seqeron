// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Analysis;

/// <summary>EMBOSS AJAX string conventions shared by the program ports of this assembly.</summary>
public static class EmbossText
{
    /// <summary>
    /// <c>ajStrAssignSubS(first, last)</c> / <c>ajStrAppendSubS</c>: inclusive positions; both ends are clipped into [0, length - 1] (<c>ajMathPosI</c>), the end not below the start, so a window that lies
    /// wholly past the end yields the last character and a reversed window (<c>last</c> &lt; <c>first</c>) one character at <c>first</c>.
    /// </summary>
    public static string Substring(string s, int first, int last)
    {
        if (s.Length == 0)
            return string.Empty;
        first = Math.Clamp(first, 0, s.Length - 1);
        last = Math.Clamp(last, first, s.Length - 1);
        return s.Substring(first, last - first + 1);
    }

    /// <summary>
    /// <c>ajStrAssignSubS</c> with <c>ajMathPosI</c>: inclusive 0-based positions, a negative position counts from the end (<c>length + position</c>), both are then clipped into [0, length - 1] and the end is not below the start
    /// (diffseq cuts its sequence at coordinates that may lie before the start: <c>-1</c>).
    /// </summary>
    public static string SubstringFromEnd(string s, long first, long last)
    {
        long len = s.Length;
        if (len == 0)
            return string.Empty;
        long b = first < 0 ? len + first : first;
        if (b >= len)
            b = len - 1;
        if (b < 0)
            b = 0;
        long e = last < 0 ? len + last : last;
        if (e >= len)
            e = len - 1;
        if (e < b)
            e = b;
        return s.Substring((int)b, (int)(e - b + 1));
    }
}
