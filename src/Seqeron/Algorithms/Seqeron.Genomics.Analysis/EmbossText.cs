namespace Seqeron.Genomics.Analysis;

/// <summary>EMBOSS AJAX string conventions shared by the program ports of this assembly.</summary>
internal static class EmbossText
{
    /// <summary>
    /// <c>ajStrAssignSubS(first, last)</c> / <c>ajStrAppendSubS</c>: inclusive positions; both ends are clipped into [0, length - 1] (<c>ajMathPosI</c>), the end not below the start, so a window that lies
    /// wholly past the end yields the last character and a reversed window (<c>last</c> &lt; <c>first</c>) one character at <c>first</c>.
    /// </summary>
    internal static string Substring(string s, int first, int last)
    {
        if (s.Length == 0)
            return string.Empty;
        first = Math.Clamp(first, 0, s.Length - 1);
        last = Math.Clamp(last, first, s.Length - 1);
        return s.Substring(first, last - first + 1);
    }
}
