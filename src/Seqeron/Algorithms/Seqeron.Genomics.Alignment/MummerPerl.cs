using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The Perl behaviours the ported MUMmer scripts (<c>mapview</c>, <c>mummerplot</c>, <c>nucmer2xfig</c>) depend on: numeric conversion of strings, number stringification, <c>split ' '</c>, line reading and <c>chomp</c>.</summary>
internal static class MummerPerl
{
    /// <summary>The numeric value of a string: the longest leading decimal number after white space (<c>12abc</c> is 12, <c>abc</c> is 0, <c>.5</c> is 0.5, <c>1e2</c> is 100); an undefined value is 0.</summary>
    public static double Num(string? s)
    {
        if (string.IsNullOrEmpty(s))
            return 0;
        int i = 0;
        while (i < s.Length && MummerFastaReader.IsSpace(s[i]))
            i++;
        int start = i;
        if (i < s.Length && (s[i] == '+' || s[i] == '-'))
            i++;
        int digits = 0;
        while (i < s.Length && s[i] >= '0' && s[i] <= '9')
        {
            i++;
            digits++;
        }

        if (i < s.Length && s[i] == '.')
        {
            i++;
            while (i < s.Length && s[i] >= '0' && s[i] <= '9')
            {
                i++;
                digits++;
            }
        }

        if (digits == 0)
        {
            string rest = s[start..].TrimStart('+', '-').ToLowerInvariant();
            if (rest.StartsWith("inf", StringComparison.Ordinal))
                return s[start] == '-' ? double.NegativeInfinity : double.PositiveInfinity;
            return rest.StartsWith("nan", StringComparison.Ordinal) ? double.NaN : 0;
        }

        if (i < s.Length && (s[i] == 'e' || s[i] == 'E'))
        {
            int j = i + 1;
            if (j < s.Length && (s[j] == '+' || s[j] == '-'))
                j++;
            int expStart = j;
            while (j < s.Length && s[j] >= '0' && s[j] <= '9')
                j++;
            if (j > expStart)
                i = j;
        }

        return double.Parse(s.AsSpan(start, i - start), NumberStyles.Float, CultureInfo.InvariantCulture);
    }

    /// <summary>The string of a number as <c>print</c> writes it (<c>%.15g</c>; a negative zero is <c>0</c>).</summary>
    public static string Str(double v)
    {
        if (v == 0)
            return "0";
        if (double.IsNaN(v))
            return "NaN";
        if (double.IsInfinity(v))
            return v > 0 ? "Inf" : "-Inf";
        return CFormat.General(v, 15);
    }

    /// <summary><c>int ($x)</c>: truncation towards zero.</summary>
    public static double Int(double v) => double.IsNaN(v) || double.IsInfinity(v) ? v : Math.Truncate(v);

    /// <summary><c>sprintf ("%.0f")</c>.</summary>
    public static string F0(double v) => CFormat.Fixed(v, 0);

    /// <summary><c>split</c> (<c>split ' ', $_</c>): fields separated by white space, leading white space ignored.</summary>
    public static List<string> Split(string? s)
    {
        var fields = new List<string>();
        if (s is null)
            return fields;
        int i = 0;
        while (i < s.Length)
        {
            while (i < s.Length && MummerFastaReader.IsSpace(s[i]))
                i++;
            int j = i;
            while (j < s.Length && !MummerFastaReader.IsSpace(s[j]))
                j++;
            if (j > i)
                fields.Add(s[i..j]);
            i = j;
        }

        return fields;
    }

    /// <summary><c>split (/regex/, $s)</c> with a one-character separator: empty leading fields are kept, trailing empty fields are removed.</summary>
    public static List<string> SplitChar(string? s, char separator)
    {
        var fields = new List<string>();
        if (string.IsNullOrEmpty(s))
            return fields;
        fields.AddRange(s.Split(separator));
        while (fields.Count > 0 && fields[^1].Length == 0)
            fields.RemoveAt(fields.Count - 1);
        return fields;
    }

    /// <summary><c>split (/\s+/, $s)</c>: like <c>split ' '</c> but an empty first field is kept when the string starts with white space; trailing empty fields are removed.</summary>
    public static List<string> SplitSpaces(string? s)
    {
        var fields = new List<string>();
        if (string.IsNullOrEmpty(s))
            return fields;
        fields.AddRange(System.Text.RegularExpressions.Regex.Split(s, @"\s+"));
        while (fields.Count > 0 && fields[^1].Length == 0)
            fields.RemoveAt(fields.Count - 1);
        return fields;
    }

    /// <summary>Is the string true in Perl (not empty and not <c>0</c>)?</summary>
    public static bool True(string? s) => !string.IsNullOrEmpty(s) && s != "0";

    /// <summary>The element <paramref name="index"/> of a Perl array (negative: from the end), null when it does not exist.</summary>
    public static string? At(IReadOnlyList<string> a, int index)
    {
        if (index < 0)
            index += a.Count;
        return index >= 0 && index < a.Count ? a[index] : null;
    }

    /// <summary>Sets the element <paramref name="index"/> of a Perl array, extending it with undefined (empty) elements.</summary>
    public static void Set(List<string> a, int index, string value)
    {
        while (a.Count <= index)
            a.Add(string.Empty);
        a[index] = value;
    }

    /// <summary>The lines of a text as <c>&lt;FH&gt;</c> reads them (each keeps its newline; the last may lack it).</summary>
    public sealed class Lines(string text)
    {
        private int _pos;

        public string? Next()
        {
            if (_pos >= text.Length)
                return null;
            int nl = text.IndexOf('\n', _pos);
            int end = nl < 0 ? text.Length : nl + 1;
            string line = text[_pos..end];
            _pos = end;
            return line;
        }
    }

    /// <summary><c>chomp</c>: removes one trailing newline.</summary>
    public static string Chomp(string s) => s.EndsWith('\n') ? s[..^1] : s;
}
