using System.Globalization;

namespace Seqeron.Genomics.Alignment;

/// <summary>The C library input routines the MUMmer programs read their text files with: <c>fgets</c> chunks and <c>sscanf</c> of <c>%ld</c>.</summary>
internal static class MummerCText
{
    /// <summary><c>fgets (Line, size, file)</c> until the end of the file: chunks of at most <c>size - 1</c> characters, each up to and including the newline.</summary>
    public static IEnumerable<string> Fgets(string text, int size)
    {
        int pos = 0;
        while (pos < text.Length)
        {
            int end = Math.Min(text.Length, pos + size - 1);
            int nl = text.IndexOf('\n', pos, end - pos);
            if (nl >= 0)
                end = nl + 1;
            yield return text.Substring(pos, end - pos);
            pos = end;
        }
    }

    /// <summary><c>sscanf (line, "%ld %ld ...")</c>: how many of the <paramref name="values"/> were converted (the others keep their value). White space (also newlines) before a number is skipped; the number is a sign and decimal digits; an overflow saturates as <c>strtol</c> does.</summary>
    public static int ScanLongs(string line, long[] values)
    {
        int pos = 0, count = 0;
        foreach (ref long slot in values.AsSpan())
        {
            while (pos < line.Length && MummerFastaReader.IsSpace(line[pos]))
                pos++;
            bool negative = false;
            if (pos < line.Length && (line[pos] == '+' || line[pos] == '-'))
            {
                negative = line[pos] == '-';
                pos++;
            }

            int digits = pos;
            System.Numerics.BigInteger value = 0;
            while (pos < line.Length && line[pos] >= '0' && line[pos] <= '9')
            {
                value = value * 10 + (line[pos] - '0');
                pos++;
            }

            if (pos == digits)
                return count;
            if (negative)
                value = -value;
            slot = Saturate(value);
            count++;
        }

        return count;
    }

    /// <summary><c>strtol (s, NULL, 10)</c> / <c>atoi</c> of the leading number of <paramref name="s"/> (0 when there is none).</summary>
    public static long StrToL(string s)
    {
        int p = 0;
        while (p < s.Length && MummerFastaReader.IsSpace(s[p]))
            p++;
        int st = p;
        if (p < s.Length && (s[p] == '+' || s[p] == '-'))
            p++;
        int d = p;
        while (p < s.Length && s[p] >= '0' && s[p] <= '9')
            p++;
        if (p == d)
            return 0;
        var v = System.Numerics.BigInteger.Parse(s.AsSpan(st, p - st), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture);
        return Saturate(v);
    }

    private static long Saturate(System.Numerics.BigInteger v)
    {
        if (v > long.MaxValue)
            return long.MaxValue;
        return v < long.MinValue ? long.MinValue : (long)v;
    }
}
