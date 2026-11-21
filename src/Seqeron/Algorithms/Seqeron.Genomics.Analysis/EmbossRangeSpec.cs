using System.Globalization;

namespace Seqeron.Genomics.Analysis;

/// <summary>One element of an EMBOSS range specification: 1-based inclusive positions and the text that follows the pair of numbers (a colour for <c>showalign -highlight</c>).</summary>
public readonly record struct EmbossRangeElement(int Start, int End, string Text);

/// <summary>
/// The EMBOSS 6.6.0 range specification of the ACD type <c>range</c> (<c>ajrange.c</c>: <c>ajRangeNewStringLimits</c>, <c>ajRangeNewFilenameLimits</c>): pairs of positions separated by any non-digit
/// characters ("24-45, 56-78", "1:45 67=99;765..888", "3+4" = three to six), each pair optionally followed by a text ("24-45 blue 56-78 orange"), or <c>@filename</c>: one range per line (a line
/// that starts with <c>#</c> is a comment). The parser keeps the program's <c>strtok</c> behaviour: the text of a pair is the non-digit run after the second number, so "3-5,7-9" gives the text "," for
/// the first pair and none for the second.
/// </summary>
public static class EmbossRangeSpec
{
    private const string NonDigit = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ \t\n\r!@#$%^&*()_-=|\\~`{[}]:;\"'<,>.?/+";
    private const string NonDigitPlus = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ \t\n\r!@#$%^&*()_-=|\\~`{[}]:;\"'<,>.?/";
    private const string Digit = "0123456789";

    /// <summary>
    /// Parses a range specification; an empty text gives no ranges. Throws <see cref="ArgumentException"/> where the program refuses the value ("Odd integer(s) in range specification",
    /// "Bad range value", "From range greater than To range", a position below <paramref name="min"/>, a missing file).
    /// </summary>
    public static IReadOnlyList<EmbossRangeElement> Parse(string? spec, int min = 1, int max = int.MaxValue)
    {
        string s = Trim(spec ?? string.Empty);
        if (s.Length > 0 && s[0] == '@')
            return ParseFile(s[1..], min, max);

        var result = new List<EmbossRangeElement>();
        var count = Tokens(s, NonDigit).Count;
        if (count == 0)
            return result;
        if (count % 2 != 0)
            throw new ArgumentException($"Odd integer(s) in range specification [{count}]", nameof(spec));

        int e = count >> 1;
        var tokens = Tokens(s, NonDigitPlus);
        int next = 0;
        var starts = new long[e];
        var ends = new long[e];
        for (int i = 0; i < e; i++)
        {
            string p = Take(tokens, ref next, spec);
            long f, t;
            if (TryUnsignedPlusUnsigned(p, out f, out long second))
            {
                t = second != 0 ? f + second - 1 : f;
            }
            else
            {
                if (!TryLeadingUnsigned(p, out f))
                    throw new ArgumentException($"Bad range value [{p}]", nameof(spec));
                p = Take(tokens, ref next, spec);
                if (p.StartsWith('+') && TryLeadingUnsigned(p[1..], out long plus))
                    t = plus != 0 ? f + plus - 1 : f;
                else if (!TryLeadingUnsigned(p, out t))
                    throw new ArgumentException($"Bad range value [{p}]", nameof(spec));
            }

            Check(f, t, min, max);
            starts[i] = f;
            ends[i] = t;
        }

        // the texts: the non-digit runs of the string (ignoring what lies between the two numbers of a pair)
        var runs = Tokens(s, Digit);
        int r = s.Length > 0 && !char.IsAsciiDigit(s[0]) ? 1 : 0;   // a leading non-digit run is read first (doneone)
        var text = new string[e];
        for (int i = 0; i < e; i++)
        {
            // runs[r] lies between the two numbers of the pair (ignored), runs[r + 1] follows them: the text
            text[i] = r + 1 < runs.Count ? Trim(runs[r + 1]) : string.Empty;
            r += 2;
        }

        for (int i = 0; i < e; i++)
            result.Add(new EmbossRangeElement((int)starts[i], (int)ends[i], text[i]));
        return result;
    }

    private static IReadOnlyList<EmbossRangeElement> ParseFile(string name, int min, int max)
    {
        name = Trim(name);
        if (!File.Exists(name))
            throw new ArgumentException($"Unable to read range file '{name}'", nameof(name));
        var result = new List<EmbossRangeElement>();
        foreach (var raw in File.ReadAllLines(name))
        {
            string line = Trim(raw);
            if (line.StartsWith('#') || line.Length == 0)
                continue;
            const string white = "- \t\n\r";
            int i = 0;
            string one = NextToken(line, white, ref i);
            string two = NextToken(line, white, ref i);
            if (two.Length == 0)
                throw new ArgumentException($"Odd integer(s) in range specification: {line}", nameof(name));
            string rest = Trim(line[Math.Min(i, line.Length)..].TrimEnd('\n', '\r'));
            if (!uint.TryParse(one, NumberStyles.None, CultureInfo.InvariantCulture, out uint a))
                throw new ArgumentException($"Bad range value [{one}]", nameof(name));
            if (!uint.TryParse(two, NumberStyles.None, CultureInfo.InvariantCulture, out uint b))
                throw new ArgumentException($"Bad range value [{two}]", nameof(name));
            Check(a, b, min, max);
            result.Add(new EmbossRangeElement((int)a, (int)b, rest));
        }

        return result;
    }

    private static string NextToken(string line, string delimiters, ref int i)
    {
        while (i < line.Length && delimiters.Contains(line[i]))
            i++;
        int start = i;
        while (i < line.Length && !delimiters.Contains(line[i]))
            i++;
        return line[start..i];
    }

    private static void Check(long f, long t, int min, int max)
    {
        if (f > t)
            throw new ArgumentException($"From range [{f}] greater than To range [{t}]");
        if (f < min)
            throw new ArgumentException($"From range [{f}] less than minimum [{min}]");
        if (t > max)
            throw new ArgumentException($"To range [{t}] greater than maximum [{max}]");
    }

    private static string Take(List<string> tokens, ref int next, string? spec)
    {
        if (next >= tokens.Count)
            throw new ArgumentException($"Bad range value [] in '{spec}'", nameof(spec));
        return tokens[next++];
    }

    // strtok(string, delimiters): the maximal runs of characters that are not delimiters
    private static List<string> Tokens(string s, string delimiters)
    {
        var tokens = new List<string>();
        int i = 0;
        while (i < s.Length)
        {
            while (i < s.Length && delimiters.Contains(s[i]))
                i++;
            if (i >= s.Length)
                break;
            int start = i;
            while (i < s.Length && !delimiters.Contains(s[i]))
                i++;
            tokens.Add(s[start..i]);
        }

        return tokens;
    }

    // sscanf(p, "%u+%u"): both numbers read
    private static bool TryUnsignedPlusUnsigned(string p, out long first, out long second)
    {
        second = 0;
        int i = 0;
        if (!ReadUnsigned(p, ref i, out first) || i >= p.Length || p[i] != '+')
            return false;
        i++;
        return ReadUnsigned(p, ref i, out second);
    }

    // sscanf(p, "%u"): the leading unsigned number (a leading '+' or '-' sign is accepted by %u; none occurs in the tokens)
    private static bool TryLeadingUnsigned(string p, out long value)
    {
        int i = 0;
        return ReadUnsigned(p, ref i, out value);
    }

    private static bool ReadUnsigned(string p, ref int i, out long value)
    {
        value = 0;
        int start = i;
        while (i < p.Length && char.IsAsciiDigit(p[i]))
        {
            value = value * 10 + (p[i] - '0');
            if (value > uint.MaxValue)
                value &= uint.MaxValue;
            i++;
        }

        return i > start;
    }

    private static string Trim(string s) => s.Trim(' ', '\t', '\n', '\r', '\f', '\v');
}
