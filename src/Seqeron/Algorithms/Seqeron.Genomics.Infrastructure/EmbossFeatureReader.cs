using System.Globalization;
using System.Numerics;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The string conversions of <c>ajstr.c</c> as the feature readers use them (<c>ajStrToInt</c>, <c>ajStrToUint</c>, <c>ajStrToFloat</c>, <c>ajStrToHex</c>: the whole string must be a number, else the value is 0 and the result false).</summary>
internal static class CStr
{
    private static bool Scan(string s, int radix, out BigInteger value)
    {
        value = BigInteger.Zero;
        int i = 0;
        while (i < s.Length && EmbossSeqFeatureTable.IsCSpace(s[i]))
            i++;
        bool neg = false;
        if (i < s.Length && s[i] is '+' or '-')
        {
            neg = s[i] == '-';
            i++;
        }

        if (radix == 16 && i + 1 < s.Length && s[i] == '0' && s[i + 1] is 'x' or 'X' && i + 2 < s.Length && Uri.IsHexDigit(s[i + 2]))
            i += 2;
        int start = i;
        BigInteger v = BigInteger.Zero;
        while (i < s.Length)
        {
            int d = s[i] switch { >= '0' and <= '9' => s[i] - '0', >= 'a' and <= 'f' => s[i] - 'a' + 10, >= 'A' and <= 'F' => s[i] - 'A' + 10, _ => 99 };
            if (d >= radix)
                break;
            v = v * radix + d;
            i++;
        }

        if (i == start || i != s.Length)
            return false;
        value = neg ? -v : v;
        return true;
    }

    /// <summary><c>ajStrToInt</c> (<c>strtol</c> base 10, clamped to the range of an <c>int</c>; a value beyond <c>long</c> fails).</summary>
    public static bool ToInt(string? s, out int value)
    {
        value = 0;
        if (string.IsNullOrEmpty(s) || !Scan(s, 10, out var v) || v > long.MaxValue || v < long.MinValue)
            return false;
        value = (int)Math.Clamp((long)v, int.MinValue, int.MaxValue);
        return true;
    }

    /// <summary><c>ajStrToHex</c>.</summary>
    public static bool ToHex(string? s, out int value)
    {
        value = 0;
        if (string.IsNullOrEmpty(s) || !Scan(s, 16, out var v) || v > long.MaxValue || v < long.MinValue)
            return false;
        value = (int)Math.Clamp((long)v, int.MinValue, int.MaxValue);
        return true;
    }

    /// <summary><c>ajStrToUint</c> (<c>strtoul</c>: a minus sign wraps).</summary>
    public static bool ToUint(string? s, out uint value)
    {
        value = 0;
        if (string.IsNullOrEmpty(s) || !Scan(s, 10, out var v))
            return false;
        BigInteger two64 = BigInteger.One << 64;
        if (v < 0)
            v += two64;
        if (v < 0 || v >= two64)
            return false;
        value = unchecked((uint)(ulong)v);
        return true;
    }

    /// <summary><c>ajStrToFloat</c> (<c>strtod</c> of the whole string; a value beyond the range of a <c>float</c> fails).</summary>
    public static bool ToFloat(string? s, out float value)
    {
        value = 0f;
        if (string.IsNullOrEmpty(s))
            return false;
        string t = s.TrimStart(' ', '\t', '\n', '\r', '\f', '\v');
        if (t.Length == 0)
            return false;
        string lower = t.ToLowerInvariant();
        string body = lower.TrimStart('+', '-');
        if (t.Length - body.Length > 1)
            return false;
        if (body == "nan" || (body.StartsWith("nan(", StringComparison.Ordinal) && body.EndsWith(')')))
        {
            value = lower.StartsWith('-') ? -float.NaN : float.NaN;
            return true;
        }

        if (body is "inf" or "infinity")
            return false;                                      // beyond FLT_MAX
        if (body.StartsWith("0x", StringComparison.Ordinal))
            return false;
        if (!Regex.IsMatch(t, @"^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$", RegexOptions.CultureInvariant))
            return false;
        if (!double.TryParse(t, NumberStyles.Float, CultureInfo.InvariantCulture, out double d))
            return false;
        if (double.IsInfinity(d))
            return false;                                      // ERANGE
        if (d > float.MaxValue || d < -float.MaxValue)
            return false;
        value = (float)d;
        return true;
    }
}

/// <summary>The line buffer of a feature file (<c>AjPFilebuff</c>): the lines are read in turn, the last one can be given back (<c>ajFilebuffClear</c> keeping one line), the whole buffer can be rewound (<c>ajFilebuffReset</c>).</summary>
internal sealed class EmbossFeatureLines(IReadOnlyList<string> lines)
{
    private int _pos;

    public int Position => _pos;

    public bool AtEnd => _pos >= lines.Count;

    /// <summary>The lines of a text: split at line feeds; a carriage return stays at the end of its line (<c>ajBuffreadLine</c> keeps the line end).</summary>
    public static EmbossFeatureLines FromText(string text)
    {
        var parts = new List<string>(text.Split('\n'));
        if (parts.Count > 0 && parts[^1].Length == 0)
            parts.RemoveAt(parts.Count - 1);
        return new EmbossFeatureLines(parts);
    }

    public string? Next() => _pos < lines.Count ? lines[_pos++] : null;

    public void Unread() => _pos--;

    public void Reset() => _pos = 0;

    public void Seek(int position) => _pos = position;
}

/// <summary>A static tokenizer of <c>ajstr.c</c> (<c>AjPStrTok</c>): the leading delimiters are skipped, a token ends at the next delimiter and ALL the delimiters after it are skipped.</summary>
internal sealed class StrTok
{
    private readonly string _s;
    private string _delim;
    private int _pos;

    public StrTok(string s, string delim)
    {
        _s = s;
        _delim = delim;
        _pos = Span(0);
    }

    private int Span(int from)
    {
        int p = from;
        while (p < _s.Length && _delim.Contains(_s[p], StringComparison.Ordinal))
            p++;
        return p;
    }

    /// <summary><c>ajStrTokenNextParse</c> (<paramref name="delim"/> null keeps the delimiters; <c>ajStrTokenNextParseC</c> sets them).</summary>
    public bool Next(out string token, string? delim = null)
    {
        if (delim is not null)
            _delim = delim;
        token = "";
        if (_pos >= _s.Length)
            return false;
        int q = _pos;
        int p = q;
        while (p < _s.Length && !_delim.Contains(_s[p], StringComparison.Ordinal))
            p++;
        token = _s[q..p];
        _pos = Span(p);
        return true;
    }

    /// <summary><c>ajStrTokenRestParse</c>.</summary>
    public bool Rest(out string rest)
    {
        rest = "";
        if (_pos >= _s.Length)
            return false;
        rest = _s[_pos..];
        _pos = _s.Length;
        return true;
    }
}
