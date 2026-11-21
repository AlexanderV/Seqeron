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

/// <summary>
/// The line buffer of a feature file (<c>AjPFilebuff</c> of <c>ajfile.c</c>: <c>ajBuffreadLine</c>, <c>ajFilebuffClear</c>, <c>ajFilebuffReset</c>, <c>ajFilebuffSetUnbuffered</c>, <c>ajFilebuffIsEmpty</c>): a line read from the file is kept in the buffer until the
/// buffer is cleared (so that another reader can try the same lines after a reset); an UNBUFFERED buffer keeps nothing but the last line of a clear that keeps one. Whether the file was read to its end (<see cref="End"/>) decides what the
/// feature input does after a reader found no feature. A5-8c3e, B09 F241.
/// </summary>
internal sealed class EmbossFeatureLines(IReadOnlyList<string> lines)
{
    private readonly List<string> _buffer = [];
    private int _filePos;
    private int _pos;
    private string? _lastRead;

    /// <summary>The file was read up to its end (<c>File-&gt;End</c>: a read found no more lines).</summary>
    public bool End { get; private set; }

    /// <summary>The buffer retains no lines (<c>Nobuff</c>).</summary>
    public bool Unbuffered { get; private set; }

    /// <summary>The lines of a text: split at line feeds; a carriage return stays at the end of its line (<c>ajBuffreadLine</c> keeps the line end).</summary>
    public static EmbossFeatureLines FromText(string text)
    {
        var parts = new List<string>(text.Split('\n'));
        if (parts.Count > 0 && parts[^1].Length == 0)
            parts.RemoveAt(parts.Count - 1);
        return new EmbossFeatureLines(parts);
    }

    /// <summary><c>ajBuffreadLine</c>.</summary>
    public string? Next()
    {
        if (_pos < _buffer.Count)
            return _buffer[_pos++];
        if (_filePos >= lines.Count)
        {
            End = true;
            return null;
        }

        string l = lines[_filePos++];
        _lastRead = l;
        if (Unbuffered)
            return l;
        _buffer.Add(l);
        _pos++;
        return l;
    }

    /// <summary><c>ajFilebuffClear</c> with <paramref name="keep"/> lines: the lines before the last <paramref name="keep"/> read ones are dropped and the next read starts at the first kept one (the last line read is "given back").</summary>
    public void Clear(int keep)
    {
        int first = keep < 0 ? _buffer.Count : _pos - keep;
        if (first < 0)
            first = 0;
        if (Unbuffered && _pos == _buffer.Count)
            first = _pos;
        _buffer.RemoveRange(0, Math.Min(first, _buffer.Count));
        _pos = 0;
        if (Unbuffered && _buffer.Count == 0 && keep == 1 && _lastRead is not null)
            _buffer.Add(_lastRead);                 // unbuffered: only the last line read can be kept
    }

    /// <summary><c>ajFilebuffReset</c>: the next read starts at the first line of the buffer.</summary>
    public void Reset() => _pos = 0;

    /// <summary><c>ajFilebuffSetUnbuffered</c>.</summary>
    public void SetUnbuffered() => Unbuffered = true;

    /// <summary><c>ajFilebuffIsEmpty</c>: nothing is left in the buffer and the file was read to its end.</summary>
    public bool IsEmpty => _pos >= _buffer.Count && End;

    /// <summary>The file is at its end with nothing buffered (no more line can be read).</summary>
    public bool AtEnd => _pos >= _buffer.Count && _filePos >= lines.Count;
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
