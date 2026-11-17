using System.Globalization;
using System.Text;

namespace SuffixTree.EmProto.Verify;

/// <summary>The content of <c>.pk.meta</c>: length, width, alphabet, counts, SHA-256 of the input (written by PH0).</summary>
internal sealed record MetaInfo(long N, int Bits, int Sigma, byte[] Letters, long[] Counts, string Sha256);

/// <summary>The content of <c>.stats</c> (written after PH5).</summary>
internal sealed record StatsInfo(long NodeCount, long Internal, long LrsDepth, long LrsRank, long N, int Sigma, int LcpWidth, long Overflow, long MaxLcp);

/// <summary>
/// Strict parsers of the two small text files of an <c>e6</c> run. Strict means one exact form (keys in the order the pipeline writes them, canonical decimal numbers, upper-case hex for the letters, lower-case hex for the hash,
/// every line ended by a line feed): a flipped byte anywhere in the file is a parse error or a different value, never a silent pass.
/// </summary>
internal static class MetaFiles
{
    internal static bool TryParseMeta(string path, out MetaInfo meta, out string error)
    {
        meta = null!;
        if (!TryLines(path, ["n", "bits", "sigma", "letters_hex", "counts", "sha256", "name"], out string[] v, out error)) return false;
        if (!Num(v[0], out long n) || !Num(v[1], out long bits) || !Num(v[2], out long sigma))
        {
            error = "n, bits or sigma is not a canonical decimal number";
            return false;
        }
        if (sigma is < 1 or > 256 || bits is < 1 or > 8 || v[3].Length != 2 * sigma || !v[3].All(c => c is >= '0' and <= '9' or >= 'A' and <= 'F'))
        {
            error = "sigma, bits or letters_hex out of range";
            return false;
        }
        string[] parts = v[4].Split(',');
        var counts = new long[parts.Length];
        for (int i = 0; i < parts.Length; i++)
        {
            if (!Num(parts[i], out counts[i]))
            {
                error = "counts is not a list of canonical decimal numbers";
                return false;
            }
        }
        if (counts.Length != sigma || v[5].Length != 64 || !v[5].All(c => c is >= '0' and <= '9' or >= 'a' and <= 'f'))
        {
            error = "counts or sha256 malformed";
            return false;
        }
        meta = new MetaInfo(n, (int)bits, (int)sigma, Convert.FromHexString(v[3]), counts, v[5]);
        return true;
    }

    internal static bool TryParseStats(string path, out StatsInfo stats, out string error)
    {
        stats = null!;
        if (!TryLines(path, ["nodecount", "internal", "lrsdepth", "lrsrank", "n", "sigma", "lcpwidth", "overflow", "maxlcp"], out string[] v, out error)) return false;
        var x = new long[v.Length];
        for (int i = 0; i < v.Length; i++)
        {
            if (!Num(v[i], out x[i]))
            {
                error = $"line {i + 1} is not a canonical decimal number";
                return false;
            }
        }
        stats = new StatsInfo(x[0], x[1], x[2], x[3], x[4], (int)Math.Min(x[5], int.MaxValue), (int)Math.Min(x[6], int.MaxValue), x[7], x[8]);
        return true;
    }

    /// <summary>The statistics file as the pipeline writes it, from recomputed values.</summary>
    internal static string StatsText(StatsInfo s) => string.Create(CultureInfo.InvariantCulture,
        $"nodecount={s.NodeCount}\ninternal={s.Internal}\nlrsdepth={s.LrsDepth}\nlrsrank={s.LrsRank}\nn={s.N}\nsigma={s.Sigma}\nlcpwidth={s.LcpWidth}\noverflow={s.Overflow}\nmaxlcp={s.MaxLcp}\n");

    internal static string ReadText(string path) => Encoding.Latin1.GetString(File.ReadAllBytes(path));

    /// <summary>The first differing line of two statistics texts, for the message of a failed check.</summary>
    internal static string FirstDifference(string expected, string found)
    {
        string[] e = expected.Split('\n'), f = found.Split('\n');
        for (int i = 0; i < Math.Max(e.Length, f.Length); i++)
        {
            string a = i < e.Length ? e[i] : "<none>", b = i < f.Length ? f[i] : "<none>";
            if (a != b) return $"line {i + 1}: expected '{a}', found '{b}'";
        }
        return "no difference";
    }

    private static bool TryLines(string path, string[] keys, out string[] values, out string error)
    {
        values = new string[keys.Length];
        error = string.Empty;
        string text = ReadText(path);
        if (!text.EndsWith('\n'))
        {
            error = "the file does not end with a line feed";
            return false;
        }
        string[] lines = text[..^1].Split('\n');
        if (lines.Length != keys.Length)
        {
            error = $"{lines.Length} lines, expected {keys.Length}";
            return false;
        }
        for (int i = 0; i < keys.Length; i++)
        {
            string prefix = keys[i] + "=";
            if (!lines[i].StartsWith(prefix, StringComparison.Ordinal))
            {
                error = $"line {i + 1} does not start with '{prefix}'";
                return false;
            }
            values[i] = lines[i][prefix.Length..];
        }
        return true;
    }

    /// <summary>A canonical non-negative decimal number: digits only, no sign, no white space, no leading zero, at most 18 digits.</summary>
    private static bool Num(string s, out long value)
    {
        value = 0;
        if (s.Length is 0 or > 18 || (s.Length > 1 && s[0] == '0') || !s.All(c => c is >= '0' and <= '9')) return false;
        value = long.Parse(s, NumberStyles.None, CultureInfo.InvariantCulture);
        return true;
    }
}
