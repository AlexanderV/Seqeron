using System.Globalization;

namespace Seqeron.Genomics.Analysis;

/// <summary><c>strtod</c> as <c>ajStrToDouble</c> uses it: the whole text must be a number (decimal or hexadecimal integer / fraction, optional exponent) and in range; otherwise the value is 0 and the call fails.</summary>
internal static class EmbossStrtod
{
    public static bool TryParse(string text, out double value)
    {
        value = 0.0;
        if (text.Length == 0)
            return false;
        string t = text;
        bool neg = false;
        if (t[0] is '+' or '-')
        {
            neg = t[0] == '-';
            t = t[1..];
        }

        double d;
        if (t.Length > 2 && t[0] == '0' && t[1] is 'x' or 'X' && long.TryParse(t[2..], NumberStyles.AllowHexSpecifier, CultureInfo.InvariantCulture, out long hex))
        {
            d = hex;
        }
        else
        {
            if (t.Length == 0 || !(char.IsAsciiDigit(t[0]) || t[0] == '.'))
                return false;
            foreach (char c in t)
            {
                if (!(char.IsAsciiDigit(c) || c is '.' or 'e' or 'E' or '+' or '-'))
                    return false;
            }

            if (!double.TryParse(t, NumberStyles.Float, CultureInfo.InvariantCulture, out d) || double.IsInfinity(d))
                return false;
            if (d == 0.0 && t.Any(c => c is >= '1' and <= '9') && !t.Contains('e', StringComparison.OrdinalIgnoreCase))
                return false;
        }

        value = neg ? -d : d;
        return true;
    }
}

/// <summary>
/// A codon usage table (<c>AjPCod</c> of <c>ajcod.c</c>, 6.6.0) as <c>makenucseq</c> uses it: the frequency per thousand of each of the 64 codons (<c>tcount</c>). The reader is <c>codReadEmboss</c>
/// (the EMBOSS / <c>cut</c> format written by <c>cusp</c>: <c>Codon AA Fraction Frequency Number</c>, comments <c>#Species:</c> ...); a file in another format is not read (B09 WP62, F281).
/// </summary>
internal sealed class EmbossCodonUsage
{
    private readonly double[] _tcount = new double[64];
    private readonly int[] _num = new int[64];
    private readonly double[] _fraction = new double[64];

    private static int Index(string codon) => Base(codon[0]) * 16 + Base(codon[1]) * 4 + Base(codon[2]);

    private static int Base(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };

    /// <summary><c>ajCodRead</c> with the default (try every format) on a file of <paramref name="files"/>; null when the file cannot be opened or no reader accepts it.</summary>
    public static EmbossCodonUsage? Read(string name, IReadOnlyDictionary<string, string> files)
    {
        if (!files.TryGetValue(name, out string? text) && !files.TryGetValue("CODONS/" + name, out text))
            return null;
        var cod = new EmbossCodonUsage();
        return cod.ReadEmboss(text) ? cod : null;
    }

    private static bool IsFraction(string t)
    {
        int len = t.Length;
        if (len == 0)
            return false;
        int p = 0;
        bool startOne = false;
        if (t[0] == '0')
        {
            p++;
            len--;
        }
        else if (t[0] == '1')
        {
            p++;
            len--;
            startOne = true;
        }

        if (p < t.Length && t[p] == '.')
        {
            p++;
            len--;
        }

        if (p >= t.Length)
            return true;
        string rest = t[p..];
        return startOne ? rest.All(c => c == '0') && len == rest.Length : rest.All(char.IsAsciiDigit) && len == rest.Length;
    }

    private static bool IsFreq(string t)
    {
        if (t.Length == 0 || !t.All(c => char.IsAsciiDigit(c) || c == '.'))
            return false;
        if (!float.TryParse(t, NumberStyles.Float, CultureInfo.InvariantCulture, out float f))
            return false;
        return f >= 0.0f && f <= 1000.0f;
    }

    private bool ReadEmboss(string text)
    {
        foreach (string raw in EmbossSequencePrograms.ReadLines(text))
        {
            if (raw[0] == '\n' || raw[0] == '#' || raw[0] == '!')
                continue;
            string[] tok = raw.Split([' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
            if (tok.Length != 5)
                return false;
            string codon = tok[0].ToUpperInvariant().Replace('U', 'T');
            if (codon.Length != 3 || !codon.All(c => c is 'A' or 'C' or 'G' or 'T'))
                return false;
            int idx = Index(codon);
            string aa = tok[1].ToUpperInvariant();
            if (aa.Length != 1 || !"ACDEFGHIKLMNPQRSTUVWY*".Contains(aa[0], StringComparison.Ordinal))
                return false;
            if (!IsFraction(tok[2]) || !IsFreq(tok[3]))
                return false;
            if (tok[4].Length == 0 || !tok[4].All(char.IsAsciiDigit))
                return false;
            EmbossStrtod.TryParse(tok[2], out double fraction);
            EmbossStrtod.TryParse(tok[3], out double tcount);
            _ = int.TryParse(tok[4], NumberStyles.Integer, CultureInfo.InvariantCulture, out int num);
            _num[idx] = num;
            _tcount[idx] = tcount;
            _fraction[idx] = fraction;
        }

        // codFix: with no frequencies the frequencies per thousand are calculated from the numbers
        double totfreq = _tcount.Sum();
        if (Math.Abs(totfreq) <= 2.22e-15)
        {
            long totnum = 0;
            foreach (int n in _num)
                totnum += n;
            for (int i = 0; i < 64; i++)
                _tcount[i] = 1000.0 * _num[i] / totnum;
        }

        return true;
    }

    /// <summary><c>ajCodGetCodonlist</c>: every codon as many times as its frequency per thousand (rounded up), in the order AAA, AAC, ... TTT.</summary>
    public void AppendCodonList(List<string> list)
    {
        const string conv = "ACGT";
        for (int i = 0; i < 64; ++i)
        {
            for (int j = 0; j < _tcount[i]; j++)
                list.Add(new string([conv[i >> 4], conv[(i & 0xf) >> 2], conv[i & 3]]));
        }
    }
}
