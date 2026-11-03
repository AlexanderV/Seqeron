using System.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The generator families of the differential fuzz harness (work package 1a.7): seeded, deterministic, every family usable at every length
/// 0..10^5 and every alphabet size 1..1000 (symbols from <see cref="StringFamilies.Alphabet"/>, so alphabets above 26 symbols use code units
/// above 255 and alphabets above 256 symbols switch SA-IS to the ushort reader). All families use <see cref="SplitMix64"/> only.
/// </summary>
internal static class FuzzTexts
{
    /// <summary>The family names, in the order used by <see cref="Pick"/>.</summary>
    internal static readonly string[] Families =
    [
        "uniform", "runs", "periodic-noise", "block-copy", "fibonacci", "thue-morse", "debruijn", "skewed", "akb", "splice",
        "abab-noise", "permutation", "unicode-mix",
    ];

    private static readonly int[] Sigmas = [1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 8, 20, 26, 27, 64, 255, 256, 257, 300, 1000];
    private static readonly int[] RunLengths = [2, 3, 8, 64, 1000];

    // ascending code units: NUL, 'a', BMP above 255, the 15-bit boundary, high and low surrogates (lone), U+FFFF
    private static readonly char[] UnicodeUnits = ['\0', 'a', 'Ā', '翿', '耀', '\uD800', '\uDBFF', '\uDC00', '\uDFFF', '￿'];

    /// <summary>Draws the family, alphabet size and parameter of a case of length <paramref name="n"/> from <paramref name="rng"/>.</summary>
    internal static FuzzCase Pick(SplitMix64 rng, int n)
    {
        string family = Families[rng.NextInt(Families.Length)];
        if (family == "permutation" && n > 65_536) family = "uniform";
        int sigma = Sigmas[rng.NextInt(Sigmas.Length)];
        int param = 0;
        switch (family)
        {
            case "runs": param = RunLengths[rng.NextInt(RunLengths.Length)]; break;
            case "periodic-noise": param = 1 + rng.NextInt(70); break;
            case "block-copy": param = 1 + rng.NextInt(Math.Max(1, n / 2)); break;
            case "fibonacci" or "thue-morse": param = rng.NextInt(5000); break;
            case "akb": param = 1 + rng.NextInt(100); break;
            case "abab-noise": param = 2 + rng.NextInt(30); break;
            case "permutation": sigma = n; break;
            case "unicode-mix": sigma = 1 + rng.NextInt(UnicodeUnits.Length); break;
        }
        return new FuzzCase(family, rng.Next(), n, sigma, param);
    }

    /// <summary>The text of <paramref name="c"/> (exactly <c>c.N</c> code units).</summary>
    internal static string Generate(FuzzCase c)
    {
        int n = c.N;
        if (n == 0) return string.Empty;
        var rng = new SplitMix64(c.Seed);
        char[] a = c.Family switch
        {
            "unicode-mix" => UnicodeSubset(rng, c.Sigma),
            "permutation" => [],
            _ => StringFamilies.Alphabet(c.Sigma),
        };
        string text = c.Family switch
        {
            "uniform" => Fill(n, () => a[rng.NextInt(a.Length)]),
            "unicode-mix" => Fill(n, () => a[rng.NextInt(a.Length)]),
            "runs" => Runs(n, c.Param, a, rng),
            "periodic-noise" => PeriodicNoise(n, c.Param, a, rng),
            "block-copy" => BlockCopy(n, c.Param, a, rng),
            "fibonacci" => Map(StringFamilies.Fibonacci(n + c.Param).Substring(c.Param, n), a),
            "thue-morse" => Map(StringFamilies.ThueMorse(n + c.Param).Substring(c.Param, n), a),
            "debruijn" => DeBruijn(n, a),
            "skewed" => Fill(n, () => Skewed(a, rng)),
            "akb" => Map(StringFamilies.AkB(n, c.Param), a),
            "splice" => Splice(n, a, rng),
            "abab-noise" => AbabNoise(n, c.Param, a, rng),
            "permutation" => StringFamilies.AllDistinct(n, c.Seed),
            _ => throw new ArgumentException($"unknown fuzz family '{c.Family}'", nameof(c)),
        };
        if (text.Length != n) throw new InvalidOperationException($"family {c.Family} produced {text.Length} units, expected {n}");
        return text;
    }

    private static string Fill(int n, Func<char> next)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = next();
        return new string(c);
    }

    /// <summary>Maps the letters 'a', 'b', ... of a family text to the symbols of the alphabet (letter index modulo the alphabet size).</summary>
    private static string Map(string text, char[] a)
    {
        var c = new char[text.Length];
        for (int i = 0; i < c.Length; i++) c[i] = a[(text[i] - 'a') % a.Length];
        return new string(c);
    }

    private static char[] UnicodeSubset(SplitMix64 rng, int sigma)
    {
        var idx = Enumerable.Range(0, UnicodeUnits.Length).ToArray();
        for (int i = 0; i < sigma; i++)
        {
            int j = i + rng.NextInt(idx.Length - i);
            (idx[i], idx[j]) = (idx[j], idx[i]);
        }
        return idx.Take(sigma).OrderBy(x => x).Select(x => UnicodeUnits[x]).ToArray();
    }

    private static string Runs(int n, int maxRun, char[] a, SplitMix64 rng)
    {
        var sb = new StringBuilder(n);
        while (sb.Length < n)
        {
            int len = 1 + rng.NextInt(maxRun);
            sb.Append(a[rng.NextInt(a.Length)], Math.Min(len, n - sb.Length));
        }
        return sb.ToString();
    }

    private static string PeriodicNoise(int n, int period, char[] a, SplitMix64 rng)
    {
        int p = Math.Clamp(period, 1, n);
        var pattern = new char[p];
        for (int i = 0; i < p; i++) pattern[i] = a[rng.NextInt(a.Length)];
        int noise = rng.NextInt(4) == 0 ? 0 : 4 + rng.NextInt(60);
        var c = new char[n];
        for (int i = 0; i < n; i++)
        {
            c[i] = pattern[i % p];
            if (noise != 0 && rng.NextInt(noise) == 0) c[i] = a[rng.NextInt(a.Length)];
        }
        return new string(c);
    }

    private static string BlockCopy(int n, int blockLength, char[] a, SplitMix64 rng)
    {
        int b = Math.Clamp(blockLength, 1, n);
        var block = new char[b];
        for (int i = 0; i < b; i++) block[i] = a[rng.NextInt(a.Length)];
        var c = new char[n];
        for (int i = 0; i < n; i++)
        {
            c[i] = block[i % b];
            if (rng.NextInt(32) == 0) c[i] = a[rng.NextInt(a.Length)]; // about 3 % substitutions, as in the rep_10m benchmark text
        }
        return new string(c);
    }

    private static string DeBruijn(int n, char[] a)
    {
        int k = Math.Min(a.Length, 4);
        if (k < 2) return new string(a[0], n);
        int order = 1;
        for (long size = k; size < n; size *= k) order++;
        return Map(StringFamilies.DeBruijn(k, order, n), a);
    }

    private static char Skewed(char[] a, SplitMix64 rng)
    {
        int i = 0;
        while (i < a.Length - 1 && (rng.Next() & 1) == 0) i++; // P(symbol i) = 2^-(i+1)
        return a[i];
    }

    private static string Splice(int n, char[] a, SplitMix64 rng)
    {
        int parts = 2 + rng.NextInt(3);
        var sb = new StringBuilder(n);
        for (int s = 0; s < parts; s++)
        {
            int len = s == parts - 1 ? n - sb.Length : rng.NextInt(n - sb.Length + 1);
            switch (rng.NextInt(4))
            {
                case 0: // uniform over at most 3 symbols
                    for (int i = 0; i < len; i++) sb.Append(a[rng.NextInt(Math.Min(3, a.Length))]);
                    break;
                case 1: // one symbol
                    sb.Append(a[rng.NextInt(a.Length)], len);
                    break;
                case 2: // short period
                {
                    int p = 1 + rng.NextInt(8);
                    var pattern = new char[p];
                    for (int i = 0; i < p; i++) pattern[i] = a[rng.NextInt(a.Length)];
                    for (int i = 0; i < len; i++) sb.Append(pattern[i % p]);
                    break;
                }
                default: // runs
                    sb.Append(Runs(len, 1 + rng.NextInt(16), a, rng));
                    break;
            }
        }
        return sb.ToString();
    }

    private static string AbabNoise(int n, int noiseDenominator, char[] a, SplitMix64 rng)
    {
        char x = a[0], y = a[Math.Min(1, a.Length - 1)];
        var c = new char[n];
        for (int i = 0; i < n; i++)
        {
            c[i] = (i & 1) == 0 ? x : y;
            if (rng.NextInt(noiseDenominator) == 0) c[i] = a[rng.NextInt(a.Length)];
        }
        return new string(c);
    }
}
