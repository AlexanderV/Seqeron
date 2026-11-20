using System.Diagnostics;
using System.Globalization;
using SuffixTree;
using SuffixTree.Esa;

namespace ReadBench;

/// <summary>
/// One process = one tree. Variants: memory (EsaSuffixTreeFactory.Create over the text), checked / trusted (EsaSuffixTreeFactory.Open at Header level), resident (when the build has the option), exp:* (the experiments of the diagnosis, built with -p:Experiments=true).
/// Output: one TSV row per process: variant, count20/s, locate8/s, absent24/s (medians of the passes), open seconds, and a few /proc counters.
/// usage: ReadBench &lt;variant&gt; &lt;container&gt; &lt;text&gt; [passes=7] [label]
/// </summary>
internal static class Program
{
    private static int Main(string[] args)
    {
        if (args.Length < 3)
        {
            Console.Error.WriteLine("usage: ReadBench <variant> <container> <text> [passes=7] [label] [alphabet=ACGT] [seed=1234]");
            return 2;
        }
        string variant = args[0], file = args[1], textPath = args[2];
        int passes = args.Length > 3 ? int.Parse(args[3], CultureInfo.InvariantCulture) : 7;
        string label = args.Length > 4 ? args[4] : "run";
        string alphabet = args.Length > 5 ? args[5] : Environment.GetEnvironmentVariable("RB_ALPHABET") ?? "ACGT";
        ulong seed = args.Length > 6 ? ulong.Parse(args[6], CultureInfo.InvariantCulture) : 1234;
        string text = File.ReadAllText(textPath, System.Text.Encoding.Latin1);

        if (variant == "open") return OpenBench(file);
#if EXPERIMENTS
        if (variant == "count") return ExperimentsCount.Run(file, textPath);
#endif
        if (variant.Contains('+', StringComparison.Ordinal)) return PairBench(variant.Split('+'), file, text, passes, label, alphabet, seed);

        long t0 = Stopwatch.GetTimestamp();
        Q tree = CreateQ(variant, file, text);
        double openSeconds = Seconds(t0);

        (string[] twenty, string[] eight, string[] absent) = Patterns(text, alphabet, seed);

        // warm-up by time: tier-1 code (and its PGO) is installed after a call-count delay, so a fixed number of calls is not enough
        WarmFor(1.5, () => RunCount(tree, twenty));
        WarmFor(0.7, () => RunLocate(tree, eight));
        WarmFor(1.0, () => RunCount(tree, absent));

        var c20 = new List<double>();
        var l8 = new List<double>();
        var a24 = new List<double>();
        long sink = 0;
        for (int r = 0; r < passes; r++)
        {
            long t = Stopwatch.GetTimestamp();
            sink += RunCount(tree, twenty);
            c20.Add(twenty.Length / Seconds(t));
            t = Stopwatch.GetTimestamp();
            sink += RunLocate(tree, eight);
            l8.Add(eight.Length / Seconds(t));
            t = Stopwatch.GetTimestamp();
            sink += RunCount(tree, absent);
            a24.Add(absent.Length / Seconds(t));
        }
        GC.KeepAlive(sink);
        string smaps = Smaps();
        Console.WriteLine(string.Join('\t', "row", label, variant + (Environment.GetEnvironmentVariable("RB_TAG") ?? ""),
            Median(c20).ToString("F0", CultureInfo.InvariantCulture), Median(l8).ToString("F0", CultureInfo.InvariantCulture), Median(a24).ToString("F0", CultureInfo.InvariantCulture),
            openSeconds.ToString("F4", CultureInfo.InvariantCulture), sink, smaps));
        return 0;
    }

    /// <summary>The queries of a pass, behind one small abstraction so that a public tree and an experiment engine are timed by the same loops.</summary>
    internal abstract class Q
    {
        internal abstract int Count(string pattern);

        internal abstract int Locate(string pattern);
    }

    private sealed class TreeQ(ISuffixTree tree) : Q
    {
        internal override int Count(string pattern) => tree.CountOccurrences(pattern);

        internal override int Locate(string pattern) => tree.FindAllOccurrences(pattern).Count;
    }

    private static Q CreateQ(string variant, string file, string text)
    {
#if EXPERIMENTS
        if (variant.StartsWith("ix-", StringComparison.Ordinal)) return ExperimentsIx.Create(variant, file, text);
#endif
        return new TreeQ(Create(variant, file, text));
    }

    private static (string[] Twenty, string[] Eight, string[] Absent) Patterns(string text, string alphabet, ulong seed)
    {
        var rng = new Mix(seed);
        int m20 = int.Parse(Environment.GetEnvironmentVariable("RB_N20") ?? "200000", CultureInfo.InvariantCulture);
        int m8 = int.Parse(Environment.GetEnvironmentVariable("RB_N8") ?? "20000", CultureInfo.InvariantCulture);
        int m24 = int.Parse(Environment.GetEnvironmentVariable("RB_N24") ?? "100000", CultureInfo.InvariantCulture);
        string[] twenty = new string[m20];
        for (int i = 0; i < twenty.Length; i++) twenty[i] = text.Substring(rng.NextInt(text.Length - 20), 20);
        string[] eight = new string[m8];
        for (int i = 0; i < eight.Length; i++) eight[i] = text.Substring(rng.NextInt(text.Length - 8), 8);
        string[] absent = new string[m24];
        for (int i = 0; i < absent.Length; i++)
        {
            var chars = new char[24];
            for (int k = 0; k < chars.Length; k++) chars[k] = alphabet[rng.NextInt(alphabet.Length)];
            absent[i] = new string(chars);
        }
        return (twenty, eight, absent);
    }

    /// <summary>
    /// In-process A/B: several trees in one process, every pass runs all of them in a rotated order (pass p starts with tree p mod k), so that the drift of the machine and the order are not a bias; the output is, per tree, the median
    /// throughput and the median over the passes of its ratio to the first tree.
    /// </summary>
    private static int PairBench(string[] variants, string file, string text, int passes, string label, string alphabet, ulong seed)
    {
        var qs = new Q[variants.Length];
        for (int i = 0; i < qs.Length; i++) qs[i] = CreateQ(variants[i], file, text);
        (string[] twenty, string[] eight, string[] absent) = Patterns(text, alphabet, seed);
        for (int i = 0; i < qs.Length; i++)
        {
            Q q = qs[i];
            WarmFor(1.5, () => RunCount(q, twenty));
            WarmFor(0.7, () => RunLocate(q, eight));
            WarmFor(1.0, () => RunCount(q, absent));
        }
        int k = qs.Length;
        var c = new List<double>[k];
        var l = new List<double>[k];
        var a = new List<double>[k];
        for (int i = 0; i < k; i++) { c[i] = []; l[i] = []; a[i] = []; }
        var sinks = new long[k];
        for (int r = 0; r < passes; r++)
        {
            for (int j = 0; j < k; j++)
            {
                int i = (r + j) % k;
                Q q = qs[i];
                long t = Stopwatch.GetTimestamp();
                sinks[i] += RunCount(q, twenty);
                c[i].Add(twenty.Length / Seconds(t));
                t = Stopwatch.GetTimestamp();
                sinks[i] += RunLocate(q, eight);
                l[i].Add(eight.Length / Seconds(t));
                t = Stopwatch.GetTimestamp();
                sinks[i] += RunCount(q, absent);
                a[i].Add(absent.Length / Seconds(t));
            }
        }
        for (int i = 0; i < k; i++)
        {
            string Ratio(List<double>[] m) => Median(m[i].Select((x, p) => x / m[0][p]).ToList()).ToString("F3", CultureInfo.InvariantCulture);
            Console.WriteLine(string.Join('\t', "pair", label, variants[i],
                Median(c[i]).ToString("F0", CultureInfo.InvariantCulture), Median(l[i]).ToString("F0", CultureInfo.InvariantCulture), Median(a[i]).ToString("F0", CultureInfo.InvariantCulture),
                "ratio_to_first:", Ratio(c), Ratio(l), Ratio(a), "sink=" + sinks[i]));
        }
        return 0;
    }

    private static ISuffixTree Create(string variant, string file, string text)
    {
        switch (variant)
        {
            case "memory": return EsaSuffixTreeFactory.Create(text);
#if !NOOPEN
            case "checked": return EsaSuffixTreeFactory.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Checked });
            case "trusted": return EsaSuffixTreeFactory.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Trusted });
#endif
            default:
#if EXPERIMENTS
                return Experiments.Create(variant, file, text);
#else
                throw new ArgumentException($"unknown variant '{variant}'");
#endif
        }
    }

    private static int OpenBench(string file)
    {
#if NOOPEN
        throw new NotSupportedException("this build has no Open");
#else
        // JIT-warm open at the Header level: three opens to warm, then the median of 25
        for (int i = 0; i < 5; i++) EsaSuffixTreeFactory.Open(file).Dispose();
        var t = new List<double>();
        for (int i = 0; i < 25; i++)
        {
            long t0 = Stopwatch.GetTimestamp();
            EsaSuffixTree tree = EsaSuffixTreeFactory.Open(file);
            t.Add(Seconds(t0));
            tree.Dispose();
        }
        Console.WriteLine(string.Join('\t', "open", "header", Median(t).ToString("F6", CultureInfo.InvariantCulture)));
        return 0;
#endif
    }

    private static long RunCount(ISuffixTree tree, string[] patterns)
    {
        long sum = 0;
        for (int i = 0; i < patterns.Length; i++) sum += tree.CountOccurrences(patterns[i]);
        return sum;
    }

    private static long RunLocate(ISuffixTree tree, string[] patterns)
    {
        long sum = 0;
        for (int i = 0; i < patterns.Length; i++) sum += tree.FindAllOccurrences(patterns[i]).Count;
        return sum;
    }

    private static long RunCount(Q q, string[] patterns)
    {
        long sum = 0;
        for (int i = 0; i < patterns.Length; i++) sum += q.Count(patterns[i]);
        return sum;
    }

    private static long RunLocate(Q q, string[] patterns)
    {
        long sum = 0;
        for (int i = 0; i < patterns.Length; i++) sum += q.Locate(patterns[i]);
        return sum;
    }

    private static void WarmFor(double seconds, Func<long> pass)
    {
        long t0 = Stopwatch.GetTimestamp();
        do { _ = pass(); } while (Seconds(t0) < seconds);
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    private static double Median(List<double> v)
    {
        double[] s = [.. v.Order()];
        return s[s.Length / 2];
    }

    private static string Smaps()
    {
        try
        {
            string rollup = File.ReadAllText("/proc/self/smaps_rollup");
            string Field(string name)
            {
                foreach (string line in rollup.Split('\n'))
                    if (line.StartsWith(name + ":", StringComparison.Ordinal)) return line[(name.Length + 1)..].Replace("kB", "", StringComparison.Ordinal).Trim();
                return "-";
            }
            string status = File.ReadAllText("/proc/self/status");
            string St(string name)
            {
                foreach (string line in status.Split('\n'))
                    if (line.StartsWith(name + ":", StringComparison.Ordinal)) return line[(name.Length + 1)..].Replace("kB", "", StringComparison.Ordinal).Trim();
                return "-";
            }
            return string.Join('\t', "rss_kb=" + Field("Rss"), "anonhuge_kb=" + Field("AnonHugePages"), "filehuge_kb=" + Field("FilePmdMapped"), "shmemhuge_kb=" + Field("ShmemPmdMapped"), "rssfile_kb=" + St("RssFile"), "rssanon_kb=" + St("RssAnon"));
        }
        catch (IOException)
        {
            return "smaps=unavailable";
        }
    }

    /// <summary>SplitMix64 with the modulo rule of the test support class (the patterns of 1b.7f are drawn the same way).</summary>
    private sealed class Mix(ulong seed)
    {
        private ulong _state = seed;

        private ulong Next()
        {
            unchecked
            {
                _state += 0x9E3779B97F4A7C15UL;
                ulong z = _state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
                return z ^ (z >> 31);
            }
        }

        internal int NextInt(int bound) => (int)(Next() % (ulong)bound);
    }
}
