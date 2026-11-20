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
        if (args.Length > 0 && args[0] == "gen") return Generate.Run(args);
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
        if (variant == "open") return OpenBench(file);
        if (variant == "opencold") return OpenCold(file, textPath);
        string text = File.ReadAllText(textPath, System.Text.Encoding.Latin1);

#if EXPERIMENTS
        if (variant == "count") return ExperimentsCount.Run(file, textPath);
#endif
        if (variant.StartsWith("mt:", StringComparison.Ordinal)) return ThreadBench(variant[3..].Split('+'), file, text, passes, label, alphabet, seed);
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

    /// <summary>
    /// A build of SuffixTree.Esa loaded into a context of its own (the dependencies of its directory, except SuffixTree.Core, which stays the one of the harness so that the ISuffixTree of both builds is one type): several
    /// commits of the in-memory tree in ONE process, with the passes in a rotated order. Variant <c>alc=&lt;path of SuffixTree.Esa.dll&gt;</c> builds the in-memory tree of the text with that build's EsaSuffixTreeFactory.Create.
    /// </summary>
    private sealed class BuildContext(string path) : System.Runtime.Loader.AssemblyLoadContext(isCollectible: false)
    {
        private readonly string _directory = Path.GetDirectoryName(path)!;

        protected override System.Reflection.Assembly? Load(System.Reflection.AssemblyName name)
        {
            if (name.Name == "SuffixTree.Core") return null;
            string candidate = Path.Combine(_directory, name.Name + ".dll");
            return File.Exists(candidate) ? LoadFromAssemblyPath(candidate) : null;
        }
    }

    private static Q InMemoryOfBuild(string path, string text)
    {
        var context = new BuildContext(path);
        System.Reflection.Assembly assembly = context.LoadFromAssemblyPath(path);
        Type factory = assembly.GetType("SuffixTree.Esa.EsaSuffixTreeFactory")!;
        object tree = factory.GetMethod("Create", [typeof(string)])!.Invoke(null, [text])!;
        return new TreeQ((ISuffixTree)tree);
    }

    /// <summary>Variant <c>alcf=&lt;path of SuffixTree.Esa.dll&gt;</c>: the file-backed tree (Checked) of the container, opened by that build.</summary>
    private static Q FileOfBuild(string path, string file)
    {
        var context = new BuildContext(path);
        System.Reflection.Assembly assembly = context.LoadFromAssemblyPath(path);
        Type factory = assembly.GetType("SuffixTree.Esa.EsaSuffixTreeFactory")!;
        Type options = assembly.GetType("SuffixTree.Esa.EsaOpenOptions")!;
        object o = Activator.CreateInstance(options)!;
        options.GetProperty("Access")!.SetValue(o, Enum.Parse(assembly.GetType("SuffixTree.Esa.EsaAccessMode")!, "Checked"));
        object tree = factory.GetMethod("Open", [typeof(string), options])!.Invoke(null, [file, o])!;
        return new TreeQ((ISuffixTree)tree);
    }

    private static Q CreateQ(string variant, string file, string text)
    {
        if (variant.StartsWith("alc=", StringComparison.Ordinal)) return InMemoryOfBuild(variant[4..], text);
        if (variant.StartsWith("alcf=", StringComparison.Ordinal)) return FileOfBuild(variant[5..], file);
#if EXPERIMENTS
        if (variant.StartsWith("ix-", StringComparison.Ordinal)) return ExperimentsIx.Create(variant, file, text);
#endif
        return new TreeQ(Create(variant, file, text));
    }

    /// <summary>
    /// The three pattern sets of the benchmark definition (BENCHMARK-1b9b.md): present long (count), present short (locate), absent (count). The lengths are 20 / 8 / 24 unless RB_L20 / RB_L8 / RB_L24 say otherwise, the numbers
    /// 200000 / 20000 / 100000 unless RB_N20 / RB_N8 / RB_N24 do; the absent patterns are drawn uniformly from the alphabet (the argument, or the distinct symbols of the text when it is "@text").
    /// </summary>
    private static (string[] Twenty, string[] Eight, string[] Absent) Patterns(string text, string alphabet, ulong seed)
    {
        var rng = new Mix(seed);
        int Env(string name, int fallback) => int.Parse(Environment.GetEnvironmentVariable(name) ?? fallback.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
        int m20 = Env("RB_N20", 200000), m8 = Env("RB_N8", 20000), m24 = Env("RB_N24", 100000);
        int l20 = Env("RB_L20", 20), l8 = Env("RB_L8", 8), l24 = Env("RB_L24", 24);
        if (alphabet == "@text") alphabet = new string(text.Distinct().Order().ToArray());
        string[] twenty = new string[m20];
        for (int i = 0; i < twenty.Length; i++) twenty[i] = text.Substring(rng.NextInt(text.Length - l20), l20);
        string[] eight = new string[m8];
        for (int i = 0; i < eight.Length; i++) eight[i] = text.Substring(rng.NextInt(text.Length - l8), l8);
        string[] absent = new string[m24];
        for (int i = 0; i < absent.Length; i++)
        {
            var chars = new char[l24];
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

    /// <summary>
    /// The four-thread instrument (work package 1b.9b, item 4): several trees (or several access paths to one tree) in ONE process, the same persistent threads for all of them. In every round, for every variant in a rotated order, the
    /// threads (RB_THREADS, default 4) are released together by a barrier and each runs its slice of the pattern set (count of the present long patterns, locate of the short ones, count of the absent ones) on the tree of the
    /// variant, the phase repeated RB_REPEAT times (default 3); the throughput is the queries of the phase over the wall time between the release and the last thread to finish. Output: per variant the median over the rounds of the
    /// three throughputs and the median over the rounds of the per-round ratio to the first variant (the pairing that cancels the drift of the machine). `passes` is the number of rounds.
    /// </summary>
    private static int ThreadBench(string[] variants, string file, string text, int rounds, string label, string alphabet, ulong seed)
    {
        int threads = int.Parse(Environment.GetEnvironmentVariable("RB_THREADS") ?? "4", CultureInfo.InvariantCulture);
        int repeat = int.Parse(Environment.GetEnvironmentVariable("RB_REPEAT") ?? "3", CultureInfo.InvariantCulture);
        var qs = new Q[variants.Length];
        for (int i = 0; i < qs.Length; i++) qs[i] = CreateQ(variants[i], file, text);
        (string[] twenty, string[] eight, string[] absent) = Patterns(text, alphabet, seed);
        string[][] sets = [twenty, eight, absent];
        int k = qs.Length;

        // the persistent threads: phase = (variant, metric) written by the main thread before the start barrier
        var start = new Barrier(threads + 1);
        var end = new Barrier(threads + 1);
        int phaseVariant = 0, phaseMetric = 0;
        bool stop = false;
        long sink = 0;
        var workers = new Thread[threads];
        for (int t = 0; t < threads; t++)
        {
            int id = t;
            workers[t] = new Thread(() =>
            {
                long local = 0;
                while (true)
                {
                    start.SignalAndWait();
                    if (Volatile.Read(ref stop))
                    {
                        Interlocked.Add(ref sink, local);
                        return;
                    }
                    Q q = qs[Volatile.Read(ref phaseVariant)];
                    int metric = Volatile.Read(ref phaseMetric);
                    string[] set = sets[metric];
                    int from = set.Length / threads * id, to = id == threads - 1 ? set.Length : set.Length / threads * (id + 1);
                    for (int rep = 0; rep < repeat; rep++)
                        for (int i = from; i < to; i++) local += metric == 1 ? q.Locate(set[i]) : q.Count(set[i]);
                    end.SignalAndWait();
                }
            }) { IsBackground = true };
            workers[t].Start();
        }

        double Phase(int variant, int metric)
        {
            phaseVariant = variant;
            phaseMetric = metric;
            start.SignalAndWait();
            long t0 = Stopwatch.GetTimestamp();
            end.SignalAndWait();
            return repeat * sets[metric].Length / Seconds(t0);
        }

        // warm-up by time: every variant, every metric, about one second each
        for (int v = 0; v < k; v++)
            for (int metric = 0; metric < 3; metric++)
            {
                long w0 = Stopwatch.GetTimestamp();
                do { _ = Phase(v, metric); } while (Seconds(w0) < 1.0);
            }
        var res = new List<double>[k, 3];
        for (int v = 0; v < k; v++)
            for (int metric = 0; metric < 3; metric++) res[v, metric] = [];
        for (int r = 0; r < rounds; r++)
            for (int j = 0; j < k; j++)
            {
                int v = (r + j) % k;
                for (int metric = 0; metric < 3; metric++) res[v, metric].Add(Phase(v, metric));
            }
        Volatile.Write(ref stop, true);
        start.SignalAndWait();
        foreach (Thread w in workers) w.Join();
        GC.KeepAlive(sink);
        for (int v = 0; v < k; v++)
        {
            string Ratio(int metric) => Median(res[v, metric].Select((x, p) => x / res[0, metric][p]).ToList()).ToString("F3", CultureInfo.InvariantCulture);
            Console.WriteLine(string.Join('\t', "mt", label, threads.ToString(CultureInfo.InvariantCulture), variants[v],
                Median(res[v, 0]).ToString("F0", CultureInfo.InvariantCulture), Median(res[v, 1]).ToString("F0", CultureInfo.InvariantCulture), Median(res[v, 2]).ToString("F0", CultureInfo.InvariantCulture),
                "ratio_to_first:", Ratio(0), Ratio(1), Ratio(2)));
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

    /// <summary>
    /// The first open of a big container in a JIT-warm process: five opens of a small container warm the code, then one open of the big one is timed (Header level; the driver drops the page cache before the process starts, so
    /// the header page, the packed text's layout and the alphabet are read from the disk). usage: ReadBench opencold &lt;big container&gt; &lt;small container&gt;
    /// </summary>
    private static int OpenCold(string big, string small)
    {
#if NOOPEN
        throw new NotSupportedException("this build has no Open");
#else
        for (int i = 0; i < 5; i++) EsaSuffixTreeFactory.Open(small).Dispose();
        long t0 = Stopwatch.GetTimestamp();
        EsaSuffixTree tree = EsaSuffixTreeFactory.Open(big);
        double seconds = Seconds(t0);
        tree.Dispose();
        Console.WriteLine(string.Join('\t', "opencold", "header", seconds.ToString("F6", CultureInfo.InvariantCulture)));
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
