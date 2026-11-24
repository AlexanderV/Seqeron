using System.Globalization;
using System.Text;

// WidthBench: see ../README.md. Commands:
//   synthetic --layouts w4,w5,w8,hyb,p34 --n 10000000,100000000 --v 6200000000 --threads 1,4 --reps 5 --ops 4000000 --out file.tsv --tag text
//   null      --layouts ... --reps 7 --ops 20000000 --out file.tsv --tag text
//   gen       --n 100000000 --out dna.txt [--seed 1]           random ACGT text (SplitMix64), one symbol per byte
//   descent   (Descent.cs)   the count descent over a real container with the width swapped
//   text      (TextBench.cs) the read cost of the packed text against the byte text
internal static class Program
{
    internal static CultureInfo Inv => CultureInfo.InvariantCulture;

    internal static Dictionary<string, string> Opts(string[] args, int from)
    {
        var d = new Dictionary<string, string>();
        for (int i = from; i < args.Length; i += 2)
        {
            if (!args[i].StartsWith("--") || i + 1 >= args.Length) throw new ArgumentException("bad argument " + args[i]);
            d[args[i][2..]] = args[i + 1];
        }
        return d;
    }

    internal static string Get(Dictionary<string, string> o, string k, string dflt = null) => o.TryGetValue(k, out string v) ? v : dflt ?? throw new ArgumentException("--" + k + " is required");

    internal static int Main(string[] args)
    {
        if (args.Length == 0) { Console.Error.WriteLine("usage: see Program.cs"); return 2; }
        var o = Opts(args, 1);
        switch (args[0])
        {
            case "synthetic": return RunSynthetic(o);
            case "null": return RunNull(o);
            case "gen": return Gen(o);
            case "selftest": return SelfTest();
            case "descent": return Descent.Run(o);
            case "text": return TextBench.Run(o);
            case "stats": return Stats.Run(o);
            default: Console.Error.WriteLine("unknown command"); return 2;
        }
    }

    private static void Emit(string path, IEnumerable<Synthetic.Row> rows)
    {
        bool header = !File.Exists(path) || new FileInfo(path).Length == 0;
        using var w = new StreamWriter(path, true, Encoding.ASCII);
        if (header) w.WriteLine("kind\tlayout\tn\tthreads\tmetric\tmedian\tmin\tmax\tunit\treps\tnotes");
        foreach (Synthetic.Row r in rows)
            w.WriteLine(string.Join('\t', r.Kind, r.Layout, r.N.ToString(Inv), r.Threads.ToString(Inv), r.Metric, r.Median.ToString("F3", Inv), r.Min.ToString("F3", Inv), r.Max.ToString("F3", Inv), r.Unit, r.Reps.ToString(Inv), r.Notes));
    }

    private static int RunSynthetic(Dictionary<string, string> o)
    {
        string[] layouts = Get(o, "layouts").Split(',');
        long[] ns = Get(o, "n").Split(',').Select(x => long.Parse(x, Inv)).ToArray();
        ulong v = ulong.Parse(Get(o, "v", "6200000000"), Inv);
        int[] threads = Get(o, "threads", "1,4").Split(',').Select(int.Parse).ToArray();
        int reps = int.Parse(Get(o, "reps", "5"), Inv);
        long ops = long.Parse(Get(o, "ops", "4000000"), Inv);
        string[] kinds = Get(o, "kinds", "chase,thr,scan").Split(',');
        int[] scans = Get(o, "scans", "8,1500,100000").Split(',').Select(int.Parse).ToArray();
        string path = Get(o, "out");
        string tag = Get(o, "tag", "");
        foreach (long n in ns)
            foreach (string layout in layouts)
            {
                var sw = System.Diagnostics.Stopwatch.StartNew();
                var rows = Synthetic.RunLayout(layout, n, v, kinds, threads, reps, ops, scans, tag);
                Emit(path, rows);
                Console.WriteLine($"{layout} n={n} done in {sw.Elapsed.TotalSeconds:F1} s");
            }
        return 0;
    }

    private static int RunNull(Dictionary<string, string> o)
    {
        string[] layouts = Get(o, "layouts").Split(',');
        int reps = int.Parse(Get(o, "reps", "7"), Inv);
        long ops = long.Parse(Get(o, "ops", "20000000"), Inv);
        string path = Get(o, "out");
        string tag = Get(o, "tag", "");
        Emit(path, Synthetic.Empty(reps, ops, tag));
        foreach (string layout in layouts) Emit(path, Synthetic.Null(layout, reps, ops, tag));
        return 0;
    }

    private static unsafe int SelfTest()
    {
        int bad = 0;
        foreach (long n in new long[] { 1, 14, 15, 16, 63, 64, 65, 959, 960, 961, 1_000_003 })
            foreach (string layout in new[] { "w4", "w5", "w8", "hyb", "p34" })
            {
                ulong v = 6_200_000_000UL;
                using Packed a = Packed.Create(layout, n, v, 99UL, 3);
                for (long i = 0; i < n; i++)
                {
                    ulong want = Rng.Value(99UL, i, v);
                    ulong got = layout switch
                    {
                        "w4" => new R4(a.Main).Get(i),
                        "w5" => new R5(a.Main).Get(i),
                        "w8" => new R8(a.Main).Get(i),
                        "p34" => new R34(a.Main).Get(i),
                        _ => new RHyb(a.Main, a.Blocks, a.Longs).Get(i),
                    };
                    if (layout == "w4") want &= 0xFFFFFFFFUL;
                    if (got != want) { if (bad++ < 10) Console.WriteLine($"MISMATCH {layout} n={n} i={i} got={got} want={want}"); }
                }
                if (layout == "hyb" && n > 1000)
                {
                    var h = new RHyb(a.Main, a.Blocks, a.Longs);
                    ulong a1 = h.SumRange(10, 500), b1 = 0;
                    for (long i = 10; i < 510; i++) b1 += h.Get(i);
                    if (a1 != b1) { bad++; Console.WriteLine("MISMATCH hyb SumRange"); }
                    Console.WriteLine($"hyb n={n} escapes={a.Escapes} ({100.0 * a.Escapes / n:F1} %) bytes/entry={(double)a.Bytes / n:F3}");
                }
            }
        Console.WriteLine(bad == 0 ? "selftest OK" : $"selftest FAILED {bad}");
        return bad == 0 ? 0 : 1;
    }

    private static int Gen(Dictionary<string, string> o)
    {
        long n = long.Parse(Get(o, "n"), Inv);
        ulong seed = ulong.Parse(Get(o, "seed", "1"), Inv);
        string path = Get(o, "out");
        byte[] map = "ACGT"u8.ToArray();
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        byte[] buf = new byte[1 << 20];
        ulong x = seed;
        long done = 0;
        while (done < n)
        {
            int m = (int)Math.Min(buf.Length, n - done);
            for (int i = 0; i < m; i += 32)
            {
                ulong r = Rng.Mix(x++);
                for (int j = 0; j < 32 && i + j < m; j++) { buf[i + j] = map[(int)(r & 3)]; r >>= 2; }
            }
            fs.Write(buf, 0, m);
            done += m;
        }
        return 0;
    }
}
