using System.Diagnostics;
using System.Reflection;
using System.Runtime.Loader;

// In-process A/B client of the suffix-tree upgrade (S1b, work package 1b.9g; the instrument of the verifications of 1b.9b): two builds of SuffixTree.Esa.dll (or two kinds of tree of one build) side by side in ONE process, each in its
// own load context, called through delegates; every measurement is an A B B A block (the order of the blocks alternates), and the result is the median over the blocks of the per-block ratio B / A. Across processes a ratio of 5-10 % cannot be
// resolved (the throughput of a process falls into one of two groups about 1.4x apart); in process it resolves 1-2 % (the null experiment, the same build twice, shows the spread).
// usage: abapp <patternfile> <passes> <threads> <seconds> <specA> <specB>
//   spec = <dir>|mem:<textfile>     the in-memory tree of the build in <dir> (a directory that holds SuffixTree.Esa.dll and its dependencies), built from the text file
//        | <dir>|file:<container>   the facade over a SESA container (the file-backed tree, Checked; "fileT:" = Trusted)
//        | <dir>|idx:<container>    the index behind that facade (no reader registration; "idxT:" = Trusted)
//   pattern file: lines "P <pattern>" (count of a present pattern), "L <pattern>" (locate: FindAllOccurrences), "A <pattern>" (count of an absent pattern); genpat.py makes them
//   passes: blocks per operation; threads: 1 = one thread, n = n threads on one tree, each from its own offset of the pattern list; seconds: length of one pass (0 with one thread = one sweep of the list)
// output: one line per operation: the median throughput of A and of B (calls per second, summed over the threads), the median of the per-block ratios B / A, its quartiles and the number of blocks.

sealed class Ctx(string dir) : AssemblyLoadContext(false)
{
    protected override Assembly? Load(AssemblyName n)
    {
        string p = Path.Combine(dir, n.Name + ".dll");
        return File.Exists(p) ? LoadFromAssemblyPath(p) : null;
    }
}

sealed record V(string Spec, Func<string, int> Count, Func<string, IReadOnlyList<int>> Locate);

static class P
{
    const BindingFlags All = BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static;

    static V Make(string spec)
    {
        string[] parts = spec.Split('|');
        string dir = parts[0], kind = parts[1];
        var ctx = new Ctx(dir);
        Assembly asm = ctx.LoadFromAssemblyPath(Path.Combine(dir, "SuffixTree.Esa.dll"));
        object tree;
        if (kind.StartsWith("mem:"))
        {
            Type f = asm.GetType("SuffixTree.Esa.EsaSuffixTreeFactory")!;
            MethodInfo m = f.GetMethod("Create", All, null, new[] { typeof(string) }, null) ?? f.GetMethod("CreateEsa", All, null, new[] { typeof(string) }, null)!;
            tree = m.Invoke(null, new object[] { File.ReadAllText(kind[4..]) })!;
        }
        else
        {
            Type f = asm.GetType("SuffixTree.Esa.SuffixIndexFactory")!;
            MethodInfo m = f.GetMethod("Open", All, null, new[] { typeof(string), asm.GetType("SuffixTree.Esa.EsaOpenOptions")! }, null)!;
            Type ot = asm.GetType("SuffixTree.Esa.EsaOpenOptions")!;
            object opts = Activator.CreateInstance(ot)!;
            string head = kind[..kind.IndexOf(':')];
            if (head.EndsWith("T")) ot.GetProperty("Access")!.SetValue(opts, Enum.ToObject(asm.GetType("SuffixTree.Esa.EsaAccessMode")!, 1));
            tree = m.Invoke(null, new object?[] { kind[(kind.IndexOf(':') + 1)..], opts })!;
            if (head.StartsWith("idx")) tree = tree.GetType().GetProperty("Index", All)!.GetValue(tree)!;
        }
        Type t = tree.GetType();
        MethodInfo c = t.GetMethod("CountOccurrences", All, null, new[] { typeof(string) }, null)!;
        MethodInfo l = t.GetMethod("FindAllOccurrences", All, null, new[] { typeof(string) }, null)!;
        return new V(spec, (Func<string, int>)Delegate.CreateDelegate(typeof(Func<string, int>), tree, c), (Func<string, IReadOnlyList<int>>)Delegate.CreateDelegate(typeof(Func<string, IReadOnlyList<int>>), tree, l));
    }

    static long sink;

    static double Pass(V v, int op, string[][] w, int th, double secs)
    {
        string[] pats = w[op];
        if (th <= 1)
        {
            long t0 = Stopwatch.GetTimestamp(); long ops = 0, sk = 0;
            if (secs <= 0) { Run(v, op, pats, 0, ref sk); ops = pats.Length; }
            else { var end = t0 + (long)(secs * Stopwatch.Frequency); while (Stopwatch.GetTimestamp() < end) { Run(v, op, pats, (int)(ops % 97), ref sk); ops += pats.Length; } }
            double r = ops / Stopwatch.GetElapsedTime(t0).TotalSeconds; Interlocked.Add(ref sink, sk); return r;
        }
        long total = 0; var go = new ManualResetEventSlim(); var sw = new Stopwatch();
        var ts = Enumerable.Range(0, th).Select(k => new Thread(() =>
        {
            go.Wait(); long ops = 0, sk = 0; int i = k * 97;
            var end = Stopwatch.GetTimestamp() + (long)(secs * Stopwatch.Frequency);
            while (Stopwatch.GetTimestamp() < end) { for (int j = 0; j < 64; j++) { sk += op == 1 ? v.Locate(pats[i++ % pats.Length]).Count : v.Count(pats[i++ % pats.Length]); ops++; } }
            Interlocked.Add(ref total, ops); Interlocked.Add(ref sink, sk);
        })).ToList();
        ts.ForEach(x => x.Start()); sw.Start(); go.Set(); ts.ForEach(x => x.Join());
        return total / sw.Elapsed.TotalSeconds;
    }

    static void Run(V v, int op, string[] pats, int off, ref long sk)
    {
        if (op == 1) foreach (string p in pats) sk += v.Locate(p).Count; else foreach (string p in pats) sk += v.Count(p);
    }

    static double Median(List<double> x) { x.Sort(); return x[x.Count / 2]; }

    static void Main(string[] a)
    {
        string[] lines = File.ReadAllLines(a[0]);
        string[][] w = { lines.Where(l => l[0] == 'P').Select(l => l[2..]).ToArray(), lines.Where(l => l[0] == 'L').Select(l => l[2..]).ToArray(), lines.Where(l => l[0] == 'A').Select(l => l[2..]).ToArray() };
        // op 0 = count present, op 1 = locate, op 2 = count absent
        int passes = int.Parse(a[1]), th = int.Parse(a[2]); double secs = double.Parse(a[3]);
        V A = Make(a[4]), B = Make(a[5]);
        string[] names = { "count_present", "locate", "count_absent" };
        for (int wu = 0; wu < 4; wu++) foreach (V v in new[] { A, B }) for (int op = 0; op < 3; op++) Pass(v, op, w, 1, 0);
        for (int op = 0; op < 3; op++)
        {
            var ra = new List<double>(); var rb = new List<double>(); var pr = new List<double>();
            for (int i = 0; i < passes; i++)
            {
                double x1, y1, y2, x2;
                if (i % 2 == 0) { x1 = Pass(A, op, w, th, secs); y1 = Pass(B, op, w, th, secs); y2 = Pass(B, op, w, th, secs); x2 = Pass(A, op, w, th, secs); }
                else { y1 = Pass(B, op, w, th, secs); x1 = Pass(A, op, w, th, secs); x2 = Pass(A, op, w, th, secs); y2 = Pass(B, op, w, th, secs); }
                double xa = (x1 + x2) / 2, yb = (y1 + y2) / 2; ra.Add(xa); rb.Add(yb); pr.Add(yb / xa);
            }
            pr.Sort();
            Console.WriteLine($"AB {names[op],-14} A {Median(ra),12:F0}  B {Median(rb),12:F0}  B/A median-of-paired {Median(new List<double>(pr)):F4}  p25 {pr[pr.Count / 4]:F4} p75 {pr[3 * pr.Count / 4]:F4}  n {pr.Count}");
        }
        Console.WriteLine("sink " + sink);
    }
}
