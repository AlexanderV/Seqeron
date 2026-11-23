using System.Diagnostics;
using System.Reflection;
using System.Runtime.Loader;

// In-process A/B client of the difference cover of the suffix-tree upgrade (S1b, work package 1b.9h, item 4): two builds of SuffixTree.Esa.dll side by side, each in its own load context, and the call that the production
// chain makes once per build (SaBuilder.CoverFor(nominal period, exact: false), which is DifferenceCover.Create(NearestSingerPeriod(nominal))) made through a delegate. Every measurement is an A B B A block; the figure is the
// median over the blocks of the ratio of the means of the B calls and the A calls, with its quartiles. The same build on both sides is the null experiment (it shows the spread of the instrument).
// usage: abcover split <dir> <nominal period> <calls>            where the time of one Create goes: the whole call, the search of the Singer set, the rest (the tables IndexOf and Witness), and the cost of the check Verify + IsPerfect
//        abcover ab <dirA> <dirB> <nominal period> <calls> <blocks>   the A/B of the two builds
//        abcover same <dirA> <dirB> <nominal period>                  the two covers compared residue by residue (exit code 1 when they differ)
//        abcover residues <dir> <nominal period>                      the period, the kind and the residues of the cover of the build (the table of the default period was made from this line)
// <dir> holds a build of SuffixTree.Esa.dll and its dependencies (dotnet build src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj -c Release -o <dir>).

sealed class Ctx(string dir) : AssemblyLoadContext(false)
{
    protected override Assembly? Load(AssemblyName n)
    {
        string p = Path.Combine(dir, n.Name + ".dll");
        return File.Exists(p) ? LoadFromAssemblyPath(p) : null;
    }
}

sealed class Build(string dir)
{
    const BindingFlags All = BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static;
    public readonly Assembly Asm = new Ctx(dir).LoadFromAssemblyPath(Path.GetFullPath(Path.Combine(dir, "SuffixTree.Esa.dll")));
    public Type Cover => Asm.GetType("SuffixTree.Esa.External.DifferenceCover", true)!;
    public Type Builder => Asm.GetType("SuffixTree.Esa.External.SaBuilder", true)!;

    public Func<int, object> CoverFor()
    {
        MethodInfo m = Builder.GetMethod("CoverFor", All)!;
        return nominal => m.Invoke(null, [nominal, false])!;
    }

    public int[] Residues(object cover) => (int[])cover.GetType().GetProperty("Residues", All)!.GetValue(cover)!;

    public Func<int, int[]>? SingerSetOnly()
    {
        MethodInfo? m = Cover.GetMethod("SingerSet", All);
        return m is null ? null : q => (int[])m.Invoke(null, [q])!;
    }

    public bool Verify(object cover) => (bool)cover.GetType().GetMethod("Verify", All)!.Invoke(cover, null)!;

    public bool IsPerfect(object cover) => (bool)cover.GetType().GetMethod("IsPerfect", All)!.Invoke(cover, null)!;

    public int NearestSinger(int nominal) => (int)Cover.GetMethod("NearestSingerPeriod", All)!.Invoke(null, [nominal])!;
}

static class P
{
    static double Time(Action a, int calls)
    {
        long t0 = Stopwatch.GetTimestamp();
        for (int i = 0; i < calls; i++) a();
        return Stopwatch.GetElapsedTime(t0).TotalMilliseconds / calls;
    }

    static double Median(List<double> x) { x.Sort(); return x[x.Count / 2]; }

    static int Main(string[] a)
    {
        switch (a[0])
        {
            case "split":
            {
                var b = new Build(a[1]);
                int nominal = int.Parse(a[2]), calls = int.Parse(a[3]);
                var cf = b.CoverFor();
                int period = b.NearestSinger(nominal);
                object cover = cf(nominal);
                Func<int, int[]>? singer = b.SingerSetOnly();
                int q = (int)Math.Round((-1 + Math.Sqrt(1 + (4.0 * (period - 1)))) / 2);
                Console.WriteLine($"period {period} (q {q}), residues {b.Residues(cover).Length}, first call (cold, includes the JIT) {Time(() => cf(nominal), 1):F2} ms");
                double whole = Time(() => cf(nominal), calls);
                Console.WriteLine($"CoverFor (nearest Singer period + Create), mean of {calls} calls: {whole:F3} ms");
                if (singer is not null)
                {
                    double search = Time(() => singer(q), calls);
                    Console.WriteLine($"SingerSet({q}) alone (the search for a primitive cubic): {search:F3} ms = {100 * search / whole:F1} % of the call");
                    Console.WriteLine($"the rest of Create (IndexOf and Witness tables, the nearest period): {whole - search:F3} ms");
                }
                double check = Time(() => { b.Verify(cover); b.IsPerfect(cover); }, Math.Max(1000, calls * 100));
                Console.WriteLine($"Verify + IsPerfect of the cover (what a table that is verified when used would cost): {check * 1000:F1} microseconds");
                return 0;
            }
            case "ab":
            {
                var A = new Build(a[1]); var B = new Build(a[2]);
                int nominal = int.Parse(a[3]), calls = int.Parse(a[4]), blocks = int.Parse(a[5]);
                var fa = A.CoverFor(); var fb = B.CoverFor();
                for (int w = 0; w < 3; w++) { fa(nominal); fb(nominal); }
                var ra = new List<double>(); var rb = new List<double>(); var pr = new List<double>();
                for (int i = 0; i < blocks; i++)
                {
                    double x1, y1, y2, x2;
                    if (i % 2 == 0) { x1 = Time(() => fa(nominal), calls); y1 = Time(() => fb(nominal), calls); y2 = Time(() => fb(nominal), calls); x2 = Time(() => fa(nominal), calls); }
                    else { y1 = Time(() => fb(nominal), calls); x1 = Time(() => fa(nominal), calls); x2 = Time(() => fa(nominal), calls); y2 = Time(() => fb(nominal), calls); }
                    double xa = (x1 + x2) / 2, yb = (y1 + y2) / 2; ra.Add(xa); rb.Add(yb); pr.Add(yb / xa);
                }
                pr.Sort();
                Console.WriteLine($"AB cover  A {Median(ra):F4} ms  B {Median(rb):F4} ms  B/A median-of-paired {Median(new List<double>(pr)):F4}  p25 {pr[pr.Count / 4]:F4} p75 {pr[3 * pr.Count / 4]:F4}  min {pr[0]:F4} max {pr[^1]:F4}  n {pr.Count}  (calls per pass {calls}; A - B = {Median(ra) - Median(rb):F4} ms)");
                return 0;
            }
            case "residues":
            {
                var b = new Build(a[1]);
                int nominal = int.Parse(a[2]);
                object cover = b.CoverFor()(nominal);
                Console.WriteLine($"period {b.NearestSinger(nominal)}, kind {cover.GetType().GetProperty("Kind", BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(cover)}, residues {b.Residues(cover).Length}, verify {b.Verify(cover)}, perfect {b.IsPerfect(cover)}");
                Console.WriteLine(string.Join(", ", b.Residues(cover)));
                return 0;
            }
            case "same":
            {
                var A = new Build(a[1]); var B = new Build(a[2]);
                int nominal = int.Parse(a[3]);
                int[] ra = A.Residues(A.CoverFor()(nominal)), rb = B.Residues(B.CoverFor()(nominal));
                bool same = ra.SequenceEqual(rb);
                Console.WriteLine($"nominal {nominal}: {ra.Length} and {rb.Length} residues, {(same ? "identical" : "DIFFERENT")}");
                return same ? 0 : 1;
            }
        }
        return 2;
    }
}
