using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.Esa.External;

// The decomposition of the per-call cost of the index at one thread (work package 1b.9g, item 2): the in-memory index of the text and the index over the container, in one process, in blocks of 0.25 s in alternating order
// (median of 7 blocks per figure), the descent alone and the materialisation of the list separately, so that a difference between the two indexes has a place.
// usage: <pattern file> <container> <text file> [C|T]      (C = Checked access, the default; T = Trusted)
// rows (ns per call): count of present patterns, of absent patterns, locate (FindAllOccurrences), count of the locate patterns, TryMatch of the locate patterns (the descent), Sa.ToArray of the intervals, Sa.CopyTo into a reused buffer, new int[count]
internal static class P
{
    static double Time(Func<long> body, double minSeconds = 0.4)
    {
        // returns ns per call
        long calls = 0; long sink = 0; long t0 = Stopwatch.GetTimestamp();
        while (Stopwatch.GetElapsedTime(t0).TotalSeconds < minSeconds) { sink += body(); calls += n; }
        double s = Stopwatch.GetElapsedTime(t0).TotalSeconds; GC.KeepAlive(sink);
        return s / calls * 1e9;
    }
    static int n;
    static void Main(string[] a)
    {
        string[] lines = File.ReadAllLines(a[0]);
        string[] pc = lines.Where(l => l[0] == 'P').Select(l => l[2..]).ToArray();
        string[] pl = lines.Where(l => l[0] == 'L').Select(l => l[2..]).ToArray();
        string[] pa = lines.Where(l => l[0] == 'A').Select(l => l[2..]).ToArray();
        string text = File.ReadAllText(a[2], System.Text.Encoding.Latin1);
        EsaIndex mem = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        EsaIndex fil = EsaIndex.FromFile(SesaFile.Open(a[1], new EsaOpenOptions { Access = a.Length > 3 && a[3] == "T" ? EsaAccessMode.Trusted : EsaAccessMode.Checked }));
        var names = new[] { "mem", "file" };
        var ix = new[] { mem, fil };
        // intervals of the locate patterns
        var loHi = new (int lo, int hi)[2][];
        for (int k = 0; k < 2; k++) { loHi[k] = pl.Select(p => { ix[k].TryMatch(p.AsSpan(), out int lo, out int hi); return (lo, hi); }).ToArray(); }
        double hits = loHi[0].Average(x => x.hi - x.lo + 1);
        Console.WriteLine($"patterns P {pc.Length} L {pl.Length} A {pa.Length}; mean locate hits {hits:F1}");
        var rows = new List<(string, double[])>();
        void Run(string label, Func<int, Func<long>> make, int count)
        {
            var r = new double[2][]; for (int k = 0; k < 2; k++) r[k] = new double[7];
            n = count;
            // warm
            for (int k = 0; k < 2; k++) { var b = make(k); Time(b, 0.3); }
            for (int rep = 0; rep < 7; rep++) for (int k = 0; k < 2; k++) { var b = make(k); r[k][rep] = Time(b, 0.25); }
            double m0 = r[0].OrderBy(x => x).ElementAt(3), m1 = r[1].OrderBy(x => x).ElementAt(3);
            Console.WriteLine($"{label,-34} mem {m0,8:F1} ns  file {m1,8:F1} ns  diff {m1 - m0,7:F1} ns  ratio file/mem {m1 / m0:F3}");
        }
        // 1: full count of present patterns (descent + facade-free)
        Run("count present (CountOccurrences)", k => () => { long s = 0; foreach (string p in pc) s += ix[k].CountOccurrences(p); return s; }, pc.Length);
        Run("count absent", k => () => { long s = 0; foreach (string p in pa) s += ix[k].CountOccurrences(p); return s; }, pa.Length);
        Run("locate (FindAllOccurrences)", k => () => { long s = 0; foreach (string p in pl) s += ix[k].FindAllOccurrences(p).Count; return s; }, pl.Length);
        Run("count of the locate patterns", k => () => { long s = 0; foreach (string p in pl) s += ix[k].CountOccurrences(p); return s; }, pl.Length);
        Run("TryMatch of the locate patterns", k => () => { long s = 0; foreach (string p in pl) { ix[k].TryMatch(p.AsSpan(), out int lo, out int hi); s += hi - lo; } return s; }, pl.Length);
        Run("Sa.ToArray of the intervals", k => () => { long s = 0; var sa = ix[k].Sa; foreach (var (lo, hi) in loHi[k]) s += sa.ToArray(lo, hi - lo + 1).Length; return s; }, pl.Length);
        Run("Sa.CopyTo into a reused buffer", k => () => { long s = 0; var sa = ix[k].Sa; var buf = new int[4096]; foreach (var (lo, hi) in loHi[k]) { int len = Math.Min(hi - lo + 1, 4096); sa.CopyTo(lo, buf.AsSpan(0, len)); s += buf[0]; } return s; }, pl.Length);
        Run("new int[count] only (no copy)", k => () => { long s = 0; foreach (var (lo, hi) in loHi[k]) s += new int[hi - lo + 1].Length; return s; }, pl.Length);
    }
}
