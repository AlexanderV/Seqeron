// allocprobe: how exact is GC.GetAllocatedBytesForCurrentThread when another thread makes garbage collections (work package 1b.9h, item 3; see README.md)
// usage: allocprobe <helper> <seconds> [region 0|1] [window n|s|l]
//   helper  0 none, 1 garbage of 2 KB flat out (gen 0 collections), 2 GC.Collect(0) + Sleep(1), 3 GC.Collect(0) tight, 4 large garbage (200 KB) flat out (gen 2 / background collections), 5 GC.Collect(2, Forced, blocking: false) + Sleep(1)
//   region  1: every window is opened with GC.TryStartNoGCRegion(16 MB, 8 MB) and closed with EndNoGCRegion (what the test project's AllocationMeter does, with a 64 MB budget)
//   window  n: the code of the window allocates nothing, s: one long[11] (112 bytes), l: one long[200001] on the large object heap (1,600,032 bytes); every window also computes for about a millisecond
// prints the windows, how many had a collection inside, and how many read more than 992 bytes above the expected figure (the slack of the 1b.3 test is 1,024 bytes)
using System.Diagnostics;
using System.Reflection;
using System.Runtime;
if (args.Length > 0 && args[0] == "real") { RealTest.Run(args); return; }
int mode = int.Parse(args[0]); int seconds = int.Parse(args[1]); bool region = args.Length > 2 && args[2] == "1"; string wtype = args.Length > 3 ? args[3] : "l"; long expected = wtype == "l" ? 1_600_032 : wtype == "s" ? 112 : 0;
var stop = false;
if (mode != 0)
{
    new Thread(() =>
    {
        var keep = new object[64]; int i = 0;
        while (!Volatile.Read(ref stop))
        {
            switch (mode)
            {
                case 1: keep[i++ & 63] = new byte[2000]; break;
                case 2: GC.Collect(0, GCCollectionMode.Forced, true); Thread.Sleep(1); break;
                case 3: GC.Collect(0, GCCollectionMode.Forced, true); break;
                case 4: keep[i++ & 7] = new byte[200_000]; break;
                case 5: GC.Collect(2, GCCollectionMode.Forced, false); Thread.Sleep(1); break;
            }
        }
    }) { IsBackground = true }.Start();
}
var large = new SortedDictionary<long, long>();
long anomaliesIntact = 0, anomaliesLost = 0; long iters = 0, regionLost = 0, regionRefused = 0, anomalies = 0, gcInWindow = 0;
var sw = Stopwatch.StartNew();
object? sink = null; long acc = 0;
while (sw.Elapsed.TotalSeconds < seconds)
{
    iters++;
    for (int k = 0; k < (iters % 7); k++) sink = new byte[100 + (int)(iters % 300)];
    bool inRegion = false;
    if (region) { try { inRegion = GC.TryStartNoGCRegion(16_000_000, 8_000_000); } catch (InvalidOperationException) { inRegion = false; } if (!inRegion) regionRefused++; }
    int g0 = GC.CollectionCount(0), g2 = GC.CollectionCount(2);
    long a = GC.GetAllocatedBytesForCurrentThread();
    if (wtype == "l") sink = new long[200001]; else if (wtype == "s") sink = new long[11];
    long s = 0; for (int k = 0; k < 1_500_000; k++) s += k ^ (s >> 3);
    acc += s;
    long b = GC.GetAllocatedBytesForCurrentThread();
    bool gc = GC.CollectionCount(0) != g0 || GC.CollectionCount(2) != g2;
    if (gc) gcInWindow++;
    bool intact = !gc;
    if (inRegion) { if (gc || GCSettings.LatencyMode != GCLatencyMode.NoGCRegion) regionLost++; try { GC.EndNoGCRegion(); } catch (InvalidOperationException) { } }
    long d = b - a;
    if (region && inRegion) { if (intact && d > expected + 992) anomaliesIntact++; if (!intact && d > expected + 992) anomaliesLost++; }
    large[d] = large.GetValueOrDefault(d) + 1;
    if (d > expected + 992) anomalies++;
}
Volatile.Write(ref stop, true);
GC.KeepAlive(sink);
Console.WriteLine($"helper {mode} region {region} window {wtype}: windows {iters}, with a GC inside {gcInWindow}, region refused {regionRefused}, region lost {regionLost}, windows over the limit of the 1b.3 test (+992 B): {anomalies} (region intact {anomaliesIntact}, region lost {anomaliesLost}); distinct deltas {large.Count}, max {large.Keys.Max()} (+{large.Keys.Max() - expected}), min {large.Keys.Min()} ({acc & 1})");


// the test of work package 1b.3 (CountDistinctSubstrings allocates nothing, CountDistinctSubstringsByLength(10) 112 bytes, ByLength(n) the result array), its body as it was before 1b.9h and as it is now, over the real
// SuffixTree.Esa.dll (internal EsaIndex, called through delegates so that the windows hold no reflection), with or without a thread that makes large garbage; usage:
//   allocprobe real <SuffixTree.Esa.dll> <helper 0|4|5> <seconds> [old|new|both]    (default both: the two bodies alternate in one process; old or new alone: only that body runs)    helper 4: 200 KB flat out (background gen 2 collections, the worst case), 5: 100 KB per millisecond (what a quiet neighbour does)
static class RealTest
{
    public static void Run(string[] a)
    {
        string dll = a[1]; int helper = int.Parse(a[2]); int seconds = int.Parse(a[3]); string which = a.Length > 4 ? a[4] : "both"; bool runOld = which != "new", runNew = which != "old";
        Assembly asm = Assembly.LoadFrom(dll);
        Type t = asm.GetType("SuffixTree.Esa.EsaIndex")!;
        var rng = new Random(11);
        const string dna = "ACGT";
        string text = string.Create(200_000, rng, (span, r) => { for (int i = 0; i < span.Length; i++) span[i] = dna[r.Next(4)]; });
        object esa = t.GetMethod("Build", BindingFlags.Static | BindingFlags.NonPublic | BindingFlags.Public)!.Invoke(null, new object?[] { text, null })!;
        var total = (Func<long>)Delegate.CreateDelegate(typeof(Func<long>), esa, t.GetMethod("CountDistinctSubstrings", BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public)!);
        var byLength = (Func<int, long[]>)Delegate.CreateDelegate(typeof(Func<int, long[]>), esa, t.GetMethod("CountDistinctSubstringsByLength", BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public)!);
        total(); byLength(10);                                   // the JIT, as in the test
        var stop = false;
        if (helper != 0)
            new Thread(() =>
            {
                var keep = new object[8]; int i = 0;
                while (!Volatile.Read(ref stop)) { keep[i++ & 7] = helper == 4 ? new byte[200_000] : new byte[100_000]; if (helper == 5) Thread.Sleep(1); }
            }) { IsBackground = true }.Start();
        long limitLarge = 8L * (text.Length + 1) + 1024;
        long runs = 0, oldFail = 0, oldTotal = 0, oldSmall = 0, oldLarge = 0, newFail = 0, newInconclusive = 0, newLost = 0, maxOld = 0;
        long[]? full = null;
        var sw = Stopwatch.StartNew();
        while (sw.Elapsed.TotalSeconds < seconds)
        {
            runs++;
            // the body of the test as it was: three windows of the plain counter of the thread
            if (runOld)
            {
            long before = GC.GetAllocatedBytesForCurrentThread(); total(); long t1 = GC.GetAllocatedBytesForCurrentThread() - before;
            before = GC.GetAllocatedBytesForCurrentThread(); byLength(10); long s1 = GC.GetAllocatedBytesForCurrentThread() - before;
            before = GC.GetAllocatedBytesForCurrentThread(); full = byLength(text.Length); long l1 = GC.GetAllocatedBytesForCurrentThread() - before;
            bool bad = false;
            if (t1 >= 1024) { oldTotal++; bad = true; }
            if (s1 >= 1024) { oldSmall++; bad = true; }
            if (l1 >= limitLarge) { oldLarge++; bad = true; }
            if (bad) oldFail++;
            maxOld = Math.Max(maxOld, Math.Max(t1, Math.Max(s1, l1 - 8L * (text.Length + 1))));
            }
            if (!runNew) continue;
            // the body as it is now: the same three windows, each in a no-GC region, a reading above the limit counts only when the window was free of collections
            (long b, bool free) t2 = Window(() => total()), s2 = Window(() => byLength(10)), l2 = Window(() => full = byLength(text.Length));
            bool fail = (t2.b >= 1024 && t2.free) || (s2.b >= 1024 && s2.free) || (l2.b >= limitLarge && l2.free);
            bool inconclusive = !fail && ((t2.b >= 1024 && !t2.free) || (s2.b >= 1024 && !s2.free) || (l2.b >= limitLarge && !l2.free));
            if (fail) newFail++;
            if (inconclusive) newInconclusive++;
            if (!t2.free || !s2.free || !l2.free) newLost++;
        }
        Volatile.Write(ref stop, true);
        GC.KeepAlive(full);
        Console.WriteLine($"real test ({which}), helper {helper}: {runs} runs of the three windows in {seconds} s; as it was (plain counter): {oldFail} runs failed ({oldFail * 100.0 / runs:F1} %: total {oldTotal}, ByLength(10) {oldSmall}, ByLength(n) {oldLarge}; the largest excess over the expected figure +{maxOld} bytes); " +
                          $"as it is now (no-GC region): {newFail} failed, {newInconclusive} inconclusive, {newLost} runs with a window that was not free of collections; gen0 {GC.CollectionCount(0)} gen1 {GC.CollectionCount(1)} gen2 {GC.CollectionCount(2)}");
    }

    static (long Bytes, bool Free) Window(Action code)
    {
        bool inRegion;
        try { inRegion = GC.TryStartNoGCRegion(64L << 20, 32L << 20); } catch (InvalidOperationException) { inRegion = false; }
        try
        {
            int c = GC.CollectionCount(0);
            long before = GC.GetAllocatedBytesForCurrentThread();
            code();
            long after = GC.GetAllocatedBytesForCurrentThread();
            return (after - before, inRegion && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == c);
        }
        finally
        {
            if (inRegion) { try { if (GCSettings.LatencyMode == GCLatencyMode.NoGCRegion) GC.EndNoGCRegion(); } catch (InvalidOperationException) { } }
        }
    }
}
