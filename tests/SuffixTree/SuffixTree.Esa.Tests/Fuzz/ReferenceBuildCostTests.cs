using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The build time and the memory of the reference trees at the sizes above the harness bands (work package 1b.8, G5: "extend the v6 comparison to the largest n at which v6 builds in reasonable time"): the v6 tree in heap mode,
/// the Ukkonen tree and the in-memory ESA on one random DNA text per size, one process, the managed heap growth after each build (a forced collection first). Category <c>Cost</c>, <c>[Explicit]</c>; the rows go to the file
/// named by <c>SESA_ROWS</c>; sizes from <c>SESA_SIZES</c> (default 1000000,4000000,10000000).
/// </summary>
[TestFixture]
[Category("Cost")]
[Explicit("build time and memory of the reference trees at 10^6 .. 10^7 symbols (G5 of 1b.8); run with --filter FullyQualifiedName~ReferenceBuildCostTests")]
public class ReferenceBuildCostTests
{
    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    private static long Heap()
    {
        GC.Collect(2, GCCollectionMode.Forced, true);
        return GC.GetTotalMemory(true);
    }

    [Test]
    public void BuildTimeAndHeap_OfV6_Ukkonen_AndTheEsa_AtLargeSizes()
    {
        string rows = Environment.GetEnvironmentVariable("SESA_ROWS") ?? throw new InvalidOperationException("SESA_ROWS");
        int[] sizes = [.. (Environment.GetEnvironmentVariable("SESA_SIZES") ?? "1000000,4000000,10000000").Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture))];
        File.AppendAllText(rows, "#refbuild: n engine build_s heap_mb_after_build note\n");
        foreach (int n in sizes)
        {
            string text = StringFamilies.RandomDna(n, 7);
            long before = Heap();
            string note = "ok";
            long t = Stopwatch.GetTimestamp();
            try
            {
                using PersistentSuffixTree v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
                double s = Seconds(t);
                File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture, $"refbuild\t{n}\tv6-heap\t{s:F2}\t{(Heap() - before) / 1048576.0:F0}\t{note}\n"));
            }
            catch (Exception e)
            {
                File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture, $"refbuild\t{n}\tv6-heap\t{Seconds(t):F2}\t0\t{e.GetType().Name}: {e.Message.Replace('\n', ' ')}\n"));
            }
            before = Heap();
            t = Stopwatch.GetTimestamp();
            try
            {
                UkkonenTree u = UkkonenTree.Build(text);
                double s = Seconds(t);
                File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture, $"refbuild\t{n}\tukkonen\t{s:F2}\t{(Heap() - before) / 1048576.0:F0}\t{note}\n"));
                GC.KeepAlive(u);
            }
            catch (Exception e)
            {
                File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture, $"refbuild\t{n}\tukkonen\t{Seconds(t):F2}\t0\t{e.GetType().Name}: {e.Message.Replace('\n', ' ')}\n"));
            }
            before = Heap();
            t = Stopwatch.GetTimestamp();
            using (EsaSuffixTree esa = EsaSuffixTree.Build(text))
            {
                double s = Seconds(t);
                File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture, $"refbuild\t{n}\tesa-memory\t{s:F2}\t{(Heap() - before) / 1048576.0:F0}\t{note}\n"));
            }
        }
        Assert.Pass();
    }
}
