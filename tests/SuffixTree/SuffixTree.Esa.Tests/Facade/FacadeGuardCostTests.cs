using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The cost of the reader count of <see cref="EsaSuffixTree"/> (work package 1b.8, the Dispose race) measured inside one process, where the machine noise of separate processes (a spread of 40 % between processes of the
/// same binary on this machine) does not enter: the same tree, the same patterns, batches of the facade call and of the call of the index behind it in alternating order, medians of the batches. The difference of the two
/// is the whole cost of the facade (the reader count, the argument check, the span of the string); the member <c>IsEmpty</c> against the read of the index length isolates the count itself. Category <c>Cost</c>,
/// <c>[Explicit]</c>, driven by <c>analysis/em/run_1b8_dispose_cost.sh</c> (variables <c>SESA_FILE</c> and <c>SESA_INPUT</c> as in <see cref="SesaReaderCostTests"/>, <c>SESA_ROWS</c> for the rows, <c>SESA_LABEL</c> for the label of a row).
/// </summary>
[TestFixture]
[Category("Cost")]
[Explicit("cost of the reader count of 1b.8: run by analysis/em/run_1b8_dispose_cost.sh, which sets SESA_FILE, SESA_INPUT and SESA_ROWS")]
public class FacadeGuardCostTests
{
    private const int Batches = 15;

    private static double Median(List<double> values)
    {
        double[] s = [.. values.Order()];
        return s[s.Length / 2];
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    [Test]
    public void TheReaderCount_CostsFewNanosecondsPerCall()
    {
        string file = Environment.GetEnvironmentVariable("SESA_FILE") ?? throw new InvalidOperationException("SESA_FILE");
        string input = Environment.GetEnvironmentVariable("SESA_INPUT") ?? throw new InvalidOperationException("SESA_INPUT");
        string rows = Environment.GetEnvironmentVariable("SESA_ROWS") ?? throw new InvalidOperationException("SESA_ROWS");
        string label = Environment.GetEnvironmentVariable("SESA_LABEL") ?? "run";
        string text = File.ReadAllText(input, System.Text.Encoding.Latin1);
        var rng = new SplitMix64(1234);
        string[] twenty = [.. Enumerable.Range(0, 100_000).Select(_ => text.Substring(rng.NextInt(text.Length - 20), 20))];
        string[] eight = [.. Enumerable.Range(0, 1_000).Select(_ => text.Substring(rng.NextInt(text.Length - 8), 8))];
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        using EsaSuffixTree checkedTree = SesaTestFile.Open(file, EsaValidationLevel.Header, EsaAccessMode.Checked);
        using EsaSuffixTree trusted = SesaTestFile.Open(file, EsaValidationLevel.Header, EsaAccessMode.Trusted);
        foreach ((string name, EsaSuffixTree tree) in new[] { ("memory", memory), ("checked", checkedTree), ("trusted", trusted) })
        {
            EsaIndex index = tree.Index;
            long sink = 0;
            for (int w = 0; w < 2; w++)
            {
                foreach (string p in twenty) sink += tree.CountOccurrences(p) + index.CountOccurrences(p.AsSpan());
                foreach (string p in eight) sink += tree.FindAllOccurrences(p).Count + index.FindAllOccurrences(p.AsSpan()).Count;
            }
            var countFacade = new List<double>();
            var countIndex = new List<double>();
            var locateFacade = new List<double>();
            var locateIndex = new List<double>();
            for (int b = 0; b < Batches; b++)
            {
                for (int order = 0; order < 2; order++)
                {
                    bool facadeFirst = ((b + order) & 1) == 0;
                    for (int part = 0; part < 2; part++)
                    {
                        bool facade = (part == 0) == facadeFirst;
                        long t = Stopwatch.GetTimestamp();
                        if (facade)
                            foreach (string p in twenty) sink += tree.CountOccurrences(p);
                        else
                            foreach (string p in twenty) sink += index.CountOccurrences(p.AsSpan());
                        (facade ? countFacade : countIndex).Add(Seconds(t) * 1e9 / twenty.Length);
                        t = Stopwatch.GetTimestamp();
                        if (facade)
                            foreach (string p in eight) sink += tree.FindAllOccurrences(p).Count;
                        else
                            foreach (string p in eight) sink += index.FindAllOccurrences(p.AsSpan()).Count;
                        (facade ? locateFacade : locateIndex).Add(Seconds(t) * 1e9 / eight.Length);
                    }
                }
            }

            // the count itself: IsEmpty (the reader count and one comparison) against the read of the length of the index (the comparison alone)
            const int loops = 20_000_000;
            var guardMember = new List<double>();
            var guardBare = new List<double>();
            for (int b = 0; b < Batches; b++)
            {
                long t = Stopwatch.GetTimestamp();
                int hits = 0;
                for (int i = 0; i < loops; i++)
                    if (tree.IsEmpty) hits++;
                guardMember.Add(Seconds(t) * 1e9 / loops);
                t = Stopwatch.GetTimestamp();
                for (int i = 0; i < loops; i++)
                    if (index.Length == 0) hits++;
                guardBare.Add(Seconds(t) * 1e9 / loops);
                sink += hits;
            }

            // four threads on the one tree (the reader count is one cache line that every reader writes): the facade call against the call of the index, per call and thread
            var contendedFacade = new List<double>();
            var contendedIndex = new List<double>();
            for (int b = 0; b < 7; b++)
            {
                for (int part = 0; part < 2; part++)
                {
                    bool facade = ((b + part) & 1) == 0;
                    long t = Stopwatch.GetTimestamp();
                    var threads = Enumerable.Range(0, 4).Select(_ => new Thread(() =>
                    {
                        long local = 0;
                        if (facade)
                            foreach (string p in twenty) local += tree.CountOccurrences(p);
                        else
                            foreach (string p in twenty) local += index.CountOccurrences(p.AsSpan());
                        Interlocked.Add(ref sink, local);
                    })).ToArray();
                    foreach (Thread thread in threads) thread.Start();
                    foreach (Thread thread in threads) thread.Join();
                    (facade ? contendedFacade : contendedIndex).Add(Seconds(t) * 1e9 / twenty.Length);
                }
            }
            GC.KeepAlive(sink);
            double kf = Median(contendedFacade), ki = Median(contendedIndex);
            double cf = Median(countFacade), ci = Median(countIndex), lf = Median(locateFacade), li = Median(locateIndex), gm = Median(guardMember), gb = Median(guardBare);
            File.AppendAllText(rows, string.Create(CultureInfo.InvariantCulture,
                $"guard\t{label}\t{name}\t{cf:F1}\t{ci:F1}\t{(cf - ci):F1}\t{(100 * (cf - ci) / cf):F2}\t{lf:F1}\t{li:F1}\t{(lf - li):F1}\t{(100 * (lf - li) / lf):F2}\t{gm:F2}\t{gb:F2}\t{(gm - gb):F2}\t{kf:F1}\t{ki:F1}\t{(kf - ki):F1}\t{(100 * (kf - ki) / kf):F2}\n"));
        }
        Assert.Pass();
    }
}
