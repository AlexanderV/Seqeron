using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Heavy checks of work package 1b.7c (gate G1). Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>TMPDIR=/dev/shm dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~LcpChildExhaustiveTests</c> (the runs msync and fsync tens of small files each: a memory
/// file system makes them cheap). Every string over sigma 1 up to length 12, sigma 2 up to 12, sigma 3 up to 6, sigma 4 up to 5 (the 9,282 strings of E6 and the unary and 4-ary lists) goes through the whole stage
/// with rotating tiny parameters (R, blocks, windows, slots, threads, both codecs); the families go through it at 10^5 and 10^6 symbols with several chunk sizes and thread counts. Every output is compared with the in-RAM structures.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter FullyQualifiedName~LcpChildExhaustiveTests (see PROGRESS-S1.md, Log 1b.7c)")]
public class LcpChildExhaustiveTests
{
    [Test]
    public void EveryString_OverSigma1To4_UpToTheLengthsOfE6_ThroughTheWholeStage()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (1, 12), (2, 12), (3, 6), (4, 5) })
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                LcpChildCase c = LcpChildCase.Of(s);
                ChildCodec codec = count % 3 == 2 ? ChildCodec.Tagged : ChildCodec.Packed;
                using var dir = new TempDir();
                LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, LcpChildRun.Options(dir.Sub("scratch"), LcpChildRun.TinyTuning(count), 1 + (count % 4), codec, c.Sigma));
                LcpChildRun.Compare(c, o, codec, $"'{s}' variant {count} {codec}");
                count++;
            }
        }
        Assert.That(count, Is.EqualTo(12 + 8190 + 1092 + 1364), "unary 1-12, binary 1-12, ternary 1-6, 4-ary 1-5");
    }

    [Test]
    public void TheFamilies_AtTenToTheFifth_AndTenToTheSixth_WithSeveralChunkSizesAndThreadCounts()
    {
        int v = 0;
        foreach (int n in new[] { 100_000, 1_000_000 })
        {
            foreach ((string name, string text) in LcpChildRun.Families(n))
            {
                LcpChildCase c = LcpChildCase.Of(text);
                foreach ((long r, int threads, int block) in new[] { (n / 7L, 4, 4096), (8_000_000L, 1, 256 * 1024), (n / 3L + 1, 2, 100_000) })
                {
                    using var dir = new TempDir();
                    var t = new LcpChildTuning { LcpChunk = r, BlockEntries = block };
                    LcpChildRun.Outcome o = LcpChildRun.Run(dir, c, new LcpChildOptions { ScratchDir = dir.Sub("scratch"), Threads = threads, Tuning = t, Sigma = c.Sigma });
                    LcpChildRun.Compare(c, o, ChildCodec.Packed, $"{name} {n} R={r} threads {threads} block {block}");
                    v++;
                }
            }
        }
        Assert.That(v, Is.EqualTo(2 * 7 * 3));
    }
}
