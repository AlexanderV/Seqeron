using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Sais;

/// <summary>
/// The stopwatch counters of <see cref="SaisTrace"/> and <see cref="LcpTrace"/> (work package 1a.9 profile): an instrumented build gives the same
/// result as an untraced one, every reached level has timed phases, and the phases account for nearly all of the build time.
/// </summary>
[TestFixture]
public class PhaseTraceTests
{
    [Test]
    public void TracedBuild_EqualsUntracedBuild_OnAllFamilies()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(1500, 3))
        {
            PreparedText p = TextPreparation.Prepare(text);
            int[] plain, traced;
            CompactLcp plainLcp, tracedLcp;
            var trace = new SaisTrace();
            var lcpTrace = new LcpTrace();
            if (p.Width == SymbolWidth.Byte)
            {
                plain = SaisReference.Build(p.ByteReader(), p.AlphabetSize);
                traced = SaisReference.Build(p.ByteReader(), p.AlphabetSize, trace);
                plainLcp = LcpReference.Build(p.ByteReader(), plain);
                tracedLcp = LcpReference.Build(p.ByteReader(), traced, lcpTrace);
            }
            else
            {
                plain = SaisReference.Build(p.UShortReader(), p.AlphabetSize);
                traced = SaisReference.Build(p.UShortReader(), p.AlphabetSize, trace);
                plainLcp = LcpReference.Build(p.UShortReader(), plain);
                tracedLcp = LcpReference.Build(p.UShortReader(), traced, lcpTrace);
            }
            Assert.That(traced, Is.EqualTo(plain), name);
            Assert.That(tracedLcp.ToArray(), Is.EqualTo(plainLcp.ToArray()), name);
            Assert.That(tracedLcp.Kind, Is.EqualTo(plainLcp.Kind), name);
            Assert.That(trace.Depth, Is.GreaterThanOrEqualTo(1), name);
            Assert.That(lcpTrace.PhiTicks, Is.GreaterThanOrEqualTo(0), name);
        }
    }

    [Test]
    public void Phases_AreTimedPerLevel_AndAccountForTheBuild()
    {
        string text = StringFamilies.RandomDna(400_000, 5);
        PreparedText p = TextPreparation.Prepare(text);
        var trace = new SaisTrace();
        var lcpTrace = new LcpTrace();

        long start = System.Diagnostics.Stopwatch.GetTimestamp();
        int[] sa = SaisReference.Build(p.ByteReader(), p.AlphabetSize, trace);
        long saisTicks = System.Diagnostics.Stopwatch.GetTimestamp() - start;
        start = System.Diagnostics.Stopwatch.GetTimestamp();
        LcpReference.Build(p.ByteReader(), sa, lcpTrace);
        long lcpTicks = System.Diagnostics.Stopwatch.GetTimestamp() - start;

        Assert.That(trace.Depth, Is.GreaterThanOrEqualTo(2), "random DNA needs at least one recursion level");
        long sum = 0;
        foreach (SaisPhase phase in Enum.GetValues<SaisPhase>())
        {
            long total = trace.TotalTicks(phase);
            Assert.That(total, Is.GreaterThanOrEqualTo(trace.Ticks(0, phase)), phase.ToString());
            sum += total;
        }
        foreach (SaisPhase phase in new[] { SaisPhase.Count, SaisPhase.Classify, SaisPhase.Round1InduceL, SaisPhase.Round1InduceS, SaisPhase.Round2InduceL, SaisPhase.Round2InduceS })
            Assert.That(trace.Ticks(0, phase), Is.GreaterThan(0), "level 0 " + phase);
        Assert.That(trace.Ticks(trace.Depth + 5, SaisPhase.Name), Is.Zero, "a level that was not reached has no time");
        Assert.That(sum, Is.LessThanOrEqualTo(saisTicks), "exclusive phases cannot exceed the whole build");
        Assert.That(sum, Is.GreaterThan(saisTicks / 2), "the phases account for the bulk of the build (the rest is array allocation and the entry points)");

        long lcpSum = lcpTrace.PhiTicks + lcpTrace.PlcpTicks + lcpTrace.CountTicks + lcpTrace.OutputTicks;
        Assert.That(lcpSum, Is.LessThanOrEqualTo(lcpTicks));
        Assert.That(lcpSum, Is.GreaterThan(lcpTicks / 2));
        Assert.That(lcpTrace.PlcpTicks, Is.GreaterThan(0));
    }
}
