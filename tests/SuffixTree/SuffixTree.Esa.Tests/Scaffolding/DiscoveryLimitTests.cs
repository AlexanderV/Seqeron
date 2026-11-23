using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Scaffolding;

/// <summary>The limit on what a test source runs at test discovery (work package 1b.9h): a value that comes back in time is returned, a failure is thrown again with its own type, and a computation that does not come back gives false at the limit.</summary>
[TestFixture]
public class DiscoveryLimitTests
{
    [Test]
    public void AComputationThatReturns_GivesItsValue()
    {
        Assert.That(DiscoveryLimit.TryRun(TimeSpan.FromSeconds(30), () => 6 * 7, out int value), Is.True);
        Assert.That(value, Is.EqualTo(42));
        Assert.That(DiscoveryLimit.TryRun(TimeSpan.FromSeconds(30), () => "x", out string? text), Is.True);
        Assert.That(text, Is.EqualTo("x"));
    }

    [Test]
    public void AComputationThatThrows_ThrowsAgainWithItsOwnTypeAndMessage()
    {
        var e = Assert.Throws<InvalidDataException>(() => DiscoveryLimit.TryRun<int>(TimeSpan.FromSeconds(30), () => throw new InvalidDataException("damaged"), out _));
        Assert.That(e!.Message, Is.EqualTo("damaged"));
    }

    [Test]
    public void AComputationThatDoesNotReturn_GivesFalse_AtTheLimit_AndLeavesNoValue()
    {
        using var release = new ManualResetEventSlim();
        long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
        bool returned = DiscoveryLimit.TryRun(TimeSpan.FromMilliseconds(300), () => { release.Wait(); return 5; }, out int value);
        double seconds = System.Diagnostics.Stopwatch.GetElapsedTime(t0).TotalSeconds;
        release.Set();
        Assert.That(returned, Is.False);
        Assert.That(value, Is.Zero);
        Assert.That(seconds, Is.GreaterThanOrEqualTo(0.25).And.LessThan(20), "it waited for the limit and not for the computation");
    }

    [Test]
    public void TheComputation_RunsOnAThreadOfItsOwn_ThatDoesNotKeepTheHostAlive()
    {
        bool background = false;
        int thread = 0;
        Assert.That(DiscoveryLimit.TryRun(TimeSpan.FromSeconds(30), () =>
        {
            background = Thread.CurrentThread.IsBackground;
            thread = Environment.CurrentManagedThreadId;
            return 0;
        }, out int _), Is.True);
        Assert.That(background, Is.True);
        Assert.That(thread, Is.Not.EqualTo(Environment.CurrentManagedThreadId));
    }
}
