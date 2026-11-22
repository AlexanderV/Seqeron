using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.Esa.Tests.Handoff;

/// <summary>
/// Work package 1b.9f, orchestrator's audit: <see cref="HandoffCompare"/> decides which answers of the suffix-array search and of the engine differ from the oracle, so every signal of it is tested alone:
/// a wrong count above the cap (the only signal there), a wrong "exists" flag with a right count, a wrong locate set at exactly the cap, a missing locate set. Written after the audit's mutants of
/// <c>HandoffCompare</c> survived the first 20 tests (a wrong count, a wrong flag and a set at exactly the cap were not seen by any test).
/// </summary>
[TestFixture]
public sealed class HandoffCompareTests
{
    private static HandoffQuery Q(int index = 0) => new(index, true, 0, "ACGT"u8.ToArray());

    private static ScanResult Oracle(long count, params long[] positions) => new([count], [[.. positions]]);

    private static List<HandoffMismatch> Compare(ScanResult oracle, QueryAnswer answer, int cap) => HandoffCompare.Compare("test", [Q()], oracle, [answer], cap);

    [Test]
    public void AnAnswerThatEqualsTheOracle_IsNotAMismatch()
    {
        Assert.That(Compare(Oracle(3, 1, 5, 9), new QueryAnswer(3, true, [1, 5, 9]), 10), Is.Empty);
        Assert.That(Compare(Oracle(0), new QueryAnswer(0, false, []), 10), Is.Empty);
        Assert.That(Compare(Oracle(50, 1, 2), new QueryAnswer(50, true, null), 2), Is.Empty, "above the cap only the count and the flag are compared");
    }

    [Test]
    public void AWrongCountAboveTheCap_IsReported_AndItIsTheOnlySignalThere()
    {
        List<HandoffMismatch> m = Compare(Oracle(5, 1, 2, 3), new QueryAnswer(4, true, null), 3);
        Assert.That(m.Select(x => x.Kind), Is.EqualTo(new[] { "count" }));
    }

    [Test]
    public void AWrongExistsFlag_WithTheRightCount_IsReported()
    {
        Assert.That(Compare(Oracle(0), new QueryAnswer(0, true, []), 10).Select(x => x.Kind), Is.EqualTo(new[] { "exists" }));
        Assert.That(Compare(Oracle(2, 4, 8), new QueryAnswer(2, false, [4, 8]), 10).Select(x => x.Kind), Is.EqualTo(new[] { "exists" }));
    }

    [Test]
    public void AWrongLocateSetAtExactlyTheCap_IsReported()
    {
        Assert.That(Compare(Oracle(3, 1, 2, 3), new QueryAnswer(3, true, [1, 2, 9]), 3).Select(x => x.Kind), Is.EqualTo(new[] { "locate" }), "a count equal to the cap is still compared");
        Assert.That(Compare(Oracle(3, 1, 2, 3), new QueryAnswer(3, true, [1, 2, 9]), 2), Is.Empty, "one above the cap it is not");
    }

    [Test]
    public void AnAnswerWithoutPositions_ForACountWithinTheCap_IsReported_WhenTheCountIsRight_AndOnlyOnce()
    {
        List<HandoffMismatch> m = Compare(Oracle(2, 4, 8), new QueryAnswer(2, true, null), 10);
        Assert.That(m.Select(x => x.Kind), Is.EqualTo(new[] { "locate" }));
        Assert.That(m[0].Detail, Does.Contain("no positions"));
        // with a wrong count too, the count is the report and the positions are not asked for twice
        Assert.That(Compare(Oracle(2, 4, 8), new QueryAnswer(3, true, null), 10).Select(x => x.Kind), Is.EqualTo(new[] { "count" }));
    }

    [Test]
    public void ALocateSetThatDiffersInOrLength_IsReportedWithTheFirstDifference()
    {
        List<HandoffMismatch> m = Compare(Oracle(3, 1, 5, 9), new QueryAnswer(3, true, [1, 5]), 10);
        Assert.That(m.Single(x => x.Kind == "locate").Detail, Does.Contain("index 2"));
        m = Compare(Oracle(3, 1, 5, 9), new QueryAnswer(3, true, [1, 6, 9]), 10);
        Assert.That(m.Single().Detail, Does.Contain("index 1 (6 against 5)"));
    }

    [Test]
    public void EveryQueryIsJudgedByItsOwnAnswer()
    {
        var qs = new List<HandoffQuery> { Q(0), Q(1), Q(2) };
        var oracle = new ScanResult([1, 0, 2], [[7], [], [3, 4]]);
        var answers = new[] { new QueryAnswer(1, true, [7]), new QueryAnswer(1, true, [7]), new QueryAnswer(2, true, [3, 4]) };
        List<HandoffMismatch> m = HandoffCompare.Compare("who", qs, oracle, answers, 10);
        Assert.That(m.Select(x => x.Index).Distinct(), Is.EqualTo(new[] { 1 }));
        Assert.That(m.Select(x => x.Answerer).Distinct(), Is.EqualTo(new[] { "who" }));
    }
}
