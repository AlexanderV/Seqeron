using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Scaffolding;

/// <summary>The bookkeeping of <see cref="HangLimitAttribute"/> (work package 1b.9h): a test that begins while another one is still inside is the sign that the other one hung, and from then on no guarded test starts.</summary>
[TestFixture]
public class HangBreakerTests
{
    [Test]
    public void TestsThatComeBack_OneAfterTheOther_AreAllStarted()
    {
        var breaker = new HangBreaker();
        foreach (string name in new[] { "A", "B", "C", "A" })
        {
            Assert.That(breaker.Begin(name), Is.Null, name);
            breaker.End(name);
        }
        Assert.That(breaker.Runaway, Is.Null);
    }

    [Test]
    public void ATestThatDidNotComeBack_IsNamed_ByTheNextOne_AndByEveryOneAfterIt()
    {
        var breaker = new HangBreaker();
        Assert.That(breaker.Begin("Hung"), Is.Null);
        Assert.That(breaker.Begin("Next"), Is.EqualTo("Hung"), "the next test is refused and told which one did not come back");
        Assert.That(breaker.Runaway, Is.EqualTo("Hung"));
        breaker.End("Next");
        Assert.That(breaker.Begin("Third"), Is.EqualTo("Hung"), "and so is every test after it");
        Assert.That(breaker.Runaway, Is.EqualTo("Hung"), "the first one stays the one that is named");
    }

    [Test]
    public void TheThreadThatIsStuck_ComingBackLater_DoesNotOpenTheBreaker()
    {
        var breaker = new HangBreaker();
        breaker.Begin("Hung");
        Assert.That(breaker.Begin("Next"), Is.EqualTo("Hung"));
        breaker.End("Hung"); // the runaway thread ended after all (its test was already reported as failed)
        Assert.That(breaker.Begin("Later"), Is.EqualTo("Hung"), "the host had a runaway thread: its timings and its memory are no longer those of a quiet run");
    }

    [Test]
    public void ARefusedTest_DoesNotMarkItselfAsInside_AndTheEndOfAnotherTestLeavesTheMarkAlone()
    {
        var breaker = new HangBreaker();
        breaker.Begin("A");
        breaker.End("Other"); // not the test that is inside: nothing changes
        Assert.That(breaker.Begin("B"), Is.EqualTo("A"), "A is still inside");
        var fresh = new HangBreaker();
        Assert.That(fresh.Begin("A"), Is.Null);
        Assert.That(fresh.Begin("A"), Is.Null, "the same test begun again (a repeat) is not a test that did not come back");
        Assert.That(fresh.Runaway, Is.Null);
    }

    [Test]
    public void TheAttribute_SetsTheLimitInMilliseconds_OnTheContextOfTheTest()
    {
        var attribute = new HangLimitAttribute(90);
        Assert.That(attribute.Seconds, Is.EqualTo(90));
        Assert.That(attribute.Targets, Is.EqualTo(ActionTargets.Test));
        var context = new NUnit.Framework.Internal.TestExecutionContext.AdhocContext();
        attribute.ApplyToContext(context);
        Assert.That(context.TestCaseTimeout, Is.EqualTo(90_000), "the limit that the tests below a class take");
        Assert.That(attribute.Properties.Get(NUnit.Framework.Internal.PropertyNames.Timeout), Is.EqualTo(90_000), "the limit that a method takes (NUnit reads the property of the test, and it wins over the limit of the class)");
        Assert.Throws<OverflowException>(() => new HangLimitAttribute(int.MaxValue).ApplyToContext(context), "a limit that does not fit milliseconds is refused, not wrapped around");
    }
}
