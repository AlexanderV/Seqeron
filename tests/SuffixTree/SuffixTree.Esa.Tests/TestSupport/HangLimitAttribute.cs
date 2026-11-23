using NUnit.Framework.Interfaces;
using NUnit.Framework.Internal;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// The second half of <see cref="HangLimitAttribute"/>, as an object of its own so that it can be tested without touching the state of the host: the test whose thread is inside its test, and the first test that did not
/// come back. A test that begins while another one is still inside is the sign that the other one hung (NUnit reported it as failed at its limit and went on; its thread is still running).
/// </summary>
internal sealed class HangBreaker
{
    private string? _inside;
    private string? _runaway;

    /// <summary>The first test that did not come back, or null.</summary>
    internal string? Runaway => Volatile.Read(ref _runaway);

    /// <summary>A test is about to start. Returns null when it may, otherwise the name of the earlier test that did not come back.</summary>
    internal string? Begin(string name)
    {
        string? before = Volatile.Read(ref _inside);
        if (before is not null && before != name) Interlocked.CompareExchange(ref _runaway, before, null);
        string? runaway = Volatile.Read(ref _runaway);
        if (runaway is null) Volatile.Write(ref _inside, name);
        return runaway;
    }

    /// <summary>A test came back (passed, failed or was ignored).</summary>
    internal void End(string name) => Interlocked.CompareExchange(ref _inside, null, name);
}

/// <summary>
/// A time limit per test that turns a hang of the code under test into a failed test (work package 1b.9h, debt D: 14 of the 28 mutants of 1b.9c stopped the test host, because a mutated loop or work pool never ends and the default
/// set has no limit per test, so the only net was the time limit of the whole job). On a class or a method. It does two things:
/// <list type="number">
/// <item>it sets the limit that NUnit's own <see cref="TimeoutAttribute"/> sets, in the two places where NUnit reads it: the property <c>Timeout</c> of the test (a limit on a method; it wins over the limit of the class, larger or smaller) and the context of the tests below a class (the test runs on a thread of its own; when the limit is over the test is reported as failed, "Test exceeded Timeout value", and the run goes on; on .NET
/// the thread that is stuck cannot be stopped, it is left running, and the host can still end);</item>
/// <item>it does not start another test of a guarded class while such a thread is still inside its test: the next test fails at once, naming the one that did not come back. A hung loop of the product code hangs every test that
/// runs it, and without this each of them would wait for its own limit (a class of 25 tests at 60 s each is 25 minutes for one bug); with it the first test pays the limit and the others fail in milliseconds, and no further
/// call is made that could hang the same way and leave one more thread behind.</item>
/// </list>
/// The limit is generous: at least ten times the slowest legitimate run of the test in the loaded default set (the figures are in ADR-004 section 5a, "1b.9h"), so a loaded machine cannot make a good test fail by it.
/// </summary>
[AttributeUsage(AttributeTargets.Class | AttributeTargets.Method, AllowMultiple = false, Inherited = true)]
public sealed class HangLimitAttribute(int seconds) : PropertyAttribute(PropertyNames.Timeout, checked(seconds * 1000)), IApplyToContext, ITestAction
{
    private static readonly HangBreaker Shared = new();

    /// <summary>The limit in seconds.</summary>
    public int Seconds { get; } = seconds;

    public ActionTargets Targets => ActionTargets.Test;

    public void ApplyToContext(TestExecutionContext context)
    {
        ArgumentNullException.ThrowIfNull(context);
        context.TestCaseTimeout = checked(Seconds * 1000);
    }

    public void BeforeTest(ITest test)
    {
        ArgumentNullException.ThrowIfNull(test);
        if (Shared.Begin(test.FullName) is { } runaway)
        {
            Assert.Fail($"not started: the earlier test {runaway} did not come back within its limit and its thread is still running in this host, so another call could hang the same way");
        }
    }

    public void AfterTest(ITest test)
    {
        ArgumentNullException.ThrowIfNull(test);
        Shared.End(test.FullName);
    }
}
