using System.Runtime.ExceptionServices;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// A time limit for the code that a test source runs at TEST DISCOVERY (work package 1b.9h): NUnit evaluates <c>[TestCaseSource]</c> members when it builds the tests of the assembly, in every run and before any test, so a source that calls
/// the code under test is outside every limit per test: one loop in the code under test that never ends stops the discovery of the whole assembly, and no test runs at all (the 14 mutants of 1b.9c that hang the builder did exactly
/// that once <c>SesaResume1b9e1Tests.GapCases</c> built a reference container at discovery). The code runs on a thread of its own; when it does not return in time the source gives one case that fails and names the cause (the thread is
/// left running, as in <see cref="HangLimitAttribute"/>).
/// </summary>
internal static class DiscoveryLimit
{
    /// <summary>Runs <paramref name="compute"/> on a background thread. True with its value when it returns within <paramref name="limit"/> (an exception that it throws is thrown again here), false when it does not.</summary>
    internal static bool TryRun<T>(TimeSpan limit, Func<T> compute, out T? value)
    {
        ArgumentNullException.ThrowIfNull(compute);
        T? result = default;
        ExceptionDispatchInfo? error = null;
        using var done = new ManualResetEventSlim();
        var thread = new Thread(() =>
        {
            try
            {
                result = compute();
            }
#pragma warning disable CA1031 // the failure of the computation is given to the caller of TryRun, with its own type and stack
            catch (Exception e)
            {
                error = ExceptionDispatchInfo.Capture(e);
            }
#pragma warning restore CA1031
            finally
            {
                done.Set();
            }
        })
        { IsBackground = true, Name = "discovery-limit" };
        thread.Start();
        if (!done.Wait(limit))
        {
            value = default;
            return false;
        }
        error?.Throw();
        value = result;
        return true;
    }
}
