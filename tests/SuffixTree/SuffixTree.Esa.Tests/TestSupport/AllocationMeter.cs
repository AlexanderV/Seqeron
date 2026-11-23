using System.Globalization;
using System.Runtime;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>What <see cref="AllocationMeter.Measure"/> saw: the bytes that the calling thread allocated while the code ran, and whether the window was free of garbage collections.</summary>
internal readonly record struct AllocationReading(long Bytes, bool GcFree)
{
    /// <summary>
    /// Asserts that the code allocated less than <paramref name="limit"/> bytes. The counter can only be too high (a collection adds to it, see <see cref="AllocationMeter"/>), never too low by more than a few
    /// bytes, so a reading below the limit proves the code is within it whatever happened in the window. A reading above the limit is a failure when no collection ran in the window; when one did, the reading
    /// proves nothing and the test is inconclusive (it says why), it is never a failure of the code under test.
    /// </summary>
    internal void AssertLessThan(long limit, string what)
    {
        if (Bytes < limit) return;
        if (GcFree) Assert.Fail(string.Create(CultureInfo.InvariantCulture, $"{what}: the thread allocated {Bytes:N0} bytes in a window in which no garbage collection ran, the limit is {limit:N0}"));
        Assert.Inconclusive(string.Create(CultureInfo.InvariantCulture, $"{what}: the counter of the thread read {Bytes:N0} bytes against the limit {limit:N0}, in a window in which a garbage collection ran (the counter is not exact then), so the reading says nothing about the code"));
    }
}

/// <summary>
/// The bytes that a piece of code allocates, as the code allocates them (work package 1b.9h; the 1b.3 test <c>NoAllocationProportionalToTheTextForTheTotal_AndOnlyTheResultForTheLengths</c> was flaky under load).
/// <para>
/// <b>Why the plain counter is not enough.</b> <see cref="GC.GetAllocatedBytesForCurrentThread"/> is the thread's own counter, so the allocations of other threads are not in it, but it is not exact when a garbage collection runs
/// on another thread inside the window: the counter holds the whole allocation context that the thread was handed (8 KB) minus the part that is still unused, and when a <b>background</b> collection (a gen 2 collection that
/// another thread's large allocations start, the default of the concurrent collector) makes the context parsable the unused part stops being subtracted: the reading jumps up by up to one context, 8 KB, whatever the window does
/// (measured: a window that allocates nothing, with another thread that allocates large arrays flat out, read more than 1 KB too high in 10-20 % of the windows, up to +8,144 bytes, and never with the concurrent collector off,
/// <c>analysis/bench/allocprobe</c>). The test of 1b.3 allowed 1 KB. Four busy processes beside the test host make a window longer and a collection inside it likelier, so the test failed now and then.
/// </para>
/// <para>
/// <b>What this does.</b> The code runs in a no-GC region (<see cref="GC.TryStartNoGCRegion(long, long)"/>: the runtime promises that no collection, background or not, runs in it as long as the allocations of all threads stay
/// below the budget, 64 MB here, 32 MB of it for large objects; the start itself collects, so the thread's allocation context is empty at the first reading), and the reading says whether the window stayed free of collections.
/// There is no retry: one window, one reading, and a reading that could not be taken cleanly is reported as such, never as a pass or as a failure of the code (<see cref="AllocationReading.AssertLessThan"/>).
/// </para>
/// </summary>
internal static class AllocationMeter
{
    /// <summary>The budget of the no-GC region: the allocations of every thread together in the window, 64 MB, of which 32 MB for large objects.</summary>
    internal const long RegionBytes = 64L << 20, RegionLargeObjectBytes = 32L << 20;

    /// <summary>
    /// Runs <paramref name="code"/> once and returns the bytes that the calling thread allocated meanwhile. Allocate what the code needs (the delegate, its captured variables) before the call: they are not in the window.
    /// A no-GC region that is already open (the caller's, or another thread's: the mode is the process's) is used as it is and left open, a failed start would end it; the window is then free of collections as long as it holds.
    /// </summary>
    internal static AllocationReading Measure(Action code)
    {
        ArgumentNullException.ThrowIfNull(code);
        bool alreadyOpen = GCSettings.LatencyMode == GCLatencyMode.NoGCRegion;
        bool opened = !alreadyOpen && TryStartRegion();
        try
        {
            int collections = GC.CollectionCount(0);                // every collection counts in generation 0, whatever its generation
            long before = GC.GetAllocatedBytesForCurrentThread();
            code();
            long after = GC.GetAllocatedBytesForCurrentThread();
            bool gcFree = (alreadyOpen || opened) && GCSettings.LatencyMode == GCLatencyMode.NoGCRegion && GC.CollectionCount(0) == collections;
            return new AllocationReading(after - before, gcFree);
        }
        finally
        {
            if (opened) EndRegion();
        }
    }

    private static bool TryStartRegion()
    {
        try
        {
            return GC.TryStartNoGCRegion(RegionBytes, RegionLargeObjectBytes);
        }
        catch (InvalidOperationException)           // another thread opened a region between the check and the start: the window is measured without the promise (and that region is gone: a failed start ends it)
        {
            return false;
        }
        catch (ArgumentOutOfRangeException)         // the runtime cannot hold the budget
        {
            return false;
        }
    }

    private static void EndRegion()
    {
        try
        {
            if (GCSettings.LatencyMode == GCLatencyMode.NoGCRegion) GC.EndNoGCRegion();
        }
        catch (InvalidOperationException)           // the region was left by a collection (the budget was exceeded or a collection was induced) between the check and the call: nothing is left to end
        {
        }
    }
}
