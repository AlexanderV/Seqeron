using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Build;

/// <summary>Bytes of scratch files alive at a time and their peak (ADR-004 Q1: peak scratch in B/char is reported next to the I/O volume).</summary>
internal static class ScratchTracker
{
    private static long s_live;
    private static long s_peak;

    internal static long Live => Interlocked.Read(ref s_live);

    internal static long Peak => Interlocked.Read(ref s_peak);

    internal static void Reset()
    {
        Interlocked.Exchange(ref s_live, 0);
        Interlocked.Exchange(ref s_peak, 0);
    }

    internal static void Add(long bytes) => Atomic.RaiseTo(ref s_peak, Interlocked.Add(ref s_live, bytes));

    internal static void Remove(long bytes) => Interlocked.Add(ref s_live, -bytes);
}
