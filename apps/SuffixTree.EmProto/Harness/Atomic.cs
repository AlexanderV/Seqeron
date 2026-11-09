namespace SuffixTree.EmProto.Harness;

/// <summary>Small lock-free helpers for the counters of the harness.</summary>
internal static class Atomic
{
    /// <summary>Sets <paramref name="target"/> to <paramref name="value"/> when that is larger (a running maximum).</summary>
    internal static void RaiseTo(ref long target, long value)
    {
        long seen = Interlocked.Read(ref target);
        while (seen < value)
        {
            long previous = Interlocked.CompareExchange(ref target, value, seen);
            if (previous == seen) return;
            seen = previous;
        }
    }
}
