#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>Small lock-free helpers for the counters of the external-memory building blocks.</summary>
internal static class Atomic
{
    /// <summary>Sets <paramref name="target"/> to <paramref name="value"/> when that is larger (a running maximum).</summary>
    internal static void RaiseTo(ref long target, long value) => throw new NotImplementedException();
}
