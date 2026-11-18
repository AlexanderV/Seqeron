#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>Where the kernel's text files (<c>/proc/self/io</c>, <c>/proc/self/status</c>) come from; injectable so that a platform without <c>/proc</c> can be simulated.</summary>
internal interface IProcSource
{
    /// <summary>Reads the start of the file at <paramref name="path"/> into <paramref name="buffer"/>; returns the bytes read, or -1 when the file does not exist or cannot be read. Never throws.</summary>
    int TryRead(string path, Span<byte> buffer);
}

/// <summary>The real file system.</summary>
internal sealed class FileProcSource : IProcSource
{
    internal static FileProcSource Instance { get; } = new();

    public int TryRead(string path, Span<byte> buffer) => throw new NotImplementedException();
}

/// <summary>Parser of the <c>Key: value</c> lines of the kernel's text files.</summary>
internal static class ProcText
{
    /// <summary>Finds the line that starts with <paramref name="key"/> (including its colon) and parses the unsigned decimal number after the blanks.</summary>
    internal static bool TryGetLong(ReadOnlySpan<byte> text, ReadOnlySpan<byte> key, out long value) => throw new NotImplementedException();
}

/// <summary>A snapshot of <c>/proc/self/io</c>: syscall-level byte counts and block-layer byte counts.</summary>
internal readonly record struct ProcIoSnapshot(long RChar, long WChar, long SysCr, long SysCw, long ReadBytes, long WriteBytes)
{
    /// <summary>Field-wise difference <c>this - earlier</c>.</summary>
    internal ProcIoSnapshot Since(ProcIoSnapshot earlier) => throw new NotImplementedException();

    /// <summary>Reads <c>/proc/self/io</c> through <paramref name="source"/> (the real files by default); null where it is unavailable.</summary>
    internal static ProcIoSnapshot? TryRead(IProcSource? source = null) => throw new NotImplementedException();
}

/// <summary>What <see cref="IoMeasurement.Finish"/> reports: the explicit counters of the interval and the <c>/proc/self/io</c> delta (null = unavailable).</summary>
internal readonly record struct IoReport(IoSnapshot Counted, ProcIoSnapshot? Proc);

/// <summary>One measured interval: the counters and <c>/proc/self/io</c> are read at construction and again at <see cref="Finish"/>.</summary>
internal sealed class IoMeasurement
{
    internal IoMeasurement(IoCounters counters, IProcSource? source = null) => throw new NotImplementedException();

    internal IoReport Finish() => throw new NotImplementedException();
}
