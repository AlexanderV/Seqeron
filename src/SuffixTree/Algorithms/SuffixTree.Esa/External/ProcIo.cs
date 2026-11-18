using System.Security;
using Microsoft.Win32.SafeHandles;

namespace SuffixTree.Esa.External;

/// <summary>
/// Where the kernel's text files (<c>/proc/self/io</c>, <c>/proc/self/status</c>) come from. Injectable so that a platform without <c>/proc</c> (macOS, Windows) can be simulated in a test, and so that the
/// counters degrade to "unavailable" there instead of throwing (the library runs on Windows too; the owner develops on macOS).
/// </summary>
internal interface IProcSource
{
    /// <summary>Reads the start of the file at <paramref name="path"/> into <paramref name="buffer"/>; returns the bytes read, or -1 when the file does not exist or cannot be read. Never throws.</summary>
    int TryRead(string path, Span<byte> buffer);
}

/// <summary>The real file system.</summary>
internal sealed class FileProcSource : IProcSource
{
    internal static FileProcSource Instance { get; } = new();

    public int TryRead(string path, Span<byte> buffer)
    {
        try
        {
            using SafeFileHandle handle = File.OpenHandle(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            int total = 0;
            while (total < buffer.Length)
            {
                int got = RandomAccess.Read(handle, buffer[total..], total);
                if (got <= 0) break;
                total += got;
            }
            return total;
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException or NotSupportedException or ArgumentException or SecurityException)
        {
            return -1;
        }
    }
}

/// <summary>Parser of the <c>Key: value</c> lines of the kernel's text files, and the guarded read of a source.</summary>
internal static class ProcText
{
    /// <summary>Reads through <paramref name="source"/>; any failure of the source (it should never throw, a hostile one might) is "unavailable" (-1); a count above the buffer is clamped.</summary>
    internal static int SafeRead(IProcSource source, string path, Span<byte> buffer)
    {
        try
        {
            int got = source.TryRead(path, buffer);
            return got < 0 ? -1 : Math.Min(got, buffer.Length);
        }
        catch (Exception ex) when (ex is IOException or UnauthorizedAccessException or NotSupportedException or ArgumentException or InvalidOperationException or SecurityException)
        {
            return -1;
        }
    }

    /// <summary>
    /// Finds the first line that starts with <paramref name="key"/> (including its colon) and parses the unsigned decimal number after the blanks that follow it. The key must stand at the start of a line, so that
    /// <c>write_bytes:</c> does not match inside <c>cancelled_write_bytes:</c>. False when there is no such line, no digit follows, or the number does not fit a <see cref="long"/>.
    /// </summary>
    internal static bool TryGetLong(ReadOnlySpan<byte> text, ReadOnlySpan<byte> key, out long value)
    {
        value = 0;
        if (key.IsEmpty) return false;
        int from = 0;
        while (from <= text.Length - key.Length)
        {
            int at = text[from..].IndexOf(key);
            if (at < 0) return false;
            int position = from + at;
            if (position == 0 || text[position - 1] == (byte)'\n') return ParseAfter(text, position + key.Length, out value);
            from = position + 1;
        }
        return false;
    }

    private static bool ParseAfter(ReadOnlySpan<byte> text, int start, out long value)
    {
        value = 0;
        int i = start;
        while (i < text.Length && (text[i] == (byte)' ' || text[i] == (byte)'\t')) i++;
        long parsed = 0;
        int digits = 0;
        while (i < text.Length && text[i] >= (byte)'0' && text[i] <= (byte)'9')
        {
            int digit = text[i] - (byte)'0';
            if (parsed > (long.MaxValue - digit) / 10) return false;
            parsed = (parsed * 10) + digit;
            digits++;
            i++;
        }
        if (digits == 0) return false;
        value = parsed;
        return true;
    }
}

/// <summary>A snapshot of <c>/proc/self/io</c>: syscall-level byte counts (<c>rchar</c>, <c>wchar</c>) and block-layer byte counts (<c>read_bytes</c>, <c>write_bytes</c>).</summary>
internal readonly record struct ProcIoSnapshot(long RChar, long WChar, long SysCr, long SysCw, long ReadBytes, long WriteBytes)
{
    private const string Path = "/proc/self/io";

    /// <summary>Field-wise difference <c>this - earlier</c>.</summary>
    internal ProcIoSnapshot Since(ProcIoSnapshot earlier)
        => new(RChar - earlier.RChar, WChar - earlier.WChar, SysCr - earlier.SysCr, SysCw - earlier.SysCw, ReadBytes - earlier.ReadBytes, WriteBytes - earlier.WriteBytes);

    /// <summary>
    /// Reads <c>/proc/self/io</c> through <paramref name="source"/> (the real files by default). Null where the file does not exist or any of the six fields is missing or unreadable ("unavailable":
    /// a platform without <c>/proc</c>, a kernel without I/O accounting, a hostile source); never throws.
    /// </summary>
    internal static ProcIoSnapshot? TryRead(IProcSource? source = null)
    {
        Span<byte> buffer = stackalloc byte[2048];
        int got = ProcText.SafeRead(source ?? FileProcSource.Instance, Path, buffer);
        if (got <= 0) return null;
        ReadOnlySpan<byte> text = buffer[..got];
        if (!ProcText.TryGetLong(text, "rchar:"u8, out long rchar)
            || !ProcText.TryGetLong(text, "wchar:"u8, out long wchar)
            || !ProcText.TryGetLong(text, "syscr:"u8, out long syscr)
            || !ProcText.TryGetLong(text, "syscw:"u8, out long syscw)
            || !ProcText.TryGetLong(text, "read_bytes:"u8, out long readBytes)
            || !ProcText.TryGetLong(text, "write_bytes:"u8, out long writeBytes))
        {
            return null;
        }
        return new ProcIoSnapshot(rchar, wchar, syscr, syscw, readBytes, writeBytes);
    }
}

/// <summary>What <see cref="IoMeasurement.Finish"/> reports: the explicit counters of the interval and the <c>/proc/self/io</c> delta (null = unavailable at either end).</summary>
internal readonly record struct IoReport(IoSnapshot Counted, ProcIoSnapshot? Proc);

/// <summary>One measured interval: the counters and <c>/proc/self/io</c> are read at construction and again at <see cref="Finish"/>.</summary>
internal sealed class IoMeasurement
{
    private readonly IoCounters _counters;
    private readonly IProcSource? _source;
    private readonly IoSnapshot _start;
    private readonly ProcIoSnapshot? _procStart;

    internal IoMeasurement(IoCounters counters, IProcSource? source = null)
    {
        ArgumentNullException.ThrowIfNull(counters);
        _counters = counters;
        _source = source;
        _start = counters.Snapshot();
        _procStart = ProcIoSnapshot.TryRead(source);
    }

    internal IoReport Finish()
    {
        IoSnapshot counted = _counters.Snapshot().Since(_start);
        ProcIoSnapshot? end = ProcIoSnapshot.TryRead(_source);
        ProcIoSnapshot? proc = _procStart is { } begin && end is { } last ? last.Since(begin) : null;
        return new IoReport(counted, proc);
    }
}
