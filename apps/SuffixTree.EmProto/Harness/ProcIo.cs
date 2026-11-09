using System.Globalization;

namespace SuffixTree.EmProto.Harness;

/// <summary>A snapshot of <c>/proc/self/io</c>: syscall-level byte counts (<c>rchar</c>, <c>wchar</c>) and block-layer byte counts (<c>read_bytes</c>, <c>write_bytes</c>).</summary>
internal readonly record struct ProcIoSnapshot(long RChar, long WChar, long SysCr, long SysCw, long ReadBytes, long WriteBytes)
{
    /// <summary>Field-wise difference <c>this - earlier</c>.</summary>
    internal ProcIoSnapshot Since(ProcIoSnapshot earlier)
        => new(RChar - earlier.RChar, WChar - earlier.WChar, SysCr - earlier.SysCr, SysCw - earlier.SysCw, ReadBytes - earlier.ReadBytes, WriteBytes - earlier.WriteBytes);

    /// <summary>Reads <c>/proc/self/io</c>; all zero where the file does not exist.</summary>
    internal static ProcIoSnapshot Read()
    {
        const string path = "/proc/self/io";
        if (!File.Exists(path)) return default;
        long rchar = 0, wchar = 0, syscr = 0, syscw = 0, rb = 0, wb = 0;
        foreach (string line in File.ReadAllLines(path))
        {
            int colon = line.IndexOf(':', StringComparison.Ordinal);
            if (colon < 0) continue;
            long value = long.Parse(line.AsSpan(colon + 1).Trim(), CultureInfo.InvariantCulture);
            switch (line[..colon])
            {
                case "rchar": rchar = value; break;
                case "wchar": wchar = value; break;
                case "syscr": syscr = value; break;
                case "syscw": syscw = value; break;
                case "read_bytes": rb = value; break;
                case "write_bytes": wb = value; break;
                default: break;
            }
        }
        return new ProcIoSnapshot(rchar, wchar, syscr, syscw, rb, wb);
    }
}
