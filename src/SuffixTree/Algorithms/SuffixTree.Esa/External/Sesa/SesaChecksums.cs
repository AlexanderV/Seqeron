using System.Diagnostics;
using System.IO.MemoryMappedFiles;

namespace SuffixTree.Esa.External;

/// <summary>
/// The CRC-32C pass of PH7 (ADR-004 section 4.2, "one extra read pass"): every section is read once, sequentially, through a read-only mapping of the finished file in views of at most <see cref="ViewBytes"/>, so that no
/// anonymous memory is used (the pages are file-backed and leave the process as the kernel drops them). The sections are done one after the other; the figure of each (bytes, seconds, hardware or table path) goes to the result.
/// </summary>
internal static class SesaChecksums
{
    /// <summary>Largest view of the pass.</summary>
    internal const int ViewBytes = 16 << 20;

    /// <summary>The CRC-32C of the payloads of <paramref name="ranges"/> of the file <paramref name="path"/>, in the order given, with the bytes and the seconds of each.</summary>
    internal static (uint Crc, SesaCrcRow Row)[] Compute(string path, IReadOnlyList<(int Kind, long Offset, long Length)> ranges, int viewBytes = ViewBytes)
    {
        ArgumentNullException.ThrowIfNull(path);
        ArgumentNullException.ThrowIfNull(ranges);
        ArgumentOutOfRangeException.ThrowIfLessThan(viewBytes, 1);
        var result = new (uint, SesaCrcRow)[ranges.Count];
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        long fileLength = stream.Length;
        MemoryMappedFile? map = fileLength > 0 ? MemoryMappedFile.CreateFromFile(stream, null, 0, MemoryMappedFileAccess.Read, HandleInheritability.None, leaveOpen: true) : null;
        try
        {
            for (int i = 0; i < ranges.Count; i++)
            {
                (int kind, long offset, long length) = ranges[i];
                if (offset < 0 || length < 0 || offset > fileLength - length) throw new InvalidDataException($"the section {kind} [{offset}, {offset + length}) lies outside the file of {fileLength} bytes");
                long t0 = Stopwatch.GetTimestamp();
                uint crc = 0;
                for (long done = 0; done < length; done += viewBytes)
                {
                    long take = Math.Min(viewBytes, length - done);
                    crc = ChecksumView(map!, offset + done, take, crc);
                }
                result[i] = (crc, new SesaCrcRow(kind, length, (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, Crc32C.HardwareAccelerated));
            }
        }
        finally
        {
            map?.Dispose();
        }
        return result;
    }

    private static unsafe uint ChecksumView(MemoryMappedFile map, long offset, long length, uint crc)
    {
        using MemoryMappedViewAccessor view = map.CreateViewAccessor(offset, length, MemoryMappedFileAccess.Read);
        byte* p = null;
        view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        try
        {
            return Crc32C.Append(crc, p + view.PointerOffset, length);
        }
        finally
        {
            view.SafeMemoryMappedViewHandle.ReleasePointer();
        }
    }
}
