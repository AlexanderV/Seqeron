using System.Diagnostics;
using System.Globalization;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>The filesystem does not make sparse files (or the volume is too full to try): a test that needs one is inconclusive, never a failure and never a dense write.</summary>
public sealed class SparseUnavailableException(string message) : Exception(message);

/// <summary>What the probe of a directory found.</summary>
/// <param name="Sparse">True when a file of <paramref name="ProbeBytes"/> bytes made with <c>SetLength</c> took less than <see cref="LargeDisk.MaxAllocatedByProbe"/> bytes of disk.</param>
/// <param name="AllocatedBytes">Bytes of disk that the probe file took (blocks times the block size).</param>
/// <param name="ProbeBytes">Logical size of the probe file.</param>
/// <param name="Message">One line that says what was found and why the answer is what it is.</param>
public sealed record SparseProbe(bool Sparse, long AllocatedBytes, long ProbeBytes, string Message);

/// <summary>
/// The disk rules of the <c>Large</c> category (work package S2.1; the volume has about 17 GB free and 252 GB nominal): the directory is <c>$SEQERON_LARGE_DIR</c> or a subdirectory of the system temporary directory; a sparse
/// file is made only after a probe has shown that the filesystem makes sparse files (<c>SetLength</c> of 1 GiB, then the blocks that the file occupies, by <c>stat</c>); no dense region above 64 MiB is ever written; the
/// nominal size of one container is capped; free space is checked before any file is made.
/// </summary>
public static class LargeDisk
{
    /// <summary>The most bytes of real data that one container may receive (the patches): never a dense region above 64 MiB.</summary>
    public const long MaxDenseBytes = 64L << 20;

    /// <summary>The largest nominal (logical) size of one sparse file: 200 GB, below the 252 GB of the volume.</summary>
    public const long MaxNominalBytes = 200L << 30;

    /// <summary>Logical size of the probe file.</summary>
    public const long ProbeBytes = 1L << 30;

    /// <summary>A probe file that takes at least this much disk is not sparse.</summary>
    public const long MaxAllocatedByProbe = 1L << 20;

    /// <summary>Free space that must remain for a sparse file to be made at all: a filesystem that turns out not to be sparse must not be able to fill the volume.</summary>
    public const long MinimumFreeBytes = 2L << 30;

    private static readonly Dictionary<string, SparseProbe> Probed = [];

    /// <summary>The directory of the sparse files: <c>$SEQERON_LARGE_DIR</c>, else <c>seqeron-large</c> under the system temporary directory.</summary>
    public static string Directory
    {
        get
        {
            string? configured = Environment.GetEnvironmentVariable("SEQERON_LARGE_DIR");
            return string.IsNullOrWhiteSpace(configured) ? Path.Combine(Path.GetTempPath(), "seqeron-large") : configured;
        }
    }

    /// <summary>A new, unused file path in <paramref name="directory"/> (the directory is created, the file is not).</summary>
    public static string NewPath(string directory, string stem)
    {
        System.IO.Directory.CreateDirectory(directory);
        return Path.Combine(directory, $"{stem}-{Guid.NewGuid():N}.sesa");
    }

    /// <summary>The decision of the probe from what it measured: sparse when the file of <paramref name="probeBytes"/> bytes took less than <see cref="MaxAllocatedByProbe"/> bytes.</summary>
    public static SparseProbe Evaluate(long allocatedBytes, long probeBytes)
    {
        bool sparse = allocatedBytes >= 0 && allocatedBytes < MaxAllocatedByProbe;
        string message = sparse
            ? $"sparse files are made here: a file of {probeBytes} bytes took {allocatedBytes} bytes of disk"
            : $"sparse files are NOT made here: a file of {probeBytes} bytes took {allocatedBytes} bytes of disk (at least {MaxAllocatedByProbe} counts as dense); no sparse container is written";
        return new SparseProbe(sparse, allocatedBytes, probeBytes, message);
    }

    /// <summary>
    /// Makes a file of <paramref name="probeBytes"/> bytes with <c>SetLength</c> in <paramref name="directory"/>, measures what it occupies and deletes it. The answer is remembered per directory and size.
    /// </summary>
    public static SparseProbe Probe(string directory, long probeBytes = ProbeBytes)
    {
        string key = Path.GetFullPath(directory) + "|" + probeBytes.ToString(CultureInfo.InvariantCulture);
        lock (Probed)
        {
            if (Probed.TryGetValue(key, out SparseProbe? known)) return known;
            System.IO.Directory.CreateDirectory(directory);
            long free = FreeBytes(directory);
            if (free < MinimumFreeBytes)
                return Remember(key, new SparseProbe(false, -1, probeBytes, $"the volume has {free} bytes free, less than the {MinimumFreeBytes} that a sparse test needs to be safe; nothing is written"));
            string path = Path.Combine(directory, $"probe-{Guid.NewGuid():N}.sparse");
            try
            {
                using (var stream = new FileStream(path, FileMode.CreateNew, FileAccess.Write, FileShare.None))
                {
                    stream.SetLength(probeBytes);
                    stream.Flush(true);
                }
                return Remember(key, Evaluate(AllocatedBytes(path), probeBytes));
            }
            finally
            {
                File.Delete(path);
            }
        }
    }

    /// <summary>Throws <see cref="SparseUnavailableException"/> unless the filesystem of <paramref name="directory"/> makes sparse files (see <see cref="Probe"/>).</summary>
    public static void RequireSparse(string directory)
    {
        SparseProbe probe = Probe(directory);
        if (!probe.Sparse) throw new SparseUnavailableException(probe.Message);
    }

    /// <summary>Bytes available to this user on the volume of <paramref name="path"/>.</summary>
    public static long FreeBytes(string path) => new DriveInfo(Path.GetFullPath(path)).AvailableFreeSpace;

    /// <summary>
    /// Bytes of disk that the file occupies (the blocks that <c>stat</c> reports times the block size), so a hole is not counted; -1 when <c>stat</c> cannot be run (the probe then says "not sparse" and the test is inconclusive).
    /// </summary>
    public static long AllocatedBytes(string path)
    {
        try
        {
            var info = new ProcessStartInfo("stat") { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
            info.ArgumentList.Add("-c");
            info.ArgumentList.Add("%b %B");
            info.ArgumentList.Add(path);
            using Process process = Process.Start(info) ?? throw new InvalidOperationException("stat did not start");
            string output = process.StandardOutput.ReadToEnd();
            process.WaitForExit();
            if (process.ExitCode != 0) return -1;
            string[] parts = output.Split(' ', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            return parts.Length == 2 ? long.Parse(parts[0], CultureInfo.InvariantCulture) * long.Parse(parts[1], CultureInfo.InvariantCulture) : -1;
        }
        catch (Exception e) when (e is System.ComponentModel.Win32Exception or InvalidOperationException or FormatException)
        {
            return -1;
        }
    }

    private static SparseProbe Remember(string key, SparseProbe probe)
    {
        Probed[key] = probe;
        return probe;
    }
}
