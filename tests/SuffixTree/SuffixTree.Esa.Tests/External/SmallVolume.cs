using System.ComponentModel;
using System.Diagnostics;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// A real small volume for the tests of a full disk (work package 1b.9d): a tmpfs of a given size mounted in a directory of the test (the sandbox runs as root; a test that cannot mount one is ignored, as in
/// <c>MappedOutputFileTests</c>). On a tmpfs the page fault of a mapped page that has no block fails with a signal exactly as on a full disk, so a mapped write that lands on an unreserved page ends the test host:
/// that is the property the tests of this package pin.
/// </summary>
internal sealed class SmallVolume : IDisposable
{
    private readonly object _gate = new();
    private bool _filled;

    private SmallVolume(string path) => Path = path;

    /// <summary>The mount point.</summary>
    internal string Path { get; }

    /// <summary>A path inside the volume.</summary>
    internal string File(string name) => System.IO.Path.Combine(Path, name);

    /// <summary>Bytes available on the volume now.</summary>
    internal long Free => SystemReservationSource.Instance.AvailableFreeBytes(Path);

    /// <summary>Mounts a tmpfs of <paramref name="bytes"/> at <paramref name="mountPoint"/>; null where that is not possible (not Linux, not root, no <c>mount</c>).</summary>
    internal static SmallVolume? TryMount(string mountPoint, long bytes)
    {
        if (!OperatingSystem.IsLinux()) return null;
        Directory.CreateDirectory(mountPoint);
        try
        {
            using Process? p = Process.Start(new ProcessStartInfo("mount", ["-t", "tmpfs", "-o", $"size={bytes}", "tmpfs", mountPoint]) { RedirectStandardError = true, RedirectStandardOutput = true });
            if (p is null || !p.WaitForExit(10_000) || p.ExitCode != 0) return null;
            return new SmallVolume(mountPoint);
        }
        catch (Win32Exception)
        {
            return null;
        }
    }

    /// <summary>A volume or an ignored test.</summary>
    internal static SmallVolume MountOrIgnore(TempDir dir, string name, long bytes)
        => TryMount(dir.Sub(name), bytes) ?? throw new IgnoreException("could not mount a tmpfs here (needs Linux and root)");

    /// <summary>Writes a file of the name <c>ballast</c> until the volume refuses the next byte: no block is free afterwards.</summary>
    internal void FillToTheBrim()
    {
        lock (_gate)
        {
            if (_filled) return; // a second call, from another thread or another hook, finds the volume full already
            _filled = true;
            Fill();
        }
    }

    private void Fill()
    {
        using Microsoft.Win32.SafeHandles.SafeFileHandle handle = System.IO.File.OpenHandle(File("ballast"), FileMode.Create, FileAccess.Write, FileShare.None);
        var chunk = new byte[1 << 16];
        long offset = 0;
        foreach (int size in new[] { 1 << 16, 4096, 1 })
        {
            try
            {
                while (true)
                {
                    RandomAccess.Write(handle, chunk.AsSpan(0, size), offset);
                    offset += size;
                }
            }
            catch (IOException)
            {
                offset = RandomAccess.GetLength(handle); // a write that was cut short moved the end by part of its length
            }
        }
    }

    /// <summary>Deletes the ballast.</summary>
    internal void RemoveBallast()
    {
        lock (_gate)
        {
            System.IO.File.Delete(File("ballast"));
            _filled = false;
        }
    }

    public void Dispose()
    {
        for (int attempt = 0; attempt < 20; attempt++)
        {
            using Process? p = Process.Start(new ProcessStartInfo("umount", [Path]) { RedirectStandardError = true, RedirectStandardOutput = true });
            if (p is not null && p.WaitForExit(10_000) && p.ExitCode == 0) return;
            Thread.Sleep(50); // a handle that a finalizer has not released yet
        }
        using Process? lazy = Process.Start(new ProcessStartInfo("umount", ["-l", Path]) { RedirectStandardError = true, RedirectStandardOutput = true });
        lazy?.WaitForExit(10_000);
        Assert.Fail($"the volume at {Path} could not be unmounted: a file or a mapping of it is still open");
    }
}
