using System.Runtime.InteropServices;

namespace SuffixTree.Persistent;

public sealed unsafe partial class MappedFileStorageProvider
{
    // ──────────────── Prefetch ────────────────

    /// <summary>
    /// Hints the OS to prefetch the mapped region into physical memory.
    /// Useful after loading a tree to avoid page faults on first access.
    /// This is advisory and fail-safe — errors are silently ignored.
    /// </summary>
    public void Prefetch()
    {
        ThrowIfDisposed();
        if (_ptr == null || _position <= 0) return;
        PrefetchRegion(_ptr, _position);
    }

    /// <summary>
    /// Hints the OS to prefetch the entire pre-allocated capacity before a build.
    /// Uses <c>PrefetchVirtualMemory</c> (Windows) / <c>posix_madvise</c> (Linux/macOS)
    /// to asynchronously start paging-in / zero-filling the mapped region.
    /// Advisory and fail-safe — returns immediately, errors are silently ignored.
    /// </summary>
    internal void PrefetchForBuild()
    {
        if (_ptr == null || _capacity <= 0) return;
        PrefetchRegion(_ptr, _capacity);
    }

#pragma warning disable CA1031 // Prefetch is advisory — must not throw
    private static void PrefetchRegion(byte* address, long length)
    {
        try
        {
            if (OperatingSystem.IsWindows())
                PrefetchWindows(address, length);
            else if (OperatingSystem.IsLinux() || OperatingSystem.IsMacOS())
                PrefetchPosix(address, length);
        }
        catch (Exception)
        {
            // Prefetch is advisory — ignore all errors
        }
    }
#pragma warning restore CA1031

    // ── Windows: PrefetchVirtualMemory ──

    private static void PrefetchWindows(byte* address, long length)
    {
        var entry = new WIN32_MEMORY_RANGE_ENTRY
        {
            VirtualAddress = address,
            NumberOfBytes = (nuint)length
        };
        PrefetchVirtualMemory(GetCurrentProcess(), 1, &entry, 0);
    }

    [StructLayout(LayoutKind.Sequential)]
    private struct WIN32_MEMORY_RANGE_ENTRY
    {
        public void* VirtualAddress;
        public nuint NumberOfBytes;
    }

    [LibraryImport("kernel32")]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    private static partial nint GetCurrentProcess();

    [LibraryImport("kernel32", SetLastError = true)]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static partial bool PrefetchVirtualMemory(
        nint hProcess,
        nuint numberOfEntries,
        WIN32_MEMORY_RANGE_ENTRY* virtualAddresses,
        uint flags);

    // ── Linux/macOS: posix_madvise ──

    private static void PrefetchPosix(byte* address, long length)
    {
        const int POSIX_MADV_WILLNEED = 3;
        _ = PosixMadvise(address, (nuint)length, POSIX_MADV_WILLNEED);
    }

    [LibraryImport("libc", EntryPoint = "posix_madvise")]
    [DefaultDllImportSearchPaths(DllImportSearchPath.System32)]
    private static partial int PosixMadvise(void* addr, nuint length, int advice);
}
