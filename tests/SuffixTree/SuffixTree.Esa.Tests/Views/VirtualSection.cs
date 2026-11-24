using System.Runtime.InteropServices;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// A range of virtual address space of any size (an anonymous mapping that is inaccessible and reserves no memory) of which a test makes a few ranges readable and writable: the pointer readers and the views over them can then be run at
/// element indexes above 2^31 and 2^32 (a section of 2^33 five-byte entries is 43 GB of address space) with the memory of the few pages that the test touches, and a read of a page that the test did not open ends the test host (as a read
/// outside a section must). Linux only (the tests that need it are inconclusive elsewhere).
/// </summary>
internal sealed unsafe class VirtualSection : IDisposable
{
    private const int ProtNone = 0, ProtRead = 1, ProtWrite = 2, MapPrivate = 2, MapAnonymous = 0x20, MapNoReserve = 0x4000;
    private const long Page = 4096;

    [DllImport("libc", SetLastError = true)]
    private static extern nint mmap(nint addr, nuint length, int prot, int flags, int fd, nint offset);

    [DllImport("libc", SetLastError = true)]
    private static extern int mprotect(nint addr, nuint length, int prot);

    [DllImport("libc", SetLastError = true)]
    private static extern int munmap(nint addr, nuint length);

    private readonly nint _base;
    private readonly nuint _size;

    /// <summary>Reserves <paramref name="bytes"/> bytes (rounded up to pages) of address space; every byte is inaccessible until <see cref="Open"/>.</summary>
    internal VirtualSection(long bytes)
    {
        if (!OperatingSystem.IsLinux()) Assert.Inconclusive("the virtual section needs mmap and mprotect (Linux)");
        _size = (nuint)(((bytes + Page - 1) / Page) * Page);
        _base = mmap(0, _size, ProtNone, MapPrivate | MapAnonymous | MapNoReserve, -1, 0);
        if (_base == -1) Assert.Inconclusive($"mmap of {bytes} bytes of address space failed: {Marshal.GetLastWin32Error()}");
        Length = bytes;
    }

    /// <summary>The first byte of the section.</summary>
    internal byte* Data => (byte*)_base;

    /// <summary>The size that was asked for.</summary>
    internal long Length { get; }

    /// <summary>Makes the bytes <c>[offset, offset + length)</c> (rounded out to pages) readable and writable; they read as zero until written.</summary>
    internal void Open(long offset, long length)
    {
        long first = offset / Page * Page;
        long last = (offset + length + Page - 1) / Page * Page;
        if (first < 0 || last > (long)_size) throw new ArgumentOutOfRangeException(nameof(offset));
        if (mprotect(_base + (nint)first, (nuint)(last - first), ProtRead | ProtWrite) != 0) throw new InvalidOperationException($"mprotect failed: {Marshal.GetLastWin32Error()}");
    }

    /// <summary>Opens the first <paramref name="head"/> and the last <paramref name="tail"/> bytes (the ends of a section that a closed form fills).</summary>
    internal void OpenEnds(long head, long tail)
    {
        Open(0, head);
        Open(Length - tail, tail);
    }

    public void Dispose() => munmap(_base, _size);
}
