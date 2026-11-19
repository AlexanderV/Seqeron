using System.Runtime.InteropServices;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The bytes of a container in memory with an inaccessible page before and after the pages that hold them (an anonymous mapping and <c>mprotect</c>), so that a read one byte outside the container ends the process instead of
/// returning a neighbour's bytes: with <c>flush</c> the last byte of the image lies right against the page behind, otherwise the first byte right behind the page before. Linux only; elsewhere the image is held in a pinned array
/// (no guard, the harness reports it).
/// </summary>
internal sealed unsafe class GuardedImage : IDisposable
{
    private const int Page = 4096, ProtNone = 0, ProtRead = 1, ProtWrite = 2, MapPrivate = 2, MapAnonymous = 0x20;

    [DllImport("libc", SetLastError = true)]
    private static extern nint mmap(nint addr, nuint length, int prot, int flags, int fd, nint offset);

    [DllImport("libc", SetLastError = true)]
    private static extern int mprotect(nint addr, nuint length, int prot);

    [DllImport("libc", SetLastError = true)]
    private static extern int munmap(nint addr, nuint length);

    private readonly nint _base;
    private readonly nuint _size;
    private readonly GCHandle _pin;

    internal GuardedImage(byte[] image, bool flush)
    {
        Length = image.Length;
        if (!OperatingSystem.IsLinux())
        {
            _pin = GCHandle.Alloc(image, GCHandleType.Pinned);
            Data = (byte*)_pin.AddrOfPinnedObject();
            Guarded = false;
            return;
        }
        int pages = (image.Length + Page - 1) / Page;
        _size = (nuint)((pages + 2) * Page);
        _base = mmap(0, _size, ProtRead | ProtWrite, MapPrivate | MapAnonymous, -1, 0);
        if (_base == -1) throw new InvalidOperationException("mmap failed");
        byte* data = (byte*)_base + Page + (flush ? (pages * Page) - image.Length : 0);
        image.CopyTo(new Span<byte>(data, image.Length));
        if (mprotect(_base, Page, ProtNone) != 0 || mprotect(_base + (nint)((pages + 1) * Page), Page, ProtNone) != 0) throw new InvalidOperationException("mprotect failed");
        Data = data;
        Guarded = true;
    }

    /// <summary>The first byte of the image.</summary>
    internal byte* Data { get; }

    /// <summary>The length of the image.</summary>
    internal int Length { get; }

    /// <summary>True when the guard pages are in place.</summary>
    internal bool Guarded { get; }

    public void Dispose()
    {
        if (_base != 0) munmap(_base, _size);
        if (_pin.IsAllocated) _pin.Free();
    }
}
