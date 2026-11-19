using System.Buffers.Binary;
using Microsoft.Win32.SafeHandles;
using SuffixTree.Persistent;

namespace SuffixTree.Esa;

/// <summary>
/// The one entry point that opens a suffix index file of either kind (ADR-003 section 8): a SESA container (this engine) or a v6 persistent suffix tree (<see cref="PersistentSuffixTreeFactory"/>), told apart by the
/// first eight bytes of the file.
/// </summary>
public static class SuffixIndexFactory
{
    /// <summary>Opens the index file at <paramref name="path"/>.</summary>
    /// <param name="path">The file.</param>
    /// <param name="options">Validation level and access mode of a SESA container (default <see cref="EsaOpenOptions"/>); a v6 file ignores them.</param>
    /// <returns>The index as <see cref="ISuffixTree"/>; it is <see cref="IDisposable"/> and holds the file mapped until it is disposed.</returns>
    /// <exception cref="ArgumentException"><paramref name="path"/> is null or empty.</exception>
    /// <exception cref="InvalidDataException">The file is of neither kind (<c>unknown suffix index kind</c>), or a SESA container is corrupt at the requested level.</exception>
    /// <exception cref="NotSupportedException">A SESA container of a larger major version, or one that needs a feature this reader does not have.</exception>
    public static ISuffixTree Open(string path, EsaOpenOptions? options = null)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        Span<byte> magic = stackalloc byte[8];
        int read;
        using (SafeFileHandle handle = File.OpenHandle(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete))
        {
            // the decision reads the first eight bytes through the handle: nothing is mapped for a file of an unknown kind (a device or a pipe that cannot be read at an offset is refused, not waited for)
            if (handle.IsInvalid) throw new IOException($"cannot open '{path}'");
            read = ReadMagic(handle, magic);
        }
        if (read == 8 && magic.SequenceEqual(External.SesaFormat.Magic)) return EsaSuffixTreeFactory.Open(path, options);
        if (read == 8 && BinaryPrimitives.ReadInt64LittleEndian(magic) == V6Magic) return PersistentSuffixTreeFactory.LoadPersistent(path);
        throw new InvalidDataException("unknown suffix index kind");
    }

    /// <summary>The v6 magic: <c>PersistentConstants.MAGIC_NUMBER</c> stored as a little-endian <c>int64</c> (the bytes <c>46 46 49 58 45 45 52 54</c>).</summary>
    private const long V6Magic = 0x5452454558494646;

    private static int ReadMagic(SafeFileHandle handle, Span<byte> magic)
    {
        try
        {
            return RandomAccess.Read(handle, magic, 0);
        }
        catch (NotSupportedException)
        {
            throw new InvalidDataException("unknown suffix index kind: the file cannot be read at an offset (not a regular file)");
        }
    }
}
