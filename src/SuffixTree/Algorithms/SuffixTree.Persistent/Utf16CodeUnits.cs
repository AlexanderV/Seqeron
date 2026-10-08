using System.Buffers.Binary;
using System.Runtime.InteropServices;

namespace SuffixTree.Persistent;

/// <summary>
/// Lossless UTF-16LE code-unit codec for persisted and exported text.
/// <para>
/// A suffix tree indexes <c>char</c> code units, so any string — including lone
/// surrogates and surrogate pairs that straddle a chunk boundary — is a valid text.
/// <see cref="System.Text.Encoding.Unicode"/> would replace unpaired surrogates with
/// U+FFFD, silently changing the stored text relative to the tree built over it.
/// These helpers copy code units verbatim (byte-identical to <c>Encoding.Unicode</c>
/// for well-formed UTF-16, so existing files and exports are unaffected).
/// </para>
/// </summary>
internal static class Utf16CodeUnits
{
    /// <summary>Writes <paramref name="chars"/> as little-endian code units; returns the byte count.</summary>
    public static int Write(ReadOnlySpan<char> chars, Span<byte> destination)
    {
        int byteCount = chars.Length * sizeof(char);
        if (destination.Length < byteCount)
            throw new ArgumentException("Destination buffer is too small.", nameof(destination));
        for (int i = 0; i < chars.Length; i++)
            BinaryPrimitives.WriteUInt16LittleEndian(destination.Slice(i * sizeof(char)), chars[i]);
        return byteCount;
    }

    /// <summary>Reads little-endian code units verbatim into a string.</summary>
    public static string Read(ReadOnlySpan<byte> bytes)
    {
        if ((bytes.Length & 1) != 0)
            throw new ArgumentException("UTF-16 payload must have an even byte length.", nameof(bytes));
        if (BitConverter.IsLittleEndian)
            return new string(MemoryMarshal.Cast<byte, char>(bytes));
        var chars = new char[bytes.Length / sizeof(char)];
        for (int i = 0; i < chars.Length; i++)
            chars[i] = (char)BinaryPrimitives.ReadUInt16LittleEndian(bytes.Slice(i * sizeof(char)));
        return new string(chars);
    }
}
