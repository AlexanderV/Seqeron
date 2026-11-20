#if EXPERIMENTS
using System.Runtime.InteropServices;
using SuffixTree;
using SuffixTree.Esa;
using SuffixTree.Esa.External;

namespace ReadBench;

/// <summary>
/// Isolated experiments of the diagnosis of 1b.9a (internals through the friend assembly name). Each variant changes one factor of the file-backed read path and nothing else.
/// resident-poh     the container read into a pinned managed byte[] (the Pinned Object Heap), the same mapped-arm code over it (SesaFile.FromMemory)
/// resident-native  the container read into NativeMemory.AlignedAlloc (4 KB aligned anonymous memory), the same code
/// memory-packed    the in-memory tree with the LCP-disambiguated child table (the layout of the container's CHILD codec 1) instead of the default tagged one: the layout alone
/// </summary>
internal static unsafe class Experiments
{
    private static byte[]? s_pinned;
    private static void* s_native;

    internal static ISuffixTree Create(string variant, string file, string text)
    {
        var options = new EsaOpenOptions { Access = EsaAccessMode.Checked };
        if (variant.EndsWith("-trusted", StringComparison.Ordinal))
        {
            options = new EsaOpenOptions { Access = EsaAccessMode.Trusted };
            variant = variant[..^"-trusted".Length];
        }
        switch (variant)
        {
            case "resident-poh":
                {
                    long length = new FileInfo(file).Length;
                    s_pinned = GC.AllocateUninitializedArray<byte>((int)length, pinned: true);
                    using (var fs = new FileStream(file, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20))
                        fs.ReadExactly(s_pinned);
                    byte* p = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(s_pinned, 0);
                    return EsaSuffixTree.FromFile(SesaFile.FromMemory(p, length, options));
                }
            case "resident-native":
                {
                    long length = new FileInfo(file).Length;
                    s_native = NativeMemory.AlignedAlloc((nuint)length, 4096);
                    using (var fs = new FileStream(file, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20))
                    {
                        var span = new Span<byte>(s_native, (int)length);
                        fs.ReadExactly(span);
                    }
                    return EsaSuffixTree.FromFile(SesaFile.FromMemory((byte*)s_native, length, options));
                }
            case "memory-packed":
                return EsaSuffixTree.Build(text, ChildLayoutKind.Packed);
            case "memory-tagged":
                return EsaSuffixTree.Build(text, ChildLayoutKind.Tagged);
            default:
                throw new ArgumentException($"unknown variant '{variant}'");
        }
    }
}
#endif
