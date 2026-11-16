using System.Runtime.InteropServices;
using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.EmProto.Ph6;

/// <summary>The ISA and PSI files (n + 1 little-endian 4-byte entries) and their relation to the in-RAM tables of the suffix-link prototype P1.</summary>
internal static class PsiFiles
{
    /// <summary>The file as int32 entries.</summary>
    internal static int[] ReadInts(string path)
    {
        if (!BitConverter.IsLittleEndian) throw new PlatformNotSupportedException("the files are little endian, as the other files of the prototype");
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20, FileOptions.SequentialScan);
        if (fs.Length % 4 != 0) throw new InvalidDataException($"{path}: {fs.Length} bytes is not a whole number of 4-byte entries");
        var values = new int[fs.Length / 4];
        fs.ReadExactly(MemoryMarshal.AsBytes(values.AsSpan()));
        return values;
    }

    /// <summary>A PSI file as P1 reads it: entry 0 (the file convention 0 or the cyclic value) becomes -1, the P1 marker of the empty suffix.</summary>
    internal static int[] LoadForP1(string psiPath)
    {
        int[] psi = ReadInts(psiPath);
        if (psi.Length > 0) psi[0] = -1;
        return psi;
    }

    /// <summary>The in-RAM ISA and Psi of P1 (<c>LinkTables</c>) from a suffix array file, written as ISA / PSI files (psi[0] = -1 is written as 0).</summary>
    internal static void WriteOracle(string saPath, string isaPath, string psiPath)
    {
        int[] sa = ReadInts(saPath);
        int[] isa = LinkTables.BuildIsa(sa);
        int[] psi = LinkTables.BuildPsi(sa, isa);
        if (psi.Length > 0) psi[0] = 0;
        Write(isaPath, isa);
        Write(psiPath, psi);
    }

    private static void Write(string path, int[] values)
    {
        using var fs = new FileStream(path, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20);
        fs.Write(MemoryMarshal.AsBytes(values.AsSpan()));
    }
}
