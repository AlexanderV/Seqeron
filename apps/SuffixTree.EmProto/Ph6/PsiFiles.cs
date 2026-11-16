namespace SuffixTree.EmProto.Ph6;

/// <summary>The ISA and PSI files (n + 1 little-endian 4-byte entries) and their relation to the in-RAM tables of the suffix-link prototype P1.</summary>
internal static class PsiFiles
{
    /// <summary>The file as int32 entries.</summary>
    internal static int[] ReadInts(string path) => throw new NotImplementedException();

    /// <summary>A PSI file as P1 reads it: entry 0 (the file convention 0 or the cyclic value) becomes -1, the P1 marker of the empty suffix.</summary>
    internal static int[] LoadForP1(string psiPath) => throw new NotImplementedException();

    /// <summary>The in-RAM ISA and Psi of P1 (<c>LinkTables</c>) from a suffix array file, written as ISA / PSI files (psi[0] = -1 is written as 0).</summary>
    internal static void WriteOracle(string saPath, string isaPath, string psiPath) => throw new NotImplementedException();
}
