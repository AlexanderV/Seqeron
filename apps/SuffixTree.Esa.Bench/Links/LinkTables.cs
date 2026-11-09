namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// The two candidate sources of the link endpoints (prototype P1): the inverse suffix array <c>ISA</c> and the successor array <c>Psi[i] = ISA[SA[i] + 1]</c>.
/// Both have n + 1 entries (4 bytes each). Built from the suffix array alone in one scattered pass (ISA) and one gather pass (Psi).
/// </summary>
internal static class LinkTables
{
    /// <summary><c>isa[sa[i]] = i</c> for every rank <c>i</c> (so <c>isa[n] = 0</c>, the rank of the empty suffix).</summary>
    internal static int[] BuildIsa(int[] sa)
    {
        ArgumentNullException.ThrowIfNull(sa);
        var isa = new int[sa.Length];
        for (int r = 0; r < sa.Length; r++) isa[sa[r]] = r;
        return isa;
    }

    /// <summary><c>psi[0] = -1</c> (the empty suffix has no successor), <c>psi[i] = isa[sa[i] + 1]</c> for <c>i &gt;= 1</c>.</summary>
    internal static int[] BuildPsi(int[] sa, int[] isa)
    {
        ArgumentNullException.ThrowIfNull(sa);
        ArgumentNullException.ThrowIfNull(isa);
        var psi = new int[sa.Length];
        if (sa.Length == 0) return psi;
        psi[0] = -1;
        for (int r = 1; r < sa.Length; r++) psi[r] = isa[sa[r] + 1];
        return psi;
    }
}
