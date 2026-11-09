#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
namespace SuffixTree.Esa.Bench.Links;

/// <summary>
/// The two candidate sources of the link endpoints (prototype P1): the inverse suffix array <c>ISA</c> and the successor array <c>Psi[i] = ISA[SA[i] + 1]</c>.
/// Both have n + 1 entries (4 bytes each).
/// </summary>
internal static class LinkTables
{
    /// <summary><c>isa[sa[i]] = i</c> for every rank <c>i</c> (so <c>isa[n] = 0</c>, the rank of the empty suffix).</summary>
    internal static int[] BuildIsa(int[] sa) => throw new NotImplementedException();

    /// <summary><c>psi[0] = -1</c> (the empty suffix has no successor), <c>psi[i] = isa[sa[i] + 1]</c> for <c>i &gt;= 1</c>.</summary>
    internal static int[] BuildPsi(int[] sa, int[] isa) => throw new NotImplementedException();
}
