#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186, S3453, CS9113 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal readonly struct FastMod32
{
    internal FastMod32(uint divisor) => throw new NotImplementedException();

    internal uint Divisor => throw new NotImplementedException();

    internal uint Div(uint a) => throw new NotImplementedException();

    internal uint Mod(uint a) => throw new NotImplementedException();
}

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal sealed unsafe class RankTable : IDisposable
{
    private RankTable()
    {
    }

    internal DifferenceCover Cover => throw new NotImplementedException();

    internal FastMod32 Mod => throw new NotImplementedException();

    internal long N => throw new NotImplementedException();

    internal long Slots => throw new NotImplementedException();

    internal static RankTable Open(string path, DifferenceCover cover, long n) => throw new NotImplementedException();

    internal uint Rank(long p) => throw new NotImplementedException();

    /// <summary>Reads every real slot once and checks that the ranks are distinct and that the unused slots of the last block are zero.</summary>
    internal void Verify() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
