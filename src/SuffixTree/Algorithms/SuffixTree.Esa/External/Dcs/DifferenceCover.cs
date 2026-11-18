#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186, S3453, CS9113 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>How a <see cref="DifferenceCover"/> was constructed.</summary>
internal enum CoverKind
{
    Singer,
    Pruned,
}

/// <summary>STUB (work package 1b.7d, tests first).</summary>
internal sealed class DifferenceCover
{
    private DifferenceCover()
    {
    }

    internal int Period => throw new NotImplementedException();

    internal int[] Residues => throw new NotImplementedException();

    internal CoverKind Kind => throw new NotImplementedException();

    internal int[] IndexOf => throw new NotImplementedException();

    internal int[] Witness => throw new NotImplementedException();

    internal double Density => throw new NotImplementedException();

    internal int PeriodBits => throw new NotImplementedException();

    internal bool Verify() => throw new NotImplementedException();

    internal bool IsPerfect() => throw new NotImplementedException();

    internal static DifferenceCover Create(int period) => throw new NotImplementedException();

    internal static int NearestSingerPeriod(int target) => throw new NotImplementedException();

    internal long SampleCount(long n) => throw new NotImplementedException();
}
