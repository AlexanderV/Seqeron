#pragma warning disable S2325, S2933, S1144, CA1823, IDE0052, CA1822, S1172, IDE0060, CA1801, S3251, S1186 // STUB: removed when the implementation lands
namespace SuffixTree.Esa.External;

/// <summary>STUB (work package 1b.7d, tests first). Counters of the work of one task.</summary>
internal sealed class SortWork
{
    internal long Compares;
    internal long TiedCompares;
    internal long LceCalls;
    internal long LceSymbols;
    internal long LceSteps;
    internal long RankLookups;
    internal long TiePositions;
    internal long LeafPositions;
    internal long ScanPositions;

    internal void Add(SortWork o) => throw new NotImplementedException();
}

/// <summary>STUB (work package 1b.7d, tests first). The suffix comparator of the bucket sort.</summary>
internal sealed class SuffixOrder(MappedTextReader text, RankTable ranks, SortWork work)
{
    internal MappedTextReader Text { get; } = text;

    internal RankTable Ranks { get; } = ranks;

    internal SortWork Work { get; } = work;

    internal int Compare(uint i, uint j) => throw new NotImplementedException();

    internal int CompareTied(uint i, uint j) => throw new NotImplementedException();
}
