#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
namespace SuffixTree.EmProto.Ph5;

/// <summary>What the backward scan did.</summary>
internal sealed class NsvResult
{
    internal long StackPeak { get; set; }

    internal SpillStats Spill { get; set; }

    internal WindowStats Windows { get; set; }

    internal double Seconds { get; set; }
}

/// <summary>STUB (tests first): the backward scan of layout F (next smaller value, written in descending rank order).</summary>
internal static class NsvScan
{
    internal static NsvResult Run(long n, CompactLcpFile lcp, string nsvPath, ChildBuildOptions options) => throw new NotImplementedException();
}
