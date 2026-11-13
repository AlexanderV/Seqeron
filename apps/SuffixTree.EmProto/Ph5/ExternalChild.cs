#pragma warning disable S2325, S1144, S1172, S3459, CS9113 // stub of the tests-first step: every member throws
namespace SuffixTree.EmProto.Ph5;

/// <summary>The files of a PH5 run (null = off): CHILD of layout A, PSV and NSV of layout F.</summary>
internal sealed record ChildBuildOutputs
{
    internal string? Child { get; init; }

    internal string? Psv { get; init; }

    internal string? Nsv { get; init; }
}

/// <summary>What a PH5 run found and did (the forward result, plus the backward pass when NSV was asked for).</summary>
internal sealed class ChildBuildStats
{
    internal required ChildScanResult Forward { get; init; }

    internal NsvResult? Backward { get; init; }
}

/// <summary>STUB (tests first): PH5 of ADR-004 on a given compact LCP file.</summary>
internal static class ExternalChild
{
    internal static ChildBuildStats BuildFromFiles(long n, string lcpBase, string lcpOverflow, ChildBuildOutputs outputs, ChildBuildOptions options) => throw new NotImplementedException();
}
