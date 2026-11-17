using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Experiments;

/// <summary>E8 of ADR-004 section 5: the hand-off package. Commands <c>e8-gen</c> (the streamed genome-like text) and <c>e8-verify</c> (the verifier of the outputs of <c>e6</c>).</summary>
internal static class E8
{
    internal static int Gen(Options o, Report r) => throw new NotImplementedException();

    internal static int Verify(Options o, Report r) => throw new NotImplementedException();
}

/// <summary>Parsing of counts with decimal suffixes (<c>k</c>, <c>m</c>, <c>g</c>), the unit of the E8 protocol (1 Gbp = 10^9 symbols).</summary>
internal static class E8Parsing
{
    internal static long Count(string text) => throw new NotImplementedException();
}
