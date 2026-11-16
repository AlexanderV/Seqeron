using SuffixTree.Esa;
using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.EmProto.Ph6;

/// <summary>One kernel x endpoint-source cell of the smoke test.</summary>
internal sealed record SmokeCell(LinkKernel Kernel, EndpointSource Source, long Nodes, long Mismatches, long GaveUp, long DefinitionMismatches);

internal sealed record SmokeResult(long NodeCount, long SampleCount, SmokeCell[] Cells)
{
    internal long Mismatches => Cells.Sum(c => c.Mismatches);

    internal long DefinitionMismatches => Cells.Sum(c => c.DefinitionMismatches);

    internal bool Ok => Mismatches == 0 && DefinitionMismatches == 0;
}

/// <summary>
/// The suffix-link smoke test of E9: every kernel of prototype P1 with every endpoint source on a sample of internal nodes, once with the in-RAM ISA / Psi of <c>LinkTables</c> and once with the tables of the streaming code
/// (ISA from DIST, Psi from the BWT, Psi from the ISA); the results (success flag, target, step count) must be identical, and, where a definition is given, equal to it.
/// </summary>
internal static class LinkSmoke
{
    internal static SmokeResult Run(EsaIndex index, int[] isaRam, int[] psiRam, int[] isaOther, int[] psiOther, int sample, ulong seed, long budget, int fanout, bool p1Sampling, Func<LinkNode, LinkTarget>? definition) => throw new NotImplementedException();
}
