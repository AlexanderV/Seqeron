#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
namespace SuffixTree.Esa.Bench.Links;

/// <summary>Counters of one <see cref="MemDriver.Run"/>.</summary>
internal sealed record MemStats(long Starts, long Links, long ExtendSteps, long ChildLookups, long RescanLookups, long Emitted);

/// <summary>
/// The throwaway matching-statistics driver of prototype P2 (work package 1b.4a, not the navigator of 1b.4): MEM enumeration on the ESA with the structure of
/// <c>SuffixTreeAlgorithms.FindMaximalExactMatches</c> (extend, follow the suffix link of the last node, rescan) where the locus is a rank interval and the suffix link
/// comes from a <see cref="LinkKernel"/>.
/// </summary>
internal sealed class MemDriver
{
    /// <param name="index">The ESA (tagged child table, at most 256 distinct symbols).</param>
    /// <param name="links">ISA / Psi / helper for the kernel.</param>
    /// <param name="kernel">The link kernel.</param>
    /// <param name="source">The endpoint source.</param>
    /// <param name="loci">1: one chain (matching statistics of the ESA); 2: the generic algorithm's two loci, the second chain only costs its suffix links (cost emulation).</param>
    internal MemDriver(EsaIndex index, LinkContext links, LinkKernel kernel, EndpointSource source, int loci) => throw new NotImplementedException();

    /// <summary>Counters of the last <see cref="Run"/>.</summary>
    internal MemStats Stats => throw new NotImplementedException();

    /// <summary>All maximal exact matches of at least <paramref name="minLength"/> symbols, ascending by query position, then text position.</summary>
    internal IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> Run(string query, int minLength) => throw new NotImplementedException();
}
