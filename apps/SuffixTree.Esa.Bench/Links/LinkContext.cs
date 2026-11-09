#pragma warning disable S2325 // stubs of the red run: the bodies are written after the tests
namespace SuffixTree.Esa.Bench.Links;

/// <summary>Everything the link kernels read: the index arrays and the optional navigation structures (ISA, Psi, the LCP pyramid).</summary>
internal sealed class LinkContext
{
    /// <exception cref="NotSupportedException">The child table of <paramref name="index"/> is not the tagged layout.</exception>
    internal LinkContext(EsaIndex index, int[]? isa, int[]? psi, LcpRangeHelper? helper) => throw new NotImplementedException();

    /// <summary>The text length n.</summary>
    internal int N => throw new NotImplementedException();

    /// <summary>The suffix array (n + 1 entries).</summary>
    internal int[] Sa => throw new NotImplementedException();

    /// <summary>The compact LCP array.</summary>
    internal CompactLcp Lcp => throw new NotImplementedException();

    /// <summary>The child table (tagged layout).</summary>
    internal TaggedChild Child => throw new NotImplementedException();

    /// <summary>The inverse suffix array, if built.</summary>
    internal int[]? Isa => throw new NotImplementedException();

    /// <summary>The successor array, if built.</summary>
    internal int[]? Psi => throw new NotImplementedException();

    /// <summary>The range-minimum / previous-smaller / next-smaller helper, if built.</summary>
    internal LcpRangeHelper? Helper => throw new NotImplementedException();
}
