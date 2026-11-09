namespace SuffixTree.Esa.Bench.Links;

/// <summary>Everything the link kernels read: the index arrays and the optional navigation structures (ISA, Psi, the LCP pyramid).</summary>
internal sealed class LinkContext
{
    /// <exception cref="NotSupportedException">The child table of <paramref name="index"/> is not the tagged layout.</exception>
    internal LinkContext(EsaIndex index, int[]? isa, int[]? psi, LcpRangeHelper? helper)
    {
        ArgumentNullException.ThrowIfNull(index);
        if (index.Child.Kind != ChildLayoutKind.Tagged)
            throw new NotSupportedException($"the link prototypes read the tagged child table, the index has {index.Child.Kind}");
        if (isa is not null && isa.Length != index.Length + 1) throw new ArgumentException("ISA needs n + 1 entries", nameof(isa));
        if (psi is not null && psi.Length != index.Length + 1) throw new ArgumentException("Psi needs n + 1 entries", nameof(psi));
        if (helper is not null && helper.Length != index.Length + 1) throw new ArgumentException("the helper must cover the n + 1 LCP entries", nameof(helper));
        N = index.Length;
        Sa = index.SuffixArray;
        Lcp = index.Lcp;
        Child = index.Child.Tagged;
        Isa = isa;
        Psi = psi;
        Helper = helper;
    }

    /// <summary>The text length n.</summary>
    internal int N { get; }

    /// <summary>The suffix array (n + 1 entries).</summary>
    internal int[] Sa { get; }

    /// <summary>The compact LCP array.</summary>
    internal CompactLcp Lcp { get; }

    /// <summary>The child table (tagged layout).</summary>
    internal TaggedChild Child { get; }

    /// <summary>The inverse suffix array, if built.</summary>
    internal int[]? Isa { get; }

    /// <summary>The successor array, if built.</summary>
    internal int[]? Psi { get; }

    /// <summary>The range-minimum / previous-smaller / next-smaller helper, if built.</summary>
    internal LcpRangeHelper? Helper { get; }
}
