#pragma warning disable S2325, S4144, S1144, S3453, CA1822 // stub of the red run (work package 1b.4): every member throws NotImplementedException; the implementation replaces this file
namespace SuffixTree.Esa;

/// <summary>
/// The navigation pack of an <see cref="EsaIndex"/> (ADR-001 R-SL, work package 1b.4): the inverse suffix array (4 bytes per rank) and the <see cref="LcpPyramid"/> (0.13 bytes per
/// rank at fanout 32), which together give the suffix link of an internal node as <c>[PrevLess(Psi(l), d - 1) .. NextLess(Psi(r) + 1, d - 1) - 1]</c> with
/// <c>Psi(i) = ISA[SA[i] + 1]</c> (kernel (c2) of the Log of 1b.4a). Built lazily by <see cref="EsaIndex.Navigation"/> (once, thread-safe, heap-resident); it also holds the pinned views
/// of the suffix array and of the LCP base array that the navigator reads through pointers. Stub of the red run.
/// </summary>
internal sealed class EsaNavigationPack : IDisposable
{
    /// <summary>Builds the pack of <paramref name="index"/> (O(n): one scatter pass for ISA, one pass for the pyramid).</summary>
    internal static EsaNavigationPack Build(EsaIndex index, int fanout = LcpPyramid.DefaultFanout)
    {
        _ = (index, fanout);
        throw new NotImplementedException();
    }

    /// <summary>The index the pack belongs to.</summary>
    internal EsaIndex Index => throw new NotImplementedException();

    /// <summary>The root node <c>(0, n, 0, 0, n, true)</c>.</summary>
    internal EsaNode Root => throw new NotImplementedException();

    /// <summary>The LCP helper.</summary>
    internal LcpPyramid Pyramid => throw new NotImplementedException();

    /// <summary><c>ISA[position]</c>, the rank of the suffix that starts at <paramref name="position"/> (<c>0..n</c>; <c>ISA[n] = 0</c>).</summary>
    internal int Isa(int position)
    {
        _ = position;
        throw new NotImplementedException();
    }

    /// <summary>The inverse suffix array as a managed array (tests and diagnostics).</summary>
    internal int[] IsaArray => throw new NotImplementedException();

    /// <summary>Bytes of ISA plus the pyramid levels (what the pack adds to the index).</summary>
    internal long PayloadBytes => throw new NotImplementedException();

    /// <summary>
    /// The suffix link of <paramref name="node"/>: for an internal node <c>aY</c> the node <c>Y</c> (the root for a node of depth 1), for a leaf the leaf of the next suffix (the root for the
    /// terminator-only leaf), for the root the root. The result carries its parent depth, last-child flag and <c>Pos</c>.
    /// </summary>
    internal EsaNode SuffixLink(in EsaNode node)
    {
        _ = node;
        throw new NotImplementedException();
    }

    /// <summary>Releases the pins; idempotent.</summary>
    public void Dispose() => throw new NotImplementedException();
}
