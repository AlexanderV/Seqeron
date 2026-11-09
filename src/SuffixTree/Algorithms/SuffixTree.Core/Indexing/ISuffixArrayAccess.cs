namespace SuffixTree.Indexing;

/// <summary>
/// Direct access to the arrays of an engine that really has a suffix array, its inverse and the LCP array. The array has <c>n + 1</c> entries:
/// rank 0 is the empty suffix (<c>SA[0] = n</c>), ranks <c>1..n</c> are the suffixes in ascending order, <c>LCP[r] = lcp(SA[r-1], SA[r])</c>, <c>LCP[0] = LCP[1] = 0</c>.
/// </summary>
public interface ISuffixArrayAccess : ISuffixIndex
{
    /// <summary><c>SA[rank]</c>: the start position of the suffix of the given rank.</summary>
    long SuffixAt(long rank);

    /// <summary><c>ISA[position]</c>: the rank of the suffix starting at <paramref name="position"/> (<c>RankOf(n) == 0</c>).</summary>
    long RankOf(long position);

    /// <summary><c>LCP[rank]</c>.</summary>
    long LcpAt(long rank);

    /// <summary>The longest common prefix of the suffixes of two ranks, a range minimum over the LCP array; <c>Lcp(r, r)</c> is the length of the suffix.</summary>
    long Lcp(long rankA, long rankB);

    /// <summary>The lowest common ancestor (an explicit node) of the leaves of two ranks; for equal ranks the leaf itself.</summary>
    IndexLocus Lca(long rankA, long rankB);
}
