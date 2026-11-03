namespace SuffixTree.Esa;

/// <summary>
/// Observer of the intermediate state of a suffix-array kernel (safety infrastructure, review section 6 T1). The kernels take it as a generic struct
/// parameter: with <see cref="NoProbe"/> (<c>Enabled</c> false) every call and the copies behind it are removed by the JIT, so the profile build and the
/// production build are one source (review MIN-8). The tests plug in a recorder and compare the states with the definitional oracle.
/// </summary>
internal interface ISaisProbe
{
    /// <summary>False for <see cref="NoProbe"/>: the kernel skips all probe work.</summary>
    static abstract bool Enabled { get; }

    /// <summary>Symbol counts of recursion <paramref name="level"/> (index = symbol).</summary>
    void Counts(int level, ReadOnlySpan<int> counts);

    /// <summary>The bucket array right after <c>SetBucketTails</c> of round 1 (exclusive bucket ends, n + 1 layout).</summary>
    void Tails(int level, ReadOnlySpan<int> tails);

    /// <summary>The bucket array right after <c>SetBucketHeads</c> of round 1 (bucket starts, n + 1 layout).</summary>
    void Heads(int level, ReadOnlySpan<int> heads);

    /// <summary>The LMS positions after induction round 1, i.e. sorted by their LMS substrings (equal substrings in an unspecified order).</summary>
    void LmsSorted(int level, ReadOnlySpan<int> positions);

    /// <summary>The number of names and the reduced string (the names in text order).</summary>
    void Reduced(int level, int names, ReadOnlySpan<int> reduced);

    /// <summary>The LMS positions in text order.</summary>
    void LmsPositions(int level, ReadOnlySpan<int> positions);
}

/// <summary>The probe of production builds: compiled out.</summary>
internal readonly struct NoProbe : ISaisProbe
{
    public static bool Enabled => false;

    public void Counts(int level, ReadOnlySpan<int> counts)
    {
    }

    public void Tails(int level, ReadOnlySpan<int> tails)
    {
    }

    public void Heads(int level, ReadOnlySpan<int> heads)
    {
    }

    public void LmsSorted(int level, ReadOnlySpan<int> positions)
    {
    }

    public void Reduced(int level, int names, ReadOnlySpan<int> reduced)
    {
    }

    public void LmsPositions(int level, ReadOnlySpan<int> positions)
    {
    }
}
