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

    /// <summary>
    /// The whole suffix array of the level (n + 1 entries) right after <paramref name="phase"/>: Round1Place, Round1InduceL, Round1InduceS, CompactLms, Name, Round2Place,
    /// Round2InduceL, Round2InduceS. The tests replay each induction and naming pass on the previous snapshot with the safe copy of the pass (1o.1c). For
    /// <see cref="SaisPhase.Reduced"/> the span is the suffix array of the REDUCED string (m + 1 entries, n + 1 layout) as it stands right after the recursion or the direct placement.
    /// </summary>
    void SaState(int level, SaisPhase phase, ReadOnlySpan<int> sa);

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

    public void SaState(int level, SaisPhase phase, ReadOnlySpan<int> sa)
    {
    }
}
