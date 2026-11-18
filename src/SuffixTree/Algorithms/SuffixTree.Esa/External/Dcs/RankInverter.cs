namespace SuffixTree.Esa.External;

/// <summary>
/// The mapper of the inverse permutation of the suffix array of a reduced string, for DIST (ADR-004 section 3.3, the rank file as a target): the entry <c>i</c> of the suffix array (the reduced position <c>SA[i]</c>) becomes the
/// tuple <c>(slot of SA[i], i)</c>; the empty suffix (entry 0) and the separators of the classes are skipped. The rank of a sample suffix is its index in the suffix array of the reduced string, separators included (they occupy the
/// first |D| ranks), the convention of the prototype's <c>DcsBuilder</c>. Counts what it emitted, so that the caller can see that every sample position got a rank.
/// </summary>
internal struct SampleRankMapper(SampleGeometry geometry) : IDistMapper<uint>
{
    /// <summary>Tuples emitted so far.</summary>
    internal long Emitted { get; private set; }

    public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
    {
        key = 0;
        value = 0;
        if (index == 0) return false; // the empty suffix
        if (sourceValue >= geometry.RedLength) throw new InvalidDataException($"SA[{index}] = {sourceValue} lies beyond the reduced string of {geometry.RedLength} symbols");
        if (!geometry.TrySlotOfRed(sourceValue, out long slot)) return false; // a separator
        key = (uint)slot;
        value = (uint)index;
        Emitted++;
        return true;
    }
}

/// <summary>The inverse permutation of the suffix array of the reduced string, straight into the rank file of the level above, by DIST (pass 1 distributes, pass 2 scatters chunk by chunk and appends through the mapped window writer).</summary>
internal static class RankInverter
{
    /// <summary>Writes the rank file <paramref name="rankPath"/> (<c>g.Slots</c> entries; the unused slots stay 0) from the suffix array file <paramref name="saPath"/> of the reduced string (<c>RedLength + 1</c> entries).</summary>
    /// <exception cref="InvalidDataException">The suffix array holds an entry beyond the reduced string, or does not rank every sample position exactly once.</exception>
    internal static DistStats Invert(string saPath, SampleGeometry g, string rankPath, DistOptions options)
    {
        var mapper = new SampleRankMapper(g);
        DistStats stats = Dist.Scatter(saPath, g.RedLength + 1, ref mapper, rankPath, new DistTarget<uint>(g.Slots, 0, false), options);
        if (mapper.Emitted != g.Samples) throw new InvalidDataException($"the suffix array of the reduced string ranks {mapper.Emitted} positions of {g.Samples}");
        return stats;
    }
}
