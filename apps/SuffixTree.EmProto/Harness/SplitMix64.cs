namespace SuffixTree.EmProto.Harness;

/// <summary>
/// SplitMix64 (Steele, Lea, Flood 2014; Vigna's public-domain constants): the same generator and the same <c>Next() % bound</c> rule as the test scaffolding of
/// the ESA project, so that families generated here equal those of <c>StringFamilies</c> and of the committed libsais vectors.
/// </summary>
internal sealed class SplitMix64(ulong seed)
{
    private ulong _state = seed;

    internal ulong Next()
    {
        unchecked
        {
            _state += 0x9E3779B97F4A7C15UL;
            ulong z = _state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
            return z ^ (z >> 31);
        }
    }

    internal int NextInt(int bound)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bound, 0);
        return (int)(Next() % (ulong)bound);
    }

    /// <summary>A double in [0, 1) from the top 53 bits.</summary>
    internal double NextDouble() => (Next() >> 11) * (1.0 / (1UL << 53));
}
