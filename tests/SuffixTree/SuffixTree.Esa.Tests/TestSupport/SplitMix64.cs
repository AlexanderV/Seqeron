namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// SplitMix64 pseudo-random generator (Steele, Lea, Flood 2014; Vigna's public-domain reference constants).
/// Used instead of <see cref="Random"/> because the sequence is a pure function of the seed and is identical in
/// every runtime version and in the out-of-repository generators (PROGRESS-S1.md, "Test strategy").
/// <see cref="NextInt"/> is <c>Next() % bound</c> (the modulo bias is irrelevant for test data and keeps the
/// rule reproducible in any language).
/// </summary>
public sealed class SplitMix64
{
    private ulong _state;

    /// <summary>Creates a generator whose whole sequence is determined by <paramref name="seed"/>.</summary>
    public SplitMix64(ulong seed) => _state = seed;

    /// <summary>Next 64-bit value.</summary>
    public ulong Next()
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

    /// <summary>A value in <c>[0, bound)</c>: <c>Next() % bound</c>.</summary>
    public int NextInt(int bound)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bound, 0);
        return (int)(Next() % (ulong)bound);
    }
}
