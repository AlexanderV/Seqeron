namespace Seqeron.Genomics.Alignment;

/// <summary>The C library's <c>rand()</c> of glibc (TYPE_3 additive feedback generator, r[i] = r[i-3] + r[i-31], default seed 1): ProbCons never calls <c>srand</c>.</summary>
public sealed class GlibcRandom
{
    private readonly int[] _state = new int[31];
    private int _f = 3, _b;

    public GlibcRandom(uint seed = 1)
    {
        // srandom_r: seed 0 means 1; the state words come from the Lehmer generator 16807 mod 2^31 - 1, then 310 outputs are discarded
        int word = seed == 0 ? 1 : unchecked((int)seed);
        _state[0] = word;
        for (int i = 1; i < 31; i++)
        {
            int hi = word / 127773, lo = word % 127773;
            word = unchecked(16807 * lo - 2836 * hi);
            if (word < 0) word += 2147483647;
            _state[i] = word;
        }

        for (int i = 0; i < 310; i++) Next();
    }

    /// <summary>The next value, 0 ... 2^31 - 1.</summary>
    public int Next()
    {
        uint val = unchecked((uint)_state[_f] + (uint)_state[_b]);
        _state[_f] = unchecked((int)val);
        int result = (int)(val >> 1);
        if (++_f >= 31) _f = 0;
        if (++_b >= 31) _b = 0;
        return result;
    }
}
