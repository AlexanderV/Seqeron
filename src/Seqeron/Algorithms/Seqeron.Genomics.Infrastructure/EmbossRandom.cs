namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The random number generator of the EMBOSS library (<c>ajRandomSeed</c>, <c>ajRandomDouble</c>, <c>ajRandomNumber</c> of <c>ajmath.c</c>, 6.6.0): the Wichmann-Hill generator fills a table of 101 terms and a
/// lagged generator modulo 1000009711 over the table is added to a second generator modulo 1 (based on dprand / sdprand by N. M. Maclaren, Cambridge 1992). One instance is the state of one program run.
/// </summary>
internal sealed class EmbossRandom
{
    private const double XMod = 1000009711.0;
    private const double YMod = 33554432.0;
    private const double Tiny = 1.0e-17;

    private readonly double[] _poly = new double[101];
    private double _other;
    private int _index;
    private bool _seeded;

    /// <summary>The seed of the run (<c>ajRandomSeed</c>: <c>tv_usec % 9999 + 1</c>, or 0 with the variable <c>timetoday</c>).</summary>
    private readonly int _seed;

    public EmbossRandom(int seed) => _seed = seed;

    /// <summary><c>ajRandomSeed</c>: the first call seeds the generator, later calls do nothing.</summary>
    public void Seed()
    {
        if (_seeded)
            return;
        _seeded = true;
        int ix = (_seed >= 0 ? _seed : -_seed) % 10000 + 1;
        int iy = 2 * ix + 1;
        int iz = 3 * ix + 1;
        double x = 0.0;
        for (int i = -11; i < 101; ++i)
        {
            if (i >= 0)
                _poly[i] = Math.Floor(XMod * x);
            ix = (171 * ix) % 30269;
            iy = (172 * iy) % 30307;
            iz = (170 * iz) % 30323;
            x = ix / 30269.0 + iy / 30307.0 + iz / 30323.0;
            x -= Math.Floor(x);
        }

        _other = Math.Floor(YMod * x) / YMod;
        _index = 0;
    }

    /// <summary><c>ajRandomDouble</c>: a number in (0, 1).</summary>
    public double Double()
    {
        const double offset = 1.0 / YMod;
        const double xmod2 = 2.0 * XMod;
        const double xmod4 = 4.0 * XMod;
        if (!_seeded)
            Seed();
        int n = _index - 64;
        if (n < 0)
            n += 101;
        double x = _poly[_index] + _poly[_index];
        x = xmod4 - _poly[n] - _poly[n] - x - x - _poly[_index];
        if (x <= 0.0)
        {
            if (x < -XMod)
                x += xmod2;
            if (x < 0.0)
                x += XMod;
        }
        else
        {
            if (x >= xmod2)
            {
                x -= xmod2;
                if (x >= XMod)
                    x -= XMod;
            }

            if (x >= XMod)
                x -= XMod;
        }

        _poly[_index] = x;
        if (++_index >= 101)
            _index = 0;
        double y;
        do
        {
            y = 37.0 * _other + offset;
            _other = y - Math.Floor(y);
        }
        while (Math.Abs(_other) <= 2.22e-15);
        if ((x = x / XMod + _other) >= 1.0)
            x -= 1.0;
        return x + Tiny;
    }

    /// <summary><c>ajRandomNumber</c>: an integer 0 to 32767.</summary>
    public int Number() => (int)Math.Floor(Double() * 32768.0);

    /// <summary><c>ajStrRandom</c>: the characters of a string in random order (a random number for every position, then the shell sort <c>ajSortUintIncI</c> of the positions by it).</summary>
    public string Shuffle(string s)
    {
        int len = s.Length;
        var na = new uint[len];
        var rn = new uint[len];
        for (int i = 0; i < len; i++)
        {
            na[i] = (uint)i;
            rn[i] = (uint)Number();
        }

        SortUintIncI(rn, na, len);
        var q = new char[len];
        for (int i = 0; i < len; i++)
            q[i] = s[(int)na[i]];
        return new string(q);
    }

    /// <summary><c>ajSortIntIncI</c>: Shell sort of the index array <paramref name="p"/> by the keys a[p[.]] (an int key array).</summary>
    public static void SortIntIncI(int[] a, uint[] p, int n)
    {
        for (uint s = (uint)n / 2; s > 0; s /= 2)
        {
            for (uint i = s; i < n; ++i)
            {
                for (long j = (long)i - s; j >= 0 && a[p[j]] >= a[p[j + s]]; j -= s)
                    (p[j], p[j + s]) = (p[j + s], p[j]);
            }
        }
    }

    // ajSortUintIncI: Shell sort of the index array p by the keys a[p[.]] (swap while a[p[j]] >= a[p[j+s]])
    private static void SortUintIncI(uint[] a, uint[] p, int n)
    {
        for (uint s = (uint)n / 2; s > 0; s /= 2)
        {
            for (uint i = s; i < n; ++i)
            {
                for (long j = (long)i - s; j >= 0 && a[p[j]] >= a[p[j + s]]; j -= s)
                    (p[j], p[j + s]) = (p[j + s], p[j]);
            }
        }
    }
}
