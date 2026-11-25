// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>
/// R's default RNG (RNG.c): Mersenne-Twister (<c>MT_genrand</c>, unsigned 32-bit state, output ·2.3283064365386963e−10,
/// <c>fixup</c> into (0, 1)), seeded as <c>set.seed(seed)</c> does (<c>Randomize</c>: 50 LCG scrambles
/// <c>seed = 69069·seed + 1</c>, then 625 LCG values fill <c>dummy[0..624]</c>, <c>dummy[0] = mti = 624</c>), with R ≥ 3.6
/// <c>sample.kind = "Rejection"</c> integer draws (<c>R_unif_index</c>: <c>rbits(ceil(log2(n)))</c> until &lt; n).
/// E.g. <c>set.seed(42); sample(5, 10, TRUE)</c> = 1 5 1 1 2 4 2 2 1 4.
/// </summary>
public sealed class RMersenneTwister
{
    private const int N = 624;
    private const int M = 397;
    private const uint MatrixA = 0x9908b0df;
    private const uint UpperMask = 0x80000000;
    private const uint LowerMask = 0x7fffffff;
    private const double TwoTo32Inverse = 2.3283064365386963e-10;
    private const double I2To32M1 = 2.328306437080797e-10;
    private readonly uint[] _mt = new uint[N];
    private int _mti;

    public RMersenneTwister(int seed)
    {
        uint s = unchecked((uint)seed);
        for (int j = 0; j < 50; j++)
        {
            s = unchecked((69069u * s) + 1u);
        }

        s = unchecked((69069u * s) + 1u); // dummy[0] (mti), overwritten by FixupSeeds
        for (int j = 0; j < N; j++)
        {
            s = unchecked((69069u * s) + 1u);
            _mt[j] = s;
        }

        _mti = N;
    }

    /// <summary>R <c>unif_rand()</c>.</summary>
    public double UnifRand()
    {
        double v = Genrand();
        if (v <= 0.0)
        {
            return 0.5 * I2To32M1;
        }

        return 1.0 - v <= 0.0 ? 1.0 - (0.5 * I2To32M1) : v;
    }

    /// <summary>R <c>R_unif_index(n)</c> (rejection sampling); 0-based index in [0, n).</summary>
    public int UnifIndex(int n)
    {
        if (n <= 0)
        {
            return 0;
        }

        int bits = (int)Math.Ceiling(Math.Log2(n));
        double dv;
        do
        {
            long v = 0;
            for (int k = 0; k <= bits; k += 16)
            {
                int v1 = (int)Math.Floor(UnifRand() * 65536);
                v = (65536 * v) + v1;
            }

            dv = v & ((1L << bits) - 1);
        }
        while (n <= dv);

        return (int)dv;
    }

    /// <summary>R <c>sample.int(n, k)</c> without replacement (random.c <c>do_sample</c>, equal probabilities,
    /// n ≤ 1e7): <c>x = 0..n−1</c>; k times <c>j = R_unif_index(n); y[i] = x[j]; x[j] = x[--n]</c>. 0-based indices in
    /// draw order.</summary>
    public int[] SampleWithoutReplacement(int n, int k)
    {
        var x = new int[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i;
        }

        var y = new int[k];
        for (int i = 0; i < k; i++)
        {
            int j = UnifIndex(n);
            y[i] = x[j];
            x[j] = x[--n];
        }

        return y;
    }

    private double Genrand()
    {
        uint y;
        if (_mti >= N)
        {
            int kk;
            for (kk = 0; kk < N - M; kk++)
            {
                y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                _mt[kk] = _mt[kk + M] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
            }

            for (; kk < N - 1; kk++)
            {
                y = (_mt[kk] & UpperMask) | (_mt[kk + 1] & LowerMask);
                _mt[kk] = _mt[kk + (M - N)] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
            }

            y = (_mt[N - 1] & UpperMask) | (_mt[0] & LowerMask);
            _mt[N - 1] = _mt[M - 1] ^ (y >> 1) ^ ((y & 1u) != 0 ? MatrixA : 0u);
            _mti = 0;
        }

        y = _mt[_mti++];
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680u;
        y ^= (y << 15) & 0xefc60000u;
        y ^= y >> 18;
        return y * TwoTo32Inverse;
    }
}
