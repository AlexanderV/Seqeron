using System.Numerics;

namespace SuffixTree.EmProto.Dcs;

/// <summary>How a <see cref="DifferenceCover"/> was constructed.</summary>
internal enum CoverKind
{
    /// <summary>Singer perfect difference set: modulus <c>q^2 + q + 1</c> (q prime), <c>q + 1</c> residues, every non-zero difference exactly once.</summary>
    Singer,

    /// <summary>Square-root construction <c>{0..m-1} + {m, 2m, .., m*m}</c> pruned greedily (any modulus, larger than Singer).</summary>
    Pruned,
}

/// <summary>
/// A difference cover D of Z_X: for every residue d there are a, b in D with <c>b - a = d (mod X)</c> (Karkkainen, Fast BWT in small space by blockwise suffix sorting, TCS 2007, section 3).
/// The ADR-004 table of covers was unverified; this class constructs a cover for any X and <see cref="Verify"/> proves the property by exhaustion
/// (every residue difference, <c>O(|D|^2 + X)</c>); the tests call it for every X in [2, 4096].
/// <see cref="Witness"/> maps a difference to one a with <c>a in D</c> and <c>a + d in D</c> (the comparator of the suffix order needs exactly this table).
/// </summary>
internal sealed class DifferenceCover
{
    private DifferenceCover(int period, int[] residues, CoverKind kind)
    {
        Period = period;
        Residues = residues;
        Kind = kind;
        IndexOf = new int[period];
        Array.Fill(IndexOf, -1);
        for (int i = 0; i < residues.Length; i++) IndexOf[residues[i]] = i;
        Witness = new int[period];
        Array.Fill(Witness, -1);
        foreach (int a in residues)
        {
            foreach (int b in residues)
            {
                int d = (b - a + period) % period;
                if (Witness[d] < 0) Witness[d] = a;
            }
        }
    }

    /// <summary>The modulus X (the period of the sample).</summary>
    internal int Period { get; }

    /// <summary>The sorted residues of D.</summary>
    internal int[] Residues { get; }

    internal CoverKind Kind { get; }

    /// <summary><c>IndexOf[r]</c> = position of r in <see cref="Residues"/>, or -1.</summary>
    internal int[] IndexOf { get; }

    /// <summary><c>Witness[d]</c> = a residue a in D with <c>(a + d) mod X</c> in D, or -1 if d is not covered (never, for a valid cover).</summary>
    internal int[] Witness { get; }

    /// <summary>Fraction of the positions that belong to the sample: <c>|D| / X</c>.</summary>
    internal double Density => (double)Residues.Length / Period;

    /// <summary>True when every residue difference modulo X occurs among D (exhaustive check, independent of the <see cref="Witness"/> table).</summary>
    internal bool Verify()
    {
        var covered = new bool[Period];
        foreach (int a in Residues)
        {
            if (a < 0 || a >= Period) return false;
            foreach (int b in Residues) covered[(b - a + Period) % Period] = true;
        }
        return covered.All(c => c);
    }

    /// <summary>True when the cover is a perfect difference set: every non-zero residue is a difference exactly once.</summary>
    internal bool IsPerfect()
    {
        var count = new int[Period];
        foreach (int a in Residues)
        {
            foreach (int b in Residues)
            {
                if (a != b) count[(b - a + Period) % Period]++;
            }
        }
        for (int d = 1; d < Period; d++)
        {
            if (count[d] != 1) return false;
        }
        return true;
    }

    /// <summary>
    /// A cover of Z_<paramref name="period"/>: the Singer set when the period is <c>q^2 + q + 1</c> with q prime, otherwise the pruned square-root construction.
    /// Deterministic.
    /// </summary>
    internal static DifferenceCover Create(int period)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(period, 2);
        if (SingerPrime(period) is int q) return new DifferenceCover(period, SingerSet(q), CoverKind.Singer);
        return new DifferenceCover(period, Pruned(period), CoverKind.Pruned);
    }

    /// <summary>The Singer modulus <c>q^2 + q + 1</c> (q prime) closest to <paramref name="target"/> (ratio-wise); for the nominal X = 512 / 1024 / 2048 of ADR-004 this gives 553 / 993 / 1893.</summary>
    internal static int NearestSingerPeriod(int target)
    {
        int best = 7;
        double bestRatio = double.MaxValue;
        for (int q = 2; q < 200; q++)
        {
            if (!IsPrime(q)) continue;
            int x = (q * q) + q + 1;
            double ratio = Math.Max((double)x / target, (double)target / x);
            if (ratio < bestRatio)
            {
                bestRatio = ratio;
                best = x;
            }
        }
        return best;
    }

    private static bool IsPrime(int v)
    {
        if (v < 2) return false;
        for (int d = 2; d * d <= v; d++)
        {
            if (v % d == 0) return false;
        }
        return true;
    }

    private static int? SingerPrime(int period)
    {
        for (int q = 2; (q * q) + q + 1 <= period; q++)
        {
            if (IsPrime(q) && (q * q) + q + 1 == period) return q;
        }
        return null;
    }

    /// <summary>
    /// Singer's construction (1938) for prime q: find a primitive cubic over GF(q), so that alpha = x generates GF(q^3)*; the residues i in [0, q^2+q+1) for which
    /// the constant coefficient of alpha^i is zero form a perfect difference set (the zero set of a non-zero linear functional is a plane through the origin of GF(q)^3, and alpha^N
    /// with N = q^2+q+1 generates GF(q)*, so membership depends on i mod N only).
    /// </summary>
    private static int[] SingerSet(int q)
    {
        int n = (q * q) + q + 1;
        long order = ((long)q * q * q) - 1;
        for (int c0 = 1; c0 < q; c0++)
        {
            for (int c1 = 0; c1 < q; c1++)
            {
                for (int c2 = 0; c2 < q; c2++)
                {
                    // x^3 = -c2 x^2 - c1 x - c0
                    int a0 = 1, a1 = 0, a2 = 0;
                    var zero = new List<int>();
                    bool primitive = true;
                    for (long k = 0; k < order; k++)
                    {
                        if (k < n && a0 == 0) zero.Add((int)k);
                        // multiply by x
                        int t = a2;
                        a2 = a1;
                        a1 = a0;
                        a0 = 0;
                        if (t != 0)
                        {
                            a0 = (((q - (t * c0 % q)) % q) + a0) % q;
                            a1 = (((q - (t * c1 % q)) % q) + a1) % q;
                            a2 = (((q - (t * c2 % q)) % q) + a2) % q;
                        }
                        if (a0 == 1 && a1 == 0 && a2 == 0 && k + 1 < order)
                        {
                            primitive = false;
                            break;
                        }
                    }
                    if (primitive && a0 == 1 && a1 == 0 && a2 == 0 && zero.Count == q + 1) return [.. zero];
                }
            }
        }
        throw new InvalidOperationException($"no primitive cubic found for q = {q}");
    }

    private static int[] Pruned(int x)
    {
        int m = 1;
        while (m * m < x) m++;
        var start = new SortedSet<int>();
        for (int j = 0; j < m; j++) start.Add(j % x);
        for (int k = 1; k <= m; k++) start.Add(m * k % x);
        List<int> best = [.. start];
        var rng = new Harness.SplitMix64(0x9E3779B97F4A7C15UL ^ (ulong)x);
        int rounds = 6;
        if (x <= 256) rounds = 40;
        else if (x <= 1024) rounds = 12;
        for (int round = 0; round < rounds; round++)
        {
            List<int> d = [.. start];
            // random order of removal attempts; keep a removal only if the cover property survives
            var order = Enumerable.Range(0, d.Count).Select(i => d[i]).OrderBy(_ => rng.Next()).ToList();
            var count = new int[x]; // number of ordered pairs (a, b) with b - a = delta
            foreach (int a in d)
            {
                foreach (int b in d) count[(b - a + x) % x]++;
            }
            foreach (int r in order)
            {
                // removing r deletes the pairs (r, b), (a, r): (r, r) counted once
                var delta = new List<int>(2 * d.Count);
                foreach (int b in d)
                {
                    delta.Add((b - r + x) % x);
                    if (b != r) delta.Add((r - b + x) % x);
                }
                bool ok = true;
                foreach (int dd in delta) count[dd]--;
                foreach (int dd in delta)
                {
                    if (count[dd] <= 0)
                    {
                        ok = false;
                        break;
                    }
                }
                if (ok) d.Remove(r);
                else foreach (int dd in delta) count[dd]++;
            }
            if (d.Count < best.Count) best = d;
        }
        best.Sort();
        return [.. best];
    }

    /// <summary>Number of positions of the sample among <c>0..n</c> (position n is part of the text: the empty suffix).</summary>
    internal long SampleCount(long n)
    {
        long blocks = n / Period;
        long rest = n % Period;
        long count = blocks * Residues.Length;
        foreach (int r in Residues)
        {
            if (r <= rest) count++;
        }
        return count;
    }

    /// <summary>log2 of the modulus rounded up (used to size buffers).</summary>
    internal int PeriodBits => BitOperations.Log2((uint)Period - 1) + 1;
}
