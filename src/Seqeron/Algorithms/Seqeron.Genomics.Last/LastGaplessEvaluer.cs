// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: LastEvaluer.cc and the ALP library of NCBI in src/alp: njn_localmaxstat.cpp, njn_localmaxstatutil.cpp, njn_dynprogproblim.cpp, njn_dynprogprob.cpp, njn_root.hpp), Copyright 2008 - 2019 Martin C. Frith, J. Spouge (NCBI); GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Collections.Concurrent;
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Last;

/// <summary>
/// The Karlin-Altschul parameters of ungapped local alignments of the scoring scheme (<c>LastEvaluer::initGapless</c>, <c>Njn::LocalMaxStat</c>): lambda by bisection and K by the dynamic program over the
/// ascending random walk (<c>descendingLadderEpoch</c>). <c>lastal</c> uses them for the default minimum score of gapless alignments (<c>calcMinScoreGapless</c>). The E-value parameters of gapped alignments (the importance-sampling
/// part of the ALP library) are not ported.
/// </summary>
internal sealed class LastGaplessEvaluer
{
    private const double RelTol = 1.0e-6;

    private static readonly ConcurrentDictionary<string, (double Lambda, double K)?> Cache = new();

    public double Lambda { get; }

    public double K { get; }

    private LastGaplessEvaluer(double lambda, double k)
    {
        Lambda = lambda;
        K = k;
    }

    /// <summary>The parameters of the scoring scheme, or null when LAST cannot get them (the scores are not "logarithmic": the expected score is not negative, or no score is positive).</summary>
    public static LastGaplessEvaluer? TryCreate(int[][] scores, int size, double[] letterProbs1, double[] letterProbs2)
    {
        var key = new StringBuilder();
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                key.Append(scores[i][j].ToString(CultureInfo.InvariantCulture)).Append(',');
        for (int i = 0; i < size; i++)
            key.Append(letterProbs1[i].ToString("R", CultureInfo.InvariantCulture)).Append(';').Append(letterProbs2[i].ToString("R", CultureInfo.InvariantCulture)).Append(';');
        var r = Cache.GetOrAdd(key.ToString(), _ => Compute(scores, size, letterProbs1, letterProbs2));
        return r is { } v ? new LastGaplessEvaluer(v.Lambda, v.K) : null;
    }

    /// <summary><c>LastEvaluer::minScore(evalue, area)</c>.</summary>
    public double MinScore(double evalue, double area)
    {
        double s = (Math.Log(K * area) - Math.Log(evalue)) / Lambda;
        return Math.Max(0.0, s);
    }

    private static (double Lambda, double K)? Compute(int[][] scores, int size, double[] p1, double[] p2)
    {
        // initGapless(alphabetSize, matrix transposed, letterFreqs2, letterFreqs1): the frequencies are normalised, the scores flattened to a distribution
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < size; i++)
        {
            sum1 += p1[i];
            sum2 += p2[i];
        }

        var byScore = new SortedDictionary<long, double>();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                long score = scores[j][i];
                double pr = p2[i] / sum2 * (p1[j] / sum1);
                byScore[score] = byScore.GetValueOrDefault(score) + pr;
            }
        }

        var score_ = byScore.Where(kv => kv.Value > 0.0).Select(kv => kv.Key).ToArray();
        var prob = byScore.Where(kv => kv.Value > 0.0).Select(kv => kv.Value).ToArray();
        if (!IsLogarithmic(score_, prob))
            return null;
        double mu = 0;
        for (int i = 0; i < score_.Length; i++)
            mu += score_[i] * prob[i];
        double TotalProbAssoc(double x) => Enumerable.Range(0, score_.Length).Sum(i => prob[i] * Math.Exp(x * score_[i]));
        double MeanAssoc(double x) => Enumerable.Range(0, score_.Length).Sum(i => score_[i] * prob[i] * Math.Exp(x * score_[i]));
        double p = -Math.Log(prob[^1]) / score_[^1];
        while (1.0 <= TotalProbAssoc(p))
            p *= 0.5;
        double q = p / 0.5;
        double lambda = Bisection(1.0, TotalProbAssoc, p, q, RelTol * Math.Abs(p - q));
        double muAssoc = MeanAssoc(lambda);
        double thetaMin = Bisection(0.0, MeanAssoc, 0.0, lambda, RelTol * Math.Abs(p - q));
        double rMin = TotalProbAssoc(thetaMin);
        long delta = 0;
        foreach (long s in score_)
            delta = Gcd(delta, s);
        double thetaMinusDelta = (1.0 - Math.Exp(-lambda * delta)) / delta;

        // descendingLadderEpoch (weak ladder epochs)
        long iterMin = (long)(Math.Log(RelTol * (1.0 - rMin)) / Math.Log(rMin));
        long iter = Math.Max(0, iterMin);
        long yMax = (long)(-Math.Log(RelTol) / lambda);
        long morgue = score_[0] - 1;
        const long entry = 0;
        long valueBegin = morgue;
        int capacity = (int)(yMax - valueBegin);
        var arrays = new[] { new double[capacity], new double[capacity] };
        arrays[0][0 - valueBegin] = 1.0;
        long valueLower = 0, valueUpper = 1;
        int step = 0;
        void Update(bool bury)
        {
            double[] old = arrays[step % 2];
            double[] next = arrays[(step + 1) % 2];
            Array.Clear(next);
            long lower = long.MaxValue, upper = long.MinValue;
            for (long i = valueLower; i < valueUpper; i++)
            {
                double po = old[i - valueBegin];
                if (po == 0.0)
                    continue;
                for (int j = 0; j < score_.Length; j++)
                {
                    if (prob[j] == 0.0)
                        continue;
                    long value;
                    if (bury)
                        value = entry < i ? i : morgue;
                    else
                        value = morgue < i ? i + score_[j] : i;
                    double pr = po * prob[j];
                    if (value < valueBegin || valueBegin + capacity <= value)
                        continue;
                    if (value < lower)
                        lower = value;
                    if (upper < value)
                        upper = value;
                    next[value - valueBegin] += pr;
                }
            }

            valueLower = lower;
            valueUpper = upper + 1;
            step++;
        }

        double ProbAt(long value) => value < valueBegin || valueBegin + capacity <= value ? 0.0 : arrays[step % 2][value - valueBegin];
        Update(false);
        double eSumAlpha = 0.0, eOneMinusExpSumAlpha = 0.0;
        for (long w = 1; w < iter; w++)
        {
            for (long value = score_[0]; value <= entry; value++)
            {
                eSumAlpha += ProbAt(value) * value;
                eOneMinusExpSumAlpha += ProbAt(value) * (1.0 - Math.Exp(lambda * value));
            }

            Update(true);
            Update(false);
        }

        for (long value = score_[0]; value <= entry; value++)
        {
            eSumAlpha += ProbAt(value) * value;
            eOneMinusExpSumAlpha += ProbAt(value) * (1.0 - Math.Exp(lambda * value));
        }

        double ratio = eOneMinusExpSumAlpha / eSumAlpha;
        double k = mu * mu / thetaMinusDelta / muAssoc * ratio * ratio;
        return (lambda, k);
    }

    private static bool IsLogarithmic(long[] score, double[] prob)
    {
        for (int i = 1; i < score.Length; i++)
            if (score[i] <= score[i - 1])
                return false;
        double sum = 0;
        foreach (double p in prob)
        {
            if (p < 0.0 || p > 1.0)
                return false;
            sum += p;
        }

        if (Math.Abs(sum - 1.0) > RelTol * 1.0)
            return false;
        double mu = 0;
        for (int i = 0; i < score.Length; i++)
            mu += score[i] * prob[i];
        return mu < 0.0 && score[^1] > 0.0;
    }

    private static long Gcd(long i, long j)
    {
        long a = Math.Max(Math.Abs(i), Math.Abs(j)), b = Math.Min(Math.Abs(i), Math.Abs(j));
        while (b > 0)
            (a, b) = (b, a % b);
        return a;
    }

    /// <summary><c>Njn::Root::bisection</c>: the root of f(x) = y in [p, q] to within the absolute tolerance (100 iterations at most).</summary>
    private static double Bisection(double y, Func<double, double> f, double p, double q, double tol)
    {
        double fp = f(p) - y, fq = f(q) - y;
        if (fp * fq > 0.0)
            throw new LastProgramException("Root::bisection : root not bracketed");
        if (fp == 0.0)
            return p;
        if (fq == 0.0)
            return q;
        if (fp > 0.0)
            (p, q) = (q, p);
        double x = 0.5 * (p + q);
        for (int iterations = 100; iterations > 0; iterations--)
        {
            double fx = f(x) - y;
            if (fx < 0.0)
                p = x;
            else
                q = x;
            x = 0.5 * (p + q);
            if (Math.Abs(p - x) <= tol)
                return x;
        }

        return double.PositiveInfinity;
    }
}
