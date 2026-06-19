using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Phylogenetics;
using static Seqeron.Genomics.Phylogenetics.PhylogeneticAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Phylogenetic area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-DIST-001 — pairwise evolutionary distance (Phylogenetic).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 39.
///
/// API under test (PhylogeneticAnalyzer.CalculatePairwiseDistance):
///   Over the columns where both aligned sequences carry a standard base (A/C/G/T; gaps and
///   ambiguous bases skipped), let p = differences / comparableSites. Then
///     Hamming         = raw difference count
///     PDistance       = p
///     JukesCantor     = −¾·ln(1 − 4p/3)                          (JC69)
///     Kimura2Parameter= −½·ln((1 − 2S − V)·√(1 − 2V))            (K80; S,V = transition,
///                                                                  transversion proportions)
///   comparableSites = 0 ⇒ distance 0.
///
/// Relations (derived from these formulas, NOT from output):
///   • INV (zero self-distance): every column of d(x,x) is identical, so differences = 0 and
///          p = 0; each formula maps p = 0 (and S = V = 0) to exactly 0.
///   • SYM (symmetry): difference, transition and transversion counts are symmetric in the two
///          sequences, so d(a,b) = d(b,a) for every method.
///   • MON (more mutations ⇒ larger distance): adding one substitution raises the difference
///          count (and S or V), and each method is strictly increasing in those — so the
///          distance strictly increases.
///   • COMP (triangle inequality): Hamming and PDistance are genuine metrics on equal-length
///          gap-free sequences (a≠c at a site implies a≠b or b≠c there), so
///          d(a,c) ≤ d(a,b) + d(b,c). The corrected distances (JukesCantor, Kimura2Parameter)
///          are CONVEX in p and designed for additivity along a tree, NOT metricity, and can
///          violate the triangle inequality — so this relation is asserted only for the two
///          raw metrics, rather than rubber-stamping a property the corrections lack.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class PhylogeneticMetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(20260619);

    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>A base guaranteed to differ from <paramref name="c"/>, so flipping yields a difference.</summary>
    private static char Flip(char c) => c == 'A' ? 'C' : 'A';

    /// <summary>Returns a copy of <paramref name="seq"/> with exactly <paramref name="count"/> spread-out positions mutated.</summary>
    private static string MutatePositions(string seq, int count)
    {
        var chars = seq.ToCharArray();
        for (int n = 0; n < count; n++)
        {
            int pos = (int)((long)n * seq.Length / count);
            chars[pos] = Flip(chars[pos]);
        }
        return new string(chars);
    }

    private static readonly DistanceMethod[] AllMethods =
    {
        DistanceMethod.PDistance,
        DistanceMethod.JukesCantor,
        DistanceMethod.Kimura2Parameter,
        DistanceMethod.Hamming,
    };

    /// <summary>Methods that are genuine metrics on equal-length gap-free sequences.</summary>
    private static readonly DistanceMethod[] MetricMethods =
    {
        DistanceMethod.PDistance,
        DistanceMethod.Hamming,
    };

    #endregion

    #region INV — d(x,x) = 0

    [Test]
    [Description("INV: a sequence has zero distance to itself under every method, because no column differs.")]
    public void PairwiseDistance_SelfDistance_IsZero()
    {
        foreach (var method in AllMethods)
        {
            foreach (int len in new[] { 10, 40, 100 })
            {
                string x = RandomDna(len);
                PhylogeneticAnalyzer.CalculatePairwiseDistance(x, x, method)
                    .Should().BeApproximately(0.0, 1e-12,
                        because: $"{method}: identical sequences differ at no comparable site, so p = 0 maps to distance 0");
            }
        }
    }

    #endregion

    #region SYM — d(a,b) = d(b,a)

    [Test]
    [Description("SYM: distance is symmetric, since difference/transition/transversion counts do not depend on argument order.")]
    public void PairwiseDistance_IsSymmetric()
    {
        foreach (var method in AllMethods)
        {
            for (int trial = 0; trial < 20; trial++)
            {
                int len = 60;
                string a = RandomDna(len);
                string b = MutatePositions(a, Rng.Next(1, 10));   // a few substitutions

                double ab = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
                double ba = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, a, method);

                ba.Should().BeApproximately(ab, 1e-12,
                    because: $"{method}: swapping the arguments leaves the per-column comparison counts unchanged");
            }
        }
    }

    #endregion

    #region MON — more substitutions strictly increase the distance

    [Test]
    [Description("MON: each additional substitution raises the difference count, and every method is strictly increasing in it, so the distance strictly increases.")]
    public void PairwiseDistance_MoreMutations_StrictlyIncreases()
    {
        const int len = 40;   // keep mutation fraction well below the JC/K80 divergence limits

        foreach (var method in AllMethods)
        {
            string baseSeq = RandomDna(len);
            double previous = double.NegativeInfinity;

            for (int d = 0; d <= 6; d++)
            {
                string variant = MutatePositions(baseSeq, d);
                double dist = PhylogeneticAnalyzer.CalculatePairwiseDistance(baseSeq, variant, method);

                dist.Should().BeGreaterThan(previous,
                    because: $"{method}: {d} substitutions differ at one more site than {d - 1}, and the distance is strictly increasing in the difference count");
                previous = dist;
            }
        }
    }

    #endregion

    #region COMP — triangle inequality for the raw metrics

    [Test]
    [Description("COMP: on equal-length gap-free sequences Hamming and PDistance satisfy the triangle inequality d(a,c) ≤ d(a,b) + d(b,c).")]
    public void PairwiseDistance_RawMetrics_SatisfyTriangleInequality()
    {
        foreach (var method in MetricMethods)
        {
            for (int trial = 0; trial < 50; trial++)
            {
                int len = 50;
                string a = RandomDna(len);
                string b = RandomDna(len);
                string c = RandomDna(len);

                double ab = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
                double bc = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, c, method);
                double ac = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, c, method);

                ac.Should().BeLessThanOrEqualTo(ab + bc + 1e-9,
                    because: $"{method}: a≠c at a site forces a≠b or b≠c there, so per-site disagreements are subadditive");
            }
        }
    }

    #endregion
}
