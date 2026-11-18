using System.Diagnostics;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-16, B09 F54): <c>AnchorBasedAligner.ChainAnchors</c> is a range-maximum sparse dynamic program
/// (Abouelhoda &amp; Ohlebusch 2005) instead of the all-pairs O(a²) scan. The O(a²) scan is kept here, verbatim, as the oracle:
/// the chain, its total score and every tie decision must be identical on random, tie-heavy and adversarial anchor sets.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class AnchorBasedAligner_ChainEquivalence_Tests
{
    private static AnchorBasedAligner.Anchor A(int r, int q, int len) => new(r, q, len);

    /// <summary>The removed O(a²) <c>ChainAnchors</c> dynamic program, verbatim (test oracle).</summary>
    private static List<AnchorBasedAligner.Anchor> ChainAnchorsQuadratic(List<AnchorBasedAligner.Anchor> anchors)
    {
        var sorted = anchors
            .Where(a => a.Length > 0 && a.RefStart >= 0 && a.QueryStart >= 0)
            .OrderBy(a => a.RefStart).ThenBy(a => a.QueryStart).ThenByDescending(a => a.Length)
            .ToList();
        int n = sorted.Count;
        if (n == 0)
            return new List<AnchorBasedAligner.Anchor>();

        var best = new long[n];
        var prev = new int[n];
        var trimOf = new int[n];
        for (int i = 0; i < n; i++)
        {
            best[i] = sorted[i].Length;
            prev[i] = -1;
            for (int j = 0; j < i; j++)
            {
                int trim = Math.Max(0, Math.Max(
                    sorted[j].RefEnd - sorted[i].RefStart,
                    sorted[j].QueryEnd - sorted[i].QueryStart));
                int gain = sorted[i].Length - trim;
                if (gain <= 0)
                    continue;
                if (best[j] + gain > best[i])
                {
                    best[i] = best[j] + gain;
                    prev[i] = j;
                    trimOf[i] = trim;
                }
            }
        }

        int bestIdx = 0;
        for (int i = 1; i < n; i++)
            if (best[i] > best[bestIdx])
                bestIdx = i;

        var chain = new List<AnchorBasedAligner.Anchor>();
        for (int idx = bestIdx; idx >= 0; idx = prev[idx])
        {
            var a = sorted[idx];
            int t = prev[idx] >= 0 ? trimOf[idx] : 0;
            chain.Add(new AnchorBasedAligner.Anchor(a.RefStart + t, a.QueryStart + t, a.Length - t));
        }

        chain.Reverse();
        return chain;
    }

    private static void AssertSameChain(List<AnchorBasedAligner.Anchor> anchors, string label)
    {
        var expected = ChainAnchorsQuadratic(anchors);
        var actual = AnchorBasedAligner.ChainAnchors(anchors);
        Assert.That(actual, Is.EqualTo(expected), label);
    }

    private static List<AnchorBasedAligner.Anchor> Random(Random rng, int count, int coord, int maxLen)
    {
        var list = new List<AnchorBasedAligner.Anchor>(count);
        for (int i = 0; i < count; i++)
            list.Add(A(rng.Next(0, coord), rng.Next(0, coord), rng.Next(1, maxLen)));
        return list;
    }

    [Test]
    public void ChainAnchors_EqualsQuadraticOracle_OnRandomSets_DenseAndSparse()
    {
        var rng = new Random(20261118);
        for (int trial = 0; trial < 1500; trial++)
        {
            // coordinate range / length scales vary from "everything overlaps and ties" to "mostly disjoint"
            int coord = new[] { 6, 15, 40, 120, 600 }[trial % 5];
            int maxLen = new[] { 4, 9, 25, 80 }[(trial / 5) % 4];
            AssertSameChain(Random(rng, rng.Next(1, 60), coord, maxLen), $"trial {trial} coord {coord} maxLen {maxLen}");
        }
    }

    [Test]
    public void ChainAnchors_EqualsQuadraticOracle_OnMemLikeSets_ManyDiagonalsAndRepeats()
    {
        // Shapes of real MEM sets: a colinear backbone plus anchors on shifted diagonals (repeats) and inverted pieces.
        var rng = new Random(77);
        for (int trial = 0; trial < 400; trial++)
        {
            var list = new List<AnchorBasedAligner.Anchor>();
            int r = 0, q = 0;
            for (int k = rng.Next(3, 40); k > 0; k--)
            {
                int len = rng.Next(8, 40);
                list.Add(A(r, q, len));
                if (rng.Next(3) == 0)
                    list.Add(A(r + rng.Next(-5, 10), q + rng.Next(-20, 20), rng.Next(8, 30)));   // repeat copy on another diagonal
                if (rng.Next(6) == 0)
                    list.Add(A(rng.Next(0, 400), rng.Next(0, 400), rng.Next(8, 30)));            // noise
                r += len + rng.Next(0, 15);
                q += len + rng.Next(0, 15);
            }

            AssertSameChain(list, $"trial {trial}");
        }
    }

    [Test]
    public void ChainAnchors_EqualsQuadraticOracle_OnAdversarialShapes()
    {
        var cases = new Dictionary<string, List<AnchorBasedAligner.Anchor>>
        {
            ["all identical"] = Enumerable.Range(0, 20).Select(_ => A(5, 5, 12)).ToList(),
            ["same start, different lengths"] = Enumerable.Range(1, 15).Select(l => A(10, 10, l)).ToList(),
            ["nested in one long anchor"] = new[] { A(0, 0, 200) }.Concat(Enumerable.Range(0, 30).Select(i => A(5 + 5 * i, 5 + 5 * i, 4))).ToList(),
            ["inverted (anti-diagonal)"] = Enumerable.Range(0, 40).Select(i => A(10 * i, 400 - 10 * i, 8)).ToList(),
            ["staircase of overlaps"] = Enumerable.Range(0, 50).Select(i => A(7 * i, 7 * i, 10)).ToList(),
            ["shifted diagonals, every pair overlaps"] = Enumerable.Range(0, 40).Select(i => A(i, 2 * i, 60)).ToList(),
            ["query-overlap only"] = Enumerable.Range(0, 30).Select(i => A(20 * i, 12 * i, 15)).ToList(),
            ["reference-overlap only"] = Enumerable.Range(0, 30).Select(i => A(12 * i, 20 * i, 15)).ToList(),
            ["duplicates of a two-anchor chain"] = new[] { A(0, 0, 9), A(20, 20, 9), A(0, 0, 9), A(20, 20, 9), A(10, 10, 9) }.ToList(),
            ["single"] = new[] { A(3, 4, 5) }.ToList(),
            ["empty"] = new List<AnchorBasedAligner.Anchor>(),
            ["invalid entries filtered"] = new[] { A(-1, 0, 5), A(0, -3, 5), A(2, 2, 0), A(4, 4, 6) }.ToList(),
        };

        foreach (var (label, list) in cases)
            AssertSameChain(list, label);
    }

    [Test]
    public void ChainAnchors_EqualsQuadraticOracle_OnLargeSets()
    {
        var rng = new Random(5);
        AssertSameChain(Random(rng, 3000, 20000, 60), "3000 anchors, sparse");
        AssertSameChain(Random(rng, 3000, 2000, 120), "3000 anchors, dense overlaps");
    }

    [Test]
    public void ChainAnchors_HundredThousandAnchors_IsFastAndStrictlyColinear()
    {
        // 10^5 anchors (the O(a²) scan would need ~5·10^9 predecessor tests): a colinear backbone with repeat copies and noise.
        var rng = new Random(2026);
        var list = new List<AnchorBasedAligner.Anchor>(100_000);
        int r = 0, q = 0;
        while (list.Count < 100_000)
        {
            int len = rng.Next(8, 60);
            list.Add(A(r, q, len));
            if (rng.Next(3) == 0)
                list.Add(A(r + rng.Next(-5, 30), q + rng.Next(-30, 30), rng.Next(8, 50)));
            r += len + rng.Next(0, 20);
            q += len + rng.Next(0, 20);
        }

        var sw = Stopwatch.StartNew();
        var chain = AnchorBasedAligner.ChainAnchors(list);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(30), $"took {sw.Elapsed.TotalSeconds:F1} s");
            Assert.That(chain.Count, Is.GreaterThan(30_000));
            for (int i = 1; i < chain.Count; i++)
            {
                Assert.That(chain[i].RefStart, Is.GreaterThanOrEqualTo(chain[i - 1].RefEnd));
                Assert.That(chain[i].QueryStart, Is.GreaterThanOrEqualTo(chain[i - 1].QueryEnd));
            }
        });
    }
}
