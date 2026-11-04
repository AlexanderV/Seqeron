using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// AnchorBasedAligner (suffix-tree seed / colinear chain / Needleman-Wunsch fill) and the exact
/// linear-space Hirschberg global aligner used for very large inter-anchor segments.
/// Test Unit: ALIGN-MULTI-001 (B09 review F5/F6; cross-batch request from B05).
///
/// Sources: Delcher et al. 1999, Nucleic Acids Res 27:2369 (MUMmer: suffix-tree matches, longest
/// colinear subset, DP between them); Abouelhoda &amp; Ohlebusch 2005, J Discrete Algorithms 3:321
/// (maximum-weight chain of colinear non-overlapping fragments, weight = length); Hirschberg 1975,
/// Commun ACM 18:341 (linear-space optimal alignment). Pair optima: Biopython 1.88 PairwiseAligner.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class AnchorBasedAligner_Tests
{
    private static AnchorBasedAligner.Anchor A(int r, int q, int len) => new(r, q, len);

    // ----- ChainAnchors -----------------------------------------------------------------------

    [Test]
    public void ChainAnchors_Empty_ReturnsEmpty()
    {
        Assert.That(AnchorBasedAligner.ChainAnchors(new List<AnchorBasedAligner.Anchor>()), Is.Empty);
    }

    [Test]
    public void ChainAnchors_DisjointColinear_AreAllKept()
    {
        var chain = AnchorBasedAligner.ChainAnchors(new() { A(0, 0, 10), A(20, 20, 10), A(40, 45, 8) });

        Assert.That(chain, Is.EqualTo(new[] { A(0, 0, 10), A(20, 20, 10), A(40, 45, 8) }));
    }

    [Test]
    public void ChainAnchors_OverlapInReference_TrimsTheLaterAnchor_KeepsBoth()
    {
        // ref [0,10) and [6,16); query [0,10) and [20,30). The second anchor overlaps the first by
        // 4 reference positions: its front is trimmed by 4 -> (10, 24, 6). Total covered = 16.
        // (Discarding one of the two, as the previous overlap filter did, would cover only 10.)
        var chain = AnchorBasedAligner.ChainAnchors(new() { A(0, 0, 10), A(6, 20, 10) });

        Assert.That(chain, Is.EqualTo(new[] { A(0, 0, 10), A(10, 24, 6) }));
    }

    [Test]
    public void ChainAnchors_OverlapInQuery_TrimsByTheLargerOverlap()
    {
        // ref [0,10) / [30,40); query [0,10) / [7,17): overlap 3 in the query only.
        var chain = AnchorBasedAligner.ChainAnchors(new() { A(0, 0, 10), A(30, 7, 10) });

        Assert.That(chain, Is.EqualTo(new[] { A(0, 0, 10), A(33, 10, 7) }));
    }

    [Test]
    public void ChainAnchors_PrefersTotalCoverageOverNumberOfAnchors()
    {
        // One 50-base anchor crosses three short colinear 3-base anchors; the weighted chain keeps
        // the long one (50 > 9). A longest-increasing-subsequence by count would keep the three.
        var chain = AnchorBasedAligner.ChainAnchors(new()
        {
            A(0, 100, 50), A(0, 0, 3), A(10, 5, 3), A(20, 10, 3)
        });

        Assert.That(chain, Is.EqualTo(new[] { A(0, 100, 50) }));
    }

    [Test]
    public void ChainAnchors_CrossingAnchors_KeepOnlyOneOfTheCrossingPair()
    {
        var chain = AnchorBasedAligner.ChainAnchors(new() { A(0, 30, 12), A(20, 0, 15) });

        Assert.That(chain, Is.EqualTo(new[] { A(20, 0, 15) }));
    }

    [Test]
    public void ChainAnchors_RandomAnchors_ChainIsStrictlyColinearNonOverlappingAndInsideInputs()
    {
        var rng = new Random(11);
        for (int trial = 0; trial < 300; trial++)
        {
            var anchors = new List<AnchorBasedAligner.Anchor>();
            for (int i = rng.Next(1, 14); i > 0; i--)
                anchors.Add(A(rng.Next(0, 60), rng.Next(0, 60), rng.Next(1, 25)));

            var chain = AnchorBasedAligner.ChainAnchors(anchors);

            Assert.That(chain, Is.Not.Empty);
            for (int i = 0; i < chain.Count; i++)
            {
                // Every chained piece is a (front-trimmed) piece of one input anchor on its diagonal.
                Assert.That(anchors.Any(o => o.RefStart - o.QueryStart == chain[i].RefStart - chain[i].QueryStart
                                              && o.RefStart <= chain[i].RefStart && chain[i].RefEnd == o.RefEnd), Is.True);
                if (i > 0)
                {
                    Assert.That(chain[i].RefStart, Is.GreaterThanOrEqualTo(chain[i - 1].RefEnd));
                    Assert.That(chain[i].QueryStart, Is.GreaterThanOrEqualTo(chain[i - 1].QueryEnd));
                }
            }

            // At least as good as the best single anchor.
            Assert.That(chain.Sum(c => c.Length), Is.GreaterThanOrEqualTo(anchors.Max(a => a.Length)));
        }
    }

    // ----- AlignWithAnchors -------------------------------------------------------------------

    private static AlignmentResult Anchored(string reference, string query, ScoringMatrix? sc = null) =>
        AnchorBasedAligner.AlignWithAnchors(reference, SuffixTree.SuffixTree.Build(reference), query, sc ?? SequenceAligner.SimpleDna);

    private static int LinearScore(AlignmentResult r, ScoringMatrix s)
    {
        int total = 0;
        for (int i = 0; i < r.AlignedSequence1.Length; i++)
        {
            char a = r.AlignedSequence1[i], b = r.AlignedSequence2[i];
            total += a == '-' || b == '-' ? s.GapExtend : a == b ? s.Match : s.Mismatch;
        }
        return total;
    }

    [Test]
    public void AlignWithAnchors_Identical_IsAllMatches()
    {
        var r = Anchored("ACGTACGTACGTTTGA", "ACGTACGTACGTTTGA");

        Assert.That(r.AlignedSequence1, Is.EqualTo("ACGTACGTACGTTTGA"));
        Assert.That(r.AlignedSequence2, Is.EqualTo("ACGTACGTACGTTTGA"));
        Assert.That(r.Score, Is.EqualTo(16));
    }

    [Test]
    public void AlignWithAnchors_EmptySide_IsAllGapAlignment_NotEmptyResult()
    {
        var r = Anchored("ACGT", "");

        Assert.That(r.AlignedSequence1, Is.EqualTo("ACGT"));
        Assert.That(r.AlignedSequence2, Is.EqualTo("----"));
        Assert.That(r.Score, Is.EqualTo(-4)); // linear model: 4 x GapExtend, GapOpen unused
    }

    [Test]
    public void AlignWithAnchors_OneSidedGapRun_ScoredLinearly_LikeGlobalAlign()
    {
        // Anchor-flanked pure deletion of 3 bases; the old code charged GapOpen + 3*GapExtend (-5)
        // for such runs but GlobalAlign's linear -3 for DP gaps, so Score disagreed with the rows.
        string reference = "ACGTACGTAGGGCATCGATCGAT";
        string query = "ACGTACGTACATCGATCGAT"; // 'GGG' deleted
        var r = Anchored(reference, query);

        Assert.That(r.Score, Is.EqualTo(LinearScore(r, SequenceAligner.SimpleDna)));
        Assert.That(r.Score, Is.EqualTo(SequenceAligner.GlobalAlign(reference, query).Score));
    }

    [Test]
    public void AlignWithAnchors_RelatedPairs_ValidAndNeverAboveGlobalOptimum()
    {
        var rng = new Random(4242);
        var sc = SequenceAligner.SimpleDna;
        for (int trial = 0; trial < 120; trial++)
        {
            string reference = RandomDna(rng, rng.Next(30, 260));
            string query = Mutate(rng, reference, 0.06 + 0.2 * rng.NextDouble());
            if (query.Length == 0) query = "A";

            var r = Anchored(reference, query, sc);
            var opt = SequenceAligner.GlobalAlign(reference, query, sc);

            Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(reference));
            Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(query));
            Assert.That(r.AlignedSequence2, Has.Length.EqualTo(r.AlignedSequence1.Length));
            Assert.That(r.Score, Is.EqualTo(LinearScore(r, sc)));
            Assert.That(r.Score, Is.LessThanOrEqualTo(opt.Score));
        }
    }

    [Test]
    public void AlignWithAnchors_RepeatedBlocks_OverlappingAnchors_StillValid()
    {
        // Tandem repeat: the query's blocks match the reference on several diagonals, so the
        // suffix-link anchors overlap each other.
        string unit = "ACGGTCATGCAA";
        string reference = unit + unit + unit + "TTGACCA" + unit;
        string query = unit + "G" + unit + "TTGACCA" + unit + unit;

        var r = Anchored(reference, query);

        Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(reference));
        Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(query));
        Assert.That(r.Score, Is.EqualTo(LinearScore(r, SequenceAligner.SimpleDna)));
        Assert.That(r.Score, Is.LessThanOrEqualTo(SequenceAligner.GlobalAlign(reference, query).Score));
    }

    // ----- GlobalAlignLinearSpace (Hirschberg) ------------------------------------------------

    [TestCase("GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT", "AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT", 24)]
    [TestCase("TAGACGCATCAAACAGCCTAGGCCGTACTGCGATTCAGCTGTCAACAGGTAAAACAGA", "TAGACGCCTAGTCAACAGCCTACGGCCCTATCAGCTGTATACATGTTAGTCCACCGTAACAGT", 27)]
    public void GlobalAlignLinearSpace_ReachesBiopythonOptimum_EvenWhenForcedToRecurse(string a, string b, int biopythonOptimum)
    {
        var r = SequenceAligner.GlobalAlignLinearSpace(a, b, SequenceAligner.SimpleDna, maxFullMatrixCells: 20);

        Assert.That(r.Score, Is.EqualTo(biopythonOptimum));
        Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(a));
        Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(b));
        Assert.That(LinearScore(r, SequenceAligner.SimpleDna), Is.EqualTo(biopythonOptimum));
    }

    [Test]
    public void GlobalAlignLinearSpace_RandomPairs_ScoreEqualsNeedlemanWunsch_AllSchemes()
    {
        var rng = new Random(99);
        var schemes = new[] { SequenceAligner.SimpleDna, SequenceAligner.BlastDna, SequenceAligner.HighIdentityDna, new ScoringMatrix(3, -2, -4, -3) };
        for (int trial = 0; trial < 200; trial++)
        {
            string a = RandomDna(rng, rng.Next(0, 70));
            string b = rng.Next(2) == 0 ? Mutate(rng, a, 0.25) : RandomDna(rng, rng.Next(0, 70));
            var sc = schemes[trial % schemes.Length];

            var nw = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b), sc);
            var h = SequenceAligner.GlobalAlignLinearSpace(a, b, sc, maxFullMatrixCells: 16);

            Assert.That(h.Score, Is.EqualTo(nw.Score), $"{a} / {b}");
            Assert.That(LinearScore(h, sc), Is.EqualTo(h.Score));
            Assert.That(h.AlignedSequence1.Replace("-", ""), Is.EqualTo(a));
            Assert.That(h.AlignedSequence2.Replace("-", ""), Is.EqualTo(b));
        }
    }

    private static string RandomDna(Random rng, int n)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[rng.Next(4)];
        return new string(c);
    }

    private static string Mutate(Random rng, string s, double p)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char ch in s)
        {
            double x = rng.NextDouble();
            if (x < p) sb.Append("ACGT"[rng.Next(4)]);
            else if (x < p * 1.5) { }
            else if (x < p * 2) { sb.Append(ch); sb.Append("ACGT"[rng.Next(4)]); }
            else sb.Append(ch);
        }
        return sb.ToString();
    }
}
