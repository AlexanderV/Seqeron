using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// Tests for GenomeAnnotator.FindPromoterMotifs — bacterial σ70 promoter boxes.
/// Test Unit: ANNOT-PROM-001
///
/// Definition (Harley &amp; Reynolds 1987, NAR 15:2343; Hawley &amp; McClure 1983, NAR 11:2237): a σ70 promoter is a
/// −35 hexamer (consensus TTGACA) followed after a 15–21 bp spacer (17 ± 1 bp in 92 %) by a −10 hexamer (consensus
/// TATAAT). FindPromoterMotifs delegates the pairing to the canonical MotifFinder.FindSigma70Promoters and reports the
/// boxes of every pair once; score = 1 − d/6 (d = Hamming distance of the box to its consensus).
/// Oracle: independent Python brute force of that definition (400 random cases, 1,131 boxes, 0 mismatches —
/// docs/Evidence/ANNOT-PROM-001-Evidence.md).
/// </summary>
[TestFixture]
public class GenomeAnnotator_PromoterMotif_Tests
{
    private static readonly string C17 = new('C', 17);

    // E. coli lac operon control region with the lacUV5 −10 (TATAAT): −35 TTTACA at 68, 18-bp spacer, −10 at 92
    // (same published fixture as MotifFinder_Sigma70Promoters_Tests).
    private const string LacUv5 =
        "TCAGCATTCGAGCTTACGGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT";

    private static List<(int position, string type, string sequence, double score)> Find(string s) =>
        GenomeAnnotator.FindPromoterMotifs(s).ToList();

    #region Consensus promoter

    /// <summary>M01: consensus −35 + 17-bp spacer + consensus −10 → exactly the two boxes, score 1.0.</summary>
    [Test]
    public void FindPromoterMotifs_ConsensusPromoter_ReportsBothBoxes()
    {
        var motifs = Find("GG" + "TTGACA" + C17 + "TATAAT" + "GG");

        Assert.That(motifs, Is.EqualTo(new[]
        {
            (2, "-35 box", "TTGACA", 1.0),
            (25, "-10 box", "TATAAT", 1.0),
        }));
    }

    /// <summary>M02/M03: an isolated −35 or −10 consensus with no partner box is not a promoter (H&amp;R pairing).</summary>
    [TestCase("GGGGGTTGACAGGGGG")]
    [TestCase("CCCCCTATAATCCCCC")]
    [TestCase("TTGACATATAAT")]
    public void FindPromoterMotifs_UnpairedBox_NotReported(string sequence)
    {
        Assert.That(Find(sequence), Is.Empty);
    }

    #endregion

    #region Spacer

    /// <summary>M04: spacer bounds 15–21 bp (Harley &amp; Reynolds 1987).</summary>
    [TestCase(14, false)]
    [TestCase(15, true)]
    [TestCase(17, true)]
    [TestCase(21, true)]
    [TestCase(22, false)]
    public void FindPromoterMotifs_SpacerBounds(int spacer, bool expected)
    {
        var motifs = Find("TTGACA" + new string('C', spacer) + "TATAAT");

        if (expected)
            Assert.That(motifs, Is.EqualTo(new[]
            {
                (0, "-35 box", "TTGACA", 1.0),
                (6 + spacer, "-10 box", "TATAAT", 1.0),
            }));
        else
            Assert.That(motifs, Is.Empty);
    }

    /// <summary>M05: overload exposes the spacer range.</summary>
    [Test]
    public void FindPromoterMotifs_CustomSpacerRange()
    {
        string seq = "TTGACA" + new string('C', 22) + "TATAAT";
        Assert.That(GenomeAnnotator.FindPromoterMotifs(seq, 0, 0, 15, 22).Select(m => m.position), Is.EqualTo(new[] { 0, 28 }));
        Assert.That(GenomeAnnotator.FindPromoterMotifs("GG" + "TTGACA" + C17 + "TATAAT", 0, 0, 15, 16), Is.Empty);
    }

    #endregion

    #region Mismatches and score

    /// <summary>M06: score = 1 − d/6 per box (1 mismatch → 5/6).</summary>
    [Test]
    public void FindPromoterMotifs_OneMismatchPerBox_ScoreFiveSixths()
    {
        var motifs = Find("GG" + "TTTACA" + C17 + "TATGAT" + "GG");

        Assert.That(motifs, Is.EqualTo(new[]
        {
            (2, "-35 box", "TTTACA", 5.0 / 6),
            (25, "-10 box", "TATGAT", 5.0 / 6),
        }));
    }

    /// <summary>M07: default limit ≤ 2 mismatches per box; the overload raises it.</summary>
    [Test]
    public void FindPromoterMotifs_MismatchLimit()
    {
        string seq = "GG" + "TGGCCA" + C17 + "TATAAT" + "GG"; // TGGCCA: 2 mismatches; GTGGCC at 1: 3 mismatches

        Assert.That(Find(seq), Is.EqualTo(new[]
        {
            (2, "-35 box", "TGGCCA", 4.0 / 6),
            (25, "-10 box", "TATAAT", 1.0),
        }));
        Assert.That(GenomeAnnotator.FindPromoterMotifs(seq, 3, 0, 15, 21).ToList(), Is.EqualTo(new[]
        {
            (1, "-35 box", "GTGGCC", 3.0 / 6),
            (2, "-35 box", "TGGCCA", 4.0 / 6),
            (25, "-10 box", "TATAAT", 1.0),
        }));
    }

    /// <summary>M08: a box paired with several partners is reported once; output ordered by position.</summary>
    [Test]
    public void FindPromoterMotifs_SharedBox_ReportedOnce()
    {
        string seq = "TTGACA" + new string('C', 15) + "TATAAT" + "TATAAT"; // −10 at spacer 15 and 21

        Assert.That(GenomeAnnotator.FindPromoterMotifs(seq, 0, 0, 15, 21).ToList(), Is.EqualTo(new[]
        {
            (0, "-35 box", "TTGACA", 1.0),
            (21, "-10 box", "TATAAT", 1.0),
            (27, "-10 box", "TATAAT", 1.0),
        }));
        Assert.That(Find(seq), Is.EqualTo(new[]
        {
            (0, "-35 box", "TTGACA", 1.0),
            (21, "-10 box", "TATAAT", 1.0),
            (24, "-10 box", "AATTAT", 4.0 / 6),
            (27, "-10 box", "TATAAT", 1.0),
        }));
    }

    #endregion

    #region Published / oracle examples

    /// <summary>M09: lacUV5 — −35 TTTACA at 68, 18-bp spacer, −10 TATAAT at 92 (published lac control region).</summary>
    [Test]
    public void FindPromoterMotifs_LacUv5_PublishedBoxes()
    {
        Assert.That(GenomeAnnotator.FindPromoterMotifs(LacUv5, 1, 0, 15, 21).ToList(), Is.EqualTo(new[]
        {
            (68, "-35 box", "TTTACA", 5.0 / 6),
            (92, "-10 box", "TATAAT", 1.0),
        }));
        // Defaults (≤ 2 mismatches) — Python oracle.
        Assert.That(Find(LacUv5), Is.EqualTo(new[]
        {
            (54, "-35 box", "TAGGCA", 4.0 / 6),
            (68, "-35 box", "TTTACA", 5.0 / 6),
            (77, "-10 box", "TATGCT", 4.0 / 6),
            (92, "-10 box", "TATAAT", 1.0),
        }));
    }

    /// <summary>M10: random AT-rich sequence, defaults — exact box set from the independent Python oracle.</summary>
    [Test]
    public void FindPromoterMotifs_RandomSequence_MatchesPythonOracle()
    {
        Assert.That(Find("CCGATATTTTAAAAGTTTTTTTAAAGACTAACTATAAATTTGACATTATG"), Is.EqualTo(new[]
        {
            (7, "-35 box", "TTTAAA", 4.0 / 6),
            (8, "-35 box", "TTAAAA", 4.0 / 6),
            (19, "-35 box", "TTTAAA", 4.0 / 6),
            (32, "-10 box", "TATAAA", 5.0 / 6),
            (34, "-10 box", "TAAATT", 4.0 / 6),
            (43, "-10 box", "CATTAT", 4.0 / 6),
        }));
    }

    /// <summary>M11: delegation — the boxes are exactly those of MotifFinder.FindSigma70Promoters' pairs.</summary>
    [Test]
    public void FindPromoterMotifs_EqualsBoxesOfFindSigma70Promoters()
    {
        var rng = new Random(11);
        for (int k = 0; k < 50; k++)
        {
            string seq = new(Enumerable.Range(0, 120).Select(_ => "AATTCG"[rng.Next(6)]).ToArray());
            var expected = MotifFinder.FindSigma70Promoters(new DnaSequence(seq))
                .SelectMany(c => new[] { (c.Minus35Start, "-35 box", c.Minus35), (c.Minus10Start, "-10 box", c.Minus10) })
                .Distinct().OrderBy(b => b.Item1).ThenBy(b => b.Item2 == "-35 box" ? 0 : 1).ToList();
            var actual = Find(seq).Select(m => (m.position, m.type, m.sequence)).ToList();
            Assert.That(actual, Is.EqualTo(expected), seq);
        }
    }

    #endregion

    #region Edge cases

    /// <summary>M12: forward strand only — a reverse-strand promoter is not reported.</summary>
    [Test]
    public void FindPromoterMotifs_ReverseStrandPromoter_NotReported()
    {
        string plus = "TTGACA" + C17 + "TATAAT";
        Assert.That(Find("AA" + DnaSequence.GetReverseComplementString(plus) + "AA"), Is.Empty);
    }

    /// <summary>M13: case-insensitive.</summary>
    [Test]
    public void FindPromoterMotifs_LowerCase_SameAsUpperCase()
    {
        string seq = "GG" + "TTTACA" + C17 + "TATGAT" + "GG";
        Assert.That(Find(seq.ToLowerInvariant()), Is.EqualTo(Find(seq)));
    }

    /// <summary>M14: empty / too-short input → empty.</summary>
    [TestCase("")]
    [TestCase("TATA")]
    [TestCase("TTGACATATAA")]
    public void FindPromoterMotifs_EmptyOrShort_ReturnsEmpty(string sequence)
    {
        Assert.That(Find(sequence), Is.Empty);
    }

    /// <summary>M15: null → ArgumentNullException; non-ACGT → ArgumentException; invalid limits → ArgumentOutOfRangeException.</summary>
    [Test]
    public void FindPromoterMotifs_InvalidInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.FindPromoterMotifs(null!));
        Assert.Throws<ArgumentException>(() => GenomeAnnotator.FindPromoterMotifs("TTGACANNN"));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAnnotator.FindPromoterMotifs("ACGT", 7, 2, 15, 21));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeAnnotator.FindPromoterMotifs("ACGT", 2, 2, 21, 15));
    }

    #endregion
}
