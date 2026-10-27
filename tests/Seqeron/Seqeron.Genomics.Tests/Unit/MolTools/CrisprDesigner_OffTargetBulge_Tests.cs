namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>Cas-OFFinder 3.0.0b3 oracle rows are recorded in CRISPR-OFF-001-Evidence.md.</summary>
[TestFixture]
public class CrisprDesigner_OffTargetBulge_Tests
{
    private const string Guide = "GGCCGACCTGTCGCTGACGC";
    private const string Toy =
        "TTTTTTTTTT" +
        "GGCCGACCATGTCGCTGACGCTGG" + "TTTTTTTTTTTTTTT" +
        "GGCCGACCGTCGCTGACGCTGG" + "TTTTTTTTTTTTTTT" +
        "CCAGCGTCAGCGACATGGTCGGCC" + "TTTTTTTTTTTTTTT" +
        "CCAGCGTCAGCGACGGTCGGCC" + "TTTTTTTTTT";

    [Test]
    public void FindOffTargetsWithBulges_CasOffinderOracle_FourExactBulgesOnBothStrands()
    {
        var hits = CrisprDesigner.FindOffTargetsWithBulges(
            Guide, new DnaSequence(Toy), 0, 2, 2).ToArray();

        Assert.That(hits, Has.Length.EqualTo(4));
        Assert.Multiple(() =>
        {
            AssertHit(hits, 31, 10, true, BulgeType.Dna,
                "GGCCGACC-TGTCGCTGACGC", "GGCCGACCATGTCGCTGACGC", 8);
            AssertHit(hits, 68, 49, true, BulgeType.Rna,
                "GGCCGACCTGTCGCTGACGC", "GGCCGACC-GTCGCTGACGC", 8);
            AssertHit(hits, 86, 89, false, BulgeType.Dna,
                "GGCCGACC-TGTCGCTGACGC", "GGCCGACCATGTCGCTGACGC", 8);
            AssertHit(hits, 125, 128, false, BulgeType.Rna,
                "GGCCGACCTGTCGCTGACGC", "GGCCGACC-GTCGCTGACGC", 8);
        });
    }

    private static void AssertHit(BulgedOffTargetSite[] hits, int pam, int targetStart,
        bool forward, BulgeType type, string guide, string target, int gap)
    {
        var hit = hits.Single(h => h.Position == pam && h.IsForwardStrand == forward);
        Assert.That(hit.TargetStart, Is.EqualTo(targetStart));
        Assert.That(hit.Bulge, Is.EqualTo(type));
        Assert.That(hit.BulgeSize, Is.EqualTo(1));
        Assert.That(hit.BulgePosition, Is.EqualTo(gap));
        Assert.That(hit.AlignedGuide, Is.EqualTo(guide));
        Assert.That(hit.AlignedTarget, Is.EqualTo(target));
        Assert.That(hit.Mismatches, Is.Zero);
    }

    [Test]
    public void FindOffTargetsWithBulges_CasOffinderOracle_OneMismatchOutsideDnaBulge()
    {
        const string sequence = "TTTTTTTTTTGGCCGACCATGTCGCAGACGCTGGTTTTTTTTTT";
        var hits = CrisprDesigner.FindOffTargetsWithBulges(
            Guide, new DnaSequence(sequence), 1, 1, 0).ToArray();

        Assert.That(hits, Has.Length.EqualTo(1));
        Assert.That(hits[0], Is.EqualTo(new BulgedOffTargetSite(
            Position: 31, TargetStart: 10,
            AlignedGuide: "GGCCGACC-TGTCGCTGACGC",
            AlignedTarget: "GGCCGACCATGTCGCAGACGC",
            Mismatches: 1, Bulge: BulgeType.Dna, BulgeSize: 1,
            BulgePosition: 8, IsForwardStrand: true, Pam: "TGG")));
    }

    [Test]
    public void FindOffTargetsWithBulges_CasOffinderOracle_TwoBaseBulges()
    {
        const string sequence =
            "TTTTTTTTTTGGCCGACCATTGTCGCTGACGCTGGTTTTTTTTTTTTTTT" +
            "GGCCGACCTCGCTGACGCTGGTTTTTTTTTT";
        var hits = CrisprDesigner.FindOffTargetsWithBulges(
            Guide, new DnaSequence(sequence), 0, 2, 2).ToArray();

        Assert.That(hits, Has.Length.EqualTo(3));
        Assert.That(hits.Select(h => (h.Position, h.TargetStart, h.Bulge,
            h.BulgeSize, h.BulgePosition, h.AlignedGuide, h.AlignedTarget, h.Mismatches)),
            Is.EquivalentTo(new[]
            {
                (32, 10, BulgeType.Dna, 2, 8,
                    "GGCCGACC--TGTCGCTGACGC", "GGCCGACCATTGTCGCTGACGC", 0),
                (68, 50, BulgeType.Rna, 2, 8,
                    "GGCCGACCTGTCGCTGACGC", "GGCCGACC--TCGCTGACGC", 0),
                (68, 50, BulgeType.Rna, 2, 9,
                    "GGCCGACCTGTCGCTGACGC", "GGCCGACCT--CGCTGACGC", 0)
            }));
    }

    // -------------------------------------------------------------------------------------------
    // Admissible gap offsets (audit round 4, H2/H3). Oracle: a source-faithful port of Cas-OFFinder
    // 3.0.0b3's compare-string enumeration (cas-offinder.cpp parseInput) that reproduces all 8 rows of the
    // real binary's recorded output (CRISPR-OFF-001-Evidence.md); the expected rows below are its output.
    //  - DNA bulge: the gap offset on the padded side of the query builds the ungapped query itself (the
    //    tool reports it as type X, size 0): offset 0 when the PAM follows the target, guide length when
    //    it precedes it. It must not come back as a bulge.
    //  - RNA bulge: the removed run starts at 0..guideLength-size-1; the base next to the end away from
    //    the padding is never removable, so a "terminal" gap run (offset == windowLength) does not exist.
    // -------------------------------------------------------------------------------------------

    private static readonly (int Position, int TargetStart, BulgeType Bulge, int Size, int Gap, int Mismatches,
        string AlignedGuide, string AlignedTarget)[] PerfectCas9SiteRows =
    {
        (40, 21, BulgeType.Rna, 1, 0, 0, "GGCCGACCTGTCGCTGACGC", "-GCCGACCTGTCGCTGACGC"),
        (40, 21, BulgeType.Rna, 1, 1, 0, "GGCCGACCTGTCGCTGACGC", "G-CCGACCTGTCGCTGACGC"),
        (40, 22, BulgeType.Rna, 2, 0, 0, "GGCCGACCTGTCGCTGACGC", "--CCGACCTGTCGCTGACGC"),
    };

    /// <summary>
    /// A perfect forward site (guide + AGG) plus a second site with the PAM-proximal guide base missing
    /// (guide[0..19) + AGG). Cas-OFFinder 3.0.0b3 reports, at zero mismatches, only the three RNA-bulge rows
    /// that remove a leading guide base from the perfect site: no DNA bulge at offset 0 (the ungapped
    /// query) and no RNA bulge at the end of the guide (the second site is a deletion of the one base no
    /// version of the tool can remove).
    /// </summary>
    [Test]
    public void FindOffTargetsWithBulges_PerfectCas9Site_ReportsOnlyTheRowsCasOffinderReports()
    {
        const string genome =
            "ACGTACTGCA" + "TCGATCGATC" + Guide + "AGG" + "CATTCAGTCA" +
            "GGCCGACCTGTCGCTGACG" + "AGG" + "TACCATGCAT" + "GTCAGTCAGT";

        var hits = CrisprDesigner.FindOffTargetsWithBulges(Guide, new DnaSequence(genome), 0, 2, 2).ToArray();

        Assert.That(hits.Select(h => (h.Position, h.TargetStart, h.Bulge, h.BulgeSize, h.BulgePosition,
                h.Mismatches, h.AlignedGuide, h.AlignedTarget)),
            Is.EquivalentTo(PerfectCas9SiteRows));
        Assert.That(hits.All(h => h.IsForwardStrand), Is.True);
    }

    /// <summary>
    /// Cas12a (PAM before the target): the collapse sits at the OTHER end — a DNA gap run after the last
    /// guide base builds the ungapped query — and the unremovable base is the last one, away from the PAM.
    /// A perfect site and a site missing its last (PAM-distal) base therefore produce no 0-mismatch row.
    /// </summary>
    [Test]
    public void FindOffTargetsWithBulges_PerfectCas12aSite_ReportsNoRowAtZeroMismatches()
    {
        const string guide = "GTGGAATTGGTGTATTCAGCTTG";
        const string genome =
            "ACGTACTGCA" + "TCGATCGATC" + "TTTA" + guide + "CATTCAGTCA" +
            "TTTC" + "GTGGAATTGGTGTATTCAGCTT" + "GTAGTCAGTA" + "CCATGCATGA";

        Assert.That(CrisprDesigner.FindOffTargetsWithBulges(
            guide, new DnaSequence(genome), 0, 2, 2, CrisprSystemType.Cas12a), Is.Empty);
    }

    /// <summary>
    /// One DNA bulge with a substitution on the forward strand, one RNA bulge on the reverse strand
    /// (a T missing from the target). The four RNA rows of the reverse site are the equivalent gap
    /// placements around the missing base that stay within one mismatch.
    /// </summary>
    [Test]
    public void FindOffTargetsWithBulges_MixedBulgesOnBothStrands_MatchCasOffinderRows()
    {
        const string genome =
            "ACGTACTGCA" + "GGCCGACCATGTCGCAGACGCTGG" + "CATTCAGTCA" +
            "CCGGCGTCAGCGACGGTCGGCC" + "GTCAGTCAGT";

        var hits = CrisprDesigner.FindOffTargetsWithBulges(Guide, new DnaSequence(genome), 1, 1, 1).ToArray();

        Assert.That(hits.Select(h => (h.IsForwardStrand, h.Position, h.Bulge, h.BulgeSize, h.BulgePosition,
                h.Mismatches, h.AlignedGuide, h.AlignedTarget)),
            Is.EquivalentTo(new[]
            {
                (false, 44, BulgeType.Rna, 1, 6, 1, "GGCCGACCTGTCGCTGACGC", "GGCCGA-CCGTCGCTGACGC"),
                (false, 44, BulgeType.Rna, 1, 7, 1, "GGCCGACCTGTCGCTGACGC", "GGCCGAC-CGTCGCTGACGC"),
                (false, 44, BulgeType.Rna, 1, 8, 0, "GGCCGACCTGTCGCTGACGC", "GGCCGACC-GTCGCTGACGC"),
                (false, 44, BulgeType.Rna, 1, 9, 1, "GGCCGACCTGTCGCTGACGC", "GGCCGACCG-TCGCTGACGC"),
                (true, 31, BulgeType.Dna, 1, 8, 1, "GGCCGACC-TGTCGCTGACGC", "GGCCGACCATGTCGCAGACGC"),
            }));
    }

    /// <summary>
    /// Structural invariants of the admissible gap offsets, over a genome dense in near-copies of the
    /// guide, for a PAM-after and a PAM-before system: DNA offsets are 1..L (after) or 0..L-1 (before),
    /// RNA offsets 0..L-size-1, and the reported target geometry follows from the aligned target.
    /// </summary>
    [TestCase(CrisprSystemType.SpCas9, "AGG", true)]
    [TestCase(CrisprSystemType.Cas12a, "TTTA", false)]
    public void FindOffTargetsWithBulges_GapOffsets_StayWithinTheCasOffinderRange(
        CrisprSystemType system, string pam, bool pamAfter)
    {
        var rng = new Random(20261027);
        string Bases(int n) => new string(Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        int len = system == CrisprSystemType.Cas12a ? 23 : 20;
        string guide = Bases(len);
        var sb = new System.Text.StringBuilder(Bases(40));
        for (int i = 0; i < 30; i++)
        {
            string body = guide;
            if (i % 3 == 0) body = body.Remove(rng.Next(1, len - 2), 1);
            if (i % 3 == 1) body = body.Insert(rng.Next(1, len - 2), "A");
            body = body.Remove(rng.Next(len - 1), 1).Insert(rng.Next(len - 1), "C");
            sb.Append(pamAfter ? body + pam : pam + body).Append(Bases(20));
        }

        string genomeText = sb.ToString();
        var hits = CrisprDesigner.FindOffTargetsWithBulges(
            guide, new DnaSequence(genomeText), 4, 2, 2, system).ToArray();

        Assert.That(hits, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            foreach (var h in hits)
            {
                switch (h.Bulge)
                {
                    case BulgeType.Dna:
                        Assert.That(h.BulgePosition, Is.InRange(pamAfter ? 1 : 0, pamAfter ? len : len - 1), h.ToString());
                        break;
                    case BulgeType.Rna:
                        Assert.That(h.BulgePosition, Is.InRange(0, len - h.BulgeSize - 1), h.ToString());
                        break;
                    default:
                        Assert.That(h.BulgePosition, Is.EqualTo(-1), h.ToString());
                        break;
                }
                string forwardPam = genomeText.Substring(h.Position, pam.Length);
                Assert.That(h.Pam, Is.EqualTo(h.IsForwardStrand
                    ? forwardPam : DnaSequence.GetReverseComplementString(forwardPam)), h.ToString());
                int targetLength = h.AlignedTarget.Replace("-", "").Length;
                Assert.That(h.TargetStart, Is.EqualTo(h.IsForwardStrand
                    ? (pamAfter ? h.Position - targetLength : h.Position + pam.Length)
                    : (pamAfter ? h.Position + pam.Length : h.Position - targetLength)), h.ToString());
            }
        });
    }

    // -------------------------------------------------------------------------------------------
    // On-target interval and PAM on the hit (audit round 4, H26). The interval has the semantics of the
    // CRISPOR-faithful aggregates (crispor.py annotateOfftargets L1920-1929): an exact ungapped match is
    // dropped only when its protospacer+PAM footprint lies inside the interval; without one, all are dropped.
    // -------------------------------------------------------------------------------------------

    private const string PerfectSiteGenome =
        "ACGTACTGCA" + "TCGATCGATC" + Guide + "AGG" + "CATTCAGTCA" + "GTCAGTCAGT";   // footprint [20, 43)

    [Test]
    public void FindOffTargetsWithBulges_OnTargetInterval_KeepsTheExactUngappedRowOutsideIt()
    {
        var genome = new DnaSequence(PerfectSiteGenome);

        var inside = CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, 0, 2, 2, onTargetStart: 20, onTargetEnd: 43);
        var outside = CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, 0, 2, 2, onTargetStart: 0, onTargetEnd: 30);
        var none = CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, 0, 2, 2);

        Assert.Multiple(() =>
        {
            Assert.That(inside.Any(h => h.Bulge == BulgeType.None), Is.False);
            Assert.That(none.Any(h => h.Bulge == BulgeType.None), Is.False);
            var exact = outside.Where(h => h.Bulge == BulgeType.None).ToArray();
            Assert.That(exact, Has.Length.EqualTo(1));
            Assert.That(exact[0], Is.EqualTo(new BulgedOffTargetSite(
                Position: 40, TargetStart: 20, AlignedGuide: Guide, AlignedTarget: Guide, Mismatches: 0,
                Bulge: BulgeType.None, BulgeSize: 0, BulgePosition: -1, IsForwardStrand: true, Pam: "AGG")));
            Assert.That(outside.Count(h => h.Bulge != BulgeType.None), Is.EqualTo(3));
        });
    }

    [Test]
    public void FindOffTargetsWithBulges_OnTargetInterval_IsValidatedLikeTheAggregates()
    {
        var genome = new DnaSequence(PerfectSiteGenome);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => CrisprDesigner.FindOffTargetsWithBulges(
                Guide, genome, 0, 2, 2, onTargetStart: 0).ToList());
            Assert.Throws<ArgumentException>(() => CrisprDesigner.FindOffTargetsWithBulges(
                Guide, genome, 0, 2, 2, onTargetEnd: 40).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.FindOffTargetsWithBulges(
                Guide, genome, 0, 2, 2, onTargetStart: -1, onTargetEnd: 40).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.FindOffTargetsWithBulges(
                Guide, genome, 0, 2, 2, onTargetStart: 40, onTargetEnd: 40).ToList());
        });
    }

    [Test]
    public void FindOffTargets_ExistingSubstitutionOnlyScan_RemainsEmptyForBulgedToy()
    {
        var genome = new DnaSequence(Toy);
        Assert.That(CrisprDesigner.FindOffTargets(Guide, genome, 0), Is.Empty);
        Assert.That(CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, 0, 2, 2), Has.Exactly(4).Items);
    }

    [TestCase(-1, 0)]
    [TestCase(3, 0)]
    [TestCase(0, -1)]
    [TestCase(0, 3)]
    public void FindOffTargetsWithBulges_InvalidBulgeBounds_Throw(int dna, int rna)
    {
        string expected = dna is < 0 or > 2 ? "maxDnaBulge" : "maxRnaBulge";
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.FindOffTargetsWithBulges(Guide, new DnaSequence(Toy), 0, dna, rna).ToList());
        Assert.That(ex!.ParamName, Is.EqualTo(expected));
    }

    [Test]
    public void FindOffTargetsWithBulges_SharedGuardOrderAndTypes_ArePreserved()
    {
        var genome = new DnaSequence(Toy);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges(null!, null!, -1, -1, -1).ToList());
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges("", null!, -1, -1, -1).ToList());
            Assert.Throws<ArgumentNullException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges(Guide, null!, -1, -1, -1).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, -1, -1, -1).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges(Guide, genome, 6, -1, -1).ToList());
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges("ACGT", genome, 0, -1, -1).ToList());
            Assert.Throws<ArgumentException>(() =>
                CrisprDesigner.FindOffTargetsWithBulges("N" + Guide[1..], genome, 0, -1, -1).ToList());
        });
    }
}
