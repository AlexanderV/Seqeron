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
            BulgePosition: 8, IsForwardStrand: true)));
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
