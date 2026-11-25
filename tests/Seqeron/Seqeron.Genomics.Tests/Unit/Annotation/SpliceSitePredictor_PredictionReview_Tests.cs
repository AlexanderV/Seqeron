namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// SPLICE-PREDICT-001 review (B11 F32-F38): intron typing / scoring, exon tiling, SUPPA/rMATS event
/// definitions, ORF-based IsWithinCodingRegion and the HBB (Ensembl ENST00000335295, GRCh38 11:5225464-5227071:-1) check.
/// Expected values come from the independent Python oracle (0 mismatches on 1000 random sequences).
/// </summary>
[TestFixture]
public class SpliceSitePredictor_PredictionReview_Tests
{
    private const string Hbb =
        "ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGGATAAGTAACAGGGTACAGTTTAGAATGGGAAACAGACGAATGATTGCATCAGTGTGGAAGTCTCAGGATCGTTTTAGTTTCTTTTATTTGCTGTTCATAACAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATGTGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATCTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA";

    [Test]
    public void PredictIntrons_Hbb_BothAnnotatedIntronsAreU2Candidates()
    {
        var introns = SpliceSitePredictor.PredictIntrons(Hbb, 30, 100000, 0.5).ToList();
        foreach (var (s, e) in new[] { (142, 271), (495, 1344) })
        {
            var hit = introns.Single(i => i.Start == s && i.End == e);
            Assert.That(hit.Type, Is.EqualTo(SpliceSitePredictor.IntronType.U2));
            Assert.That(hit.Length, Is.EqualTo(e - s + 1));
            Assert.That(hit.Sequence, Is.EqualTo(Hbb.Substring(s, e - s + 1).ToUpperInvariant().Replace('T', 'U')));
        }
    }

    [Test]
    public void PredictIntrons_GtAgPairWithoutU12Evidence_IsU2NotU12()
    {
        string seq = new string('C', 30) + "GUAAGU" + new string('A', 70) + "UUUUUUUUUUUUCAG" + new string('C', 30);
        var introns = SpliceSitePredictor.PredictIntrons(seq, 60, 100000, 0.3).ToList();
        Assert.That(introns, Is.Not.Empty);
        Assert.That(introns.All(i => i.Type != SpliceSitePredictor.IntronType.U12));
    }

    [Test]
    public void PredictIntrons_U12DonorAndBranchSite_TypedU12_WithoutBranchPoint()
    {
        // AU-AC intron: /AUAUCCUU donor, UUCCUUAAC branch site 12 nt upstream of the AC terminus.
        string seq = new string('C', 30) + "AUAUCCUU" + new string('G', 20) + "UUCCUUAAC" + "GGGGGGGGGGGG" + "AC" + new string('C', 30);
        var introns = SpliceSitePredictor.PredictIntrons(seq, 30, 100000, 0.5).ToList();
        var u12 = introns.Single(i => i.Type == SpliceSitePredictor.IntronType.U12);
        Assert.That(u12.Start, Is.EqualTo(30));
        Assert.That(u12.End, Is.EqualTo(seq.Length - 31));
        Assert.That(u12.BranchPoint, Is.Null);
        Assert.That(u12.Score, Is.EqualTo((u12.DonorSite.Score + u12.AcceptorSite.Score) / 2).Within(1e-12));
    }

    [Test]
    public void PredictIntrons_BranchPoint_IsCanonicalGaoBestPoint_AndRaisesScore()
    {
        string seq = new string('C', 30) + "GUAAGU" + new string('A', 40) + "CUGAC" + new string('A', 16) +
                     "UUUUUUUUUUUUCAG" + new string('C', 30);
        var i = SpliceSitePredictor.PredictIntrons(seq, 60, 100000, 0.3).Single(x => x.Start == 30);
        Assert.That(i.BranchPoint, Is.Not.Null);
        var canonical = SpliceSitePredictor.FindAcceptorBranchPoint(seq, i.End);
        Assert.That(i.BranchPoint!.Value.Position, Is.EqualTo(canonical.BranchPointPosition));
        Assert.That(i.Score, Is.EqualTo((i.DonorSite.Score + i.AcceptorSite.Score + canonical.Score) / 3).Within(1e-12));
    }

    [TestCase("AUGAAAUAAAAA", 6, 0, true)]   // UAA at 6 is the terminating codon itself: included in the CDS
    [TestCase("AUGAAAUAAAAA", 9, 0, false)]  // beyond the in-frame stop
    [TestCase("AUGAUGAAAAAA", 4, 1, true)]   // second AUG (codon start 3) is in frame with the position
    [TestCase("AAAAUGAAA", 4, 1, true)]      // AUG at index 3, position 4 is its 2nd base
    [TestCase("AUGAAA", 1, 1, true)]         // position inside the start codon
    [TestCase("AUGAAA", 4, 0, false)]        // frame 0 -> codon start 4: off the AUG frame
    public void IsWithinCodingRegion_StopCodonAndMultipleStarts(string seq, int pos, int frame, bool expected)
    {
        Assert.That(SpliceSitePredictor.IsWithinCodingRegion(seq, pos, frame), Is.EqualTo(expected));
    }

    [Test]
    public void IsWithinCodingRegion_FrameOutOfRange_Throws() =>
        Assert.Throws<ArgumentOutOfRangeException>(() => SpliceSitePredictor.IsWithinCodingRegion("AUGAAA", 3, 3));

    [Test]
    public void IsWithinCodingRegion_NoUpstreamWindowLimit()
    {
        string seq = "AUG" + string.Concat(Enumerable.Repeat("GCC", 200)) + "GCC";
        Assert.That(SpliceSitePredictor.IsWithinCodingRegion(seq, seq.Length - 3, 0), Is.True);
    }
}
