namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// MOTIF-REGULATORY-001 (review 2026-09, B05 follow-up): strand-annotated scan of the regulatory library.
/// Orientation-independent elements: CCAAT (Mantovani 1998, NAR 26:1135), GC box (Gidoni et al. 1985,
/// Science 230:511), enhancer sites AP-1 / NF-κB / E-box / CREB (Banerji et al. 1981, Cell 27:299).
/// Reference: Biopython 1.88 <c>Bio.SeqUtils.nt_search(seq, pattern)</c> on the sequence and on
/// <c>seq.reverse_complement()</c> (minus start = n − p − m).
/// </summary>
[TestFixture]
[Category("MOTIF-REGULATORY-001")]
public class MotifFinder_RegulatoryStrands_Tests
{
    private static IEnumerable<TestCaseData> BiopythonBothStrandCases()
    {
        yield return new TestCaseData("ATTGGTTTATAAACCGCCCATCCAATGGAAAGTCCCTGACTCAGGGCGGA", new[]
        {
            ("TATA Box", 7, '+', "TATAAA"), ("CAAT Box", 0, '-', "CCAAT"), ("CAAT Box", 21, '+', "CCAAT"),
            ("GC Box", 13, '-', "GGGCGG"), ("GC Box", 43, '+', "GGGCGG"), ("AP-1", 36, '+', "TGACTCA"),
            ("NF-κB", 26, '-', "GGGACTTTCC"),
        }).SetName("MixedOrientation_Case1");
        yield return new TestCaseData("GGGACTTTCCATTGGCCAATTTATTTAGGCACGTGCCGCCCGGGCGG", new[]
        {
            ("CAAT Box", 10, '-', "CCAAT"), ("CAAT Box", 15, '+', "CCAAT"), ("GC Box", 35, '-', "GGGCGG"),
            ("GC Box", 41, '+', "GGGCGG"), ("E-box", 29, '+', "CACGTG"), ("NF-κB", 0, '+', "GGGACTTTCC"),
        }).SetName("MixedOrientation_Case2");
    }

    [TestCaseSource(nameof(BiopythonBothStrandCases))]
    [Description("bothStrands=true equals Biopython nt_search on both strands for orientation-independent, non-self-RC elements")]
    public void FindRegulatoryElements_BothStrands_EqualsBiopythonNtSearch(string sequence, (string, int, char, string)[] expected)
    {
        var hits = MotifFinder.FindRegulatoryElements(new DnaSequence(sequence), bothStrands: true)
            .Select(e => (e.Name, e.Position, e.Strand, e.Sequence)).ToArray();
        Assert.That(hits, Is.EqualTo(expected));
    }

    [Test]
    [Description("bothStrands=false is exactly the legacy FindRegulatoryElements output with Strand '+'")]
    public void FindRegulatoryElements_SingleStrand_EqualsLegacyOverload()
    {
        var seq = new DnaSequence("ATTGGTTTATAAACCGCCCATCCAATGGAAAGTCCCTGACTCAGGGCGGACACGTGTGACGTCA");
        var legacy = MotifFinder.FindRegulatoryElements(seq)
            .Select(e => (e.Name, e.Position, e.Sequence, e.Pattern, e.Description));
        var stranded = MotifFinder.FindRegulatoryElements(seq, bothStrands: false).ToList();
        Assert.That(stranded.Select(e => (e.Name, e.Position, e.Sequence, e.Pattern, e.Description)), Is.EqualTo(legacy));
        Assert.That(stranded.All(e => e.Strand == '+'), Is.True);
    }

    [Test]
    [Description("Self-reverse-complementary AP-1/E-box/CREB are not duplicated on '-'; strand-specific TATA/polyA not scanned on '-'")]
    public void FindRegulatoryElements_BothStrands_NoPalindromeDuplicates_StrandSpecificPlusOnly()
    {
        // TGACTCA (AP-1), CACGTG (E-box), TGACGTCA (CREB) are palindromic; TTTATA = revcomp(TATAAA), TTTATT = revcomp(AATAAA).
        var hits = MotifFinder.FindRegulatoryElements(new DnaSequence("TGACTCAGGCACGTGGGTGACGTCAGGTTTATAGGTTTATT"), true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(hits.Where(h => h.Strand == '-'), Is.Empty);
            Assert.That(hits.Count(h => h.Name == "AP-1"), Is.EqualTo(1));
            Assert.That(hits.Count(h => h.Name == "CREB"), Is.EqualTo(1));
            Assert.That(hits.Any(h => h.Name is "TATA Box" or "Poly(A) Signal"), Is.False);
        });
    }

    [Test]
    public void OrientationIndependentRegulatoryElements_AreTheSourcedSet()
    {
        Assert.That(MotifFinder.OrientationIndependentRegulatoryElements,
            Is.EqualTo(new[] { "CAAT Box", "GC Box", "E-box", "AP-1", "NF-κB", "CREB" }));
    }

    [Test]
    public void FindRegulatoryElements_BothStrands_NullSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => MotifFinder.FindRegulatoryElements(null!, true));
    }
}
