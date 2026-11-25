namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssemblerTests
{
    // NOTE: FindOverlap and FindAllOverlaps tests were moved to the canonical
    // ASSEMBLY-OLC-001 fixture: SequenceAssembler_AssembleOLC_Tests.cs.

    #region CalculateIdentity Tests

    [Test]
    public void CalculateIdentity_IdenticalSequences_ReturnsOne()
    {
        double identity = SequenceAssembler.CalculateIdentity("ACGT", "ACGT");
        Assert.That(identity, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateIdentity_CompletelyDifferent_ReturnsZero()
    {
        double identity = SequenceAssembler.CalculateIdentity("AAAA", "TTTT");
        Assert.That(identity, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateIdentity_HalfMatch_ReturnsFifty()
    {
        double identity = SequenceAssembler.CalculateIdentity("AATT", "AAGG");
        Assert.That(identity, Is.EqualTo(0.5));
    }

    [Test]
    public void CalculateIdentity_CaseInsensitive()
    {
        double identity = SequenceAssembler.CalculateIdentity("acgt", "ACGT");
        Assert.That(identity, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateIdentity_EmptySequences_ReturnsOne()
    {
        double identity = SequenceAssembler.CalculateIdentity("", "");
        Assert.That(identity, Is.EqualTo(1.0));
    }

    [Test]
    public void CalculateIdentity_DifferentLengths_ReturnsZero()
    {
        double identity = SequenceAssembler.CalculateIdentity("ACGT", "ACG");
        Assert.That(identity, Is.EqualTo(0.0));
    }

    #endregion

    // NOTE: AssembleOLC tests were moved to the canonical ASSEMBLY-OLC-001 fixture:
    // SequenceAssembler_AssembleOLC_Tests.cs.

    #region AssembleDeBruijn Tests

    [Test]
    public void AssembleDeBruijn_SimpleReads_ProducesContigs()
    {
        var reads = new List<string>
        {
            "ACGTACGTACGT",
            "ACGTACGTTTTT"
        };

        var result = SequenceAssembler.AssembleDeBruijn(reads, new SequenceAssembler.AssemblyParameters(
            KmerSize: 5, MinContigLength: 5));

        Assert.That(result.Contigs.Count, Is.GreaterThan(0));
    }

    [Test]
    public void AssembleDeBruijn_EmptyReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleDeBruijn(new List<string>());

        Assert.That(result.Contigs.Count, Is.EqualTo(0));
    }

    #endregion

    #region MergeContigs Tests

    [Test]
    public void MergeContigs_ValidOverlap_MergesCorrectly()
    {
        string contig1 = "ACGTACGT";
        string contig2 = "ACGTTTTT";

        string merged = SequenceAssembler.MergeContigs(contig1, contig2, 4);

        Assert.That(merged, Is.EqualTo("ACGTACGTTTTT"));
    }

    [Test]
    public void MergeContigs_ZeroOverlap_Concatenates()
    {
        string merged = SequenceAssembler.MergeContigs("AAAA", "TTTT", 0);
        Assert.That(merged, Is.EqualTo("AAAATTTT"));
    }

    [Test]
    public void MergeContigs_OverlapTooLarge_Throws()
    {
        Assert.That(() => SequenceAssembler.MergeContigs("AAA", "TTT", 10),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
    }

    #endregion

    #region CalculateStats Tests

    [Test]
    public void CalculateStats_ValidContigs_CalculatesN50()
    {
        var contigs = new List<string>
        {
            "AAAAAAAAAA",     // 10
            "CCCCCCCCCCCCCC", // 14
            "GGGGG"           // 5
        };

        var stats = SequenceAssembler.CalculateStats(contigs, 10);

        Assert.That(stats.TotalLength, Is.EqualTo(29));
        Assert.That(stats.LongestContig, Is.EqualTo(14));
        // QUAST N50([14,10,5]) = 10: 14 covers 14/29 < 50 %, 14+10 = 24 >= 14.5 (was 14 via 29/2 = 14).
        Assert.That(stats.N50, Is.EqualTo(10));
    }

    // ASSEMBLY-STATS-001 F: N50 used integer halving (total/2), so an odd total stopped one contig early.
    // Expected values: QUAST quast_libs/N50.py N50() on the same lengths.
    [TestCase(new[] { 4, 3, 2 }, 3.0)]
    [TestCase(new[] { 2, 1, 1, 1 }, 1.0)]
    [TestCase(new[] { 10, 20, 30 }, 30.0)]
    [TestCase(new[] { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 }, 7.0)]
    public void CalculateStats_N50_MatchesQuast(int[] lengths, double expectedN50)
    {
        var contigs = lengths.Select(l => new string('A', l)).ToList();

        var stats = SequenceAssembler.CalculateStats(contigs, 0);

        Assert.That(stats.N50, Is.EqualTo(expectedN50));
    }

    [Test]
    public void CalculateStats_N50_EqualsCanonicalGenomeAssemblyAnalyzerN50()
    {
        var lengths = new[] { 17, 3, 99, 42, 42, 8, 1, 65 };
        var stats = SequenceAssembler.CalculateStats(lengths.Select(l => new string('C', l)).ToList(), 0);

        Assert.That(stats.N50, Is.EqualTo(GenomeAssemblyAnalyzer.CalculateN50(lengths)),
            "One Nx implementation: both statistics paths agree.");
    }

    [Test]
    public void CalculateStats_AssembledReads_DefaultsToTotalAndAcceptsExplicitCount()
    {
        var contigs = new List<string> { "ACGTACGT" };

        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.CalculateStats(contigs, 7).AssembledReads, Is.EqualTo(7));
            Assert.That(SequenceAssembler.CalculateStats(contigs, 7, 4).AssembledReads, Is.EqualTo(4));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAssembler.CalculateStats(contigs, 7, 8));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAssembler.CalculateStats(contigs, 7, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAssembler.CalculateStats(contigs, -1));
            Assert.Throws<ArgumentNullException>(() => SequenceAssembler.CalculateStats(null!, 1));
        });
    }

    // Velvet usedReads (graphStats.c): a read is used when it lies on a reported contig (node of at
    // least the minimum contig length). r1+r2 overlap by 6 -> contig of 14 (kept); r3 alone (10) is
    // filtered by MinContigLength 12 -> 2 of 3 reads assembled.
    [Test]
    public void AssembleOLC_AssembledReads_CountsOnlyReadsInReportedContigs()
    {
        var reads = new List<string> { "ACGTTGCAAG", "TGCAAGCTTA", "GGGGGGGGGG" };
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 5, MinIdentity: 1.0, KmerSize: 31, MinContigLength: 12);

        var result = SequenceAssembler.AssembleOLC(reads, param);

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTGCAAGCTTA" }));
            Assert.That(result.TotalReads, Is.EqualTo(3));
            Assert.That(result.AssembledReads, Is.EqualTo(2));
        });
    }

    // Velvet usedReads for the de Bruijn graph: reads with a k-mer in a reported contig. k=4:
    // ACGTTGCA spells the 8-bp unitig (kept); CGTTG lies inside it; TTTTT is an isolated TTT cycle
    // spelled TTTT (4 < MinContigLength 5, dropped); AC is shorter than k -> 2 of 4 reads used.
    [Test]
    public void AssembleDeBruijn_AssembledReads_CountsReadsWithAKmerInAReportedContig()
    {
        var reads = new List<string> { "ACGTTGCA", "CGTTG", "TTTTT", "AC" };
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 20, MinIdentity: 0.9, KmerSize: 4, MinContigLength: 5);

        var result = SequenceAssembler.AssembleDeBruijn(reads, param);

        Assert.Multiple(() =>
        {
            Assert.That(result.Contigs, Is.EqualTo(new[] { "ACGTTGCA" }));
            Assert.That(result.TotalReads, Is.EqualTo(4));
            Assert.That(result.AssembledReads, Is.EqualTo(2));
        });
    }

    // Double-stranded DBG: a read given on the reverse strand is used through the twin contig.
    [Test]
    public void AssembleDeBruijn_ReverseComplementOption_CountsOppositeStrandReadAsAssembled()
    {
        // TGCAACGT is the reverse complement of ACGTTGCA.
        var reads = new List<string> { "ACGTTGCA", "TGCAACGT" };
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 20, MinIdentity: 0.9, KmerSize: 5, MinContigLength: 1);

        var result = SequenceAssembler.AssembleDeBruijn(reads, param, includeReverseComplement: true);

        Assert.That(result.AssembledReads, Is.EqualTo(2));
    }

    [Test]
    public void CalculateStats_EmptyContigs_ReturnsZeros()
    {
        var stats = SequenceAssembler.CalculateStats(new List<string>(), 10);

        Assert.That(stats.TotalLength, Is.EqualTo(0));
        Assert.That(stats.N50, Is.EqualTo(0));
    }

    #endregion

    #region Scaffold Tests

    [Test]
    public void Scaffold_WithLinks_JoinsContigs()
    {
        var contigs = new List<string>
        {
            "AAAAAAA",
            "CCCCCCC",
            "GGGGGGG"
        };

        var links = new List<(int, int, int)>
        {
            (0, 1, 5),
            (1, 2, 10)
        };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.That(scaffolds.Count, Is.EqualTo(1));
        Assert.That(scaffolds[0], Does.Contain("NNNNN")); // Gap
    }

    [Test]
    public void Scaffold_NoLinks_ReturnsOriginalContigs()
    {
        var contigs = new List<string> { "AAAA", "CCCC" };
        var scaffolds = SequenceAssembler.Scaffold(contigs, new List<(int, int, int)>());

        Assert.That(scaffolds.Count, Is.EqualTo(2));
    }

    [Test]
    public void Scaffold_CustomGapCharacter()
    {
        var contigs = new List<string> { "AAA", "CCC" };
        var links = new List<(int, int, int)> { (0, 1, 3) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links, gapCharacter: 'X');

        Assert.That(scaffolds[0], Does.Contain("XXX"));
    }

    #endregion

    #region CalculateCoverage Tests

    [Test]
    public void CalculateCoverage_ReadsMapToReference()
    {
        string reference = "ACGTACGTACGT";
        var reads = new List<string> { "ACGTACGT", "CGTACGTA" };

        var coverage = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(coverage.Length, Is.EqualTo(reference.Length));
        Assert.That(coverage.Any(c => c > 0), Is.True);
    }

    [Test]
    public void CalculateCoverage_NoMatchingReads_ZeroCoverage()
    {
        string reference = "AAAAAAAAAA";
        var reads = new List<string> { "CCCCCCCC" };

        var coverage = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(coverage.All(c => c == 0), Is.True);
    }

    #endregion

    // ComputeConsensus tests consolidated into SequenceAssembler_ComputeConsensus_Tests.cs
    // (canonical file for ASSEMBLY-CONSENSUS-001).

    // QualityTrimReads tests consolidated into SequenceAssembler_QualityTrimReads_Tests.cs
    // (canonical file for ASSEMBLY-TRIM-001).

    #region ErrorCorrectReads Tests

    [Test]
    public void ErrorCorrectReads_CorrectsSingleErrors()
    {
        var reads = new List<string>
        {
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTNCGTACGTACGT"  // Error: N instead of A at position 8
        };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 7, minKmerFrequency: 3);

        // The corrected reads should be returned
        Assert.That(corrected.Count, Is.EqualTo(5));
    }

    [Test]
    public void ErrorCorrectReads_NoChangesNeeded()
    {
        var reads = new List<string>
        {
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT"
        };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 5);

        Assert.That(corrected, Is.EqualTo(reads));
    }

    #endregion

}
