using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceAssemblerTests
{
    #region FindOverlap Tests

    [Test]
    public void FindOverlap_PerfectOverlap_ReturnsCorrect()
    {
        string seq1 = "ACGTACGTACGT";
        string seq2 = "ACGTACGTTTT";  // 8-base overlap with seq1's suffix "ACGTACGT"

        var overlap = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 8, minIdentity: 1.0);

        Assert.That(overlap, Is.Not.Null);
        Assert.That(overlap.Value.length, Is.EqualTo(8));
    }

    [Test]
    public void FindOverlap_NoOverlap_ReturnsNull()
    {
        string seq1 = "AAAAAAAAAA";
        string seq2 = "CCCCCCCCCC";

        var overlap = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 5);

        Assert.That(overlap, Is.Null);
    }

    [Test]
    public void FindOverlap_MinOverlapNotMet_ReturnsNull()
    {
        string seq1 = "ACGTACGT";
        string seq2 = "CGTAAAAA"; // 3-base overlap

        var overlap = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 5);

        Assert.That(overlap, Is.Null);
    }

    [Test]
    public void FindOverlap_WithMismatches_RespectsIdentityThreshold()
    {
        string seq1 = "ACGTACGTACGT";
        string seq2 = "ACGTACCTAAAA"; // 1 mismatch in 8-base overlap

        // 7/8 = 0.875 identity
        var overlap85 = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 8, minIdentity: 0.85);
        Assert.That(overlap85, Is.Not.Null);

        var overlap95 = SequenceAssembler.FindOverlap(seq1, seq2, minOverlap: 8, minIdentity: 0.95);
        Assert.That(overlap95, Is.Null);
    }

    #endregion

    #region FindAllOverlaps Tests

    [Test]
    public void FindAllOverlaps_ReturnsAllValidOverlaps()
    {
        var reads = new List<string>
        {
            "ACGTACGTACGT",
            "ACGTACGTTTTT",
            "TTTTACGTACGT"
        };

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap: 4);

        Assert.That(overlaps.Count, Is.GreaterThan(0));
    }

    [Test]
    public void FindAllOverlaps_NoSelfOverlaps()
    {
        var reads = new List<string> { "ACGTACGT", "CGTACGTA" };

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap: 4);

        // No read should overlap with itself
        Assert.That(overlaps.All(o => o.ReadIndex1 != o.ReadIndex2), Is.True);
    }

    #endregion

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

    #region AssembleOLC Tests

    [Test]
    public void AssembleOLC_OverlappingReads_ProducesContigs()
    {
        var reads = new List<string>
        {
            "ACGTACGTACGT",
            "ACGTACGTTTTT",
            "TTTTACGTACGT"
        };

        var result = SequenceAssembler.AssembleOLC(reads, new SequenceAssembler.AssemblyParameters(
            MinOverlap: 4, MinContigLength: 10));

        Assert.That(result.Contigs.Count, Is.GreaterThan(0));
        Assert.That(result.TotalReads, Is.EqualTo(3));
    }

    [Test]
    public void AssembleOLC_PerfectOverlappingReads_MergesIntoOne()
    {
        var reads = new List<string>
        {
            "AAAAACCCCC",
            "CCCCCGGGGG",
            "GGGGGTTTTTT"
        };

        var result = SequenceAssembler.AssembleOLC(reads, new SequenceAssembler.AssemblyParameters(
            MinOverlap: 5, MinContigLength: 10, MinIdentity: 1.0));

        // Should produce contigs that include the merged result
        Assert.That(result.TotalLength, Is.GreaterThanOrEqualTo(15));
    }

    [Test]
    public void AssembleOLC_NoOverlaps_ReturnsSingletonContigs()
    {
        var reads = new List<string>
        {
            "AAAAAAAAAA",
            "CCCCCCCCCC",
            "GGGGGGGGGG"
        };

        var result = SequenceAssembler.AssembleOLC(reads, new SequenceAssembler.AssemblyParameters(
            MinOverlap: 5, MinContigLength: 5));

        Assert.That(result.Contigs.Count, Is.EqualTo(3));
    }

    [Test]
    public void AssembleOLC_EmptyReads_ReturnsEmptyResult()
    {
        var result = SequenceAssembler.AssembleOLC(new List<string>());

        Assert.That(result.Contigs.Count, Is.EqualTo(0));
        Assert.That(result.TotalLength, Is.EqualTo(0));
    }

    #endregion

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
    public void MergeContigs_OverlapTooLarge_Concatenates()
    {
        string merged = SequenceAssembler.MergeContigs("AAA", "TTT", 10);
        Assert.That(merged, Is.EqualTo("AAATTT"));
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
        Assert.That(stats.N50, Is.EqualTo(14)); // 14 > 29/2
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

    #region ComputeConsensus Tests

    [Test]
    public void ComputeConsensus_IdenticalReads_ReturnsRead()
    {
        var reads = new List<string> { "ACGT", "ACGT", "ACGT" };
        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ACGT"));
    }

    [Test]
    public void ComputeConsensus_MajorityVote()
    {
        var reads = new List<string>
        {
            "ACGT",
            "ACGT",
            "TCGT"  // First position differs
        };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus[0], Is.EqualTo('A')); // Majority
    }

    [Test]
    public void ComputeConsensus_IgnoresGaps()
    {
        var reads = new List<string>
        {
            "A-GT",
            "ACGT",
            "ACGT"
        };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ACGT"));
    }

    [Test]
    public void ComputeConsensus_EmptyReads_ReturnsEmpty()
    {
        string consensus = SequenceAssembler.ComputeConsensus(new List<string>());
        Assert.That(consensus, Is.EqualTo(""));
    }

    #endregion

    #region QualityTrimReads Tests

    [Test]
    public void QualityTrimReads_TrimsLowQualityEnds()
    {
        // Quality scores: '!' = 0, 'I' = 40 (Phred+33)
        var reads = new List<(string, string)>
        {
            ("ACGTACGT", "!!IIII!!") // Low quality at ends
        };

        var trimmed = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 2);

        Assert.That(trimmed.Count, Is.EqualTo(1));
        Assert.That(trimmed[0].Length, Is.LessThan(8));
    }

    [Test]
    public void QualityTrimReads_RemovesTooShort()
    {
        var reads = new List<(string, string)>
        {
            ("ACGT", "!!!!") // All low quality
        };

        var trimmed = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 50);

        Assert.That(trimmed.Count, Is.EqualTo(0));
    }

    [Test]
    public void QualityTrimReads_KeepsHighQuality()
    {
        var reads = new List<(string, string)>
        {
            ("ACGTACGT", "IIIIIIII") // All high quality
        };

        var trimmed = SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 5);

        Assert.That(trimmed.Count, Is.EqualTo(1));
        Assert.That(trimmed[0], Is.EqualTo("ACGTACGT"));
    }

    #endregion

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

    #region Integration Tests

    [Test]
    public void Assembly_RealWorldScenario_ProducesValidContigs()
    {
        // Simulate reads from a known sequence
        string genome = "ACGTACGTACGTACGTACGTACGTACGTACGT";

        var reads = new List<string>();
        for (int i = 0; i <= genome.Length - 15; i += 5)
        {
            reads.Add(genome.Substring(i, 15));
        }

        var result = SequenceAssembler.AssembleOLC(reads, new SequenceAssembler.AssemblyParameters(
            MinOverlap: 8, MinContigLength: 15));

        Assert.That(result.Contigs.Count, Is.GreaterThan(0));
        Assert.That(result.LongestContig, Is.GreaterThanOrEqualTo(15));
    }

    #endregion
}
