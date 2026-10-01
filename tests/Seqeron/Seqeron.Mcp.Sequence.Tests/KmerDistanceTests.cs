using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class KmerDistanceTests
{
    [Test]
    public void KmerDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerDistance("ATGCGATCG", "ATGCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("", "ATGC", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", "", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance(null!, "ATGC", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance("ATGC", "ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerDistance_Binding_InvokesSuccessfully()
    {
        // Identical sequences should have distance 0
        var identical = SequenceTools.KmerDistance("ATGCATGC", "ATGCATGC", 3);
        Assert.That(identical.Distance, Is.EqualTo(0).Within(0.0001));
        Assert.That(identical.K, Is.EqualTo(3));

        // Different sequences should have distance > 0
        var different = SequenceTools.KmerDistance("ATGCATGC", "CCCCCCCC", 3);
        Assert.That(different.Distance, Is.GreaterThan(0));
    }

    /// <summary>
    /// Audit round 1 WP4: parity with the Analysis server's kmer_distance — optional metric / markovOrder delegate to
    /// the same library overload (defaults unchanged).
    /// </summary>
    [Test]
    public void KmerDistance_Metric_ParityWithAnalysisServer()
    {
        const string s1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
        const string s2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.KmerDistance(s1, s2, 3).Distance,
                Is.EqualTo(Seqeron.Genomics.Analysis.KmerAnalyzer.KmerDistance(s1, s2, 3)));
            // Zielezinski et al. 2017 Fig. 1: Blaisdell d_E = 3, D2 = 5 (WP2 reference values).
            Assert.That(SequenceTools.KmerDistance("ATGTGTG", "CATGTG", 3, "squared_euclidean_counts").Distance, Is.EqualTo(3.0));
            Assert.That(SequenceTools.KmerDistance("ATGTGTG", "CATGTG", 3, "d2").Distance, Is.EqualTo(5.0));
            Assert.That(SequenceTools.KmerDistance(s1, s2, 3, "d2star").Distance, Is.EqualTo(0.44457941706964565).Within(1e-12));
            Assert.That(SequenceTools.KmerDistance(s1, s2, 3, "d2shepherd", 1).Distance, Is.EqualTo(0.5226793773737308).Within(1e-12));
        });
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance(s1, s2, 3, "bogus"));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerDistance(s1, s2, 3, "cosine", 2));
    }
}
