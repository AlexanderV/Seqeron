using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>most_frequent_kmers</c> MCP tool.
/// Expected values from the Frequent Words Problem (Rosalind BA1B sample) and hand
/// counts, NOT the wrapper's output.
/// </summary>
[TestFixture]
public class MostFrequentKmersTests
{
    [Test]
    public void MostFrequentKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MostFrequentKmers("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MostFrequentKmers("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MostFrequentKmers(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MostFrequentKmers("AAAA", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MostFrequentKmers("AAAA", -2));
    }

    [Test]
    public void MostFrequentKmers_Binding_InvokesSuccessfully()
    {
        // "ATGATG" k=3 -> ATG (count 2) is the unique winner.
        var single = AnalysisTools.MostFrequentKmers("ATGATG", 3).Kmers;
        Assert.That(single, Is.EquivalentTo(new[] { "ATG" }));

        // Rosalind BA1B sample: CATG and GCAT both occur 3 times.
        var tie = AnalysisTools.MostFrequentKmers("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4).Kmers;
        Assert.That(tie, Is.EquivalentTo(new[] { "CATG", "GCAT" }));
    }

    [Test]
    public void MostFrequentKmers_CanonicalAndAcgtOnly_EqualJellyfishDump()
    {
        // jellyfish 2.3.1 count -m 4 -C + dump -c on the BA1B sample: ATGC = 4 is the maximum (CATG 3);
        // without -C (ACGT-only) CATG = GCAT = 3.
        const string ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
        Assert.That(AnalysisTools.MostFrequentKmers(ba1b, 4, canonical: true).Kmers, Is.EquivalentTo(new[] { "ATGC" }));
        // AAAANTTTTGGGGuCCCC k=2: -C -> AA = CC = 6; ACGT-only -> AA = CC = GG = TT = 3.
        Assert.That(AnalysisTools.MostFrequentKmers("AAAANTTTTGGGGuCCCC", 2, canonical: true).Kmers, Is.EquivalentTo(new[] { "AA", "CC" }));
        Assert.That(AnalysisTools.MostFrequentKmers("AAAANTTTTGGGGuCCCC", 2, acgtOnly: true).Kmers, Is.EquivalentTo(new[] { "AA", "CC", "GG", "TT" }));
    }
}
