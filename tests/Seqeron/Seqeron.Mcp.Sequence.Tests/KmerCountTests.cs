using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class KmerCountTests
{
    [Test]
    public void KmerCount_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.KmerCount("ATGCGATCGATCG", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount("", 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount(null!, 3));
        Assert.Throws<ArgumentException>(() => SequenceTools.KmerCount("ATGC", 0)); // k < 1
    }

    [Test]
    public void KmerCount_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.KmerCount("ATGCATGC", 3);
        Assert.That(result.K, Is.EqualTo(3));
        Assert.That(result.Counts, Is.Not.Null);
        Assert.That(result.UniqueKmers, Is.GreaterThan(0));
        Assert.That(result.TotalKmers, Is.GreaterThan(0));

        // ATG appears twice
        Assert.That(result.Counts.ContainsKey("ATG"));
        Assert.That(result.Counts["ATG"], Is.EqualTo(2));
    }

    [Test]
    public void KmerCount_CanonicalAndAcgtOnly_MatchAnalysisCountKmersAndJellyfish()
    {
        // jellyfish count -m 3 [-C] + dump -c on ACGTNACGTAAcgtRTT: -C -> ACG 6, AAC 1, GTA 1, TAA 1;
        // without -C (ACGT-only) -> ACG 3, CGT 3, AAC 1, GTA 1, TAA 1.
        var canonical = SequenceTools.KmerCount("ACGTNACGTAAcgtRTT", 3, canonical: true);
        Assert.That(canonical.Counts, Is.EquivalentTo(new Dictionary<string, int> { ["ACG"] = 6, ["AAC"] = 1, ["GTA"] = 1, ["TAA"] = 1 }));
        Assert.That(canonical.UniqueKmers, Is.EqualTo(4));
        Assert.That(canonical.TotalKmers, Is.EqualTo(9));
        var acgt = SequenceTools.KmerCount("ACGTNACGTAAcgtRTT", 3, acgtOnly: true);
        Assert.That(acgt.Counts, Is.EquivalentTo(new Dictionary<string, int> { ["ACG"] = 3, ["CGT"] = 3, ["AAC"] = 1, ["GTA"] = 1, ["TAA"] = 1 }));
    }
}
