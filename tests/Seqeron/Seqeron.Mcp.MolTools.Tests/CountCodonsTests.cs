using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CountCodonsTests
{
    [Test]
    public void CountCodons_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.count_codons("ATGATGTTT"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.count_codons(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.count_codons(null!));
    }

    [Test]
    public void CountCodons_Binding_InvokesSuccessfully()
    {
        // "ATGATGTTT" -> ATG,ATG,TTT.
        var counts = MolToolsTools.count_codons("ATGATGTTT").Counts;
        Assert.Multiple(() =>
        {
            Assert.That(counts["ATG"], Is.EqualTo(2));
            Assert.That(counts["TTT"], Is.EqualTo(1));
            Assert.That(counts.Count, Is.EqualTo(2));
        });

        // Non-ACGT codon (ANG) and trailing partial (last 'C') are skipped.
        // "ATGANGTTTC" -> ATG(valid), ANG(invalid), TTT(valid), trailing "C" ignored.
        var filtered = MolToolsTools.count_codons("ATGANGTTTC").Counts;
        Assert.Multiple(() =>
        {
            Assert.That(filtered["ATG"], Is.EqualTo(1));
            Assert.That(filtered["TTT"], Is.EqualTo(1));
            Assert.That(filtered.ContainsKey("ANG"), Is.False);
            Assert.That(filtered.Count, Is.EqualTo(2));
        });
    }
}
