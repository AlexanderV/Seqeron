using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CalculateTpmTests
{
    [Test]
    public void CalculateTpm_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CalculateTpm(
            new List<GeneCountInputDto> { new("G1", 10, 1000) }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.CalculateTpm(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CalculateTpm(new List<GeneCountInputDto>()));
    }

    [Test]
    public void CalculateTpm_TwoEqualLengthGenes_TpmProportionalToCounts()
    {
        // Rates: G1 = 10/1000 = 0.01, G2 = 30/1000 = 0.03; sumRates = 0.04.
        // TPM = rate/sumRates * 1e6 -> G1 = 250000, G2 = 750000 (sum = 1e6).
        // FPKM = count*1e9/(length*totalReads), totalReads = 40 -> G1 = 250000, G2 = 750000.
        var result = AnnotationTools.CalculateTpm(new List<GeneCountInputDto>
        {
            new("G1", 10, 1000),
            new("G2", 30, 1000)
        });

        Assert.That(result.Expressions, Has.Count.EqualTo(2));
        var g1 = result.Expressions.Single(e => e.GeneId == "G1");
        var g2 = result.Expressions.Single(e => e.GeneId == "G2");
        Assert.Multiple(() =>
        {
            Assert.That(g1.Tpm, Is.EqualTo(250_000).Within(1e-6));
            Assert.That(g2.Tpm, Is.EqualTo(750_000).Within(1e-6));
            Assert.That(g1.Tpm + g2.Tpm, Is.EqualTo(1_000_000).Within(1e-6));
            Assert.That(g1.Fpkm, Is.EqualTo(250_000).Within(1e-6));
            Assert.That(g2.Fpkm, Is.EqualTo(750_000).Within(1e-6));
            Assert.That(g1.RawCount, Is.EqualTo(10));
            Assert.That(g1.Length, Is.EqualTo(1000));
        });
    }

    [Test]
    public void CalculateTpm_ZeroCounts_AllTpmZero()
    {
        var result = AnnotationTools.CalculateTpm(new List<GeneCountInputDto>
        {
            new("G1", 0, 1000),
            new("G2", 0, 1000)
        });
        Assert.That(result.Expressions.All(e => e.Tpm == 0 && e.Fpkm == 0), Is.True);
    }
}
