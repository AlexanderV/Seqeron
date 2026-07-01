using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class AlleleFrequenciesTests
{
    [Test]
    public void AlleleFrequencies_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.AlleleFrequencies(49, 42, 9));
        Assert.DoesNotThrow(() => PopulationTools.AlleleFrequencies(0, 0, 0));

        Assert.Throws<ArgumentException>(() => PopulationTools.AlleleFrequencies(-1, 0, 0));
        Assert.Throws<ArgumentException>(() => PopulationTools.AlleleFrequencies(0, -1, 0));
        Assert.Throws<ArgumentException>(() => PopulationTools.AlleleFrequencies(0, 0, -1));
    }

    [Test]
    public void AlleleFrequencies_Binding_InvokesSuccessfully()
    {
        // Wikipedia four-o'clock flower example: 49 AA, 42 Aa, 9 aa
        //   p = (2×49 + 42)/200 = 140/200 = 0.70
        //   q = (2×9  + 42)/200 =  60/200 = 0.30
        var flower = PopulationTools.AlleleFrequencies(49, 42, 9);
        Assert.Multiple(() =>
        {
            Assert.That(flower.MajorFreq, Is.EqualTo(0.70).Within(1e-10));
            Assert.That(flower.MinorFreq, Is.EqualTo(0.30).Within(1e-10));
            Assert.That(flower.MajorFreq + flower.MinorFreq, Is.EqualTo(1.0).Within(1e-10));
        });

        // Wikipedia diploid example: 6 AA, 3 AB, 1 BB → p = 15/20 = 0.75, q = 5/20 = 0.25
        var diploid = PopulationTools.AlleleFrequencies(6, 3, 1);
        Assert.Multiple(() =>
        {
            Assert.That(diploid.MajorFreq, Is.EqualTo(0.75).Within(1e-10));
            Assert.That(diploid.MinorFreq, Is.EqualTo(0.25).Within(1e-10));
        });
    }

    [Test]
    public void AlleleFrequencies_Binding_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            // All homozygous major → (1.0, 0.0)
            var allMajor = PopulationTools.AlleleFrequencies(100, 0, 0);
            Assert.That(allMajor.MajorFreq, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(allMajor.MinorFreq, Is.EqualTo(0.0).Within(1e-10));

            // All heterozygous → (0.5, 0.5)
            var allHet = PopulationTools.AlleleFrequencies(0, 100, 0);
            Assert.That(allHet.MajorFreq, Is.EqualTo(0.5).Within(1e-10));
            Assert.That(allHet.MinorFreq, Is.EqualTo(0.5).Within(1e-10));

            // Zero samples → (0, 0)
            var empty = PopulationTools.AlleleFrequencies(0, 0, 0);
            Assert.That(empty.MajorFreq, Is.EqualTo(0.0));
            Assert.That(empty.MinorFreq, Is.EqualTo(0.0));
        });
    }
}
