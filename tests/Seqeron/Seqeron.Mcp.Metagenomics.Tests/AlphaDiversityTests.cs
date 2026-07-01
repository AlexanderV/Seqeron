using System;
using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.CalculateAlphaDiversity.
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_AlphaDiversity_Tests
// (Shannon 1948; Simpson 1949; Hill 1973; Chao 1984).
[TestFixture]
public class AlphaDiversityTests
{
    private const double Ln2 = 0.6931471805599453;

    private static AbundanceItem[] Ab(params (string Name, double Fraction)[] items)
        => Array.ConvertAll(items, t => new AbundanceItem(t.Name, t.Fraction));

    [Test]
    public void AlphaDiversity_Schema_ValidatesCorrectly()
    {
        // Valid vector — no throw.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.AlphaDiversity(Ab(("s1", 0.5), ("s2", 0.5))));

        // Empty vector is defined input (all-zero metrics), not an error.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.AlphaDiversity(Array.Empty<AbundanceItem>()));
    }

    [Test]
    public void AlphaDiversity_Binding_InvokesSuccessfully()
    {
        // Two equal species (0.5, 0.5): Shannon = ln2, Simpson = 0.5, InvSimpson = 2,
        // Pielou = 1, ObservedSpecies = 2, Chao1 = 2 (proportional data -> S_obs).
        var even = MetagenomicsTools.AlphaDiversity(Ab(("s1", 0.5), ("s2", 0.5)));

        Assert.Multiple(() =>
        {
            Assert.That(even.ObservedSpecies, Is.EqualTo(2));
            Assert.That(even.ShannonIndex, Is.EqualTo(Ln2).Within(1e-10));
            Assert.That(even.SimpsonIndex, Is.EqualTo(0.5).Within(1e-10));
            Assert.That(even.InverseSimpson, Is.EqualTo(2.0).Within(1e-10));
            Assert.That(even.PielouEvenness, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(even.Chao1Estimate, Is.EqualTo(2.0).Within(1e-10));
        });

        // Chao1 with count data {50,30,1,1,2}: S_obs=5, f1=2, f2=1 -> 5 + 2^2/(2*1) = 7.
        var chao = MetagenomicsTools.AlphaDiversity(
            Ab(("d1", 50), ("d2", 30), ("s1", 1), ("s2", 1), ("db", 2)));
        Assert.Multiple(() =>
        {
            Assert.That(chao.ObservedSpecies, Is.EqualTo(5));
            Assert.That(chao.Chao1Estimate, Is.EqualTo(7.0).Within(1e-10),
                "Chao1 = 5 + 2^2/(2*1) = 7 (Chao 1984).");
        });
    }
}
