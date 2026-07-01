using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.CalculateFunctionalDiversity: functional richness = distinct
// function count; functional diversity = Shannon entropy of function counts; per-pathway hit
// counts. Reference values hand-derived from the algorithm contract.
[TestFixture]
public class FunctionalDiversityTests
{
    private static MetagenomicsAnalyzer.FunctionalAnnotation Ann(string gene, string function, string pathway)
        => new(gene, function, pathway, KoNumber: "K0", CogCategory: "C", EValue: 1e-5, BitScore: 50);

    private static MetagenomicsAnalyzer.FunctionalAnnotation[] Annotations() => new[]
    {
        Ann("g1", "F1", "P1"),
        Ann("g2", "F1", "P1"),
        Ann("g3", "F2", "P2"),
    };

    [Test]
    public void FunctionalDiversity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.FunctionalDiversity(Annotations()));

        // Empty annotation set is defined (zero richness), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.FunctionalDiversity(
            Array.Empty<MetagenomicsAnalyzer.FunctionalAnnotation>()));
    }

    [Test]
    public void FunctionalDiversity_Binding_InvokesSuccessfully()
    {
        // 3 annotations: F1 x2, F2 x1 (richness 2). Shannon = -(2/3 ln 2/3 + 1/3 ln 1/3).
        // Pathways: P1 x2, P2 x1.
        var result = MetagenomicsTools.FunctionalDiversity(Annotations());

        double expectedShannon = -(2.0 / 3 * Math.Log(2.0 / 3) + 1.0 / 3 * Math.Log(1.0 / 3));

        Assert.Multiple(() =>
        {
            Assert.That(result.FunctionalRichness, Is.EqualTo(2.0),
                "Two distinct functions -> richness 2.");
            Assert.That(result.FunctionalDiversity, Is.EqualTo(expectedShannon).Within(1e-10),
                "Shannon entropy of function counts {2,1}.");
            var p1 = result.PathwayCounts.Single(p => p.Pathway == "P1");
            var p2 = result.PathwayCounts.Single(p => p.Pathway == "P2");
            Assert.That(p1.Count, Is.EqualTo(2));
            Assert.That(p2.Count, Is.EqualTo(1));
        });
    }
}
