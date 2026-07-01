using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CodingPotentialTests
{
    private static readonly Dictionary<string, double> Coding = new() { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
    private static readonly Dictionary<string, double> Noncoding = new() { ["ATGAAA"] = 2, ["AAACCC"] = 4 };

    [Test]
    public void CodingPotential_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CodingPotential("ATGAAACCC", Coding, Noncoding));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodingPotential("", Coding, Noncoding));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodingPotential(null!, Coding, Noncoding));
    }

    [Test]
    public void CodingPotential_TwoInFrameHexamers_ReturnsMeanLogRatio()
    {
        // Mirrors GenomeAnnotator_CalculateCodingPotential_Tests M1.
        // ATGAAA: ln(8/2)=1.3862943611...; AAACCC: ln(2/4)=-0.6931471805...; mean=0.34657359027997264.
        var result = AnnotationTools.CodingPotential("ATGAAACCC", Coding, Noncoding);
        Assert.That(result.Score, Is.EqualTo(0.34657359027997264).Within(1e-10));
    }

    [Test]
    public void CodingPotential_SingleHexamerBothTables_ReturnsNaturalLogRatio()
    {
        // Mirrors M2: ln(4/1) = ln 4 = 1.3862943611198906.
        var coding = new Dictionary<string, double> { ["ATGAAA"] = 4 };
        var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 1 };
        var result = AnnotationTools.CodingPotential("ATGAAA", coding, noncoding);
        Assert.That(result.Score, Is.EqualTo(1.3862943611198906).Within(1e-10));
    }

    [Test]
    public void CodingPotential_CodingAndNoncodingOnly_ContributePlusMinusOne()
    {
        // Mirrors M3/M4: coding-only -> +1, noncoding-only -> -1 for a single hexamer.
        var codingOnly = AnnotationTools.CodingPotential(
            "ATGAAA",
            new Dictionary<string, double> { ["ATGAAA"] = 5 },
            new Dictionary<string, double> { ["ATGAAA"] = 0 });
        var noncodingOnly = AnnotationTools.CodingPotential(
            "ATGAAA",
            new Dictionary<string, double> { ["ATGAAA"] = 0 },
            new Dictionary<string, double> { ["ATGAAA"] = 5 });

        Assert.Multiple(() =>
        {
            Assert.That(codingOnly.Score, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(noncodingOnly.Score, Is.EqualTo(-1.0).Within(1e-10));
        });
    }

    [Test]
    public void CodingPotential_SequenceShorterThanHexamer_ReturnsZero()
    {
        // Mirrors M7: len < wordSize -> 0.
        var result = AnnotationTools.CodingPotential(
            "ATGAA",
            new Dictionary<string, double> { ["ATGAA"] = 5 },
            new Dictionary<string, double> { ["ATGAA"] = 5 });
        Assert.That(result.Score, Is.EqualTo(0.0));
    }
}
