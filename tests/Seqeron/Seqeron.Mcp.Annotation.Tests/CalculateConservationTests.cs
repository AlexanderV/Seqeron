using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CalculateConservationTests
{
    private static ConservationPositionInputDto Pos(string alleles) => new("chr1", 100, alleles);

    [Test]
    public void CalculateConservation_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CalculateConservation(
            new List<ConservationPositionInputDto> { Pos("AAAA") }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CalculateConservation(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CalculateConservation(
            new List<ConservationPositionInputDto>()));
    }

    [Test]
    public void CalculateConservation_FullyConserved_MaxScores()
    {
        // 10 identical alleles: fraction = 1.0 -> PhyloP = (1-0.5)*12 = 6, PhastCons = 1,
        // GERP = clamp((10 - 2.5)/max(1,7.5) * 6) = 6, conserved = 10.
        var result = AnnotationTools.CalculateConservation(
            new List<ConservationPositionInputDto> { Pos("AAAAAAAAAA") });

        Assert.That(result.Scores, Has.Count.EqualTo(1));
        var s = result.Scores[0];
        Assert.Multiple(() =>
        {
            Assert.That(s.PhyloP, Is.EqualTo(6.0).Within(1e-9));
            Assert.That(s.PhastCons, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(s.Gerp, Is.EqualTo(6.0).Within(1e-9));
            Assert.That(s.ConservedSpeciesCount, Is.EqualTo(10));
        });
    }

    [Test]
    public void CalculateConservation_HalfConserved_MidScores()
    {
        // "AATT": ref = A, conserved = 2/4 = 0.5 -> PhyloP = 0, PhastCons = 0.5,
        // GERP = clamp((2 - 1)/max(1,3) * 6) = 2, conserved = 2.
        var s = AnnotationTools.CalculateConservation(
            new List<ConservationPositionInputDto> { Pos("AATT") }).Scores[0];
        Assert.Multiple(() =>
        {
            Assert.That(s.PhyloP, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(s.PhastCons, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(s.Gerp, Is.EqualTo(2.0).Within(1e-9));
            Assert.That(s.ConservedSpeciesCount, Is.EqualTo(2));
        });
    }

    [Test]
    public void CalculateConservation_EmptyAlleles_ZeroScores()
    {
        var s = AnnotationTools.CalculateConservation(
            new List<ConservationPositionInputDto> { Pos("") }).Scores[0];
        Assert.Multiple(() =>
        {
            Assert.That(s.PhyloP, Is.EqualTo(0.0));
            Assert.That(s.PhastCons, Is.EqualTo(0.0));
            Assert.That(s.Gerp, Is.EqualTo(0.0));
            Assert.That(s.ConservedSpeciesCount, Is.EqualTo(0));
        });
    }
}
