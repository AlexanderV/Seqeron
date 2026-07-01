using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindConservedElementsTests
{
    private static IReadOnlyList<ConservationScoreDto> Run(int start, int count, double phastCons)
    {
        var list = new List<ConservationScoreDto>();
        for (int i = 0; i < count; i++)
            list.Add(new ConservationScoreDto("chr1", start + i, 0.0, phastCons, 0.0, 0));
        return list;
    }

    [Test]
    public void FindConservedElements_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindConservedElements(Run(100, 30, 0.9)));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindConservedElements(Array.Empty<ConservationScoreDto>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindConservedElements(null!));
    }

    [Test]
    public void FindConservedElements_ContiguousHighScores_ReturnsOneElement()
    {
        // 30 contiguous positions [100..129] with PhastCons 0.9 >= threshold 0.8, minLength 20.
        var result = AnnotationTools.FindConservedElements(Run(100, 30, 0.9), threshold: 0.8, minLength: 20);

        Assert.That(result.Elements, Has.Count.EqualTo(1));
        var e = result.Elements[0];
        Assert.Multiple(() =>
        {
            Assert.That(e.Chromosome, Is.EqualTo("chr1"));
            Assert.That(e.Start, Is.EqualTo(100));
            Assert.That(e.End, Is.EqualTo(129));
            Assert.That(e.Score, Is.EqualTo(0.9).Within(1e-10));
        });
    }

    [Test]
    public void FindConservedElements_BelowThreshold_ReturnsNone()
    {
        // All PhastCons 0.5 < threshold 0.8 -> no conserved element.
        var result = AnnotationTools.FindConservedElements(Run(100, 30, 0.5), threshold: 0.8, minLength: 20);
        Assert.That(result.Elements, Is.Empty);
    }

    [Test]
    public void FindConservedElements_TooShort_ReturnsNone()
    {
        // Only 10 conserved positions but minLength 20 -> excluded.
        var result = AnnotationTools.FindConservedElements(Run(100, 10, 0.9), threshold: 0.8, minLength: 20);
        Assert.That(result.Elements, Is.Empty);
    }
}
