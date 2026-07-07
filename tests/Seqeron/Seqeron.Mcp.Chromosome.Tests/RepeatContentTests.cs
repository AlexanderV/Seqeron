using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;
using SourceGa = Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>repeat_content</c>. GenomeAssemblyAnalyzer.CalculateRepeatContent sums inclusive
/// (End - Start + 1) lengths, computes total*100/genomeLength, and groups lengths by RepeatClass.
/// </summary>
[TestFixture]
public class RepeatContentTests
{
    private static SourceGa.RepeatAnnotation R(int start, int end, string cls) =>
        new("s", start, end, cls, "fam", 0.0, '+');

    [Test]
    public void RepeatContent_Schema_ValidatesCorrectly()
    {
        var repeats = new List<SourceGa.RepeatAnnotation> { R(0, 9, "LINE") };
        Assert.DoesNotThrow(() => ChromosomeTools.RepeatContent(repeats, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.RepeatContent(null!, 100));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.RepeatContent(repeats, 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.RepeatContent(repeats, -1));
    }

    [Test]
    public void RepeatContent_Binding_InvokesSuccessfully()
    {
        // LINE 0..99 (100 bp), SINE 100..149 (50 bp); genome 1000 -> 150 total, 15%.
        var repeats = new List<SourceGa.RepeatAnnotation> { R(0, 99, "LINE"), R(100, 149, "SINE") };
        var rc = ChromosomeTools.RepeatContent(repeats, 1000);

        Assert.Multiple(() =>
        {
            Assert.That(rc.TotalRepeatLength, Is.EqualTo(150));
            Assert.That(rc.RepeatPercentage, Is.EqualTo(15.0).Within(1e-9));
            Assert.That(rc.RepeatClassLengths["LINE"], Is.EqualTo(100));
            Assert.That(rc.RepeatClassLengths["SINE"], Is.EqualTo(50));
        });
    }
}
