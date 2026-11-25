using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DinucleotideShuffleTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    private static Dictionary<string, int> Doublets(string s)
    {
        var d = new Dictionary<string, int>();
        for (int i = 0; i < s.Length - 1; i++)
            d[s.Substring(i, 2)] = d.GetValueOrDefault(s.Substring(i, 2)) + 1;
        return d;
    }

    [Test]
    public void DinucleotideShuffle_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.DinucleotideShuffle(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.DinucleotideShuffle(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.DinucleotideShuffle(Mir21, 1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.DinucleotideShuffle(Mir21, 1, 1001));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.DinucleotideShuffle(new string('A', 10001)));
    }

    [Test]
    public void DinucleotideShuffle_ManySeeds_PreservesDinucleotideCountsFirstLastLength()
    {
        // Altschul & Erickson 1985 invariant (library CL6), over many seeds and shuffles.
        var orig = Doublets(Mir21);
        for (int seed = 0; seed < 50; seed++)
        {
            var r = AnnotationTools.DinucleotideShuffle(Mir21, seed, 5);
            Assert.That(r.Shuffles, Has.Count.EqualTo(5));
            foreach (var s in r.Shuffles)
            {
                Assert.That(s.Length, Is.EqualTo(Mir21.Length));
                Assert.That(s[0], Is.EqualTo(Mir21[0]));
                Assert.That(s[^1], Is.EqualTo(Mir21[^1]));
                Assert.That(Doublets(s), Is.EquivalentTo(orig));
            }
        }
    }

    [Test]
    public void DinucleotideShuffle_Binding_EqualsLibrarySequentialDraws()
    {
        var r = AnnotationTools.DinucleotideShuffle(Mir21, 7, 3);
        var rng = new Random(7);
        for (int i = 0; i < 3; i++)
            Assert.That(r.Shuffles[i], Is.EqualTo(MiRnaAnalyzer.DinucleotideShuffle(Mir21, rng)));
        Assert.That(r.Seed, Is.EqualTo(7));
        Assert.That(r.Count, Is.EqualTo(3));
        Assert.That(AnnotationTools.DinucleotideShuffle(Mir21, 7, 3).Shuffles, Is.EqualTo(r.Shuffles)); // deterministic
        Assert.That(AnnotationTools.DinucleotideShuffle(Mir21, 8, 3).Shuffles, Is.Not.EqualTo(r.Shuffles));
    }

    [Test]
    public void DinucleotideShuffle_DnaInputNormalisedAndShortReturnedUnchanged()
    {
        Assert.That(AnnotationTools.DinucleotideShuffle("acgt", 1).Shuffles[0], Has.Length.EqualTo(4).And.Matches("^[ACGU]+$"));
        Assert.That(AnnotationTools.DinucleotideShuffle("AU", 1).Shuffles[0], Is.EqualTo("AU"));
    }
}
