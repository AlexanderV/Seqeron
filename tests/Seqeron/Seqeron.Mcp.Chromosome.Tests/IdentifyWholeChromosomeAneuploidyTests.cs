using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;
using SourceCa = Seqeron.Genomics.Chromosome.ChromosomeAnalyzer;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>identify_whole_chromosome_aneuploidy</c>. ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy
/// flags a chromosome whose dominant copy-number bin covers >= minFraction of bins and is != 2, mapping
/// the copy number to an ISCN term (3 -> Trisomy).
/// </summary>
[TestFixture]
public class IdentifyWholeChromosomeAneuploidyTests
{
    private static SourceCa.CopyNumberState S(int cn) => new("chr1", 0, 999, cn, 0.0, 0.9);

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_Schema_ValidatesCorrectly()
    {
        var states = new List<SourceCa.CopyNumberState> { S(3) };
        Assert.DoesNotThrow(() => ChromosomeTools.IdentifyWholeChromosomeAneuploidy(states));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.IdentifyWholeChromosomeAneuploidy(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.IdentifyWholeChromosomeAneuploidy(states, minFraction: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.IdentifyWholeChromosomeAneuploidy(states, minFraction: 1.5));
    }

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_Binding_InvokesSuccessfully()
    {
        // All three bins copy number 3 -> Trisomy chr1.
        var states = new List<SourceCa.CopyNumberState> { S(3), S(3), S(3) };
        var result = ChromosomeTools.IdentifyWholeChromosomeAneuploidy(states, 0.8);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var a = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.Chromosome, Is.EqualTo("chr1"));
            Assert.That(a.CopyNumber, Is.EqualTo(3));
            Assert.That(a.Type, Is.EqualTo("Trisomy"));
        });
    }

    [Test]
    public void IdentifyWholeChromosomeAneuploidy_DiploidNotReported()
    {
        // Dominant copy number 2 is normal -> no aneuploidy reported.
        var states = new List<SourceCa.CopyNumberState> { S(2), S(2), S(2) };
        var result = ChromosomeTools.IdentifyWholeChromosomeAneuploidy(states, 0.8);
        Assert.That(result.Items, Is.Empty);
    }
}
