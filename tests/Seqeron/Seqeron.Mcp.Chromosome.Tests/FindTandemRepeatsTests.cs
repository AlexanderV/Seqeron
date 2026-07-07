using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_tandem_repeats</c>. GenomeAssemblyAnalyzer.FindTandemRepeats reports tandem arrays
/// with unit, copy number and purity. A pure "CAG" x 30 array is one repeat of unit "CAG", 30 copies.
/// </summary>
[TestFixture]
public class FindTandemRepeatsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void FindTandemRepeats_Schema_ValidatesCorrectly()
    {
        var seqs = new List<NamedSequence> { new("s", Rep("CAG", 30)) };
        Assert.DoesNotThrow(() => ChromosomeTools.FindTandemRepeats(seqs));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindTandemRepeats(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindTandemRepeats(seqs, minUnitLength: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindTandemRepeats(seqs, minUnitLength: 5, maxUnitLength: 2));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindTandemRepeats(seqs, minCopies: 0));
    }

    [Test]
    public void FindTandemRepeats_Binding_InvokesSuccessfully()
    {
        // 30 tandem copies of "CAG" (90 bp), purity 1.0.
        var result = ChromosomeTools.FindTandemRepeats(new List<NamedSequence> { new("s", Rep("CAG", 30)) }, 2, 50, 3);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var r = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceId, Is.EqualTo("s"));
            Assert.That(r.Unit, Is.EqualTo("CAG"));
            Assert.That(r.Copies, Is.EqualTo(30));
            Assert.That(r.Start, Is.EqualTo(0));
            Assert.That(r.End, Is.EqualTo(89));
            Assert.That(r.Purity, Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
