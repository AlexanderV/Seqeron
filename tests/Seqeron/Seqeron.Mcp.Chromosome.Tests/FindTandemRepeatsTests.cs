using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_tandem_repeats</c>. GenomeAssemblyAnalyzer.FindTandemRepeats runs Tandem Repeats Finder
/// (Benson 1999) and reports consensus unit, copy number and purity (match fraction between adjacent copies).
/// A pure "CAG" x 30 array is one repeat of unit "CAG", 30 copies.
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
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindTandemRepeats(seqs, minScore: 0));
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

    [Test]
    public void FindTandemRepeats_SubstitutedArray_ReportsTrfFractionalPurityAndMinScore()
    {
        // (ATG)x7 · ATC · (ATG)x7: TRF score 81, 40 of 42 adjacent-copy columns match.
        string seq = Rep("ATG", 7) + "ATC" + Rep("ATG", 7);
        var seqs = new List<NamedSequence> { new("s", seq) };

        var r = ChromosomeTools.FindTandemRepeats(seqs, 2, 50, 3).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That((r.Start, r.End, r.Unit), Is.EqualTo((0, 44, "ATG")));
            Assert.That(r.Copies, Is.EqualTo(15.0).Within(1e-9));
            Assert.That(r.Purity, Is.EqualTo(40.0 / 42.0).Within(1e-9));
        });

        // Minscore above the alignment score (81) rejects it.
        Assert.That(ChromosomeTools.FindTandemRepeats(seqs, 2, 50, 3, minScore: 82).Items, Is.Empty);
    }
}
