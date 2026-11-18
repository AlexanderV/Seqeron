using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_repeats</c> MCP tool.
/// Expected repeats derived from the "substring occurring >= 2 times" definition and
/// overlapping-occurrence counting, NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindRepeatsTests
{
    [Test]
    public void FindRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindRepeats("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeats("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeats(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeats("XYZ", 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindRepeats("AAAA", 0));
    }

    [Test]
    public void FindRepeats_Binding_InvokesSuccessfully()
    {
        // "AAAA" minLength 2 -> AA (0,1,2; count 3) and AAA (0,1; count 2).
        var homo = AnalysisTools.FindRepeats("AAAA", 2).Items;
        var byseq = homo.ToDictionary(r => r.Sequence);
        Assert.Multiple(() =>
        {
            Assert.That(byseq.ContainsKey("AA"), Is.True);
            Assert.That(byseq["AA"].Count, Is.EqualTo(3));
            Assert.That(byseq["AA"].Length, Is.EqualTo(2));
            Assert.That(byseq["AA"].Positions, Is.EquivalentTo(new[] { 0, 1, 2 }));
            Assert.That(byseq.ContainsKey("AAA"), Is.True);
            Assert.That(byseq["AAA"].Count, Is.EqualTo(2));
            Assert.That(byseq["AAA"].Positions, Is.EquivalentTo(new[] { 0, 1 }));
        });

        // "ATGATG" minLength 3 -> only ATG (0,3).
        var tandem = AnalysisTools.FindRepeats("ATGATG", 3).Items;
        Assert.Multiple(() =>
        {
            Assert.That(tandem, Has.Length.EqualTo(1));
            Assert.That(tandem[0].Sequence, Is.EqualTo("ATG"));
            Assert.That(tandem[0].Count, Is.EqualTo(2));
            Assert.That(tandem[0].Positions, Is.EquivalentTo(new[] { 0, 3 }));
        });
    }

    [Test]
    public void FindRepeats_MaximalOnly_ReturnsGusfieldMaximalRepeats()
    {
        // ACGTACGT: ACGT@{0,4} is the only left- and right-diverse repeat (brute force / generate_repeat_oracle.py).
        var maximal = AnalysisTools.FindRepeats("ACGTACGT", 1, maximalOnly: true).Items;
        Assert.Multiple(() =>
        {
            Assert.That(maximal, Has.Length.EqualTo(1));
            Assert.That(maximal[0].Sequence, Is.EqualTo("ACGT"));
            Assert.That(maximal[0].Positions, Is.EqualTo(new[] { 0, 4 }));
        });

        // Default (all repeated substrings) lists the prefixes and suffixes too.
        var all = AnalysisTools.FindRepeats("ACGTACGT", 1).Items.Select(r => r.Sequence).ToArray();
        Assert.That(all, Is.EqualTo(new[] { "A", "AC", "ACG", "ACGT", "C", "CG", "CGT", "G", "GT", "T" }));

        // ATCGATCGA: A@{0,4,8} is left-diverse (start, T, T) and right-diverse (T, T, end); ATCGA@{0,4} too.
        var atcga = AnalysisTools.FindRepeats("ATCGATCGA", 1, maximalOnly: true).Items;
        Assert.That(atcga.Select(r => r.Sequence), Is.EqualTo(new[] { "A", "ATCGA" }));
    }
}
