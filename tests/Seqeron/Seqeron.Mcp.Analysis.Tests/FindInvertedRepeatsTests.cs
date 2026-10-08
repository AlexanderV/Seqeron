using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_inverted_repeats</c> MCP tool.
/// Expected values from the reverse-complement-arm definition on a GGGG-AAA-CCCC
/// hairpin (totalLength = 2*arm + loop), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindInvertedRepeatsTests
{
    [Test]
    public void FindInvertedRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 4, 50, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindInvertedRepeats("", 4, 50, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindInvertedRepeats(null!, 4, 50, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 1, 50, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 4, 50, -1));
    }

    [Test]
    public void FindInvertedRepeats_Binding_InvokesSuccessfully()
    {
        // GGGG (0) ... loop AAA ... CCCC (7): CCCC = revcomp(GGGG).
        var items = AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 4, 50, 3).Items;
        var hp = items.Single(i => i is { LeftArmStart: 0, ArmLength: 4 });
        Assert.Multiple(() =>
        {
            Assert.That(hp.RightArmStart, Is.EqualTo(7));
            Assert.That(hp.LoopLength, Is.EqualTo(3));
            Assert.That(hp.LeftArm, Is.EqualTo("GGGG"));
            Assert.That(hp.RightArm, Is.EqualTo("CCCC"));
            Assert.That(hp.Loop, Is.EqualTo("AAA"));
            Assert.That(hp.CanFormHairpin, Is.True);
            Assert.That(hp.TotalLength, Is.EqualTo(11));
        });

        // Poly-A tract has no reverse-complement arms.
        var none = AnalysisTools.FindInvertedRepeats("AAAAAAAA", 4, 50, 3).Items;
        Assert.That(none, Is.Empty);
    }

    [Test]
    public void FindInvertedRepeats_ReportsMaximalStemOnly()
    {
        // EMBOSS palindrome (-nummismatches 0 -overlap Y, minpallen 4, gaplimit 50): [(0,8,5)] —
        // the arm-4 sub-stems (1,8,4) and (0,9,4) lie inside the 5-bp stem and are not reported.
        var items = AnalysisTools.FindInvertedRepeats("GGGGGAAACCCCC", 4, 50, 3).Items;
        Assert.That(items.Select(i => (i.LeftArmStart, i.RightArmStart, i.ArmLength)),
            Is.EqualTo(new[] { (0, 8, 5) }));

        // A loop window with maxLoopLength < minLoopLength is rejected by the core.
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 4, 2, 3));
    }

    [Test]
    public void FindInvertedRepeats_OptionalMismatchesArmLengthWobble()
    {
        // EMBOSS palindrome -minpallen 4 -gaplimit 10 -nummismatches 1: gaattcagg / cttaactcc (1 interior mismatch).
        var mm = AnalysisTools.FindInvertedRepeats("GAATTCAGGAAAACCTCAATTC", 4, 10, 0, maxMismatches: 1).Items;
        // EMBOSS palindrome -maxpallen 5: the arm-10 G/C stem is neither reported nor split.
        var capped = AnalysisTools.FindInvertedRepeats("GGGGGGGGGGAAACCCCCCCCCC", 4, 10, 0, maxArmLength: 5).Items;
        // G-U wobble completes the arm-5 stem G-C, G.U, U-A, G-C, C-G.
        var wobble = AnalysisTools.FindInvertedRepeats("GGUGCAAAAGCAUC", 4, 10, 3, allowWobble: true).Items;

        Assert.Multiple(() =>
        {
            Assert.That(mm.Select(i => (i.LeftArmStart, i.RightArmStart, i.ArmLength, i.Mismatches)),
                Is.EqualTo(new[] { (0, 13, 9, 1) }));
            Assert.That(capped, Is.Empty);
            Assert.That(wobble.Select(i => (i.LeftArmStart, i.RightArmStart, i.ArmLength, i.Mismatches)),
                Is.EqualTo(new[] { (0, 9, 5, 0) }));
            Assert.That(AnalysisTools.FindInvertedRepeats("GGUGCAAAAGCAUC", 4, 10, 3).Items, Is.Empty);
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", maxMismatches: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindInvertedRepeats("GGGGAAACCCC", 4, maxArmLength: 3));
        });
    }
}
