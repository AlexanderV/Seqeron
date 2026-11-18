using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// B09 F46 (audit A1-6): the pairwise MCP tools align inputs above 25 000 000 DP cells in exact linear space automatically
/// (same output as the full-matrix aligner); <c>general_gap_align</c> rejects inputs a non-affine gap function cannot handle
/// with a clear error instead of an OutOfMemoryException.
/// </summary>
[TestFixture]
public class LinearSpaceAlignTests
{
    private static string Random(int seed, int length)
    {
        var rnd = new System.Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++) chars[i] = "ACGT"[rnd.Next(4)];
        return new string(chars);
    }

    private static string Mutate(int seed, string s)
    {
        var rnd = new System.Random(seed);
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            double r = rnd.NextDouble();
            if (r < 0.04) sb.Append("ACGT"[rnd.Next(4)]);
            else if (r < 0.06) { }
            else if (r < 0.08) { sb.Append(c); sb.Append("ACGT"[rnd.Next(4)]); }
            else sb.Append(c);
        }
        return sb.ToString();
    }

    [Test]
    public void GlobalAlignAffine_AboveThreshold_EqualsExplicitLinearSpaceAndKeepsTheContract()
    {
        string a = Random(11, 5300);
        string b = Mutate(12, a); // 5 300 x ~5 300 = 28 M cells > 25 M: automatic linear-space switch
        Assert.That((long)(a.Length + 1) * (b.Length + 1), Is.GreaterThan(SequenceAligner.LinearSpaceCellThreshold));

        var tool = AlignmentTools.GlobalAlignAffine(a, b, 2, -3, -5, -2);
        var library = SequenceAligner.GlobalAlignAffineLinearSpace(a, b, SequenceAligner.BlastDna, 1_000_000);
        Assert.Multiple(() =>
        {
            Assert.That(tool.Score, Is.EqualTo(library.Score));
            Assert.That(tool.AlignedSequence1, Is.EqualTo(library.AlignedSequence1));
            Assert.That(tool.AlignedSequence2, Is.EqualTo(library.AlignedSequence2));
            Assert.That(tool.AlignmentType, Is.EqualTo("Global"));
            Assert.That(tool.AlignedSequence1.Replace("-", ""), Is.EqualTo(a));
            Assert.That(tool.AlignedSequence2.Replace("-", ""), Is.EqualTo(b));
        });
    }

    [Test]
    public void LocalAlignAffine_AboveThreshold_ReturnsTheLocalAlignmentCoordinates()
    {
        string a = Random(21, 5300);
        string b = Mutate(22, a);
        var tool = AlignmentTools.LocalAlignAffine(a, b, 2, -3, -5, -2);
        var library = SequenceAligner.LocalAlignAffineLinearSpace(a, b, SequenceAligner.BlastDna, 1_000_000);
        Assert.Multiple(() =>
        {
            Assert.That(tool.Score, Is.EqualTo(library.Score));
            Assert.That(tool.AlignedSequence1, Is.EqualTo(library.AlignedSequence1));
            Assert.That((tool.StartPosition1, tool.StartPosition2, tool.EndPosition1, tool.EndPosition2),
                Is.EqualTo((library.StartPosition1, library.StartPosition2, library.EndPosition1, library.EndPosition2)));
            Assert.That(tool.AlignmentType, Is.EqualTo("Local"));
        });
    }

    [Test]
    public void GeneralGapAlign_NonAffineGapFunction_AboveCellLimit_IsRejectedWithAClearError()
    {
        string a = Random(31, 2100);
        string b = Random(32, 2100); // 4.4 M cells > 4 M
        var table = Enumerable.Range(1, 2100).Select(k => -(2 + 3 * Math.Log(k))).ToArray();
        var ex = Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign(a, b, table, "global", 2, -3));
        Assert.That(ex!.Message, Does.Contain("exceeds the limit"));
    }
}
