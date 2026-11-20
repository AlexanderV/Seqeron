using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>merger (audit A3-13, B09 F139): EMBOSS merger 6.6.0. The tool delegates to SequenceAligner.Merger; the 450-row oracle comparison lives in the library tests.</summary>
[TestFixture]
public class MergerTests
{
    [Test]
    public void Merger_MergesADovetailOverlap_AndListsTheConflict()
    {
        // run on merger 6.6.0 (default EDNAFULL, gapopen 50, gapextend 5): the 3' end of a and the 5' end of b overlap over 20 bases with one mismatch
        const string a = "GATTACAGATTACACCGGTTAACGTACGTTAGC";
        const string b = "CCGGTTCACGTACGTTAGCAAATTTGGGCCCAAA";
        var r = AlignmentTools.Merger(a, b);
        Assert.Multiple(() =>
        {
            // merger prints "# Score: 86.0", the conflict "21 'A' 7 'c' 'A'" and the merged sequence below
            Assert.That(r.MergedSequence, Is.EqualTo("gattacagattacaccggttAacgtacgttagcaaatttgggcccaaa"));
            Assert.That(r.Conflicts.Select(c => (c.PositionA, c.BaseA, c.PositionB, c.BaseB, c.Using)), Is.EqualTo(new[] { (21, "A", 7, "c", "A") }));
            Assert.That(r.Score, Is.EqualTo(86.0));
        });
        Assert.Throws<ArgumentException>(() => AlignmentTools.Merger("", "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.Merger("ACGT", "ACGT", gapOpen: 200));
    }
}
