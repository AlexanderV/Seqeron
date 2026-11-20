using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>A4-3 (B09 F153): the MCP EMBOSS tools read plain sequences like the programs: gap characters, digits, blanks and every other non-letter except <c>*</c> are dropped first; lengths and offsets of the formatted output count the cleaned sequence.</summary>
[TestFixture]
public class EmbossReaderGapCharactersTests
{
    [Test]
    public void NeedleAlign_OfGappedText_EqualsTheCleanedText_Score58AndLength11()
    {
        // needle 6.6.0: `AC-GTTAGC.TT*` vs `ACGTTAGCTT` prints Length 11 and Score 58.0
        var r = AlignmentTools.NeedleAlign("AC-GTTAGC.TT*", "ACGTTAGCTT");
        Assert.That(r.ExactScore, Is.EqualTo(58.0));
        Assert.That(r.AlignedSequence1, Is.EqualTo("ACGTTAGCTT*"));
        Assert.That(r.AlignedSequence2, Is.EqualTo("ACGTTAGCTT-"));
        Assert.That(r.EndPosition1, Is.EqualTo(11));
        var formatted = AlignmentTools.NeedleAlign("AC-GTTAGC.TT*", "ACGTTAGCTT", format: "pair").Formatted;
        Assert.That(formatted, Does.Contain("# Length: 11").And.Contain("# Score: 58.0"));
        Assert.That(AlignmentTools.WaterAlign("AC-GTTAGC.TT*", "ACGTTAGCTT").ExactScore, Is.EqualTo(58.0));
    }

    [Test]
    public void DotmatcherAndMerger_DropGapCharacters()
    {
        Assert.That(AlignmentTools.Dotmatcher("ACGTACGT-AACCGGTT", "ACGT.ACGTAACCGGTT", 6, 20).Segments,
            Is.EqualTo(AlignmentTools.Dotmatcher("ACGTACGTAACCGGTT", "ACGTACGTAACCGGTT", 6, 20).Segments));
        var gapped = AlignmentTools.Merger("acgtacgtac-gtacgtt1ggcc", "gtacgtt ggccaa.ttacg");
        var clean = AlignmentTools.Merger("acgtacgtacgtacgttggcc", "gtacgttggccaattacg");
        Assert.That(gapped.MergedSequence, Is.EqualTo(clean.MergedSequence));
    }
}
