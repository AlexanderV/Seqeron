using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FormatAlignmentTests
{
    [Test]
    public void FormatAlignment_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FormatAlignment("ACGT", "ACGA", 80));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FormatAlignment("", "ACGA", 80));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FormatAlignment("ACGT", null!, 80));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FormatAlignment("ACGT", "ACGA", 0));
    }

    [Test]
    public void FormatAlignment_Binding_InvokesSuccessfully()
    {
        // EMBOSS srspair legend: '|' identity, ' ' mismatch/gap. "ACGT" vs "ACGA" -> "||| ".
        // A single block (width 80) is seq1 line, markup line, seq2 line, then a blank line.
        var r = AlignmentTools.FormatAlignment("ACGT", "ACGA", 80);
        Assert.That(r.Formatted, Is.EqualTo("ACGT\n||| \nACGA\n\n"));
    }
}
