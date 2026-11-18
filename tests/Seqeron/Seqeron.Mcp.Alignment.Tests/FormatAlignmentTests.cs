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
        // EMBOSS needle 6.6.0 srspair: '|' identity, '.' mismatch, ' ' gap. "ACGT" vs "ACGA" -> "|||.".
        // A single block (width 80) is seq1 line, markup line, seq2 line, then a blank line.
        var r = AlignmentTools.FormatAlignment("ACGT", "ACGA", 80);
        Assert.That(r.Formatted, Is.EqualTo("ACGT\n|||.\nACGA\n\n"));
    }

    [Test]
    public void FormatAlignment_Matrix_MarksSimilarColumns()
    {
        // BLOSUM62: L/I +2 and A/S +1 are similar (':'), A/T 0 is a mismatch ('.'), K/R +2 similar, gap ' '.
        var r = AlignmentTools.FormatAlignment("LAAK-A", "ISTRAA", 80, "BLOSUM62");
        Assert.That(r.Formatted, Is.EqualTo("LAAK-A\n::.: |\nISTRAA\n\n"));
        // Scalar model with positive mismatch marks every substitution ':'.
        Assert.That(AlignmentTools.FormatAlignment("ACGT", "ACGA", 80, null, 1, 1).Formatted, Is.EqualTo("ACGT\n|||:\nACGA\n\n"));
    }
}
