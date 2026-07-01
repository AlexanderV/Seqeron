using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AlignmentStatisticsTests
{
    [Test]
    public void AlignmentStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AlignmentStatistics("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("ACGT", null!));
        // Unequal aligned lengths are invalid.
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("ACGT", "ACG"));
    }

    [Test]
    public void AlignmentStatistics_Binding_InvokesSuccessfully()
    {
        // EMBOSS needle convention: denominator is the full alignment length (incl. gap columns).
        // "ACGT-A" vs "ACGTGA": 5 identical columns, 0 mismatch, 1 gap column, length 6.
        var r = AlignmentTools.AlignmentStatistics("ACGT-A", "ACGTGA");
        Assert.Multiple(() =>
        {
            Assert.That(r.Matches, Is.EqualTo(5));
            Assert.That(r.Mismatches, Is.EqualTo(0));
            Assert.That(r.Gaps, Is.EqualTo(1));
            Assert.That(r.AlignmentLength, Is.EqualTo(6));
            Assert.That(r.Identity, Is.EqualTo(500.0 / 6).Within(1e-9));
            Assert.That(r.Similarity, Is.EqualTo(500.0 / 6).Within(1e-9));
            Assert.That(r.GapPercent, Is.EqualTo(100.0 / 6).Within(1e-9));
        });

        // Perfect identity.
        var id = AlignmentTools.AlignmentStatistics("ACGT", "ACGT");
        Assert.Multiple(() =>
        {
            Assert.That(id.Matches, Is.EqualTo(4));
            Assert.That(id.Gaps, Is.EqualTo(0));
            Assert.That(id.Identity, Is.EqualTo(100.0));
            Assert.That(id.GapPercent, Is.EqualTo(0.0));
        });
    }
}
