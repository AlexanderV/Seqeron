using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FindWithEditsTests
{
    [Test]
    public void FindWithEdits_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindWithEdits("ACGTACGT", "ACG", 1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindWithEdits("", "ACG", 1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindWithEdits("ACGT", null!, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FindWithEdits("ACGT", "ACG", -1));
    }

    [Test]
    public void FindWithEdits_Binding_InvokesSuccessfully()
    {
        // Levenshtein with sliding windows of length 3±1 = [2,4]. Pattern "AGT" in "ACGT":
        // window "ACGT" (edit 1, deletion), window "CGT" (edit 1, substitution),
        // window "GT" (edit 1, deletion) all qualify at maxEdits=1.
        var r = AlignmentTools.FindWithEdits("ACGT", "AGT", 1);
        Assert.Multiple(() =>
        {
            Assert.That(r.Items, Has.Length.EqualTo(3));
            Assert.That(r.Items[0].MatchedSequence, Is.EqualTo("ACGT"));
            Assert.That(r.Items[0].Distance, Is.EqualTo(1));
            Assert.That(r.Items[1].MatchedSequence, Is.EqualTo("CGT"));
            Assert.That(r.Items[1].Distance, Is.EqualTo(1));
        });
    }
}
