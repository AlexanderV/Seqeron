using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplementBaseTests
{
    [Test]
    public void ComplementBase_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplementBase("A"));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplementBase(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplementBase(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplementBase("AT")); // Too long
    }

    [Test]
    public void ComplementBase_Binding_InvokesSuccessfully()
    {
        // DNA complements
        Assert.That(SequenceTools.ComplementBase("A").Complement, Is.EqualTo("T"));
        Assert.That(SequenceTools.ComplementBase("T").Complement, Is.EqualTo("A"));
        Assert.That(SequenceTools.ComplementBase("G").Complement, Is.EqualTo("C"));
        Assert.That(SequenceTools.ComplementBase("C").Complement, Is.EqualTo("G"));

        // RNA complement
        Assert.That(SequenceTools.ComplementBase("U").Complement, Is.EqualTo("A"));

        // Case insensitive
        Assert.That(SequenceTools.ComplementBase("a").Complement, Is.EqualTo("T"));
    }

    /// <summary>
    /// A1-12 / F25: default (DNA alphabet) equals Biopython 1.88 <c>complement</c>:
    /// complement("ACGTURYSWKMBDHVN") = "TGCAAYRSWMKVHDBN" (U → A, not "A↔U").
    /// </summary>
    [Test]
    public void ComplementBase_Default_MatchesBiopythonComplement()
    {
        const string input = "ACGTURYSWKMBDHVN";
        const string expected = "TGCAAYRSWMKVHDBN";
        for (int i = 0; i < input.Length; i++)
            Assert.That(SequenceTools.ComplementBase(input[i].ToString()).Complement, Is.EqualTo(expected[i].ToString()), $"base {input[i]}");
    }

    /// <summary>
    /// A1-12 / F25: rna=true equals Biopython 1.88 <c>complement_rna</c>:
    /// complement_rna("ACGTURYSWKMBDHVN") = "UGCAAYRSWMKVHDBN" (A → U, T → A); lower case → upper case; '-' passes through.
    /// </summary>
    [Test]
    public void ComplementBase_Rna_MatchesBiopythonComplementRna()
    {
        const string input = "ACGTURYSWKMBDHVN";
        const string expected = "UGCAAYRSWMKVHDBN";
        for (int i = 0; i < input.Length; i++)
            Assert.That(SequenceTools.ComplementBase(input[i].ToString(), rna: true).Complement, Is.EqualTo(expected[i].ToString()), $"base {input[i]}");
        Assert.Multiple(() =>
        {
            Assert.That(SequenceTools.ComplementBase("a", rna: true).Complement, Is.EqualTo("U"));
            Assert.That(SequenceTools.ComplementBase("-", rna: true).Complement, Is.EqualTo("-"));
            Assert.That(SequenceTools.ComplementBase("A", rna: true).Original, Is.EqualTo("A"));
        });
    }
}
