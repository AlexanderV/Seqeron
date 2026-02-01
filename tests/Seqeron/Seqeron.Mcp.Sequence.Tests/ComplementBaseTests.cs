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
}
