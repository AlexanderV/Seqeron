using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RestrictionDigestTests
{
    [Test]
    public void RestrictionDigest_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", null!));
    }

    [Test]
    public void RestrictionDigest_Binding_InvokesSuccessfully()
    {
        // "AAAGAATTCAAA" (len 12), EcoRI cuts forward at position 4.
        // Fragments: [0,4)="AAAG", [4,12)="AATTCAAA".
        var fragments = MolToolsTools.restriction_digest("AAAGAATTCAAA", new[] { "EcoRI" }).Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(2));

            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAG"));
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(4));
            Assert.That(fragments[0].LeftEnzyme, Is.Null);
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));

            Assert.That(fragments[1].Sequence, Is.EqualTo("AATTCAAA"));
            Assert.That(fragments[1].StartPosition, Is.EqualTo(4));
            Assert.That(fragments[1].Length, Is.EqualTo(8));
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.Null);
        });

        // No cut site -> single fragment equal to the whole sequence.
        var uncut = MolToolsTools.restriction_digest("AAAAAAAA", new[] { "EcoRI" }).Fragments;
        Assert.Multiple(() =>
        {
            Assert.That(uncut, Has.Count.EqualTo(1));
            Assert.That(uncut[0].Sequence, Is.EqualTo("AAAAAAAA"));
        });
    }
}
