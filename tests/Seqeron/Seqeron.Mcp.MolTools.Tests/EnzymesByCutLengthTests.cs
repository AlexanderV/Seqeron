using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EnzymesByCutLengthTests
{
    [Test]
    public void EnzymesByCutLength_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.enzymes_by_cut_length(6));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_by_cut_length(0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_by_cut_length(-4));
    }

    [Test]
    public void EnzymesByCutLength_Binding_InvokesSuccessfully()
    {
        var four = MolToolsTools.enzymes_by_cut_length(4);
        var six = MolToolsTools.enzymes_by_cut_length(6);
        var eight = MolToolsTools.enzymes_by_cut_length(8);
        var five = MolToolsTools.enzymes_by_cut_length(5);

        Assert.Multiple(() =>
        {
            // Built-in database: 9 four-cutters, 24 six-cutters, 5 eight-cutters.
            Assert.That(four.Enzymes, Has.Count.EqualTo(9));
            Assert.That(six.Enzymes, Has.Count.EqualTo(24));
            Assert.That(eight.Enzymes, Has.Count.EqualTo(5));
            // No 5-cutters in the database.
            Assert.That(five.Enzymes, Is.Empty);

            // Every returned enzyme really has the requested recognition length.
            Assert.That(four.Enzymes.All(e => e.RecognitionSequence.Length == 4), Is.True);
            Assert.That(six.Enzymes.All(e => e.RecognitionSequence.Length == 6), Is.True);

            // Spot-check membership.
            Assert.That(six.Enzymes.Select(e => e.Name), Does.Contain("EcoRI"));
            Assert.That(eight.Enzymes.Select(e => e.Name), Does.Contain("NotI"));
        });
    }
}
