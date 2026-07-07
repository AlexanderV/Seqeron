using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class StickyCuttersTests
{
    [Test]
    public void StickyCutters_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.sticky_cutters());
        var result = MolToolsTools.sticky_cutters();
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Enzymes, Is.Not.Null);
    }

    [Test]
    public void StickyCutters_Binding_InvokesSuccessfully()
    {
        var result = MolToolsTools.sticky_cutters();
        var names = result.Enzymes.Select(e => e.Name).ToHashSet();

        Assert.Multiple(() =>
        {
            // 39 total enzymes - 10 blunt cutters = 29 sticky cutters.
            Assert.That(result.Enzymes, Has.Count.EqualTo(29));

            // Every returned enzyme is NOT blunt (staggered cut).
            Assert.That(result.Enzymes.All(e => !e.IsBluntEnd), Is.True);
            Assert.That(result.Enzymes.All(e => e.CutPositionForward != e.CutPositionReverse), Is.True);

            // Known sticky cutters present; known blunt cutters absent.
            Assert.That(names, Does.Contain("EcoRI"));
            Assert.That(names, Does.Contain("BamHI"));
            Assert.That(names, Does.Not.Contain("EcoRV"));
            Assert.That(names, Does.Not.Contain("SmaI"));
        });
    }
}
