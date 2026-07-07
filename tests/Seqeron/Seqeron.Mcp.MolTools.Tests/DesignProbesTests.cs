using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignProbesTests
{
    // 81-nt ACGT-repeat, the Microarray target used by Seqeron.Genomics.Tests (PROBE-DESIGN-001).
    private const string Target =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    [Test]
    public void DesignProbes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_probes(Target));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes(Target, max_probes: 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes(Target, max_probes: -1));
    }

    [Test]
    public void DesignProbes_TooShort_ReturnsEmpty()
    {
        // Target shorter than the Microarray minimum length (50) yields no probes (algorithm yields break).
        var probes = MolToolsTools.design_probes("ACGTACGTACGT").Probes;
        Assert.That(probes, Is.Empty);
    }

    [Test]
    public void DesignProbes_Binding_InvokesSuccessfully()
    {
        // maxProbes limits the count exactly; the 81-nt ACGT repeat yields exactly maxProbes probes.
        var probes = MolToolsTools.design_probes(Target, max_probes: 3).Probes;

        Assert.Multiple(() =>
        {
            Assert.That(probes.Count, Is.EqualTo(3));
            // Sorted by score, descending.
            for (int i = 0; i + 1 < probes.Count; i++)
                Assert.That(probes[i].Score, Is.GreaterThanOrEqualTo(probes[i + 1].Score));
            foreach (var p in probes)
            {
                // Microarray length window is 50..60.
                Assert.That(p.Sequence.Length, Is.InRange(50, 60));
                // Probe sequence equals the target substring at [Start, End].
                Assert.That(p.Sequence, Is.EqualTo(Target.Substring(p.Start, p.End - p.Start + 1)));
                Assert.That(p.Tm, Is.GreaterThan(0));
            }
        });
    }
}
