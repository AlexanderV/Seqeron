using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignTilingProbesTests
{
    // 208-nt target from Seqeron.Genomics.Tests (PROBE-DESIGN-001, M9): A*100 + GCGCGCGC + T*100.
    private static readonly string Target = new string('A', 100) + "GCGCGCGC" + new string('T', 100);

    [Test]
    public void DesignTilingProbes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_tiling_probes(Target, 50, 10));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes("", 50, 10));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes(null!, 50, 10));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes(Target, 0, 10));
        // Overlap must be < probe length.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes(Target, 50, 50));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes(Target, 50, -1));
    }

    [Test]
    public void DesignTilingProbes_Binding_InvokesSuccessfully()
    {
        // probeLength=50, overlap=10 -> step=40. Starts {0,40,80,120}; coverage = 170 positions.
        var set = MolToolsTools.design_tiling_probes(Target, probe_length: 50, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(set.Probes.Count, Is.EqualTo(4));
            Assert.That(set.Probes.Select(p => p.Start), Is.EqualTo(new[] { 0, 40, 80, 120 }));
            Assert.That(set.Coverage, Is.EqualTo(170));
            Assert.That(set.Probes.All(p => p.Type == ProbeDesigner.ProbeType.Tiling), Is.True);
            // MeanTm / TmRange consistency with the individual probes.
            Assert.That(set.MeanTm, Is.EqualTo(set.Probes.Average(p => p.Tm)).Within(1e-6));
            Assert.That(set.TmRange, Is.EqualTo(set.Probes.Max(p => p.Tm) - set.Probes.Min(p => p.Tm)).Within(1e-6));
        });
    }
}
