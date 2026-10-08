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
        // B07 audit round 7, A7-2: a target shorter than probe_length (30 nt, default 60) is an argument error, not
        // InvalidOperationException "Sequence contains no elements".
        var ex = Assert.Throws<ArgumentException>(() => MolToolsTools.design_tiling_probes(new string('A', 30)));
        Assert.That(ex!.ParamName, Is.EqualTo("probe_length"));
    }

    [Test]
    public void DesignTilingProbes_Binding_InvokesSuccessfully()
    {
        // probeLength=50, overlap=10 -> step=40. Grid starts {0,40,80,120} end at 169; (208 − 50) mod 40 = 38 ≠ 0, so a
        // final window is anchored at the end (CATCH): start 208 − 50 = 158 → every one of the 208 positions covered
        // (B07 audit round 7, A7-2; previously 4 probes, coverage 170).
        var set = MolToolsTools.design_tiling_probes(Target, probe_length: 50, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(set.Probes.Count, Is.EqualTo(5));
            Assert.That(set.Probes.Select(p => p.Start), Is.EqualTo(new[] { 0, 40, 80, 120, 158 }));
            Assert.That(set.Probes[^1].End, Is.EqualTo(207));
            Assert.That(set.Coverage, Is.EqualTo(208));
            Assert.That(set.Probes.All(p => p.Type == ProbeDesigner.ProbeType.Tiling), Is.True);
            // MeanTm / TmRange consistency with the individual probes.
            Assert.That(set.MeanTm, Is.EqualTo(set.Probes.Average(p => p.Tm)).Within(1e-6));
            Assert.That(set.TmRange, Is.EqualTo(set.Probes.Max(p => p.Tm) - set.Probes.Min(p => p.Tm)).Within(1e-6));
        });
    }
}
