using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignAntisenseProbesTests
{
    // 81-nt mRNA-sense ACGT repeat; the designer works on its reverse complement.
    private const string Mrna =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    [Test]
    public void DesignAntisenseProbes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_antisense_probes(Mrna));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_antisense_probes(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_antisense_probes(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_antisense_probes(Mrna, max_probes: 0));
    }

    [Test]
    public void DesignAntisenseProbes_Binding_InvokesSuccessfully()
    {
        var probes = MolToolsTools.design_antisense_probes(Mrna, max_probes: 5).Probes;
        string antisense = DnaSequence.GetReverseComplementString(Mrna.ToUpperInvariant());

        Assert.Multiple(() =>
        {
            Assert.That(probes.Count, Is.EqualTo(5));
            // Every probe is tagged Antisense and is a substring of the reverse complement of the mRNA.
            Assert.That(probes.All(p => p.Type == ProbeDesigner.ProbeType.Antisense), Is.True);
            foreach (var p in probes)
                Assert.That(p.Sequence, Is.EqualTo(antisense.Substring(p.Start, p.End - p.Start + 1)));
        });
    }
}
