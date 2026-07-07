using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignMolecularBeaconTests
{
    private const string Target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    [Test]
    public void DesignMolecularBeacon_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_molecular_beacon(Target, 20, 5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_molecular_beacon("", 20, 5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_molecular_beacon(null!, 20, 5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_molecular_beacon(Target, 0, 5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_molecular_beacon(Target, 20, 0));
    }

    [Test]
    public void DesignMolecularBeacon_ShortTarget_ReturnsNull()
    {
        // Target shorter than probe_length -> probe is null (not an error).
        var result = MolToolsTools.design_molecular_beacon("ACGT", probe_length: 20, stem_length: 5);
        Assert.That(result.Probe, Is.Null);
    }

    [Test]
    public void DesignMolecularBeacon_Binding_InvokesSuccessfully()
    {
        // stem_length=5 -> stem5 = "GG" + "CCC" = "GGCCC"; stem3 = revcomp("GGCCC") = "GGGCC".
        // Beacon = stem5(5) + loop(20) + stem3(5) = 30; loop is target[0..19].
        var result = MolToolsTools.design_molecular_beacon(Target, probe_length: 20, stem_length: 5);

        Assert.That(result.Probe, Is.Not.Null);
        var b = result.Probe!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(b.Type, Is.EqualTo(ProbeDesigner.ProbeType.MolecularBeacon));
            Assert.That(b.Sequence.Length, Is.EqualTo(30));
            Assert.That(b.Sequence.Substring(0, 5), Is.EqualTo("GGCCC"));
            Assert.That(b.Sequence.Substring(b.Sequence.Length - 5, 5), Is.EqualTo("GGGCC"));
            Assert.That(b.Sequence, Is.EqualTo(string.Concat("GGCCC", Target.AsSpan(0, 20), "GGGCC")));
            Assert.That(b.Start, Is.EqualTo(0));
            Assert.That(b.End, Is.EqualTo(19));
        });
    }
}
