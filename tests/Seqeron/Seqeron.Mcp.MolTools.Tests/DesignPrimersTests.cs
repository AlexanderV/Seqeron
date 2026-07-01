using System;
using System.Text;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignPrimersTests
{
    // Standard template from Seqeron.Genomics.Tests PrimerDesigner_PrimerDesign_Tests:
    // 100 bp GAACTCGT-unit forward region + 50 bp poly-T target + 100 bp TCCGAAGT-unit
    // reverse region = 258 bp. Documented to yield a valid primer pair for target [100,150).
    private static string StandardTemplate()
    {
        var sb = new StringBuilder();
        while (sb.Length < 100) sb.Append("GAACTCGT");
        sb.Append(new string('T', 50));
        int revStart = sb.Length;
        while (sb.Length - revStart < 100) sb.Append("TCCGAAGT");
        return sb.ToString();
    }

    [Test]
    public void DesignPrimers_Schema_ValidatesCorrectly()
    {
        var t = StandardTemplate();
        Assert.DoesNotThrow(() => MolToolsTools.design_primers(t, 100, 150));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers("", 100, 150));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(null!, 100, 150));
        // Negative start.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, -1, 150));
        // End at/after template length.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, t.Length));
        // start >= end.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 150, 100));
    }

    [Test]
    public void DesignPrimers_Binding_InvokesSuccessfully()
    {
        var t = StandardTemplate();
        Assert.That(t.Length, Is.EqualTo(258));

        var r = MolToolsTools.design_primers(t, 100, 150);

        Assert.Multiple(() =>
        {
            // Documented behavior of PrimerDesigner.DesignPrimers on this fixture.
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.Forward, Is.Not.Null);
            Assert.That(r.Reverse, Is.Not.Null);

            // Forward is upstream of the target and ends at/before target start.
            Assert.That(r.Forward!.Position, Is.LessThan(100));
            Assert.That(r.Forward.Position + r.Forward.Length, Is.LessThanOrEqualTo(100));
            // Reverse is downstream of target end.
            Assert.That(r.Reverse!.Position, Is.GreaterThanOrEqualTo(150));

            // Exact selected pair (highest-scoring valid candidates).
            Assert.That(r.Forward.Position, Is.EqualTo(0));
            Assert.That(r.Forward.Length, Is.EqualTo(25));
            Assert.That(r.Reverse.Position, Is.EqualTo(155));
            Assert.That(r.Reverse.Length, Is.EqualTo(25));

            // Product size = reverse.Position + reverse.Length - forward.Position = 155 + 25 - 0.
            Assert.That(r.ProductSize, Is.EqualTo(180));
            Assert.That(r.ProductSize,
                Is.EqualTo(r.Reverse.Position + r.Reverse.Length - r.Forward.Position));

            // Both primers in the documented GC (40-60%) and Tm (57-63°C) windows.
            Assert.That(r.Forward.GcContent, Is.InRange(40.0, 60.0));
            Assert.That(r.Reverse.GcContent, Is.InRange(40.0, 60.0));
            Assert.That(r.Forward.MeltingTemperature, Is.InRange(57.0, 63.0));
            Assert.That(r.Reverse.MeltingTemperature, Is.InRange(57.0, 63.0));
            // Pair Tm difference within 5°C (compatibility rule).
            Assert.That(Math.Abs(r.Forward.MeltingTemperature - r.Reverse.MeltingTemperature),
                Is.LessThanOrEqualTo(5.0));
        });
    }
}
