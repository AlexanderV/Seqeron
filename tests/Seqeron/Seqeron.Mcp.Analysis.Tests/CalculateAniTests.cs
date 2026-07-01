using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>calculate_ani</c> MCP tool.
/// Expected values taken from COMPGEN-ANI-001 / the algorithm's own unit tests
/// (ComparativeGenomics_CalculateANI_Tests, Goris et al. 2007), NOT the wrapper's output.
/// The wrapper exposes minFragmentIdentity (algorithm's minIdentity); tests set it
/// explicitly to reproduce the algorithm's documented cut-off behaviour.
/// </summary>
[TestFixture]
public class CalculateAniTests
{
    private const string Reference = "AAAACCCCGGGGTTTT";

    [Test]
    public void CalculateAni_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CalculateAni(Reference, Reference, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CalculateAni("", Reference, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CalculateAni(Reference, null!, 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CalculateAni(Reference, Reference, 0));
    }

    [Test]
    public void CalculateAni_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // M1 identical genomes: every 4-nt fragment perfect -> ANI 1.0.
            Assert.That(AnalysisTools.CalculateAni(Reference, Reference, 4, 0.30).Ani,
                Is.EqualTo(1.0).Within(1e-10));
            // M2 one substitution: last fragment 0.75 -> (1+1+1+0.75)/4 = 0.9375.
            Assert.That(AnalysisTools.CalculateAni("AAAACCCCGGGGTTTA", Reference, 4, 0.30).Ani,
                Is.EqualTo(0.9375).Within(1e-10));
            // M3 half-identity qualifying fragment (0.5 > 0.30) -> (1+1+1+0.5)/4 = 0.875.
            Assert.That(AnalysisTools.CalculateAni("AAAACCCCGGGGAATT", Reference, 4, 0.30).Ani,
                Is.EqualTo(0.875).Within(1e-10));
        });
    }

    [Test]
    public void CalculateAni_IdentityCutoff_ExcludesLowFragment()
    {
        // M4: "CGTC" identity 0.0 vs all-A ref is <= cutoff -> excluded; only "AAAA" kept -> 1.0.
        Assert.That(AnalysisTools.CalculateAni("AAAACGTC", "AAAAAAAA", 4, 0.30).Ani,
            Is.EqualTo(1.0).Within(1e-10));
    }
}
