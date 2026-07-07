using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class ValidateProbeTests
{
    // From Seqeron.Genomics.Tests PROBE-VALID-001.
    private const string UniqueProbe = "ATCGATCGATCGATCGATCG";
    private static readonly string[] SingleMatchReference = { "NNNNNATCGATCGATCGATCGATCGNNNN" };

    [Test]
    public void ValidateProbe_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(null!, SingleMatchReference));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(UniqueProbe, null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, max_mismatches: -1));
    }

    [Test]
    public void ValidateProbe_UniqueProbe_SpecificityOne()
    {
        // 1 exact hit -> specificity 1.0.
        var v = MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, max_mismatches: 0);
        Assert.Multiple(() =>
        {
            Assert.That(v.OffTargetHits, Is.EqualTo(1));
            Assert.That(v.SpecificityScore, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    public void ValidateProbe_MultipleHits_SpecificityIsOneOverCount()
    {
        // 10-mer poly-A in a 34-mer poly-A: 34-10+1 = 25 exact-match positions -> specificity 1/25.
        var v = MolToolsTools.validate_probe("AAAAAAAAAA",
            new[] { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" }, max_mismatches: 0);
        Assert.Multiple(() =>
        {
            Assert.That(v.OffTargetHits, Is.EqualTo(25));
            Assert.That(v.SpecificityScore, Is.EqualTo(1.0 / 25).Within(1e-9));
        });
    }

    [Test]
    public void ValidateProbe_EmptyProbe_InvalidZeroSpecificity()
    {
        // Empty probe is a degenerate input: invalid, zero specificity, no hits.
        var v = MolToolsTools.validate_probe("", SingleMatchReference);
        Assert.Multiple(() =>
        {
            Assert.That(v.IsValid, Is.False);
            Assert.That(v.SpecificityScore, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(v.OffTargetHits, Is.EqualTo(0));
        });
    }
}
