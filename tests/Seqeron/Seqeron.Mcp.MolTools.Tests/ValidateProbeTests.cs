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

    [Test]
    public void ValidateProbe_ThermodynamicScreen_ReportsPrimer3NtthalTm()
    {
        // primer3-py 2.3.1 calc_homodimer / calc_hairpin at mv 50, dv 0, dntp 0, dna 50 (Primer3 probe conditions).
        var v = MolToolsTools.validate_probe("GCGCGCGCGCGCGCGCGCGC", Array.Empty<string>());
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.True);
            Assert.That(v.SelfDimerTm!.Value, Is.EqualTo(78.85652531616256).Within(1e-9));
            Assert.That(v.HairpinTm!.Value, Is.EqualTo(87.30265612393043).Within(1e-9));
            Assert.That(v.IsValid, Is.False);
        });
    }

    [Test]
    public void ValidateProbe_NonTargets_KaneCriteria()
    {
        // Kane et al. (2000): identity > 75 % or > 15 contiguous identical bases. Non-target = probe with a
        // substitution at every 5th position (Biopython local alignment: 40/50 identical → 80 %).
        const string probe = "TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA";
        const string nonTarget = "CCGCACCATGAGACTGTTTCTATGGCTCCTGTACCTCAAGTACATTTAGGCTTACGAGACAAATGCCAACCACATCGGCTTCGCACGTCT";
        var v = MolToolsTools.validate_probe(probe, new[] { probe }, non_target_sequences: new[] { nonTarget });
        Assert.Multiple(() =>
        {
            Assert.That(v.CrossHybridization, Has.Count.EqualTo(2));
            Assert.That(v.CrossHybridization[0].Identity, Is.EqualTo(0.8).Within(1e-12));
            Assert.That(v.CrossHybridization[0].CrossHybridizes, Is.True);
            Assert.That(v.IsValid, Is.False);
        });
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(probe, new[] { probe }, max_non_target_identity: 1.5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(probe, new[] { probe }, max_contiguous_match: -1));
    }
}
