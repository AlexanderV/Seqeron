using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CrisprSpecificityScoreTests
{
    // 20-nt guide with no internal GG/CC so the only PAM in the genome is the on-target's.
    private const string Guide = "ATATATATATATATATATAT";

    // Genome = guide + "AGG": the single forward NGG PAM is the trailing "AGG", whose
    // 20-nt upstream target equals the guide exactly (0 mismatches -> not an off-target).
    private const string Genome = "ATATATATATATATATATATAGG";

    [Test]
    public void CrisprSpecificityScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.crispr_specificity_score(Guide, Genome));
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score("", Genome));
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score(null!, Genome));
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score(Guide, ""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score(Guide, null!));
        // Wrong guide length for the system must throw from the underlying method.
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score("ACGT", Genome, CrisprSystemType.SpCas9));
    }

    [Test]
    public void CrisprSpecificityScore_Binding_InvokesSuccessfully()
    {
        // Only the on-target (0 mismatches) exists -> no off-targets -> specificity = 100.
        var result = MolToolsTools.crispr_specificity_score(Guide, Genome);
        Assert.That(result.Specificity, Is.EqualTo(100.0).Within(1e-9));
    }
}
