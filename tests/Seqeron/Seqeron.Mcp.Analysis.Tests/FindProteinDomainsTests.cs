using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_protein_domains</c> MCP tool.
/// Expected values from ProteinMotifFinder's own unit test
/// (ProteinMotifFinder_DomainPrediction_Tests, C2H2 zinc finger PS00028/PF00096 on
/// "AAAACAACAAALEEEEEEEEHAAAHAAAA" spanning 4..24), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindProteinDomainsTests
{
    [Test]
    public void FindProteinDomains_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindProteinDomains("AAAACAACAAALEEEEEEEEHAAAHAAAA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinDomains(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinDomains(null!));
    }

    [Test]
    public void FindProteinDomains_Binding_InvokesSuccessfully()
    {
        // C2H2 zinc finger domain: one hit spanning residues 4..24.
        var domains = AnalysisTools.FindProteinDomains("AAAACAACAAALEEEEEEEEHAAAHAAAA").Items;
        Assert.Multiple(() =>
        {
            Assert.That(domains, Has.Length.EqualTo(1));
            Assert.That(domains[0].Name, Is.EqualTo("Zinc Finger C2H2"));
            Assert.That(domains[0].Accession, Is.EqualTo("PF00096"));
            Assert.That(domains[0].Start, Is.EqualTo(4));
            Assert.That(domains[0].End, Is.EqualTo(24));
        });

        // Poly-alanine has no detectable domain.
        var none = AnalysisTools.FindProteinDomains("AAAAAAAA").Items;
        Assert.That(none, Is.Empty);
    }
}
